# sdustrs

[![crates.io](https://img.shields.io/crates/v/sdustrs.svg)][crates-io]

[Symmetric DUST][paper] 算法的 Rust 实现，用于在 DNA 序列中查找**低复杂度区域**
（low-complexity region）。本项目是 [lh3/sdust][sdust] C 版本的逐行移植：对同样的
参数和输入，输出与 C 版本**逐字节一致**，因此也与 NCBI 的 [dustmasker][dm] 一致
（组装 gap 除外），同时保持单遍扫描、低内存占用的特性。

[paper]: http://www.ncbi.nlm.nih.gov/pubmed/16796549
[dm]: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/
[sdust]: https://github.com/lh3/sdust
[crates-io]: https://crates.io/crates/sdustrs

## 构建

```sh
cargo build --release          # 产物：target/release/sdust
cargo install sdustrs          # 从 crates.io 安装，命令名仍为 sdust
cargo install --path .         # 可选：从本地源码安装到 ~/.cargo/bin
```

crate 名为 `sdustrs`（`sdust` 在 crates.io 上已被占用），但命令行程序与库内的
部分 API 沿用上游命名，安装后的可执行文件叫 `sdust`。

只需要 Rust 1.70 及以上版本；第三方依赖只有 `clap`（命令行解析）与 `flate2`（读取 gzip）。

## 命令行用法

```
Finds low-complexity regions in FASTA/FASTQ sequences and prints one
tab-separated `name start finish` line per region

Usage: sdust [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input FASTA/FASTQ file, or `-` for standard input

Options:
  -w, --window <INT>  Window size [default: 64]
  -t, --total <INT>   Score threshold [default: 20]
  -h, --help          Print help
  -V, --version       Print version

The input may be gzip compressed (including BGZF); use `-` to read from
standard input.
```

* `<in.fa>` 可以是 FASTA 或 FASTQ；`-` 表示从标准输入读取。
* 输入会自动识别 gzip / BGZF（按魔数判断，无需关心扩展名），
  因此 `sdust in.fa.gz` 与 `zcat in.fa.gz | sdust -` 等价。
* 序列可以跨多行，支持 CRLF 换行、大小写混合；FASTQ 的质量值会被校验但不会使用。

输出为每行一条区间，制表符分隔，坐标为 0-based 半开区间 `[start, finish)`，
与 BED 格式一致：

```
$ sdust test.fa
polyA_1	120	200
ca_rep_1	59	160
at_rep_1	60	181
acgt_rep_1	60	220
polyA25_1	100	125
n_gap_1	260	320
lower_polyA	60	122
```

需要注意，区间可能略微超出实际的低复杂度片段（例如向前多包含几个碱基，或越过
序列末尾），这是原算法按窗口整体报告区间的结果，C 版本行为完全相同。

## 库用法

```toml
[dependencies]
sdustrs = "0.2"
```

```rust
use sdustrs::{Masker, Params, sdust};

// 一次性调用
let regions = sdust(b"ACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGT", Params::default());
for r in &regions {
    println!("{}\t{}", r.start, r.finish);
}

// 频繁调用时复用缓冲区（对应 C API 的 sdust_buf_t），不产生额外堆分配
let mut masker = Masker::new();
for record in sequences {
    let regions = masker.mask(&record, Params::new(20, 64));
    // regions 借用自 masker，下一次调用前需要处理完
}
```

`Params::new(total, window)` 对应命令行的 `-t` 与 `-w`。参数非法时
`Params::validate()` 会返回错误说明；`Masker::mask` 本身不会 panic。

另外提供 `sdustrs::fastx` 模块，它是一个与 `kseq.h` 行为一致的流式 FASTA/FASTQ
解析器（`fastx::Reader`），读写缓冲区复用，可直接用于其他工具。

## 与 C 版本的一致性

移植过程中做了逐字节的差分验证：

* 仓库自带测试数据的 90 组参数组合（`-w ∈ {3,4,8,16,32,64}` × `-t ∈ {1,5,10,20,30}`
  × FASTA / FASTQ / gzip）输出完全一致，其中一部分已固化为 `tests/data/expected/`
  下的夹具，由 `cargo test` 端到端校验。
* 随机差分测试：随机生成包含串联重复、poly-A、`N` gap、大小写混合、多行 FASTQ
  的序列（长度 0–600 bp，7 种窗口 × 5 种阈值），与 C 版本逐字节对比，无差异。
* 5 Mbp 的随机序列、基因组风格序列（3% 重复）和高重复序列（20% 重复）输出全部一致。
* 边界情况：多成员（BGZF）gzip、CRLF 换行、多行质量值、缺少结尾换行、前后混入
  垃圾字符的输入，结果均与 C 版本一致。

## 与 C 版本的有意差异

1. **FASTQ 数据损坏时报错**：质量值长度不符时，Rust 版本向 stderr 输出
   `truncated quality string in FASTQ input` 并以退出码 1 结束；C 版本是静默停止
   （`kseq_read` 返回 -2，主循环直接退出）。
2. **参数校验**：`-w` 必须在 `3..=65535` 之间、`-t` 必须为正，否则报错退出。
   C 版本对非法值不检查（`-w 0` 会读越界，`-w` 过大时 `int` 计数溢出）。
3. **命令行是 C 版本的超集**：额外支持长选项、`-h/--help`、`-V/--version`；
   缺少输入文件或给出多余的位置参数会报错，而不是静默忽略。
4. **库接口是 Rust 风格**：区间用 `Interval { start, finish }` 结构体表示，
   而不是把两个 32 位数打包进 `u64`；`Masker` 取代 `sdust_buf_t`/`sdust_core`；
   序列长度上限为 `u32::MAX`（C 版本内部使用 `int`，上限 2^31-1）。

## 性能

在同一台机器上用 `-w 64 -t 20` 处理 5 Mbp 序列（单条记录，取三次最小值）：

| 输入 | C 版本 | Rust 版本 | 相对速度 |
| --- | --- | --- | --- |
| 纯随机序列 5 Mbp | 0.064 s | 0.057 s | 1.12× |
| 基因组风格 5 Mbp（3% 重复） | 0.524 s | 0.430 s | 1.22× |
| 高重复 1 Mbp（20% 重复） | 2.004 s | 1.678 s | 1.19× |
| 高重复 5 Mbp（20% 重复） | 15.366 s | 12.809 s | 1.20× |

高重复输入的时间主要花在 `perfect` 区间的插入排序与扫描上。Rust 版本在此处做了
两点优化：区间结构体压缩到 16 字节（与 C 的 `perf_intv_t` 相同），并用二分查找
（`partition_point`）定位插入位置、以「哨兵」初值消去内层循环的额外分支。

## 代码结构

| 文件 | 作用 |
| --- | --- |
| `src/lib.rs` | 核心算法：3-mer 滚动窗口（`Ring`/`Sliding`）、`find_perfect`、`save_masked_regions`、`Masker` |
| `src/fastx.rs` | FASTA/FASTQ 流式解析（`kseq.h` 的移植） |
| `src/main.rs` | 命令行解析、输入解压、输出 |
| `tests/cli.rs` | 端到端测试，与 C 版本生成的夹具对比 |
| `tests/data/` | 测试数据及期望输出 |

`src/lib.rs` 中的变量命名与注释尽量保持与 C 代码一一对应，便于对照阅读：
`wl`（窗口长度）、`l`（当前 A/C/G/T 连续片段长度）、`r`（相同三联体对数）、
`P`/`perfect`（当前窗口内的完美区间列表）等。

## 测试

```sh
cargo test              # 单元测试 + 端到端测试 + 文档测试
cargo clippy --all-targets
```

`tests/data/expected/*.out` 由原 C 版本生成，重新生成方式见 `tests/cli.rs` 顶部注释。

## 许可与致谢

算法与原始实现来自 Heng Li 的 [sdust][sdust]（MIT 许可），以及
Morgulis 等人发表的对称 DUST 论文。本移植同样以 MIT 许可发布。
