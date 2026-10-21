//! End-to-end tests of the `sdust` binary.
//!
//! The expected outputs in `tests/data/expected/` were produced by the original
//! C implementation (https://github.com/lh3/sdust), e.g.
//!
//! ```text
//! ./sdust test.fa > expected/test.fa.default.out
//! ./sdust -w 4 -t 5 test.fa > expected/test.fa.w4t5.out
//! ```
//!
//! so these tests double as a consistency check against the reference.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_sdust");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn expected(name: &str) -> Vec<u8> {
    std::fs::read(data(&format!("expected/{name}"))).expect("expected output")
}

/// Runs the binary with `args`, feeding `input` on standard input.
fn run(args: &[&str], input: Option<&[u8]>) -> Output {
    let mut child = Command::new(BIN)
        .args(args)
        .stdin(if input.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to run sdust");
    if let Some(input) = input {
        child
            .stdin
            .take()
            .expect("stdin")
            .write_all(input)
            .expect("writing stdin");
    }
    child.wait_with_output().expect("waiting for sdust")
}

#[test]
fn default_parameters_on_fasta() {
    let out = run(&[data("test.fa").to_str().unwrap()], None);
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fa.default.out"));
}

#[test]
fn default_parameters_on_fastq() {
    let out = run(&[data("test.fq").to_str().unwrap()], None);
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fq.default.out"));
}

#[test]
fn default_parameters_on_gzip() {
    let out = run(&[data("test.fa.gz").to_str().unwrap()], None);
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fa.gz.default.out"));
}

#[test]
fn custom_parameters_match_the_c_implementation() {
    let fa = data("test.fa");
    let fa = fa.to_str().unwrap();
    let fq = data("test.fq");
    let fq = fq.to_str().unwrap();
    let cases: &[(&[&str], &str)] = &[
        (&["-w", "4", "-t", "5", fa], "test.fa.w4t5.out"),
        (&["-w", "3", "-t", "1", fa], "test.fa.w3t1.out"),
        (&["-w", "16", "-t", "30", fa], "test.fa.w16t30.out"),
        (&["-w", "100", "-t", "25", fa], "test.fa.w100t25.out"),
        (&["-w", "8", "-t", "5", fq], "test.fq.w8t5.out"),
    ];
    for (args, name) in cases {
        let out = run(args, None);
        assert!(out.status.success(), "args: {args:?}");
        assert_eq!(out.stdout, expected(name), "args: {args:?}");
    }
}

#[test]
fn option_forms_are_equivalent() {
    let fa = data("test.fa");
    let fa = fa.to_str().unwrap();
    let reference = expected("test.fa.w16t30.out");
    for args in [
        vec!["-w", "16", "-t", "30", fa],
        vec!["-w16", "-t30", fa],
        vec!["--window", "16", "--total=30", fa],
        vec![fa, "-t", "30", "-w", "16"],
    ] {
        let out = run(&args, None);
        assert!(out.status.success(), "args: {args:?}");
        assert_eq!(out.stdout, reference, "args: {args:?}");
    }
}

#[test]
fn reads_plain_and_gzip_from_stdin() {
    let plain = std::fs::read(data("test.fa")).unwrap();
    let out = run(&["-"], Some(&plain));
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fa.default.out"));

    let out = run(&["-w", "16", "-t", "30", "-"], Some(&plain));
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fa.w16t30.out"));

    let gzipped = std::fs::read(data("test.fa.gz")).unwrap();
    let out = run(&["-"], Some(&gzipped));
    assert!(out.status.success());
    assert_eq!(out.stdout, expected("test.fa.gz.default.out"));
}

#[test]
fn help_and_version() {
    let out = run(&["--help"], None);
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("Usage: sdust"), "{text}");
    assert!(text.contains("--window"), "{text}");

    let out = run(&["-V"], None);
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).starts_with("sdust "));
}

#[test]
fn bad_usage_fails_with_a_message() {
    for args in [
        vec!["--nonsense"],
        vec!["-w"],
        vec!["-w", "zero"],
        vec!["-w", "0"],
        vec!["-t", "0"],
        vec!["-w", "2"],
        vec!["-w", "70000"],
        vec!["a.fa", "b.fa"],
        vec![],
    ] {
        let out = run(&args, None);
        assert!(!out.status.success(), "args: {args:?} should fail");
        let text = String::from_utf8_lossy(&out.stderr);
        assert!(text.starts_with("sdust: "), "args: {args:?}: {text}");
        assert!(text.contains("Usage: sdust"), "args: {args:?}: {text}");
    }
}

#[test]
fn missing_input_file_fails() {
    let out = run(&["/nonexistent/file.fa"], None);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot open"));
}

#[test]
fn truncated_fastq_is_reported() {
    let out = run(&["-"], Some(b"@r1\nACGTACGT\n+\nIIII\n"));
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("truncated quality"));
}

#[test]
fn empty_input_produces_no_output() {
    let out = run(&["-"], Some(b""));
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
}
