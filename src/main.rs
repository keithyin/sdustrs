//! Command line interface of `sdust`.
//!
//! ```text
//! Usage: sdust [-w 64] [-t 20] <in.fa>
//! ```

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

use clap::{error::ErrorKind, Parser};
use flate2::read::MultiGzDecoder;

use sdust::fastx;
use sdust::{Interval, Masker, Params};

/// Size of the buffers used for input and output.
const BUF_SIZE: usize = 64 * 1024;

const USAGE: &str = "Usage: sdust [-w 64] [-t 20] <in.fa>";

/// Command line arguments.
#[derive(Debug, Parser)]
#[command(
    name = "sdust",
    version,
    about = "Finds low-complexity regions in FASTA/FASTQ sequences and prints one \
             tab-separated `name start finish` line per region.",
    after_help = "The input may be gzip compressed (including BGZF); use `-` to read from \
                 standard input."
)]
struct Cli {
    /// Window size [64]
    #[arg(short = 'w', long = "window", value_name = "INT", default_value_t = Params::DEFAULT_WINDOW)]
    window: u32,

    /// Score threshold [20]
    #[arg(short = 't', long = "total", value_name = "INT", default_value_t = Params::DEFAULT_TOTAL)]
    total: u32,

    /// Input FASTA/FASTQ file, or `-` for standard input
    input: String,
}

impl Cli {
    fn params(&self) -> Params {
        Params::new(self.total, self.window)
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        // A closed pipe (`sdust big.fa | head`) is not an error.
        Err(Error::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("sdust: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), Error> {
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        // `--help` and `--version` are reported as errors by `clap`; they are
        // printed on standard output and are not failures.
        Err(e) if matches!(e.kind(), ErrorKind::DisplayHelp | ErrorKind::DisplayVersion) => {
            e.print().map_err(Error::Io)?;
            return Ok(());
        }
        // A parse error carries its own message; make sure a usage line is
        // always reported alongside it (`clap` omits it for some kinds).
        Err(e) => {
            let msg = e.to_string();
            return if msg.contains("Usage:") {
                Err(Error::Failure(msg))
            } else {
                Err(Error::Usage(msg))
            };
        }
    };
    cli.params()
        .validate()
        .map_err(|msg| Error::Usage(msg.to_string()))?;

    let input = open_input(&cli.input)?;
    let mut reader = fastx::Reader::new(input);
    let mut masker = Masker::new();
    let mut out = BufWriter::with_capacity(BUF_SIZE, io::stdout().lock());
    let mut line = Vec::with_capacity(128);

    while let Some(record) = reader.read()? {
        let regions = masker.mask(&record.seq, cli.params());
        for region in regions {
            write_record(&mut out, &mut line, &record.name, *region)?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes one `name<TAB>start<TAB>finish<LF>` line.
fn write_record(
    out: &mut impl Write,
    line: &mut Vec<u8>,
    name: &[u8],
    region: Interval,
) -> io::Result<()> {
    line.clear();
    line.extend_from_slice(name);
    line.push(b'\t');
    push_u32(line, region.start);
    line.push(b'\t');
    push_u32(line, region.finish);
    line.push(b'\n');
    out.write_all(line)
}

/// Appends the decimal representation of `value` to `buf`.
fn push_u32(buf: &mut Vec<u8>, mut value: u32) {
    let mut digits = [0u8; 10];
    let mut len = 0;
    loop {
        digits[len] = b'0' + (value % 10) as u8;
        len += 1;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    buf.extend(digits[..len].iter().rev());
}

/// An error reported by the CLI.
#[derive(Debug)]
enum Error {
    /// Bad parameters; the usage line is appended.
    Usage(String),
    /// Input/output failure.
    Io(io::Error),
    /// Malformed FASTA/FASTQ input.
    Fastx(fastx::Error),
    /// Any other failure, reported without the usage line.
    Failure(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Usage(msg) => write!(f, "{msg}\n{USAGE}"),
            Self::Io(e) => write!(f, "{e}"),
            Self::Fastx(e) => write!(f, "{e}"),
            Self::Failure(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<fastx::Error> for Error {
    fn from(e: fastx::Error) -> Self {
        match e {
            fastx::Error::Io(e) => Self::Io(e),
            other => Self::Fastx(other),
        }
    }
}

/// Opens `path` (`-` for standard input), transparently decompressing gzip.
fn open_input(path: &str) -> Result<Box<dyn BufRead>, Error> {
    let raw: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin().lock())
    } else {
        Box::new(
            File::open(path)
                .map_err(|e| Error::Failure(format!("cannot open '{}': {}", path, e)))?,
        )
    };
    let mut buffered = BufReader::with_capacity(BUF_SIZE, raw);
    let is_gzip = {
        let head = buffered.fill_buf()?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(
            BUF_SIZE,
            MultiGzDecoder::new(buffered),
        )))
    } else {
        Ok(Box::new(buffered))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn itoa_matches_std() {
        for value in [0u32, 1, 9, 10, 99, 100, 1234, 65_535, u32::MAX] {
            let mut buf = Vec::new();
            push_u32(&mut buf, value);
            assert_eq!(String::from_utf8(buf).unwrap(), value.to_string());
        }
    }

    /// Parses `list` as if the binary had been invoked with those arguments.
    fn try_args(list: &[&str]) -> Result<Cli, clap::Error> {
        Cli::try_parse_from(std::iter::once("sdust").chain(list.iter().copied()))
    }

    #[test]
    fn defaults() {
        let cli = try_args(&["in.fa"]).unwrap();
        assert_eq!(cli.params(), Params::default());
        assert_eq!(cli.input, "in.fa");
    }

    #[test]
    fn attached_and_separated_values() {
        for list in [
            &["-w", "32", "-t", "10", "in.fa"][..],
            &["-w32", "-t10", "in.fa"],
        ] {
            let cli = try_args(list).unwrap();
            assert_eq!(cli.params(), Params::new(10, 32));
            assert_eq!(cli.input, "in.fa");
        }
    }

    #[test]
    fn long_options_and_ordering() {
        let cli = try_args(&["in.fa", "--window=32", "--total", "10"]).unwrap();
        assert_eq!(cli.params(), Params::new(10, 32));
        assert_eq!(cli.input, "in.fa");
    }

    #[test]
    fn dash_is_stdin() {
        let cli = try_args(&["-"]).unwrap();
        assert_eq!(cli.input, "-");
    }

    #[test]
    fn help_and_version() {
        assert!(matches!(
            try_args(&["-h"]).unwrap_err().kind(),
            ErrorKind::DisplayHelp
        ));
        assert!(matches!(
            try_args(&["--help"]).unwrap_err().kind(),
            ErrorKind::DisplayHelp
        ));
        assert!(matches!(
            try_args(&["-V"]).unwrap_err().kind(),
            ErrorKind::DisplayVersion
        ));
    }

    #[test]
    fn bad_arguments_are_rejected() {
        for list in [
            &["-x"][..],
            &["-w"],
            &["-w", "abc", "in.fa"],
            &["-w", "-3", "in.fa"],
            &["a.fa", "b.fa"],
            &[],
        ] {
            assert!(try_args(list).is_err(), "args: {list:?} should fail");
        }
    }

}
