//! A streaming FASTA/FASTQ reader that mirrors the behaviour of `kseq.h`.
//!
//! The parser is intentionally as permissive as `kseq` (and therefore as the C
//! tool): a record starts at any `>` or `@` found in the input, the sequence is
//! read line by line and the quality string of a FASTQ record is validated but
//! otherwise ignored.

use std::error::Error as StdError;
use std::fmt;
use std::io::{self, BufRead};

/// A FASTA/FASTQ record.
///
/// Buffers are reused between records, so the contents are only valid until the
/// next call to [`Reader::read`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Record {
    /// Sequence name, i.e. the header without the leading `>`/`@` and without
    /// the comment.
    pub name: Vec<u8>,
    /// Header comment, i.e. everything after the name on the header line.
    pub comment: Vec<u8>,
    /// Sequence bases, concatenated across lines and without whitespace.
    pub seq: Vec<u8>,
    /// Quality string, empty for FASTA records.
    pub qual: Vec<u8>,
}

impl Record {
    /// Returns the sequence name as a string, with invalid UTF-8 replaced.
    pub fn name_lossy(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.name)
    }

    /// Returns the sequence as a string, with invalid UTF-8 replaced.
    pub fn seq_lossy(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(&self.seq)
    }
}

/// Errors that can be reported while reading records.
#[derive(Debug)]
pub enum Error {
    /// The underlying stream failed.
    Io(io::Error),
    /// A FASTQ record ended before a complete quality string was read.
    TruncatedQuality,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::TruncatedQuality => write!(f, "truncated quality string in FASTQ input"),
        }
    }
}

impl StdError for Error {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::TruncatedQuality => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

/// A streaming FASTA/FASTQ reader.
pub struct Reader<R> {
    inner: R,
    last_char: Option<u8>,
    record: Record,
}

impl<R: BufRead> Reader<R> {
    /// Creates a reader on top of an already buffered stream.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            last_char: None,
            record: Record::default(),
        }
    }

    /// Reads the next record.
    ///
    /// Returns `Ok(None)` at end of input. The returned record borrows internal
    /// buffers and is invalidated by the next call.
    pub fn read(&mut self) -> Result<Option<&Record>, Error> {
        let Self {
            inner,
            last_char,
            record,
        } = self;

        // Jump to the next header line unless the previous record already
        // consumed its first character.
        if last_char.is_none() {
            loop {
                match getc(inner)? {
                    None => return Ok(None),
                    Some(c @ (b'>' | b'@')) => {
                        *last_char = Some(c);
                        break;
                    }
                    Some(_) => {}
                }
            }
        }
        // The header character has been consumed either way.
        *last_char = None;

        record.name.clear();
        record.comment.clear();
        record.seq.clear();
        record.qual.clear();

        // Header: name up to the first whitespace, then the rest of the line as
        // the comment.
        match read_until(inner, &mut record.name, is_space)? {
            None if record.name.is_empty() => return Ok(None),
            Some(b'\n') => {}
            Some(_) => {
                read_line(inner, &mut record.comment)?;
            }
            None => {}
        }

        // Sequence: lines are read whole, so a header marker only terminates the
        // sequence when it is the first character of a line.
        let terminator = loop {
            let c = match getc(inner)? {
                None => break None,
                Some(c) => c,
            };
            match c {
                b'>' | b'@' | b'+' => break Some(c),
                b'\n' => {}
                _ => {
                    record.seq.push(c);
                    read_line(inner, &mut record.seq)?;
                }
            }
        };
        if let Some(c @ (b'>' | b'@')) = terminator {
            *last_char = Some(c);
        }
        if terminator != Some(b'+') {
            return Ok(Some(record));
        }

        // FASTQ: skip the rest of the '+' line, then read the quality string.
        if !skip_line(inner)? {
            return Err(Error::TruncatedQuality);
        }
        while record.qual.len() < record.seq.len() {
            if read_line(inner, &mut record.qual)?.is_none() {
                break;
            }
        }
        if record.qual.len() != record.seq.len() {
            return Err(Error::TruncatedQuality);
        }
        Ok(Some(record))
    }
}

/// Matches `isspace()` in the C locale.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
}

/// Reads a single byte.
fn getc<R: BufRead>(reader: &mut R) -> io::Result<Option<u8>> {
    let byte = {
        let buf = reader.fill_buf()?;
        buf.first().copied()
    };
    if byte.is_some() {
        reader.consume(1);
    }
    Ok(byte)
}

/// Appends bytes to `out` until `stop` matches, consuming the matching byte.
///
/// Returns the terminating byte, or `None` at end of input.
fn read_until<R, F>(reader: &mut R, out: &mut Vec<u8>, stop: F) -> io::Result<Option<u8>>
where
    R: BufRead,
    F: Fn(u8) -> bool,
{
    loop {
        let (consumed, terminator) = {
            let buf = reader.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            match buf.iter().position(|&c| stop(c)) {
                Some(pos) => {
                    out.extend_from_slice(&buf[..pos]);
                    (pos + 1, Some(buf[pos]))
                }
                None => {
                    out.extend_from_slice(buf);
                    (buf.len(), None)
                }
            }
        };
        reader.consume(consumed);
        if terminator.is_some() {
            return Ok(terminator);
        }
    }
    Ok(None)
}

/// Appends the rest of the current line (without the newline) to `out`.
///
/// A trailing carriage return is dropped, like `ks_getuntil2()` does, so that
/// CRLF files are parsed correctly.
fn read_line<R: BufRead>(reader: &mut R, out: &mut Vec<u8>) -> io::Result<Option<u8>> {
    let terminator = read_until(reader, out, |c| c == b'\n')?;
    if out.len() > 1 && out.last() == Some(&b'\r') {
        out.pop();
    }
    Ok(terminator)
}

/// Consumes the rest of the current line, including the newline.
///
/// Returns `false` if the input ended before a newline was found.
fn skip_line<R: BufRead>(reader: &mut R) -> io::Result<bool> {
    loop {
        let (consumed, found) = {
            let buf = reader.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            match buf.iter().position(|&c| c == b'\n') {
                Some(pos) => (pos + 1, true),
                None => (buf.len(), false),
            }
        };
        reader.consume(consumed);
        if found {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(input: &str) -> Result<Vec<Record>, Error> {
        let mut reader = Reader::new(input.as_bytes());
        let mut records = Vec::new();
        while let Some(record) = reader.read()? {
            records.push(record.clone());
        }
        Ok(records)
    }

    #[test]
    fn fasta_with_comment() {
        let records = parse(">seq1 some comment\nACGT\nACGT\n").unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, b"seq1");
        assert_eq!(records[0].comment, b"some comment");
        assert_eq!(records[0].seq, b"ACGTACGT");
        assert!(records[0].qual.is_empty());
    }

    #[test]
    fn fasta_without_trailing_newline() {
        let records = parse(">seq1\nACGT").unwrap();
        assert_eq!(records[0].seq, b"ACGT");
    }

    #[test]
    fn fasta_with_crlf() {
        let records = parse(">seq1\r\nACGT\r\nACGT\r\n>seq2\r\nTTTT\r\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"seq1");
        assert_eq!(records[0].seq, b"ACGTACGT");
        assert_eq!(records[1].name, b"seq2");
        assert_eq!(records[1].seq, b"TTTT");
    }

    #[test]
    fn fastq_record() {
        let records = parse("@r1 comment\nACGT\n+\nIIII\n@r2\nTTTT\n+r2\n!!!!\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].comment, b"comment");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[0].qual, b"IIII");
        assert_eq!(records[1].name, b"r2");
        assert_eq!(records[1].seq, b"TTTT");
        assert_eq!(records[1].qual, b"!!!!");
    }

    #[test]
    fn fastq_multiline_quality() {
        let records = parse("@r1\nACGTACGT\n+\nIIII\nIIII\n").unwrap();
        assert_eq!(records[0].qual, b"IIIIIIII");
    }

    #[test]
    fn truncated_quality_is_an_error() {
        assert!(matches!(
            parse("@r1\nACGTACGT\n+\nIIII\n"),
            Err(Error::TruncatedQuality)
        ));
        assert!(matches!(
            parse("@r1\nACGT\n+\n"),
            Err(Error::TruncatedQuality)
        ));
    }

    #[test]
    fn leading_garbage_is_skipped() {
        let records = parse("garbage\n>seq1\nACGT\n").unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, b"seq1");
    }

    #[test]
    fn empty_input() {
        assert!(parse("").unwrap().is_empty());
        assert!(parse("ACGT\n").unwrap().is_empty());
    }

    #[test]
    fn empty_sequence_record() {
        let records = parse(">seq1\n>seq2\nACGT\n").unwrap();
        assert_eq!(records.len(), 2);
        assert!(records[0].seq.is_empty());
        assert_eq!(records[1].seq, b"ACGT");
    }

    #[test]
    fn header_marker_in_the_middle_of_a_line_is_part_of_the_sequence() {
        let records = parse(">seq1\nAC>GT\n").unwrap();
        assert_eq!(records[0].seq, b"AC>GT");
    }
}
