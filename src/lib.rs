//! Symmetric DUST: detection of low-complexity regions in DNA sequences.
//!
//! This crate is a Rust port of the C tool [`sdust`](https://github.com/lh3/sdust),
//! which is itself a reimplementation of NCBI's symmetric DUST algorithm
//! (Morgulis et al., 2006). For the same parameters it produces exactly the same
//! intervals as the C implementation, and therefore the same output as
//! `dustmasker` except in assembly gaps.
//!
//! # Quick start
//!
//! ```
//! use sdust::{Params, sdust};
//!
//! let seq = b"ACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGT";
//! let regions = sdust(seq, Params::default());
//! for r in &regions {
//!     println!("{}\t{}", r.start, r.finish);
//! }
//! ```
//!
//! When masking many sequences, reuse a [`Masker`] to avoid re-allocating the
//! internal buffers (this is the counterpart of `sdust_buf_t` in the C API):
//!
//! ```
//! use sdust::{Masker, Params};
//!
//! let mut masker = Masker::new();
//! for seq in [b"ACGTACGT".as_slice(), b"AAAAAAAAAAAAAAAAAAAAAAAA".as_slice()] {
//!     let regions = masker.mask(seq, Params::default()).to_vec();
//!     assert!(regions.is_empty() || regions[0].start < seq.len() as u32);
//! }
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod fastx;

/// Length of the words the algorithm is built on (3-mers).
const WORD_LEN: u32 = 3;
/// Number of distinct words (`4^3`).
const WORD_TOT: usize = 1 << (WORD_LEN << 1);
/// Mask to keep the low `2 * WLEN` bits of the rolling word.
const WORD_MASK: u32 = WORD_TOT as u32 - 1;

/// Maps a byte to a 2-bit nucleotide code: `A=0`, `C=1`, `G=2`, `T=3`, anything
/// else (including `N` and lowercase non-`acgt` bytes) to `4`.
///
/// Uppercase and lowercase `A/C/G/T` are equivalent, exactly like the table in
/// the C implementation.
static NT4: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// A half-open interval `[start, finish)` of a low-complexity region.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Interval {
    /// Start coordinate (0-based, inclusive).
    pub start: u32,
    /// End coordinate (0-based, exclusive).
    pub finish: u32,
}

impl Interval {
    /// Creates an interval.
    pub fn new(start: u32, finish: u32) -> Self {
        Self { start, finish }
    }

    /// Length of the interval in bases.
    pub fn len(&self) -> u32 {
        self.finish.saturating_sub(self.start)
    }

    /// Returns `true` if the interval is empty.
    pub fn is_empty(&self) -> bool {
        self.finish <= self.start
    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}", self.start, self.finish)
    }
}

/// Parameters of the symmetric DUST algorithm.
///
/// The defaults match the command line defaults of the C tool.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Params {
    /// `T`: the score threshold. A window of `n` triplet positions is considered
    /// low-complexity when the number of repeated triplet *pairs* `r` it contains
    /// satisfies `10 * r > T * n`. Default: 20.
    pub total: u32,
    /// `W`: the size of the sliding window in bases. Default: 64.
    pub window: u32,
}

impl Params {
    /// Default score threshold (`-t`).
    pub const DEFAULT_TOTAL: u32 = 20;
    /// Default window size (`-w`).
    pub const DEFAULT_WINDOW: u32 = 64;
    /// Largest supported window size, bounded by the 32-bit pair counters.
    pub const MAX_WINDOW: u32 = 65535;

    /// Creates parameters from a score threshold and a window size.
    pub const fn new(total: u32, window: u32) -> Self {
        Self { total, window }
    }

    /// Checks that the parameters are usable.
    ///
    /// The window has to be at least 3 bases long for the algorithm to be well
    /// defined (the C implementation silently misbehaves below that), and it has
    /// to stay below [`Params::MAX_WINDOW`] so that the pair counts of a window
    /// fit in 32 bits (the C implementation silently overflows instead).
    pub fn validate(&self) -> Result<(), &'static str> {
        if self.total == 0 {
            return Err("the score threshold (-t) must be positive");
        }
        if self.window < WORD_LEN {
            return Err("the window size (-w) must be at least 3");
        }
        if self.window > Self::MAX_WINDOW {
            return Err("the window size (-w) must be at most 65535");
        }
        Ok(())
    }
}

impl Default for Params {
    fn default() -> Self {
        Self::new(Self::DEFAULT_TOTAL, Self::DEFAULT_WINDOW)
    }
}

/// A FIFO of 3-mer codes with O(1) random access by position.
///
/// This is the equivalent of `kdq_t(int)` in the C code: a ring buffer that
/// grows by doubling and keeps its contents contiguous in logical order.
#[derive(Debug)]
struct Ring {
    buf: Vec<u8>,
    head: usize,
    len: usize,
}

impl Ring {
    fn new() -> Self {
        Self {
            buf: vec![0; 16],
            head: 0,
            len: 0,
        }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn clear(&mut self) {
        self.head = 0;
        self.len = 0;
    }

    /// Returns the `i`-th element, counting from the oldest one.
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    fn get(&self, i: usize) -> u8 {
        assert!(i < self.len, "ring index out of bounds");
        self.buf[(self.head + i) & (self.buf.len() - 1)]
    }

    fn push(&mut self, value: u8) {
        if self.len == self.buf.len() {
            self.grow();
        }
        let mask = self.buf.len() - 1;
        self.buf[(self.head + self.len) & mask] = value;
        self.len += 1;
    }

    /// Removes and returns the oldest element.
    fn pop_front(&mut self) -> Option<u8> {
        if self.len == 0 {
            return None;
        }
        let value = self.buf[self.head];
        self.head = (self.head + 1) & (self.buf.len() - 1);
        self.len -= 1;
        Some(value)
    }

    fn grow(&mut self) {
        let mut buf = vec![0; self.buf.len() * 2];
        let mask = self.buf.len() - 1;
        for (i, slot) in buf.iter_mut().enumerate().take(self.len) {
            *slot = self.buf[(self.head + i) & mask];
        }
        self.buf = buf;
        self.head = 0;
    }
}

/// The sliding window of 3-mers plus its incrementally maintained statistics.
///
/// Two nested windows are tracked, mirroring the C code:
///
/// * the ring buffer, which keeps up to `W - 2` words and grows one word per
///   base until it is full;
/// * the *active* window of `len` words, which is shorter than the ring when
///   the most recent word is over-represented.
///
/// For both, the number of *pairs* of identical 3-mers is maintained, since the
/// symmetric DUST score is based on pair counts rather than on word counts.
#[derive(Debug)]
struct Sliding {
    ring: Ring,
    /// Occurrences of each 3-mer in the ring.
    ring_count: [u32; WORD_TOT],
    /// Occurrences of each 3-mer in the active window.
    win_count: [u32; WORD_TOT],
    /// Pairs of identical 3-mers in the ring.
    ring_pairs: i64,
    /// Pairs of identical 3-mers in the active window.
    win_pairs: i64,
    /// Length of the active window, in words.
    len: i64,
}

impl Sliding {
    fn new() -> Self {
        Self {
            ring: Ring::new(),
            ring_count: [0; WORD_TOT],
            win_count: [0; WORD_TOT],
            ring_pairs: 0,
            win_pairs: 0,
            len: 0,
        }
    }

    fn reset(&mut self) {
        self.ring.clear();
        self.ring_count = [0; WORD_TOT];
        self.win_count = [0; WORD_TOT];
        self.ring_pairs = 0;
        self.win_pairs = 0;
        self.len = 0;
    }

    fn ring_len(&self) -> i64 {
        self.ring.len() as i64
    }

    /// Appends one 3-mer, sliding the window forward by one base.
    ///
    /// `ring_pairs` and `win_pairs` are the `rw` and `rv` counters of the C
    /// implementation: they count pairs of identical 3-mers, and are updated
    /// incrementally so that no window is ever rescanned.
    fn push(&mut self, word: u8, total: i64, window: i64) {
        let word = word as usize;
        // Drop the oldest word once the ring holds a full window.
        if !self.ring.is_empty() && self.ring_len() > window - WORD_LEN as i64 {
            let oldest = self.ring.pop_front().expect("ring is not empty") as usize;
            self.ring_count[oldest] -= 1;
            self.ring_pairs -= self.ring_count[oldest] as i64;
            if self.len > self.ring_len() {
                self.len -= 1;
                self.win_count[oldest] -= 1;
                self.win_pairs -= self.win_count[oldest] as i64;
            }
        }
        self.ring.push(word as u8);
        self.len += 1;
        self.ring_pairs += self.ring_count[word] as i64;
        self.ring_count[word] += 1;
        self.win_pairs += self.win_count[word] as i64;
        self.win_count[word] += 1;

        // The new word is over-represented: shrink the window from the left
        // until that word has been consumed.
        if self.win_count[word] as i64 * 10 > total * 2 {
            loop {
                let offset = self.ring_len() - self.len;
                debug_assert!(offset >= 0);
                let oldest = self.ring.get(offset as usize) as usize;
                self.win_count[oldest] -= 1;
                self.win_pairs -= self.win_count[oldest] as i64;
                self.len -= 1;
                if oldest == word {
                    break;
                }
            }
        }
    }
}

/// A perfect (maximal-scoring) interval found inside the current window.
///
/// Intervals are ordered by descending `start` and then by ascending `finish`.
///
/// The fields are kept narrow (16 bytes in total) because this list is the hot
/// spot of the algorithm: inserting an interval shifts the elements to its
/// right.
#[derive(Clone, Copy, Debug)]
struct Perfect {
    start: u32,
    finish: u32,
    /// Number of triplet pairs counted so far.
    pairs: u32,
    /// Number of triplet positions covered so far.
    len: u32,
}

/// Reusable buffers for masking an unbounded number of sequences.
///
/// This is the counterpart of `sdust_buf_t` in the C API: allocate one masker
/// and call [`Masker::mask`] for every sequence to avoid repeated allocations.
#[derive(Debug)]
pub struct Masker {
    sliding: Sliding,
    /// Perfect intervals of the current window, sorted by descending `start`
    /// and then by ascending `finish`.
    perfect: Vec<Perfect>,
    /// Low-complexity intervals found so far, sorted and merged.
    result: Vec<Interval>,
}

impl Masker {
    /// Creates an empty masker.
    pub fn new() -> Self {
        Self {
            sliding: Sliding::new(),
            perfect: Vec::new(),
            result: Vec::new(),
        }
    }

    /// Drops all state, keeping the allocated memory.
    pub fn reset(&mut self) {
        self.sliding.reset();
        self.perfect.clear();
        self.result.clear();
    }

    /// Finds the low-complexity intervals of `seq`.
    ///
    /// `seq` is interpreted as a nucleotide string: `A/C/G/T` in either case
    /// are valid bases, every other byte (notably `N`) breaks the sequence into
    /// independent pieces. The returned intervals are sorted, merged when they
    /// overlap or are adjacent, and borrowed from the masker, so they have to be
    /// consumed (or copied) before the next call.
    ///
    /// The parameters are not validated; call [`Params::validate`] first if they
    /// come from an untrusted source. Invalid parameters never panic, but the
    /// output is then unspecified.
    ///
    /// Coordinates are 32 bits wide, so sequences must be shorter than
    /// `u32::MAX` bases (longer than the C implementation supports).
    pub fn mask(&mut self, seq: &[u8], params: Params) -> &[Interval] {
        self.reset();
        let total = params.total as i64;
        let window = params.window as i64;
        // `perfect` never holds more intervals than there are words in a window.
        self.perfect.reserve(params.window as usize);

        let mut word: u32 = 0;
        // Length of the current contiguous A/C/G/T stretch.
        let mut length_of_cur_contiguous_stretch: i64 = 0;

        // The iteration runs one past the end so that the trailing piece of
        // sequence is flushed (a sentinel non-ACGT byte).
        for cur_base_position in 0..=seq.len() {
            let code = if cur_base_position < seq.len() {
                NT4[seq[cur_base_position] as usize]
            } else {
                4
            } as u32;
            if code < 4 {
                length_of_cur_contiguous_stretch += 1;
                word = (word << 2 | code) & WORD_MASK;
                if length_of_cur_contiguous_stretch >= WORD_LEN as i64 {
                    // Start of the current window.
                    let start_of_current_window = (length_of_cur_contiguous_stretch - window)
                        .max(0)
                        + (cur_base_position as i64 + 1 - length_of_cur_contiguous_stretch);

                    save_masked_regions(
                        &mut self.result,
                        &mut self.perfect,
                        start_of_current_window,
                    );

                    self.sliding.push(word as u8, total, window);
                    // Only scan when the window holds enough repeated words to
                    // possibly pass the threshold.
                    // ring_pairs / len > total / 10
                    if self.sliding.ring_pairs * 10 > self.sliding.len * total {
                        find_perfect(
                            &mut self.perfect,
                            &self.sliding,
                            total,
                            start_of_current_window,
                        );
                    }
                }
            } else {
                // `N` or the end of the sequence: flush the perfect intervals
                // that are still buffered, scanning the window to the right.
                let mut start_of_current_window = (length_of_cur_contiguous_stretch - window + 1)
                    .max(0)
                    + (cur_base_position as i64 + 1 - length_of_cur_contiguous_stretch);
                while !self.perfect.is_empty() {
                    save_masked_regions(
                        &mut self.result,
                        &mut self.perfect,
                        start_of_current_window,
                    );
                    start_of_current_window += 1;
                }
                length_of_cur_contiguous_stretch = 0;
                word = 0;
            }
        }
        &self.result
    }
}

impl Default for Masker {
    fn default() -> Self {
        Self::new()
    }
}

/// Saves the leftmost perfect interval holding the highest score.
///
/// Only intervals that start before `start` are eligible; they are dropped from
/// `perfect` as they leave the window and appended to (or merged into) `result`.
fn save_masked_regions(result: &mut Vec<Interval>, perfect: &mut Vec<Perfect>, start: i64) {
    let start = start as u32;
    let leftmost = match perfect.last() {
        Some(last) => *last,
        None => return,
    };
    if leftmost.start >= start {
        return;
    }
    let interval = Interval {
        start: leftmost.start,
        finish: leftmost.finish,
    };
    // Merge with the previous interval if they overlap or are adjacent.
    match result.last_mut() {
        Some(prev) if interval.start <= prev.finish => {
            prev.finish = prev.finish.max(interval.finish);
        }
        _ => result.push(interval),
    }
    // Forget the perfect intervals that have fallen out of the window.
    let mut keep = perfect.len();
    while keep > 0 && perfect[keep - 1].start < start {
        keep -= 1;
    }
    perfect.truncate(keep);
}

/// Scans the window backwards for perfect intervals and inserts them into
/// `perfect`, keeping it sorted by descending start and ascending finish.
fn find_perfect(perfect: &mut Vec<Perfect>, sliding: &Sliding, total: i64, start: i64) {
    let size = sliding.ring_len();
    let window = sliding.len;
    if size <= window {
        return;
    }
    let mut win_counts = sliding.win_count;
    // Number of repeated pairs of the suffix ending at the window edge.
    let mut win_pairs = sliding.win_pairs;
    // Best score among the intervals that start to the right of the current
    // position, as the fraction `max_r / max_l`. The initial value `0 / 1`
    // behaves like "no interval seen yet": any real interval beats it in the
    // loop below, and the comparison after the loop always succeeds.
    let mut max_pairs = 0u64;
    let mut max_len = 1u64;
    let finish = (size + (WORD_LEN as i64 - 1) + start) as u32;

    for i in (0..=size - window - 1).rev() {
        let word = sliding.ring.get(i as usize) as usize;
        win_pairs += win_counts[word] as i64;
        win_counts[word] += 1;
        let new_win_pairs = win_pairs as u64;
        let new_len = (size - i - 1) as u64;
        // new_win_pairs / new_len < total / 10
        // new_win_pairs / new_len 这个值越大 说明这块区域 越是 low-complexity
        if (new_win_pairs as i64) * 10 <= total * new_len as i64 {
            continue;
        }
        // Find the insertion position, remembering the best score seen so far.
        // The intervals are sorted by descending start, so a binary search finds
        // the first one that starts before `start_i`.
        let start_i = (i + start) as u32;
        let j = perfect.partition_point(|p| p.start >= start_i);
        for cur_perfect in &perfect[..j] {
            // cur_perfect.pairs / cur_perfect.len > max_pairs / max_len
            if cur_perfect.pairs as u64 * max_len > max_pairs * cur_perfect.len as u64 {
                max_pairs = cur_perfect.pairs as u64;
                max_len = cur_perfect.len as u64;
            }
        }
        // Only keep the new interval if it is at least as good as those that
        // start later, which keeps `perfect` free of redundant intervals.
        // new_win_pairs / new_len > max_pairs / max_len
        if new_win_pairs * max_len >= max_pairs * new_len {
            max_pairs = new_win_pairs;
            max_len = new_len;
            perfect.insert(
                j,
                Perfect {
                    start: start_i,
                    finish,
                    pairs: new_win_pairs as u32,
                    len: new_len as u32,
                },
            );
        }
    }
}

/// Masks a single sequence and returns its low-complexity intervals.
///
/// This is the counterpart of `sdust()` in the C API; prefer [`Masker`] when
/// many sequences are processed.
pub fn sdust(seq: &[u8], params: Params) -> Vec<Interval> {
    let mut masker = Masker::new();
    masker.mask(seq, params).to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn intervals(seq: &[u8], params: Params) -> Vec<(u32, u32)> {
        sdust(seq, params)
            .into_iter()
            .map(|i| (i.start, i.finish))
            .collect()
    }

    #[test]
    fn nt4_table() {
        assert_eq!(NT4[b'A' as usize], 0);
        assert_eq!(NT4[b'c' as usize], 1);
        assert_eq!(NT4[b'G' as usize], 2);
        assert_eq!(NT4[b't' as usize], 3);
        for b in 0..=255u8 {
            if !b"ACGTacgt".contains(&b) {
                assert_eq!(NT4[b as usize], 4, "byte {b} should be non-ACGT");
            }
        }
    }

    #[test]
    fn empty_and_short_sequences_have_no_regions() {
        assert!(intervals(b"", Params::default()).is_empty());
        assert!(intervals(b"AC", Params::default()).is_empty());
        assert!(intervals(b"NNNNNNNNNNNNNNNNNNNN", Params::default()).is_empty());
    }

    #[test]
    fn random_sequence_has_no_regions() {
        assert!(intervals(
            b"TGCTATACAGCACTACCCAACAGTTGCTAGCAACTCTGCGCAAGACATCATTCCTTCCGT",
            Params::default()
        )
        .is_empty());
    }

    #[test]
    fn poly_a_is_masked() {
        // Note: a masked interval may start before the low-complexity stretch,
        // and may end past the end of the sequence, because the window covering
        // the stretch is always reported as a whole. This matches the C output.
        let seq = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let got = intervals(seq, Params::default());
        assert_eq!(got, vec![(0, 129)]);
    }

    #[test]
    fn n_breaks_the_sequence_into_independent_pieces() {
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let got = intervals(seq, Params::default());
        assert_eq!(got, vec![(0, 30), (31, 89)]);
    }

    #[test]
    fn masker_state_is_independent_between_sequences() {
        let a = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGT";
        let b = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT";
        let params = Params::default();
        let expected_a = intervals(a, params);
        let expected_b = intervals(b, params);

        let mut masker = Masker::new();
        assert_eq!(masker.mask(a, params).len(), expected_a.len());
        assert_eq!(masker.mask(b, params).len(), expected_b.len());
        // Masking `a` again after `b` must give the same result as the first time.
        let again: Vec<_> = masker
            .mask(a, params)
            .iter()
            .map(|i| (i.start, i.finish))
            .collect();
        assert_eq!(again, expected_a);
    }

    #[test]
    fn lowercase_is_masked_like_uppercase() {
        let mut upper = b"G".repeat(20);
        upper.extend(b"A".repeat(80));
        let mut lower = b"g".repeat(20);
        lower.extend(b"a".repeat(80));
        let expected = intervals(&upper, Params::default());
        assert!(!expected.is_empty());
        assert_eq!(intervals(&lower, Params::default()), expected);
    }

    #[test]
    fn params_validation() {
        assert!(Params::default().validate().is_ok());
        assert!(Params::new(0, 64).validate().is_err());
        assert!(Params::new(20, 2).validate().is_err());
    }

    #[test]
    fn ring_keeps_logical_order() {
        let mut ring = Ring::new();
        for i in 0..100u8 {
            ring.push(i);
        }
        assert_eq!(ring.len(), 100);
        for i in 0..100 {
            assert_eq!(ring.get(i), i as u8);
        }
        for i in 0..50u8 {
            assert_eq!(ring.pop_front(), Some(i));
        }
        assert_eq!(ring.len(), 50);
        assert_eq!(ring.get(0), 50);
        for i in 100..200u8 {
            ring.push(i);
        }
        assert_eq!(ring.get(0), 50);
        assert_eq!(ring.get(149), 199);
    }

    #[test]
    fn tiny_windows_do_not_panic() {
        for window in 3..10 {
            for total in [1, 5, 20] {
                let _ = intervals(
                    b"ACGTACGTAAAAAAAAAAAAAAAAACGTACGTACGTAAAAAAAAAAAAAAAAAACGT",
                    Params::new(total, window),
                );
            }
        }
    }
}
