//! Shared `.clstr` model + parser, plus ports of the CD-HIT Perl
//! post-processing scripts as subcommands.
//!
//! CD-HIT's cluster file looks like:
//! ```text
//! >Cluster 0
//! 0	150nt, >seqA... *
//! 1	150nt, >seqB... at 1:150:1:150/+/98.00%
//! ```
//! Each subcommand mirrors one `clstr_*.pl` script and reproduces its exact
//! output (bit-for-bit), reusing this parser. The parser keeps each member's
//! raw line so scripts that re-emit clusters verbatim (e.g. `clstr_sort_by`)
//! stay lossless.

pub mod ops;

/// A single member line within a cluster.
#[derive(Clone)]
pub struct Member {
    /// The full original line, without its trailing newline.
    pub raw: Vec<u8>,
    /// Length (the `<n>aa`/`<n>nt` number).
    pub len: i64,
    /// Identifier (text between `>` and `...`).
    pub id: Vec<u8>,
    /// Whether this is the representative (line ends in `*`).
    pub is_rep: bool,
    /// Identity string as it appears, e.g. `100%` or `98.00%` (`100%` for the
    /// representative). Empty if the line had no parseable identity.
    pub iden: Vec<u8>,
}

/// One cluster: its header line and member lines.
#[derive(Clone)]
pub struct Cluster {
    /// Cluster number parsed from `>Cluster N` (−1 if unparented).
    pub num: i64,
    /// The full `>Cluster ...` header line, without trailing newline.
    pub header_raw: Vec<u8>,
    pub members: Vec<Member>,
}

impl Cluster {
    /// Maximum member length (the representative's length in practice).
    pub fn max_len(&self) -> i64 {
        self.members.iter().map(|m| m.len).max().unwrap_or(0)
    }
    /// Number of members.
    pub fn size(&self) -> usize {
        self.members.len()
    }
}

fn strip_newline(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

/// Parse a member line, mirroring the regexes in the Perl scripts:
/// `\d+\t(\d+)[a-z]{2}, >(.+)\.\.\. \*` (representative) and
/// `\d+\t(\d+)[a-z]{2}, >(.+)\.\.\.` with trailing `(\d+%|\d+\.\d+%)` identity.
fn parse_member(raw: &[u8]) -> Member {
    let s = raw;
    let mut len = 0i64;
    let mut id = Vec::new();

    // find the tab, then the length digits after it, then "aa"/"nt", ", >".
    if let Some(tab) = s.iter().position(|&b| b == b'\t') {
        let after = &s[tab + 1..];
        // length digits
        let mut i = 0;
        while i < after.len() && after[i].is_ascii_digit() {
            i += 1;
        }
        if i > 0 {
            len = std::str::from_utf8(&after[..i]).unwrap().parse().unwrap_or(0);
        }
        // find ">" then id up to the last "..."
        if let Some(gt) = after.iter().position(|&b| b == b'>') {
            let rest = &after[gt + 1..];
            if let Some(dots) = find_last(rest, b"...") {
                id = rest[..dots].to_vec();
            }
        }
    }
    // representative if the line ends with "*" (after "... ").
    let (is_rep, iden) = if strip_newline(s).last() == Some(&b'*') {
        (true, b"100%".to_vec())
    } else {
        (false, trailing_identity(strip_newline(s)))
    };
    Member {
        raw: strip_newline(s).to_vec(),
        len,
        id,
        is_rep,
        iden,
    }
}

/// Last occurrence of `needle` in `hay`.
fn find_last(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    (0..=hay.len() - needle.len())
        .rev()
        .find(|&i| &hay[i..i + needle.len()] == needle)
}

/// Extract a trailing `\d+%` or `\d+\.\d+%` from the end of the line.
fn trailing_identity(s: &[u8]) -> Vec<u8> {
    if s.last() != Some(&b'%') {
        return Vec::new();
    }
    let mut i = s.len() - 1; // at '%'
    let mut j = i; // walk back over digits and one dot
    let mut seen_digit = false;
    let mut seen_dot = false;
    while j > 0 {
        let c = s[j - 1];
        if c.is_ascii_digit() {
            seen_digit = true;
            j -= 1;
        } else if c == b'.' && !seen_dot && seen_digit {
            seen_dot = true;
            j -= 1;
        } else {
            break;
        }
    }
    if seen_digit {
        i += 1; // include '%'
        s[j..i].to_vec()
    } else {
        Vec::new()
    }
}

/// Parse cluster number from a `>Cluster N ...` header (−1 if absent).
fn parse_cluster_num(header: &[u8]) -> i64 {
    // Skip ">Cluster" and whitespace, then read digits.
    let prefix = b">Cluster";
    if header.len() < prefix.len() || &header[..prefix.len()] != prefix {
        return -1;
    }
    let mut i = prefix.len();
    while i < header.len() && (header[i] == b' ' || header[i] == b'\t') {
        i += 1;
    }
    let start = i;
    while i < header.len() && header[i].is_ascii_digit() {
        i += 1;
    }
    if i > start {
        std::str::from_utf8(&header[start..i])
            .unwrap()
            .parse()
            .unwrap_or(-1)
    } else {
        -1
    }
}

/// Parse a `.clstr` buffer into clusters.
pub fn parse_clstr(input: &[u8]) -> Vec<Cluster> {
    let mut clusters: Vec<Cluster> = Vec::new();
    let mut cur: Option<Cluster> = None;
    for line in input.split_inclusive(|&b| b == b'\n') {
        if line.first() == Some(&b'>') {
            if let Some(c) = cur.take() {
                clusters.push(c);
            }
            let header = strip_newline(line).to_vec();
            let num = parse_cluster_num(&header);
            cur = Some(Cluster {
                num,
                header_raw: header,
                members: Vec::new(),
            });
        } else if !strip_newline(line).is_empty() {
            if let Some(c) = cur.as_mut() {
                c.members.push(parse_member(line));
            }
        }
    }
    if let Some(c) = cur.take() {
        clusters.push(c);
    }
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_clstr() {
        let c = b">Cluster 0\n0\t150nt, >seqA... *\n1\t150nt, >seqB... at 1:150:1:150/+/98.00%\n>Cluster 1\n0\t99aa, >seqC... *\n";
        let cs = parse_clstr(c);
        assert_eq!(cs.len(), 2);
        assert_eq!(cs[0].num, 0);
        assert_eq!(cs[0].members.len(), 2);
        assert_eq!(cs[0].members[0].id, b"seqA");
        assert!(cs[0].members[0].is_rep);
        assert_eq!(cs[0].members[0].len, 150);
        assert_eq!(cs[0].members[0].iden, b"100%");
        assert!(!cs[0].members[1].is_rep);
        assert_eq!(cs[0].members[1].id, b"seqB");
        assert_eq!(cs[0].members[1].iden, b"98.00%");
        assert_eq!(cs[1].members[0].id, b"seqC");
        assert_eq!(cs[1].members[0].len, 99);
    }
}
