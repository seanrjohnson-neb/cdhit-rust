//! `cd-hit-lap` — cluster reads that overlap end-to-end.
//!
//! Port of `cd-hit-auxtools/cdhit-lap.cxx`. Sequences are sorted longest-first
//! and each is tested for a perfect (0-mismatch) suffix→prefix or prefix→suffix
//! overlap against existing cluster representatives, on both strands. A `WORD`
//! of length `-m` is hashed (2-bit packed) to index candidate representatives
//! by their head and tail words. Output (`.clstr` + representative FASTA) is
//! bit-for-bit compatible with the reference.
//!
//! Quirks reproduced: the reverse-complement table maps every non-`ACGT` base
//! (e.g. `N`) to `A`; the reported identity is always `100.00%` (`mm` is never
//! set in the C++); residues of the file's last record are not upper-cased (see
//! [`super::bioseq`]).

use std::collections::HashMap;

use super::bioseq::{parse_fastaq, Sequence};

/// Parameters mirroring the `cd-hit-lap` command line.
pub struct LapParams {
    /// `-m`: minimum overlap length (also the hash word length). Default 20.
    pub minlen: i32,
    /// `-p`: minimum overlap as a fraction of the read length. Default 0.
    pub minper: f32,
    /// `-d`: description length (0 = truncate at first whitespace).
    pub deslen: i32,
    /// `-s`: shuffle seed (0 = no shuffle). Non-zero seeds are NOT bit-identical
    /// to the C++ (which depends on libc `rand`); only the default path matches.
    pub seed: u32,
}

/// Result of a `cd-hit-lap` run: representative FASTA/FASTQ, the `.clstr` file,
/// and the log text.
pub struct LapOutput {
    pub rep: Vec<u8>,
    pub clstr: Vec<u8>,
    pub log: String,
}

/// Base code table (`a/A`→0, `c/C`→1, `g/G`→2, `t/T`→3, else 0).
fn base_code(c: u8) -> u64 {
    match c {
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' => 3,
        _ => 0,
    }
}

/// Reverse-complement table: `memset('A')` then a/c/g/t overrides — so any
/// non-`ACGT` base maps to `A`.
fn rev_comp(c: u8) -> u8 {
    match c {
        b'a' => b't',
        b'A' => b'T',
        b'c' => b'g',
        b'C' => b'G',
        b'g' => b'c',
        b'G' => b'C',
        b't' => b'a',
        b'T' => b'A',
        _ => b'A',
    }
}

/// Encode `word[0..w]` as a 2-bit-packed integer (port of `EncodeWord`).
fn encode_word(word: &[u8], w: usize) -> u64 {
    let mut code = 0u64;
    let mut power = 1u64;
    for &c in &word[..w] {
        code += base_code(c) * power;
        power = power.wrapping_mul(4);
    }
    code
}

#[derive(Clone)]
struct ClusterMember {
    seq: usize, // index into the sorted sequence list
    offset: i32,
    offset2: i32,
    revcomp: bool,
}

struct SequenceCluster {
    members: Vec<ClusterMember>,
}

/// Counting sort by descending length (port of `SequenceList::SortByLength`).
fn sort_by_length(seqs: &mut Vec<Sequence>, max_len: usize, min_len: usize) {
    if max_len == min_len {
        return;
    }
    let n = seqs.len();
    let m = max_len - min_len + 1;
    let mut count = vec![0usize; m];
    let mut accum = vec![0usize; m];
    let mut offset = vec![0usize; m];
    for s in seqs.iter() {
        count[max_len - s.len()] += 1;
    }
    for i in 1..m {
        accum[i] = accum[i - 1] + count[i - 1];
    }
    let mut sorted: Vec<Option<Sequence>> = (0..n).map(|_| None).collect();
    for s in seqs.drain(..) {
        let len = max_len - s.len();
        let id = accum[len] + offset[len];
        offset[len] += 1;
        sorted[id] = Some(s);
    }
    *seqs = sorted.into_iter().map(|s| s.unwrap()).collect();
}

/// Trim a description to its identifier token (port of
/// `SequenceCluster::GetDescription`, equivalent to `Sequence::get_description`).
fn get_description(seq: &Sequence, deslen: i32) -> Vec<u8> {
    seq.get_description(if deslen < 0 { 0 } else { deslen as usize })
}

fn write_cluster(
    out: &mut Vec<u8>,
    seqs: &[Sequence],
    cluster: &SequenceCluster,
    id: i32,
    deslen: i32,
) {
    let members = &cluster.members;
    let rep = &seqs[members[0].seq];
    let des = get_description(rep, deslen);
    let len2 = rep.len() as i32;

    out.extend_from_slice(format!(">Cluster {}\n", id).as_bytes());
    out.extend_from_slice(b"0\t");
    out.extend_from_slice(len2.to_string().as_bytes());
    out.extend_from_slice(b"nt, >");
    out.extend_from_slice(&des);
    out.extend_from_slice(b"... *\n");

    for i in 1..members.len() {
        let m = &members[i];
        let s = &seqs[m.seq];
        let des_i = get_description(s, deslen);
        let len = s.len() as i32;
        let mut start_q = m.offset + 1;
        let start_r = m.offset2 + 1;
        let mut end_q = len;
        let mut end_r = len2;
        if (end_q - start_q) < (end_r - start_r) {
            end_r = start_r + end_q - start_q;
        } else {
            end_q = start_q + end_r - start_r;
        }
        if m.revcomp {
            start_q = len - start_q + 1;
            end_q = len - end_q + 1;
        }
        let pct = 100.0f32 * (len as f32) / (len as f32); // mm always 0 -> 100.00
        out.extend_from_slice(
            format!(
                "{}\t{}nt, >{}... at {}:{}:{}:{}/{}/{:.2}%\n",
                i,
                len,
                String::from_utf8_lossy(&des_i),
                start_q,
                end_q,
                start_r,
                end_r,
                if m.revcomp { '-' } else { '+' },
                pct
            )
            .as_bytes(),
        );
    }
}

/// Port of `ClusterOverlap`.
fn cluster_overlap(
    seqs: &[Sequence],
    min: usize,
    minper: f32,
    log: &mut String,
) -> Vec<SequenceCluster> {
    let word = min;
    let mut clusters: Vec<SequenceCluster> = Vec::new();
    let mut head_hashes: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut tail_hashes: HashMap<u64, Vec<usize>> = HashMap::new();

    for i in 0..seqs.len() {
        let seq = &seqs[i];
        let len = seq.len();
        let mut min2 = (len as f32 * minper) as usize;
        if min2 < min {
            min2 = min;
        }
        let mut clustered = false;

        // reverse-complement of the read (N -> A via the table).
        let mut reverse = vec![0u8; len];
        for j in 0..len {
            reverse[j] = rev_comp(seq.seq[len - j - 1]);
        }

        // Phase 1: suffix words of this read against representative tail words.
        if len >= word + 1 && min2 >= word {
            let mut j = len as i32 - word as i32 - 1;
            let low = min2 as i32 - word as i32;
            while j >= low {
                for i2 in 0..2 {
                    let ss: &[u8] = if i2 == 1 { &reverse } else { &seq.seq };
                    let hash = encode_word(&ss[j as usize..], word);
                    if let Some(bucket) = tail_hashes.get(&hash) {
                        for &ci in bucket {
                            let rep = &seqs[clusters[ci].members[0].seq].seq;
                            if rep.len() < (j as usize + word) {
                                continue;
                            }
                            let off = rep.len() - j as usize - word;
                            let big_m = j as usize + word;
                            if ss[..big_m] == rep[off..off + big_m] {
                                clusters[ci].members.push(ClusterMember {
                                    seq: i,
                                    offset: 0,
                                    offset2: (rep.len() - j as usize - word) as i32,
                                    revcomp: i2 == 1,
                                });
                                clustered = true;
                                break;
                            }
                        }
                    }
                    if clustered {
                        break;
                    }
                }
                if clustered {
                    break;
                }
                j -= 1;
            }
        }

        // Phase 2: prefix words of this read against representative head words.
        if !clustered && len >= min2 {
            let mut j = 0usize;
            while j <= len - min2 {
                for i2 in 0..2 {
                    let ss: &[u8] = if i2 == 1 { &reverse } else { &seq.seq };
                    let hash = encode_word(&ss[j..], word);
                    if let Some(bucket) = head_hashes.get(&hash) {
                        for &ci in bucket {
                            let rep = &seqs[clusters[ci].members[0].seq].seq;
                            if rep.len() < (len - j) {
                                continue;
                            }
                            let big_m = len - j;
                            if ss[j..j + big_m] == rep[..big_m] {
                                clusters[ci].members.push(ClusterMember {
                                    seq: i,
                                    offset: j as i32,
                                    offset2: 0,
                                    revcomp: i2 == 1,
                                });
                                clustered = true;
                                break;
                            }
                        }
                    }
                    if clustered {
                        break;
                    }
                }
                if clustered {
                    break;
                }
                j += 1;
            }
        }

        if !clustered {
            let ci = clusters.len();
            clusters.push(SequenceCluster {
                members: vec![ClusterMember {
                    seq: i,
                    offset: 0,
                    offset2: 0,
                    revcomp: false,
                }],
            });
            // Reads shorter than the word length cannot be indexed (the C++
            // reads out of bounds here); such a read just becomes a singleton.
            if len >= word {
                let head = encode_word(&seq.seq, word);
                let tail = encode_word(&seq.seq[len - word..], word);
                head_hashes.entry(head).or_default().push(ci);
                tail_hashes.entry(tail).or_default().push(ci);
            }
        }

        if (i + 1) % 100000 == 0 {
            log.push_str(&format!(
                "Clustered {:9} sequences with {:9} clusters ...\n",
                i + 1,
                clusters.len()
            ));
        }
    }
    clusters
}

/// Run `cd-hit-lap` over an in-memory FASTA/FASTQ buffer.
pub fn cd_hit_lap(input: &[u8], params: &LapParams) -> LapOutput {
    let mut log = String::new();
    let mut seqs = parse_fastaq(input, false);

    let count = seqs.len();
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);

    log.push_str(&format!("Total number of sequences: {}\n", count));
    log.push_str(&format!("Longest: {}\n", max_len));
    log.push_str(&format!("Shortest: {}\n", min_len));

    if params.seed != 0 {
        shuffle(&mut seqs, params.seed);
    }
    if max_len != min_len {
        sort_by_length(&mut seqs, max_len, min_len);
        log.push_str("Sorted by length ...\n");
    }

    log.push_str("Start clustering duplicated sequences ...\n");
    let clusters = cluster_overlap(&seqs, params.minlen as usize, params.minper, &mut log);
    log.push_str(&format!("Number of clusters found: {}\n", clusters.len()));

    // WriteClusters: rep FASTA + .clstr, renumbering non-empty clusters.
    let mut rep = Vec::new();
    let mut clstr = Vec::new();
    let mut k1 = 0i32;
    for cluster in &clusters {
        if cluster.members.is_empty() {
            continue;
        }
        seqs[cluster.members[0].seq].print(&mut rep);
        write_cluster(&mut clstr, &seqs, cluster, k1, params.deslen);
        k1 += 1;
    }
    log.push_str("Done!\n");

    LapOutput { rep, clstr, log }
}

/// Port of `SequenceList::Shuffle`. NOTE: uses a local PRNG, not libc `rand`,
/// so output under `-s` is not bit-identical to the C++.
fn shuffle(seqs: &mut [Sequence], seed: u32) {
    let n = seqs.len();
    if n <= 1 {
        return;
    }
    let mut state = seed as u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    for i in (0..n).rev() {
        let rnd = (next() as usize) % (i + 1);
        seqs.swap(i, rnd);
    }
}
