//! `cd-hit-dup` — detect and remove duplicate / near-duplicate reads.
//!
//! Port of `cd-hit-auxtools/cdhit-dup.cxx`. Supports single-end and paired-end
//! (`-i2`/`-o2`) FASTA/FASTQ, exact and up-to-`-e`-mismatch duplicate
//! detection, match-length control (`-m`) and prefix analysis (`-u`).
//!
//! Bucketing must match the C++ exactly, so [`murmur_hash2`] is reproduced
//! byte-for-byte (the C++ maps key on the MurmurHash value, and near-duplicate
//! discovery depends on the shared-subword hash tables). The compare/mismatch
//! primitives mirror `mintlib/minString.cxx`.
//!
//! NOTE: chimera filtering (`-f`/`-s`) is handled in [`super::dup_chimera`] and
//! wired in here; the re-read-by-index behaviour of the C++ (used to restore
//! full-length reads for output after `-u`/PE modification) is reproduced
//! faithfully, including its dependence on whether length-sorting occurred.

use std::collections::HashMap;

use super::bioseq::{parse_fastaq, Sequence};

const HASH_SEED: u32 = 0xda0;

/// MurmurHash2 (Austin Appleby), matching `mintlib/minMap.cxx` including its
/// little-endian 4-byte reads (the C++ is explicitly endian-dependent; the
/// reference was built on x86, so little-endian is correct).
pub fn murmur_hash2(key: &[u8], seed: u32) -> u32 {
    let m: u32 = 0x5bd1_e995;
    let r = 24;
    let mut len = key.len();
    let mut h: u32 = seed ^ (len as u32);
    let mut i = 0usize;
    while len >= 4 {
        let mut k = u32::from_le_bytes([key[i], key[i + 1], key[i + 2], key[i + 3]]);
        k = k.wrapping_mul(m);
        k ^= k >> r;
        k = k.wrapping_mul(m);
        h = h.wrapping_mul(m);
        h ^= k;
        i += 4;
        len -= 4;
    }
    if len == 3 {
        h ^= (key[i + 2] as u32) << 16;
    }
    if len >= 2 {
        h ^= (key[i + 1] as u32) << 8;
    }
    if len >= 1 {
        h ^= key[i] as u32;
        h = h.wrapping_mul(m);
    }
    h ^= h >> 13;
    h = h.wrapping_mul(m);
    h ^= h >> 15;
    h
}

/// `MakeHash(s, length)`: hash the first `length` bytes (`length<=0` = full).
fn make_hash(s: &[u8], length: i32) -> u32 {
    let mut length = length;
    let n = s.len() as i32;
    if length <= 0 {
        length = n;
    }
    if length > n {
        length = n;
    }
    murmur_hash2(&s[..length as usize], HASH_SEED)
}

/// `MakeHash(s, offset, length)`: hash `length` bytes from `offset`.
pub(crate) fn make_hash_off(s: &[u8], offset: i32, length: i32) -> u32 {
    let n = s.len() as i32;
    let mut offset = offset;
    let mut length = length;
    if offset < 0 {
        offset = 0;
    }
    if length <= 0 {
        length = n - offset;
    }
    if length + offset > n {
        length = n - offset;
    }
    let o = offset as usize;
    murmur_hash2(&s[o..o + length as usize], HASH_SEED)
}

/// `Compare`: lexicographic, length-aware (equivalent to the C++ `Compare`).
fn compare(left: &[u8], right: &[u8]) -> i32 {
    let m = left.len();
    let n = right.len();
    let mut i = 0;
    while i < m && i < n && left[i] == right[i] {
        i += 1;
    }
    if m == n && i == n {
        return 0;
    }
    if i == m {
        return -1;
    }
    if i == n {
        return 1;
    }
    if left[i] < right[i] {
        -1
    } else {
        1
    }
}

/// `ComparePrefix(left, right, maxmm)`: number of leading positions matched
/// within `maxmm` mismatches.
pub(crate) fn compare_prefix(left: &[u8], right: &[u8], maxmm: i32) -> usize {
    let k = left.len().min(right.len());
    let mut i = 0usize;
    let mut mm = 0i32;
    while i < k {
        mm += (left[i] != right[i]) as i32;
        if mm > maxmm {
            break;
        }
        i += 1;
    }
    i
}

/// 6-arg `ComparePrefix`: returns the number of aligned positions matched from
/// (`lstart`,`rstart`) within `maxmm` mismatches, capped by `max` (0 = no cap).
pub(crate) fn compare_prefix2(
    left: &[u8],
    lstart: i32,
    right: &[u8],
    rstart: i32,
    max: i32,
    maxmm: i32,
) -> i32 {
    let mut m = left.len() as i32;
    let mut n = right.len() as i32;
    if max != 0 {
        if lstart + max < m {
            m = lstart + max;
        }
        if rstart + max < n {
            n = rstart + max;
        }
    }
    let mut i = lstart;
    let mut j = rstart;
    let mut mm = 0i32;
    while i < m && j < n {
        mm += (left[i as usize] != right[j as usize]) as i32;
        if mm > maxmm {
            break;
        }
        i += 1;
        j += 1;
    }
    i - lstart
}

/// `CountMismatch(left, right, from, to)`: mismatches over `[from,to)`
/// (`to<0` = to end), clamped to both lengths.
pub(crate) fn count_mismatch(left: &[u8], right: &[u8], from: i32, to: i32) -> i32 {
    let m = left.len() as i32;
    let n = right.len() as i32;
    let mut to = if to < 0 { m } else { to };
    if to > m {
        to = m;
    }
    if to > n {
        to = n;
    }
    let mut mm = 0i32;
    let mut i = from;
    while i < to {
        mm += (left[i as usize] != right[i as usize]) as i32;
        i += 1;
    }
    mm
}

/// `CountMismatch2`: like [`count_mismatch`] but stops once `maxmm` is exceeded
/// (may return `maxmm+1`).
pub(crate) fn count_mismatch2(left: &[u8], right: &[u8], from: i32, to: i32, maxmm: i32) -> i32 {
    let m = left.len() as i32;
    let n = right.len() as i32;
    let mut to = if to < 0 { m } else { to };
    if to > m {
        to = m;
    }
    if to > n {
        to = n;
    }
    let mut mm = 0i32;
    let mut i = from;
    while i < to && mm <= maxmm {
        mm += (left[i as usize] != right[i as usize]) as i32;
        i += 1;
    }
    mm
}

/// `CountMismatch3`: mismatches from (`lstart`,`rstart`) for up to `max`
/// aligned positions, stopping once `maxmm` is exceeded.
pub(crate) fn count_mismatch3(
    left: &[u8],
    lstart: i32,
    right: &[u8],
    rstart: i32,
    max: i32,
    maxmm: i32,
) -> i32 {
    let mut m = left.len() as i32;
    let mut n = right.len() as i32;
    if max != 0 {
        if lstart + max < m {
            m = lstart + max;
        }
        if rstart + max < n {
            n = rstart + max;
        }
    }
    let mut i = lstart;
    let mut j = rstart;
    let mut mm = 0i32;
    while i < m && j < n {
        mm += (left[i as usize] != right[j as usize]) as i32;
        if mm > maxmm {
            break;
        }
        i += 1;
        j += 1;
    }
    mm
}

/// `CompareSuffix`: number of positions matched walking backwards from
/// (`lstart`,`rstart`) within `maxmm` mismatches (`max` bounds the window).
pub(crate) fn compare_suffix(
    left: &[u8],
    lstart: i32,
    right: &[u8],
    rstart: i32,
    max: i32,
    maxmm: i32,
) -> i32 {
    let mut m = 0i32;
    let mut n = 0i32;
    if max != 0 {
        if lstart + 1 - max >= 0 {
            m = lstart + 1 - max;
        }
        if rstart + 1 - max >= 0 {
            n = rstart + 1 - max;
        }
    }
    let mut i = lstart;
    let mut j = rstart;
    let mut mm = 0i32;
    while i >= m && j >= n {
        mm += (left[i as usize] != right[j as usize]) as i32;
        if mm > maxmm {
            break;
        }
        i -= 1;
        j -= 1;
    }
    lstart - i
}

fn hashing_depth(len: i32, min: i32) -> i32 {
    // (int)sqrt( (len - min) / 10 ) with integer division; negative -> 0.
    let v = (len - min) / 10;
    if v < 0 {
        0
    } else {
        (v as f64).sqrt() as i32
    }
}

fn hashing_length(dep: i32, min: i32) -> i32 {
    min + 10 * dep * dep
}

/// A cluster of duplicate reads (indices into the working sequence list).
pub struct DupCluster {
    pub members: Vec<usize>,
    pub id: i32,
    pub abundance: i32,
    pub chi_head: i32,
    pub chi_tail: i32,
}

impl DupCluster {
    fn new(rep: usize) -> Self {
        DupCluster {
            members: vec![rep],
            id: 0,
            abundance: 0,
            chi_head: 0,
            chi_tail: 0,
        }
    }
}

type HashTable = HashMap<u32, Vec<i32>>;

fn update_hash_tables(
    middles: &mut [Vec<HashTable>],
    seq: &[u8],
    id: i32,
    shared: i32,
    min: i32,
    primer: i32,
) {
    let size = seq.len() as i32;
    let mut j = primer;
    while (j + shared) <= size {
        let dep = if j != 0 { 0 } else { hashing_depth(size - j, min) };
        for k in 0..=dep {
            let hash = make_hash_off(seq, j, hashing_length(k, min));
            middles[j as usize][k as usize]
                .entry(hash)
                .or_default()
                .push(id);
        }
        j += shared;
    }
}

/// Port of `ClusterDuplicate`. Greedy: each read joins the first existing
/// cluster whose representative it exactly (or within `maxmm`) matches,
/// otherwise starts a new cluster.
pub fn cluster_duplicate(
    seqlist: &[Sequence],
    mlen: bool,
    errors: i32,
    errors2: f32,
    log: &mut String,
) -> Vec<DupCluster> {
    let n = seqlist.len();
    let max = seqlist.iter().map(|s| s.len()).max().unwrap_or(0) as i32;
    let minlen = seqlist.iter().map(|s| s.len()).min().unwrap_or(0) as i32;

    let mut clusters: Vec<DupCluster> = Vec::new();
    let mut wholes: HashTable = HashMap::new();

    // primer = longest exact common prefix across all reads (>=15, else 0).
    let mut primer = seqlist[0].seq.len() as i32;
    for i in 1..n {
        primer = compare_prefix2(&seqlist[0].seq, 0, &seqlist[i].seq, 0, primer, 0);
        if primer < 15 {
            primer = 0;
            break;
        }
    }
    let mut maxmm0 = errors;
    if errors == 0 {
        maxmm0 = (max as f32 * errors2) as i32;
    }
    log.push_str(&format!("primer = {}\n", primer));
    let mut shared = (minlen - primer) / (maxmm0 + 1);
    if shared < 30 {
        shared = 30;
    }
    let min = shared;

    // middles[j][k]: shared-subword hash table at position j, depth k.
    let mut middles: Vec<Vec<HashTable>> = (0..max.max(0) as usize).map(|_| Vec::new()).collect();
    let mut i = 0i32;
    while (i + shared) <= max {
        let depth = hashing_depth(max - i, min) + 1;
        middles[i as usize] = (0..depth).map(|_| HashMap::new()).collect();
        i += 1;
    }

    for i in 0..n {
        let seq = &seqlist[i].seq;
        let qlen = seq.len() as i32;
        let mut clustered = false;
        let mut maxmm = errors;
        if errors == 0 {
            maxmm = (qlen as f32 * errors2) as i32;
        }
        if i > 0 && i % 10000 == 0 {
            log.push_str(&format!(
                "Clustered {:9} sequences with {:9} clusters ...\n",
                i,
                clusters.len()
            ));
        }

        // whole-sequence table: exact duplicate, then within-maxmm prefix.
        let hash = make_hash(seq, 0);
        if let Some(hits) = wholes.get(&hash) {
            for &h in hits {
                let rep = &seqlist[clusters[h as usize].members[0]].seq;
                if qlen as usize == rep.len() && compare(seq, rep) == 0 {
                    clusters[h as usize].members.push(i);
                    clustered = true;
                    break;
                }
            }
            if !clustered {
                for &h in hits {
                    let rep = &seqlist[clusters[h as usize].members[0]].seq;
                    if mlen && rep.len() as i32 != qlen {
                        continue;
                    }
                    if compare_prefix(seq, rep, maxmm) as i32 == qlen {
                        clusters[h as usize].members.push(i);
                        clustered = true;
                        break;
                    }
                }
            }
        }
        if clustered {
            continue;
        }

        // shared-subword table at the primer offset.
        let dep = hashing_depth(qlen - primer, min);
        let hashlen = hashing_length(dep, min);
        let hash = make_hash_off(seq, primer, hashlen);
        if let Some(bucket) = middles
            .get(primer as usize)
            .and_then(|v| v.get(dep as usize))
        {
            if let Some(hits) = bucket.get(&hash) {
                for &h in hits {
                    let rep = &seqlist[clusters[h as usize].members[0]].seq;
                    if mlen && rep.len() as i32 != qlen {
                        continue;
                    }
                    if compare_prefix(seq, rep, maxmm) as i32 == qlen {
                        clusters[h as usize].members.push(i);
                        clustered = true;
                        break;
                    }
                }
            }
        }

        // sliding shared-subword probes (only when maxmm == 0).
        let mut start = primer;
        while maxmm == 0
            && !clustered
            && (start + 2 * shared) <= qlen
            && start <= primer + (maxmm + 1) * shared
        {
            let hash = make_hash_off(seq, start, shared);
            let found = middles
                .get(start as usize)
                .and_then(|v| v.first())
                .and_then(|b| b.get(&hash))
                .cloned();
            start += shared;
            if let Some(hits) = found {
                for h in hits {
                    let rep = &seqlist[clusters[h as usize].members[0]].seq;
                    if mlen && rep.len() as i32 != qlen {
                        continue;
                    }
                    if compare_prefix(seq, rep, maxmm) as i32 == qlen {
                        clusters[h as usize].members.push(i);
                        clustered = true;
                        break;
                    }
                }
            }
        }

        if !clustered {
            let hash = make_hash(seq, 0);
            update_hash_tables(&mut middles, seq, clusters.len() as i32, shared, min, primer);
            wholes.entry(hash).or_default().push(clusters.len() as i32);
            clusters.push(DupCluster::new(i));
        }
    }
    let _ = maxmm0;
    clusters
}

/// Counting sort by descending length (port of `SequenceList::SortByLength`).
pub(crate) fn sort_by_length(seqs: &mut Vec<Sequence>) {
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
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

fn write_cluster(
    out: &mut Vec<u8>,
    seqlist: &[Sequence],
    cluster: &DupCluster,
    id: i32,
    deslen: i32,
    cdes: Option<&str>,
) {
    let deslen = deslen.max(0) as usize;
    let rep = &seqlist[cluster.members[0]];
    let des = rep.get_description(deslen);
    out.extend_from_slice(format!(">Cluster {}{}\n", id, cdes.unwrap_or("")).as_bytes());
    out.extend_from_slice(b"0\t");
    out.extend_from_slice(rep.len().to_string().as_bytes());
    out.extend_from_slice(b"nt, >");
    out.extend_from_slice(&des);
    out.extend_from_slice(b"... *\n");
    for idx in 1..cluster.members.len() {
        let m = &seqlist[cluster.members[idx]];
        let des_m = m.get_description(deslen);
        let len = m.len() as i32;
        let mm = count_mismatch(&rep.seq, &m.seq, 0, -1);
        let pct = ((100 * (len - mm)) as f32 / len as f32) as f64;
        out.extend_from_slice(format!("{}\t{}nt, >", idx, len).as_bytes());
        out.extend_from_slice(&des_m);
        out.extend_from_slice(format!("... at 1:{}:1:{}/+/{:.2}%\n", len, len, pct).as_bytes());
    }
}

/// Emit the `.clstr`/`2.clstr` files (and, when `with_reps`, representative
/// records) — port of `WriteClusters` / `WriteClusters_clstronly`.
fn write_clusters(
    seqlist: &[Sequence],
    clusters: &mut [DupCluster],
    deslen: i32,
    with_reps: bool,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let mut reps = Vec::new();
    let mut clstr = Vec::new();
    let mut clstr2 = Vec::new();
    let mut k1 = 0i32;
    let mut k2 = 0i32;
    for i in 0..clusters.len() {
        if clusters[i].members.is_empty() {
            continue;
        }
        let head = clusters[i].chi_head;
        let tail = clusters[i].chi_tail;
        if head == tail {
            if with_reps {
                seqlist[clusters[i].members[0]].print(&mut reps);
            }
            clusters[i].id = k1;
            write_cluster(&mut clstr, seqlist, &clusters[i], k1, deslen, None);
            k1 += 1;
        } else {
            let head_id = clusters[head as usize].id;
            let tail_id = clusters[tail as usize].id;
            let cdes = format!(
                " chimeric_parent1={},chimeric_parent2={}",
                head_id, tail_id
            );
            write_cluster(&mut clstr2, seqlist, &clusters[i], k2, deslen, Some(&cdes));
            k2 += 1;
        }
    }
    (reps, clstr, clstr2)
}

/// Representative-only output (port of `WriteClusters_seqonly`).
fn write_reps_only(seqlist: &[Sequence], clusters: &[DupCluster]) -> Vec<u8> {
    let mut reps = Vec::new();
    for c in clusters {
        if c.members.is_empty() {
            continue;
        }
        if c.chi_head == c.chi_tail {
            seqlist[c.members[0]].print(&mut reps);
        }
    }
    reps
}

/// Parameters mirroring the `cd-hit-dup` command line.
pub struct DupParams<'a> {
    pub input_name: &'a str,
    /// Second (R2) input: (bytes, name).
    pub input2: Option<(&'a [u8], &'a str)>,
    pub match_length: bool, // -m
    pub abundance: i32,     // -a (-1 = auto)
    pub deslen: i32,        // -d
    pub uselen: i32,        // -u
    pub errors: i32,        // -e (int form)
    pub errors2: f32,       // -e (float form)
    pub nochimeric: bool,   // -f / -s
    pub shared: i32,        // -s
    pub abratio: f32,       // -b
    pub percent: f32,       // -p
}

/// Result of a `cd-hit-dup` run.
pub struct DupOutput {
    pub log: String,
    pub reps_r1: Vec<u8>,
    pub reps_r2: Option<Vec<u8>>,
    pub clstr: Vec<u8>,
    pub clstr2: Vec<u8>,
}

fn max_len(seqs: &[Sequence]) -> i32 {
    seqs.iter().map(|s| s.len()).max().unwrap_or(0) as i32
}
fn min_len(seqs: &[Sequence]) -> i32 {
    seqs.iter().map(|s| s.len()).min().unwrap_or(0) as i32
}

/// Run `cd-hit-dup`. Returns an error string (matching the C++ stderr message)
/// for the fatal conditions the C++ reports with a non-zero exit.
pub fn cd_hit_dup(input: &[u8], params: &DupParams) -> Result<DupOutput, String> {
    let mut log = String::new();
    let mut seqlist = parse_fastaq(input, false);
    let mut seqlist_modified = false;

    let seqlist2_src = params.input2.map(|(b, _)| parse_fastaq(b, false));

    if params.nochimeric && params.input2.is_some() {
        return Err("ERROR: Chimeric filtering can only be done with single (-i) input file!".into());
    }

    log.push_str(&format!("From input: {}\n", params.input_name));
    log.push_str(&format!("Total number of sequences: {}\n", seqlist.len()));
    log.push_str(&format!("Longest: {}\n", max_len(&seqlist)));
    log.push_str(&format!("Shortest: {}\n", min_len(&seqlist)));

    let mut uselen = params.uselen;
    if let Some(seqlist2) = &seqlist2_src {
        let (_, name2) = params.input2.unwrap();
        if seqlist2.len() != seqlist.len() {
            return Err("Error: the pair end files contain different number of reads!".into());
        }
        log.push_str(&format!("\nFrom input: {}\n", name2));
        log.push_str(&format!("Total number of sequences: {}\n", seqlist2.len()));
        log.push_str(&format!("Longest: {}\n", max_len(seqlist2)));
        log.push_str(&format!("Shortest: {}\n", min_len(seqlist2)));

        let max = max_len(&seqlist).max(max_len(seqlist2));
        if uselen <= 0 {
            uselen = max;
        }
        if uselen > max {
            uselen = max;
        }
        let ul = uselen as usize;
        for k in 0..seqlist.len() {
            let has_qs = !seqlist[k].qs.is_empty();
            let s1 = seqlist[k].seq.clone();
            let s2 = &seqlist2[k].seq;
            let qs1 = seqlist[k].qs.clone();
            let qs2 = &seqlist2[k].qs;
            let mut ns = Vec::with_capacity(2 * ul);
            for j in 0..ul {
                ns.push(if j < s1.len() { s1[j] } else { b'N' });
            }
            for j in 0..ul {
                ns.push(if j < s2.len() { s2[j] } else { b'N' });
            }
            seqlist[k].seq = ns;
            if has_qs {
                let mut nq = Vec::with_capacity(2 * ul);
                for j in 0..ul {
                    nq.push(if j < qs1.len() { qs1[j] } else { 0 });
                }
                for j in 0..ul {
                    nq.push(if j < qs2.len() { qs2[j] } else { 0 });
                }
                seqlist[k].qs = nq;
            }
        }
        seqlist_modified = true;
    } else if uselen > 0 {
        let ul = uselen as usize;
        for k in 0..seqlist.len() {
            let n1 = seqlist[k].seq.len();
            if n1 < ul {
                continue;
            }
            seqlist_modified = true;
            seqlist[k].seq.truncate(ul);
            if !seqlist[k].qs.is_empty() {
                seqlist[k].qs.truncate(ul);
            }
        }
    }

    if seqlist.is_empty() {
        return Err("Abort: no sequence available for the analysis!".into());
    }

    if max_len(&seqlist) != min_len(&seqlist) {
        sort_by_length(&mut seqlist);
        log.push_str("Sorted by length ...\n");
    }

    log.push_str("Start clustering duplicated sequences ...\n");
    let mut clusters = cluster_duplicate(
        &seqlist,
        params.match_length,
        params.errors,
        params.errors2,
        &mut log,
    );
    log.push_str(&format!("Number of reads: {}\n", seqlist.len()));
    log.push_str(&format!("Number of clusters found: {}\n", clusters.len()));

    for c in clusters.iter_mut() {
        c.abundance = c.members.len() as i32;
    }

    let mut abundance = params.abundance;
    if params.nochimeric {
        if abundance < 0 {
            abundance = 2;
        }
        if seqlist.len() == clusters.len() {
            for c in clusters.iter_mut() {
                let des = &seqlist[c.members[0]].des;
                let mut ab = 1i32;
                if let Some(pos) = find_subslice(des, b"_abundance_") {
                    ab = parse_leading_int(&des[pos + 11..]);
                }
                c.abundance = ab;
            }
        }
        super::dup_chimera::sort_by_abundance(&mut clusters);
        let chistat = super::dup_chimera::detect_chimeric(
            &seqlist,
            &clusters,
            max_len(&seqlist),
            params.shared,
            params.percent,
            params.abratio,
            abundance,
            &mut log,
        );
        for cs in &chistat {
            clusters[cs.index as usize].chi_head = cs.head;
            clusters[cs.index as usize].chi_tail = cs.tail;
        }
        log.push_str(&format!(
            "Number of chimeric clusters found: {}\n",
            chistat.len()
        ));
    }
    if abundance < 0 {
        abundance = 1;
    }
    let mut above = 0i32;
    let mut below = 0i32;
    for c in &clusters {
        if c.chi_head == c.chi_tail {
            if c.abundance < abundance {
                below += 1;
            } else {
                above += 1;
            }
        }
    }
    log.push_str(&format!(
        "Number of clusters with abundance above the cutoff (={}): {}\n",
        abundance, above
    ));
    log.push_str(&format!(
        "Number of clusters with abundance below the cutoff (={}): {}\n",
        abundance, below
    ));
    log.push_str("Writing clusters to files ...\n");

    let (reps_from_full, clstr, clstr2) =
        write_clusters(&seqlist, &mut clusters, params.deslen, !seqlist_modified);

    let mut reps_r1 = reps_from_full;
    let mut reps_r2 = None;

    // R2 output: restore full-length R2 reads by index, then emit reps.
    if let Some((b2, _)) = params.input2 {
        log.push_str("Write R2 reads\n");
        let r2 = parse_fastaq(b2, false);
        for k in 0..seqlist.len() {
            seqlist[k].seq = r2[k].seq.clone();
            seqlist[k].des = r2[k].des.clone();
            seqlist[k].qs = r2[k].qs.clone();
        }
        reps_r2 = Some(write_reps_only(&seqlist, &clusters));
    }
    // R1 output: if modified, restore full-length R1 reads by index, then emit.
    if seqlist_modified {
        log.push_str("Write R1 reads\n");
        let r1 = parse_fastaq(input, false);
        for k in 0..seqlist.len() {
            seqlist[k].seq = r1[k].seq.clone();
            seqlist[k].des = r1[k].des.clone();
            seqlist[k].qs = r1[k].qs.clone();
        }
        reps_r1 = write_reps_only(&seqlist, &clusters);
    }

    log.push_str("Done!\n");
    Ok(DupOutput {
        log,
        reps_r1,
        reps_r2,
        clstr,
        clstr2,
    })
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.len() < needle.len() {
        return None;
    }
    (0..=haystack.len() - needle.len()).find(|&i| &haystack[i..i + needle.len()] == needle)
}

fn parse_leading_int(s: &[u8]) -> i32 {
    let mut i = 0;
    let mut sign = 1i32;
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        if s[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    let mut v = 0i32;
    while i < s.len() && s[i].is_ascii_digit() {
        v = v * 10 + (s[i] - b'0') as i32;
        i += 1;
    }
    sign * v
}
