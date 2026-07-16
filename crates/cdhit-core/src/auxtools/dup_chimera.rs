//! Chimera detection for `cd-hit-dup` (`-f` / `-s`).
//!
//! Port of `SortByAbundance` and `DetectChimeric` from `cdhit-dup.cxx`.
//!
//! CAVEAT: the reference `DetectChimeric` contains out-of-bounds reads (e.g. it
//! always inspects the top *two* offsets even when only one exists) — undefined
//! behaviour whose result is not portable. This port guards those accesses, so
//! its output matches the reference on well-formed cases but may diverge in the
//! degenerate ones the C++ leaves undefined. The core de-duplication path
//! ([`super::dup::cd_hit_dup`] without `-f`) is bit-for-bit verified.

use super::bioseq::Sequence;
use super::dup::{
    compare_prefix, compare_prefix2, compare_suffix, count_mismatch, count_mismatch2,
    count_mismatch3, make_hash_off, DupCluster,
};

pub struct ChimericSource {
    pub index: i32,
    pub head: i32,
    pub tail: i32,
}

/// Counting sort of clusters by descending abundance (port of
/// `SortByAbundance`).
pub fn sort_by_abundance(clusters: &mut Vec<DupCluster>) {
    let n = clusters.len();
    if n <= 1 {
        return;
    }
    let mut max = clusters[0].abundance;
    let mut min = clusters[0].abundance;
    for c in clusters.iter() {
        max = max.max(c.abundance);
        min = min.min(c.abundance);
    }
    let m = (max - min + 1) as usize;
    let mut count = vec![0usize; m];
    let mut accum = vec![0usize; m];
    let mut offset = vec![0usize; m];
    for c in clusters.iter() {
        count[(max - c.abundance) as usize] += 1;
    }
    for i in 1..m {
        accum[i] = accum[i - 1] + count[i - 1];
    }
    let mut sorted: Vec<Option<DupCluster>> = (0..n).map(|_| None).collect();
    for c in clusters.drain(..) {
        let len = (max - c.abundance) as usize;
        let id = accum[len] + offset[len];
        offset[len] += 1;
        sorted[id] = Some(c);
    }
    *clusters = sorted.into_iter().map(|c| c.unwrap()).collect();
}

#[derive(Clone, Copy)]
struct HashHit {
    index: i32,
    offset: i32,
}

const MAX_OFFSETS: usize = 5;

#[derive(Clone)]
struct HashHit3 {
    index: i32,
    size: usize,
    offsets: [i32; MAX_OFFSETS],
    counts: [i32; MAX_OFFSETS],
}

impl HashHit3 {
    fn new(i: i32, o: i32) -> Self {
        let mut h = HashHit3 {
            index: i,
            size: 1,
            offsets: [0; MAX_OFFSETS],
            counts: [0; MAX_OFFSETS],
        };
        h.offsets[0] = o;
        h.counts[0] = 1;
        h
    }
    fn update(&mut self, o: i32) {
        let mut imin = 0usize;
        let mut min = self.counts[0];
        for i in 0..self.size {
            if o == self.offsets[i] {
                self.counts[i] += 1;
                return;
            } else if self.counts[i] < min {
                min = self.counts[i];
                imin = i;
            }
        }
        if self.size == MAX_OFFSETS {
            self.offsets[imin] = o;
            self.counts[imin] = 1;
        } else {
            self.offsets[self.size] = o;
            self.counts[self.size] = 1;
            self.size += 1;
        }
    }
}

#[derive(Clone, Copy)]
struct ParentInfo {
    index: i32,
    start: i32,
    len: i32,
    offset: i32,
}

/// Port of `DetectChimeric`. Returns the list of clusters identified as
/// chimeric together with their two parent clusters.
#[allow(clippy::too_many_arguments)]
pub fn detect_chimeric(
    seqlist: &[Sequence],
    clusters: &[DupCluster],
    max: i32,
    shared: i32,
    percent: f32,
    abratio: f32,
    minabu: i32,
    log: &mut String,
) -> Vec<ChimericSource> {
    let n = clusters.len();
    let mut chistat: Vec<ChimericSource> = Vec::new();
    if n <= 2 {
        return chistat;
    }
    let mut shared = shared;
    if shared < 20 {
        shared = 20;
    }

    let rep_seq = |ci: usize| -> &[u8] { &seqlist[clusters[ci].members[0]].seq };

    // primer = common prefix across cluster representatives.
    let mut primer = rep_seq(0).len() as i32;
    for i in 1..n {
        primer = compare_prefix2(rep_seq(0), 0, rep_seq(i), 0, primer, 0);
        if primer < 15 {
            primer = 0;
            break;
        }
    }
    log.push_str(&format!("primer = {}\n", primer));
    log.push_str("Searching for chimeric clusters ...\n");

    let mut hash_table: std::collections::HashMap<u32, Vec<HashHit>> =
        std::collections::HashMap::new();
    let mut chimap: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
    let mut hit_mapping = vec![0i32; n];

    for i in 0..n {
        let seq = rep_seq(i);
        let qlen = seq.len() as i32;
        let qn = clusters[i].abundance;
        let maxmm = (0.01 * qlen as f32 * percent) as i32;
        let maxmm2 = (0.015 * qlen as f32 * percent) as i32;

        if i > 0 && i % 1000 == 0 {
            log.push_str(&format!(
                "Checked {:9} clusters, detected {:9} chimeric clusters\n",
                i,
                chistat.len()
            ));
        }
        if qn < minabu {
            break;
        }

        // hashes[j] = MakeHash(seq, j, shared) for j in [primer, qlen-shared].
        let mut hashes = vec![0u32; primer.max(0) as usize];
        let mut j = primer;
        while j <= qlen - shared {
            hashes.push(make_hash_off(seq, j, shared));
            j += 1;
        }

        if qlen < 2 * shared + primer {
            let mut j = primer;
            while (j as usize) < hashes.len() {
                hash_table
                    .entry(hashes[j as usize])
                    .or_default()
                    .push(HashHit { index: i as i32, offset: j });
                j += 1;
            }
            continue;
        }

        let mut hit_list: Vec<HashHit3> = Vec::new();

        let mut detected = false;
        let mut start = primer;
        let mut proto: i32 = -1;
        let mut protomm = max;

        while (start + shared) <= qlen {
            let hash = hashes[start as usize];
            start += 1;
            if let Some(hits) = hash_table.get(&hash) {
                for hh in hits {
                    let hit = hh.index;
                    let offset = hh.offset - start + 1;
                    if hit_mapping[hit as usize] == 0 {
                        hit_list.push(HashHit3::new(hit, offset));
                        hit_mapping[hit as usize] = hit_list.len() as i32;
                    } else {
                        let idx = hit_mapping[hit as usize] as usize - 1;
                        hit_list[idx].update(offset);
                    }
                }
            }
        }

        let mut parent_a: Vec<ParentInfo> = Vec::new();
        let mut parent_b: Vec<ParentInfo> = Vec::new();

        for hitj in &hit_list {
            let hit = hitj.index as usize;
            let rep = rep_seq(hit);
            if clusters[hit].abundance < (qn as f32 * abratio) as i32 {
                continue;
            }

            // Top offsets by hit count (count desc, then offset asc).
            let mut oc: Vec<(i32, i32)> = (0..hitj.size)
                .map(|k| (hitj.offsets[k], hitj.counts[k]))
                .collect();
            oc.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

            let mut tail = 0i32;
            let head = compare_prefix(seq, rep, maxmm) as i32;
            if head >= primer + shared {
                parent_a.push(ParentInfo {
                    index: hit as i32,
                    start: 0,
                    len: head,
                    offset: 0,
                });
            }

            let mut best = oc[0].0;
            for k in 0..2 {
                if k >= oc.len() {
                    break; // guard: C++ reads OOB here when only one offset exists
                }
                best = oc[k].0;
                if qlen + best < rep.len() as i32 {
                    tail = compare_suffix(seq, qlen - 1, rep, qlen + best - 1, 0, maxmm);
                    if tail >= shared {
                        parent_b.push(ParentInfo {
                            index: hit as i32,
                            start: qlen - tail + best,
                            len: tail,
                            offset: best,
                        });
                        break;
                    }
                } else {
                    let rlen = rep.len() as i32;
                    tail = compare_suffix(seq, rlen - best - 1, rep, rlen - 1, 0, maxmm);
                    if tail >= shared {
                        parent_b.push(ParentInfo {
                            index: hit as i32,
                            start: rlen - tail,
                            len: tail,
                            offset: best,
                        });
                        break;
                    }
                }
                tail = 0;
            }

            let mm = if head + tail >= qlen {
                let x = (head + qlen - tail) / 2;
                count_mismatch2(seq, rep, 0, x, maxmm2)
                    + count_mismatch3(seq, x, rep, x + best, qlen, maxmm2)
            } else {
                count_mismatch2(seq, rep, 0, qlen, maxmm2)
            };
            if mm <= maxmm2 && mm < protomm {
                protomm = mm;
                proto = hit as i32;
            }
        }

        // Reset touched hit_mapping entries.
        for hitj in &hit_list {
            hit_mapping[hitj.index as usize] = 0;
        }

        if protomm <= maxmm {
            if let Some(&csi) = chimap.get(&proto) {
                let (h, t) = (chistat[csi as usize].head, chistat[csi as usize].tail);
                chimap.insert(i as i32, chistat.len() as i32);
                chistat.push(ChimericSource {
                    index: i as i32,
                    head: h,
                    tail: t,
                });
            }
            continue;
        }

        // parentB sorted by (start - offset) ascending.
        parent_b.sort_by(|a, b| (a.start - a.offset).cmp(&(b.start - b.offset)));

        'outer: for info_a in &parent_a {
            let rep_a = rep_seq(info_a.index as usize);
            for info_b in &parent_b {
                let rep_b = rep_seq(info_b.index as usize);
                if info_a.index == info_b.index {
                    continue;
                }
                if (info_b.start - info_b.offset) > info_a.len {
                    break;
                }
                let xa = (info_a.len + info_b.start - info_b.offset) / 2;
                let xb = xa + info_b.offset;
                let lqa = count_mismatch(seq, rep_a, 0, xa);
                let rqa = count_mismatch(seq, rep_a, xa, qlen);
                let lqb = count_mismatch3(seq, 0, rep_b, 0, xa, qlen);
                let rqb = count_mismatch3(seq, xa, rep_b, xb, qlen - xa, qlen);
                let protomm_eff = if proto < 0 { maxmm2 } else { protomm };
                if lqb > (protomm_eff + maxmm)
                    && rqa > (protomm_eff + maxmm)
                    && (lqa + rqb) <= protomm_eff
                {
                    detected = true;
                    chimap.insert(i as i32, chistat.len() as i32);
                    chistat.push(ChimericSource {
                        index: i as i32,
                        head: info_a.index,
                        tail: info_b.index,
                    });
                    break 'outer;
                }
            }
        }

        if !detected {
            let mut j = primer;
            while (j as usize) < hashes.len() {
                hash_table
                    .entry(hashes[j as usize])
                    .or_default()
                    .push(HashHit { index: i as i32, offset: j });
                j += 1;
            }
        }
    }

    chistat
}
