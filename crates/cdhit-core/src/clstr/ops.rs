//! Ports of individual CD-HIT `clstr_*` Perl post-processing scripts. Each
//! function takes the `.clstr` bytes and returns the exact stdout bytes the
//! corresponding Perl script would produce.

use super::{parse_clstr, Cluster};

/// `clstr_sort_by.pl`: re-emit clusters sorted by member count (`"no"`,
/// default) or by representative length (`"len"`), renumbered from 0.
pub fn sort_by(input: &[u8], sort_by_what: &str) -> Vec<u8> {
    let mut clusters = parse_clstr(input);
    // Perl sort is stable; sort keys are descending.
    match sort_by_what {
        "no" => clusters.sort_by(|a, b| b.size().cmp(&a.size())),
        "len" => clusters.sort_by(|a, b| b.max_len().cmp(&a.max_len())),
        _ => {}
    }
    let mut out = Vec::new();
    for (i, c) in clusters.iter().enumerate() {
        out.extend_from_slice(format!(">Cluster {}\n", i).as_bytes());
        for m in &c.members {
            out.extend_from_slice(&m.raw);
            out.push(b'\n');
        }
    }
    out
}

/// `clstr_size_stat.pl`: distribution of cluster sizes.
pub fn size_stat(input: &[u8]) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let max_size = clusters.iter().map(|c| c.size()).max().unwrap_or(0);
    let mut counts = vec![0i64; max_size + 1];
    for c in &clusters {
        counts[c.size()] += 1;
    }
    let mut out = Vec::new();
    out.extend_from_slice(b"size\tNo.clstr\tNo.seq\n");
    for i in 1..=max_size {
        let noc = counts[i];
        if noc == 0 {
            continue;
        }
        let nos = noc * i as i64;
        out.extend_from_slice(format!("{}\t{}\t{}\n", i, noc, nos).as_bytes());
    }
    out
}

/// `clstr_size_histogram.pl`: histogram of cluster sizes binned by `step`
/// (default 100). Bin index is `int((size-1)/step)`.
pub fn size_histogram(input: &[u8], step: i64) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut bins: Vec<i64> = Vec::new();
    for c in &clusters {
        let this_no = c.size() as i64;
        let bin = ((this_no - 1) / step) as usize;
        if bin >= bins.len() {
            bins.resize(bin + 1, 0);
        }
        bins[bin] += 1;
    }
    let mut out = Vec::new();
    out.extend_from_slice(b"bin_size\tNo_of_clusters\n");
    for (i, count) in bins.iter().enumerate() {
        let i = i as i64;
        out.extend_from_slice(
            format!("{}-{}\t{}\n", i * step + 1, i * step + step, count).as_bytes(),
        );
    }
    out
}

/// `clstr2txt.pl`: tabular per-sequence view of the cluster file.
pub fn to_txt(input: &[u8]) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    out.extend_from_slice(b"id\tclstr\tclstr_size\tlength\tclstr_rep\tclstr_iden\tclstr_cov\n");
    for c in &clusters {
        process_cluster_txt(&mut out, c);
    }
    out
}

fn process_cluster_txt(out: &mut Vec<u8>, c: &Cluster) {
    let no = c.members.len();
    // sort by (is_rep desc, len desc) — stable, as in Perl.
    let mut t: Vec<&super::Member> = c.members.iter().collect();
    t.sort_by(|a, b| {
        (b.is_rep as i32)
            .cmp(&(a.is_rep as i32))
            .then(b.len.cmp(&a.len))
    });
    let longest = t.iter().filter(|m| m.is_rep).map(|m| m.len).max().unwrap_or(0);
    for m in &t {
        let cov = if longest != 0 {
            (m.len as f64 / longest as f64 * 100.0) as i64
        } else {
            0
        };
        let rep_flag = if m.is_rep { 1 } else { 0 };
        // rep identity is "100" (no %), members carry the parsed "...%" string.
        let iden: &[u8] = if m.is_rep { b"100" } else { &m.iden };
        out.extend_from_slice(&m.id);
        out.extend_from_slice(format!("\t{}\t{}\t{}\t{}\t", c.num, no, m.len, rep_flag).as_bytes());
        out.extend_from_slice(iden);
        out.extend_from_slice(format!("\t{}%\n", cov).as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &[u8] = b">Cluster 0\n0\t150nt, >seqA... *\n1\t150nt, >seqB... at 1:150:1:150/+/98.00%\n>Cluster 1\n0\t99nt, >seqC... *\n";

    #[test]
    fn size_stat_counts() {
        let out = size_stat(SAMPLE);
        let s = String::from_utf8(out).unwrap();
        assert!(s.starts_with("size\tNo.clstr\tNo.seq\n"));
        assert!(s.contains("1\t1\t1\n")); // one cluster of size 1
        assert!(s.contains("2\t1\t2\n")); // one cluster of size 2
    }

    #[test]
    fn sort_by_no_orders_and_renumbers() {
        let out = sort_by(SAMPLE, "no");
        let s = String::from_utf8(out).unwrap();
        // the 2-member cluster comes first, renumbered to 0.
        assert!(s.starts_with(">Cluster 0\n0\t150nt, >seqA... *\n"));
    }
}
