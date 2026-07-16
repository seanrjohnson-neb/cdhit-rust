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

/// Replace the leading run of digits in a member line with `n`
/// (port of the Perl `s/^\d+/$n/`).
fn replace_leading_int(raw: &[u8], n: usize) -> Vec<u8> {
    let mut i = 0;
    while i < raw.len() && raw[i].is_ascii_digit() {
        i += 1;
    }
    let mut out = n.to_string().into_bytes();
    out.extend_from_slice(&raw[i..]);
    out
}

/// `clstr_renumber.pl`: renumber clusters (from 0) and each cluster's member
/// indices (from 0), leaving the rest of each line intact.
pub fn renumber(input: &[u8]) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for (i, c) in clusters.iter().enumerate() {
        out.extend_from_slice(format!(">Cluster {}\n", i).as_bytes());
        for (j, m) in c.members.iter().enumerate() {
            out.extend_from_slice(&replace_leading_int(&m.raw, j));
            out.push(b'\n');
        }
    }
    out
}

/// `clstr_select.pl min max`: print (verbatim) clusters whose size is within
/// `[min, max]`.
pub fn select(input: &[u8], min: usize, max: usize) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let no = c.size();
        if no >= min && no <= max {
            out.extend_from_slice(&c.header_raw);
            out.push(b'\n');
            for m in &c.members {
                out.extend_from_slice(&m.raw);
                out.push(b'\n');
            }
        }
    }
    out
}

/// `clstr_cut.pl N`: keep only the top `N` members of each cluster (the
/// representative is always kept). Clusters are printed verbatim (original
/// header and member lines).
pub fn cut(input: &[u8], no_cutoff: i64) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let mut kept: Vec<&super::Member> = Vec::new();
        let mut no = 0i64;
        let mut repout = 1i64;
        for m in &c.members {
            if no < no_cutoff - repout || m.is_rep {
                if m.is_rep {
                    repout = 0;
                }
                kept.push(m);
                no += 1;
            }
        }
        if no > 0 {
            out.extend_from_slice(&c.header_raw);
            out.push(b'\n');
            for m in kept {
                out.extend_from_slice(&m.raw);
                out.push(b'\n');
            }
        }
    }
    out
}

/// `clstr_rep.pl`: print `cluster_id<TAB>rep_id<TAB>size` for each cluster.
/// The reference only recognises protein (`aa`) representative lines; a
/// representative line that is not `aa` is a format error.
pub fn rep(input: &[u8]) -> Result<Vec<u8>, String> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let mut rep_id: Vec<u8> = Vec::new();
        for m in &c.members {
            if m.is_rep {
                if find(&m.raw, b"aa, >").is_none() {
                    return Err(format!(
                        "format error {}",
                        String::from_utf8_lossy(&m.raw)
                    ));
                }
                rep_id = m.id.clone();
            }
        }
        out.extend_from_slice(format!("{}\t", c.num).as_bytes());
        out.extend_from_slice(&rep_id);
        out.extend_from_slice(format!("\t{}\n", c.size()).as_bytes());
    }
    Ok(out)
}

/// `clstr_merge.pl <master> <div1> [div2 ...]`: merge divided cluster files
/// into the master. The cluster order must be identical across files (each
/// master cluster is matched to the div cluster with the same representative);
/// the div clusters' extra members are appended, renumbered continuously.
pub fn merge(master: &[u8], divs: &[&[u8]]) -> Vec<u8> {
    let master_clusters = parse_clstr(master);
    let div_clusters: Vec<Vec<Cluster>> = divs.iter().map(|d| parse_clstr(d)).collect();
    let mut div_idx = vec![0usize; divs.len()];
    let mut out = Vec::new();

    let rep_id = |c: &Cluster| -> Option<Vec<u8>> {
        c.members.iter().find(|m| m.is_rep).map(|m| m.id.clone())
    };

    for mc in &master_clusters {
        let master_rep = match rep_id(mc) {
            Some(r) => r,
            None => continue,
        };
        // print the master cluster verbatim.
        out.extend_from_slice(&mc.header_raw);
        out.push(b'\n');
        for m in &mc.members {
            out.extend_from_slice(&m.raw);
            out.push(b'\n');
        }
        let mut rep_no = mc.members.len();

        for (i, dc_list) in div_clusters.iter().enumerate() {
            // advance this div file until a cluster with the matching rep.
            while div_idx[i] < dc_list.len() {
                let dc = &dc_list[div_idx[i]];
                div_idx[i] += 1;
                if rep_id(dc).as_deref() == Some(master_rep.as_slice()) {
                    if dc.members.len() > 1 {
                        let mut j = rep_no;
                        for m in &dc.members {
                            if m.is_rep {
                                continue;
                            }
                            out.extend_from_slice(&replace_leading_int(&m.raw, j));
                            out.push(b'\n');
                            j += 1;
                        }
                    }
                    rep_no += dc.members.len() - 1;
                    break;
                }
            }
        }
    }
    out
}

/// `clstr_reps_faa_rev.pl <clstr> <fasta> <cutoff>`: output the FASTA, keeping
/// for each cluster only its first `cutoff` sequences (the representative plus
/// up to `cutoff-1` others); the rest are dropped.
pub fn reps_faa_rev(clstr: &[u8], fasta: &[u8], cutoff: usize) -> Vec<u8> {
    // Build the set of ids to skip.
    let mut skip: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    for c in parse_clstr(clstr) {
        // gis: representative first, then the others in order.
        let mut gis: Vec<Vec<u8>> = Vec::new();
        for m in &c.members {
            if m.is_rep {
                gis.insert(0, m.id.clone());
            } else {
                gis.push(m.id.clone());
            }
        }
        for id in gis.iter().skip(cutoff) {
            skip.insert(id.clone());
        }
    }

    let mut out = Vec::new();
    let mut flag = false;
    for line in fasta.split_inclusive(|&b| b == b'\n') {
        let l = strip_line(line);
        if l.first() == Some(&b'>') {
            // id = substring after '>' up to first whitespace.
            let after = &l[1..];
            let end = after
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .unwrap_or(after.len());
            let gi = after[..end].to_vec();
            flag = !skip.contains(&gi);
        }
        if flag {
            out.extend_from_slice(l);
            out.push(b'\n');
        }
    }
    out
}

/// `clstr_rev.pl <file90> <file80>`: rebuild a coarse cluster file (`file80`,
/// clustered from `file90`'s representatives) so its members are the original
/// finer members from `file90`. Used to flatten a two-level hierarchical
/// clustering back onto the original sequences.
pub fn rev(file90: &[u8], file80: &[u8]) -> Vec<u8> {
    // gi2clstr[rep_id] = the finer cluster's member lines (only for clusters
    // with more than one member and a representative).
    let mut gi2clstr: std::collections::HashMap<Vec<u8>, Vec<Vec<u8>>> =
        std::collections::HashMap::new();
    for c in parse_clstr(file90) {
        if c.members.len() > 1 {
            if let Some(rep) = c.members.iter().find(|m| m.is_rep) {
                gi2clstr.insert(
                    rep.id.clone(),
                    c.members.iter().map(|m| m.raw.clone()).collect(),
                );
            }
        }
    }

    let mut out = Vec::new();
    let mut no = 0usize;
    for line in file80.split_inclusive(|&b| b == b'\n') {
        let l = strip_line(line);
        if l.first() == Some(&b'>') {
            out.extend_from_slice(line);
            no = 0;
        } else if let Some(id) = member_id(l) {
            let rep = l.last() == Some(&b'*');
            // identity = text after " at ", else "100%".
            let iden: Vec<u8> = match find_sub(l, b" at ") {
                Some(p) => l[p + 4..].to_vec(),
                None => b"100%".to_vec(),
            };
            if let Some(child_lines) = gi2clstr.get(&id) {
                for a in child_lines {
                    let mut a2 = replace_leading_int(a, no);
                    if !rep {
                        if a2.last() == Some(&b'*') {
                            // s/\*/at $iden/  (first '*')
                            if let Some(star) = a2.iter().position(|&b| b == b'*') {
                                let mut rep_line = a2[..star].to_vec();
                                rep_line.extend_from_slice(b"at ");
                                rep_line.extend_from_slice(&iden);
                                rep_line.extend_from_slice(&a2[star + 1..]);
                                a2 = rep_line;
                            }
                        } else if let Some(p) = find_sub(&a2, b"at ") {
                            // s/at (.+)$/at $iden,$1/  (first "at ")
                            let tail = a2[p + 3..].to_vec();
                            let mut nl = a2[..p].to_vec();
                            nl.extend_from_slice(b"at ");
                            nl.extend_from_slice(&iden);
                            nl.push(b',');
                            nl.extend_from_slice(&tail);
                            a2 = nl;
                        }
                    }
                    out.extend_from_slice(&a2);
                    out.push(b'\n');
                    no += 1;
                }
            } else {
                out.extend_from_slice(&replace_leading_int(l, no));
                out.push(b'\n');
                no += 1;
            }
        } else {
            out.extend_from_slice(line);
        }
    }
    out
}

fn strip_line(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

/// The id between `>` and `...` in a member line (mirrors `>(.+)\.\.\.`).
fn member_id(l: &[u8]) -> Option<Vec<u8>> {
    let gt = l.iter().position(|&b| b == b'>')?;
    let rest = &l[gt + 1..];
    // greedy up to the last "..."
    let mut last = None;
    let mut i = 0;
    while i + 3 <= rest.len() {
        if &rest[i..i + 3] == b"..." {
            last = Some(i);
        }
        i += 1;
    }
    last.map(|d| rest[..d].to_vec())
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    (0..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
}

/// `clstr_reduce.pl <segs> <rate>`: partition clusters into size-segments
/// (`segs` like `"1-1,2-5,6-100"`) and keep every `rate`-th cluster within each
/// segment. Clusters are emitted verbatim.
pub fn reduce(input: &[u8], segs: &str, rate: i64) -> Vec<u8> {
    // no2seg_idx[size] = segment index; uncovered sizes map to 0 (Perl undef).
    let mut no2seg: Vec<usize> = Vec::new();
    let seg_list: Vec<&str> = segs.split(',').collect();
    let nseg = seg_list.len();
    for (i, seg) in seg_list.iter().enumerate() {
        let (b, e) = if let Some((b, e)) = seg.split_once('-') {
            (b.parse().unwrap_or(0), e.parse().unwrap_or(0))
        } else {
            let v = seg.parse().unwrap_or(0);
            (v, v)
        };
        for j in b..=e {
            let j = j as usize;
            if j >= no2seg.len() {
                no2seg.resize(j + 1, 0);
            }
            no2seg[j] = i;
        }
    }
    let mut segs_no = vec![0i64; nseg.max(1)];
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let this_no = c.size();
        let seg = no2seg.get(this_no).copied().unwrap_or(0);
        if segs_no[seg] % rate == 0 {
            out.extend_from_slice(&c.header_raw);
            out.push(b'\n');
            for m in &c.members {
                out.extend_from_slice(&m.raw);
                out.push(b'\n');
            }
        }
        segs_no[seg] += 1;
    }
    out
}

/// Format a float the way Perl stringifies numbers (`%.15g`): 15 significant
/// digits, trailing zeros stripped. Cleans up float error (e.g. `1-0.8` prints
/// as `0.2`).
fn perl_g(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let p = 15i32;
    let neg = x < 0.0;
    let ax = x.abs();
    let exp = ax.log10().floor() as i32;
    let mut s = if exp >= -4 && exp < p {
        let prec = (p - 1 - exp).max(0) as usize;
        format!("{:.*}", prec, ax)
    } else {
        format!("{:.*e}", (p - 1) as usize, ax)
    };
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    if neg {
        format!("-{}", s)
    } else {
        s
    }
}

/// `clstr2tree.pl <clstr> <fr>`: emit a Newick-like tree. `fr` is the identity
/// fraction of this level (e.g. 0.8), the branch length inside a cluster is
/// `1-fr`.
pub fn to_tree(input: &[u8], fr: &str) -> Vec<u8> {
    let fr_val: f64 = fr.parse().unwrap_or(0.0);
    let fra = perl_g(1.0 - fr_val);
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    out.extend_from_slice(b"(\n");
    for (ci, c) in clusters.iter().enumerate() {
        let last = ci + 1 == clusters.len();
        let comma = if last { "" } else { "," };
        if c.members.len() == 1 {
            let rep = c
                .members
                .iter()
                .find(|m| m.is_rep)
                .map(|m| m.id.clone())
                .unwrap_or_default();
            out.extend_from_slice(&rep);
            out.extend_from_slice(format!(":1.0{}\n", comma).as_bytes());
        } else {
            out.extend_from_slice(b"(\n");
            let mms: Vec<String> = c
                .members
                .iter()
                .map(|m| format!("{}:{}", String::from_utf8_lossy(&m.id), fra))
                .collect();
            out.extend_from_slice(mms.join(",\n").as_bytes());
            out.extend_from_slice(format!("\n):{}{}\n", fr, comma).as_bytes());
        }
    }
    out.extend_from_slice(b");\n");
    out
}

/// Parse the leading numeric prefix of a byte string the way Perl coerces a
/// string to a number (e.g. `"97.50%"` → 97.5, `"+"` → 0, `""` → 0).
fn perl_num(s: &[u8]) -> f64 {
    let t = std::str::from_utf8(s).unwrap_or("");
    let b = t.trim_start().as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut has_digits = false;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        has_digits = true;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            has_digits = true;
        }
    }
    if has_digits && i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let save = i;
        i += 1;
        if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            i += 1;
        }
        let mut edig = false;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            edig = true;
        }
        if !edig {
            i = save;
        }
    }
    if !has_digits {
        return 0.0;
    }
    std::str::from_utf8(&b[..i]).unwrap().parse().unwrap_or(0.0)
}

/// `cd-hit-clstr_2_blm8.pl`: convert a `.clstr` into BLAST tabular (m8) rows,
/// one per non-representative member vs its representative.
///
/// This reproduces the reference behaviour *including its quirks*: protein
/// (`aa`) clusters and strand-only nucleotide lines lack alignment coordinates,
/// so the C++ Perl produces deterministic "garbage" (negative alignment
/// lengths, a raw string in the `q_b` column) via Perl's numeric coercion —
/// which is matched here exactly.
pub fn to_blm8(input: &[u8]) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let rep = c
            .members
            .iter()
            .find(|m| m.is_rep)
            .map(|m| m.id.clone())
            .unwrap_or_default();
        for m in &c.members {
            if m.is_rep {
                continue;
            }
            // last whitespace-delimited field of the line.
            let last_field = m
                .raw
                .split(|&b| b == b' ' || b == b'\t')
                .filter(|s| !s.is_empty())
                .next_back()
                .unwrap_or(b"");
            let mms: Vec<&[u8]> = last_field.split(|&b| b == b'/').collect();
            let a = mms.first().copied().unwrap_or(b"");
            // identity = last "/"-field with its trailing char chopped (the %).
            let iden_raw = mms.last().copied().unwrap_or(b"");
            let iden_str: &[u8] = if iden_raw.is_empty() {
                iden_raw
            } else {
                &iden_raw[..iden_raw.len() - 1]
            };
            let iden = perl_num(iden_str);
            // split coordinates by ':'; missing fields print as empty.
            let coords: Vec<&[u8]> = a.split(|&b| b == b':').collect();
            let field = |i: usize| coords.get(i).copied().unwrap_or(b"");
            let (qb, qe, sb, se) = (field(0), field(1), field(2), field(3));
            let alnln = perl_num(qe) - perl_num(qb) + 1.0;
            let mis = (alnln * (100.0 - iden) / 100.0).trunc();
            let bit = alnln * 3.0 - mis * 6.0;
            out.extend_from_slice(&m.id);
            out.push(b'\t');
            out.extend_from_slice(&rep);
            out.push(b'\t');
            out.extend_from_slice(iden_str);
            out.extend_from_slice(format!("\t{}\t{}\t0\t", perl_g(alnln), perl_g(mis)).as_bytes());
            out.extend_from_slice(qb);
            out.push(b'\t');
            out.extend_from_slice(qe);
            out.push(b'\t');
            out.extend_from_slice(sb);
            out.push(b'\t');
            out.extend_from_slice(se);
            out.extend_from_slice(format!("\t0\t{}\n", perl_g(bit)).as_bytes());
        }
    }
    out
}

fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    (0..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
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
