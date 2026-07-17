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

/// True if a member line contains the `(aa|nt), >` form the Perl regexes
/// require (`(\d+)(aa|nt), >(.+)\.\.\.`). Lines that don't match are skipped
/// by scripts like `clstr_sort_prot_by.pl`.
fn has_len_id(raw: &[u8]) -> bool {
    find(raw, b"aa, >").is_some() || find(raw, b"nt, >").is_some()
}

/// Drop a leading run of digits followed by a single tab (Perl `s/^\d+\t//`).
fn strip_leading_index_tab(raw: &[u8]) -> Vec<u8> {
    let mut i = 0;
    while i < raw.len() && raw[i].is_ascii_digit() {
        i += 1;
    }
    if i > 0 && i < raw.len() && raw[i] == b'\t' {
        raw[i + 1..].to_vec()
    } else {
        raw.to_vec()
    }
}

/// `clstr_select_rep.pl <min> <max>`: for each cluster whose member count is
/// within `[min, max]`, print the representative id (an empty line if the
/// cluster has no representative). A representative line that is not `aa`/`nt`
/// is a format error, as upstream.
pub fn select_rep(input: &[u8], min: usize, max: usize) -> Result<Vec<u8>, String> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        let no = c.size();
        let mut rep: Vec<u8> = Vec::new();
        for m in &c.members {
            if m.is_rep {
                // /\*$/ then require /\d+(aa|nt), >(.+)\.\.\./
                if !has_len_id(&m.raw) || m.id.is_empty() {
                    return Err(format!("format error {}", String::from_utf8_lossy(&m.raw)));
                }
                rep = m.id.clone();
            }
        }
        if no >= min && no <= max {
            out.extend_from_slice(&rep);
            out.push(b'\n');
        }
    }
    Ok(out)
}

/// `clstr_sort_prot_by.pl [len|id|<other>]`: sort the members *within* each
/// cluster, re-emitting the cluster header verbatim then `<i>\t<line>` for each
/// member (the leading `\d+\t` index is stripped, then re-added as the new
/// rank). Representatives sort first (their length is treated as 99999999).
///
/// Sort keys (all stable, matching Perl's stable mergesort):
/// - `"len"` (default): length descending.
/// - `"id"`: id ascending, then length descending.
/// - anything else: length descending, then id ascending.
pub fn sort_prot_by(input: &[u8], sort_by: &str) -> Vec<u8> {
    let clusters = parse_clstr(input);
    let mut out = Vec::new();
    for c in &clusters {
        // (line-without-index, len, id)
        let mut items: Vec<(Vec<u8>, i64, &[u8])> = Vec::new();
        for m in &c.members {
            if !has_len_id(&m.raw) || m.id.is_empty() {
                continue;
            }
            // Perl sets len to 99999999 when the line matches /\*/ (anywhere);
            // in practice that is the representative line.
            let len = if find(&m.raw, b"*").is_some() {
                99999999
            } else {
                m.len
            };
            items.push((strip_leading_index_tab(&m.raw), len, m.id.as_slice()));
        }
        match sort_by {
            "len" => items.sort_by(|a, b| b.1.cmp(&a.1)),
            "id" => items.sort_by(|a, b| a.2.cmp(b.2).then(b.1.cmp(&a.1))),
            _ => items.sort_by(|a, b| b.1.cmp(&a.1).then(a.2.cmp(b.2))),
        }
        out.extend_from_slice(&c.header_raw);
        out.push(b'\n');
        for (i, (line, _, _)) in items.iter().enumerate() {
            out.extend_from_slice(format!("{}\t", i).as_bytes());
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

/// `clstr_merge_noorder.pl <master> <div1> [div2 ...]`: like `clstr_merge` but
/// the div clusters may appear in any order — each div file is read fully into
/// a `rep_id -> [member lines]` map, then master clusters are printed verbatim
/// with the matching div members appended and renumbered continuously.
pub fn merge_noorder(master: &[u8], divs: &[&[u8]]) -> Vec<u8> {
    // slave[rep_id] = accumulated non-representative member lines (verbatim).
    let mut slave: std::collections::HashMap<Vec<u8>, Vec<Vec<u8>>> =
        std::collections::HashMap::new();
    for d in divs {
        for c in parse_clstr(d) {
            let rep = c.members.iter().find(|m| m.is_rep).map(|m| m.id.clone());
            let members: Vec<Vec<u8>> = c
                .members
                .iter()
                .filter(|m| !m.is_rep)
                .map(|m| m.raw.clone())
                .collect();
            if !members.is_empty() {
                // Perl dies here if there is no rep; we skip unmatched members.
                if let Some(rep) = rep {
                    slave.entry(rep).or_default().extend(members);
                }
            }
        }
    }

    let mut out = Vec::new();
    for mc in parse_clstr(master) {
        let master_rep = match mc.members.iter().find(|m| m.is_rep) {
            Some(m) => m.id.clone(),
            None => continue, // Perl only prints clusters that have a rep.
        };
        out.extend_from_slice(&mc.header_raw);
        out.push(b'\n');
        for m in &mc.members {
            out.extend_from_slice(&m.raw);
            out.push(b'\n');
        }
        let mut rep_no = mc.members.len();
        if let Some(list) = slave.get(&master_rep) {
            for line in list {
                out.extend_from_slice(&replace_leading_int(line, rep_no));
                out.push(b'\n');
                rep_no += 1;
            }
        }
    }
    out
}

/// Split an id on the first-and-subsequent `||` separators (Perl
/// `split(/\|\|/, $id)`).
fn split_bars(id: &[u8]) -> Vec<&[u8]> {
    let mut parts = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i + 1 < id.len() {
        if id[i] == b'|' && id[i + 1] == b'|' {
            parts.push(&id[start..i]);
            i += 2;
            start = i;
        } else {
            i += 1;
        }
    }
    parts.push(&id[start..]);
    parts
}

/// `clstr_quality_eval_by_link.pl`: sensitivity/specificity of the clustering
/// vs a benchmark encoded in each id as `>seq_id||benchmark_id`. Uses
/// independent links (n-1 per group of n) rather than all pairs. Output is
/// deterministic (only sums over hash values are used).
pub fn quality_eval_by_link(input: &[u8]) -> Result<Vec<u8>, String> {
    let mut bench_count: std::collections::HashMap<Vec<u8>, i64> =
        std::collections::HashMap::new();
    let mut total_cdhit_links = 0i64;
    let mut correct_links = 0i64;

    for c in parse_clstr(input) {
        let mut clstr_by_ben: std::collections::HashMap<Vec<u8>, i64> =
            std::collections::HashMap::new();
        let mut t_no = 0i64;
        for m in &c.members {
            if m.id.is_empty() || !has_len_id(&m.raw) {
                continue;
            }
            let parts = split_bars(&m.id);
            let ben_id = parts.get(1).copied().unwrap_or(b"");
            // Perl truthiness: skip if undef/empty/"0".
            if ben_id.is_empty() || ben_id == b"0" {
                continue;
            }
            *bench_count.entry(ben_id.to_vec()).or_insert(0) += 1;
            *clstr_by_ben.entry(ben_id.to_vec()).or_insert(0) += 1;
            t_no += 1;
        }
        if t_no > 1 {
            for cnt in clstr_by_ben.values() {
                correct_links += cnt - 1;
            }
            total_cdhit_links += t_no - 1;
        }
    }

    let mut total_bench_links = 0i64;
    for n in bench_count.values() {
        total_bench_links += n - 1;
    }

    if total_bench_links == 0 || total_cdhit_links == 0 {
        return Err("Illegal division by zero".to_string());
    }
    let sen = correct_links as f64 / total_bench_links as f64;
    let spe = correct_links as f64 / total_cdhit_links as f64;

    let mut out = Vec::new();
    out.extend_from_slice(format!("Total benchmark links\t{}\n", total_bench_links).as_bytes());
    out.extend_from_slice(format!("Total cd-hit links\t{}\n", total_cdhit_links).as_bytes());
    out.extend_from_slice(format!("Total correct links\t{}\n", correct_links).as_bytes());
    out.extend_from_slice(format!("Sensitivity\t{}\n", perl_g(sen)).as_bytes());
    out.extend_from_slice(format!("Specificity\t {}\n", perl_g(spe)).as_bytes());
    Ok(out)
}

/// Parse a decimal integer prefix the way Perl coerces a string to a number,
/// truncated to an integer (used for the segment bounds of `plot_len1`).
fn perl_int(s: &[u8]) -> i64 {
    perl_num(s) as i64
}

/// `plot_len1.pl <clstr> <segs> <len_segs>`: a tabular count of sequences and
/// clusters, cross-tabulated by cluster-size segment (`segs`, e.g.
/// `"1,2-5,6-up"`) against representative-length segment (`len_segs`, e.g.
/// `"1-100,101-200"`). Despite the name it emits no plot, just the table.
pub fn plot_len1(input: &[u8], segs: &str, len_segs: &str) -> Vec<u8> {
    // clstr_nos[size] = #clusters of that size; clstr_len[size] = rep lengths.
    let mut clstr_nos: std::collections::HashMap<usize, i64> = std::collections::HashMap::new();
    let mut clstr_len: std::collections::HashMap<usize, Vec<i64>> =
        std::collections::HashMap::new();
    let mut max_no = 0usize;
    // `this_len` persists across clusters in the Perl (only reset by a rep line).
    let mut this_len = 0i64;
    for c in parse_clstr(input) {
        let this_no = c.members.len();
        for m in &c.members {
            if m.is_rep && has_len_id(&m.raw) {
                this_len = m.len;
            }
        }
        *clstr_nos.entry(this_no).or_insert(0) += 1;
        if this_no > max_no {
            max_no = this_no;
        }
        clstr_len.entry(this_no).or_default().push(this_len);
    }

    let seg_list: Vec<&str> = if segs.is_empty() {
        Vec::new()
    } else {
        segs.split(',').collect()
    };
    let len_seg_list: Vec<&str> = if len_segs.is_empty() {
        Vec::new()
    } else {
        len_segs.split(',').collect()
    };

    let nos = |j: usize| clstr_nos.get(&j).copied().unwrap_or(0);
    let empty: Vec<i64> = Vec::new();
    let lens_at = |j: usize| clstr_len.get(&j).unwrap_or(&empty);

    let mut out = Vec::new();
    out.extend_from_slice(b"Size\tNo. seq\tNo. clstr");
    let mut tlen_nos = vec![0i64; len_seg_list.len()];
    for ls in &len_seg_list {
        out.extend_from_slice(format!("\t{}", ls).as_bytes());
    }
    out.push(b'\n');

    // Parse a length segment into inclusive [b, e] bounds.
    let len_bounds = |ls: &str| -> (i64, i64) {
        if let Some((b, e)) = ls.split_once('-') {
            (perl_int(b.as_bytes()), perl_int(e.as_bytes()))
        } else {
            let v = perl_int(ls.as_bytes());
            (v, v)
        }
    };

    let mut tno = 0i64;
    let mut tno1 = 0i64;
    for seg in &seg_list {
        let mut lens: Vec<i64> = Vec::new();
        if seg.contains('-') {
            let (b_str, e_str) = seg.split_once('-').unwrap();
            let b = perl_int(b_str.as_bytes());
            let e = if e_str.to_ascii_lowercase().contains("up") {
                max_no as i64
            } else {
                perl_int(e_str.as_bytes())
            };
            let mut no = 0i64;
            let mut no1 = 0i64;
            let mut j = b;
            while j <= e {
                if j >= 0 {
                    let ju = j as usize;
                    no += j * nos(ju);
                    no1 += nos(ju);
                    lens.extend_from_slice(lens_at(ju));
                }
                j += 1;
            }
            tno += no;
            tno1 += no1;
            out.extend_from_slice(format!("{}\t{}\t{}", seg, no, no1).as_bytes());
        } else {
            let s = perl_int(seg.as_bytes());
            let su = if s >= 0 { s as usize } else { usize::MAX };
            let count = if su == usize::MAX { 0 } else { nos(su) };
            if su != usize::MAX {
                lens.extend_from_slice(lens_at(su));
            }
            tno += s * count;
            tno1 += count;
            // Third column is interpolated `$clstr_nos[$seg]`: empty if unseen.
            let field3 = match clstr_nos.get(&su) {
                Some(v) => v.to_string(),
                None => String::new(),
            };
            out.extend_from_slice(format!("{}\t{}\t{}", seg, s * count, field3).as_bytes());
        }
        for (j, ls) in len_seg_list.iter().enumerate() {
            let (lb, le) = len_bounds(ls);
            let cnt = lens.iter().filter(|&&t| t >= lb && t <= le).count() as i64;
            out.extend_from_slice(format!("\t{}", cnt).as_bytes());
            tlen_nos[j] += cnt;
        }
        out.push(b'\n');
    }
    out.extend_from_slice(format!("Total\t{}\t{}", tno, tno1).as_bytes());
    for v in &tlen_nos {
        out.extend_from_slice(format!("\t{}", v).as_bytes());
    }
    out.push(b'\n');
    out
}

/// Split a line on tabs, dropping trailing empty fields (Perl's default
/// `split(/\t/, $line)` behaviour).
fn split_tab_trim(line: &[u8]) -> Vec<&[u8]> {
    let mut parts: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    while parts.last() == Some(&(b"".as_slice())) {
        parts.pop();
    }
    parts
}

/// `clstr_sql_tbl_sort.pl <table_file> <level>`: stable-sort a tab-separated
/// table by numeric columns counted from the end (last-2 ascending, last
/// descending; for level 2/3 the preceding column pairs break ties). Returns
/// an error ("error level") if the first row has fewer than `level*2+2`
/// columns.
pub fn sql_tbl_sort(input: &[u8], level: i64) -> Result<Vec<u8>, String> {
    let mut lls: Vec<Vec<u8>> = input
        .split_inclusive(|&b| b == b'\n')
        .map(|l| strip_line(l).to_vec())
        .collect();

    let first_cols = lls.first().map(|l| split_tab_trim(l).len()).unwrap_or(0);
    if (first_cols as i64) < level * 2 + 2 {
        return Err("error level".to_string());
    }

    // Numeric value of the column `neg` positions from the end (Perl $x[-neg]);
    // undef (out of range) coerces to 0.
    let field = |line: &[u8], neg: usize| -> f64 {
        let cols = split_tab_trim(line);
        if neg <= cols.len() {
            perl_num(cols[cols.len() - neg])
        } else {
            0.0
        }
    };
    let numcmp = |x: f64, y: f64| x.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Equal);

    lls.sort_by(|a, b| {
        // level >= 1: a[-2] asc, then b[-1] desc.
        let mut ord = numcmp(field(a, 2), field(b, 2));
        if ord == std::cmp::Ordering::Equal {
            ord = numcmp(field(b, 1), field(a, 1));
        }
        if level >= 2 && ord == std::cmp::Ordering::Equal {
            ord = numcmp(field(a, 4), field(b, 4));
            if ord == std::cmp::Ordering::Equal {
                ord = numcmp(field(b, 3), field(a, 3));
            }
        }
        if level >= 3 && ord == std::cmp::Ordering::Equal {
            ord = numcmp(field(a, 6), field(b, 6));
            if ord == std::cmp::Ordering::Equal {
                ord = numcmp(field(b, 5), field(a, 5));
            }
        }
        ord
    });

    let mut out = Vec::new();
    for l in &lls {
        out.extend_from_slice(l);
        out.push(b'\n');
    }
    Ok(out)
}

/// Output of [`dup_pe_out`]: the two filtered paired-end files.
pub struct PeOut {
    pub out1: Vec<u8>,
    pub out2: Vec<u8>,
}

/// Extract an id from a FASTA/FASTQ header line the way `cd-hit-dup-PE-out.pl`
/// does: drop the leading `>`/`@`, Perl `chop` (remove the final byte, normally
/// the newline), then truncate at the first whitespace.
fn pe_header_id(line: &[u8]) -> &[u8] {
    if line.is_empty() {
        return b"";
    }
    let after = &line[1..];
    let after = if after.is_empty() {
        after
    } else {
        &after[..after.len() - 1] // chop
    };
    let end = after
        .iter()
        .position(|&c| c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' || c == 0x0b || c == 0x0c)
        .unwrap_or(after.len());
    &after[..end]
}

/// `cd-hit-dup-PE-out.pl -i R1 -j R2 -c clstr -o out1 -p out2`: export the
/// representative paired-end reads (those whose R1 *or* R2 id is a cluster
/// representative in the `.clstr`) into two parallel output files. Detects
/// FASTA vs FASTQ from the first byte of `in1`.
pub fn dup_pe_out(clstr: &[u8], in1: &[u8], in2: &[u8]) -> PeOut {
    // Representative ids from `*` lines (Perl requires `\s(\d+)(aa|nt), >`).
    let mut rep_ids: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    for c in parse_clstr(clstr) {
        for m in &c.members {
            if m.is_rep && has_len_id(&m.raw) && !m.id.is_empty() {
                rep_ids.insert(m.id.clone());
            }
        }
    }

    let lines1: Vec<&[u8]> = in1.split_inclusive(|&b| b == b'\n').collect();
    let lines2: Vec<&[u8]> = in2.split_inclusive(|&b| b == b'\n').collect();

    let is_fasta = in1.first() == Some(&b'>');
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();

    let is_rep = |line: &[u8]| -> bool {
        let ida = pe_header_id(line);
        rep_ids.contains(ida)
    };

    if is_fasta {
        let mut flag = false;
        let n = lines1.len().min(lines2.len());
        for i in 0..n {
            let a = lines1[i];
            let b = lines2[i];
            if a.first() == Some(&b'>') && b.first() == Some(&b'>') {
                flag = is_rep(a) || is_rep(b);
            }
            if flag {
                out1.extend_from_slice(a);
                out2.extend_from_slice(b);
            }
        }
    } else {
        let mut ia = 0usize;
        let mut ib = 0usize;
        while ia < lines1.len() && ib < lines2.len() {
            let a = lines1[ia];
            ia += 1;
            let b = lines2[ib];
            ib += 1;
            if a.first() == Some(&b'@') && b.first() == Some(&b'@') {
                let flag = is_rep(a) || is_rep(b);
                if flag {
                    out1.extend_from_slice(a);
                    out2.extend_from_slice(b);
                    for _ in 0..3 {
                        if ia < lines1.len() {
                            out1.extend_from_slice(lines1[ia]);
                            ia += 1;
                        }
                    }
                    for _ in 0..3 {
                        if ib < lines2.len() {
                            out2.extend_from_slice(lines2[ib]);
                            ib += 1;
                        }
                    }
                }
            }
        }
    }
    PeOut { out1, out2 }
}

/// `clstr_sql_tbl.pl <clstr> <tbl>`: build (or extend) a hierarchical cluster
/// table. `existing` is `None` to create a fresh table from `clstr` (columns
/// `id  len  cid  rep`), or `Some(table)` to append two columns (`cid  rep`) at
/// the next clustering level — resolving each row's new cluster via its own id,
/// or, failing that, via the representative of its last-level cluster (as in
/// hierarchical `db90 -> db60 -> db30` runs).
pub fn sql_tbl(clstr: &[u8], existing: Option<&[u8]>) -> Result<Vec<u8>, String> {
    match existing {
        None => {
            // Create mode: one row per member, cluster index counted from 0.
            let mut out = Vec::new();
            for (cid, c) in parse_clstr(clstr).iter().enumerate() {
                for m in &c.members {
                    if !has_len_id(&m.raw) || m.id.is_empty() {
                        return Err(format!("format error {}", String::from_utf8_lossy(&m.raw)));
                    }
                    let rep = if m.is_rep { 1 } else { 0 };
                    out.extend_from_slice(&m.id);
                    out.extend_from_slice(format!("\t{}\t{}\t{}\n", m.len, cid, rep).as_bytes());
                }
            }
            Ok(out)
        }
        Some(table) => {
            // Append mode: map each id to its new-level cluster and rep flag.
            let mut id2cid: std::collections::HashMap<Vec<u8>, i64> =
                std::collections::HashMap::new();
            let mut idisrep: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
            for (cid, c) in parse_clstr(clstr).iter().enumerate() {
                for m in &c.members {
                    if !has_len_id(&m.raw) || m.id.is_empty() {
                        return Err(format!("format error {}", String::from_utf8_lossy(&m.raw)));
                    }
                    if m.is_rep {
                        idisrep.insert(m.id.clone());
                    }
                    id2cid.insert(m.id.clone(), cid as i64);
                }
            }

            // Table rows, newline-stripped (Perl `chop`).
            let rows: Vec<Vec<u8>> = table
                .split_inclusive(|&b| b == b'\n')
                .map(|l| strip_line(l).to_vec())
                .filter(|l| !l.is_empty())
                .collect();

            // last_cid_2_id[last_level_cid] = the row id whose last-level rep==1.
            let mut last_cid_2_id: std::collections::HashMap<i64, Vec<u8>> =
                std::collections::HashMap::new();
            for row in &rows {
                let cols = split_tab_trim(row);
                if cols.len() < 2 {
                    continue;
                }
                let last_cid = perl_num(cols[cols.len() - 2]) as i64;
                let last_rep = perl_num(cols[cols.len() - 1]);
                if last_rep == 1.0 {
                    last_cid_2_id.insert(last_cid, cols[0].to_vec());
                }
            }

            let mut out = Vec::new();
            for row in &rows {
                let cols = split_tab_trim(row);
                let id = cols.first().copied().unwrap_or(b"");
                let last_cid = if cols.len() >= 2 {
                    perl_num(cols[cols.len() - 2]) as i64
                } else {
                    0
                };
                let this_rep = if idisrep.contains(id) { 1 } else { 0 };
                let this_cid = if let Some(&c) = id2cid.get(id) {
                    c
                } else {
                    // Fall back to the representative of the last-level cluster.
                    match last_cid_2_id.get(&last_cid).and_then(|r| id2cid.get(r)) {
                        Some(&c) => c,
                        None => return Err(format!("at {}", String::from_utf8_lossy(row))),
                    }
                };
                out.extend_from_slice(row);
                out.extend_from_slice(format!("\t{}\t{}\n", this_cid, this_rep).as_bytes());
            }
            Ok(out)
        }
    }
}

/// One per-cluster FASTA file produced by [`make_multi_seq`]: the cluster id
/// (used as the file name) and the file's contents.
pub struct MultiSeqFile {
    pub cid: Vec<u8>,
    pub content: Vec<u8>,
}

/// The first non-whitespace token after a leading `>` (Perl `/^>(\S+)/`).
/// Returns `None` if the line is not a `>`-header with a token.
fn fasta_defline_id(line: &[u8]) -> Option<&[u8]> {
    if line.first() != Some(&b'>') {
        return None;
    }
    let after = &line[1..];
    let end = after
        .iter()
        .position(|&c| c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' || c == 0x0b || c == 0x0c)
        .unwrap_or(after.len());
    if end == 0 {
        None
    } else {
        Some(&after[..end])
    }
}

/// `make_multi_seq.pl <fasta> <clstr> <out_dir> <size_cutoff>`: emit one FASTA
/// file per cluster whose member count is `>= size_cutoff`, named by the
/// cluster number, containing that cluster's sequences drawn from `fasta`.
/// Files are returned in the order their clusters first appear in `fasta`
/// (native-only: it produces many files rather than a single stream).
pub fn make_multi_seq(
    fasta: &[u8],
    clstr: &[u8],
    size_cutoff: usize,
) -> Result<Vec<MultiSeqFile>, String> {
    // id -> cluster number, only for clusters at or above the size cutoff.
    let mut id2cid: std::collections::HashMap<Vec<u8>, i64> = std::collections::HashMap::new();
    for c in parse_clstr(clstr) {
        // Header must be `>Cluster <n>`; parse_cluster_num yields -1 otherwise.
        if c.num < 0 {
            return Err(format!(
                "Wrong format {}",
                String::from_utf8_lossy(&c.header_raw)
            ));
        }
        for m in &c.members {
            if !has_len_id(&m.raw) || m.id.is_empty() {
                return Err(format!("Wrong format {}", String::from_utf8_lossy(&m.raw)));
            }
        }
        if c.size() >= size_cutoff {
            for m in &c.members {
                id2cid.insert(m.id.clone(), c.num);
            }
        }
    }

    // Stream the FASTA, routing each record to its cluster's buffer.
    let mut order: Vec<i64> = Vec::new();
    let mut buffers: std::collections::HashMap<i64, Vec<u8>> = std::collections::HashMap::new();
    let mut cur: Option<i64> = None;
    let mut flag = false;
    for line in fasta.split_inclusive(|&b| b == b'\n') {
        if let Some(id) = fasta_defline_id(line) {
            if let Some(&cid) = id2cid.get(id) {
                buffers.entry(cid).or_insert_with(|| {
                    order.push(cid);
                    Vec::new()
                });
                cur = Some(cid);
                flag = true;
            } else {
                flag = false;
            }
        }
        if flag {
            if let Some(cid) = cur {
                buffers.get_mut(&cid).unwrap().extend_from_slice(line);
            }
        }
    }

    Ok(order
        .into_iter()
        .map(|cid| MultiSeqFile {
            cid: cid.to_string().into_bytes(),
            content: buffers.remove(&cid).unwrap(),
        })
        .collect())
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
