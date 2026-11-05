//! The rounds' inputs: the uncorrected tokens, the q-symbol counts and a partner table's donor rows.

use super::*;

fn int(field: &str) -> Result<i64> {
    field.trim().parse::<i64>().map_err(|e| format!("bad integer {:?}: {}", field, e).into())
}

/// The uncorrected tokens (a folder of token arrays): read ids, source lengths and tokens, through the partner
/// search's loader.
pub(super) fn load_prepared(dir: &Path) -> Result<(Vec<String>, Vec<i64>, Vec<Vec<Token>>)> {
    crate::partners::load_prepared(dir)
}

/// Each read's id and q-symbol count from a q-symbols file.
pub(super) fn symbol_counts(path: &Path) -> Result<(Vec<String>, Vec<i64>)> {
    let mut lines = open_text(path)?.lines();
    if lines.next().transpose()?.as_deref() != Some(crate::io::QSYMBOLS_TAG) {
        return Err(format!("bad q-symbols tag in {}", path.display()).into());
    }
    lines.next();
    match lines.next().transpose()?.as_deref() {
        Some("read_id\tsymbols") | Some("#columns=read_id\tsymbols") => {}
        _ => return Err(format!("bad q-symbols header in {}", path.display()).into()),
    }
    let (mut ids, mut counts) = (Vec::new(), Vec::new());
    for line in lines {
        let line = line?;
        let (id, field) = line.split_once('\t').ok_or("bad q-symbols row")?;
        ids.push(id.to_string());
        counts.push(if field.is_empty() { 0 } else { field.bytes().filter(|&b| b == b',').count() as i64 + 1 });
    }
    Ok((ids, counts))
}

/// A partner table's donor rows grouped by target ordinal, each target's sorted by (rank, donor_id), and the number of
/// rows. The table is read a line at a time, the donor's id compared through the catalog, and each target's rows kept
/// at their exact length.
pub(super) fn load_donors(path: &Path, index: &HashMap<&str, usize>, reads: usize) -> Result<(Vec<Vec<DonorRow>>, usize)> {
    let mut lines = open_text(path)?.lines();
    let header = lines.next().ok_or_else(|| format!("{}: empty table", path.display()))??;
    let columns: HashMap<&str, usize> = header.split('\t').enumerate().map(|(i, c)| (c, i)).collect();
    let c = |name: &str| columns.get(name).copied().ok_or_else(|| format!("missing column {}", name));
    let (c_t, c_d, c_rank) = (c("target_id")?, c("donor_id")?, c("rank")?);
    let (c_ts, c_te, c_ds, c_de) = (c("target_start")?, c("target_end")?, c("donor_start")?, c("donor_end")?);
    let mut grouped: Vec<Vec<(i64, DonorRow)>> = vec![Vec::new(); reads];
    let mut selected = 0usize;
    let need = *[c_t, c_d, c_rank, c_ts, c_te, c_ds, c_de].iter().max().unwrap();
    for line in lines {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= need {
            return Err(format!("{}: short row", path.display()).into());
        }
        let (target, donor) = match (index.get(f[c_t]), index.get(f[c_d])) {
            (Some(&t), Some(&d)) => (t, d),
            _ => return Err("selected donor endpoint is outside the token catalog".into()),
        };
        grouped[target].push((
            int(f[c_rank])?,
            DonorRow { donor, target_start: int(f[c_ts])?, target_end: int(f[c_te])?, donor_start: int(f[c_ds])?, donor_end: int(f[c_de])? },
        ));
        selected += 1;
    }
    // (rank, donor_id): the donor's id from the catalog (`index` maps each id to its ordinal, so ids[ordinal] is it)
    let mut ids: Vec<&str> = vec![""; reads];
    for (&id, &r) in index {
        ids[r] = id;
    }
    let grouped = grouped
        .into_iter()
        .map(|mut rows| {
            rows.sort_by(|a, b| (a.0, ids[a.1.donor]).cmp(&(b.0, ids[b.1.donor])));
            let mut out = Vec::with_capacity(rows.len());
            out.extend(rows.into_iter().map(|x| x.1));
            out
        })
        .collect();
    Ok((grouped, selected))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The donor table read a line at a time gives the rows of a reference reader (grouped by target, sorted by rank,
    /// then donor id as text), in the same order, each target's rows at their exact length.
    #[test]
    fn streamed_donor_table_equals_the_table_reader() {
        let ids: Vec<String> = ["r10", "r2", "r1", "a7", "r20", "b3"].iter().map(|x| x.to_string()).collect();
        let index: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, x)| (x.as_str(), i)).collect();
        let mut state = 5u64;
        let mut next = move |m: usize| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize % m
        };
        let mut text = String::from(crate::partners::COLUMNS);
        text.push('\n');
        let mut rows: Vec<Vec<String>> = Vec::new();
        for _ in 0..200 {
            let (t, d) = (next(ids.len()), next(ids.len()));
            let row: Vec<String> = vec![ids[t].clone(), ids[d].clone(), (1 + next(4)).to_string(), "+".into(), "5".into(), "5".into(),
                                        next(900).to_string(), next(900).to_string(), next(900).to_string(), next(900).to_string(),
                                        "5.0".into(), "0.5".into()];
            text.push_str(&row.join("\t"));
            text.push('\n');
            rows.push(row);
        }
        let path = std::env::temp_dir().join(format!("donor_table_{}.tsv", std::process::id()));
        std::fs::write(&path, &text).unwrap();
        let (got, n) = load_donors(&path, &index, ids.len()).unwrap();
        std::fs::remove_file(&path).unwrap();
        // the reference: grouping and (rank, donor id) sort on the parsed fields
        let header: Vec<&str> = crate::partners::COLUMNS.split('\t').collect();
        let col = |name: &str| header.iter().position(|c| *c == name).unwrap();
        let mut want: Vec<Vec<(i64, String, (usize, i64, i64, i64, i64))>> = vec![Vec::new(); ids.len()];
        for r in &rows {
            let v = |name: &str| r[col(name)].parse::<i64>().unwrap();
            want[index[r[col("target_id")].as_str()]].push((v("rank"), r[col("donor_id")].clone(),
                (index[r[col("donor_id")].as_str()], v("target_start"), v("target_end"), v("donor_start"), v("donor_end"))));
        }
        for w in want.iter_mut() {
            w.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
        }
        assert_eq!(n, rows.len());
        for (g, w) in got.iter().zip(&want) {
            let g2: Vec<_> = g.iter().map(|x| (x.donor, x.target_start, x.target_end, x.donor_start, x.donor_end)).collect();
            assert_eq!(g2, w.iter().map(|x| x.2).collect::<Vec<_>>());
            assert_eq!(g.capacity(), g.len());
        }
    }
}
