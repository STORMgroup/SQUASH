//! The raw-evidence guard of the rounds' partner searches: its view of the uncorrected tokens (kept in a
//! node-local file between searches) and its rule.

use super::*;

/// The guard's view of an uncorrected token: value and start only (8 bytes; the guard never reads a
/// token's end), and only tokens whose value is kept (`kept_values`): the guard looks every donor token up among the
/// target's kept values, so the others never count (`guard_tokens`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GuardTok {
    pub value: u32,
    pub start: i32,
}

/// The guard's view of uncorrected reads: each read's tokens with a kept value, as (value, start), in order. Each
/// read's view is stored at its exact length (counted first) as its values and the steps between consecutive starts
/// (u16; a step that does not fit, or goes back, is `u16::MAX` with the start in `far`), 6 bytes a token; it is read in
/// order by `iter`, as every guard helper reads it.
pub fn guard_tokens(tokens: &[Vec<Token>], kept: &HashSet<u32>) -> Vec<GuardRead> {
    tokens.par_iter().map(|read| GuardRead::new(read.iter().filter(|t| kept.contains(&t.value)))).collect()
}

/// One read's guard view (`guard_tokens`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GuardRead {
    pub(super) values: Vec<u32>,
    pub(super) steps: Vec<u16>,
    far: Vec<i32>,
}

impl GuardRead {
    pub(super) fn new<'a>(kept: impl Iterator<Item = &'a Token> + Clone) -> GuardRead {
        let n = kept.clone().count();
        let (mut values, mut steps, mut far) = (Vec::with_capacity(n), Vec::with_capacity(n), Vec::new());
        let mut previous = 0i64;
        for t in kept {
            let step = t.s() - previous;
            if (0..u16::MAX as i64).contains(&step) {
                steps.push(step as u16);
            } else {
                steps.push(u16::MAX);
                far.push(<i32 as std::convert::TryFrom<i64>>::try_from(t.s()).expect("token start outside i32"));
            }
            values.push(t.value);
            previous = t.s();
        }
        far.shrink_to_fit();
        GuardRead { values, steps, far }
    }

    /// The read's kept tokens, in order.
    pub fn iter(&self) -> impl Iterator<Item = GuardTok> + '_ {
        let mut start = 0i64;
        let mut far = self.far.iter();
        self.values.iter().zip(&self.steps).map(move |(&value, &step)| {
            start = if step == u16::MAX { *far.next().expect("a far start") as i64 } else { start + step as i64 };
            GuardTok { value, start: start as i32 }
        })
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Heap bytes (memory report).
    pub fn bytes(&self) -> usize {
        self.values.capacity() * 4 + self.steps.capacity() * 2 + self.far.capacity() * 4 + 3 * std::mem::size_of::<Vec<u8>>()
    }
}

/// The correction rounds need the guard's view (`guard_tokens`) only while they search partners, and it never changes
/// (uncorrected tokens), so it is written once to a node-local file and read back before each search
/// (`load_guard_view`). Layout: the read count, then per read its value / far counts and the three arrays,
/// little-endian.
pub fn spill_guard_view(path: &Path, view: &[GuardRead]) -> Result<()> {
    let tmp = path.with_extension("tmp");
    let mut out = BufWriter::with_capacity(1 << 24, File::create(&tmp)?);
    out.write_all(&(view.len() as u64).to_le_bytes())?;
    for read in view {
        out.write_all(&(read.values.len() as u64).to_le_bytes())?;
        out.write_all(&(read.far.len() as u64).to_le_bytes())?;
        for v in &read.values {
            out.write_all(&v.to_le_bytes())?;
        }
        for s in &read.steps {
            out.write_all(&s.to_le_bytes())?;
        }
        for f in &read.far {
            out.write_all(&f.to_le_bytes())?;
        }
    }
    out.flush()?;
    drop(out);
    std::fs::rename(&tmp, path)?;
    Ok(())
}

/// The guard's view from `spill_guard_view`'s file, each read's arrays at their exact length (as `guard_tokens` makes
/// them).
pub fn load_guard_view(path: &Path) -> Result<Vec<GuardRead>> {
    let mut f = BufReader::with_capacity(1 << 24, File::open(path)?);
    let mut word = [0u8; 8];
    let mut read_u64 = |f: &mut BufReader<File>| -> Result<usize> {
        f.read_exact(&mut word)?;
        Ok(u64::from_le_bytes(word) as usize)
    };
    let n = read_u64(&mut f)?;
    let mut view = Vec::with_capacity(n);
    let mut buf: Vec<u8> = Vec::new();
    for _ in 0..n {
        let (len, n_far) = (read_u64(&mut f)?, read_u64(&mut f)?);
        buf.resize(4 * len, 0);
        f.read_exact(&mut buf)?;
        let values: Vec<u32> = buf.chunks_exact(4).map(|b| u32::from_le_bytes(b.try_into().unwrap())).collect();
        buf.resize(2 * len, 0);
        f.read_exact(&mut buf)?;
        let steps: Vec<u16> = buf.chunks_exact(2).map(|b| u16::from_le_bytes(b.try_into().unwrap())).collect();
        buf.resize(4 * n_far, 0);
        f.read_exact(&mut buf)?;
        let far: Vec<i32> = buf.chunks_exact(4).map(|b| i32::from_le_bytes(b.try_into().unwrap())).collect();
        view.push(GuardRead { values, steps, far });
    }
    Ok(view)
}


/// The raw-evidence guard's rule (`guard_supported`): at least `min_anchors` uncorrected anchors in the implied
/// overlap, within max(`band`, `band_fraction` x the overlap) q-symbols of the chain's offset.
#[derive(Clone, Copy, Debug)]
pub struct Guard {
    pub min_anchors: usize,
    pub band: i64,
    pub band_fraction: f64,
}

impl Default for Guard {
    fn default() -> Self {
        Guard { min_anchors: 2, band: 150, band_fraction: 0.15 }
    }
}


/// The guard's uncorrected seed starts of one target: kept raw value -> its q starts.
pub(super) fn guard_starts(raw: impl IntoIterator<Item = GuardTok>, guard_kept: &HashSet<u32>) -> HashMap<u32, Vec<i64>> {
    let mut starts: HashMap<u32, Vec<i64>> = HashMap::new();
    for t in raw {
        if guard_kept.contains(&t.value) {
            starts.entry(t.value).or_default().push(t.start as i64);
        }
    }
    starts
}

/// The raw-evidence guard of one chain: true when at least `g.min_anchors` uncorrected anchors (every pair of
/// equal kept uncorrected values of the two reads) lie in the implied overlap [ts, te) on the target, on a diagonal
/// within max(`g.band`, `g.band_fraction` x the overlap) of the chain's `offset` (q-symbols). Computed for this donor
/// only, and only for chains that passed the other filters.
pub(super) fn guard_supported(target_starts: &HashMap<u32, Vec<i64>>, donor_raw: impl IntoIterator<Item = GuardTok>, ts: i64, te: i64, offset: i64, g: &Guard) -> bool {
    let band = f64::max(g.band as f64, g.band_fraction * ((te - ts) as f64));
    let mut count = 0usize;
    for t in donor_raw {
        if let Some(starts) = target_starts.get(&t.value) {
            for &tq in starts {
                if ts <= tq && tq < te && (((t.start as i64 - tq) - offset).abs() as f64) <= band {
                    count += 1;
                    if count >= g.min_anchors {
                        return true;
                    }
                }
            }
        }
    }
    count >= g.min_anchors
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The per-pair guard (`guard_supported` on one donor) keeps the chains of a per-target guard that collects every
    /// donor's uncorrected anchors before chaining.
    #[test]
    fn per_pair_guard_keeps_the_chains_of_the_per_target_guard() {
        let (lengths, tokens) = noisy_reads(10, 1);
        let (_, raw) = noisy_reads(60, 1);
        let p = Params { max_reads: 0, min_anchors: 4, min_overlap: 500, min_cover: 0.1, repeat_fraction: 0.01, ..Params::default() };
        let g = Guard { min_anchors: 3, band: 20, band_fraction: 0.01 };
        let index = build_index(&tokens, p.max_reads, p.repeat_fraction);
        let raw_index = build_index(&raw, p.max_reads, p.repeat_fraction);
        let guard_kept = kept_values(&raw, p.max_reads, p.repeat_fraction);
        let guard_view = guard_tokens(&raw, &guard_kept);   // the guard's view
        let key = |c: &Chain| (c.donor, c.anchors, c.t0, c.t1, c.ts, c.te, c.ds, c.de);
        let mut scratch = Scratch::new(tokens.len());
        let (mut unguarded_total, mut guarded_total) = (0, 0);
        for r in 0..tokens.len() {
            let unguarded = target_chains(r, &lengths, &tokens, &index, &p, None, &mut scratch, None);
            let mut support: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
            for t in &raw[r] {
                for (d, j) in raw_index.get(t.value) {
                    if d as usize != r {
                        support.entry(d as usize).or_default().push((t.s(), raw[d as usize][j as usize].s()));
                    }
                }
            }
            let mut expected: Vec<_> = unguarded
                .iter()
                .filter(|c| {
                    let band = f64::max(g.band as f64, g.band_fraction * ((c.te - c.ts) as f64));
                    let offset = c.ds - c.ts;
                    support.get(&c.donor).map_or(0, |a| {
                        a.iter().filter(|&&(tq, dq)| c.ts <= tq && tq < c.te && (((dq - tq) - offset).abs() as f64) <= band).count()
                    }) >= g.min_anchors
                })
                .map(key)
                .collect();
            let mut observed: Vec<_> =
                target_chains(r, &lengths, &tokens, &index, &p, Some((&guard_view, &guard_kept, &g)), &mut scratch, None).iter().map(key).collect();
            expected.sort();
            observed.sort();
            assert_eq!(observed, expected, "target {}", r);
            unguarded_total += unguarded.len();
            guarded_total += observed.len();
        }
        assert!(0 < guarded_total && guarded_total < unguarded_total, "{} of {}", guarded_total, unguarded_total);
        let mut keys: Vec<u32> = raw_index.values().collect();
        let mut kept: Vec<u32> = guard_kept.into_iter().collect();
        keys.sort_unstable();
        kept.sort_unstable();
        assert_eq!(keys, kept);
    }

    /// The guard's view written to the local file comes back the same, far starts included.
    #[test]
    fn guard_view_file_round_trips() {
        let (_, mut tokens) = noisy_reads(10, 6);
        // far starts: steps that do not fit u16 and a step back are stored apart
        tokens[3] = (0..40)
            .map(|i| {
                let start = if i == 20 { 5 } else { 70_000 * (i % 25) as i64 };
                Token::new(i as u32 % 9, start, start + 9)
            })
            .collect();
        let kept = kept_values(&tokens, 0, 0.0);
        let view = guard_tokens(&tokens, &kept);
        assert!(view.iter().any(|v| !v.far.is_empty()));
        let path = std::env::temp_dir().join(format!("guard_view_round_trip_{}.bin", std::process::id()));
        spill_guard_view(&path, &view).unwrap();
        let back = load_guard_view(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(back, view);
        assert!(back.iter().zip(&view).all(|(a, b)| a.bytes() == b.bytes()));
    }
}
