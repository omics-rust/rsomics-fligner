//! Fligner-Killeen test for equal variances — the rank-based statistic and its
//! chi-squared p-value.
//!
//! Each group is centered (median by default — the distribution-free variant —
//! or mean / trimmed mean). With `aᵢⱼ = |xᵢⱼ − centerᵢ|` pooled and jointly
//! tie-averaged ranked to rank `rᵢⱼ`, the normal scores are
//! `qᵢⱼ = Φ⁻¹(rᵢⱼ/(2(N+1)) + 0.5)`. With per-group mean score `Āᵢ`, grand mean
//! `Ā`, and `V = var(q, ddof=1)`,
//!
//! ```text
//! X²  =  Σ nᵢ (Āᵢ − Ā)² / V
//! ```
//!
//! and `p = P(χ²_{k−1} > X²)`. The reductions mirror NumPy's pairwise summation
//! and `scipy.stats.fligner`'s arithmetic so the result is value-exact.

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::igamc::chi2_sf;
use crate::ndtri::ndtri;
use crate::rank::rankdata_average;

/// How each group is centered before the absolute-deviation transform.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Center {
    Median,
    Mean,
    Trimmed,
}

/// The Fligner-Killeen result: the statistic, its degrees of freedom, p-value.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct FlignerResult {
    pub statistic: f64,
    pub df: f64,
    pub pvalue: f64,
}

/// NumPy `np.add.reduce` pairwise summation: blocks of ≤128 summed with an
/// 8-accumulator unrolled loop, larger spans split in half. Matching this is
/// what makes `np.mean` / `np.var` bit-identical at large N.
fn pairwise_sum(x: &[f64]) -> f64 {
    let n = x.len();
    if n <= 128 {
        if n == 0 {
            return 0.0;
        }
        if n < 8 {
            let mut s = x[0];
            for &v in &x[1..] {
                s += v;
            }
            return s;
        }
        let mut acc = [x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]];
        let mut i = 8;
        while i + 8 <= n {
            for k in 0..8 {
                acc[k] += x[i + k];
            }
            i += 8;
        }
        let mut s =
            ((acc[0] + acc[1]) + (acc[2] + acc[3])) + ((acc[4] + acc[5]) + (acc[6] + acc[7]));
        while i < n {
            s += x[i];
            i += 1;
        }
        s
    } else {
        // Split at the largest multiple of 8 not exceeding n/2.
        let mut half = n / 2;
        half -= half % 8;
        pairwise_sum(&x[..half]) + pairwise_sum(&x[half..])
    }
}

fn mean(x: &[f64]) -> f64 {
    pairwise_sum(x) / x.len() as f64
}

/// NumPy `np.var(x, ddof=1)`: divisor (n−1), the squared deviations summed with
/// the same pairwise reduction as the mean.
fn var_ddof1(x: &[f64], xmean: f64) -> f64 {
    let sq: Vec<f64> = x.iter().map(|&v| (v - xmean) * (v - xmean)).collect();
    pairwise_sum(&sq) / (x.len() as f64 - 1.0)
}

/// numpy-compatible median: average of the two middle order statistics for even
/// length, the single middle for odd. `select_nth_unstable` keeps it O(n).
fn median(values: &mut [f64]) -> f64 {
    let n = values.len();
    let mid = n / 2;
    let (_, hi, _) = values.select_nth_unstable_by(mid, f64::total_cmp);
    let upper = *hi;
    if n % 2 == 1 {
        upper
    } else {
        let lower = values[..mid]
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        (lower + upper) / 2.0
    }
}

/// scipy's `trim_mean`: drop `floor(p·n)` order statistics from each sorted end,
/// mean the rest. scipy rejects `lowercut == uppercut` upstream.
fn trimmed_mean(values: &mut [f64], proportiontocut: f64) -> Result<f64> {
    let n = values.len();
    let lowercut = (proportiontocut * n as f64) as usize;
    let uppercut = n - lowercut;
    if lowercut >= uppercut {
        return Err(RsomicsError::InvalidInput(
            "proportiontocut too large: nothing left after trimming".into(),
        ));
    }
    values.sort_unstable_by(f64::total_cmp);
    let kept = &values[lowercut..uppercut];
    Ok(mean(kept))
}

fn center_of(group: &[f64], center: Center, proportiontocut: f64) -> Result<f64> {
    match center {
        Center::Mean => Ok(mean(group)),
        Center::Median => {
            let mut scratch = group.to_vec();
            Ok(median(&mut scratch))
        }
        Center::Trimmed => {
            let mut scratch = group.to_vec();
            trimmed_mean(&mut scratch, proportiontocut)
        }
    }
}

/// Run the Fligner-Killeen test across `groups` (each a sample). Mirrors the
/// transform order of `scipy.stats.fligner` (`_morestats.py`).
pub fn fligner(groups: &[Vec<f64>], center: Center, proportiontocut: f64) -> Result<FlignerResult> {
    let k = groups.len();
    if k < 2 {
        return Err(RsomicsError::InvalidInput(
            "Fligner-Killeen test needs at least two groups".into(),
        ));
    }
    for (i, g) in groups.iter().enumerate() {
        if g.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "group {} is empty",
                i + 1
            )));
        }
    }

    let ni: Vec<usize> = groups.iter().map(Vec::len).collect();
    let ntot: usize = ni.iter().sum();

    // Pooled |xij - center_i| in group order, then jointly tie-averaged ranked.
    let mut pooled: Vec<f64> = Vec::with_capacity(ntot);
    for g in groups {
        let c = center_of(g, center, proportiontocut)?;
        pooled.extend(g.iter().map(|&y| (y - c).abs()));
    }

    let ranks = rankdata_average(&pooled);
    let denom = 2.0 * (ntot as f64 + 1.0);
    let scores: Vec<f64> = ranks.iter().map(|&r| ndtri(r / denom + 0.5)).collect();

    let abar = mean(&scores);
    let v2 = var_ddof1(&scores, abar);

    let mut numer = 0.0;
    let mut start = 0;
    for &n in &ni {
        let aibar = mean(&scores[start..start + n]);
        let d = aibar - abar;
        numer += n as f64 * (d * d);
        start += n;
    }

    let statistic = numer / v2;
    let df = (k - 1) as f64;
    let pvalue = chi2_sf(df, statistic);

    Ok(FlignerResult {
        statistic,
        df,
        pvalue,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) {
        let r = (a - b).abs() / b.abs().max(f64::MIN_POSITIVE);
        assert!(r <= tol, "got {a}, want {b}, rel {r:e}");
    }

    #[test]
    fn docstring_example_median() {
        // scipy.stats.fligner(a, b, c) from the fligner docstring.
        let a = vec![8.88, 9.12, 9.04, 8.98, 9.00, 9.08, 9.01, 8.85, 9.06, 8.99];
        let b = vec![8.88, 8.95, 9.29, 9.44, 9.15, 9.58, 8.36, 9.18, 8.67, 9.05];
        let c = vec![8.95, 9.12, 8.95, 8.85, 9.03, 8.84, 9.07, 8.98, 8.86, 8.98];
        let r = fligner(&[a, b, c], Center::Median, 0.05).unwrap();
        approx(r.pvalue, 0.004_508_260_800_047_75, 1e-12);
        assert_eq!(r.df, 2.0);
    }

    #[test]
    fn two_group_median() {
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let b = vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0];
        let r = fligner(&[a, b], Center::Median, 0.05).unwrap();
        assert_eq!(r.df, 1.0);
        assert!(r.statistic > 0.0 && r.pvalue >= 0.0 && r.pvalue <= 1.0);
    }
}
