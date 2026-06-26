//! Tie-averaged ranking of pooled values, matching SciPy `_rankdata` with
//! `method='average'`.
//!
//! SciPy sorts the pooled data with a stable argsort and assigns every member
//! of a tie group the average of the ordinal ranks it spans (`first + (count-1)/2`,
//! 1-based). Fligner-Killeen needs the ranks back in pooled (input) order so each
//! can be turned into a normal score, so this returns the full per-element rank
//! vector rather than per-group aggregates.

/// Tie-averaged ranks of `values`, returned in the input order. Sorting is
/// stable on the value (NaN-free domain) to mirror SciPy's argsort.
#[must_use]
pub fn rankdata_average(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));

    let mut ranks = vec![0.0_f64; n];
    let mut i = 0;
    while i < n {
        let v = values[order[i]];
        let mut j = i + 1;
        while j < n && values[order[j]] == v {
            j += 1;
        }
        let count = j - i;
        let avg_rank = (i + 1) as f64 + (count as f64 - 1.0) / 2.0;
        for &idx in &order[i..j] {
            ranks[idx] = avg_rank;
        }
        i = j;
    }
    ranks
}

#[cfg(test)]
mod tests {
    use super::rankdata_average;

    #[test]
    fn distinct_values() {
        let r = rankdata_average(&[3.0, 1.0, 2.0]);
        assert_eq!(r, vec![3.0, 1.0, 2.0]);
    }

    #[test]
    fn tie_averaging() {
        // scipy.stats.rankdata([1, 2, 2, 2, 1, 3]) == [1.5, 4, 4, 4, 1.5, 6]
        let r = rankdata_average(&[1.0, 2.0, 2.0, 2.0, 1.0, 3.0]);
        assert_eq!(r, vec![1.5, 4.0, 4.0, 4.0, 1.5, 6.0]);
    }
}
