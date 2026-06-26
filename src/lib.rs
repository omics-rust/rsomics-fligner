//! Fligner-Killeen test for equality of variances (median-centered by default;
//! also mean and trimmed mean) — a value-exact, faster `scipy.stats.fligner`.
//!
//! The test is non-parametric: each group is centered, the pooled absolute
//! deviations are jointly tie-averaged ranked, the ranks are turned into normal
//! scores, and the between-group score dispersion is compared to a χ² with
//! `k − 1` degrees of freedom. Robust to non-normality, unlike Bartlett's test.

mod fligner;
mod igamc;
mod ndtri;
mod rank;

use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

pub use fligner::{Center, FlignerResult, fligner};

/// Parse one single-column file into a group: one numeric value per non-empty
/// line. Used for the default mode (one file per group).
pub fn parse_column<R: BufRead>(reader: R) -> Result<Vec<f64>> {
    let mut group = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let s = line.trim();
        if s.is_empty() {
            continue;
        }
        let v: f64 = s.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {}: '{s}' is not a number", i + 1))
        })?;
        group.push(v);
    }
    if group.is_empty() {
        return Err(RsomicsError::InvalidInput("no data in input".into()));
    }
    Ok(group)
}

/// Parse a long-format file: `value<TAB>group` per non-empty line. Groups are
/// returned in first-appearance order.
pub fn parse_long<R: BufRead>(reader: R) -> Result<Vec<Vec<f64>>> {
    let mut labels: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<f64>> = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = if trimmed.contains('\t') {
            trimmed.split('\t').map(str::trim).collect()
        } else {
            trimmed.split_whitespace().collect()
        };
        if f.len() != 2 {
            return Err(RsomicsError::InvalidInput(format!(
                "line {}: expected 'value<TAB>group'",
                i + 1
            )));
        }
        let value: f64 = f[0].parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {}: '{}' is not a number", i + 1, f[0]))
        })?;
        let label = f[1];
        let idx = match labels.iter().position(|l| l == label) {
            Some(j) => j,
            None => {
                labels.push(label.to_string());
                groups.push(Vec::new());
                groups.len() - 1
            }
        };
        groups[idx].push(value);
    }
    if groups.len() < 2 {
        return Err(RsomicsError::InvalidInput(
            "long input needs at least two distinct groups".into(),
        ));
    }
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parses_single_column() {
        let g = parse_column(Cursor::new("1\n2\n3\n")).unwrap();
        assert_eq!(g, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn column_skips_blank_lines() {
        let g = parse_column(Cursor::new("1\n\n2\n")).unwrap();
        assert_eq!(g, vec![1.0, 2.0]);
    }

    #[test]
    fn parses_long_layout() {
        let data = "1\ta\n2\tb\n3\ta\n4\tb\n";
        let g = parse_long(Cursor::new(data)).unwrap();
        assert_eq!(g, vec![vec![1.0, 3.0], vec![2.0, 4.0]]);
    }

    #[test]
    fn long_rejects_single_group() {
        let data = "1\ta\n2\ta\n";
        assert!(parse_long(Cursor::new(data)).is_err());
    }
}
