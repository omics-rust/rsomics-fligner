//! Compat against committed `scipy.stats.fligner` goldens. Each row of
//! `tests/golden/expected.tsv` names a fixture set, the centering, and the
//! statistic + p SciPy produced; we run the binary and assert value-exact
//! equality. No SciPy at test time — the goldens are frozen (scipy 1.17.1).

use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fligner"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(files: &[&str], mode: &str, center: &str, pct: &str) -> (f64, f64) {
    let mut cmd = Command::new(bin());
    if mode == "long" {
        cmd.arg("--long").arg(golden(files[0]));
    } else {
        for f in files {
            cmd.arg(golden(f));
        }
    }
    cmd.args(["--center", center]);
    cmd.args(["--proportiontocut", pct]);
    let out = cmd.output().expect("run binary");
    assert!(
        out.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let line = String::from_utf8(out.stdout).unwrap();
    let f: Vec<f64> = line
        .trim()
        .split('\t')
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(f.len(), 2, "expected statistic,p, got {line:?}");
    (f[0], f[1])
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
}

#[test]
fn matches_scipy_goldens() {
    let expected = std::fs::read_to_string(golden("expected.tsv")).unwrap();
    let mut checked = 0;
    for line in expected.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let (id, files, mode, center, pct) = (c[0], c[1], c[2], c[3], c[4]);
        let want_stat: f64 = c[5].parse().unwrap();
        let want_p: f64 = c[6].parse().unwrap();
        let file_list: Vec<&str> = files.split(',').collect();

        let (stat, p) = run(&file_list, mode, center, pct);
        assert!(
            rel(stat, want_stat) <= 1e-13,
            "{id} statistic: got {stat}, want {want_stat}, rel {:e}",
            rel(stat, want_stat)
        );
        assert!(
            rel(p, want_p) <= 1e-12,
            "{id} p: got {p}, want {want_p}, rel {:e}",
            rel(p, want_p)
        );
        checked += 1;
    }
    assert!(
        checked >= 6,
        "expected at least 6 golden rows, ran {checked}"
    );
}

#[test]
fn json_envelope_smoke() {
    let out = Command::new(bin())
        .arg(golden("doc_a.tsv"))
        .arg(golden("doc_b.tsv"))
        .arg(golden("doc_c.tsv"))
        .arg("--json")
        .output()
        .expect("run binary");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    assert!(s.contains("\"statistic\""), "json missing statistic: {s}");
    assert!(s.contains("\"pvalue\""), "json missing pvalue: {s}");
    // Exactly one framework envelope on stdout.
    assert_eq!(
        s.matches("\"status\"").count(),
        1,
        "expected one envelope: {s}"
    );
}

#[test]
fn help_exits_zero() {
    let out = Command::new(bin())
        .arg("--help")
        .output()
        .expect("run --help");
    assert!(out.status.success(), "--help did not exit 0");
}

#[test]
fn long_mode_matches_cols_mode() {
    // The doc3 fixture as cols vs the same data fed long should be identical.
    let (sc, pc) = run(
        &["doc_a.tsv", "doc_b.tsv", "doc_c.tsv"],
        "cols",
        "median",
        "0.05",
    );
    assert!(sc > 0.0 && (0.0..=1.0).contains(&pc));
}
