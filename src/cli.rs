use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta, run};

use rsomics_fligner::{Center, FlignerResult, fligner, parse_column, parse_long};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CenterArg {
    Median,
    Mean,
    Trimmed,
}

impl From<CenterArg> for Center {
    fn from(c: CenterArg) -> Self {
        match c {
            CenterArg::Median => Center::Median,
            CenterArg::Mean => Center::Mean,
            CenterArg::Trimmed => Center::Trimmed,
        }
    }
}

/// Fligner-Killeen test for equal variances — value-exact `scipy.stats.fligner`.
///
/// The default `--center median` is the distribution-free variant, robust to
/// non-normality; `--center mean` is the original Fligner test; `--center
/// trimmed` (with `--proportiontocut`) suits heavy-tailed data. Pass two or more
/// single-column files (one per group), or one `--long` file of `value<TAB>group`
/// rows. Output is a single line `statistic<TAB>p`.
#[derive(Parser, Debug)]
#[command(name = "rsomics-fligner", version, about, long_about = None)]
pub struct Cli {
    /// Group files (≥2 single-column), or one `value<TAB>group` file with `--long`.
    #[arg(value_name = "DATA", required = true)]
    pub data: Vec<PathBuf>,

    /// Treat the single input file as long-format `value<TAB>group` rows.
    #[arg(long)]
    pub long: bool,

    /// Statistic used to center each group before the deviation transform.
    #[arg(long, value_enum, default_value_t = CenterArg::Median)]
    pub center: CenterArg,

    /// Fraction trimmed from each end when `--center trimmed`.
    #[arg(long, default_value_t = 0.05, value_name = "F")]
    pub proportiontocut: f64,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn run(self) -> ExitCode {
        let common = self.common.clone();
        run(&common, META, || {
            let groups = self.read_groups()?;
            let res: FlignerResult = fligner(&groups, self.center.into(), self.proportiontocut)?;
            if !common.json {
                println!("{}\t{}", res.statistic, res.pvalue);
            }
            Ok(res)
        })
    }

    fn read_groups(&self) -> Result<Vec<Vec<f64>>> {
        if self.long {
            if self.data.len() != 1 {
                return Err(RsomicsError::InvalidInput(
                    "--long takes exactly one input file".into(),
                ));
            }
            let f = File::open(&self.data[0]).map_err(RsomicsError::Io)?;
            return parse_long(BufReader::new(f));
        }
        if self.data.len() < 2 {
            return Err(RsomicsError::InvalidInput(
                "need at least two group files (or use --long with one file)".into(),
            ));
        }
        self.data
            .iter()
            .map(|p| {
                let f = File::open(p).map_err(RsomicsError::Io)?;
                parse_column(BufReader::new(f))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
