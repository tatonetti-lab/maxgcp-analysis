import argparse
import pathlib

import polars as pl  # type: ignore
import sklearn.metrics  # type: ignore


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pheno", type=pathlib.Path, nargs="+", required=True)
    parser.add_argument("--g", type=pathlib.Path, nargs="+", required=True)
    parser.add_argument("--causal", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()

    pheno_df = read_full_gwas_df(args.pheno)
    g_df = read_full_gwas_df(args.g).drop("method")

    n_chars = len(str(len(args.pheno)))

    causal_df = (
        pl.read_csv(args.causal, separator="\t")
        .melt(id_vars="phenotype")
        .with_columns(
            phenotype=pl.lit("Trait_")
            + pl.col("phenotype").str.extract("([0-9]+)$").str.zfill(n_chars),
            ID=pl.col("variable").str.extract("([0-9:]+)$"),
            causal=pl.col("value").ne(0),
        )
        .select(["phenotype", "ID", "causal"])
    )

    summary_df = summarize(pheno_df, g_df, causal_df)
    summary_df.write_csv(args.output, separator="\t")


def read_full_gwas_df(paths: list[pathlib.Path]) -> pl.DataFrame:
    full_df = pl.DataFrame()
    for path in paths:
        df = read_gwas_df(path)
        full_df = pl.concat([full_df, df])

    return full_df


def read_gwas_df(path: pathlib.Path) -> pl.DataFrame:
    method, phenotype = (
        path.with_suffix("").with_suffix("").with_suffix("").name.rsplit(".", 1)
    )
    return pl.read_csv(path, separator="\t", columns=["ID", "P", "BETA"]).with_columns(
        phenotype=pl.lit(phenotype),
        method=pl.lit(method),
    )


def summarize(
    pheno_df: pl.DataFrame, g_df: pl.DataFrame, causal_df: pl.DataFrame
) -> pl.DataFrame:
    df = (
        pheno_df.join(g_df, on=["ID", "phenotype"], suffix="_true")
        .join(causal_df, on=["ID", "phenotype"])
        .with_columns(causal=pl.col("causal").ne(0))
        .with_columns(
            positive=pl.col("P_true").lt(0.05) & pl.col("causal"),
            negative=pl.col("P_true").gt(0.05) & ~pl.col("causal"),
        )
        .filter(pl.col("positive") | pl.col("negative"))
    )

    if df["positive"].sum() == 0 or df["negative"].sum() == 0:
        auroc = lambda _: -1
    else:
        auroc = lambda list_of_series: sklearn.metrics.roc_auc_score(
            list_of_series[0], list_of_series[1]
        )

    return (
        df.with_columns(
            true_positive=pl.col("P").lt(0.05) & pl.col("positive"),
            true_negative=pl.col("P").gt(0.05) & pl.col("negative"),
            false_positive=pl.col("P").lt(0.05) & pl.col("negative"),
            false_negative=pl.col("P").gt(0.05) & pl.col("positive"),
        )
        .group_by(["phenotype", "method"])
        .agg(
            P=pl.sum("positive"),
            N=pl.sum("negative"),
            TP=pl.sum("true_positive"),
            TN=pl.sum("true_negative"),
            FP=pl.sum("false_positive"),
            FN=pl.sum("false_negative"),
            total=pl.len(),
            bias=(pl.col("BETA") - pl.col("BETA_true")).mean(),
            rmse=(pl.col("BETA") - pl.col("BETA_true")).pow(2).mean().sqrt(),
            AUROC=pl.map_groups(
                exprs=["positive", -pl.col("P").log10()],
                function=auroc,
            ),
        )
        .with_columns(
            sensitivity=pl.col("TP") / pl.col("P"),
            specificity=pl.col("TN") / pl.col("N"),
        )
    )


if __name__ == "__main__":
    main()
