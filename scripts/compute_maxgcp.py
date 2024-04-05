import argparse
import pathlib

import maxgcp  # type: ignore
import pandas as pd  # type: ignore


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gcov", type=pathlib.Path)
    parser.add_argument("--pcov", type=pathlib.Path)
    parser.add_argument("--pheno", type=pathlib.Path, required=False)
    parser.add_argument("--output-coef", type=pathlib.Path)
    parser.add_argument("--output-h2", type=pathlib.Path, required=False)
    parser.add_argument("--output-pheno", type=pathlib.Path, required=False)
    args = parser.parse_args()

    gcov_df = pd.read_csv(args.gcov, sep="\t", index_col=0)
    assert list(gcov_df.index) == list(gcov_df.columns), (
        gcov_df.index,
        gcov_df.columns,
    )
    feature_names = gcov_df.index

    pcov_df = pd.read_csv(args.pcov, sep="\t", index_col=0)
    assert list(pcov_df.index) == list(pcov_df.columns), (
        pcov_df.index,
        pcov_df.columns,
    )
    pcov_df = pcov_df.loc[feature_names, feature_names]  # type: ignore

    gcov = gcov_df.values
    pcov = pcov_df.values

    coef = maxgcp.fit_coheritability(gcov, pcov)
    coef_df = pd.DataFrame(coef, index=feature_names, columns=feature_names)
    coef_df.to_csv(args.output_coef, sep="\t", index=True)

    if args.output_h2 is not None:
        idx_to_name = {i: name for i, name in enumerate(feature_names)}
        summary_df = (
            maxgcp.summary_metrics(coef, gcov, pcov)
            .assign(phenotype=lambda df: df["phenotype_idx"].map(idx_to_name))
            .loc[:, ["phenotype", "h2", "rg", "coher", "loss"]]
        )
        summary_df.to_csv(args.output_h2, sep="\t", index=False)

    if args.pheno is not None:
        phenotypes_df = pd.read_csv(args.pheno, sep="\t", index_col=[0, 1]).loc[
            :, feature_names  # type: ignore
        ]
        maxgcp_phenotypes_df = phenotypes_df @ coef_df
        maxgcp_phenotypes_df.to_csv(args.output_pheno, sep="\t", index=True)


if __name__ == "__main__":
    main()
