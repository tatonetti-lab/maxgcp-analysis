import argparse
import pathlib

import maxgcp  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import statsmodels.api as sm  # type: ignore


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--g", type=pathlib.Path, required=True)
    parser.add_argument("--p", type=pathlib.Path, required=True)
    parser.add_argument("--output-coef", type=pathlib.Path, required=False)
    parser.add_argument("--output-h2", type=pathlib.Path, required=False)
    parser.add_argument("--output-pheno", type=pathlib.Path, required=False)
    args = parser.parse_args()

    g_df = pd.read_csv(args.g, sep="\t", index_col=[0, 1])
    p_df = pd.read_csv(args.p, sep="\t", index_col=[0, 1])

    X = sm.add_constant(p_df)
    assert isinstance(X, pd.DataFrame)

    # Save the optimal coefficients
    coef, _, _, _ = np.linalg.lstsq(X, g_df, rcond=None)
    coef_df = pd.DataFrame(coef, index=X.columns, columns=g_df.columns)
    coef_df = coef_df.drop(index=["const"])
    coef = coef_df.values
    coef_df.to_csv(args.output_coef, sep="\t")

    # Save h2, rg, coher, loss
    gcov = np.cov(g_df.T)
    pcov = np.cov(p_df.T)
    feature_names = p_df.columns
    idx_to_name = {i: name for i, name in enumerate(feature_names)}
    summary_df = (
        maxgcp.summary_metrics(coef, gcov, pcov)
        .assign(phenotype=lambda df: df["phenotype_idx"].map(idx_to_name))
        .loc[:, ["phenotype", "h2", "rg", "coher", "loss"]]
    )
    summary_df.to_csv(args.output_h2, sep="\t", index=False)

    # Save phenotypes
    maxgcp_phenotypes_df = p_df @ coef_df
    maxgcp_phenotypes_df.to_csv(args.output_pheno, sep="\t", index=True)


if __name__ == "__main__":
    main()
