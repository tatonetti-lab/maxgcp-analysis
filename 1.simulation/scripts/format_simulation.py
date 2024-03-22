import argparse
import pathlib

import numpy as np  # type: ignore
import pandas as pd  # type: ignore


def main():
    """
    PhenotypeSimulator scales the phenotypes to mean zero, unit variance
    as the last step. Therefore, to scale the genetic fixed effects properly,
    we have to sum up all the relevant phenotypic components, compute their
    mean and standard deviation, and then scale the genetic fixed effects
    accordingly.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--ysim", type=pathlib.Path)
    parser.add_argument("--ygfixed", type=pathlib.Path)
    parser.add_argument("--ygbg", type=pathlib.Path, required=False)
    parser.add_argument("--ynbg", type=pathlib.Path)
    parser.add_argument("--output-pheno", type=pathlib.Path)
    parser.add_argument("--output-g", type=pathlib.Path)
    args = parser.parse_args()

    y_df = pd.read_csv(args.ysim, index_col=0).pipe(format_raw_df)
    gfixed_df = pd.read_csv(args.ygfixed, index_col=0).pipe(format_raw_df)
    nbg_df = pd.read_csv(args.ynbg, index_col=0).pipe(format_raw_df)
    if args.ygbg is not None:
        gbg_df = pd.read_csv(args.ygbg, index_col=0).pipe(format_raw_df)
    else:
        gbg_df = None

    gfixed_df = scale_gfixed(y_df, gfixed_df, nbg_df, gbg_df)

    y_df.to_csv(args.output_pheno, sep="\t")
    gfixed_df.to_csv(args.output_g, sep="\t")


def format_raw_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Format the phenotype names to have the same number of digits, set the index
    """
    n_items = len(str(df.shape[1]))
    old_to_new = dict()
    for col in df.columns:
        if "Trait" in col:
            trait, num = col.split("_")
            old_to_new[col] = f"{trait}_{num.zfill(n_items)}"
        else:
            old_to_new[col] = col

    return (
        df.rename_axis(index="IID")
        .reset_index()
        .assign(FID=lambda df: df["IID"])
        .set_index(["FID", "IID"])
        .rename(columns=old_to_new)
    )


def scale_gfixed(
    y_df: pd.DataFrame,
    gfixed_df: pd.DataFrame,
    nbg_df: pd.DataFrame,
    gpg_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Scale the genetic fixed effects to the appropriate liability scale.
    """
    if gpg_df is not None:
        y_tot_df = gfixed_df + gpg_df + nbg_df
    else:
        y_tot_df = gfixed_df + nbg_df

    mean = y_tot_df.mean(axis=0)
    std = y_tot_df.std(axis=0)

    y_scaled = (y_tot_df - mean) / std
    assert np.allclose(y_scaled, y_df), "There were additional phenotypic components!"

    gfixed_df = (gfixed_df - mean) / std
    return gfixed_df


def flip_phenotypes(
    seed: int,
    y_df: pd.DataFrame,
    gfixed_df: pd.DataFrame,
    nbg_df: pd.DataFrame,
    gbg_df: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """
    Randomly flip the sign of half the phenotypes. Simulated phenotype covariances
    are all positive. This is a simple way to introduce negative covariances.
    """
    np.random.seed(seed)
    flip_mask = np.random.choice([1, -1], size=y_df.shape[1])

    y_df *= flip_mask
    gfixed_df *= flip_mask
    nbg_df *= flip_mask

    if gbg_df is not None:
        gbg_df *= flip_mask

    return y_df, gfixed_df, gbg_df


if __name__ == "__main__":
    main()
