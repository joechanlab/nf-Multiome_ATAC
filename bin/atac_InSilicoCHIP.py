#!/usr/bin/env python3

import argparse
import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.io import mmread
from scipy.sparse import csc_matrix, csr_matrix
from sklearn.preprocessing import MinMaxScaler


def parse_arguments():
    """
    Parse command line arguments for the in silico ChIP-seq calculation.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Calculate in silico ChIP-seq")
    parser.add_argument("sample", help="Sample name")
    parser.add_argument("atac_dir", help="ATAC directory")
    parser.add_argument("rna_h5ad", help="RNA path (after SEACells)")
    parser.add_argument("rna_dir", help="RNA directory (SEACells)")
    parser.add_argument("output_csv", help="Path to output csv")
    parser.add_argument("output_pkl", help="Path to output pkl")
    parser.add_argument(
        "--mode",
        choices=["pearson", "spearman", "spearman.include_repress"],
        default="pearson",
        help="Correlation mode",
    )
    return parser.parse_args()


def load_input_files(args):
    """
    Load input files for the analysis.

    Args:
        args (argparse.Namespace): Parsed command line arguments

    Returns:
        tuple: Loaded data (df, rna_adata, mc_assignments, bc, mat, peak_df)
    """
    motif_path = os.path.join(args.atac_dir, "PeaksOverlapMotifs.csv")
    df = pd.read_csv(motif_path, sep="\t", index_col="tf_name")
    rna_adata = sc.read_h5ad(args.rna_h5ad)
    if rna_adata.raw is None:
        obs_names = args.sample + "#" + rna_adata.obs_names
        mat = pd.DataFrame(
            rna_adata.X.todense(), index=obs_names, columns=rna_adata.var_names
        )
    else:
        obs_names = args.sample + "#" + rna_adata.raw.obs_names
        mat = pd.DataFrame(
            rna_adata.raw.X.todense(), index=obs_names, columns=rna_adata.raw.var_names
        )
    rna_adata = sc.AnnData(mat)
    rna_adata.raw = rna_adata
    sc.pp.filter_genes(rna_adata, min_cells=1)
    sc.pp.calculate_qc_metrics(rna_adata, inplace=True)
    rna_adata.obs.loc[:, "log10GenesPerUMI"] = (
        np.log10(rna_adata.obs.n_genes_by_counts)
    ).div(np.log10(rna_adata.obs.total_counts))
    rna_adata.obs["sample"] = args.sample
    rna_adata.obs["original_total_counts"] = rna_adata.obs["total_counts"]
    rna_adata.obs["log10_original_total_counts"] = np.log10(
        rna_adata.obs["original_total_counts"]
    )
    mito_genes = rna_adata.var_names.str.upper().str.startswith("MT-")
    rna_adata.obs["mito_frac"] = np.sum(rna_adata[:, mito_genes].X, axis=1) / np.sum(
        rna_adata.X, axis=1
    )
    ribo_genes = rna_adata.var_names.str.upper().str.startswith("RPS", "RPL")
    rna_adata.obs["RBP_frac"] = np.sum(rna_adata[:, ribo_genes].X, axis=1) / np.sum(
        rna_adata.X, axis=1
    )
    sc.pp.filter_genes(rna_adata, min_cells=1)
    sc.pp.normalize_total(rna_adata)
    sc.pp.log1p(rna_adata, base=2)

    mc_assignments = pd.read_csv(
        os.path.join(args.rna_dir, "mc_assignments.csv"), sep=",", index_col=0
    )["SEACell"]
    mc_assignments.index = f"{args.sample}#" + mc_assignments.index.str.replace(
        f"{args.sample}_", ""
    )
    mc_assignments = f"{args.sample}#" + mc_assignments

    bc = pd.read_csv(os.path.join(args.atac_dir, "peak_counts/cells.csv"), index_col=0)
    bc.index = bc.index.rename("single_cell")

    mat = mmread(os.path.join(args.atac_dir, "peak_counts/counts.mtx"))
    mat = csc_matrix(mat)

    peak_df = pd.read_csv(os.path.join(args.atac_dir, "peak_counts/peaks.csv"), sep=",")
    peak_names = (
        peak_df.seqnames
        + ":"
        + peak_df.start.astype(str)
        + "-"
        + peak_df.end.astype(str)
    )
    peak_df.index = peak_names

    return df, rna_adata, mc_assignments, bc, mat, peak_df


def preprocess_data(rna_adata, bc, mat, mc_assignments, peak_df):
    """
    Preprocess the loaded data for analysis.

    Args:
        rna_adata (anndata.AnnData): RNA expression data
        bc (pandas.DataFrame): Barcode information
        mat (scipy.sparse.csc_matrix): Peak count matrix
        mc_assignments (pandas.Series): Metacell assignments
        peak_df (pandas.DataFrame): Peak information

    Returns:
        tuple: Preprocessed data (rna_adata, atac_adata_peakscores, peak_names)
    """
    ind = bc.index.isin(mc_assignments.index)
    bc = bc[ind]
    mat = mat[:, ind]

    sc.pp.filter_genes(rna_adata, min_cells=1)
    sc.pp.calculate_qc_metrics(rna_adata, inplace=True)
    rna_adata.obs.loc[:, "log10GenesPerUMI"] = (
        np.log10(rna_adata.obs.n_genes_by_counts)
    ).div(np.log10(rna_adata.obs.total_counts))

    mat2 = (
        pd.DataFrame(mat.todense().T, index=bc.index)
        .groupby(mc_assignments.loc[bc.index])
        .sum()
    )
    atac_adata_peakscores = sc.AnnData(
        csr_matrix(mat2), obs=mat2.loc[:, []], var=peak_df
    )
    bc_intersect = rna_adata.obs_names.intersection(atac_adata_peakscores.obs_names)
    atac_adata_peakscores = atac_adata_peakscores[bc_intersect, :]
    rna_adata = rna_adata[bc_intersect, :]
    sc.pp.filter_genes(atac_adata_peakscores, min_cells=1)
    peak_names = atac_adata_peakscores.var_names

    rna_adata.obs.loc[:, "SEACell"] = rna_adata.obs.index.str.replace(
        ".*#", "", regex=True
    )
    atac_adata_peakscores.obs["SEACell"] = rna_adata.obs["SEACell"]

    return rna_adata, atac_adata_peakscores, peak_names


def calculate_peak_acc_tf_idf(atac_adata_peakscores):
    """
    Calculate peak accessibility TF-IDF scores.

    Args:
        atac_adata_peakscores (anndata.AnnData): ATAC peak scores

    Returns:
        numpy.ndarray: Peak accessibility TF-IDF scores
    """
    peak_acc_tf_idf_1 = (
        atac_adata_peakscores.X.toarray().T / atac_adata_peakscores.X.toarray().sum(1)
    ).T
    sums = (atac_adata_peakscores.X > 0).sum(0)
    idf = atac_adata_peakscores.shape[0] / (sums + np.random.randn(len(sums)) * 0.0001)
    idf = np.asarray(idf).ravel()
    return peak_acc_tf_idf_1 * idf


def calculate_tf_peak_correlations(TFExp, PeakAcc, mode):
    """
    Calculate correlations between TF expression and peak accessibility.

    Args:
        TFExp (numpy.ndarray): TF expression matrix
        PeakAcc (numpy.ndarray): Peak accessibility matrix
        mode (str): Correlation mode ('pearson' or 'spearman')

    Returns:
        numpy.ndarray: TF-peak correlations
    """
    if mode == "pearson":
        return 1.0 - scipy.spatial.distance.cdist(
            TFExp.T, PeakAcc.T, metric="correlation"
        )
    else:
        return 1.0 - scipy.spatial.distance.cdist(
            scipy.stats.rankdata(TFExp.T, axis=1),
            scipy.stats.rankdata(PeakAcc.T, axis=1),
            metric="correlation",
        )


def calculate_binding_scores(tf_peak_correlations, max_peak_acc, motif_scores):
    """
    Calculate binding scores for TF-peak pairs.

    Args:
        tf_peak_correlations (numpy.ndarray): TF-peak correlations
        max_peak_acc (numpy.ndarray): Maximum peak accessibility
        motif_scores (numpy.ndarray): Motif scores

    Returns:
        numpy.ndarray: Binding scores
    """
    return (
        tf_peak_correlations
        * MinMaxScaler()
        .fit_transform((max_peak_acc * motif_scores).reshape(-1, 1))
        .ravel()
    )


def main():
    """
    Main function to run the in silico ChIP-seq analysis.
    """
    # Parse command line arguments and load input files
    args = parse_arguments()
    df, rna_adata, mc_assignments, bc, mat, peak_df = load_input_files(args)

    # Preprocess RNA and ATAC data
    rna_adata, atac_adata_peakscores, peak_names = preprocess_data(
        rna_adata, bc, mat, mc_assignments, peak_df
    )

    # Calculate peak accessibility TF-IDF
    peak_acc_tf_idf = calculate_peak_acc_tf_idf(atac_adata_peakscores)

    # Filter transcription factors
    tfs = np.intersect1d(rna_adata.var_names, df.index)
    df = df.loc[tfs]

    # Extract TF expression and peak accessibility matrices
    TFExp = rna_adata[:, tfs].X.toarray()
    PeakAcc = peak_acc_tf_idf

    # Calculate TF-peak correlations
    tf_peak_correlations = calculate_tf_peak_correlations(TFExp, PeakAcc, args.mode)
    max_peak_acc = np.max(PeakAcc, axis=0)

    # Prepare peak information dataframe
    new_peaks_df = atac_adata_peakscores.var.loc[:, ["seqnames", "start", "end"]]
    new_peaks_df.loc[:, "peak_name"] = (
        new_peaks_df.seqnames
        + ":"
        + new_peaks_df.start.astype(str)
        + "-"
        + new_peaks_df.end.astype(str)
    )

    # Sort and deduplicate motif scores
    df = df.sort_values("motifscore", ascending=False).drop_duplicates(
        ["tf", "peak_name"]
    )

    # Initialize arrays for binding scores and peaks
    tf_binding_scores = np.zeros((tfs.shape[0], peak_names.shape[0]))
    tf_peaks = {}

    # Calculate binding scores for each TF
    for i, tf in enumerate(tfs):
        print(tf)
        peak_correlations = tf_peak_correlations[i, :]
        try:
            motif_scores = (
                df.loc[tf]
                .merge(
                    new_peaks_df, left_on="peak_name", right_on="peak_name", how="right"
                )["motifscore"]
                .fillna(0)
                .values
            )
        except AttributeError:
            print("\toops")
            motif_scores = (
                pd.DataFrame(df.loc[tf])
                .T.merge(
                    new_peaks_df, left_on="peak_name", right_on="peak_name", how="right"
                )["motifscore"]
                .fillna(0)
                .values
            )
        binding_scores = calculate_binding_scores(
            peak_correlations, max_peak_acc, motif_scores
        )
        tf_binding_scores[i, :] = binding_scores
        if args.mode == "spearman.include_repress":
            abs_scores = np.abs(binding_scores)
            tf_peaks[tf] = np.where(abs_scores > np.percentile(abs_scores, 95))[0]
        else:
            tf_peaks[tf] = np.where(binding_scores > np.percentile(binding_scores, 95))[
                0
            ]

    # Save results to files
    pd.DataFrame(tf_binding_scores, index=tfs, columns=new_peaks_df.peak_name).to_csv(
        args.output_csv, sep="\t"
    )

    with open(args.output_pkl, "wb") as handle:
        pickle.dump(tf_peaks, handle)


if __name__ == "__main__":
    main()
