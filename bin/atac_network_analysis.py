#!/usr/bin/env python3

import argparse
import pickle
import snapatac2 as snap
import pandas as pd
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(
        description="Get TF binding and peak lineage information using SnapATAC2"
    )
    parser.add_argument(
        "peak_mtx_h5ad", required=True, help="Path to peak matrix h5ad file"
    )
    parser.add_argument(
        "gene_mtx_h5ad", required=True, help="Path to gene matrix h5ad file"
    )
    parser.add_argument(
        "diff_peak_csv", required=True, help="Path to peak matrix csv file"
    )
    parser.add_argument("motif_pkl", required=True, help="Path to motif pkl file")
    parser.add_argument("--output", required=True, help="Path to output network file")
    parser.add_argument(
        "--genome", default=snap.genome.hg38, help="Path to genome fasta file"
    )
    args = parser.parse_args()

    # Load the AnnData object
    diff_peaks = pd.read_csv(args.diff_peak_csv)
    with open(args.motif_pkl, "rb") as f:
        motifs = pickle.load(f)

    # Merge motifs into a single DataFrame and find unique values of the id column
    merged_motifs = pd.concat(
        [df.to_pandas() for df in motifs.values()], ignore_index=True
    )
    filtered_motifs = merged_motifs[merged_motifs["adjusted p-value"] < 0.05]
    unique_motif_ids = filtered_motifs["id"].unique()
    print(f"Number of unique motif IDs after filtering: {len(unique_motif_ids)}")

    peak_mat = snap.read(args.peak_mtx_h5ad, backed=None)
    peak_mat = peak_mat[:, diff_peaks.loc[:, "feature name"].unique()]
    sc.pp.normalize_total(peak_mat)
    sc.pp.log1p(peak_mat)

    gene_mat = snap.read(args.gene_mtx_h5ad, backed=None)
    sc.pp.filter_genes(gene_mat, min_cells=5)
    sc.pp.normalize_total(gene_mat)
    sc.pp.log1p(gene_mat)

    # Obtain the actual list cis-bp motifs
    cis_bp_motifs = snap.datasets.cis_bp(unique=True)
    # Filter cis_bp_motifs based on unique_motif_ids
    filtered_cis_bp_motifs = [
        motif for motif in cis_bp_motifs if motif.id in unique_motif_ids
    ]
    motif_dict = {motif.id: motif for motif in filtered_cis_bp_motifs}
    selected_motifs = [
        motif_dict[motif_id] for motif_id in unique_motif_ids if motif_id in motif_dict
    ]
    print(f"Number of selected motifs: {len(selected_motifs)}")

    # Initialize the network from gene annotations
    regions = diff_peaks.loc[:, "feature name"].unique()
    network = snap.tl.init_network_from_annotation(
        regions=regions, anno_file=args.genome, upstream=250000, downstream=250000
    )

    # Add TF binding information
    snap.tl.add_tf_binding(
        network, motifs=selected_motifs, genome_fasta=args.genome, pvalue=1e-5
    )

    # Add correlation scores
    snap.tl.add_cor_scores(network, peak_mat=peak_mat, gene_mat=gene_mat)

    # Link TFs to genes
    genetic_network = snap.tl.link_tf_to_gene(network)

    # Prune the network
    pruned_network = snap.tl.prune_network(
        genetic_network,
        edge_filter=lambda _, __, data: abs(data.cor_score) > 0.3
        if data.cor_score is not None
        else False,
    )

    # Save the network
    with open(args.output, "wb") as f:
        pickle.dump(pruned_network, f)

    print(f"Network saved to {args.output}")


if __name__ == "__main__":
    main()
