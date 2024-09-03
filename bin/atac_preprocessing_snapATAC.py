#!/usr/bin/env python3

import os
import argparse
import snapatac2 as snap
import anndata as ad
import pickle


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="ATAC-seq preprocessing using snapATAC2"
    )
    parser.add_argument(
        "-f", "--input_fragment", required=True, help="Path to input fragment file"
    )
    parser.add_argument(
        "-a", "--input_h5ad", required=True, help="Path to input RNA h5ad file"
    )
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument(
        "-o", "--output_dir", default="snapatac2_out", help="Output path for h5ad file"
    )
    parser.add_argument("-g", "--genome", default=snap.genome.hg38, help="Genome name")
    return parser.parse_args()


def setup_output_paths(args):
    return {
        "output_h5ad": os.path.join(
            args.output_dir, f"{args.sample_name}_snapatac2.h5ad"
        ),
        "gene_mtx_h5ad": os.path.join(
            args.output_dir, f"{args.sample_name}_gene_mtx.h5ad"
        ),
        "peak_mtx_h5ad": os.path.join(
            args.output_dir, f"{args.sample_name}_peak_mtx.h5ad"
        ),
        "motifs_output": os.path.join(
            args.output_dir, f"{args.sample_name}_motif_enrichment.pkl"
        ),
    }


def process_atac_data(args, paths):
    filterFrags = 5000 if "RU581_LN" in args.sample_name else 2000
    data = snap.pp.import_data(
        args.input_fragment,
        chrom_sizes=args.genome,
        file=paths["output_h5ad"],
        sorted_by_barcode=False,
        min_num_fragments=filterFrags,
    )
    snap.metrics.tsse(data, args.genome)
    snap.pp.filter_cells(data, min_tsse=9, max_counts=1e20)
    snap.pp.add_tile_matrix(data)
    return data


def subset_multiome_cells(data, rna_h5ad_path):
    print("Reading RNA data and subsetting cells...")
    adata = ad.read_h5ad(rna_h5ad_path)
    multiome_cells = [cell for cell in adata.obs_names if cell in data.obs_names]
    data.subset(multiome_cells)


def perform_analysis(data):
    print("Performing doublet removal, dimensionality reduction and clustering...")
    snap.pp.select_features(data, n_features=25000)
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)
    snap.tl.spectral(data)
    snap.tl.umap(data)
    snap.pp.knn(data)
    snap.tl.leiden(data)


def create_gene_matrix(data, genome, output_path):
    gene_matrix = snap.pp.make_gene_matrix(data, genome)
    gene_matrix.write(output_path, compression="gzip")


def perform_peak_analysis(data, genome, peak_mtx_path, motifs_output_path):
    print("Peak calling and motif annotations...")
    snap.tl.macs3(data, groupby="leiden")
    peaks = snap.tl.merge_peaks(data.uns["macs3"], genome)
    peak_mat = snap.pp.make_peak_matrix(data, use_rep=peaks["Peaks"])
    peak_mat.write(peak_mtx_path, compression="gzip")

    marker_peaks = snap.tl.marker_regions(peak_mat, groupby="leiden", pvalue=0.01)
    motifs = snap.tl.motif_enrichment(
        motifs=snap.datasets.cis_bp(unique=True),
        regions=marker_peaks,
        genome_fasta=genome,
    )
    with open(motifs_output_path, "wb") as f:
        pickle.dump(motifs, f)


def main():
    args = parse_arguments()
    paths = setup_output_paths(args)

    data = process_atac_data(args, paths)
    subset_multiome_cells(data, args.input_h5ad)
    perform_analysis(data)
    create_gene_matrix(data, args.genome, paths["gene_mtx_h5ad"])
    perform_peak_analysis(
        data, args.genome, paths["peak_mtx_h5ad"], paths["motifs_output"]
    )

    data.close()


if __name__ == "__main__":
    main()
