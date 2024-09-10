#!/usr/bin/env python3

import argparse
import os
import snapatac2 as snap
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create joint embedding from ATAC and RNA h5ad files using snapATAC2"
    )
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("atac_h5ad", help="Path to ATAC h5ad file")
    parser.add_argument("rna_h5ad", help="Path to RNA h5ad file")
    parser.add_argument(
        "--n_comps",
        type=int,
        default=30,
        help="Number of components for dimensionality reduction",
    )
    return parser.parse_args()


def setup_output_paths(args):
    return {
        "output_rna_h5ad": os.path.join(
            args.sample_name, f"{args.sample_name}_rna.h5ad"
        ),
        "output_atac_h5ad": os.path.join(
            args.sample_name, f"{args.sample_name}_atac.h5ad"
        ),
        "output_joint_embedding": os.path.join(
            args.sample_name, f"{args.sample_name}_joint_embedding.npy"
        ),
    }


def load_data(atac_path, rna_path):
    print("Loading ATAC and RNA data...")
    atac_data = snap.read(atac_path, backed=None)
    rna_data = snap.read(rna_path, backed=None)
    rna_data.X = rna_data.X.tocsr()
    intersect_cells = atac_data.obs_names.intersection(rna_data.obs_names)
    atac_data = atac_data[intersect_cells]
    rna_data = rna_data[intersect_cells]
    return atac_data, rna_data


def create_joint_embedding(atac_data, rna_data, n_comps):
    print("Creating joint embedding...")
    joint_embedding = snap.tl.multi_spectral(
        [atac_data, rna_data], n_comps=n_comps, features=None
    )[1]
    return joint_embedding


def save_joint_embedding(joint_embedding, atac_data, rna_data, paths):
    print("Saving joint embedding...")
    np.save(paths["output_joint_embedding"], joint_embedding)
    atac_data.write(paths["output_atac_h5ad"])
    rna_data.write(paths["output_rna_h5ad"])


def main():
    args = parse_arguments()
    os.makedirs(args.sample_name, exist_ok=True)
    paths = setup_output_paths(args)
    atac_data, rna_data = load_data(args.atac_h5ad, args.rna_h5ad)
    joint_embedding = create_joint_embedding(atac_data, rna_data, args.n_comps)
    save_joint_embedding(joint_embedding, atac_data, rna_data, paths)


if __name__ == "__main__":
    main()
