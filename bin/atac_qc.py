import argparse
import numpy as np

import scanpy as sc
import muon as mu
from muon import atac as ac

parser = argparse.ArgumentParser(description="Preprocessing ATAC data by muon.")
parser.add_argument("input_h5ad", help="Path to the input 10x h5ad format file.")
parser.add_argument("input_fragment", help="Path to the input 10x fragment file.")
parser.add_argument("output_h5ad", help="The path to the output 10x h5ad file")
args = parser.parse_args()

mdata = mu.read_10x_h5(args.input_h5ad)
mdata.var_names_make_unique()

atac = mdata.mod["atac"]

# Calculate general qc metrics using scanpy
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

# Rename columns
atac.obs.rename(
    columns={
        "n_genes_by_counts": "n_features_per_cell",
        "total_counts": "total_fragment_counts",
    },
    inplace=True,
)

# log-transform total counts and add as column
atac.obs["log_total_fragment_counts"] = np.log10(atac.obs["total_fragment_counts"])

# Calculate the nucleosome signal across cells
# set n=10e3*atac.n_obs for rough estimate but faster run time
ac.tl.locate_fragments(atac, args.input_fragment)
ac.tl.nucleosome_signal(atac, n=10e3 * atac.n_obs)

# Add group labels for above and below the nucleosome signal threshold
nuc_signal_threshold = 2
atac.obs["nuc_signal_filter"] = [
    "NS_FAIL" if ns > nuc_signal_threshold else "NS_PASS"
    for ns in atac.obs["nucleosome_signal"]
]

atac.write_h5ad(args.output_h5ad)
