import argparse
import os
import shutil
import pickle
import subprocess
import shlex
from collections import defaultdict

import scanpy as sc
import muon
import pandas as pd
import numpy as np
from tqdm import tqdm
import synapseclient

import json


def get_promoter_name(row_info):
    # Determine the promoter name based on type
    if row_info["type"] == "non-targeting":
        return "non-targeting"
    else:
        promoter_name = f"{row_info['intended_target_name']}|"
        promoter_name += (
            f"{row_info['intended_target_chr']}:"
            f"{row_info['intended_target_start']:.0f}-"
            f"{row_info['intended_target_end']:.0f}"
        )
        return promoter_name


def infer_data_label(args):
    """
    Decide the output file stem.

    Synapse mode:
        --synapse_id syn70753570
        -> syn70753570.h5mu

    gcloud mode:
        --synapse_id syn70753570 --gcloud_uri gs://bucket/file.h5mu
        -> syn70753570.h5mu

        --gcloud_uri gs://bucket/path/file.h5mu
        -> file.h5mu
    """
    if args.synapse_id:
        return args.synapse_id

    basename = os.path.basename(args.gcloud_uri.rstrip("/"))
    stem, _ = os.path.splitext(basename)

    if not stem:
        raise ValueError("Could not infer output file name from --gcloud_uri.")

    return stem


def download_from_gcloud(gcloud_uri, final_path_name):
    """
    Download a file from Google Cloud Storage by executing the gcloud command.

    Example command:
        gcloud --quiet storage cp gs://bucket/path/file.h5mu ./data/file.h5mu
    """
    if not gcloud_uri.startswith("gs://"):
        raise ValueError("--gcloud_uri must start with 'gs://'.")

    cmd = [
        "gcloud",
        "--quiet",
        "storage",
        "cp",
        gcloud_uri,
        final_path_name,
    ]

    print(f"Downloading from Google Cloud Storage...")
    print(f"Running command: {shlex.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as e:
        raise RuntimeError(
            "gcloud command was not found. Please install Google Cloud CLI "
            "and make sure 'gcloud' is available in PATH."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"gcloud download failed with exit code {e.returncode}."
        ) from e


def download_from_synapse(synapse_id, auth_token, final_path_name):
    """
    Download a file from Synapse.
    """
    print("Logging in to Synapse...")
    syn = synapseclient.Synapse()
    syn.login(authToken=auth_token)

    print(f"Downloading {synapse_id} from Synapse...")
    entity = syn.get(synapse_id)

    print("Moving downloaded file to target directory...")
    if os.path.abspath(entity.path) != os.path.abspath(final_path_name):
        shutil.move(entity.path, final_path_name)


def main():
    parser = argparse.ArgumentParser(
        description="Process single-cell RNA and sgRNA data from Synapse or Google Cloud Storage."
    )

    parser.add_argument(
        "--synapse_id",
        required=False,
        help="Synapse ID to download, e.g. syn70753570. In gcloud mode, this is used as the output file stem if provided.",
    )
    parser.add_argument(
        "--auth_token",
        required=False,
        help="Synapse personal access token for login. Required only when not using --gcloud.",
    )
    parser.add_argument(
        "--out_dir",
        default="./data",
        help="Output directory for processed files.",
    )
    parser.add_argument(
        "--gcloud",
        action="store_true",
        help="Use gcloud command to download data instead of Synapse.",
    )
    parser.add_argument(
        "--gcloud_uri",
        required=False,
        help="Google Cloud Storage URI to download, e.g. gs://bucket/path/file.h5mu. Required when --gcloud is specified.",
    )

    args = parser.parse_args()

    if args.gcloud:
        if not args.gcloud_uri:
            parser.error("--gcloud_uri is required when --gcloud is specified.")
    else:
        if not args.synapse_id:
            parser.error("--synapse_id is required when --gcloud is not specified.")
        if not args.auth_token:
            parser.error("--auth_token is required when --gcloud is not specified.")

    out_dir = args.out_dir

    # Prepare data folder
    os.makedirs(out_dir, exist_ok=True)

    data_label = infer_data_label(args)

    config1_2_path = os.path.join(out_dir, "config1_2.json")
    config3_path = os.path.join(out_dir, "config3.json")

    # Define output file paths
    final_path_name = os.path.join(out_dir, f"{data_label}.h5mu")
    dict_file = os.path.join(out_dir, "gRNA_dict.pickle")
    pca_file = os.path.join(out_dir, "pca_dataframe.pickle")
    annotation_file_path = os.path.join(out_dir, "annotation_table.csv")

    # Download mudata
    if args.gcloud:
        download_from_gcloud(args.gcloud_uri, final_path_name)
    else:
        download_from_synapse(args.synapse_id, args.auth_token, final_path_name)

    # Load mudata
    print("Loading MuData...")
    adata_mu = muon.read_h5mu(final_path_name)
    print(adata_mu)

    ### Process gene expression
    print("Processing gene expression...")
    adata_exp = adata_mu["gene"].copy()

    # only consider genes with more than 1 count
    sc.pp.filter_genes(adata_exp, min_counts=1)
    sc.pp.normalize_total(adata_exp)

    # Normalize and scale the matrix
    sc.pp.log1p(adata_exp)
    sc.pp.scale(adata_exp)

    # PCA
    print("Running PCA...")
    sc.tl.pca(adata_exp, random_state=0, n_comps=50)

    X = pd.DataFrame(
        adata_exp.obsm["X_pca"].copy(),
        index=adata_exp.obs.index,
    )

    print(f"Saving PCA data to '{pca_file}'...")
    X.to_pickle(pca_file)

    ### Process gRNA
    print("Processing gRNA...")
    adata_guide = adata_mu["guide"].copy()

    gRNA_info_df = adata_guide.var.copy()
    promoter_intended_names = gRNA_info_df.apply(
        lambda x: get_promoter_name(x),
        axis=1,
    )
    gRNA_info_df["intended_target_promoter"] = promoter_intended_names

    print(f"Saving annotation table to '{annotation_file_path}'...")
    gRNA_info_df.to_csv(annotation_file_path)

    # Modify gRNA dict
    print("Extracting non-zero counts for gRNA dictionary...")

    if hasattr(adata_guide.X, "nonzero"):
        non_zero_rows, non_zero_cols = adata_guide.X.nonzero()
    else:
        non_zero_rows, non_zero_cols = np.where(adata_guide.X != 0)

    if len(non_zero_rows) == 0:
        print("Warning: No non-zero counts found in the sgRNA data. Returning an empty dictionary.")
        gRNA_dict = {}
    else:
        gRNA_name_list = adata_guide.var.index.to_list()
        cell_name_list = adata_guide.obs.index.to_list()

        print("Creating gRNA dictionary...")
        gRNA_dict_default = defaultdict(list)
        num_pairs = len(non_zero_rows)

        try:
            for i in tqdm(range(num_pairs), desc="Processing gRNA pairs"):
                row_idx = non_zero_rows[i]
                col_idx = non_zero_cols[i]

                cell_name = cell_name_list[row_idx]
                gRNA_name = gRNA_name_list[col_idx]

                gRNA_dict_default[gRNA_name].append(cell_name)

            gRNA_dict = dict(gRNA_dict_default)
            print(f"gRNA dictionary creation complete. Found {len(gRNA_dict)} types of gRNAs.")

        except Exception as e:
            print(f"An error occurred during gRNA dictionary creation loop: {e}")
            gRNA_dict = {}

    # Save dictionary
    print(f"Saving gRNA dictionary to '{dict_file}'...")

    try:
        with open(dict_file, mode="wb") as fo:
            pickle.dump(gRNA_dict, fo)
        print("gRNA dictionary saved.")

    except IOError as e:
        print(f"Error: Failed to write gRNA dictionary file '{dict_file}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred while saving the gRNA dictionary: {e}")

    # Prepare config file for pipeline
    with open("./energy_dist_pipeline/config.json", "r") as f:
        step12_config = json.load(f)

    with open("./energy_dist_pipeline/config_clustering.json", "r") as f:
        step3_config = json.load(f)

    # Modify config file for this pipeline
    step12_config["output_file_name_list"]["OUTPUT_FOLDER"] = out_dir
    step12_config["output_file_name_list"]["pca_table"] = "pca_dataframe.pickle"
    step12_config["output_file_name_list"]["gRNA_dict"] = "gRNA_dict.pickle"
    step12_config["output_file_name_list"]["OVERWRITE_PCA_DICT"] = False

    step12_config["input_data"]["annotation_file"]["file_path"] = annotation_file_path
    step12_config["input_data"]["annotation_file"]["concatenate_key"] = "intended_target_promoter"

    step12_config["gRNA_filtering"]["perform_targeting_filtering"] = True
    step12_config["gRNA_filtering"]["perform_nontargeting_filtering"] = True

    # As this pipeline doesn't use h5ad or sgRNA dataframe, fill with dummy
    step12_config["input_data"]["h5ad_file"]["file_path"] = "dummy"
    step12_config["input_data"]["sgRNA_file"]["file_path"] = "dummy"

    with open(config1_2_path, "w") as f:
        json.dump(step12_config, f, indent=2)

    with open(config3_path, "w") as f:
        json.dump(step3_config, f, indent=2)

    print("Processing complete!")


if __name__ == "__main__":
    main()