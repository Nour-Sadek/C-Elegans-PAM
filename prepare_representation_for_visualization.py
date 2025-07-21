import pandas as pd
import requests
import json
from Bio import Cluster

REPRESENTATION_PATH = "./caenorhabditis_only_3UTR_representation_for_caenorhabditis_elegans.json"
REPRESENTATION_NAME = "caenorhabditis_only_3UTR"
ENSEMBL_REST = "https://rest.ensembl.org"
SOURCE_SPECIES = "caenorhabditis_elegans"
GENE_INFO_DICT_AVAILABLE = False
GENE_INFO_PATH = "./gene_ids_names_descriptions.json"


def get_gene_descriptions_names(gene_ids: list[str]) -> dict[str: list[str]]:
    """Return a dictionary where there are three different key-value pairs which are parallel lists where the values at
    each position relate to each other; it is in this format to be able to easily be loaded as a pandas dataframe.

    first key is "gene_id" where its value is a list of WormBase gene ids, where those are taken from the <gene_ids>
    list argument
    second key is "gene_name" and third key is "description" where their values are lists whose contents are filled by
    fetch requests using Ensembl's API through the URL: <ENSEMBL_REST>/lookup/id/<gene_id>?expand=1;species=<SOURCE_SPECIES>.
    If that url doesn't contain corresponding values for a gene_id of either, they are replaced with an empty string."""

    gene_info = {"gene_id": [], "gene_name": [], "description": []}
    for gene_id in gene_ids:
        # Get gene info for source species
        url = f"{ENSEMBL_REST}/lookup/id/{gene_id}?expand=1;species={SOURCE_SPECIES}"
        r = requests.get(url, headers={"content-Type": "application/json"})
        if r.ok:
            # Get the common gene name, if it is available, else the canonical transcript name
            request = r.json()
            if "display_name" in request:
                gene_name = request["display_name"]
            else:
                if "canonical_transcript" in request:
                    gene_name = request["canonical_transcript"]
                else:
                    gene_name = ""
            # Get the gene description, if it is available
            if "description" in request:
                description = request["description"]
            else:
                description = ""
            # Add the information to the dictionary
            gene_info["gene_id"].append(gene_id)
            gene_info["gene_name"].append(gene_name)
            gene_info["description"].append(description.split("[")[0].strip())
    return gene_info


if __name__ == "__main__":

    # ==== Fetch the given name and description for each WormBase GeneID ====

    if not GENE_INFO_DICT_AVAILABLE:

        # Get all the gene ids
        gene_ids = pd.read_excel("celegans_gene_ids.xlsx")
        gene_ids = list(gene_ids["gene_id"])

        # Determine the given name and description for each WormBase gene id
        gene_info = get_gene_descriptions_names(gene_ids)

        # Save the dictionary as a json file
        with open(GENE_INFO_PATH, "w") as file:
            json.dump(gene_info, file, indent=4)

    else:

        # load the gene ids, names, and descriptions dictionary
        with open(GENE_INFO_PATH, "r") as file:
            gene_info = json.load(file)

    # ==== Load in the representation json file and prepare it for loading into the Cluster 3 software ====
    
    with open(REPRESENTATION_PATH, "r") as file:
        representation = json.load(file)

    representation_df = pd.DataFrame.from_dict(representation, orient='index')

    # median center the columns (motifs)
    representation_df = representation_df - representation_df.median()
    
    # Add the gene name and description to each gene id
    # Merge <representation_df> with <gene_info>'s gene_id
    representation_df_reset = representation_df.reset_index().rename(columns={"index": "gene_id"})
    merged_df = pd.merge(representation_df_reset, gene_info, on="gene_id", how="left")

    # Combine the <gene_id>, <gene_name>, and <description> columns into one where the values are separated by |
    merged_df["combined_id"] = merged_df["gene_id"] + " | " + merged_df["gene_name"] + " | " + merged_df["description"]

    # Drop the <gene_id>, <gene_name>, and <description> columns and assign <combined_id> as the new index
    final_df = merged_df.drop(columns=["gene_id", "gene_name", "description"])
    final_df = final_df.set_index("combined_id")
    final_df.index.name = "FILE"

    # Save as a .txt file
    final_df.to_csv(f"{REPRESENTATION_NAME}_representation.txt", sep="\t")
    
    # ==== Use Cluster 3 to cluster both the genes (rows) and motifs (columns) ====
    
    handle = open(f"{REPRESENTATION_NAME}_representation.txt")
    record = Cluster.read(handle)

    # Un-centered correlation on genes and arrays, average linkage
    gene_tree = record.treecluster(transpose=False, method="a", dist="u")  # for genes
    gene_tree.scale()
    motif_tree = record.treecluster(transpose=True, method="a", dist="u")  # for motifs
    motif_tree.scale()

    # Save to input into Java TreeView
    record.save(f"{REPRESENTATION_NAME}_representation", gene_tree, motif_tree)
