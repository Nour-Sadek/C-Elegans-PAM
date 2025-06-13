import pandas as pd
import requests
import json

# Constants

ENSEMBL_REST = "https://rest.ensembl.org"
SOURCE_SPECIES = "saccharomyces_cerevisiae"
COMPARA = "fungi"
PROMOTER_LENGTH = 500
UTR_LENGTH = 200


def get_gene_info_from_request(all_gene_info: dict) -> dict:
    """Return a dictionary that stores the required information from a gene in order to be able to request from
    Ensembl's REST API the promoter and 3'UTR sequences in relation to that gene.

    <all_gene_info> is a dictionary that is the result of calling .json() on a response object generated using Ensembl's
    REST API. The url to get the response was of the format:
        <ENSEMBL_REST>/lookup/id/<gene_id>?expand=1;species=<species>

    This is the information that is of interest that can be obtained from the response from this fetch request:
    - start and end coordinates of the transcript and coding region of the gene
    - strand, either +1 (forward) or -1 (reverse)
    - region that the gene is in (e.g. chromosome number)
    """

    filtered_gene_info = {}

    # Getting the information
    translation_start_site = all_gene_info["Transcript"][0]["Translation"]["start"]
    translation_end_site = all_gene_info["Transcript"][0]["Translation"]["end"]
    transcription_start_site = all_gene_info["start"]
    transcription_end_site = all_gene_info["end"]
    strand = all_gene_info["strand"]
    chromosome = all_gene_info["seq_region_name"]
    # Saving the information in a dictionary specific to <target_species>
    filtered_gene_info["translation_start_site"] = translation_start_site
    filtered_gene_info["translation_end_site"] = translation_end_site
    filtered_gene_info["transcription_start_site"] = transcription_start_site
    filtered_gene_info["transcription_end_site"] = transcription_end_site
    filtered_gene_info["strand"] = strand
    filtered_gene_info["chromosome"] = chromosome

    return filtered_gene_info


def fetch_ranged_sequence(start, end, strand, chromosome, length, species, seq_type) -> str:
    """Return a string of length <length> that represents the <seq_type> associated with the gene whose position in
    <chromosome> for <species> is bounded by the <start> and <end> coordinates, with strand <strand>.

    <seq_type> can be either 'promoter' or '3UTR'
    <strand> can be either +1 (forward) or -1 (reverse)
    """
    # Finding out the coordinate bounds for <seq_type> based on the <start> and <end> coordinates of the gene
    if (seq_type == "3UTR" and strand == -1) or (seq_type == "promoter" and strand == 1):
        seq_start = start - length
        seq_end = start - 1
    elif (seq_type == "3UTR" and strand == 1) or (seq_type == "promoter" and strand == -1):
        seq_start = end + 1
        seq_end = end + length

    # Request for the region
    region = f"{chromosome}:{seq_start}..{seq_end}:{strand}"
    region_url = f"{ENSEMBL_REST}/sequence/region/{species}/{region}"
    r = requests.get(region_url, headers={"Content-Type": "text/plain"})
    if not r.ok:
        return None

    # Extract only the sequence portion from the returned FASTA
    lines = r.text.splitlines()
    sequence = "".join(line for line in lines if not line.startswith(">"))
    return sequence


def accumulate_source_orthologs_info(gene_ids: list[str], print_progress=True) -> dict:
    """Return a dictionary of the format:
    key = gene_id associated to <SOURCE_SPECIES> (string)
    value = a dictionary of the format:
        key = species name (string)
        value = a dictionary, which is the return value of calling the function <get_gene_info_from_request> which
            contains enough information about the gene_id in order to be able to fetch the DNA sequences of the
            promoter and 3'UTR associated to it.
    """
    source_orthologs_info = {}
    num_genes = len(gene_ids)
    i = 1
    for gene_id in gene_ids:
        curr_gene_info = {}
        # Get gene info for source species
        url = f"{ENSEMBL_REST}/lookup/id/{gene_id}?expand=1;species={SOURCE_SPECIES}"
        r = requests.get(url, headers={"content-Type": "application/json"})
        if r.ok:
            curr_gene_info[SOURCE_SPECIES] = get_gene_info_from_request(r.json())
        # Get gene info for ortholog species
        url = f"{ENSEMBL_REST}/homology/id/{SOURCE_SPECIES}/{gene_id}?type=orthologues;compara={COMPARA};content-type=application/json"
        r = requests.get(url, headers={"content-Type": "application/json"})
        if r.ok:
            curr_orthologs_info = r.json()
            curr_orthologs_info = curr_orthologs_info["data"][0]["homologies"]
            for ortholog in curr_orthologs_info:
                if ortholog["type"] == "ortholog_one2one":
                    target_species = ortholog["target"]["species"]
                    target_gene_id = ortholog["target"]["id"]
                    url = f"{ENSEMBL_REST}/lookup/id/{target_gene_id}?expand=1;species={target_species}"
                    r = requests.get(url, headers={"content-Type": "application/json"})
                    if r.ok:
                        curr_gene_info[target_species] = get_gene_info_from_request(r.json())
        source_orthologs_info[gene_id] = curr_gene_info
        if print_progress:
            if i % 100 == 0:
                print(f"Information of {i} {SOURCE_SPECIES} genes and their orthologs out of {num_genes} have been "
                      f"collected.")
            i = i + 1
    return source_orthologs_info


def get_source_orthologs_sequences(source_orthologs_info, print_progress=True) -> dict:
    """Return a dictionary of the format:
    key = gene_id associated to <SOURCE_SPECIES> (string)
    value = a dictionary of the format:
        key = species name (string)
        value = a dictionary which contains two key-value pairs:
            promoter_sequence (string): DNA sequence of the promoter associated to the gene_id of length
                <PROMOTER_LENGTH>.
            3UTR_sequence (string): DNA sequence of the 3'UTR associated to the gene_id of length <UTR_LENGTH>.

    The parameter <source_orthologs_info> is the return value of calling the function <accumulate_source_orthologs_info>
    on a list of Ensembl protein-coding genes for <SOURCE_SPECIES>.
    """
    num_genes = len(source_orthologs_info)
    i = 1
    source_orthologs_sequences = {}
    for gene_id in source_orthologs_info:
        curr_gene_sequences_info = {}
        curr_gene_species_info = source_orthologs_info[gene_id]
        for species, gene_info in curr_gene_species_info.items():
            promoter_sequence = fetch_ranged_sequence(gene_info["transcription_start_site"],
                                                      gene_info["transcription_end_site"],
                                                      gene_info["strand"], gene_info["chromosome"],
                                                      PROMOTER_LENGTH, species, "promoter")
            UTR_sequence = fetch_ranged_sequence(gene_info["translation_start_site"],
                                                 gene_info["translation_end_site"],
                                                 gene_info["strand"], gene_info["chromosome"],
                                                 UTR_LENGTH, species, "3UTR")
            curr_gene_sequences_info[species] = {"promoter_sequence": promoter_sequence,
                                                 "3UTR_sequence": UTR_sequence}
        source_orthologs_sequences[gene_id] = curr_gene_sequences_info
        if print_progress:
            if i % 100 == 0:
                print(f"Promoter and 3'UTR sequences of {i} {SOURCE_SPECIES} genes and their orthologs out of "
                      f"{num_genes} have been processed.")
            i = i + 1
    return source_orthologs_sequences


if __name__ == "__main__":

    # A list of gene ids was extracted from the gtf file in R for <SOURCE_SPECIES> and is being imported here
    gene_ids = pd.read_excel("./sgd_gene_ids.xlsx")
    gene_ids = list(gene_ids["gene_id"])

    print("This program runs in a two step process.")
    print(f"Step 1 of 2: Create a dictionary that contains the required information for each gene of {SOURCE_SPECIES} "
          f"and its orthologs from Ensembl.")

    source_orthologs_info = accumulate_source_orthologs_info(gene_ids, print_progress=True)
    # Save the dictionary as a json file
    with open("source_orthologs_info.json", "w") as file:
        json.dump(source_orthologs_info, file, indent=4)

    print("Finished generating the dictionary, and it was successfully saved as a json object. Now on to step 2.")
    print("Step 2 of 2: Create a dictionary that uses the saved information to determine the promoter and 3'UTR "
          "sequences of each gene.")

    source_orthologs_sequences = get_source_orthologs_sequences(source_orthologs_info, print_progress=True)
    # Save the dictionary as a json file
    with open("source_orthologs_sequences.json", "w") as file:
        json.dump(source_orthologs_sequences, file, indent=4)

    print("Program Finished!")
