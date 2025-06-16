import pandas as pd
import requests
import json
import os
import time

# Constants

ENSEMBL_REST = "https://rest.ensembl.org"
SOURCE_SPECIES = "caenorhabditis_elegans"
GENE_IDS_PATH = "./celegans_gene_ids.xlsx"
COMPARA = "metazoa"
PROMOTER_LENGTH = 600
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

    canonical_transcript_gene_info = {}

    # Getting the information from the canonical transcript
    for transcript in all_gene_info["Transcript"]:
        if transcript["is_canonical"] == 1:
            canonical_transcript_gene_info = transcript
            break
    translation_start_site = canonical_transcript_gene_info["Translation"]["start"]
    translation_end_site = canonical_transcript_gene_info["Translation"]["end"]
    transcription_start_site = canonical_transcript_gene_info["start"]
    transcription_end_site = canonical_transcript_gene_info["end"]
    strand = canonical_transcript_gene_info["strand"]
    chromosome = canonical_transcript_gene_info["seq_region_name"]

    # Saving the information in a dictionary specific to <target_species>
    filtered_gene_info = {"translation_start_site": translation_start_site,
                          "translation_end_site": translation_end_site,
                          "transcription_start_site": transcription_start_site,
                          "transcription_end_site": transcription_end_site,
                          "strand": strand, "chromosome": chromosome}

    return filtered_gene_info


def fetch_ranged_sequence(start, end, strand, chromosome, length, species, seq_type):
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


def accumulate_source_orthologs_info(gene_ids, print_progress=True):
    """For every protein-coding gene in <gene_ids> save a json file in <general_info> directory that stores enough
    information about the gene in order to be able to fetch the DNA sequences of the promoter and 3'UTR associated
    to it in the format:
    key: species name (string)
    value: dictionary that contains information about the start and end coordinates, strand, and chromosome.
    """
    num_genes = len(gene_ids)
    i = 1
    for gene_id in gene_ids:
        curr_gene_info = {}
        # Get gene info for source species
        url = f"{ENSEMBL_REST}/lookup/id/{gene_id}?expand=1;species={SOURCE_SPECIES}"
        r = requests.get(url, headers={"content-Type": "application/json"})
        if r.ok:
            curr_gene_info[SOURCE_SPECIES] = get_gene_info_from_request(r.json())
        else:
            print(f"Error fetching gene info for {gene_id} in {SOURCE_SPECIES}: {r.status_code} - {r.text}")

        # Get gene info for ortholog species
        print(f"Fetching orthologs for {gene_id} in {SOURCE_SPECIES}...")
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
                    time.sleep(0.10)
                    if r.ok:
                        curr_gene_info[target_species] = get_gene_info_from_request(r.json())
                    else:
                        print(f"Error fetching gene info for {target_gene_id} in {target_species}: "
                              f"{r.status_code} - {r.text}")
        else:
            print(f"Error fetching orthologs for {gene_id} in {SOURCE_SPECIES}: {r.status_code} - {r.text}")
        # Save the info for this gene as a json file
        with open(f"./general_info/{SOURCE_SPECIES}_{gene_id}.json", "w") as file:
            json.dump(curr_gene_info, file, indent=4)
        if print_progress:
            if i % 1 == 0:
                print(f"Information of {i} {SOURCE_SPECIES} genes and their orthologs out of {num_genes} have been "
                      f"collected.")
            i = i + 1


def get_source_orthologs_sequences(gene_ids, print_progress=True) -> None:
    """For every protein-coding gene in <gene_ids> save a json file in <promoter_sequences> or <3UTR_sequences>
    directories that stores the promoter and 3'UTR DNA sequences respectively associated to each ortholog gene in the
    format:
    key: species name (string)
    value: DNA sequence, promoter of length <PROMOTER_LENGTH> or 3'UTR of length <UTR_LENGTH> (string)

    Within this function, it reads in the files in directory <general_info> which were created after calling the
    function <accumulate_source_orthologs_info>, which store the required information for each ortholog species to
    extract the associated promoter and 3'UTR sequences.
    """
    num_genes = len(gene_ids)
    i = 1
    for gene_id in gene_ids:
        # Read the general info for this gene_id
        file_path = f"./general_info/{SOURCE_SPECIES}_{gene_id}.json"
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                curr_gene_species_info = json.load(f)
        else:
            continue
        curr_gene_promoters_info = {}
        curr_gene_UTRs_info = {}
        for species, gene_info in curr_gene_species_info.items():
            promoter_sequence = fetch_ranged_sequence(gene_info["transcription_start_site"],
                                                      gene_info["transcription_end_site"],
                                                      gene_info["strand"], gene_info["chromosome"],
                                                      PROMOTER_LENGTH, species, "promoter")
            UTR_sequence = fetch_ranged_sequence(gene_info["translation_start_site"],
                                                 gene_info["translation_end_site"],
                                                 gene_info["strand"], gene_info["chromosome"],
                                                 UTR_LENGTH, species, "3UTR")
            curr_gene_promoters_info[species] = promoter_sequence
            curr_gene_UTRs_info[species] = UTR_sequence
        # Save the promoter and 3'UTR sequences info for this gene as a json file
        with open(f"./promoter_sequences/{SOURCE_SPECIES}_{gene_id}.json", "w") as file:
            json.dump(curr_gene_promoters_info, file, indent=4)
        with open(f"./3UTR_sequences/{SOURCE_SPECIES}_{gene_id}.json", "w") as file:
            json.dump(curr_gene_UTRs_info, file, indent=4)
        if print_progress:
            if i % 1 == 0:
                print(f"Promoter and 3'UTR sequences of {i} {SOURCE_SPECIES} genes and their orthologs out of "
                      f"{num_genes} have been processed.")
            i = i + 1


if __name__ == "__main__":
    # A list of gene ids was extracted from the gtf file in R for <SOURCE_SPECIES> and is being imported here
    gene_ids = pd.read_excel(GENE_IDS_PATH)
    gene_ids = list(gene_ids["gene_id"])

    print("This program runs in a two step process.")
    print(f"Step 1 of 2: Save the required information for each gene of {SOURCE_SPECIES} "
          f"and its orthologs from Ensembl.")

    accumulate_source_orthologs_info(gene_ids, print_progress=True)

    print("Finished fetching, and general information for the gene ids were successfully saved as json objects. "
          "Now on to step 2.")
    print("Step 2 of 2: Use the saved information to determine the promoter and 3'UTR sequences of each gene and save "
          "as json objects.")

    get_source_orthologs_sequences(gene_ids, print_progress=True)

    print("Program Finished!")
