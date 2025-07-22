import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
import pickle
import json
import os

# url to metazoa ensembl organisms
# https://metazoa.ensembl.org/index.html

# Constants

SOURCE_SPECIES = "caenorhabditis_elegans"
MOTIF_TYPE = "TF"  # "TF" for promoter or "RBP" for 3UTR
PWM_FILE_PATH = f"./PWMs_of_{MOTIF_TYPE}_for_{SOURCE_SPECIES}.pkl"
PROMOTERS_DIR = "./promoter_sequences"
UTRs_DIR = "./3UTR_sequences"
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
COMPLEMENT_ORDER = [3, 2, 1, 0]
PROMOTER_LENGTH = 600
UTR_LENGTH = 200
# HOMOLOGOUS_SPECIES = [
#     "ascaris_summ", "brugia_malayi", "caenorhabditis_brenneri", "caenorhabditis_briggsae", "caenorhabditis_elegans",
#     "caenorhabditis_japonica", "caenorhabditis_remanei", "loa_loa", "necator_americanus", "onchocerca_volvulus",
#     "pristionchus_pacificus", "strongyloides_ratti", "trichinella_spiralis", "trichuris_muris"
# ]  # only consider nemotoda species
# HOMOLOGOUS_SPECIES = []  # consider all metazoan species
HOMOLOGOUS_SPECIES = ["caenorhabditis_brenneri", "caenorhabditis_briggsae", "caenorhabditis_elegans",
                      "caenorhabditis_japonica", "caenorhabditis_remanei"]  # only consider caenorhabditis species


def one_hot_encode(sequence: str) -> torch.tensor:
    """Return a one-hot encoded representation of the DNA sequence <sequence> os a tensor of shape
    (1 sequence, 4 nucleotides, len(>str>)).
    Which row corresponds to which nucleotide is determined by the global constant <BASE_TO_INDEX>."""

    L = len(sequence)
    tensor = torch.zeros(4, L, dtype=torch.float32)
    for index, base in enumerate(sequence.upper()):
        tensor[BASE_TO_INDEX[base], index] = 1
    return tensor.unsqueeze(0)  # Shape: (1, 4, L) as in 1 channel, 4 rows (nucleotides), L columns


def get_pwm_reverse_complement(pwm: torch.tensor) -> torch.tensor:
    """Return the reverse complement of <pwm> where <pwm> is first reversed, then the complement is determined by the
    global constant <COMPLEMENT_ORDER> which is based on the order of the nucleotides in <BASE_TO_INDEX>.
    Both <pwm> and the returned output are of shape (4 nucleotides, length of motif)."""

    # Reverse the sequence (columns)
    rev = torch.flip(pwm, dims=[1])

    # Get the complement of the bases (rows)
    complement_indices = torch.tensor(COMPLEMENT_ORDER)
    rev_comp = rev[complement_indices]

    return rev_comp


def calculate_z_score(PAM_scores_scrambled, PAM_score_real):
    """Return the z-score between <PAM_score_real> and the distribution <PAM_scores_scrambled>.
    <PAM_scores_scrambled> is a tensor of shape (n,), and PAM_score_real is a float."""

    PAM_scores_scrambled = PAM_scores_scrambled.detach().numpy()

    # Compute mean and std of scrambled scores
    mean = np.mean(PAM_scores_scrambled)
    std = np.std(PAM_scores_scrambled)

    # Compute z-score between real score and distribution
    z_score = (PAM_score_real - mean) / std

    return z_score


def correct_seq_type(seq_type: str) -> None:
    """Raise an error if <seq_type> isn't either "promoter" or "3UTR", do nothing otherwise."""
    if seq_type not in ["promoter", "3UTR"]:
        raise ValueError(f"Invalid seq_type: {seq_type}. Expected 'promoter' or '3UTR'.")


def valid_sequence(sequence: str, seq_type) -> bool:
    """Return True if <sequence> is a valid sequence where it only contains these 4 letters {"A", "T", "G", "C"} and
    it is the required length (600bp if <seq_type> is "promoter" or 200bp if <seq_type> is "3UTR"), return False
    otherwise.

    Also, raise an error if <seq_type> isn't an allowed string value of either "promoter" or "3UTR"."""

    # If an error isn't raised, then <seq_type> is either "promoter" or "3UTR"
    correct_seq_type(seq_type)

    if seq_type == "promoter":
        seq_length = PROMOTER_LENGTH
    else:  # seq_type == "3UTR"
        seq_length = UTR_LENGTH

    return sequence is not None and set(sequence.upper()) <= {"A", "T", "G", "C"} and len(sequence) == seq_length


def scan_sequence(pwm: torch.tensor, seq_tensor: torch.tensor) -> torch.tensor:
    """Return the result of performing a convolution of <pwm> over the one-hot encoded <seq_tensor>.
    <seq_tensor> is of shape (1, 4 nucleotides, length of sequence) and <pwm> is of shape (4, length of motif).
    The result of the convolution will output a tensor of shape (1, length_of_sequence - length_of_motif + 1)."""

    conv_filter = pwm.unsqueeze(0)  # (4, length of motif) => (1, 4, length of motif)
    conv_result = F.conv1d(seq_tensor, conv_filter)
    return conv_result


def scramble_pwm(pwm: torch.tensor, n: int = 100) -> torch.tensor:
    """Return <n> stacked tensors that are scrambled versions of <pwm> row-wise. The way <pwm> is scrambled is that for
    every column, the rows (4 nucleotides) are randomly switched so that the probability of each nucleotide at each
    position varies. Of course, this is one out of many ways that <pwm> could be scrambled.

    <pwm> is a tensor of shape (4, length of motif). The output scrambled tensor is of shape (n, 4, length of motif)."""

    L = pwm.shape[1]  # length of the motif (number of positions)
    scrambled_pwms = []

    for _ in range(n):  # n is the number of randomly scrambled motifs to generate
        scrambled = np.zeros_like(pwm)
        # For each column (position), randomly permute the probabilities of the bases A, T, G, C
        for i in range(L):
            scrambled[:, i] = np.random.permutation(pwm[:, i])
        scrambled_pwms.append(scrambled)

    return np.stack(scrambled_pwms)  # shape: (n, 4, L)


def scan_scrambled(conv: nn.Conv1d, filters: torch.tensor, seq_tensor: torch.tensor) -> torch.tensor:
    """Return the result of performing a convolution using lhe layer <conv> and weights <filters> over the one-hot
    encoded <seq_tensor>.
    <seq_tensor> is of shape (1, 4, length of sequence) and <filters> are of shape <n, 4, length of motif> where n is
    the number of filters (PWMs). <filters> is the output of calling the <scramble_pwm> function.
    <conv> is a 1d convolutional layer with 4 inputs channels (4 nucleotides) and n filters of size length of motif,
    with the bias term disabled.
    The result of the convolution will output a tensor of shape (n, length_of_sequence - length_of_motif + 1)."""

    with torch.no_grad():  # Assign weights
        conv.weight.copy_(filters)
    conv_result = conv(seq_tensor)
    return conv_result


def determine_eligible_genes(sequences_dir: str, seq_type: str, min_homologous_species: int) -> tuple[int, int]:
    """Return a tuple of integers that specify the number of genes and the total number of sequences (from
    <SOURCE_SPECIES> and others) that satisfy two criteria:
    - The orthologous sequences of the gene only contain the four letters {"A", "T", "G", "C"}.
    - The orthologous sequences are of the length determined for sequence types of type <seq_type>. <seq_type> can only
    be either "promoter" or "3UTR" where that length would be 600bp for "promoter" and 200bp for "3UTR".
    - There are at least <min_homologous_species> orthologous sequences available for the gene. If the global constant
    <HOMOLOGOUS_SPECIES> is not an empty list, then only the species specified in that global variable would be
    considered in the count.

    <sequence_dir> is a string that specifies the path to the directory that contains a json file for each gene and its
    orthologous sequences.
    <seq_type> specifies that type of sequence those json files contain, either "promoter" or "3UTR" sequences."""

    eligible_genes = 0
    total_num_seqs = 0
    for file_name in os.listdir(sequences_dir):
        # Get the sequences for <SOURCE_SPECIES> and the homologous species
        path = os.path.join(sequences_dir, file_name)
        with open(path, "r") as file:
            sequences_per_species = json.load(file)

        if len(HOMOLOGOUS_SPECIES) != 0:
            # Only keep the species that are part of <HOMOLOGOUS_SPECIES>
            species_of_interest = list(set(sequences_per_species.keys()) & set(HOMOLOGOUS_SPECIES))
        else:
            species_of_interest = sequences_per_species.keys()

        # Only keep those species whose corresponding sequence is a valid one
        # (only contain A, T, G, C and are of valid length)
        species_of_interest_keep = [species for species in species_of_interest if
                                    valid_sequence(sequences_per_species[species], seq_type)]

        # If the remaining homologous species is less than <min_homologous_species>, then this gene is not eligible
        if len(species_of_interest_keep) >= min_homologous_species:
            eligible_genes = eligible_genes + 1
            total_num_seqs = total_num_seqs + len(species_of_interest_keep)
    return eligible_genes, total_num_seqs


def build_representation(pwms: dict[str, torch.tensor], seq_type: str, min_homologous_species: int) -> dict[str, dict[str, float]]:
    """Return the representation as a dictionary of the form:
    key: gene id and the value:
        dictionary of the form:
        key: motif name and the value: Phylogenetically Averaged Motif (PAM) score.

    <pwms> is a dictionary where the key is a string representing the motif name and the value is a tensor of shape
    (4, length of motif). This dictionary can be created by running the <collect_motif_data.py> script which will save
    it as a pickle file and can later be loaded in as a dictionary.

    <seq_type> is a string that can take the values "promoter" or "3UTR", which specifies the directory to access the
    sequences from. These two directories can be created by running the <fetching_promoters_and_3UTRs.py> script.

    <min_homologous_species> is a positive integer that specifies the minimum number of orthologous sequences (sequences
    from homologous species plus the sequence from <SOURCE_SPECIES>) that is needed for the gene to be included in the
    final representation. A global constant <HOMOLOGOUS_SPECIES>, if it is not an empty list, can also be incorporated
    where a minimum of <min_homologous_species> sequences from the species specified in <HOMOLOGOUS_SPECIES> need to be
    included for the gene to be included in the final representation.

    The function goes through every gene file in the directory specified by <seq_type> and, if the gene specifies the
    conditions for inclusion, will further go through every PWM in <pwms> and calculate the z-score between the real
    PAM score and the scrambled PAM scores for that gene for every pwm and store those values in a dictionary which
    will be the final return output of the function after every json file has been processed."""

    # Raise Value Error if incorrect <seq_type> is provided
    correct_seq_type(seq_type)

    if seq_type == "promoter":
        sequences_dir = PROMOTERS_DIR
    else:  # seq_type == "3UTR"
        sequences_dir = UTRs_DIR

    representation = {}

    i = 1

    for file_name in os.listdir(sequences_dir):

        # get the gene_id of <SOURCE_SPECIES> associated with this file
        # remove extension ("<SOURCE-SPECIES>_<gene_id>")
        base = file_name.rsplit('.', 1)[0]
        # remove <SOURCE_SPECIES> and underscore (+1 to skip the underscore)
        gene_id = base[len(SOURCE_SPECIES) + 1:]

        # Get the sequences for <SOURCE_SPECIES> and the homologous species
        path = os.path.join(sequences_dir, file_name)
        with open(path, "r") as file:
            sequences_per_species = json.load(file)

        if len(HOMOLOGOUS_SPECIES) != 0:
            # Only keep the species that are part of <HOMOLOGOUS_SPECIES>
            species_of_interest = list(set(sequences_per_species.keys()) & set(HOMOLOGOUS_SPECIES))
        else:
            # Use all species
            species_of_interest = sequences_per_species.keys()

        # Only keep those species whose corresponding sequence is a valid one
        # (only contain A, T, G, C and are of valid length)
        species_of_interest_keep = [species for species in species_of_interest if
                                    valid_sequence(sequences_per_species[species], seq_type)]

        # If the remaining homologous species is less than <min_homologous_species>, don't process
        if len(species_of_interest_keep) < min_homologous_species:
            print(f"{i}: {gene_id} won't be processed; minimum homologous species of {min_homologous_species} was not met.")
            i = i + 1
            continue

        # One-hot encode the sequences
        one_hot_encoded_sequences_per_species = {k: one_hot_encode(v) for k, v in sequences_per_species.items() if
                                                 k in species_of_interest_keep}

        # To accumulate the z-scores for each motif
        curr_gene_z_scores = {}

        for motif_name, pwm in pwms.items():

            pwm = torch.from_numpy(pwm).float()

            L = pwm.shape[1]

            # Get the best match score for this <pwm> over each orthologous sequences of this <gene_id>
            max_conv_scores_per_seq = []
            if seq_type == "promoter":
                rev_pwm = get_pwm_reverse_complement(pwm)
            for one_hot_encoded_sequence in one_hot_encoded_sequences_per_species.values():
                fwd_conv_result = scan_sequence(pwm, one_hot_encoded_sequence)
                if seq_type == "promoter":
                    rev_conv_result = scan_sequence(rev_pwm, one_hot_encoded_sequence)
                    # Get the maximum of the rev/fwd at each position
                    best_conv_of_both = torch.max(fwd_conv_result, rev_conv_result)
                    # Get the max score of this convolution (max pooling) after considering both the fwd and rev PWM
                    max_conv_scores_per_seq.append(torch.max(best_conv_of_both).item())
                else:  # seq_type == "3UTR"
                    # Get the max score of this convolution after considering the fwd PWM only
                    max_conv_scores_per_seq.append(torch.max(fwd_conv_result).item())
            PAM_score_real = np.mean(max_conv_scores_per_seq)

            # Normalize this phylogenetically averaged score
            n = 100
            fwd_filters = torch.tensor(scramble_pwm(pwm, n))  # shape: [n, 4, L]
            conv = nn.Conv1d(in_channels=4, out_channels=n, kernel_size=L, bias=False)
            if seq_type == "promoter":
                rev_filters = torch.tensor(scramble_pwm(get_pwm_reverse_complement(pwm), n))
            max_PWMs_scores_per_seq = []
            for one_hot_encoded_sequence in one_hot_encoded_sequences_per_species.values():
                fwd_conv_results = scan_scrambled(conv, fwd_filters,
                                                  one_hot_encoded_sequence)  # shape (1, n, len_after_conv)
                if seq_type == "promoter":
                    rev_conv_results = scan_scrambled(conv, rev_filters,
                                                      one_hot_encoded_sequence)  # shape (1, n, len_after_conv)
                    # Get the maximum of the rev/fwd at each position
                    best_convs_of_both = torch.max(fwd_conv_results, rev_conv_results)  # shape (1, n, len_after_conv)
                    # Get the max match score for each PWM
                    max_PWMs_scores_per_seq.append(torch.max(best_convs_of_both, dim=2).values)  # shape (1, n)
                else:  # seq_type == "3UTR"
                    max_PWMs_scores_per_seq.append(torch.max(fwd_conv_results, dim=2).values)  # shape (1, n)
            PAM_scores_scrambled = torch.stack(max_PWMs_scores_per_seq).mean(dim=0).squeeze(0)  # shape (n,)

            # Determine the z-score between the real match score and the distribution of the scrambled match scores
            phylogenetically_averaged_best_match_z_score = calculate_z_score(PAM_scores_scrambled, PAM_score_real)

            # Save the final z-score for this <pwm> and this <gene_id>
            curr_gene_z_scores[motif_name] = phylogenetically_averaged_best_match_z_score.item()

        # Save the z-scores of all pwms for this <gene_id>
        representation[gene_id] = curr_gene_z_scores
        print(f"PAM scores for {i}: {gene_id} have been processed")
        i = i + 1

    return representation


if __name__ == "__main__":

    # Load the PWM dictionary
    with open(PWM_FILE_PATH, "rb") as file:
        pwms = pickle.load(file)

    seq_type = "promoter"
    # Build the representation (heatmap of genes and motifs but returned as a dictionary)
    representation = build_representation(pwms, seq_type, 5)

    # Save the representation as a json file
    # file_name = f"{seq_type}_representation_for_{SOURCE_SPECIES}.json"
    file_name = f"{seq_type}_representation_using_caenorhabditis_species_only.json"
    with open(file_name, "w") as file:
        json.dump(representation, file, indent=4)

    """
    # <HOMOLOGOUS_SPECIES> is an empty list
    eligible_genes, total_num_genes = determine_eligible_genes(PROMOTERS_DIR, "promoter", 10)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total promoter sequences")
    # There are 5596 genes, with 322077 total promoter sequences (min_homologous_species = 10)
    # There are 19985 genes, with 363477 total promoter sequences (min_homologous_species = 1)

    eligible_genes, total_num_genes = determine_eligible_genes(UTRs_DIR, "3UTR", 10)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total 3'UTR sequences")
    # There are 5705 genes, with 331660 total 3'UTR sequences (min_homologous_species = 10)
    # There are 19985 genes, with 372990 total 3'UTR sequences (min_homologous_species = 1)

    # <HOMOLOGOUS_SPECIES> includes the nematoda species
    eligible_genes, total_num_genes = determine_eligible_genes(PROMOTERS_DIR, "promoter", 5)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total promoter sequences")
    # There are 8233 genes, with 67938 total promoter sequences (min_homologous_species = 5)

    eligible_genes, total_num_genes = determine_eligible_genes(UTRs_DIR, "3UTR", 5)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total 3'UTR sequences")
    # There are 8467 genes, with 70623 total 3'UTR sequences (min_homologous_species = 5)

    # <HOMOLOGOUS_SPECIES> includes the 5 caenorhabditis species
    eligible_genes, total_num_genes = determine_eligible_genes(PROMOTERS_DIR, "promoter", 5)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total promoter sequences")
    # There are 3515 genes, with 17575 total promoter sequences

    eligible_genes, total_num_genes = determine_eligible_genes(UTRs_DIR, "3UTR", 5)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total 3'UTR sequences")
    # There are 3908 genes, with 19540 total 3'UTR sequences
    """
