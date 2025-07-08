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
REVERSE_ORDER = [3, 2, 1, 0]
PROMOTER_LENGTH = 600
UTR_LENGTH = 200


def one_hot_encode(sequence: str):
    L = len(sequence)
    tensor = torch.zeros(4, L, dtype=torch.float32)
    for index, base in enumerate(sequence.upper()):
        tensor[BASE_TO_INDEX[base], index] = 1
    return tensor.unsqueeze(0)  # Shape: (1, 4, L) as in 1 channel, 4 rows (nucleotides), L columns


def get_pwm_reverse_complement(pwm):
    # Reverse the sequence (columns)
    rev = torch.flip(pwm, dims=[1])

    # Get the complement of the bases (rows)
    complement_indices = torch.tensor(REVERSE_ORDER)
    rev_comp = rev[complement_indices]

    return rev_comp


def pwm_to_filter(pwm):
    # PWM shape: [4, length of motif] => Convolutional filter shape: [1, 4, length of motif]
    return pwm.unsqueeze(0)


def calculate_z_score(PAM_scores_scrambled, PAM_score_real):
    """<PAM_scores_scrambled> is a tensor of shape (n,) and PAM_score_real is a float."""
    PAM_scores_scrambled = PAM_scores_scrambled.detach().numpy()

    # Compute mean and std of scrambled scores
    mean = np.mean(PAM_scores_scrambled)
    std = np.std(PAM_scores_scrambled)

    # Compute z-score between real score and distribution
    z_score = (PAM_score_real - mean) / std

    return z_score


def correct_seq_type(seq_type):
    if seq_type not in ["promoter", "3UTR"]:
        raise ValueError(f"Invalid seq_type: {seq_type}. Expected 'promoter' or '3UTR'.")


def valid_sequence(sequence: str, seq_type) -> bool:
    """seq_type can be either 'promoter' or '3UTR'"""
    correct_seq_type(seq_type)

    if seq_type == "promoter":
        seq_length = PROMOTER_LENGTH
    else:  # seq_type == "3UTR"
        seq_length = UTR_LENGTH

    return sequence is not None and set(sequence.upper()) <= {"A", "T", "G", "C"} and len(sequence) == seq_length


def scan_sequence(pwm, seq_tensor):
    conv_filter = pwm_to_filter(pwm)
    conv_result = F.conv1d(seq_tensor, conv_filter)
    # max_score = torch.max(conv_result).item()
    return conv_result


def scramble_pwm(pwm, n=100):
    L = pwm.shape[1]  # length of the motif (number of positions)
    scrambled_pwms = []

    for _ in range(n):  # n is the number of randomly scrambled motifs to generate
        scrambled = np.zeros_like(pwm)
        # For each column (position), randomly permute the probabilities of the bases A, T, G, C
        for i in range(L):
            scrambled[:, i] = np.random.permutation(pwm[:, i])
        scrambled_pwms.append(scrambled)

    return np.stack(scrambled_pwms)  # shape: (n, 4, L)


def scan_scrambled(conv, filters, seq_tensor):
    with torch.no_grad():  # Assign weights
        conv.weight.copy_(filters)
    conv_result = conv(seq_tensor)
    # Max pool over the sequence positions
    # max_scores, _ = torch.max(conv_result, dim=2)  # shape: (1, n)
    # max_scores = max_scores.squeeze(0)  # flatten to shape: (n)

    # return max_scores.tolist()
    return conv_result


def determine_eligible_genes(sequences_dir, seq_type, min_orthologous_species):
    eligible_genes = 0
    total_num_genes = 0
    for file_name in os.listdir(sequences_dir):
        # Get the sequences for <SOURCE_SPECIES> and the orthologous species
        path = os.path.join(sequences_dir, file_name)
        with open(path, "r") as file:
            sequences_per_species = json.load(file)

        # Only keep the species that are part of <ORTHOLOGOUS_SPECIES>
        # species_of_interest = list(set(sequences_per_species.keys()) & set(ORTHOLOGOUS_SPECIES))

        # Only keep those species whose corresponding sequence is a valid one
        # (only contain A, T, G, C and are of valid length)
        species_of_interest_keep = [species for species in sequences_per_species.keys() if
                                    valid_sequence(sequences_per_species[species], seq_type)]

        # If the remaining orthologous species is less than <min_orthologous_species>, then this gene is not eligible
        if len(species_of_interest_keep) >= min_orthologous_species:
            eligible_genes = eligible_genes + 1
            total_num_genes = total_num_genes + len(species_of_interest_keep)
    return eligible_genes, total_num_genes


def build_representation(pwms, seq_type, min_orthologous_species):
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

        # Get the sequences for <SOURCE_SPECIES> and the orthologous species
        path = os.path.join(sequences_dir, file_name)
        with open(path, "r") as file:
            sequences_per_species = json.load(file)

        # Only keep the species that are part of <ORTHOLOGOUS_SPECIES>
        # species_of_interest = list(set(sequences_per_species.keys()) & set(ORTHOLOGOUS_SPECIES))

        # Only keep those species whose corresponding sequence is a valid one
        # (only contain A, T, G, C and are of valid length)
        species_of_interest_keep = [species for species in sequences_per_species.keys() if
                                    valid_sequence(sequences_per_species[species], seq_type)]

        # If the remaining orthologous species is less than <min_orthologous_species>, don't process
        if len(species_of_interest_keep) < min_orthologous_species:
            print(f"{i}: {gene_id} won't be processed; minimum orthologous species of {min_orthologous_species} was not met.")
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
    representation = build_representation(pwms, seq_type, 10)

    # Save the representation as a json file
    file_name = f"{seq_type}_representation_for_{SOURCE_SPECIES}.json"
    with open(file_name, "w") as file:
        json.dump(representation, file, indent=4)

    """
    eligible_genes, total_num_genes = determine_eligible_genes(PROMOTERS_DIR, "promoter", 10)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total promoter sequences")
    # There are 5596 genes, with 322077 total promoter sequences (min_orthologous_species = 10)
    # There are 19985 genes, with 363477 total promoter sequences (min_orthologous_species = 1)

    eligible_genes, total_num_genes = determine_eligible_genes(UTRs_DIR, "3UTR", 10)
    print(f"There are {eligible_genes} genes, with {total_num_genes} total 3'UTR sequences")
    # There are 5705 genes, with 331660 total 3'UTR sequences (min_orthologous_species = 10)
    # There are 19985 genes, with 372990 total 3'UTR sequences (min_orthologous_species = 1)
    """
