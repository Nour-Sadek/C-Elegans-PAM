import numpy as np
import pandas as pd
import os
import pickle
import sys

# Constants

SPECIES = "caenorhabditis_elegans"  # or "saccharomyces_cerevisiae" but only with TF
MOTIF_TYPE = "RBP"  # "TF" or "RBP"


# Read in PFMs from JASPAR
def read_jaspar_pfm(filepath):
    with open(filepath, "r") as file:
        lines = file.readlines()

    # Extract motif name from header
    header = lines[0]
    identifier, identifier_dot_motif_name = header[1:].strip().split('\t')  # remove '>' and split
    motif_name = identifier_dot_motif_name.replace(f"{identifier}.", "", 1)  # remove prefix only once

    # Extract the PFM
    content = lines[1:]
    pfm = [list(map(float, line.strip().split()[2:-1])) for line in content]

    pfm = np.array(pfm)  # shape (4 x motif_length)

    return motif_name, pfm


# Read in PFMs from YetFaSCo
def read_yetfasco_pfm(file_path):
    # Read in the file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Get the motif name
    motif_name = os.path.splitext(file_name)[0]

    # Strip whitespace and split values
    pfm = [list(map(float, line.strip().split()[1:])) for line in lines]

    # Convert to numpy array of shape (4, motif_length), 4 because A, T, G, C
    pfm = np.array(pfm)

    return motif_name, pfm


def read_rbp_pfm(file_path):
    # Read in the file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Get the matrix values
    content = lines[1:]
    pfm = [list(map(float, line.strip().split()[1:])) for line in content]
    pfm = np.array(pfm)  # shape (motif_length x 4)
    pfm = pfm.T  # shape (4 x motif_length)

    return pfm


def pfm_to_pwm(pfm, given_as_counts=True, pseudo_counts=1):
    """pfm_matrix is of the shape (4, motif_length) where the rows represent the nucleotides and the columns represent
    the position of the nucleotides in the motif.
    Background frequency of 0.25 per nucleotide is used."""

    if given_as_counts:
        pfm_counts = pfm
    else:
        # Find the smallest non-zero number in the PFM probability matrix
        min_p = np.min(pfm[pfm > 0])
        # Determine the effective number of observations Ne
        Ne = 1 / min_p
        # Estimate counts from probabilities in pfm_matrix
        pfm_counts = Ne * pfm

    # Add pseudo-counts to total and per base (0.25 * pseudo-counts)
    smoothed_counts = pfm_counts + (0.25 * pseudo_counts)
    total_counts = smoothed_counts.sum(axis=0)  # compute total counts per column, len = motif_length

    # Determine PWM using a log2 probability ratio
    pwm = np.log2((smoothed_counts / total_counts) / (0.25 / pseudo_counts))

    return pwm


# Go through the PFMs file by file and convert them to PWMs. Save the PWMs in a dictionary of the form:
# Key: name of motif (string)
# Value: PWM (numpy matrix)

if __name__ == "__main__":

    pwms = {}

    if SPECIES == "caenorhabditis_elegans" and MOTIF_TYPE == "RBP":

        # Reading in the information
        rbp_info = pd.read_excel("./Caenorhabditis_elegans_2025_07_08_11_58_am/RBP_Information.xlsx",
                                 usecols=["Motif_ID", "RBP_Name"])
        rbp_info = rbp_info[rbp_info["Motif_ID"] != "."]
        rbp_info = rbp_info.drop_duplicates()

        for motif_id in rbp_info["Motif_ID"].unique():
            file_name = f"{motif_id}.txt"
            dir_name = "./Caenorhabditis_elegans_2025_07_08_11_58_am/pwms_all_motifs"
            if file_name in os.listdir(dir_name):
                path = os.path.join(dir_name, file_name)
                # Get the motif_name
                motifs = rbp_info[rbp_info["Motif_ID"] == motif_id]["RBP_Name"]
                if len(motifs) > 1:
                    motif_name = "_".join(list(motifs))
                else:
                    motif_name = motifs.item()
                # Get the PWM
                pfm = read_rbp_pfm(path)
                pwms[motif_name] = pfm_to_pwm(pfm, given_as_counts=False)

    elif SPECIES == "caenorhabditis_elegans" and MOTIF_TYPE == "TF":

        dir_name = "./20250615154045_JASPAR2024_individual_matrices_1110029_jaspar"
        for file_name in os.listdir(dir_name):
            if file_name.endswith(".jaspar"):
                path = os.path.join(dir_name, file_name)
                motif_name, pfm = read_jaspar_pfm(path)
                pwms[motif_name] = pfm_to_pwm(pfm)

    elif SPECIES == "saccharomyces_cerevisiae" and MOTIF_TYPE == "TF":
        dir_name = "./saccharomyces_cerevisiae/Expert_PFMs/1.02/ALIGNED_ENOLOGO_FORMAT_PFMS"
        for file_name in os.listdir(dir_name):
            if file_name.endswith(".pfm"):
                path = os.path.join(dir_name, file_name)
                motif_name, pfm = read_jaspar_pfm(path)
                pwms[motif_name] = pfm_to_pwm(pfm)

    else:
        print("No PWMs to generate. Please check the values of the SPECIES and MOTIF_TYPE constants")
        sys.exit()

    # Save the dictionary of motif data as a pickle file
    pwms_file_name = f"PWMs_of_{MOTIF_TYPE}_for_{SPECIES}.pkl"
    with open(pwms_file_name, "wb") as file:
        pickle.dump(pwms, file)
