import numpy as np
import pandas as pd
import os
import pickle
import sys

# Constants

SPECIES = "caenorhabditis_elegans"
MOTIF_TYPE = "RBP"  # "TF" or "RBP"
TF_PFM_DIR = "./20250615154045_JASPAR2024_individual_matrices_1110029_jaspar"
RBP_PFM_DIR = "Caenorhabditis_elegans_2025_07_08_11_58_am"


def read_jaspar_pfm(filepath: str) -> tuple[str, np.ndarray]:
    """Return a tuple of the motif name and PFM associated with it as a numpy array, which are described in the file
    <filepath> which is a string that specifies the path to the file.

    This function parses through a file with a .jaspar extension to extract the motif name and PFM which contains the
    nucleotides in the order ACGT. These files are from the JASPAR database."""

    # Read in the file
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


def read_rbp_pfm(file_path: str) -> np.ndarray:
    """Return a 2d numpy array of the PFM that is described in the file at the path <file_path>.

    This function parses through a file with the .txt extension to extract the PFM which contains the nucleotides in the
    order ACGT. These files are from the cisBP database."""

    # Read in the file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Get the matrix values
    content = lines[1:]
    pfm = [list(map(float, line.strip().split()[1:])) for line in content]
    pfm = np.array(pfm)  # shape (motif_length x 4)
    pfm = pfm.T  # shape (4 x motif_length)

    return pfm


def pfm_to_pwm(pfm: np.ndarray, given_as_counts: bool = True, pseudo_counts: int = 1) -> np.ndarray:
    """Return a 2d numpy array where it is te result of converting <pfm> into a PWM matrix using a log ratio of motifs
    over a background frequency of 0.25 per nucleotide.

    How this conversion is done depends on whether <pfm> is given as nucleotide counts or probabilities between 0 and 1.
    <given_as_counts> is True for the former, and False for the latter. <pseudo_counts> is used for the conversion and a
    default value of 1 is used, however other values can be used as well. The conversion follows the same protocol as
    the one outlined in Alam et al. 2023 paper."""

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

    return pwm  # shape (4 x motif_length)


# Go through the PFMs file by file and convert them to PWMs. Save the PWMs in a dictionary of the form:
# Key: name of motif (string)
# Value: PWM (2d numpy matrix of shape: 4 x length of motif)

if __name__ == "__main__":

    pwms = {}

    if MOTIF_TYPE == "RBP":
        # Reading in the information about the motif ids and RBP names associated with the motifs
        rbp_info = pd.read_excel(f"./{RBP_PFM_DIR}/RBP_Information.xlsx", usecols=["Motif_ID", "RBP_Name"])
        rbp_info = rbp_info[rbp_info["Motif_ID"] != "."]
        rbp_info = rbp_info.drop_duplicates()

        for motif_id in rbp_info["Motif_ID"].unique():
            file_name = f"{motif_id}.txt"
            dir_name = f"./{RBP_PFM_DIR}/pwms_all_motifs"
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

    elif MOTIF_TYPE == "TF":
        for file_name in os.listdir(TF_PFM_DIR):
            if file_name.endswith(".jaspar"):
                path = os.path.join(TF_PFM_DIR, file_name)
                motif_name, pfm = read_jaspar_pfm(path)
                pwms[motif_name] = pfm_to_pwm(pfm, given_as_counts=True)

    else:
        print("No PWMs to generate. Please check the value of the MOTIF_TYPE constant.")
        sys.exit()

    # Save the dictionary of motif data as a pickle file
    pwms_file_name = f"PWMs_of_{MOTIF_TYPE}_for_{SPECIES}.pkl"
    with open(pwms_file_name, "wb") as file:
        pickle.dump(pwms, file)
