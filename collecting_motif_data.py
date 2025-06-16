import numpy as np
import os
import pickle

# Constants

PFM_DIR = "./20250615154045_JASPAR2024_individual_matrices_1110029_jaspar"
FILES_EXTENSION = ".jaspar"
SOURCE_SPECIES = "caenorhabditis_elegans"


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


pwms = {}

for file_name in os.listdir(PFM_DIR):
    if file_name.endswith(FILES_EXTENSION):
        path = os.path.join(PFM_DIR, file_name)
        motif_name, pfm = read_jaspar_pfm(path)
        pwms[motif_name] = pfm_to_pwm(pfm)

# Save the dictionary of motif data as a pickle file
pwms_file_name = f"PWMs_of_TFs_{SOURCE_SPECIES}.pkl"
with open(pwms_file_name, "wb") as file:
    pickle.dump(pwms, file)
