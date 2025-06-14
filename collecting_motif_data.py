import numpy as np
import os
import pickle

# Constants

PFM_DIR = "./Expert_PFMs/1.02/ALIGNED_ENOLOGO_FORMAT_PFMS"
FILES_EXTENSION = ".pfm"
SOURCE_SPECIES = "saccharomyces_cerevisiae"


def read_pfm(file_path):
    # Read in the file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Strip whitespace and split values
    pfm = [list(map(float, line.strip().split()[1:])) for line in lines]

    # Convert to numpy array of shape (4, motif_length), 4 because A, T, G, C
    pfm = np.array(pfm)

    return pfm


def pfm_to_pwm(pfm, pseudo_counts=1):
    """pfm_matrix is of the shape (4, motif_length) where the nucleotides are of the order A, T, G, C
    background frequency of 0.25 per nucleotide is used."""

    # Find the smallest non-zero number in the PFM probability matrix
    min_p = np.min(pfm[pfm > 0])
    # Determine the effective number of observations Ne
    Ne = 1 / min_p

    # Estimate counts from probabilities in pfm_matrix
    estimated_counts = Ne * pfm

    # Add pseudo-counts to total and per base (0.25 * pseudo-counts)
    smoothed_counts = estimated_counts + (0.25 * pseudo_counts)
    total_counts = Ne + pseudo_counts

    # Determine PWM using a log2 probability ratio
    pwm = np.log2((smoothed_counts / total_counts) / (0.25 / pseudo_counts))

    return pwm


# Go through the PFMs file by file and convert them to PWMs. Save the PWMs in a dictionary of the form:
# Key: name of motif (string)
# Value: PWM (numpy matrix)

pwms = {}

for file_name in os.listdir(PFM_DIR):
    if file_name.endswith(FILES_EXTENSION):
        motif_name = os.path.splitext(file_name)[0]
        path = os.path.join(PFM_DIR, file_name)
        pwms[motif_name] = pfm_to_pwm(read_pfm(path))

# Save the dictionary of motif data as a pickle file
pwms_file_name = f"PWMs_of_{SOURCE_SPECIES}.pkl"
with open(pwms_file_name, "wb") as file:
    pickle.dump(pwms, file)
