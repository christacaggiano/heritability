import numpy as np
import pandas as pd
from pandas_plink import read_grm
import matplotlib.pyplot as plt 

def read_phenotype_file(file_path):
    # Read the phenotype file
    data = pd.read_csv(file_path, delim_whitespace=True, header=None)
    sample_ids = data.iloc[:, 0].values  # First column: Sample IDs
    phenotypes = data.iloc[:, 1].values  # Second column: Phenotype values
    
    # Handle missing phenotypes (assume missing values are represented by NaN)
    na = ~np.isnan(phenotypes)  # True if phenotype is present, False if missing
    
    return sample_ids, phenotypes, na


def read_grm_file(prefix, sample_ids): 

    grm = read_grm(f"{prefix}.grm.bin", f"{prefix}.grm.id")
    grm_df = pd.DataFrame(grm[0].values)
    grm_df.index = sample_ids 
    grm_df.columns = sample_ids

    return grm_df


def calculate_covariance(grm_df, sample_ids, phenotypes, na, bins):
    # Initialize variables
    covarianceN = np.zeros(len(bins) - 1)
    covariance = np.zeros(len(bins) - 1)
    avgBins = np.zeros(len(bins) - 1)

    # Iterate over pairs of samples
    for i in range(len(sample_ids)):
        if not na[i]:
            continue

        for j in range(i):
            if not na[j]:
                continue

            # Get the GRM value for the (i, j) pair
            tmp = grm_df.loc[sample_ids[i], sample_ids[j]]

            # Start k at 0 to find the correct bin
            k = 0
            while k < len(bins) - 1:
                if tmp < bins[k + 1]:
                    covarianceN[k] += 1.0
                    covariance[k] += phenotypes[i] * phenotypes[j]
                    avgBins[k] += tmp
                    break
                k += 1

    return covarianceN, covariance, avgBins

def plot_average_covariance(bins, covariance, covarianceN, output_file):
    # Calculate average covariance per bin
    average_covariance = np.divide(covariance, covarianceN, where=covarianceN != 0)

    # Plot the average covariance for each bin
    mid_bins = 0.5 * (bins[:-1] + bins[1:])  # Midpoints of each bin for plotting

    plt.figure(figsize=(10, 6))
    plt.scatter(mid_bins, average_covariance, edgecolor="black", alpha=0.7)
    plt.xlabel("GRM Value Bins")
    plt.ylabel("Average Phenotypic Covariance")
    plt.title("Average Phenotypic Covariance Across GRM Value Bins")
    plt.savefig(output_file, format='png', dpi=300)
    plt.show()


if __name__ == "__main__": 

    phenotype_file = "bmi_sample.txt"
    sample_ids, phenotypes, na = read_phenotype_file(phenotype_file)

    # Read GRM data using pandas-plink
    grm_prefix = "grm/european_irish_sample"  # Prefix of the GRM files (e.g., "example_grm.grm.bin")
    grm_df = read_grm_file(grm_prefix, sample_ids)

    output_file = "irish_test.png"

    bins = np.linspace(0, 1, 21)  # Example bin edges, e.g., 20 bins between 0 and 1

    covarianceN, covariance, avgBins = calculate_covariance(grm_df, sample_ids, phenotypes, na, bins)
    print(covarianceN)
    print(covariance)
    print(avgBins)
    print(bins)
    # plot_average_covariance(bins, covariance, covarianceN, output_file)