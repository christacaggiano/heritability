import sys
import numpy as np
import pandas as pd 

def read_phenotype_file(file_path):
    # Read the phenotype file
    data = pd.read_csv(file_path, sep="\t", header=None).dropna().values
    print(data) 
    sample_ids = data[:, 0]  # First column: Sample IDs
    phenotypes = data[:, 1]  # Second column: Phenotype values
    
    # Handle missing phenotypes (assume missing values are represented by NaN)
    na = ~np.isnan(phenotypes)  # True if phenotype is present, False if missing
    
    return sample_ids, phenotypes, na

def calculate_covariance(fileNameGRM, sample_ids, phenotypes, na, bins):
    # Initialize variables
    covarianceN = np.zeros(len(bins) - 1)
    covariance = np.zeros(len(bins) - 1)
    avgBins = np.zeros(len(bins) - 1)

    # Read binary GRM file
    with open(fileNameGRM + ".grm.bin", "rb") as grm_bin:
        n = len(sample_ids)
        for i in range(1, n + 1):
            if not na[i - 1]:
                continue

            i2 = 0.5 * i
            if i % 2 == 0:
                i2 *= (i - 1)
            else:
                i2 *= i

            for j in range(1, i):
                if not na[j - 1]:
                    continue

                cell = int(i2 + j)
                m = 4 * (cell - 1)

                # Read the value from the binary file
                grm_bin.seek(m)
                tmp = np.frombuffer(grm_bin.read(4), dtype=np.float32)[0]

                k = 20 if tmp > 0 else 0
                while k < len(bins) - 1:
                    if tmp < bins[k + 1]:
                        covarianceN[k] += 1.0
                        covariance[k] += phenotypes[i - 1] * phenotypes[j - 1]
                        avgBins[k] += tmp
                        break
                    k += 1

    print('...a total of', np.sum(covarianceN), 'comparisons')

    return covarianceN, covariance, avgBins


phenotype_file = "bmi_sample.txt"
sample_ids, phenotypes, na = read_phenotype_file(phenotype_file)

fileNameGRM = "grm/european_irish_sample"  # Base name for the .grm.bin file

bins = np.linspace(0, 1, 21)  # Example bin edges, e.g., 20 bins between 0 and 1

covarianceN, covariance, avgBins = calculate_covariance(fileNameGRM, sample_ids, phenotypes, na, bins)
print(covariance)
# grm = read_grm("grm/european_irish_sample.grm.bin", "grm/european_irish_sample.grm.id")
 