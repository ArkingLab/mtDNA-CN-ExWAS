# Power Analysis for Quantitative Trait
# This script calculates and visualizes power for a quantitative trait association study
# Author: Vamsee Pillalamarri, adapted from various sources
# Sources:
# https://web.pdx.edu/~newsomj/uvclass/ho_power.pdf
# https://www.sfu.ca/~lockhart/richard/350/08_2/lectures/PowerSampleSize/web.pdf
# etc.

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Constants and parameters
GENOME_WIDE_P_VALUE = 5.0e-8
SAMPLE_SIZE = 10500
MAF_RANGE = np.arange(0.05, 0.51, 0.05)
EFFECT_SIZE_RANGE = np.arange(0, 0.301, 0.005)

def calculate_ncp(n, maf, effect_size):
    """
    Calculate non-centrality parameter.

    Args:
    n (int): Sample size
    maf (float): Minor allele frequency
    effect_size (float): Effect size in standard deviation units

    Returns:
    float: Non-centrality parameter
    """
    return 2 * maf * (1 - maf) * n * effect_size**2

def calculate_power(critical_chisq, ncp):
    """
    Calculate power.

    Args:
    critical_chisq (float): Critical chi-square value
    ncp (float): Non-centrality parameter

    Returns:
    float: Power
    """
    return 1 - stats.ncx2.cdf(critical_chisq, df=1, nc=ncp)

def generate_power_matrix(maf_range, effect_size_range, n, critical_chisq):
    """
    Generate power matrix.

    Args:
    maf_range (np.array): Range of minor allele frequencies
    effect_size_range (np.array): Range of effect sizes
    n (int): Sample size
    critical_chisq (float): Critical chi-square value

    Returns:
    np.array: Matrix of power values
    """
    power_matrix = np.zeros((len(maf_range), len(effect_size_range)))

    for i, maf in enumerate(maf_range):
        for j, effect_size in enumerate(effect_size_range):
            ncp = calculate_ncp(n, maf, effect_size)
            power_matrix[i, j] = calculate_power(critical_chisq, ncp)

    return power_matrix

def create_power_heatmap(maf_range, effect_size_range, power_matrix):
    """
    Create power heatmap.

    Args:
    maf_range (np.array): Range of minor allele frequencies
    effect_size_range (np.array): Range of effect sizes
    power_matrix (np.array): Matrix of power values

    Returns:
    tuple: Figure and Axes objects
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(power_matrix, cmap='Greys', aspect='auto', origin='lower',
                   extent=[maf_range[0], maf_range[-1], effect_size_range[0], effect_size_range[-1]])

    ax.set_xlabel("Cumulative Minor Allele Frequency")
    ax.set_ylabel("Effect Size (sd)")
    ax.set_title("Power Analysis for Quantitative Trait")

    cbar = fig.colorbar(im)
    cbar.set_label("Power")

    return fig, ax

def main():
    critical_chisq = stats.chi2.isf(GENOME_WIDE_P_VALUE, df=1)
    power_matrix = generate_power_matrix(MAF_RANGE, EFFECT_SIZE_RANGE, SAMPLE_SIZE, critical_chisq)

    fig, ax = create_power_heatmap(MAF_RANGE, EFFECT_SIZE_RANGE, power_matrix)

    # Display the plot
    plt.show()

    # Optionally save the plot
    # plt.savefig("quantitative_power_analysis.pdf")

if __name__ == "__main__":
    main()
