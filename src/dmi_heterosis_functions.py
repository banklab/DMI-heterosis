# ------------------------------------
# title: dmi_heterosis_functions.py
# author: Julio Ayala Lopez
# Created on: 08.09.2025
# Functions to simulate population dynamics in discrete time for overdominance and Dobzhansky-Muller incompatibilities in two and three loci
# ------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import poisson
import pandas as pd
import seaborn as sns
import warnings

#from IPython.display import set_matplotlib_formats
#%config InlineBackend.figure_formats = ['svg', 'pdf']

np.set_printoptions(precision=4)

# Standarize font
plt.rcParams['font.sans-serif'] 
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

# Exclude warnings
warnings.filterwarnings('ignore')


def haplo_to_allele(haplo_freq):
    '''Function to transform from haplotype to allele frequencies (2 loci)
    returns frequencies and LD as: [A,B,LD]
    Parameters: 
    float[4] haplo_freq: haplotype frequencies
    '''
    return [haplo_freq[1]+haplo_freq[3], haplo_freq[2]+haplo_freq[3], haplo_freq[0]*haplo_freq[3]-haplo_freq[1]*haplo_freq[2]]


def get_alleles(haplotypes):
    '''Function to transform from haplotype to allele frequencies (3 loci)
    returns frequencies and LD as: [A,B,LD]
    Parameters: 
    float[8] haplo_freq: haplotype frequencies
    '''

    ## TODO: Potentially add LD between pairs (AB-C, A-BC)
    pA = haplotypes[[1,3,5,7]]
    pB = haplotypes[[2,3,6,7]]
    pC = haplotypes[[4,5,6,7]]
    return [sum(pA),sum(pB),sum(pC)]


def get_heterozygosity(frequency):
    '''Function to get heterozygosity from a diallelic locus'''
    return 1 - (((frequency)**2) + ((1-frequency)**2))



def rec_to_distance(r):
    '''Function to convert between recombinaiton rate and distance'''
    return (-0.5*np.log(1-2*r))

def distance_to_rec(d):
    '''Function to convert between recombinaiton rate and distance'''
    return (0.5*(1-np.exp((-2*d))))


def get_ld(haplodf, neutb=False):
    '''
    Gets LD of 
    x1*x4 - x2*x3
    x1 = bn
    x2 = bN
    x3 = Bn
    x4 = BN
    
    pBN - pBpN
    pB = Bn + BN
    pN = bN + BN
    '''
    #x1*x4 - x2*x3
    #x1 = bn
    #x2 = bN
    #x3 = Bn
    #x4 = BN
    
    #pBN - pBpN
    # pB = Bn + BN
    # pN = bN + BN
    if neutb:
        haplodf['x_bn'] = haplodf.Anb + haplodf.anb
        haplodf['x_bN'] = haplodf.ANb + haplodf.aNb
        haplodf['x_Bn'] = haplodf.AnB + haplodf.anB
        haplodf['x_BN'] = haplodf.ANB + haplodf.aNB
    else:
        haplodf['x_bn'] = haplodf.Abn + haplodf.abn
        haplodf['x_bN'] = haplodf.AbN + haplodf.abN
        haplodf['x_Bn'] = haplodf.ABn + haplodf.aBn
        haplodf['x_BN'] = haplodf.ABN + haplodf.aBN
    
    
    haplodf['ld_bn'] = (haplodf.x_bn*haplodf.x_BN) - (haplodf.x_bN*haplodf.x_Bn)
    
    return haplodf


def simulate_one_gen(haplo_freq, selection_mat, mig_vec, rec):
    '''Function to simulate one generation for two loci.
    Returns a vector of haplotype frequencies ab, aB, Ab, AB
    Parameters:
    float[4] haplo_freq: initial frequencies
    float[4][4] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from aB and Ab.
    float rec: recombination rate
    '''
    # Life cycle: 1. Parental gametes to newborn genotypes w random mating; 2. Selection, 3. Migration, 4. Recombination
    
    # Determine genotype frequencies assuming random mating for all 16 genotypes
    geno_mat = np.outer(haplo_freq, haplo_freq)

    # Selection step
    selection = selection_mat * geno_mat
    fit_mean = np.sum(selection)
    geno_mat = selection / fit_mean

    # Migration step (if adult genotypes arrive after selection) -> Should the negative term apply to all?
    geno_mat *= (1 - np.sum(mig_vec))
    geno_mat[1, 1] += mig_vec[0]
    geno_mat[2, 2] += mig_vec[1]

    # Return to gametes/haplotypes with recombination
    haplo_freq = np.array([
        np.sum(geno_mat[0, :] / 2 + geno_mat[:, 0] / 2) - rec * (geno_mat[0, 3] + geno_mat[3, 0]) / 2 + rec * (geno_mat[1, 2] + geno_mat[2, 1]) / 2,
        np.sum(geno_mat[1, :] / 2 + geno_mat[:, 1] / 2) - rec * (geno_mat[1, 2] + geno_mat[2, 1]) / 2 + rec * (geno_mat[3, 0] + geno_mat[0, 3]) / 2,
        np.sum(geno_mat[2, :] / 2 + geno_mat[:, 2] / 2) - rec * (geno_mat[1, 2] + geno_mat[2, 1]) / 2 + rec * (geno_mat[3, 0] + geno_mat[0, 3]) / 2,
        np.sum(geno_mat[3, :] / 2 + geno_mat[:, 3] / 2) - rec * (geno_mat[0, 3] + geno_mat[3, 0]) / 2 + rec * (geno_mat[1, 2] + geno_mat[2, 1]) / 2
    ])

    return haplo_freq

# Test the function
#print(simulate_one_gen(np.array([0, 0.5, 0.5, 0]), selection_mat, mig_vec, rec))

#
def simulate_one_gen_stoch(haplo_freq, selection_mat, mig_vec, popsize, rec_rate):
    '''Function to simulate one generation for two loci, with multinomial sampling for the population composition.
    Returns a vector of haplotype frequencies ab, aB, Ab, AB
    Parameters:
    float[4] haplo_freq: initial frequencies
    float[4][4] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from aB and Ab.
    integer popsize: population size
    float rec: recombination rate
    '''
    # Convert migration vector from rate to absolute number of migrants
    migno_vec = mig_vec * popsize

    # Determine genotype frequencies assuming random mating
    geno_mat = np.outer(haplo_freq, haplo_freq) # Get all 16 genotype frequencies

    # Selection step
    selection = selection_mat * geno_mat 
    fit_mean = np.sum(selection)
    geno_mat = selection / fit_mean

    # Determine number of migrants that arrive this generation
    migrants = [poisson.rvs(migno_vec[0]), poisson.rvs(migno_vec[1])]

    # Determine population composition after selection (leaving space for migrants)
    flat_genfreq = geno_mat.flatten()
    
    popsize_matrix = np.reshape(np.random.multinomial(popsize - sum(migrants), flat_genfreq), (4, 4))
        
    # Introduce migrant genotypes
    popsize_matrix[1, 1] += migrants[0]
    popsize_matrix[2, 2] += migrants[1]

    # Return to genotype frequencies
    geno_mat = popsize_matrix / np.sum(popsize_matrix)

    # Return to gametes/haplotypes with recombination
    haplo_freq = np.array([
        np.sum(geno_mat[0, :] / 2 + geno_mat[:, 0] / 2) - rec_rate * (geno_mat[0, 3] + geno_mat[3, 0]) / 2 + rec_rate * (geno_mat[1, 2] + geno_mat[2, 1]) / 2,
        np.sum(geno_mat[1, :] / 2 + geno_mat[:, 1] / 2) - rec_rate * (geno_mat[1, 2] + geno_mat[2, 1]) / 2 + rec_rate * (geno_mat[3, 0] + geno_mat[0, 3]) / 2,
        np.sum(geno_mat[2, :] / 2 + geno_mat[:, 2] / 2) - rec_rate * (geno_mat[1, 2] + geno_mat[2, 1]) / 2 + rec_rate * (geno_mat[3, 0] + geno_mat[0, 3]) / 2,
        np.sum(geno_mat[3, :] / 2 + geno_mat[:, 3] / 2) - rec_rate * (geno_mat[0, 3] + geno_mat[3, 0]) / 2 + rec_rate * (geno_mat[1, 2] + geno_mat[2, 1]) / 2
    ])

    return haplo_freq

# Example parameters
#selection_mat = np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])  # Replace with actual selection matrix
#mig_vec = np.array([0.01, 0.01])  # Replace with actual migration vector
#popsize = 10000
#rec_rate = 0.01  # Replace with actual recombination rate

# Test the function
#initial_haplo_freq = np.array([0.885, 0.075, 0.030, 0.010])
#print(simulate_one_gen_stoch(initial_haplo_freq, selection_mat, mig_vec, popsize, rec_rate))



def simulate_pop(haplo_freq0, t_max, selection_mat, mig_vec, rec_rate, popsize=0):
    '''
    Simulate a population up to t_max timesteps
    Parameters:
    float[4] haplo_freq0: initial frequencies
    int t_max: maximum number of generations
    float[4][4] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from aB and Ab.
    integer popsize: population size. 
    float rec: recombination rate
    '''
    # Create the matrix to store results
    pop_vector = [haplo_freq0]
    pop_new = haplo_freq0

    # Run the simulation until generation t_max
    for i in range(1, t_max + 1):
        if popsize==0: # If the parameter is 0, deterministic
            # Update the population one generation later
            pop_new = simulate_one_gen(pop_new, selection_mat, mig_vec, rec_rate)
        else:
            pop_new = simulate_one_gen_stoch(pop_new, selection_mat, mig_vec, popsize, rec_rate)
        pop_vector.append(pop_new)
        # Stop simulation early if allele frequency change is very small
        if np.sum(np.abs(pop_vector[i] - pop_vector[i - 1])) < 0.000001:
            break

    # Convert the result list to a numpy array
    pop_vector = np.array(pop_vector)
    
    # Define row and column names (as labels in the output)
    return pop_vector


def simulate_one_gen_3(haplo_freq, selection_mat, mig_vec, rec):
    '''
    Simulate one generation for three loci
    Parameters:
    float[8] haplo_freq0: initial frequencies
    float[8][8] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from Abc and aBC.
    integer popsize: population size. 
    float[2] rec: recombination rate, [ab, bc]
    '''
    # Life cycle: 1. Parental gametes to newborn genotypes w random mating; 2. Selection, 3. Migration, 4. Recombination

    geno_mat = np.outer(haplo_freq, haplo_freq)

    # Determine genotype frequencies assuming random mating for all 16 genotypes
    geno_mat = np.outer(haplo_freq, haplo_freq)

    # Selection step
    selection = selection_mat * geno_mat
    fit_mean = np.sum(selection)
    geno_mat = selection / fit_mean

    # Migration step (if adult genotypes arrive after selection) -> Should the negative term apply to all?
    geno_mat *= (1 - np.sum(mig_vec))
    geno_mat[1, 1] += mig_vec[0]
    geno_mat[6, 6] += mig_vec[1]


    # Recombination terms
    rec_vec = np.zeros(8)
    r_ab = rec[0]
    r_bc = rec[1]

    ## abc
    # Loses
    rec_vec[0] = (-r_bc*(geno_mat[0,6] + geno_mat[6,0])/2         # abc/aBC -> abC + aBc
    - r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC -> abC + Abc
    - r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC -> abC + Abc
    - r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # abc/ABc -> Abc + aBc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[0,7] + geno_mat[7,0])/2 # abc/ABC -> aBc + AbC lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # Abc/aBc -> ABc + abc
    + r_ab*(1-r_bc)*(geno_mat[1,6] + geno_mat[6,1])/2             # Abc/aBC -> ABC + abc
    + r_bc*(geno_mat[2,4] + geno_mat[4,2])/2                      # abC/aBc -> abc + aBC
    + r_bc*(1-r_ab)*(geno_mat[3,4] + geno_mat[4,3])/2             # abC/ABc -> abc + ABC
    + r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # abC/Abc -> abc + AbC
    + r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # abC/Abc -> abc + AbC
    + r_ab*r_bc*(geno_mat[2,5] + geno_mat[5,2])/2)                # aBc/AbC -> abc + ABC
    
    
    ## Abc
    # Loses
    rec_vec[1] = (-r_bc*(geno_mat[1,7] + geno_mat[7,1])/2         # Abc/ABC -> AbC + ABc
    - r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC -> AbC + abc
    - r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC -> AbC + abc
    - r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # Abc/aBc -> ABc + abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[1,6] + geno_mat[6,1])/2 # Abc/aBC lose when rec a, c and b
    # Gains 
    + r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # abc/ABc
    + r_ab*(1-r_bc)*(geno_mat[0,7] + geno_mat[7,0])/2             # abc/ABC 
    + r_bc*(geno_mat[5,3] + geno_mat[3,5])/2                      # AbC/ABc
    + r_bc*(1-r_ab)*(geno_mat[2,5] + geno_mat[5,2])/2             # AbC/aBc
    + r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*r_bc*(geno_mat[3,4] + geno_mat[4,3])/2)                # ABc/abC
    
    
    ## aBc
    # Loses
    rec_vec[2] = (-r_bc*(geno_mat[2,4] + geno_mat[4,2])/2         # aBc/abC
    - r_bc*(1-r_ab)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    - r_ab*(1-r_bc)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    - r_ab*(geno_mat[2,1] + geno_mat[1,2])/2                      # aBc/Abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[2,5] + geno_mat[5,2])/2 # aBc/AbC lose when rec a, c and b
    # Gains 
    + r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # ABc/abc
    + r_ab*(1-r_bc)*(geno_mat[3,4] + geno_mat[4,3])/2             # ABc/abC 
    + r_bc*(geno_mat[6,0] + geno_mat[0,6])/2                      # aBC/abc
    + r_bc*(1-r_ab)*(geno_mat[6,1] + geno_mat[1,6])/2             # aBC/Abc
    + r_bc*(1-r_ab)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc
    + r_ab*(1-r_bc)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc
    + r_ab*r_bc*(geno_mat[0,7] + geno_mat[7,0])/2)                # abc/ABC
    
    ## ABc
    # Loses
    rec_vec[3] = (-r_bc*(geno_mat[3,5] + geno_mat[5,3])/2         # ABc/AbC
    - r_bc*(1-r_ab)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    - r_ab*(1-r_bc)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    - r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # ABc/abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[3,4] + geno_mat[4,3])/2 # ABc/abC lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # aBc/Abc
    + r_ab*(1-r_bc)*(geno_mat[2,5] + geno_mat[5,2])/2             # aBc/AbC
    + r_bc*(geno_mat[7,1] + geno_mat[1,7])/2                      # ABC/Abc
    + r_bc*(1-r_ab)*(geno_mat[7,0] + geno_mat[0,7])/2             # ABC/abc
    + r_bc*(1-r_ab)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc
    + r_ab*(1-r_bc)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc
    + r_ab*r_bc*(geno_mat[1,6] + geno_mat[6,1])/2)                # Abc/aBC
    
    ## abC
    # Loses
    rec_vec[4] = (-r_bc*(geno_mat[4,2] + geno_mat[2,4])/2         # abC/aBc
    - r_bc*(1-r_ab)*(geno_mat[4,1] + geno_mat[1,4])/2             # abC/Abc
    - r_ab*(1-r_bc)*(geno_mat[4,1] + geno_mat[1,4])/2             # abC/Abc
    - r_ab*(geno_mat[4,7] + geno_mat[7,4])/2                      # abC/ABC
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[4,3] + geno_mat[3,4])/2 # abC/ABc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[5,6] + geno_mat[6,5])/2                      # AbC/aBC
    + r_ab*(1-r_bc)*(geno_mat[5,2] + geno_mat[2,5])/2             # AbC/aBc
    + r_bc*(geno_mat[0,6] + geno_mat[6,0])/2                      # abc/aBC
    + r_bc*(1-r_ab)*(geno_mat[0,7] + geno_mat[7,0])/2             # abc/ABC
    + r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*r_bc*(geno_mat[6,1] + geno_mat[1,6])/2)                # aBC/Abc
    
    ## AbC
    # Loses
    rec_vec[5] = (-r_bc*(geno_mat[5,3] + geno_mat[3,5])/2         # AbC/ABc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[5,0] + geno_mat[0,5])/2             # AbC/abc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[5,0] + geno_mat[0,5])/2             # AbC/abc "
    - r_ab*(geno_mat[5,6] + geno_mat[6,5])/2                      # AbC/aBC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[5,2] + geno_mat[2,5])/2 # AbC/aBc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[4,7] + geno_mat[7,4])/2                      # abC/ABC
    + r_ab*(1-r_bc)*(geno_mat[4,3] + geno_mat[3,4])/2             # abC/ABc
    + r_bc*(geno_mat[1,7] + geno_mat[7,1])/2                      # Abc/ABC
    + r_bc*(1-r_ab)*(geno_mat[1,6] + geno_mat[6,1])/2             # Abc/aBC
    + r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC
    + r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC
    + r_ab*r_bc*(geno_mat[7,0] + geno_mat[0,7])/2)                # ABC/abc
    
    
    ## aBC
    # Loses
    rec_vec[6] = (-r_bc*(geno_mat[6,0] + geno_mat[0,6])/2         # aBC/abc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc "
    - r_ab*(geno_mat[6,5] + geno_mat[5,6])/2                      # aBC/AbC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[6,1] + geno_mat[1,6])/2 # aBC/Abc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[7,4] + geno_mat[4,7])/2                      # ABC/abC
    + r_ab*(1-r_bc)*(geno_mat[7,0] + geno_mat[0,7])/2             # ABC/abc
    + r_bc*(geno_mat[2,4] + geno_mat[4,2])/2                      # aBc/abC
    + r_bc*(1-r_ab)*(geno_mat[2,5] + geno_mat[5,2])/2             # aBc/AbC
    + r_bc*(1-r_ab)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    + r_ab*(1-r_bc)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    + r_ab*r_bc*(geno_mat[4,3] + geno_mat[3,4])/2)                # abC/ABc
    
    
    ## ABC
    # Loses
    rec_vec[7] = (-r_bc*(geno_mat[7,1] + geno_mat[1,7])/2         # ABC/Abc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc "
    - r_ab*(geno_mat[7,4] + geno_mat[4,7])/2                      # ABC/abC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[7,0] + geno_mat[0,7])/2 # ABC/abc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[6,5] + geno_mat[5,6])/2                      # aBC/AbC
    + r_ab*(1-r_bc)*(geno_mat[6,1] + geno_mat[1,6])/2             # aBC/Abc
    + r_bc*(geno_mat[3,5] + geno_mat[5,3])/2                      # ABc/AbC
    + r_bc*(1-r_ab)*(geno_mat[3,4] + geno_mat[4,3])/2             # ABc/abC
    + r_bc*(1-r_ab)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    + r_ab*(1-r_bc)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    + r_ab*r_bc*(geno_mat[5,2] + geno_mat[2,5])/2)                # AbC/aBc

    
    # Return to gametes/haplotypes with recombination
    haplo_freq = np.array([
        np.sum(geno_mat[0, :] / 2 + geno_mat[:, 0] / 2) + rec_vec[0],
        np.sum(geno_mat[1, :] / 2 + geno_mat[:, 1] / 2) + rec_vec[1],
        np.sum(geno_mat[2, :] / 2 + geno_mat[:, 2] / 2) + rec_vec[2],
        np.sum(geno_mat[3, :] / 2 + geno_mat[:, 3] / 2) + rec_vec[3],
        np.sum(geno_mat[4, :] / 2 + geno_mat[:, 4] / 2) + rec_vec[4],
        np.sum(geno_mat[5, :] / 2 + geno_mat[:, 5] / 2) + rec_vec[5],
        np.sum(geno_mat[6, :] / 2 + geno_mat[:, 6] / 2) + rec_vec[6],
        np.sum(geno_mat[7, :] / 2 + geno_mat[:, 7] / 2) + rec_vec[7]
    ])

    return haplo_freq

#Stochastic simulation
def simulate_one_gen_stoch_3(haplo_freq, selection_mat, mig_vec, popsize, rec_rate):
    # Convert migration vector from rate to absolute number of migrants
    migno_vec = mig_vec * popsize

    # Determine genotype frequencies assuming random mating
    geno_mat = np.outer(haplo_freq, haplo_freq) # Get all 16 genotype frequencies

    # Selection step
    selection = selection_mat * geno_mat 
    fit_mean = np.sum(selection)
    geno_mat = selection / fit_mean

    # Determine number of migrants that arrive this generation
    migrants = [poisson.rvs(migno_vec[0]), poisson.rvs(migno_vec[1])]

    # Determine population composition after selection (leaving space for migrants)
    flat_genfreq = geno_mat.flatten()
    
    popsize_matrix = np.reshape(np.random.multinomial(popsize - sum(migrants), flat_genfreq), (8, 8))
        
    # Introduce migrant genotypes
    popsize_matrix[1, 1] += migrants[0]
    popsize_matrix[6, 6] += migrants[1]

    # Return to genotype frequencies
    geno_mat = popsize_matrix / np.sum(popsize_matrix)

    # Recombination terms
    # Recombination terms
    rec_vec = np.zeros(8)
    r_ab = rec_rate[0]
    r_bc = rec_rate[1]

    ## abc
    # Loses
    rec_vec[0] = (-r_bc*(geno_mat[0,6] + geno_mat[6,0])/2         # abc/aBC -> abC + aBc
    - r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC -> abC + Abc
    - r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC -> abC + Abc
    - r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # abc/ABc -> Abc + aBc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[0,7] + geno_mat[7,0])/2 # abc/ABC -> aBc + AbC lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # Abc/aBc -> ABc + abc
    + r_ab*(1-r_bc)*(geno_mat[1,6] + geno_mat[6,1])/2             # Abc/aBC -> ABC + abc
    + r_bc*(geno_mat[2,4] + geno_mat[4,2])/2                      # abC/aBc -> abc + aBC
    + r_bc*(1-r_ab)*(geno_mat[3,4] + geno_mat[4,3])/2             # abC/ABc -> abc + ABC
    + r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # abC/Abc -> abc + AbC
    + r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # abC/Abc -> abc + AbC
    + r_ab*r_bc*(geno_mat[2,5] + geno_mat[5,2])/2)                # aBc/AbC -> abc + ABC
    
    
    ## Abc
    # Loses
    rec_vec[1] = (-r_bc*(geno_mat[1,7] + geno_mat[7,1])/2         # Abc/ABC -> AbC + ABc
    - r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC -> AbC + abc
    - r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC -> AbC + abc
    - r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # Abc/aBc -> ABc + abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[1,6] + geno_mat[6,1])/2 # Abc/aBC lose when rec a, c and b
    # Gains 
    + r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # abc/ABc
    + r_ab*(1-r_bc)*(geno_mat[0,7] + geno_mat[7,0])/2             # abc/ABC 
    + r_bc*(geno_mat[5,3] + geno_mat[3,5])/2                      # AbC/ABc
    + r_bc*(1-r_ab)*(geno_mat[2,5] + geno_mat[5,2])/2             # AbC/aBc
    + r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*r_bc*(geno_mat[3,4] + geno_mat[4,3])/2)                # ABc/abC
    
    
    ## aBc
    # Loses
    rec_vec[2] = (-r_bc*(geno_mat[2,4] + geno_mat[4,2])/2         # aBc/abC
    - r_bc*(1-r_ab)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    - r_ab*(1-r_bc)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    - r_ab*(geno_mat[2,1] + geno_mat[1,2])/2                      # aBc/Abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[2,5] + geno_mat[5,2])/2 # aBc/AbC lose when rec a, c and b
    # Gains 
    + r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # ABc/abc
    + r_ab*(1-r_bc)*(geno_mat[3,4] + geno_mat[4,3])/2             # ABc/abC 
    + r_bc*(geno_mat[6,0] + geno_mat[0,6])/2                      # aBC/abc
    + r_bc*(1-r_ab)*(geno_mat[6,1] + geno_mat[1,6])/2             # aBC/Abc
    + r_bc*(1-r_ab)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc
    + r_ab*(1-r_bc)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc
    + r_ab*r_bc*(geno_mat[0,7] + geno_mat[7,0])/2)                # abc/ABC
    
    ## ABc
    # Loses
    rec_vec[3] = (-r_bc*(geno_mat[3,5] + geno_mat[5,3])/2         # ABc/AbC
    - r_bc*(1-r_ab)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    - r_ab*(1-r_bc)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    - r_ab*(geno_mat[0,3] + geno_mat[3,0])/2                      # ABc/abc
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[3,4] + geno_mat[4,3])/2 # ABc/abC lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[1,2] + geno_mat[2,1])/2                      # aBc/Abc
    + r_ab*(1-r_bc)*(geno_mat[2,5] + geno_mat[5,2])/2             # aBc/AbC
    + r_bc*(geno_mat[7,1] + geno_mat[1,7])/2                      # ABC/Abc
    + r_bc*(1-r_ab)*(geno_mat[7,0] + geno_mat[0,7])/2             # ABC/abc
    + r_bc*(1-r_ab)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc
    + r_ab*(1-r_bc)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc
    + r_ab*r_bc*(geno_mat[1,6] + geno_mat[6,1])/2)                # Abc/aBC
    
    ## abC
    # Loses
    rec_vec[4] = (-r_bc*(geno_mat[4,2] + geno_mat[2,4])/2         # abC/aBc
    - r_bc*(1-r_ab)*(geno_mat[4,1] + geno_mat[1,4])/2             # abC/Abc
    - r_ab*(1-r_bc)*(geno_mat[4,1] + geno_mat[1,4])/2             # abC/Abc
    - r_ab*(geno_mat[4,7] + geno_mat[7,4])/2                      # abC/ABC
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[4,3] + geno_mat[3,4])/2 # abC/ABc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[5,6] + geno_mat[6,5])/2                      # AbC/aBC
    + r_ab*(1-r_bc)*(geno_mat[5,2] + geno_mat[2,5])/2             # AbC/aBc
    + r_bc*(geno_mat[0,6] + geno_mat[6,0])/2                      # abc/aBC
    + r_bc*(1-r_ab)*(geno_mat[0,7] + geno_mat[7,0])/2             # abc/ABC
    + r_bc*(1-r_ab)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*(1-r_bc)*(geno_mat[0,5] + geno_mat[5,0])/2             # abc/AbC
    + r_ab*r_bc*(geno_mat[6,1] + geno_mat[1,6])/2)                # aBC/Abc
    
    ## AbC
    # Loses
    rec_vec[5] = (-r_bc*(geno_mat[5,3] + geno_mat[3,5])/2         # AbC/ABc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[5,0] + geno_mat[0,5])/2             # AbC/abc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[5,0] + geno_mat[0,5])/2             # AbC/abc "
    - r_ab*(geno_mat[5,6] + geno_mat[6,5])/2                      # AbC/aBC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[5,2] + geno_mat[2,5])/2 # AbC/aBc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[4,7] + geno_mat[7,4])/2                      # abC/ABC
    + r_ab*(1-r_bc)*(geno_mat[4,3] + geno_mat[3,4])/2             # abC/ABc
    + r_bc*(geno_mat[1,7] + geno_mat[7,1])/2                      # Abc/ABC
    + r_bc*(1-r_ab)*(geno_mat[1,6] + geno_mat[6,1])/2             # Abc/aBC
    + r_bc*(1-r_ab)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC
    + r_ab*(1-r_bc)*(geno_mat[1,4] + geno_mat[4,1])/2             # Abc/abC
    + r_ab*r_bc*(geno_mat[7,0] + geno_mat[0,7])/2)                # ABC/abc
    
    
    ## aBC
    # Loses
    rec_vec[6] = (-r_bc*(geno_mat[6,0] + geno_mat[0,6])/2         # aBC/abc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[6,3] + geno_mat[3,6])/2             # aBC/ABc "
    - r_ab*(geno_mat[6,5] + geno_mat[5,6])/2                      # aBC/AbC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[6,1] + geno_mat[1,6])/2 # aBC/Abc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[7,4] + geno_mat[4,7])/2                      # ABC/abC
    + r_ab*(1-r_bc)*(geno_mat[7,0] + geno_mat[0,7])/2             # ABC/abc
    + r_bc*(geno_mat[2,4] + geno_mat[4,2])/2                      # aBc/abC
    + r_bc*(1-r_ab)*(geno_mat[2,5] + geno_mat[5,2])/2             # aBc/AbC
    + r_bc*(1-r_ab)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    + r_ab*(1-r_bc)*(geno_mat[2,7] + geno_mat[7,2])/2             # aBc/ABC
    + r_ab*r_bc*(geno_mat[4,3] + geno_mat[3,4])/2)                # abC/ABc
    
    
    ## ABC
    # Loses
    rec_vec[7] = (-r_bc*(geno_mat[7,1] + geno_mat[1,7])/2         # ABC/Abc switch b and c
    - r_bc*(1-r_ab)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc switch a and c
    - r_ab*(1-r_bc)*(geno_mat[7,2] + geno_mat[2,7])/2             # ABC/aBc "
    - r_ab*(geno_mat[7,4] + geno_mat[4,7])/2                      # ABC/abC Switch a and b
    - (r_ab + r_bc - r_ab*r_bc)*(geno_mat[7,0] + geno_mat[0,7])/2 # ABC/abc lose when rec a, c and b
    # Gains
    + r_ab*(geno_mat[6,5] + geno_mat[5,6])/2                      # aBC/AbC
    + r_ab*(1-r_bc)*(geno_mat[6,1] + geno_mat[1,6])/2             # aBC/Abc
    + r_bc*(geno_mat[3,5] + geno_mat[5,3])/2                      # ABc/AbC
    + r_bc*(1-r_ab)*(geno_mat[3,4] + geno_mat[4,3])/2             # ABc/abC
    + r_bc*(1-r_ab)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    + r_ab*(1-r_bc)*(geno_mat[3,6] + geno_mat[6,3])/2             # ABc/aBC
    + r_ab*r_bc*(geno_mat[5,2] + geno_mat[2,5])/2)                # AbC/aBc

    # Return to gametes/haplotypes with recombination
    haplo_freq = np.array([
        np.sum(geno_mat[0, :] / 2 + geno_mat[:, 0] / 2) + rec_vec[0],
        np.sum(geno_mat[1, :] / 2 + geno_mat[:, 1] / 2) + rec_vec[1],
        np.sum(geno_mat[2, :] / 2 + geno_mat[:, 2] / 2) + rec_vec[2],
        np.sum(geno_mat[3, :] / 2 + geno_mat[:, 3] / 2) + rec_vec[3],
        np.sum(geno_mat[4, :] / 2 + geno_mat[:, 4] / 2) + rec_vec[4],
        np.sum(geno_mat[5, :] / 2 + geno_mat[:, 5] / 2) + rec_vec[5],
        np.sum(geno_mat[6, :] / 2 + geno_mat[:, 6] / 2) + rec_vec[6],
        np.sum(geno_mat[7, :] / 2 + geno_mat[:, 7] / 2) + rec_vec[7]
    ])

    return haplo_freq


# max_gen = 25000
# haplotype_freq = np.array([0.07, 0.05, 0.22, 0.41, 0.06, 0.04, 0.13, 0.02])
# popsize=10000
# print(simulate_one_gen_stoch_3(haplotype_freq, selection_mat_3, mig_vec, popsize, rec))


# Test the function
#print(simulate_one_gen_3(np.array([0.07, 0.05, 0.12, 0.51, 0.06, 0.04, 0.13, 0.02]), selection_mat_3, mig_vec, rec))



def simulate_pop_3(haplo_freq0, t_max, selection_mat, mig_vec, rec_rate, popsize=0):
    '''
    Simulate a 3 locus population until t_max
    Parameters:
    float[8] haplo_freq0: initial frequencies
    int t_max: maximum number of generations
    float[8][8] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from aB and Ab.
    integer popsize: population size. If non-0, runs stochastically
    float rec: recombination rate
    
    '''
    # Create the matrix to store results
    pop_vector = [haplo_freq0]
    pop_new = haplo_freq0

    # Run the simulation until generation t_max
    for i in range(1, t_max + 1):
        if popsize==0: # If the parameter is 0, do the deterministic case
            # Update the population one generation later
            pop_new = simulate_one_gen_3(pop_new, selection_mat, mig_vec, rec_rate)
        else:
            pop_new = simulate_one_gen_stoch_3(pop_new, selection_mat, mig_vec, popsize, rec_rate)
        pop_vector.append(pop_new)
        # Stop simulation early if allele frequency change is in the order of magnitude of 1e-10
        if np.sum(np.abs(pop_vector[i] - pop_vector[i - 1])) < 1e-10:
            break

    # Convert the result list to a numpy array
    pop_vector = np.array(pop_vector)
    
    # Define row and column names (as labels in the output)
    return pop_vector


def get_selection_mat(overdom_vec, epi_vec, neutb=False):
    '''
    Function to get the selection matrix given overdominance and epistatic coefficients
    Parameters:
    float[2] overdom_vec: Vector with alpha and beta coefficients
    float[4] epi_vec: Vector with epistatic coefficients
    boolean neutb: Flag to determine the architecture. If false, ABN; else ANB.
    '''
    alpha = overdom_vec[0]
    beta = overdom_vec[1]
    gamma1 = epi_vec[0]
    gamma2 = epi_vec[1]
    gamma3 = epi_vec[2]
    gamma4 = epi_vec[3]
        #ab                                       Ab                                       aB                                       AB
    selection_mat = np.array([
        [1,                                        1 + alpha,                               1 + beta,                               (1 + alpha) * (1 + beta) * (1 - gamma1)],
        [1 + alpha,                                1,                                      (1 + alpha) * (1 + beta) * (1 - gamma1), (1 + alpha) * (1 - gamma3)],
        [1 + beta,                                (1 + alpha) * (1 + beta) * (1 - gamma1),  1,                                      (1 + beta) * (1 - gamma2)],
        [(1 + alpha) * (1 + beta) * (1 - gamma1), (1 + alpha) * (1 - gamma3),              (1 + beta) * (1 - gamma2),                1 - gamma4]
    ])

    selection_mat_3 = np.tile(selection_mat,(2,2)) # The selection matrix is a repetition of the AB matrix 
    
    if neutb: # Neutral B has a different selection matrix
        selection_mat_3 = np.array([
            [1, 1+alpha, 1, 1+alpha, (1+beta), (1+alpha)*(1+beta)*(1-gamma1), 1+beta, (1+alpha)*(1+beta)*(1-gamma1)],
            [1+alpha, 1, 1+alpha, 1, (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2)],
            [1, 1+alpha, 1, 1+alpha, (1+beta), (1+alpha)*(1+beta)*(1-gamma1), 1+beta, (1+alpha)*(1+beta)*(1-gamma1)],
            [1+alpha, 1, 1+alpha, 1, (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2)],
            [1+beta, (1+alpha)*(1+beta)*(1-gamma1), 1+beta, (1+alpha)*(1+beta)*(1-gamma1), 1, (1+alpha)*(1-gamma2), 1, (1+alpha)*(1-gamma2)],
            [(1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1-gamma2), 1-gamma4, (1+alpha)*(1-gamma2), 1-gamma4],
            [1+beta, (1+alpha)*(1+beta)*(1-gamma1), 1+beta, (1+alpha)*(1+beta)*(1-gamma1), 1, (1+alpha)*(1-gamma2), 1, (1+alpha)*(1-gamma2)],
            [(1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1+beta)*(1-gamma1), (1+beta)*(1-gamma2), (1+alpha)*(1-gamma2), 1-gamma4, (1+alpha)*(1-gamma2), 1-gamma4]    
        ]) 
    return selection_mat_3



def simulate_populations(haplo_freq, 
                         max_gen,
                         selection_mat,
                         mig,      
                         recombination_list,
                         neutb=False,
                         outtype='alleles',
                         **kwargs, #rmin=-4,rmax=0,steps=20,steps=20, npoints=25
                        ):

    '''
    Simulate several populations with different recombination rates between A and B.
    For the ABN architecture, it takes a number of equally spaced 'steps' on a logistic scale.
    Parameters: 
    float[8] haplo_freq0: initial frequencies
    int t_max: maximum number of generations
    float[8][8] selection_mat: Selection coeficients for all haplotype combinations
    float[2] mig_vec: migration rate from aB and Ab.
    integer popsize: population size. If non-0, runs stochastically
    float rec: recombination rate

    '''

    defaultkwargs = {'rmin':-5, 'rmax':0, 'steps':25, 'npoints':25}    # Defaulting kwargs
    kwargs = {**kwargs, **defaultkwargs}                             # Unpack
    out = []
    outhaplo = []
    outdata = []
    alleles = ['A', 'B', 'C']

    if neutb: # Neutral B (A-N-B)
        everyn = int(100/kwargs['npoints'])
        for totr in recombination_list:
            for rab in range(0,101,everyn):
                freqs = []
                #for i in np.linspace(-4,0,20):
                r_ab = totr*rab/100
                r_bc = totr-(totr*rab/100)
                output = simulate_pop_3(haplo_freq, max_gen, selection_mat, [0,0], [r_ab, r_bc])
                t_max = len(output)
                freqs.append([t_max, r_ab, r_bc, totr, get_alleles(output[-1])])
                out.append(freqs)
                outhaplo.append([t_max, r_ab, r_bc])
                outhaplo[-1].extend(output[-1])
    
        #Extract allele data and create df
        haplodf = pd.DataFrame(outhaplo, columns=['tmax','r_ab', 'r_bc', 'anb', 'Anb', 'aNb', 'ANb', 'anB', 'AnB', 'aNB', 'ANB'])
        
        for feature_idx, feature in enumerate(out):
            for t_max, rab, rbc, totr, y_values in feature:
                for allele, y in zip(alleles, y_values):
                    outdata.append({'t_max': t_max,'r_ab': totr, 'r_bc': rbc, 'allele': allele, 'freq': y, 'pos': rab/totr, 'r_an':rab})


    else: # Neutral C (A-B-N)
        for rab in recombination_list:
            freqs = []
            for i in np.linspace(kwargs['rmin'],kwargs['rmax'],kwargs['steps']):
                output = simulate_pop_3(haplo_freq, max_gen, selection_mat, mig, [rab, 0.5*10**(i)])
                t_max = len(output)
                freqs.append([t_max, 0.5*10**(i), get_alleles(output[-1])])
                outhaplo.append([t_max, rab, 0.5*10**(i)])
                outhaplo[-1].extend(output[-1])
            out.append(freqs)
        haplodf = pd.DataFrame(outhaplo, columns=['tmax','r_ab', 'r_bn', 'abn',  'Abn', 'aBn', 'ABn', 'abN', 'AbN', 'aBN', 'ABN'])

    
        for feature_idx, feature in enumerate(out):
            for t_max, x_value, y_values in feature:
                for allele, y in zip(alleles, y_values):
                    outdata.append({'t_max': t_max, 'r_ab': recombination_list[feature_idx], 'r_bc': x_value, 'allele': allele, 'freq': y })


    # Create a DataFrame
    df = pd.DataFrame(outdata)
    df['het'] = 1 - (((df.freq)**2) + ((1-df.freq)**2)) # Calculate heterozygosity
    if neutb:
        df.loc[:, 'allele'] = df['allele'].replace({'B': 'N'})
        df.loc[:, 'allele'] = df['allele'].replace({'C': 'B'})
    else:
        df.loc[:, 'allele'] = df['allele'].replace({'C': 'N'})       # Replace 
    df['r_ab'] = [f"{i:.0e}" for i in df['r_ab']]        # Change to scientific notation      

    return df, haplodf
    
    
#####################################
# Plotting functions
#####################################


def plot_rec_freqs(haplo_freq, max_gen, selection_mat, mig, rab=0.00001):
    '''
    # UNUSED?
    Plot equilibrium frequencies at different recombination distances
    '''

    freqs = []
    for i in np.linspace(-4,0.5,20):
        output = simulate_pop_3(haplo_freq, max_gen, selection_mat, mig, [rab, 0.5*10**(i)])
        freqs.append([0.5*10**(i), get_alleles(output[-1])])
        
    plt.plot([i[0] for i in freqs], [i[1] for i in freqs], '--s', label=["A","B","C"], alpha=0.5)
    plt.ylim(0, 1)
    plt.colormaps
    plt.xscale('log')
    plt.xlabel("r_bc")
    plt.ylabel("Allele frequency")
    plt.legend()
    plt.show()
    return(freqs[-1][1])
    

def plot_dyn_haplo(haplo_freq, 
                    max_gen, 
                    selection_mat, 
                    mig, 
                    r, 
                    title="",
                    neut="C", 
                    exclude=[],
                    haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
                    **kwargs):
    '''
    Plot haplotype dynamics over time
    '''
    # Simulate population dynamics
    output = simulate_pop_3(haplo_freq, max_gen, selection_mat, mig, r)
    
    # Create a new figure explicitly
    fig, ax = plt.subplots(figsize=(8, 5))    
    x_range = np.arange(len(output))
    
    
    for i in range(len(haplotypes)):
        ax.plot(x_range, output[:, i], label=haplotypes[i])
        
    ax.set_ylim(0, np.max(output) + 0.05)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Haplotype frequency")
    ax.set_title(title)
    ax.legend()
    
    # Save the figure if a filename is provided
    if "filename" in kwargs:
        fig.savefig(kwargs["filename"], dpi=300, bbox_inches='tight')
    
    plt.show()
    
    return output[-1]



def plot_dyn_freqs(haplo_freq, 
                   max_gen, 
                   selection_mat, 
                   mig, 
                   r, 
                   title="", 
                   neut="C", 
                   palette = {'A':'#5C43AF','B':'#FE842F', 'N':'#194122', 'A,B':'#FE842F'},
                   exclude=[],
                   **kwargs):

    '''
    Plot allele frequency dynamics over time
    Parameters:
    haplo_freq: Initial haplotype frequencies
    selection_mat: Selection matrix
    mig: migration vector
    r: recombination rate vector
    [title=""]: tile
    [neut="C"]: Neutral locus (B/C)
    [palette = {'A':'#5C43AF','B':'#FE842F', 'N':'#194122', 'A,B':'#FE842F'}]: Colour palette
    [exclude=[]]: Loci to exclude 
    kwargs: 
        filename: File to save plot to
    '''
                          # abc,  Abc, aBc, ABc, abC,  AbC,  aBC,  ABC 
    output = simulate_pop_3(haplo_freq, max_gen, selection_mat, mig, r)
    
    # Show the last few lines of the data table
    #print(output[-5:])
    
    
    # Plot the output
    x_range = np.arange(len(output))
    allele_freq = np.array(list(map(get_alleles, output)))
    alleles = ["A", "B", "C"]
    alleles[alleles.index(neut)] = "N"
    
    df_alleles = pd.DataFrame(allele_freq, columns=alleles)
    df_alleles['t'] = df_alleles.index
    df_alleles = df_alleles.melt(id_vars='t', value_vars=alleles, 
                        var_name='allele', value_name='freq')

    if exclude:
        df_alleles = df_alleles[~df_alleles['allele'].isin(exclude)]
    
    palette = {'A':'#5C43AF','B':'#FE842F', 'N':'#194122', 'A,B':'#FE842F'}

    fig, ax = plt.subplots()
    figfreqs = sns.lineplot(data=df_alleles, x='t', y='freq', hue='allele',
            palette=palette, alpha=0.8)
    
    plt.legend(title='Allele')
    plt.xlabel("Generation")
    plt.ylabel("Allele frequency")
    sns.despine()
    plt.title(title)
    plt.legend(title="Allele", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, frameon=False)
    plt.show()
    if "filename" in kwargs.keys():
        fig.savefig(kwargs["filename"], dpi=300, bbox_inches='tight')
    return(get_alleles(output[-1]))



## Updated function of above, taking a dataframe as a parameter
def plot_equilibrium_over_dist(alleledf,
                            neutb=False,
                             palette = {'A':'#5C43AF','B':'#FE842F', 'N':'#194122', 'A,B':'#FE842F'},
                             yaxis ='het',
                             allelelist = ['B', 'N'],
                             **kwargs):

    
    ## Get parameters for customized plot
    legendpos = kwargs.get("legendpos", "below")
    fontsize = kwargs.get("fontsize", 18)
    markersize = kwargs.get("markersize", 16)
    linewidth = kwargs.get("linewidth", 1.5)
    dasheslist=kwargs.get("dasheslist", [(4,1),(4,1.5),(4,1.25)])
    plt.rcParams.update({'font.size': fontsize}) 
    height = kwargs.get("height", 6)
    width = kwargs.get("width", height * 1.3)  # default aspect ratio if width not provided
    aspect = width / height
    
    alleledf = alleledf[alleledf['allele'].isin(allelelist)]
    if ('A' not in allelelist and yaxis=='het'): # Change the label of B to A,B for clarity
        alleledf.loc[:, 'allele'] = alleledf['allele'].replace({'B': 'A,B'})

    #dashesdict = {'A':(2, 2),'B':(2, 3), 'N':(2, 5), 'A,B':(2, 3)}
    markerlist = ["X", "P", "^", "o"] if len(pd.unique(alleledf['r_ab']))==3 else ["o"]
    
    figfreqs = sns.relplot(data=alleledf, x='pos' if neutb else 'r_bc', y=yaxis, hue='allele',
        style='r_ab', kind="line", markers=markerlist if markersize>0 else False,
        palette=palette, alpha=0.8, linewidth=linewidth,
        dashes=dasheslist,
        #dashes=False,
        height=height, aspect=aspect)


    
    legend = figfreqs._legend  
    if legend:
        for text in legend.get_texts():
            if text.get_text() == "allele":
                if yaxis=='het':
                    text.set_text("Locus")  # Change label
                elif yaxis=='freq':
                    text.set_text("Allele")  
            elif text.get_text() == "r_ab":
                text.set_text(r"$r_{ab}$")             


        if legendpos=="below":
            sns.move_legend(
                figfreqs, "lower center",
                bbox_to_anchor=(.5, -0.05), ncol=8, title=None, frameon=False,
            )
        elif legendpos=="none":
            # Move the legend to the upper right
            #legend.set_bbox_to_anchor((0.8, 0.78))  # Position it in the upper right corner
            #legend.set_frame_on(True)
            figfreqs._legend.remove()
    

    
    
    for ax in figfreqs.axes.flat:
        for line in ax.lines:
            line.set_markersize(markersize) 

            # Change line color alpha only (line.get_color() returns RGB or RGBA string)
            lc = line.get_color()  # e.g. '#FE842F' or (r,g,b,a)
            # Convert color to RGBA tuple:
            rgba = mcolors.to_rgba(lc)
            # Replace alpha for line only, e.g. 0.3
            new_line_color = (rgba[0], rgba[1], rgba[2], 0.4)
            line.set_color(new_line_color)
            
            mfc = line.get_markerfacecolor() # Change markers to 0.8
            if isinstance(mfc, tuple) and len(mfc) == 4:
                new_mfc = (mfc[0], mfc[1], mfc[2], 0.8)  
                line.set_markerfacecolor(new_mfc)

            mec = line.get_markeredgecolor()
            if isinstance(mec, tuple) and len(mec) == 4:
                new_mec = (mec[0], mec[1], mec[2], 0.8)
                line.set_markeredgecolor(new_mec)

    if neutb:
        plt.xlabel("Relative position of N")
        if yaxis=='het':
            plt.ylabel("Heterozygosity (H)")
        elif yaxis=='freq':
            plt.ylabel("Frequency")
    
        # Get X-axis limits
        xmin, xmax = ax.get_xlim()
    
        # Add text "A" at the beginning and "B" at the end of the X-axis
        ax.text(0.005, -0.03, "A", ha='center', va='top', fontsize=fontsize-2, transform=ax.transAxes)
        ax.text(0.985, -0.03, "B", ha='center', va='top', fontsize=fontsize-2, transform=ax.transAxes)

    else:
        plt.xlabel(r"$r_{bn}$")
        if yaxis=='het':
            plt.ylabel("Heterozygosity (H)")
        elif yaxis=='freq':
            plt.ylabel("Frequency")
    
        # Get X-axis limits
        ax = figfreqs.ax
        xmin, xmax = ax.get_xlim()
    
        # Set log scale on x-axis
        ax.set_xscale('log')
    
        # ax.axhline(0.5, ls=':', linewidth=1, color='red', alpha = 0.5) # Red line at 0.5, for testing
    
        ax.text(0.01, -0.03, "AB", ha='center', va='top', fontsize=fontsize-2, transform=ax.transAxes)
        ax.text(0.98, -0.03, "→N", ha='center', va='top', fontsize=fontsize-2, transform=ax.transAxes)
    
    
    if yaxis=='het':
        plt.ylim(-0.05, 0.55)
    else:
        plt.ylim(-0.05, 1.05)
    if "title" in kwargs.keys():
        plt.title(kwargs["title"])

    plt.grid(linewidth=0.5, alpha=0.5)
    
    if "skip_show" not in kwargs.keys():
    	plt.show()
    

    if "filename" in kwargs.keys():
        figfreqs.savefig(kwargs["filename"])

    plt.rcParams.update({'font.size': 10}) ## Return font size to original
    #plt.savefig('fighet_ppt.svg')-


def load_dataset(haplotype_filename, allele_filename):
	return 1

