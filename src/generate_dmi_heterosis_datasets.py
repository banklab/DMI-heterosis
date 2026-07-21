# ------------------------------------
# title: generate_dmi_heterosis_datasets.py
# author: Julio Ayala Lopez
# Created on: 02.10.2025
# Script to generate different datasets, for codominant DMIs and heterosis
# ------------------------------------

from dmi_heterosis_functions import *
from pathlib import Path


# Setting up sumulation parameters
out_file_prefix = '../Output_new/'
max_gen = 10000
rab = [1e-6,1e-2,1e-1]
mig_vec = [0,0]

            
for pA in [0.8,0.2,0.55,0.45,0.6,0.4,0.5,0.7,0.3,0.51,0.49]:
                  # Haplotype frequencies
                  # abc,  Abc, aBc, ABc, abC,  AbC,  aBC,             ABC 
    haplo_freq = [0. ,    pA,  0. , 0. , 0. ,  0. ,  round(1 - pA,2), 0.]
    pA0 = haplo_freq[1] + haplo_freq[3] + haplo_freq[5] + haplo_freq[7]
    pB0 = haplo_freq[2] + haplo_freq[3] + haplo_freq[6] + haplo_freq[7]
    pC0 = haplo_freq[4] + haplo_freq[5] + haplo_freq[6] + haplo_freq[7]

    # Loop over combinations of parameters and generate datasets
    for base_alpha in [0.01, 0.02]: # Overdominance coefficients
        for i in range(0,4):               # relative strength of DMIs
            overdom = [base_alpha, base_alpha]
            epi_vec = [base_alpha*i/2,base_alpha*i,base_alpha*i,base_alpha*i*2]
            
            # Change pattern for overdominance only
            if i==0:
                out_file_suffix = '_g='+str(int(i))+'a='+str(base_alpha)+'_pA0='+str(pA0)+'_pB0='+str(pB0)+'_pC0='+str(pC0)
            else:
                out_file_suffix = '_g='+str(int(i))+'a='+str(epi_vec[1])+'_pA0='+str(pA0)+'_pB0='+str(pB0)+'_pC0='+str(pC0)
            filename_alleles_ABN = out_file_prefix+'allele_ABN'+out_file_suffix+'.csv'
            filename_alleles_ANB = out_file_prefix+'allele_ANB'+out_file_suffix+'.csv'
            
            # Create file if it does not exist yet
            print(Path(filename_alleles_ANB))
            if not Path(filename_alleles_ANB).is_file():
                print('Generating '+out_file_suffix)

                # Simulate for flanking regions
                resall = simulate_populations(haplo_freq,
                                         max_gen,
                                         get_selection_mat(overdom, epi_vec),
                                         mig_vec,
                                         rab,
                                         outtype='alleles', rmin=-5, rmax=0, steps=25)

                # Simulate for neutral between the DMI
                resall_neutb = simulate_populations(haplo_freq, 
                                         max_gen,
                                         get_selection_mat(overdom, epi_vec, neutb=True),
                                         mig_vec,
                                         rab,
                                         neutb=True,
                                         outtype='alleles', npoints=25)


                # Get output dataframes for alleles and haplotypes (+ calculating pairwise LD) 
                # TODO: Move LD function call to simulate_populations
                resall_allele = resall[0]
                resall_haplo = get_ld(resall[1])
                resall_neutb_allele = resall_neutb[0]
                resall_neutb_haplo = get_ld(resall_neutb[1], neutb=True)

                # Outfiles
                resall_allele.to_csv(out_file_prefix+'allele_ABN'+out_file_suffix+'.csv')
                resall_haplo.to_csv(out_file_prefix+'haplotypes_ABN'+out_file_suffix+'.csv')
                
                resall_neutb_allele.to_csv(out_file_prefix+'allele_ANB'+out_file_suffix+'.csv')
                resall_neutb_haplo.to_csv(out_file_prefix+'haplotypes_ANB'+out_file_suffix+'.csv')
            

                # Plotting heterozygosity and equilibrium frequencies
                plot_equilibrium_over_dist(resall_allele, 
                                            markersize=9, fontsize=12, linewidth=1,
                                            allelelist=['A', 'B', 'N'],
                                            yaxis='het',
                                            filename=out_file_prefix+'het_ABN'+out_file_suffix+'.png', skip_show=True)

                plot_equilibrium_over_dist(resall_neutb_allele, 
                                            markersize=9, fontsize=12, linewidth=1,
                                            neutb=True, 
                                            allelelist=['A', 'B', 'N'],
                                            yaxis='het',
                                            filename=out_file_prefix+'het_ANB'+out_file_suffix+'.png', skip_show=True)

                plot_equilibrium_over_dist(resall_allele, 
                                            markersize=9, fontsize=12, linewidth=1,
                                            allelelist=['A', 'B', 'N'],
                                            yaxis='freq',
                                            filename=out_file_prefix+'freq_ABN'+out_file_suffix+'.png', skip_show=True)
                
                plot_equilibrium_over_dist(resall_neutb_allele, 
                                            markersize=9, fontsize=12, linewidth=1,
                                            neutb=True, 
                                            allelelist=['A', 'B', 'N'],
                                            yaxis='freq',
                                            filename=out_file_prefix+'freq_ANB'+out_file_suffix+'.png', skip_show=True)



