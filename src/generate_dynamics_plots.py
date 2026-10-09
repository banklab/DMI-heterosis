# ------------------------------------
# title: generate_dynamics_plots.py
# author: Julio Ayala Lopez
# Created on: 23.10.2025
# Script to generate trajectory plots, for codominant BDMIs and BDMIs+heterosis
# ------------------------------------

from dmi_heterosis_functions import *
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext', fontset='stix')
#rc('text', usetex=True)


## Codominant BDMI + Overdominance
# Decrease in heterozygosity
haplo_freq = [0. , 0.7, 0. , 0. , 0. , 0. , 0.3, 0.]
max_gen=100000*3
overdom=[0.02,0.02]
epi_vec=[0.03,0.06,0.06,0.12]
r = [1e-2,1e-3]
mig_vec=[0,0]


plot_dyn_freqs(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title=("Allele frequency over time for a codominant BDMI and heterosis. \n"r"$\alpha$=0.02; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/allele_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title=("Haplotype frequency over time for a codominant BDMI and heterosis. \n"r"$\alpha$=0.02; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.7.png")


## Increase in neutral heterozygosity
haplo_freq = [0. , 0.3, 0. , 0. , 0. , 0. , 0.7, 0.]
max_gen=100000*3
overdom=[0.02,0.02]
epi_vec=[0.03,0.06,0.06,0.12]
r = [1e-2,1e-3]
mig_vec=[0,0]


plot_dyn_freqs(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	palette = {'B':'#5C43AF','A':'#FE842F', 'N':'#194122', 'A,B':'#FE842F'},
	title=("Allele frequency over time for a codominant BDMI and heterosis. \n"r"$\alpha$=0.02; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/allele_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.3.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title=("Haplotype frequency over time for a codominant BDMI and heterosis. \n"r"$\alpha$=0.02; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.3.png")



plot_dyn_haplo([0. , 0.3, 0.7 , 0. , 0. , 0. , 0., 0.], 
	max_gen, 
	get_selection_mat(overdom, [0.005,0.01,0.01,0.02]), 
	mig_vec, 
	[1e-6,1e-3], 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title=("Haplotype frequency over time for a codominant BDMI and heterosis. \n"r"$\alpha$=0.02; $\gamma$=0.01; $r_{ab}$=$10^-6$; $r_{bn}$=0.001"),
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.01;r_b=1e-6,1e-3;pA0=0.3.png")

#############################################
## Codominant BDMI
haplo_freq = [0. , 0.7, 0. , 0. , 0. , 0. , 0.3, 0.]
max_gen=100000*3
overdom=[0.0,0.0]
epi_vec=[0.01,0.02,0.02,0.04]
r = [1e-2,1e-3]
mig_vec=[0,0]


plot_dyn_freqs(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title=("Allele frequency over time for a codominant BDMI. \n"r"$\alpha$=0; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/allele_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-2,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title=("Haplotype frequency over time for a codominant BDMI. \n"r"$\alpha$=0; $\gamma$=0.06; $r_{ab}$=0.01; $r_{bn}$=0.001"),
	filename="../Fig/haplotype_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-2,1e-3;pA0=0.7.png")



#### With r_ab=1e-6
haplo_freq = [0. , 0.7, 0. , 0. , 0. , 0. , 0.3, 0.]
max_gen=100000*3
overdom=[0.0,0.0]
epi_vec=[0.01,0.02,0.02,0.04]
r = [1e-6,1e-3]
mig_vec=[0,0]


plot_dyn_freqs(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title=("Allele frequency over time for a codominant BDMI. \n"r"$\alpha$=0; $\gamma$=0.06; $r_{ab}$=$10^-6$; $r_{bn}$=0.001"),
	filename="../Fig/allele_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-6,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title=("Haplotype frequency over time for a codominant BDMI. \n"r"$\alpha$=0; $\gamma$=0.06; $r_{ab}$=$10^-6$; $r_{bn}$=0.001"),
	filename="../Fig/haplotype_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-6,1e-3;pA0=0.7.png")



###############################################
## Codominant BDMI + Overdominance haplotype at limit in complete linkage
# Decrease in heterozygosity
haplo_freq = [0.1 , 0.2, 0. , 0. , 0. , 0. , 0.7, 0.]
max_gen=100000*3
overdom=[0.01,0.01]
epi_vec=[0.039/2,0.039,0.039,0.039*2]
r = [0,0]
mig_vec=[0,0]


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r,
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title=("Haplotype frequency over time for a codominant BDMI with overdominance. \n"r"$\alpha$=0.01; $\gamma$=0.039; $r_{ab}$=0; $r_{bn}$=0"),
	filename="../Fig/haplotype_trajectory_codominant_DMIs_a=0.01;g=0.039;r_b=0,0;pAbn0=0.2; paBN0=0.7; pabn0=0.1.png")




###############################################
## Codominant BDMI + Overdominance frequency on extinction
haplo_freq = [0. , 0.7, 0.3 , 0. , 0. , 0. , 0., 0.]
max_gen=100000*3
overdom=[0.02,0.02]
epi_vec=[0.035,0.07,0.07,0.14]
r = [1e-2,0]
mig_vec=[0,0]


plot_dyn_freqs(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title=("Allele frequency over time for a codominant BDMI with overdominance. \n"r"$\alpha$=0.02; $\gamma$=0.07; $r_{ab}$=0.01"),
	filename="../Fig/allele_trajectory_overdominance+DMIs_a=0.02;g=0.07;r_b=1e-2,0;pA0=0.7.png",
	exclude=['N'])
