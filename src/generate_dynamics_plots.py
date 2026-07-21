# ------------------------------------
# title: generate_dynamics_plots.py
# author: Julio Ayala Lopez
# Created on: 23.10.2025
# Script to generate trajectory plots, for codominant DMIs and DMIs+heterosis
# ------------------------------------

from dmi_heterosis_functions import *
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)


## Codominant DMI + Overdominance
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
	title="Allele frequency over time for a codominant DMI and heterosis. \na=0.02; g=0.06; r_ab=0.01; r_bn=0.001",
	filename="../Fig/allele_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	title="Haplotype frequency over time for a codominant DMI and heterosis. \na=0.02; g=0.06; r_ab=0.01; r_bn=0.001",
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.7.png")


## Increase in heterozygosity
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
	title="Allele frequency over time for a codominant DMI and heterosis. \na=0.02; g=0.06; r_ab=0.01; r_bn=0.001",
	filename="../Fig/allele_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.3.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title="Haplotype frequency over time for a codominant DMI and heterosis. \na=0.02; g=0.06; r_ab=0.01; r_bn=0.001",
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.06;r_b=1e-2,1e-3;pA0=0.3.png")



plot_dyn_haplo([0. , 0.3, 0.7 , 0. , 0. , 0. , 0., 0.], 
	max_gen, 
	get_selection_mat(overdom, [0.005,0.01,0.01,0.02]), 
	mig_vec, 
	[1e-6,1e-3], 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title="Haplotype frequency over time for a codominant DMI and heterosis. \na=0.02; g=0.01; r_ab=1e-6; r_bn=0.001",
	filename="../Fig/haplotype_trajectory_overdominance+DMIs_a=0.02;g=0.01;r_b=1e-6,1e-3;pA0=0.3.png")

#############################################
## Codominant DMI
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
	title="Allele frequency over time for a codominant DMI. \na=0; g=0.06; r_ab=0.01; r_bn=0.001",
	filename="../Fig/allele_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-2,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title="Haplotype frequency over time for a codominant DMI. \na=0; g=0.06; r_ab=0.01; r_bn=0.001",
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
	title="Allele frequency over time for a codominant DMI. \na=0; g=0.06; r_ab=1e-6; r_bn=0.001",
	filename="../Fig/allele_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-6,1e-3;pA0=0.7.png")


plot_dyn_haplo(haplo_freq, 
	max_gen, 
	get_selection_mat(overdom, epi_vec), 
	mig_vec, 
	r, 
	haplotypes = ["abn", "Abn", "aBn", "ABn", "abN", "AbN", "aBN", "ABN"],
	title="Haplotype frequency over time for a codominant DMI. \na=0; g=0.06; r_ab=1e-6; r_bn=0.001",
	filename="../Fig/haplotype_trajectory_codominant_DMIs_a=0;g=0.02;r_b=1e-6,1e-3;pA0=0.7.png")



###############################################
## Codominant DMI + Overdominance haplotype at limit in complete linkage
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
	title="Haplotype frequency over time for a codominant DMI with overdominance. \na=0.01; g=0.039; r_ab=0; r_bn=0",
	filename="../Fig/haplotype_trajectory_codominant_DMIs_a=0.01;g=0.039;r_b=0,0;pAbn0=0.2; paBN0=0.7; pabn0=0.1.png")




