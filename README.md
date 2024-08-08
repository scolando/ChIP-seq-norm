# ChIP-Seq-norm

Repository with code associated with the paper, "Selecting ChIP-Seq Normalization Methods from the Perspective of their Technical Conditions"

1.  The **simulation_code** folder contains the code to simulate the ChIP-Seq read counts and conduct downstream normalization and differential binding analysis (dba) using MAnorm2 an DiffBind. *A brief description of each file within the simulation_code folder can be found below*

2.  The **sim_visualization_code** folder contains the code to analyze the outputted .csvs from the simulated data, comparing average false discovery rate, power, and absolute size factor ratios across the eight distinct simulation conditions.

3.  The **experiment_code** folder contains the code to conduct PCA and K-means clustering on actual CutandRun data, which is from the paper "Genomic Occupancy of the Bromodomain Protein Bdf3 Is Dynamic during Differentiation of African Trypanosomes from Bloodstream to Procyclic Forms" by Ashby et al. (2023). *The Fastq files for RNA-seq used in this paper are available from the NCBI SRA under accession number [PRJNA924324](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA924324).*

4.  The **png_outputs** folder contains the visualizations created by analyzing the simulated and experimental data.

## Simulation_code Files

*A brief description of each the files within the simulation_code folder.*

1.  **HPCrun.sh** is the shell script to run the HPCrun.R via slurm 

2.  **HPCrun.R**

3.  **big_function.R**

4.  **simulate_reads_proportions.R**

5.  **addBackground_new.R**

6.  **macs2.sh**

7.  **read_analysis.R**

8.  **dba.normalize.sim.R**

9.  **MAnorm2.sh**

10. **MAnorm2.sim.R**
