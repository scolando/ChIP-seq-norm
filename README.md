
# ChIP-Seq-norm

Repository with code associated with the paper, "Selecting ChIP-Seq Normalization Methods from the Perspective of their Technical Conditions"

1.  The **simulation_code** folder contains the code to simulate the ChIP-Seq read counts and conduct downstream normalization and differential binding analysis (dba) using MAnorm2 an DiffBind. *A brief description of the files within the simulation_code folder can be found below*

2.  The **sim_visualization_code** folder contains the code to analyze the outputted .csvs from the simulated data, comparing average false discovery rate, power, and absolute size factor ratios across our eight distinct simulation conditions.

3.  The **experiment_code** folder contains the code to conduct PCA on experimental CUT&RUN data, which is from the paper "Genomic Occupancy of the Bromodomain Protein Bdf3 Is Dynamic during Differentiation of African Trypanosomes from Bloodstream to Procyclic Forms" by Ashby et al. (2023). *The Fastq files for RNA-seq used in this paper are available from the NCBI SRA under accession number [PRJNA924324](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA924324).*

4.  The **png_outputs** folder contains the visualizations created by analyzing the simulated and experimental data.


## Files in simulation_code Folder
*A brief description of each the files within the simulation_code folder.*

| File Name                        | Purpose                                                                                                                                                                                                                                                                                                                                                                                                               |
|:---|:---|
| **HPCrun.sh**                    | Shell script to run the HPCrun.R via a slurm job                                                                                                                                                                                                                                                                                                                                                                      |
| **HPCrun.R**                     | Rscript which sets the simulation parameters, executes the purrr parallelization of the big_function.R, and writes the output to a CSV based on the simulations specs                                                                                                                                                                                                                                                 |
| **big_function.R**               | Function which contains the two simulation sub-functions: simulate_reads_proportions.R and read_analysis.R. It also recursively deletes previous simulation files to free up the system's memory.                                                                                                                                                                                                                     |
| **simulate_reads_proportions.R** | Function which creates .bam and .bam.bai files based on simulation parameters, as well as preforms peak-calling via [MACS2](https://pypi.org/project/MACS2/) (by calling macs2.sh).                                                                                                                                                                                                                                   |
| **addBackground_new.R**          | Function called within the simulate_reads_proportions.R function to add background binding to each sample replicate. This function is an adapted version of the addBackground function in the [xscss package](https://bioinf.wehi.edu.au/csaw/), which accompanies the [csaw package](https://www.bioconductor.org/packages/release/bioc/html/csaw.html).                                                             |
| **macs2.sh**                     | Shell script which performs MACS2 peak-calling generating the relevant .narrowpeak files for each replicate.                                                                                                                                                                                                                                                                                                          |
| **read_analysis.R**              | Function which analyzes the files outputted from simulate_reads_proportions.R (i.e., .bam, .bam.bai, and .narrowpeak files for each replicate), using the [DiffBind package](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to identify consensus peaks, conduct between-sample normalization (using dba.normalize.sim.R), and differential binding analysis (using the method native to DESeq2). |                                                                                                                                                                                                                                                                                                                                                                                                                     |To perform MAnorm2 normalization and differential binding analysis, the function calls MAnorm2.sh and MAnorm2.sim.R since MAnorm2 is not available in DiffBind.                                                                                                                                                                                                                                                       |
| **dba.normalize.sim.R**          | Wrapper function for the dba.normalize function in the DiffBind package, which allows the normalization method's name to be input as a string to conduct the proper between-sample normalization. Supports the seven between-sample normalization methods outlined on pg. 46 of DiffBind's vignette (Section 7.5).                                                                                                    |
| **MAnorm2.sh**                   | Shell script which generates the .csv needed for MAnorm2 normalization and differential binding analysis using the [MAnorm2_utils package](https://github.com/tushiqi/MAnorm2_utils).                                                                                                                                                                                                                                 |
| **MAnorm2.sim.R**                | Function which performs MANorm2 within and between-sample normalization, as well as differential binding analysis as suggested in the [MAnorm2 vignette](https://cran.r-project.org/web/packages/MAnorm2/vignettes/MAnorm2_vignette.html). Adds the MAnorm2 results to the dataset which contains the results from all other between-sample normalization methods tested.                                             |

## Required Packages

**Bioconductor Packages**

-   [DiffBind](https://www.bioconductor.org/packages/release/bioc/html/DiffBind.html), [csaw](https://www.bioconductor.org/packages/release/bioc/html/csaw.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [GreyListChIP](https://www.bioconductor.org/packages/release/bioc/html/GreyListChIP.html)

**CRAN Packages**

-   [MAnorm2](https://cran.r-project.org/web/packages/MAnorm2/index.html), [Tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html), [psych](https://cran.r-project.org/web/packages/psych/index.html), [truncnorm](https://cran.r-project.org/web/packages/truncnorm/index.html), [fuzzyjoin](https://cran.r-project.org/web/packages/fuzzyjoin/index.html), [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html), [UpSetR](https://cran.r-project.org/web/packages/UpSetR/index.html), [broom](https://cran.r-project.org/web/packages/broom/index.html), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)

**Other Packages**

-   [MACS2](https://pypi.org/project/MACS2/), [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html), [MAnorm2_utils](https://pypi.org/project/MAnorm2-utils/)
