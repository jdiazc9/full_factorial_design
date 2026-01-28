# Full factorial construction of synthetic microbial communities

This repository includes all data and code for the paper:

Diaz-Colunga J, Catalan P, San Roman M, Arrabal A, Sanchez A (2024). Full factorial construction of synthetic microbial communities. *eLife* **13**:RP101906. DOI: [10.7554/eLife.101906.1](https://doi.org/10.7554/eLife.101906.1)

The ./protocol directory contains the detailed step-by-step protocol (protocol_8species_25uL.txt) followed for the assembly of colorant combinations and *Pseudomonas* consortia. It also contains an R script (protocol_generator.R) for users to produce tailored step-by-step protocols with customizable parameters.

The ./data directory contains data for the experiment mixing synthetic colorants (colorants.txt), and for the experiment mixing *Pseudomonas* strains (pseudo.txt).

The ./scripts directory contains the R scripts to produce the panels in the main and supplementary figures. Plots are saved under the ./plots directory. Files description:

**colorants.R** produces the panels in Figure 2 & Supplementary Figures 2, 3, and 4.

**pseudo.R** produces the panels in Figures 3 and 4.

**pseudo_intxns.R** produces the panels in Figure 5 (except for panel 5H) and in Supplementary Figure 6.

**pseudo_taylor_vs_fourier.R** produces panel 5H (functional variance by interaction order) and Supplementary Figure 5.

The ./pics directory contains the pictures of the 96-well plates shown in Supplementary Figure 1.

</br></br>

Please direct any questions to J.D.-C. ([juan.diaz\@ipla.csic.es](mailto:juan.diaz@ipla.csic.es))
