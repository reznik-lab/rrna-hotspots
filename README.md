# Mitochondrial ribosomal RNA is the target of functionally dominant hotspot mutations in cancer

Clone the directory:
```
git clone https://github.com/reznik-lab/rrna-hotspots.git
```
All packages can be installed from
```
./analysis/prerequisites.R
```
## Whole-Genome-Sequencing Analysis (GEL)
1. mtDNA mutation calling across the GEL cohort was conducted within the genomics england research environment.

2. To calculate hotspot mutations, a master maf file is passed as input. 
```
./analysis/run_hotspots_snv.R
```
3. To produce figures 1 and 2, the code is available here:
```
./analysis/generate_all_figures_01_05.Rmd
```

## Single-Cell Multiome Analysis
1. To call mtDNA mutations, mgatk v.7.0 was used on the 10x output to produce a parseable Seurat object:
```
git clone https://github.com/caleblareau/mgatk
```
Then specific variants were parsed using:
```
./analysis/get_mutations_atac_data.R
```
The remaining analysis to produce figure 3 is available at
```
./analysis/generate_fig3scRNAseqfigs.Rmd
```

## Proteomics and Metabolomics Analysis 
Code to generate figures is available at:
```
./analysis/generate_functional_figs.Rmd
```

