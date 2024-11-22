# Mitochondrial ribosomal RNA is the target of functionally dominant hotspot mutations in cancer

Clone the directory:
```
git clone https://github.com/reznik-lab/rrna-hotspots.git
```
1. All packages can be installed from
```
./analysis/prerequisites.R
```
2. mtDNA mutation calling across the GEL cohort was conducted within the genomics england research environment.

3. To calculate hotspot mutations, a master maf file is passed as input.
```
./analysis/run_hotspots_snv.R
```
4. To produce figures 1 and 2, the code is available here:
```
./analysis/generate_all_figures_01_05.Rmd
```
5. 
