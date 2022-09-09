# lwf_LM_popgen
 
# How to cite:
Yue Shi, Jared J. Homola, Peter T. Euclide, Daniel A. Isermann, David C. Caroffino, Megan V. McPhee, Wesley A. Larson. High-density genomic data reveal fine-scale population structure and pronounced islands of adaptive divergence in lake whitefish (*Coregonus clupeaformis*) from Lake Michigan. *Evolutionary Applications* (2022.08.26). https://doi.org/10.1111/eva.13475

# Data

Input and intermediate data for various analyses in the study. 

**Note**:

Demultiplexed initial RAD sequencing data (N=30) and Rapture data (N=951) used in this study are archived in the NCBI Sequence Read Archive with a BioProject ID, PRJNA687593. The Rapture data also included additional 175 individuals with mixed stock origins for a different study, and these 175 samples were sequenced and processed together with 951 baseline samples. Sample meta information along with sequence accession numbers for all samples can be found in Table S5. 

A fasta file for the Rapture panel with 100K baited loci and a vcf file after filtering (829 individuals and 197,588 SNPs) are archived on DRYAD [https://doi.org/10.5061/dryad.r4xgxd2gq](https://doi.org/10.5061/dryad.r4xgxd2gq). PLINK files used in various analyses are also archived on DRYAD due to the file size limit of Github.

# Scripts

Scripts are for reference only. Please adjust accordingly for your computing environment and working directory. For certain analyses, scripts are demonstrated using one region-specific dataset as an example. 

`./scripts/0_popgen/`: 
 - 00_fig1A_map.r: Figure 1A.
 - 01_fig1BC_pca.r: Figure 1B & 1C.
 - 02_admixture.r: Figure S1.
 - 03_admixture_east.r: Figure S2.
 - 04_pwfst.r: Figure S3.
 - 05_IBDF_waterDist.R
 - 06_IBD_mantel.r: Figure 2. 
 
`./scripts/1_neuetral_vs_outlier/`: pcadapt analyses
  - 10_pcadapt_SNPthinning.r
  - 11_pcadapt_ol.r
  - 12_pcadapt_ol_summary.r

`./scripts/2_islands/`: 
  - 20_slidingwindow.r
  - 21_islands_ol.r
  - 22_islands_manhattan_plot.r: Figure 3. 
  - 23_region_ol_pwFst.r
  - 24_region_ol_pwFst_heatmap.r: Figure 4.
  - 25_pwFst_all_vs_islands.r: Figure S4.
  - 26_islands_afreqBypop.r: Figure S5.
  - 27_LD_ol_bychr.r: Figure S6.
  - 28_LD_ol_bychr_permutation.r: Table S4. 
  - 29_LDheatmap_bychr.r: Figure S7.
  - 210_lostruct.r
  - 211_chr20_inv.r: Figure 5. 
  
 

  


 


