# Design and model choices shape inference of age-varying genetic effects on complex traits


## 1. Preprocessing of phenotype data

- [extractPheno.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/extractPheno.R)
- Extract and recode phenotype data
- Generate slopes of change for genome-wide tests on change


## 2. Perform genome-wide scans using REGENIE

### GWA in REGENIE - Step 1:

```
$myprog \
--step 1 \
--bed $HOME/data/ukbb/_001_ukb_cal_allChrs \
--phenoFile $HOME/projects/$project/data/gwa \
--keep $HOME/data/ukbb/qc_pass_MAC100.id \
--extract $HOME/data/ukbb/qc_pass_MAC100.snplist \
--force-qt \
--minINFO 0.9 \
--phenoCol $run \
--covarFile $HOME/projects/$project/data/gwaCovar \
--covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--catCovarList SEX,batch \
--maxCatLevels 3000 \
--bsize 1000 \
--lowmem \
--force-qt \
--lowmem-prefix $HOME/projects/$project/output/gwas/${run} \
--out ${HOME}/projects/$project/output/gwas/${run}_s1 \
--threads $SLURM_CPUS_PER_TASK
```

#### GWA in REGENIE - Step 2:

```
$myprog \
--step 2 \
--pgen $UKBB/pgen/ukb22828_c"$a"_b0_v3 \
--ref-first \
--phenoFile $HOME/projects/$project/data/gwa \
--covarFile $HOME/projects/$project/data/gwaCovar \
--catCovarList SEX,batch \
--maxCatLevels 3000 \
--phenoCol $run \
--approx \
--pred $HOME/projects/$project/output/gwas/${run}_s1_pred.list \
--bsize 400 \
--out $HOME/projects/$project/output/gwas/chr${a}_${run} \
--force-qt \
--threads $SLURM_CPUS_PER_TASK
```

- Process the results using [processGWA.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/processGWA.R)

## 8. Summary of the results

- Script to summarize the results and generate manuscript plots: [summaryLongitudinal.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/summaryLongitudinal.R)

- Script to generate supplement tables: [createTables.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/createTables.R)

