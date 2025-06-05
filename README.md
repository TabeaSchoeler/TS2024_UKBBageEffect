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
--bed $DATA/ukbb/_001_ukb_cal_allChrs \
--phenoFile $OUT/output/regenie/phenofiles/$run \
--keep $DATA/ukbb/qc_pass_MAC100.id \
--extract $DATA/ukbb/qc_pass_MAC100.snplist \
--force-qt \
--minINFO 0.9 \
--phenoCol $run \
--covarFile $OUT/output/regenie/phenofiles/${run}_Covar \
--catCovarList SEX,batch \
--maxCatLevels 3000 \
--bsize 1000 \
--lowmem \
--lowmem-prefix $OUT/output/regenie/genofiles/${run} \
--out $OUT/output/regenie/genofiles/${run}_s1 \
--threads $SLURM_CPUS_PER_TASK 
```

#### GWA in REGENIE - Step 2:

```
$myprog \
--step 2 \
--pgen $UKBB/genotypes/imp/v3/pgen/ukb22828_c"$a"_b0_v3 \
--ref-first \
--phenoFile $OUT/output/regenie/phenofiles/$run \
--covarFile $OUT/output/regenie/phenofiles/${run}_Covar \
--catCovarList SEX,batch \
--maxCatLevels 3000 \
--phenoCol $run \
--approx \
--pred $OUT/output/regenie/genofiles/${run}_s1_pred.list \
--bsize 400 \
--out $OUT/output/regenie/genofiles/chr${a}_${run} \
--force-qt \
--threads $SLURM_CPUS_PER_TASK \
$extra_args
```

Parameter specifications

- gene-by-age-interaction model (cross-sectional): ``extra_args="--interaction age --no-condtl"``
- phenotype outcome: ``$run`` (e.g, run=bmi_0 for baseline bmi, run=bmi_change for change in bmi)
- REGENIE software: ``myprog=$PROGS/regenie/regenie_v3.2.6.gz_x86_64_Linux``
- chromosome number: ``$a``

#### Process the results 
- [processGWA.R](https://github.com/TabeaSchoeler/TS2024_UKBBageEffect/blob/main/analysis/processGWA.R)

    - generates summary statistic files
    - performs clumping to select LD-independent variants

- [sumGWA.R](https://github.com/TabeaSchoeler/TS2024_UKBBageEffect/blob/main/analysis/sumGWA.R)

    - summarize clumping output for all traits and models
 
## 3. Processing of variant-specific data

- [gene.R](https://github.com/TabeaSchoeler/TS2024_UKBBageEffect/blob/main/analysis/gene.R)

  - Extract defined genetic variants for UKB participants
  - Estimate genetic effects across age-heterogeneous subgroups per variant
 

## 4. Obtained age estimates and age-varying genetic effects for candidate genetic variants
 
- [lme.R](https://github.com/TabeaSchoeler/TS2024_UKBBageEffect/blob/main/analysis/lme.R)


## 5. Summary of the results

- Script to summarize the results and generate manuscript plots: [summaryLongitudinal.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/summaryLongitudinal.R)

- Script to generate supplement tables: [createTables.R](https://github.com/TabeaSchoeler/TS2023_UKBBlongitudinal/blob/main/analysis/createTables.R)

