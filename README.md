# bgen_utils

Programme to allow basic simulation of phenotypes with a heritable component from `bgen` file format. Available for beta testing.

## Install
Dependencies:
- installed version of the BGEN library
- boost (specifically, boost_iostreams to allow read/write to gzipped files)

In the CMakeLists.txt file I currently set paths to these dependencies in lines 9,17,25 (BGEN) and lines 13,21,29 (boost_iostreams).

After editing CMakeLists.txt you can compile with:
```
mkdir bin
cd bin
cmake ..
cd ..
cmake --build bin -- -j 4
```

## Phenotype simulation
Basic usage:
```
bin/bgen_utils --sim_pheno \
--bgen unit/data/n50_p100.bgen \
--coeffs unit/data/coeffs_w_snpid.txt \
--out tmp.txt
```

Here we simulate
```
Y = X \beta + \epsilon
```
where columns of the dosage matrix X are normalised to have mean zero and variance one.

There are two possible formats allowed for the text file given to --coeffs.
1. Txt file with header "beta" and M rows below each containing a coefficient, where M is the number of SNPs in the file passed to --bgen.
2. Txt file with header "SNPID beta". Any number of rows are allowed, and SNPs whose SNPID matches that given in the --coeffs file will use the corresponding coefficient (and zero effect otherwise).

To simulate a phenotype with multiplicative GxE effects with a single environmental variable use:
```
bin/bgen_utils --sim_pheno \
--bgen unit/data/n50_p100.bgen \
--coeffs unit/data/coeffs_w_snpid_w_gxe_effects.txt \
--environment unit/data/n50_p100_single_env.txt \
--out tmp.txt
```

This simulates according to the model
```
Y = X \beta + E \odot X \gamma + \epsilon
```
where `E` is `Nx1` and `\gamma` is `Mx1`.

To simulate a phenotype with multiplicative GxE effects with multiple environmental variables use:
```
bin/bgen_utils --sim_pheno \
--bgen unit/data/n50_p100.bgen \
--coeffs unit/data/coeffs_w_snpid_w_gxe_effects.txt \
--environment unit/data/n50_p100_env.txt \
--environment_weights unit/data/n50_p100_env_weights.txt \
--out tmp.txt
```

This simulates according to the model
```
Y = X \beta + (Ew) \odot X \gamma + \epsilon
```
where `E` is `NxL` and `w` is `Lx1`.

Finally, we can also simulate GxE effects with multiple gxe latent components using
```
bin/bgen_utils --sim_pheno \
--bgen unit/data/n50_p100.bgen \
--coeffs unit/data/coeffs_w_snpid_w_gxe_2comp.txt \
--environment unit/data/n50_p100_env.txt \
--environment_weights unit/data/n50_p100_env_weights_2comp.txt \
--out tmp.txt
```

This simulates according to the model
```
Y = X \beta + rowSum{ (E w) \odot (X \gamma) } + \epsilon
```
where `E` is `NxL`, `w` is `LxK` and `gamma` is `MxK`. Using eigen notation this is expressed as
```
Y = X_ij \beta_j + X_ij E_il w_lk \gamma_jk + \epsilon.
```

## Compute SNP-Env correlations

Given a BGEN file containing the `NxM` dosage matrix `X` and `NxL` matrix of environmental variables `E`, this programme computes the following correlations
```
\sum_i \sum_{l_1, l_2} X_{ij}^2 E_{il_1} E_{il_2}.
```
This output is used by LEMMA. Both `X` and `E` are by default normalised to have column mean zero and column variance one.

Basic usage:
```
bin/bgen_utils \
--compute-env-snp-correlations --mode_low_mem \
--incl_sample_ids unit/data/sample_ids.txt \
--bgen unit/data/n50_p100.bgen \
--environment unit/data/n50_p100_env.txt \
--out n50_p100_dxteex.txt.gz
```
Note that the sample filter above is optional.

This computation is computationally intensive. To compute in parallel the following can be used:
```
while read -r n chr st en; do \
  bin/bgen_utils \
  --compute-env-snp-correlations --mode_low_mem \
  --incl_sample_ids unit/data/sample_ids.txt \
  --range ${chr}:${st}-${en} \
  --bgen unit/data/n50_p100.bgen \
  --environment unit/data/n50_p100_env.txt \
  --out n50_p100_dxteex_chunk${n}.txt.gz;
done < unit/data/n50_p100_snp_chunks.txt

cat n50_p100_dxteex_chunk*.txt.gz > n50_p100_dxteex_stitched.txt.gz
zdiff n50_p100_dxteex_stitched.txt.gz n50_p100_dxteex.txt.gz
```
