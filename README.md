# bgen_utils

Programme to allow basic simulation of phenotypes with a heritabile component from `bgen` file format. 

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

Basic usage:
```
cmake-build-debug/bgen_utils --sim_pheno \
--bgen unit/data/n50_p100.bgen \
--coeffs unit/data/coeffs_w_snpkey.txt \
--out tmp.txt
```

There are two possible formats allowed for the text file given to --coeffs.
1. Txt file with header "beta" and M rows below each containing a coefficient, where M is the number of SNPs in the file passed to --bgen.
2. Txt file with header "SNPID beta". Any number of rows are allowed, and SNPs whose SNPid matches that given in the --coeffs file will use the corresponding coefficient (and zero effect otherwise).
