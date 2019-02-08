// File of classes for use with src/bgen_prog.cpp
#ifndef CLASS_H
#define CLASS_H

#include <iostream>
#include <string>
#include <vector>

class parameters {
	public :
		std::string bgen_file, bgi_file, pheno_file, out_file, covar_file, incl_rsids_file, incl_sids_file;
		std::string excl_rsids_file, coeffs_file, coeffs2_file, chr, env_file, env_profile_file;
		std::vector< std::string > rsid;
		int chunk_size, missing_code;
		uint32_t start, end;
		bool range, maf_lim, info_lim, mode_gen_pheno, mode_pred_pheno, mode_gen2_pheno, mode_ssv, incl_snps, excl_snps;
		bool geno_check, keep_constant_variants, sim_w_noise, mode_print_keys, flip_high_maf_variants;
		bool select_rsid, mode_compute_correlations, use_raw_covars, use_raw_env, mode_low_mem;
		bool print_causal_rsids;
		double min_maf, min_info, sigma, hb, hg, hb2, hg2, hc, he, rescale_coeffs;
		std::vector < std::string > incl_sample_ids;
		unsigned int random_seed;

	// constructors/destructors
	parameters() {
		bgen_file = "NULL";
		pheno_file = "NULL";
		bgi_file = "NULL";
		env_file = "NULL";
		covar_file = "NULL";
		coeffs_file = "NULL";
		coeffs2_file = "NULL";
		out_file = "NULL";
		incl_sids_file = "NULL";
		incl_rsids_file = "NULL";
		excl_rsids_file = "NULL";
		chunk_size = 256;
		missing_code = -999;
		range = false;
		maf_lim = false;
		info_lim = false;
		mode_ssv = false;
		mode_gen_pheno = false;
		mode_pred_pheno = false;
		mode_gen2_pheno = false;
		mode_print_keys = false;
		print_causal_rsids = false;
		mode_compute_correlations = false;
		flip_high_maf_variants = false;
		incl_snps = false;
		excl_snps = false;
		select_rsid = false;
		keep_constant_variants = false;
		rescale_coeffs = false;
		sim_w_noise = false;
		mode_low_mem = false;
		use_raw_covars = false;
		use_raw_env = false;
		geno_check = true; // check allele probs sum to 1 by default
		sigma = 1;
		hb = 0;  // trait variance explained by additive genetics
		hg = 0;  // trait variance explained by additive GxE
		hb2 = 0;  // trait variance explained by additive genetics
		hg2 = 0;  // trait variance explained by additive GxE
		hc = 0;  // trait variance explained by additive covar
		he = 0;  // trait variance explained by additive env
		random_seed = -1;
	}

	~parameters() {
	}
};

#endif
