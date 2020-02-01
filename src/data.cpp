//
// Created by kerin on 2019-02-26.
//

#include "parameters.hpp"
#include "utils.hpp"
#include "data.hpp"
#include "file_utils.hpp"
#include "eigen_utils.hpp"

#include "tools/eigen3.3/Dense"
#include "bgen_parser.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <regex>
#include <string>
#include <set>
#include <stdexcept>

namespace boost_io = boost::iostreams;

Data::Data(const parameters &my_params) : params(my_params){
	bgenView = genfile::bgen::View::create(my_params.bgen_file);
	bgen_pass = true;
	n_samples = bgenView->number_of_samples();
	n_var_parsed = 0;
	n_var = 0;
	n_constant_variance = 0;
	n_covar = 0;
	n_env = 0;
	n_gxe_components = 0;
	n_gxe_components2 = 0;

	// system time at start
	start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout << "Starting analysis at " << std::ctime(&start_time) << std::endl;
}

Data::~Data() {
	// system time at end
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << std::endl << "Analysis finished at " << std::ctime(&end_time);
	std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

	boost_io::close(outf);
}

void Data::gen_genetic_effects() {

	// Step 1; Read in raw covariates
	if(params.covar_file != "NULL") {
		read_covar();
	}
	if(params.env_file != "NULL") {
		read_environment();
	}
	read_coeffs();
	assert(n_gxe_components <= n_env);
	if(n_env > 1 && params.env_profile_file != "NULL") {
		read_env_profile();
		assert(env_profile.cols() == n_gxe_components);
	} else if (n_env > 1) {
		throw std::runtime_error("--environment_weights must be provided to "
		                         "generate a GxE effect with multiple environments.");
	} else if (n_env == 1) {
		env_profile = EigenDataMatrix::Constant(1, 1, 1);
	}

	if(params.coeffs2_file != "NULL") {
		read_coeffs2();
		assert(n_gxe_components2 <= n_env);
		if(n_env > 1 && params.env_profile_file2 != "NULL") {
			read_env_profile();
			assert(env_profile2.cols() == n_gxe_components2);
		} else if (n_env > 1) {
			assert(B2.cols() == B.cols());
			env_profile2 = env_profile;
		} else if (n_env == 1) {
			env_profile = EigenDataMatrix::Constant(1, 1, 1);
		}
	}


	// Step 2; Reduce raw covariates and phenotypes to complete cases
	// - may change value of n_samples
	// - will also skip these cases when reading bgen later
	reduce_to_complete_cases();

	// Step 3; Normalise covars
	if(n_covar > 0 && !params.use_raw_covars) {
		EigenUtils::center_matrix(W, n_covar);
		EigenUtils::scale_matrix(W, n_covar, covar_names);
	}
	if(n_env > 0 && !params.use_raw_env) {
		EigenUtils::center_matrix(E, n_env);
		EigenUtils::scale_matrix(E, n_env, env_names);
	}
	if(n_env > 0) {
		// Problem if scaling matrix removes columns from E
		eta = E * env_profile;
	}

	// Output to user
	std::cout << std::endl << "Generating predicted phenotype" << std::endl;
	if (params.normalise_genotypes) {
		std::cout << "- using normalised expected dosage (mean zero, variance one)" << std::endl;
	} else {
		std::cout << "- using raw expected dosage"<< std::endl;
	}
	if(params.mode_low_mem) {
		std::cout << "- simulating low memory compression used by LEMMA" << std::endl;
	}
	if(n_env > 0) {
		if (params.use_raw_env) {
			std::cout << "- env variables unprocessed" << std::endl;
		} else {
			std::cout << "- env variables normalised to mean zero, variance one" << std::endl;
		}
	}
	std::cout << std::endl;

	// S2; genetic or interaction effects
	Xb  = EigenDataArrayX::Zero(n_samples);
	Xb2 = EigenDataArrayX::Zero(n_samples);
	Xg  = EigenDataArrayXX::Zero(n_samples, n_gxe_components);
	Xg2 = EigenDataArrayXX::Zero(n_samples, n_gxe_components2);
	int ch = 0;
	long coeff_index, n_matched = 0;
	while (read_bgen_chunk()) {
		if(ch % 10 == 0) {
			std::cout << "Chunk " << ch+1 << " read (size " << n_var_chunk;
			std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
			std::cout << " variants parsed)" << std::endl;
		}

		if(B.rows() < bgenView->number_of_variants()) {
			assert(B.rows() >= n_var + n_var_chunk);
		}

		// Add effects to Y
		for (int kk = 0; kk < n_var_chunk; kk++) {
			if(B.rows() < bgenView->number_of_variants()) {
				coeff_index = kk + n_var;
			} else {
				coeff_index = valid_var_index_chunk[kk];
			}

			Xb += G.col(kk).array() * B(coeff_index, 0);
			if(n_env > 0) {
//				Xg += G.col(kk).array() * B.block(coeff_index, 1, 1, n_gxe_components);
				for (long cc = 0; cc < n_gxe_components; cc++) {
					Xg.col(cc) += B(coeff_index, cc + 1) * G.col(kk).array();
				}
			}

			if(params.coeffs2_file != "NULL") {
				Xb2 += G.col(kk).array() * B2(coeff_index, 0);
				if(n_env > 0) {
					for (long cc = 0; cc < n_gxe_components2; cc++) {
						Xg2.col(cc) += B2(coeff_index, cc + 1) * G.col(kk).array();
					}
				}
			}
		}
		n_var += n_var_chunk;
		ch++;
	}

	// Error checking
	if(B.rows() != n_var_parsed && B.rows() != n_var) {
		std::cout << "ERROR: n var read in = " << n_var << std::endl;
		std::cout << "ERROR: n coeffs read in = " << B.rows() << std::endl;
		assert(n_var == B.rows() || n_var_parsed == B.rows());
	}

	// Output
	std::cout << "Used " << n_var << " SNPs to construct genetic effects" << std::endl;
	if(n_var < n_var_parsed) {
		std::cout << "- Removed " << n_constant_variance  << " SNP(s) with zero variance" << std::endl;
		if(params.maf_lim || params.info_lim) {
			std::cout << "- Removed " << n_var_parsed - n_var - n_constant_variance;
			std::cout << " SNP(s) due to MAF or INFO limits" << std::endl;
		}
	}

	// Subset to coefficients corresponding to valid SNPs
	if(n_var < n_var_parsed && B.rows() == n_var_parsed) {
		B = EigenUtils::subset_rows(B, valid_var_index);
		if(params.coeffs2_file != "NULL") {
			B2 = EigenUtils::subset_rows(B2, valid_var_index);
		}
	}

	Y = Xb + Xb2;
	Zg = EigenDataArrayX::Zero(n_samples);
	Zg2 = EigenDataArrayX::Zero(n_samples);
	if(n_env > 0) {
		Zg = (Xg * eta.array()).rowwise().sum();
		if(params.coeffs2_file != "NULL") {
			Zg2 = (Xg2 * eta.array()).rowwise().sum();
		}
		Y += Zg + Zg2;
	}
}

void Data::pred_pheno() {
	gen_genetic_effects();
	Wtau = EigenDataVector::Zero(n_samples);
	Ealpha = EigenDataVector::Zero(n_samples);
	Y = Wtau + Ealpha + Xb + Zg + Xb2 + Zg2;
}

void Data::sim_pheno() {

	// Get predicted effects
	gen_genetic_effects();
	Wtau = EigenDataVector::Zero(n_samples);
	Ealpha = EigenDataVector::Zero(n_samples);

	// Generate random effects
	// http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
	if(params.random_seed == -1) {
		std::random_device rd;
		params.random_seed = rd();
	}
	std::cout << std::endl << "Initialising random sample generator with seed " << params.random_seed << std::endl;
	std::mt19937 generator{params.random_seed};

	// additive noise
	noise.resize(n_samples);
	sim_gaussian_noise(noise, 1, generator);
	noise *= std::sqrt(params.sigma / var(noise));

	// target variances
	double tol = 0.000001;
	if(params.rescale_coeffs) {
		std::cout << "Rescaling components to ensure that sample ";
		std::cout << "heritability matches expected heritability" << std::endl;
		double resid_pve = 1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2;

		// additive covar effects
		if(params.hc > tol) {
			EigenDataVector tau(n_covar);
			sim_gaussian_noise(tau, 1, generator);
			Wtau = W * tau;
			double sf = std::sqrt(params.sigma * params.hc / resid_pve / var(Wtau));
			Wtau     *= sf;
		}

		// additive env effects
		if(params.he > tol) {
			EigenDataVector alpha(n_env);
			sim_gaussian_noise(alpha, 1, generator);
			Ealpha = E * alpha;
			double sf  = std::sqrt(params.sigma * params.he / resid_pve / var(Ealpha));
			Ealpha    *= sf;
		}

		// rescaling for correct heritability
		scalarData var_xb = var(Xb);
		Xb       *= std::sqrt(params.sigma * params.hb / resid_pve / var_xb);
		B.col(0) *= std::sqrt(params.sigma * params.hb / resid_pve / var_xb);

		if(n_env > 0) {
			scalarData var_zg = var(Zg);
			double sf = std::sqrt(params.sigma * params.hg / resid_pve / var_zg);
			Zg       *= sf;
			Xg       *= sf;
			B.block(0, 1, B.rows(), n_gxe_components) *= sf;
		}

		if(params.coeffs2_file != "NULL") {
			scalarData var_xb = var(Xb2);
			double sf  = std::sqrt(params.sigma * params.hb2 / resid_pve / var_xb);
			Xb2       *= sf;
			B2.col(0) *= sf;

			if(n_env > 0) {
				scalarData var_zg = var(Zg2);
				sf         = std::sqrt(params.sigma * params.hg2 / resid_pve / var_xb);
				Zg2       *= sf;
				Xg2       *= sf;
				B2.block(0, 1, B2.rows(), n_gxe_components2) *= sf;
			}
		}
	}

	Y = Wtau + Ealpha + Xb + Zg + Xb2 + Zg2 + noise;
	std::cout << "Empirical PVE from each effect [computed as Var(<effect>) / Var(y)]:" << std::endl;
	if(n_env > 0 && params.he > tol)   std::cout << "- PVE-Env = "    << var(Ealpha) / var(Y) << std::endl;
	if(n_covar > 0 && params.hc > tol) std::cout << "- PVE-Covars = " << var(Wtau) / var(Y) << std::endl;
	std::cout << "- PVE-G = " << var(Xb) / var(Y) << std::endl;
	if(n_env > 0)                      std::cout << "- PVE-GxE = "    << var(Zg) / var(Y) << std::endl;
	if(params.coeffs2_file != "NULL") {
		                               std::cout << "- PVE-G2 = "     << var(Xb2) / var(Y) << std::endl;
		if(n_env > 0)                  std::cout << "- PVE-GxE2 = "   << var(Zg2) / var(Y) << std::endl;
	}
	std::cout << std::endl;
}

void Data::calc_ssv() {
	// Returns sum of sample variances from columns of X and Z

	// Step 1; Read in raw covariates
	// - also makes a record of missing values
	if(params.covar_file != "NULL") {
		read_covar();
	}

	// Step 2; Reduce raw covariates and phenotypes to complete cases
	// - may change value of n_samples
	// - will also skip these cases when reading bgen later
	reduce_to_complete_cases();

	// Step 3; Normalise covars
	if(n_covar > 0) {
		EigenUtils::center_matrix( W, n_covar );
		EigenUtils::scale_matrix( W, n_covar, covar_names );
	}

	int ch = 0;
	double tmp_x, tmp_z, mean_z;
	s_x = 0.0, s_z = 0.0;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var_chunk;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;
		assert(n_var_chunk == G.cols());
		n_var += G.cols();

		// Compute sum of column variances
		for (int kk = 0; kk < n_var_chunk; kk++) {
			tmp_x = G.col(kk).dot(G.col(kk));
			s_x += tmp_x / ((double) (n_samples - 1.0));

			if(n_covar > 0) {
				mean_z = G.col(kk).dot(W.col(0)) / (double) n_samples;
				tmp_z = 0.0;
				for (std::size_t ii = 0; ii < n_samples; ii++) {
					tmp_z += (W(ii,0)*G(ii,kk) - mean_z)*(W(ii,0)*G(ii,kk) - mean_z);
				}
				s_z += tmp_z / ((double) (n_samples - 1.0));
			}
		}
		ch++;
	}

	// Output
	std::cout << "Used " << n_var << " SNPs to compute sum of column variances" << std::endl;
	if(n_var < n_var_parsed) {
		std::cout << "- Removed " << n_constant_variance  << " SNP(s) with zero variance" << std::endl;
		if(params.maf_lim || params.info_lim) {
			std::cout << "- Removed " << n_var_parsed - n_var - n_constant_variance;
			std::cout << " SNP(s) due to MAF or INFO limits" << std::endl;
		}
	}
}

void Data::compute_correlations() {
	// Gen vector of fitted effects Y_hat = age + X beta + Z gamma
	// Gaussian noise added when writing to file

	// Step 1; Read in raw covariates
	// - also makes a record of missing values
	read_environment();
	if(params.covar_file != "NULL") {
		read_covar();
	}

	// Step 2; Reduce raw covariates and phenotypes to complete cases
	// - may change value of n_samples
	// - will also skip these cases when reading bgen later
	reduce_to_complete_cases();

	// Step 3; Normalise covars
	EigenUtils::center_matrix( E, n_env );
	EigenUtils::scale_matrix( E, n_env, env_names );

	if(n_covar > 0) {
		EigenUtils::center_matrix(W, n_covar);
		EigenUtils::scale_matrix(W, n_covar, covar_names);
	}
	if(n_covar > 0) {
		regress_first_mat_from_second(W, "covars", covar_names, E, "env");
	}

	// Output to user
	std::cout << std::endl << "Generating SNP-Env correlations" << std::endl;
	if (params.normalise_genotypes) {
		std::cout << "- using normalised expected dosage (mean zero, variance one)" << std::endl;
	} else {
		std::cout << "- using raw expected dosage"<< std::endl;
	}
	if(params.mode_low_mem) {
		std::cout << "- simulating low memory compression used by LEMMA" << std::endl;
	}
	if(n_env > 0) {
		if (params.use_raw_env) {
			std::cout << "- env variables unprocessed" << std::endl;
		} else {
			std::cout << "- env variables normalised to mean zero, variance one" << std::endl;
		}
	}
	std::cout << std::endl;

	// Step 4; compute correlations
	int ch = 0;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var_chunk;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;

		assert(n_var_chunk == G.cols());

		EigenDataArrayXX dXtEEX_chunk(n_var_chunk, n_env * n_env);
		compute_correlations_chunk(dXtEEX_chunk);
		output_correlations(dXtEEX_chunk);

		n_var += n_var_chunk;
		ch++;
	}

	// Output
	std::cout << "Used " << n_var << " SNPs to compute env-snp correlations" << std::endl;
	if(n_var < n_var_parsed) {
		std::cout << "- Removed " << n_constant_variance  << " SNP(s) with zero variance" << std::endl;
		if(params.maf_lim || params.info_lim) {
			std::cout << "- Removed " << n_var_parsed - n_var - n_constant_variance;
			std::cout << " SNP(s) due to MAF or INFO limits" << std::endl;
		}
	}
}

void Data::print_keys() {
	// Step 4; compute correlations
	int ch = 0;
	reduce_to_complete_cases();
	// std::cout << "num samples: " << n_samples << std::endl;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var_chunk;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;

		for (std::size_t ii = 0; ii < n_var_chunk; ii++) {
			outf << SNPID_chunk[ii] << " " << chromosome_chunk[ii] << " " << rsid_chunk[ii];
			outf << " " << position_chunk[ii];
			outf << " " << alleles_chunk[ii][0] << " " << alleles_chunk[ii][1];
			outf << " " << maf_chunk[ii] << " " << info_chunk[ii] << std::endl;
		}

		n_var += n_var_chunk;
		ch++;
	}

	// Output
	std::cout << "Used " << n_var << " SNPs" << std::endl;
	if(n_var < n_var_parsed) {
		std::cout << "- Removed " << n_constant_variance  << " SNP(s) with zero variance" << std::endl;
		if(params.maf_lim || params.info_lim) {
			std::cout << "- Removed " << n_var_parsed - n_var - n_constant_variance;
			std::cout << " SNP(s) due to MAF or INFO limits" << std::endl;
		}
	}
}

void Data::output_init() {
	std::string ofile = fileUtils::fstream_init(outf, params.out_file, "", "");


	// Output header for ssv file
	if(params.mode_ssv) {
		std::cout << "Writing ssv stats to " << ofile << std::endl;
		outf << "s_x\ts_z\tn\tp" << std::endl;
	} else if (params.mode_compute_correlations) {
		std::cout << "Writing snp-environment correlations to " << ofile << std::endl;
	} else if (params.mode_print_keys) {
		std::cout << "Writing snp-keys to " << ofile << std::endl;
		outf << "SNPID chr rsid pos a0 a1 maf info" << std::endl;
	}
}

void Data::output_results() {
	if(params.mode_ssv) {
		outf << s_x << "\t" << s_z << "\t" << n_samples << "\t";
		outf << n_var << std::endl;
	} else if (params.mode_pred_pheno || params.mode_gen_pheno) {
		{
			std::string ofile = fileUtils::filepath_format(params.out_file, "", "");
			std::cout << "Writing predicted phenotype to " << ofile << std::endl;
			std::vector<std::string> header = {"y"};

			fileUtils::write_matrix(Y, ofile, header);
		}

		{
			std::string ofile   = fileUtils::filepath_format(params.out_file, "", "_predicted_effects");
			std::cout << "Writing predicted effects to " << ofile << std::endl;
			std::vector<std::string> header;
			Eigen::MatrixXd mat(n_samples, 4 + 2 * n_gxe_components);
			long cc = 0;
			if(n_covar > 0) {
				header.emplace_back("Ctau");
				mat.col(cc) = Wtau;
				cc++;
			}
			if(n_env > 0) {
				header.emplace_back("Ealpha");
				mat.col(cc) = Ealpha;
				cc++;
			}
			{
				header.emplace_back("Xbeta");
				mat.col(cc) = Xb;
				cc++;
			}
			if(n_env > 0) {
				for (long ii = 0; ii < n_gxe_components; ii++) {
					header.emplace_back("eta" + std::to_string(ii));
					mat.col(cc + 2 * ii) = eta.col(ii);
					header.emplace_back("Xgamma" + std::to_string(ii));
					mat.col(cc + 2 * ii + 1) = Xg.col(ii);
				}
				cc += 2 * n_gxe_components;
				header.emplace_back("Zgamma");
				mat.col(cc) = Zg;
				cc++;
			}
			mat.conservativeResize(mat.rows(), header.size());
			fileUtils::write_matrix(mat, ofile, header);
		}

		{
			std::string ofile   = fileUtils::filepath_format(params.out_file, "", "_true_rescaled_coeffs");
			std::cout << "Writing rescaled coeffs to " << ofile << std::endl;
			std::vector<std::string> header = {"SNPID", "chr", "rsid", "pos",
				                               "a0", "a1", "af", "beta"};
			for (long cc = 0; cc < n_gxe_components; cc++) {
				header.emplace_back("gamma" + to_string(cc));
			}
			if(params.coeffs2_file != "NULL") {
				header.emplace_back("beta_v2");
				for (long cc = 0; cc < n_gxe_components2; cc++) {
					header.emplace_back("gamma" + to_string(cc) + "_v2");
				}
			}

			dump_coeffs(ofile, header, B, B2);
		}
	}
}

bool Data::read_bgen_chunk() {
	// Wrapper around BgenView to read in a 'chunk' of data. Remembers
	// if last call hit the EOF, and returns false if so.
	// Assumed that:
	// - commandline args parsed and passed to params
	// - bgenView initialised with correct filename
	// - scale + centering happening internally

	// Exit function if last call hit EOF.
	if (!bgen_pass) return false;

	// Temporary variables to store info from read_variant()
	uint32_t pos_j;
	std::string rsid_j, chr_j, SNPID_j;
	std::vector< std::string > alleles_j;

	long nInvalid = sample_is_invalid.size() - n_samples;
	DosageSetter setter_v2(sample_is_invalid, nInvalid);

	double chunk_missingness = 0;
	long n_var_incomplete = 0;

	// Wipe variant context from last chunk
	maf_chunk.clear();
	info_chunk.clear();
	rsid_chunk.clear();
	chromosome_chunk.clear();
	position_chunk.clear();
	alleles_chunk.clear();
	SNPKEY_chunk.clear();
	SNPID_chunk.clear();
	valid_var_index_chunk.clear();

	// Resize genotype matrix
	G.resize(n_samples, params.chunk_size);

	long jj = 0;
	valid_var_index_chunk.clear();
	while ( jj < params.chunk_size && bgen_pass ) {
		bgen_pass = bgenView->read_variant(&SNPID_j, &rsid_j, &chr_j, &pos_j, &alleles_j);
		if (!bgen_pass) break;
		n_var_parsed++;

		// Read probs + check maf filter
		bgenView->read_genotype_data_block( setter_v2 );

		double d1     = setter_v2.m_sum_eij;
		double maf_j  = setter_v2.m_maf;
		double info_j = setter_v2.m_info;
		double mu     = setter_v2.m_mean;
		double missingness_j    = setter_v2.m_missingness;
		double sigma = std::sqrt(setter_v2.m_sigma2);
		long valid_count = n_samples;
		EigenDataArrayX dosage_j = setter_v2.m_dosage;

		// Filters
		if (params.maf_lim && (maf_j < params.min_maf || maf_j > 1 - params.min_maf)) {
			continue;
		}
		if (params.info_lim && info_j < params.min_info) {
			continue;
		}
		if(!params.keep_constant_variants && d1 < 5.0) {
			n_constant_variance++;
			continue;
		}
		if(!params.keep_constant_variants && sigma <= 1e-12) {
			n_constant_variance++;
			continue;
		}
		if(!std::isfinite(mu) || !std::isfinite(sigma)) {
			std::cout << "WARNING: non-finite mean / variance detected for SNP with:" <<std::endl;
			std::cout << "SNPID: " << SNPID_j << std::endl;
			std::cout << "RSID: " << rsid_j << std::endl;
		}

		// filters passed; write contextual info
		chunk_missingness += missingness_j;
		if(missingness_j > 0) n_var_incomplete++;

		maf_chunk.push_back(maf_j);
		info_chunk.push_back(info_j);
		rsid_chunk.push_back(rsid_j);
		chromosome_chunk.push_back(chr_j);
		position_chunk.push_back(pos_j);
		alleles_chunk.push_back(alleles_j);
		SNPID_chunk.push_back(SNPID_j);

		std::string key_j = gen_snpkey(chr_j, pos_j, alleles_j);
		SNPKEY_chunk.push_back(key_j);
		if(params.mode_print_keys) {
			jj++;
			continue;
		}

		if(params.mode_low_mem) {
			double L = 2.0;
			double intervalWidth = L / 256.0;
			double invIntervalWidth = 256.0 / L;

			// Compress
			Eigen::Array<unsigned char, Eigen::Dynamic, 1> M_j(n_samples);
			for (std::uint32_t ii = 0; ii < n_samples; ii++) {
				M_j[ii] = std::floor(std::min(dosage_j[ii], (scalarData) (L - 1e-6)) * invIntervalWidth);
			}
			// Decompress
			dosage_j = (M_j.cast<scalarData>() + 0.5) * intervalWidth;

			mu    = dosage_j.mean();
			sigma = (dosage_j - mu).square().sum();
			sigma = std::sqrt(sigma/(valid_count - 1.0));
		}

		if(params.normalise_genotypes) {
			dosage_j -= mu;
			dosage_j /= sigma;
		}
		G.col(jj) = dosage_j;

		jj++;
		valid_var_index_chunk.push_back(n_var_parsed - 1);
	}

	// need to resize G whilst retaining existing coefficients if while
	// loop exits early due to EOF.
	G.conservativeResize(n_samples, jj);
	assert( rsid_chunk.size() == jj );
	assert( chromosome_chunk.size() == jj );
	assert( position_chunk.size() == jj );
	assert( alleles_chunk.size() == jj );
	n_var_chunk = jj;

	chunk_missingness /= jj;
	if(chunk_missingness > 0.0) {
		std::cout << "Mean chunk missingness " << chunk_missingness << "(";
		std::cout << n_var_incomplete << "/" << n_var_chunk;
		std::cout << " variants contain >=1 imputed entry)" << std::endl;
	}

	maf.insert(maf.end(),               maf_chunk.begin(), maf_chunk.end());
	rsid.insert(rsid.end(),             rsid_chunk.begin(), rsid_chunk.end());
	chromosome.insert(chromosome.end(), chromosome_chunk.begin(), chromosome_chunk.end());
	SNPID.insert(SNPID.end(), SNPID_chunk.begin(), SNPID_chunk.end());
	position.insert(position.end(),     position_chunk.begin(), position_chunk.end());
	alleles.insert(alleles.end(),       alleles_chunk.begin(), alleles_chunk.end());
	valid_var_index.insert(valid_var_index.end(), valid_var_index_chunk.begin(), valid_var_index_chunk.end());

	if(jj == 0) {
		// Immediate EOF
		return false;
	} else {
		return true;
	}
}

template <typename EigenMat>
EigenMat reduce_mat_to_complete_cases( EigenMat& M,
                                       int n_cols,
                                       std::map<long, bool > incomplete_cases ) {
	// Remove rows contained in incomplete_cases
	long nn = M.rows();
	int n_incomplete;
	EigenMat M_tmp;

	// Create temporary matrix of complete cases
	n_incomplete = incomplete_cases.size();
	M_tmp.resize(nn - n_incomplete, n_cols);

	// Fill M_tmp with non-missing entries of M
	int ii_tmp = 0;
	for (int ii = 0; ii < nn; ii++) {
		if (incomplete_cases.count(ii) == 0) {
			for (int kk = 0; kk < n_cols; kk++) {
				M_tmp(ii_tmp, kk) = M(ii, kk);
			}
			ii_tmp++;
		}
	}

	return M_tmp;
}


inline Eigen::MatrixXd solve(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b) {
	Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
	if (fabs((double)((A * x - b).norm()/b.norm())) > 1e-8) {
		std::cout << "ERROR: could not solve covariate scatter matrix." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return x;
}

inline scalarData var(const EigenDataArrayX& vv){
	scalarData res = (vv - vv.mean()).square().sum();
	res /= ((scalarData) vv.rows() - 1.0);
	return res;
}

void sim_gaussian_noise(EigenRefDataArrayX vv, const double& sigma_sq, std::mt19937& generator){
	std::normal_distribution<scalarData> gaussian(0.0, std::sqrt(sigma_sq));
	long nrows = vv.rows();
	for (long ii = 0; ii < nrows; ii++) {
		vv(ii) = gaussian(generator);
	}
}

void Data::read_txt_file_w_context(const std::string &filename, const int &col_offset, EigenDataMatrix &M,
                                   std::vector<std::string> &M_snpids, std::vector<std::string> &col_names) {
	/*
	   Txt file where the first column is snp ids, then x-1 contextual,
	   then a matrix to be read into memory.

	   Reads file twice to ascertain number of lines.

	   col_offset - how many contextual columns to skip
	 */
	M_snpids.clear();

	// Reading from file
	boost_io::filtering_istream fg;
	std::string gz_str = ".gz";
	if (filename.find(gz_str) != std::string::npos) {
		fg.push(boost_io::gzip_decompressor());
	}
	fg.push(boost_io::file_source(filename));
	if (!fg) {
		std::cout << "ERROR: " << filename << " not opened." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// Read file twice to ascertain number of lines
	int n_lines = 0;
	std::string line;
	getline(fg, line);
	while (getline(fg, line)) {
		n_lines++;
	}
	fg.reset();
	if (filename.find(gz_str) != std::string::npos) {
		fg.push(boost_io::gzip_decompressor());
	}
	fg.push(boost_io::file_source(filename));

	// Reading column names
	if (!getline(fg, line)) {
		std::cout << "ERROR: " << filename << " contains zero lines" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::stringstream ss;
	std::string s;
	int n_cols = 0;
	col_names.clear();
	ss.clear();
	ss.str(line);
	while (ss >> s) {
		++n_cols;
		col_names.push_back(s);
	}
	assert(n_cols > col_offset);

	// Write remainder of file to Eigen matrix M
	M.resize(n_lines, n_cols - col_offset);
	int i = 0;
	double tmp_d;
	while (getline(fg, line)) {
		ss.clear();
		ss.str(line);
		for (int k = 0; k < n_cols; k++) {
			std::string sss;
			ss >> sss;
			if (k == 0) {
				M_snpids.push_back(sss);
			}
			if (k >= col_offset) {
				try{
					M(i, k-col_offset) = stod(sss);
				} catch (const std::invalid_argument &exc) {
					std::cout << "Found value " << sss << " on line " << i;
					std::cout << " of file " << filename << std::endl;
					throw std::runtime_error("Unexpected value");
				}
				if(!std::isfinite(M(i, k-col_offset))) {
					std::cout << "WARNING: " << std::endl;
					std::cout << "Found non-finite value " << sss << " on line " << i;
					std::cout << " of file " << filename << std::endl;
				}
			}
		}
		i++;
	}
	std::cout << "Reading matrix of size " << n_lines << " x " << n_cols - col_offset << " from " << filename << std::endl;
}

void Data::reduce_to_complete_cases() {
	// Remove any samples with incomplete covariates or phenotypes from
	// Y and W.
	// Note; other functions (eg. read_incl_sids) may add to incomplete_cases
	// Note; during unit testing sometimes only phenos or covars present.

	incomplete_cases.insert(missing_covars.begin(), missing_covars.end());
	incomplete_cases.insert(missing_envs.begin(), missing_envs.end());

	sample_is_invalid.clear();
	for (long ii = 0; ii < n_samples; ii++) {
		if (incomplete_cases.find(ii) == incomplete_cases.end()) {
			sample_is_invalid[ii] = false;
		} else {
			sample_is_invalid[ii] = true;
		}
	}

	if(n_covar > 0) {
		W = reduce_mat_to_complete_cases(W, n_covar, incomplete_cases);
	}
	if(n_env > 0) {
		E = reduce_mat_to_complete_cases(E, n_env, incomplete_cases);
	}
	n_samples -= incomplete_cases.size();
	missing_phenos.clear();
	missing_covars.clear();
	missing_envs.clear();

	if(!incomplete_cases.empty()) {
		std::cout << std::endl << "Reduced to " << n_samples << " samples with complete data." << std::endl;
	}
}

void Data::compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk) {
	assert(dXtEEX_chunk.rows() == n_var_chunk);
	assert(dXtEEX_chunk.cols() == n_env * n_env);
	EigenDataArrayX cl_j;
	for (int jj = 0; jj < n_var_chunk; jj++) {
		cl_j = G.col(jj);
		for (int ll = 0; ll < n_env; ll++) {
			for (int mm = 0; mm <= ll; mm++) {
				double x = (cl_j * E.array().col(ll) * E.array().col(mm) * cl_j).sum();
				dXtEEX_chunk(jj, ll*n_env + mm) = x;
				dXtEEX_chunk(jj, mm*n_env + ll) = x;
			}
		}
	}
}

void Data::regress_first_mat_from_second(const EigenDataMatrix &A, const std::string &Astring,
                                         const std::vector<std::string> &A_names, EigenDataMatrix &yy,
                                         const std::string &yy_string) {
	//
	std::cout << "Regressing " << Astring << " from " << yy_string << ":" << std::endl;
	unsigned long nnn = A_names.size();
	for(int cc = 0; cc < std::min(nnn, (unsigned long) 10); cc++) {
		std::cout << ( cc > 0 ? ", " : "" ) << A_names[cc];
	}
	if (nnn > 10) {
		std::cout << "... (" << nnn << " variables)";
	}
	std::cout << std::endl;

	Eigen::MatrixXd AtA = (A.transpose() * A).cast<double>();
	Eigen::MatrixXd Aty = (A.transpose() * yy).cast<double>();

	Eigen::MatrixXd bb = solve(AtA, Aty);
	yy -= A * bb.cast<scalarData>();
}

void Data::read_rsids(const std::string &filename, std::vector<std::string> &rsid_list) {
	rsid_list.clear();
	boost_io::filtering_istream fg;
	fg.push(boost_io::file_source(params.incl_rsids_file));
	if (!fg) {
		std::cout << "ERROR: " << params.incl_rsids_file << " not opened." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	std::stringstream ss;
	std::string line;
	while (getline(fg, line)) {
		ss.clear();
		ss.str(line);
		std::string s;
		while(ss >> s) {
			rsid_list.push_back(s);
		}
	}

	std::sort(rsid_list.begin(), rsid_list.end());
	rsid_list.erase(std::unique(rsid_list.begin(), rsid_list.end()), rsid_list.end());
}

void Data::output_correlations(EigenRefDataArrayXX vec) {
	outf << std::scientific << std::setprecision(8);
	for (long jj = 0; jj < n_var_chunk; jj++) {
		outf << SNPID_chunk[jj] << " " << chromosome_chunk[jj] << " " << rsid_chunk[jj] << " " << position_chunk[jj];
		outf << " " << alleles_chunk[jj][0] << " " << alleles_chunk[jj][1];
		for (int ll = 0; ll < n_env * n_env; ll++) {
			outf << " " << vec(jj, ll);
		}
		outf << std::endl;
	}
}

void Data::read_covar() {
	// Read covariates to Eigen matrix W
	if ( params.covar_file != "NULL" ) {
		fileUtils::read_matrix(params.covar_file, n_samples, W, covar_names, missing_covars);
		n_covar = covar_names.size();
	} else {
		throw std::logic_error( "Tried to read NULL covar file." );
	}
}

EigenDataMatrix Data::sort_coeffs(EigenRefDataMatrix orig,
                                  const std::unordered_map<std::string, long>& key_map,
                                  std::string key_type){
	assert(key_type == "SNPKEY" || key_type == "SNPID");

	// Create query
	genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(params.bgi_file);
	if (params.range) {
		if (params.chr.length() == 1) params.chr = "0" + params.chr;
		genfile::bgen::IndexQuery::GenomicRange rr1(params.chr, params.start, params.end);
		query->include_range( rr1 );
	}
	if(params.incl_snps) {
		read_incl_rsids();
		query->include_rsids(incl_rsid_list);
	}
	if(params.excl_snps) {
		read_excl_rsids();
		query->exclude_rsids(excl_rsid_list);
	}
	if(params.select_rsid) {
		std::sort(params.rsid.begin(), params.rsid.end());
		query->include_rsids(params.rsid);
	}
	query->initialise();

	genfile::bgen::View::UniquePtr view = genfile::bgen::View::create(params.bgen_file);
	view->set_query(query);

	if(key_type == "SNPKEY") {
		std::cout << "- Matching SNPKEYs given in --coeffs to bgen file..." << std::endl;
	} else {
		std::cout << "- Matching SNPIDs given in --coeffs to bgen file..." << std::endl;
	}
	EigenDataMatrix res = EigenDataMatrix::Zero(view->number_of_variants(), orig.cols());

	uint32_t pos_j;
	std::string rsid_j, chr_j, SNPID_j;
	std::vector< std::string > alleles_j;
	long n_matched = 0, jj = 0;
	while(view->read_variant(&SNPID_j, &rsid_j, &chr_j, &pos_j, &alleles_j)) {

		if(key_type == "SNPKEY") {
			std::string SNPKEY_j = gen_snpkey(chr_j, pos_j, alleles_j);
			auto it = key_map.find(SNPKEY_j);
			if (it != key_map.end()) {
				res.row(jj) = orig.row(it->second);
				n_matched++;
			}
		} else if(key_type == "SNPID") {
			auto it = key_map.find(SNPID_j);
			if (it != key_map.end()) {
				res.row(jj) = orig.row(it->second);
				n_matched++;
			}
		}

		view->ignore_genotype_data_block();
		jj++;
	}
	std::cout << "- " << n_matched << " matching SNPs found" << std::endl;
	return res;
}

void Data::read_coeffs_file(std::string filename,
                            EigenDataMatrix &coeffs_mat,
                            long& n_gxe_components) {
	// Read coefficients to eigen matrix B
	std::vector<std::string> coeff_names;
	read_file_header(filename, coeff_names);

	std::regex base("gamma.*");
	n_gxe_components = std::count_if(coeff_names.begin(), coeff_names.end(),
	                                 [&](const std::string& s){
		bool res = std::regex_match(s.begin(), s.end(), base);
		return res;
	});

	std::vector<std::string> case1a = {"beta"};
	std::vector<std::string> case2a = {"SNPKEY", "beta"};
	std::vector<std::string> case3a = {"SNPID", "beta"};

	if(coeff_names.size() == 1 + n_gxe_components &&
	   coeff_names[0] == case1a[0]) {
		fileUtils::read_matrix(filename, coeffs_mat, coeff_names);
	} else if (coeff_names.size() == 2 + n_gxe_components &&
	           coeff_names[0] == case2a[0] &&
	           coeff_names[1] == case2a[1]) {
		std::unordered_map<std::string, long> B_SNP_map;
		std::vector<std::string> B_SNPS;

		read_txt_file_w_context(filename, 1, coeffs_mat, B_SNPS, coeff_names);
		for (long jj = 0; jj < B_SNPS.size(); jj++) {
			B_SNP_map[B_SNPS[jj]] = jj;
		}

		coeffs_mat = Data::sort_coeffs(coeffs_mat, B_SNP_map, "SNPKEY");
	} else if (coeff_names.size() == 2 + n_gxe_components &&
	           coeff_names[0] == case3a[0] &&
	           coeff_names[1] == case3a[1]) {
		std::unordered_map<std::string, long> B_SNP_map;
		std::vector<std::string> B_SNPS;

		read_txt_file_w_context(filename, 1, coeffs_mat, B_SNPS, coeff_names);
		for (long jj = 0; jj < B_SNPS.size(); jj++) {
			B_SNP_map[B_SNPS[jj]] = jj;
		}

		coeffs_mat = Data::sort_coeffs(coeffs_mat, B_SNP_map, "SNPID");
	} else {
		throw std::logic_error("Unexpected header in " + filename);
	}

	// n_gxe_components = std::max(n_gxe_components, 1l);
}

void Data::read_environment() {
	fileUtils::read_matrix(params.env_file, n_samples, E, env_names, missing_envs);
	n_env = env_names.size();
}

void Data::read_env_profile() {
	std::vector<std::string> placeholder;
	fileUtils::read_matrix(params.env_profile_file, env_profile, placeholder);
	assert(env_profile.rows() == n_env || env_profile.cols() == n_env);
	if(env_profile.cols() == n_env && env_profile.rows() != n_env) {
		env_profile.transposeInPlace();
	}
}

void Data::read_env_profile2() {
	std::vector<std::string> placeholder;
	fileUtils::read_matrix(params.env_profile_file2, env_profile2, placeholder);
	assert(env_profile2.rows() == n_env || env_profile2.cols() == n_env);
	if(env_profile2.cols() == n_env && env_profile2.rows() != n_env) {
		env_profile2.transposeInPlace();
	}
}

void Data::read_incl_sids() {
	boost_io::filtering_istream fg;
	std::string gz_str = ".gz";
	if (params.incl_sids_file.find(gz_str) != std::string::npos) {
		fg.push(boost_io::gzip_decompressor());
	}
	fg.push(boost_io::file_source(params.incl_sids_file));
	if (!fg) {
		std::cout << "ERROR: " << params.incl_sids_file << " not opened." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	std::vector<std::string> bgen_ids;
	bgenView->get_sample_ids(
		[&]( std::string const& id ) {
		bgen_ids.push_back(id);
	}
		);

	std::stringstream ss;
	std::string line;
	std::set<std::string> user_sample_ids;
	while (getline(fg, line)) {
		ss.clear();
		ss.str(line);
		std::string s;
		ss >> s;
		user_sample_ids.insert(s);
	}

	std::set<std::string>::iterator it;
	for (long ii = 0; ii < n_samples; ii++) {
		it = user_sample_ids.find(bgen_ids[ii]);
		if (it == user_sample_ids.end()) {
			incomplete_cases[ii] = true;
		}
	}

	std::cout << "Subsetted down to " << user_sample_ids.size() << " ids from --incl_sample_ids";
	std::cout << std::endl;
}

void Data::dump_coeffs(const std::string& filename,
                       const std::vector<std::string>& header,
                       Eigen::Ref<Eigen::MatrixXd> coeffs,
                       Eigen::Ref<Eigen::MatrixXd> coeffs2){
	long n_cols = coeffs.cols();
	long n_rows = coeffs.rows();
	assert(header.size() == 7 + coeffs.cols() + coeffs2.cols());

	boost_io::filtering_ostream outf;
	std::string gz_str = ".gz";
	if (filename.find(gz_str) != std::string::npos) {
		outf.push(boost_io::gzip_compressor());
	}
	outf.push(boost_io::file_sink(filename));

	for (long cc = 0; cc < header.size(); cc++) {
		outf << header[cc];
		outf << (cc != header.size() - 1 ? " " : "");
	}
	outf << std::endl;
	for (long ii = 0; ii < n_rows; ii++) {
		outf << SNPID[ii] << chromosome[ii] << " " << rsid[ii] << " " << position[ii];
		outf << " " << alleles[ii][0] << " " << alleles[ii][1];
		outf << " " << maf[ii];

		for (long cc = 0; cc < n_cols; cc++) {
			outf << " " << coeffs(ii, cc);
		}
		if(coeffs2.rows() > 0) {
			for (long cc = 0; cc < n_cols; cc++) {
				outf << " " << coeffs2(ii, cc);
			}
		}
		outf << std::endl;
	}

	outf.pop();
	boost_io::close(outf);
}
