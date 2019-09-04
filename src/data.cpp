//
// Created by kerin on 2019-02-26.
//

#include "parameters.hpp"
#include "utils.hpp"
#include "data.hpp"

#include "tools/eigen3.3/Dense"
#include "bgen_parser.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstddef>     // for ptrdiff_t class
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <set>
#include <stdexcept>


Data::Data(const parameters &my_params) : params(my_params){
	bgenView = genfile::bgen::View::create(my_params.bgen_file);
	bgen_pass = true;
	n_samples = bgenView->number_of_samples();
	n_var_parsed = 0;
	n_total_var = 0;
	n_constant_variance = 0;
	n_covar = 0;
	n_env = 0;

	match_snpkeys = false;
	match_snpids = false;

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
	std::cout << "Analysis finished at " << std::ctime(&end_time);
	std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

	boost_io::close(outf);
	boost_io::close(outf_pred);
	boost_io::close(outf_coeffs);
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
		center_matrix( W, n_covar );
		scale_matrix( W, n_covar, covar_names );
	}

	int ch = 0;
	double tmp_x, tmp_z, mean_z;
	s_x = 0.0, s_z = 0.0;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;
		assert(n_var == G.cols());
		n_total_var += G.cols();

		// Compute sum of column variances
		for (int kk = 0; kk < n_var; kk++) {
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
	if(n_constant_variance > 0) {
		std::cout << " Removed " << n_constant_variance  << " column(s) variance < 1e-12 or mac < 5:" << std::endl;
	}
}

void Data::pred_pheno() {
	// Gen vector of fitted effects Y_hat = age + X beta + Z gamma
	// Gaussian noise added when writing to file

	// Step 1; Read in raw covariates
	// - also makes a record of missing values
	if(params.covar_file != "NULL") {
		read_covar();
	}
	if(params.env_file != "NULL") {
		read_environment();
	}
	read_coeffs();
	if(params.coeffs2_file != "NULL") {
		read_coeffs2();
	}

	if(n_env > 1) {
		std::cout << "WARNING: " << n_env << " cols detected in " << params.env_file << ". Just using first" << std::endl;
		E = E.col(0);
		n_env = 1;
		env_names.resize(1);
	}

	// Step 2; Reduce raw covariates and phenotypes to complete cases
	// - may change value of n_samples
	// - will also skip these cases when reading bgen later
	reduce_to_complete_cases();

	// Step 3; Normalise covars
	if(n_covar > 0 && !params.use_raw_covars) {
		center_matrix( W, n_covar );
		scale_matrix( W, n_covar, covar_names );
	}
	if(n_env > 0 && !params.use_raw_env) {
		center_matrix( E, n_env );
		scale_matrix( E, n_env, env_names );
	}

	// S2; genetic or interaction effects
	Xb = EigenDataArrayX::Zero(n_samples);
	Xg = EigenDataArrayX::Zero(n_samples);
	Zg = EigenDataArrayX::Zero(n_samples);
	Xb2 = EigenDataArrayX::Zero(n_samples);
	Xg2 = EigenDataArrayX::Zero(n_samples);
	Zg2 = EigenDataArrayX::Zero(n_samples);
	Wtau = EigenDataVector::Zero(n_samples);
	Ealpha = EigenDataVector::Zero(n_samples);
	int ch = 0;
	long n_matched = 0;
	beta_vec1.reserve(bgenView->number_of_variants());
	beta_vec2.reserve(bgenView->number_of_variants());
	gamma_vec1.reserve(bgenView->number_of_variants());
	gamma_vec2.reserve(bgenView->number_of_variants());

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
			std::cout << "- env-vector unprocessed" << std::endl;
		} else {
			std::cout << "- env-vector normalised to mean zero, variance one" << std::endl;
		}
	}
	std::cout << std::endl;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		if(ch % 10 == 0) {
			std::cout << "Chunk " << ch+1 << " read (size " << n_var;
			std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
			std::cout << " variants parsed)" << std::endl;
		}

		assert(G.cols() == n_var);
		assert(G.rows() == n_samples);
		if(!(match_snpkeys || match_snpids)) {
			assert(B.rows() >= n_total_var + n_var);
		}

		// Add effects to Y
		long coeff_index = 0;
		for (int kk = 0; kk < n_var; kk++) {
			if(match_snpkeys) {
				auto it = B_SNPKEYS_map.find(SNPKEYS[kk]);
				if (it != B_SNPKEYS_map.end()) {
					coeff_index = it->second;
					n_matched++;
				} else {
					beta_vec1.push_back(0);
					beta_vec2.push_back(0);
					gamma_vec1.push_back(0);
					gamma_vec2.push_back(0);
					continue;
				}
			} else if(match_snpids) {
				auto it = B_SNPIDS_map.find(SNPID[kk]);
				if (it != B_SNPIDS_map.end()) {
					coeff_index = it->second;
					n_matched++;
				} else {
					beta_vec1.push_back(0);
					beta_vec2.push_back(0);
					gamma_vec1.push_back(0);
					gamma_vec2.push_back(0);
					continue;
				}
			} else {
				coeff_index = kk + n_total_var;
			}

			Xb += G.col(kk).array() * B(coeff_index, 0);
			beta_vec1.push_back(B(coeff_index, 0));
			if(n_env > 0) {
				Xg += G.col(kk).array() * B(coeff_index, 1);
				Zg += G.col(kk).array() * E.array() * B(coeff_index, 1);
				gamma_vec1.push_back(B(coeff_index, 1));
			} else {
				gamma_vec1.push_back(0);
			}

			if(params.coeffs2_file != "NULL") {
				// WARNING: ASSUMED THAT SNPKEYS IN COEFFS and COEFFS2 ARE IDENTICAL
				Xb2 += G.col(kk).array() * B2(coeff_index, 0);
				beta_vec2.push_back(B2(coeff_index, 0));
				if(n_env > 0) {
					Xg2 += G.col(kk).array() * B2(coeff_index, 1);
					Zg2 += G.col(kk).array() * E.array() * B2(coeff_index, 1);
					gamma_vec2.push_back(B2(coeff_index, 1));
				} else {
					gamma_vec2.push_back(0);
				}
			}
		}

		n_total_var += n_var;
		ch++;
	}
	if(n_total_var != B.rows() & !(match_snpkeys || match_snpids)) {
		std::cout << "ERROR: n var read in = " << n_total_var << std::endl;
		std::cout << "ERROR: n coeffs read in = " << B.rows() << std::endl;
		assert(n_total_var == B.rows());
	}
	if(n_constant_variance > 0) {
		std::cout << "- Removed " << n_constant_variance  << " SNP(s) with zero variance" << std::endl;
	}
	if(match_snpkeys || match_snpids) {
		std::cout << "- " << n_matched << " SNPs found matching those given in --coeffs" << std::endl;
	} else {
		std::cout << "- " << n_total_var << " SNPs used to construct genetic effects" << std::endl;
	}
	if(n_env > 0) {
		Zg = Xg * E.col(0).array();
		Zg2 = Xg2 * E.col(0).array();
	}
	Y = Xb + Zg + Xb2 + Zg2;
}

void Data::sim_pheno() {

	// Get predicted effects
	pred_pheno();

	// Generate random effects
	// http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
	// Random seed
	if(params.random_seed == -1) {
		std::random_device rd;
		params.random_seed = rd();
	}
	std::cout << std::endl << "Initialising random sample generator with seed " << params.random_seed << std::endl;
	std::mt19937 generator{params.random_seed};

	// additive noise
	noise.resize(n_samples);
	sim_gaussian_noise(noise, 1, generator);
	noise    *= std::sqrt(params.sigma / var(noise));

	// target variances
	if(params.rescale_coeffs) {
		std::cout << "Rescaling components to ensure that sample ";
		std::cout << "heritability matches expected heritability" << std::endl;
		double sigma_b, sigma_g = 0, sigma_b2, sigma_g2 = 0, sigma_c = 0, sigma_e = 0;
		sigma_b = params.hb / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
		sigma_b2 = params.hb2 / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
		if(n_env > 0) {
			sigma_g = params.hg / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
			sigma_g2 = params.hg2 / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
			sigma_e = params.he / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
		}
		if(n_covar > 0) {
			sigma_c = params.hc / (1.0 - params.hc - params.he - params.hb - params.hg - params.hb2 - params.hg2);
		}

		// additive covar effects
		if(params.hc > 0) {
			EigenDataVector tau(n_covar);
			sim_gaussian_noise(tau, 1, generator);
			Wtau = W * tau;
			double sf = params.sigma * sigma_c / var(Wtau);
			Wtau     *= std::sqrt(sf);
		}

		// additive env effects
		if(params.he > 0) {
			EigenDataVector alpha(n_env);
			sim_gaussian_noise(alpha, 1, generator);
			Ealpha = E * alpha;
			double sf = params.sigma * sigma_e / var(Ealpha);
			Ealpha    *= std::sqrt(sf);
		}

		// rescaling for correct heritability
		scalarData var_xb = var(Xb);
		Xb       *= std::sqrt(params.sigma * sigma_b / var_xb);
		B.col(0) *= std::sqrt(params.sigma * sigma_b / var_xb);

		if(n_env > 0) {
			scalarData var_zg = var(Zg);
			Zg       *= std::sqrt(params.sigma * sigma_g / var_zg);
			Xg       *= std::sqrt(params.sigma * sigma_g / var_zg);
			B.col(1) *= std::sqrt(params.sigma * sigma_g / var_zg);
		}

		if(params.coeffs2_file != "NULL") {
			scalarData var_xb = var(Xb2);
			Xb2       *= std::sqrt(params.sigma * sigma_b2 / var_xb);
			B2.col(0) *= std::sqrt(params.sigma * sigma_b2 / var_xb);

			if(n_env > 0) {
				scalarData var_zg = var(Zg2);
				Zg2       *= std::sqrt(params.sigma * sigma_g2 / var_zg);
				Xg2       *= std::sqrt(params.sigma * sigma_g2 / var_zg);
				B2.col(1) *= std::sqrt(params.sigma * sigma_g2 / var_zg);
			}
		}
	}

	Y = Wtau + Ealpha + Xb + Zg + Xb2 + Zg2 + noise;
	std::cout << "Heritability partition (ignoring covariance):" << std::endl;
	std::cout << "h2-G = " << var(Xb) / var(Y) << std::endl;
	if(n_env > 0) std::cout << "h2-GxE = " << var(Zg) / var(Y) << std::endl;
	if(params.coeffs2_file != "NULL") std::cout << "hb2 = " << var(Xb2) / var(Y) << std::endl;
	if(params.coeffs2_file != "NULL") std::cout << "hg2 = " << var(Zg2) / var(Y) << std::endl;
	if(n_env > 0) std::cout << "h2-Env = " << var(Ealpha) / var(Y) << std::endl;
	if(n_covar > 0) std::cout << "h2-Covars = " << var(Wtau) / var(Y) << std::endl;
	std::cout << "Noise has variance: " << params.sigma << std::endl;
}

void Data::output_init() {
	std::string ofile      = fstream_init(outf, params.out_file, "");


	// Output header for ssv file
	if(params.mode_ssv) {
		std::cout << "Writing ssv stats to " << ofile << std::endl;
		outf << "s_x\ts_z\tn\tp" << std::endl;
	} else if (params.mode_pred_pheno) {
		std::cout << "Writing predicted phenotype to " << ofile << std::endl;
		outf << "y" << std::endl;

		std::string ofile_pred   = fstream_init(outf_pred, params.out_file, "_predicted_effects");
		std::cout << "Writing predicted effects to " << ofile_pred << std::endl;
		outf_pred << "Wtau\tEalpha\tXbeta\teta\tXgamma\tZgamma" << std::endl;

		std::string ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_true_rescaled_coeffs");
		std::cout << "Writing rescaled coeffs to " << ofile_coeffs << std::endl;
		outf_coeffs << "chr rsid pos a0 a1 af beta gamma";
		if(params.coeffs2_file != "NULL"){
			outf_coeffs << " beta2 gamma2";
		}
		outf_coeffs << std::endl;

	} else if(params.mode_gen_pheno || params.mode_gen2_pheno) {
		std::string ofile_pred   = fstream_init(outf_pred, params.out_file, "_predicted_effects");

		std::cout << "Writing predicted effects to " << ofile_pred << std::endl;
		outf_pred << "Wtau\tEalpha\tXbeta\teta\tXgamma\tZgamma\tnoise" << std::endl;

		if(params.sim_w_noise) {
			std::cout << "Writing simulated phenotype to " << ofile << std::endl;
			outf << "y" << std::endl;
		}

		std::string ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_true_rescaled_coeffs");
		std::cout << "Writing rescaled coeffs to " << ofile_coeffs << std::endl;
		outf_coeffs << "chr rsid pos a0 a1 af beta gamma";
		if(params.coeffs2_file != "NULL"){
			outf_coeffs << " beta2 gamma2";
		}
		outf_coeffs << std::endl;

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
		outf << n_total_var << std::endl;
	} else if (params.mode_pred_pheno) {
		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf_pred << Wtau(ii) <<"\t" << Ealpha(ii);
			outf_pred << "\t" << Xb(ii) << "\t" << E(ii, 0) << "\t" << Xg(ii);
			outf_pred << "\t" << Xg(ii) * E(ii, 0) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf << Y(ii) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_total_var; ii++) {
			outf_coeffs << chromosome_cum[ii] << " " << rsid_cum[ii] << " " << position_cum[ii];
			outf_coeffs << " " << alleles_cum[ii][0] << " " << alleles_cum[ii][1];
			outf_coeffs << " " << maf_cum[ii] << " " << beta_vec1[ii] << " " << gamma_vec1[ii];
			if(params.coeffs2_file != "NULL"){
				outf_coeffs << " " << beta_vec2[ii] << " " << gamma_vec2[ii];
			}
			outf_coeffs << std::endl;
		}
	} else if(params.mode_gen_pheno || params.mode_gen2_pheno) {
		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf_pred << Wtau(ii) <<"\t" << Ealpha(ii) <<"\t" << Xb(ii) << "\t" << E(ii, 0);
			outf_pred << "\t" << Xg(ii) << "\t" << Xg(ii) * E(ii, 0) << "\t" << noise(ii) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf << Y(ii) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_total_var; ii++) {
			outf_coeffs << chromosome_cum[ii] << " " << rsid_cum[ii] << " " << position_cum[ii];
			outf_coeffs << " " << alleles_cum[ii][0] << " " << alleles_cum[ii][1];
			outf_coeffs << " " << maf_cum[ii] << " " << beta_vec1[ii] << " " << gamma_vec1[ii];
			if(params.coeffs2_file != "NULL"){
				outf_coeffs << " " << beta_vec2[ii] << " " << gamma_vec2[ii];
			}
			outf_coeffs << std::endl;
		}
	}

	if(params.print_causal_rsids) {
		std::string ofile_coeffs;
		ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_nonzero_beta_rsids");
		std::cout << "Writing rsids of nonzero betas to " << ofile_coeffs << std::endl;
		for (std::size_t ii = 0; ii < n_total_var; ii++) {
			if(std::abs(B(ii, 0)) > 1e-6) {
				outf_coeffs << rsid_cum[ii] << std::endl;
			}
		}

		ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_nonzero_gamma_rsids");
		std::cout << "Writing rsids of nonzero gammas to " << ofile_coeffs << std::endl;
		for (std::size_t ii = 0; ii < n_total_var; ii++) {
			if(std::abs(B(ii, 1)) > 1e-6) {
				outf_coeffs << rsid_cum[ii] << std::endl;
			}
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
	std::string chr_j;
	uint32_t pos_j;
	std::string rsid_j;
	std::vector< std::string > alleles_j;
	std::string SNPID_j;

	long nInvalid = sample_is_invalid.size() - n_samples;
	DosageSetter setter_v2(sample_is_invalid, nInvalid);

	double chunk_missingness = 0;
	long n_var_incomplete = 0;

	// Wipe variant context from last chunk
	maf.clear();
	info.clear();
	rsid.clear();
	chromosome.clear();
	position.clear();
	alleles.clear();
	SNPKEYS.clear();
	SNPID.clear();

	// Resize genotype matrix
	G.resize(n_samples, params.chunk_size);

	long jj = 0;
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

		maf.push_back(maf_j);
		info.push_back(info_j);
		rsid.push_back(rsid_j);
		chromosome.push_back(chr_j);
		position.push_back(pos_j);
		alleles.push_back(alleles_j);
		SNPID.push_back(SNPID_j);

		std::string key_j = gen_snpkey(chr_j, pos_j, alleles_j);
		SNPKEYS.push_back(key_j);
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
	}

	// need to resize G whilst retaining existing coefficients if while
	// loop exits early due to EOF.
	G.conservativeResize(n_samples, jj);
	assert( rsid.size() == jj );
	assert( chromosome.size() == jj );
	assert( position.size() == jj );
	assert( alleles.size() == jj );
	n_var = jj;

	chunk_missingness /= jj;
	if(chunk_missingness > 0.0) {
		std::cout << "Mean chunk missingness " << chunk_missingness << "(";
		std::cout << n_var_incomplete << "/" << n_var;
		std::cout << " variants contain >=1 imputed entry)" << std::endl;
	}

	maf_cum.insert(maf_cum.end(),               maf.begin(), maf.end());
	rsid_cum.insert(rsid_cum.end(),             rsid.begin(), rsid.end());
	chromosome_cum.insert(chromosome_cum.end(), chromosome.begin(), chromosome.end());
	position_cum.insert(position_cum.end(),     position.begin(), position.end());
	alleles_cum.insert(alleles_cum.end(),       alleles.begin(), alleles.end());

	if(jj == 0) {
		// Immediate EOF
		return false;
	} else {
		return true;
	}
}

template <typename EigenMat>
EigenMat reduce_mat_to_complete_cases( EigenMat& M,
                                       bool& matrix_reduced,
                                       int n_cols,
                                       std::map<long, bool > incomplete_cases ) {
	// Remove rows contained in incomplete_cases
	long nn = M.rows();
	int n_incomplete;
	EigenMat M_tmp;
	if (matrix_reduced) {
		throw std::runtime_error("ERROR: Trying to remove incomplete cases twice...");
	}

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
		} else {
			// std::cout << "Deleting row " << ii << std::endl;
		}
	}

	// Assign new values to reference variables
	matrix_reduced = true;
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

void Data::read_grid_file(const std::string &filename, EigenDataMatrix &M, std::vector<std::string> &col_names) {
	// Used in mode_vb only.
	// Slightly different from read_txt_file in that I don't know
	// how many rows there will be and we can assume no missing values.

	boost_io::filtering_istream fg;
	fg.push(boost_io::file_source(filename));
	if (!fg) {
		std::cout << "ERROR: " << filename << " not opened." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// Read file twice to acertain number of lines
	std::string line;
	int n_grid = 0;
	getline(fg, line);
	while (getline(fg, line)) {
		n_grid++;
	}
	fg.reset();
	fg.push(boost_io::file_source(filename));

	// Reading column names
	if (!getline(fg, line)) {
		std::cout << "ERROR: " << filename << " not read." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::stringstream ss;
	std::string s;
	int n_cols = 0;
	ss.clear();
	ss.str(line);
	while (ss >> s) {
		++n_cols;
		col_names.push_back(s);
	}
	std::cout << " Reading matrix of size " << n_grid << " x " << n_cols << " from " << filename << std::endl;

	// Write remainder of file to Eigen matrix M
	M.resize(n_grid, n_cols);
	int i = 0;
	double tmp_d;
	try {
		while (getline(fg, line)) {
			if (i >= n_grid) {
				throw std::runtime_error("ERROR: could not convert txt file (too many lines).");
			}
			ss.clear();
			ss.str(line);
			for (int k = 0; k < n_cols; k++) {
				std::string sss;
				ss >> sss;
				try{
					tmp_d = stod(sss);
				} catch (const std::invalid_argument &exc) {
					std::cout << sss << " on line " << i << std::endl;
					throw;
				}
				if(!std::isfinite(tmp_d)) {
					std::cout << "WARNING: " << std::endl;
					std::cout << "Found non-finite value " << sss << " on line " << i;
					std::cout << " of file " << filename << std::endl;
				}

				M(i, k) = tmp_d;
			}
			i++;
		}
		if (i < n_grid) {
			throw std::runtime_error("ERROR: could not convert txt file (too few lines).");
		}
	} catch (const std::exception &exc) {
		throw;
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
	std::cout << n_lines << " rows found in " << filename << std::endl;
}

void Data::center_matrix(EigenDataMatrix &M, int &n_cols) {
	// Center eigen matrix passed by reference.
	// Only call on matrixes which have been reduced to complete cases,
	// as no check for incomplete rows.

	std::vector<size_t> keep;
	for (int k = 0; k < n_cols; k++) {
		double mu = 0.0;
		double count = 0;
		for (int i = 0; i < n_samples; i++) {
			mu += M(i, k);
			count += 1;
		}

		mu = mu / count;
		for (int i = 0; i < n_samples; i++) {
			M(i, k) -= mu;
		}
		// std::cout << "Mean centered matrix:" << std::endl << M << std::endl;
	}
}

void Data::scale_matrix(EigenDataMatrix &M, int &n_cols, std::vector<std::string> &col_names) {
	// Scale eigen matrix passed by reference.
	// Removes columns with zero variance + updates col_names.
	// Only call on matrixes which have been reduced to complete cases,
	// as no check for incomplete rows.

	std::vector<size_t> keep;
	std::vector<std::string> keep_names;
	std::vector<std::string> reject_names;
	for (int k = 0; k < n_cols; k++) {
		double sigma = 0.0;
		double count = 0;
		for (int i = 0; i < n_samples; i++) {
			double val = M(i, k);
			sigma += val * val;
			count += 1;
		}

		sigma = sqrt(sigma/(count - 1));
		if (sigma > 1e-12) {
			for (int i = 0; i < n_samples; i++) {
				M(i, k) /= sigma;
			}
			keep.push_back(k);
			keep_names.push_back(col_names[k]);
		} else {
			reject_names.push_back(col_names[k]);
		}
	}

	if (keep.size() != n_cols) {
		n_constant_variance += (n_cols - keep.size());
		// std::cout << " Removing " << (n_cols - keep.size())  << " column(s) with zero variance:" << std::endl;
		// for(int kk = 0; kk < (n_cols - keep.size()); kk++){
		//  std::cout << reject_names[kk] << std::endl;
		// }
//			M = getCols(M, keep);
		for (std::size_t i = 0; i < keep.size(); i++) {
			M.col(i) = M.col(keep[i]);
		}
		M.conservativeResize(M.rows(), keep.size());

		n_cols = keep.size();
		col_names = keep_names;
	}

	if (n_cols == 0) {
		std::cout << "WARNING: No columns left with nonzero variance after scale_matrix()" << std::endl;
	}
}

void Data::scale_matrix_conserved(EigenDataMatrix &M, int &n_cols) {
	// Scale eigen matrix passed by reference.
	// Does not remove columns with zero variance.
	// Only call on matrixes which have been reduced to complete cases,
	// as no check for incomplete rows.

	for (int k = 0; k < n_cols; k++) {
		double sigma = 0.0;
		double count = 0;
		for (int i = 0; i < n_samples; i++) {
			double val = M(i, k);
			sigma += val * val;
			count += 1;
		}

		sigma = sqrt(sigma/(count - 1));
		if (sigma > 1e-12) {
			for (int i = 0; i < n_samples; i++) {
				M(i, k) /= sigma;
			}
		}
	}
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
		W = reduce_mat_to_complete_cases( W, W_reduced, n_covar, incomplete_cases );
	}
	if(n_env > 0) {
		E = reduce_mat_to_complete_cases( E, E_reduced, n_env, incomplete_cases );
	}
	n_samples -= incomplete_cases.size();
	missing_phenos.clear();
	missing_covars.clear();
	missing_envs.clear();

	if(!incomplete_cases.empty()) {
		std::cout << "Reduced to " << n_samples << " samples with complete data." << std::endl;
	}
}

void Data::compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk) {
	assert(dXtEEX_chunk.rows() == n_var);
	assert(dXtEEX_chunk.cols() == n_env * n_env);
	EigenDataArrayX cl_j;
	for (int jj = 0; jj < n_var; jj++) {
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

void Data::print_keys() {
	// Step 4; compute correlations
	int ch = 0;
	reduce_to_complete_cases();
	// std::cout << "num samples: " << n_samples << std::endl;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;

		for (std::size_t ii = 0; ii < n_var; ii++) {
			outf << SNPID[ii] << " " << chromosome[ii] << " " << rsid[ii];
			outf << " " << position[ii];
			outf << " " << alleles[ii][0] << " " << alleles[ii][1];
			outf << " " << maf[ii] << " " << info[ii] << std::endl;
		}

		n_total_var += n_var;
		ch++;
	}

	if(n_constant_variance > 0) {
		std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
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
	center_matrix( E, n_env );
	scale_matrix( E, n_env, env_names );

	if(n_covar > 0) {
		center_matrix(W, n_covar);
		scale_matrix(W, n_covar, covar_names);
	}
	if(n_covar > 0) {
		regress_first_mat_from_second(W, "covars", covar_names, E, "env");
	}

	// Step 4; compute correlations
	int ch = 0;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;

		assert(n_var == G.cols());

		EigenDataArrayXX dXtEEX_chunk(n_var, n_env * n_env);
		compute_correlations_chunk(dXtEEX_chunk);
		output_correlations(dXtEEX_chunk);

		n_total_var += n_var;
		ch++;
	}

	if(n_constant_variance > 0) {
		std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
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

std::string Data::fstream_init(boost_io::filtering_ostream &my_outf, const std::string &orig_file_path,
                               const std::string &file_suffix) {

	std::string filepath   = orig_file_path;
	std::string dir        = filepath.substr(0, filepath.rfind('/')+1);
	std::string stem_w_dir = filepath.substr(0, filepath.find('.'));
	std::string stem       = stem_w_dir.substr(stem_w_dir.rfind('/')+1, stem_w_dir.size());
	std::string ext        = filepath.substr(filepath.find('.'), filepath.size());

	std::string ofile      = dir + stem + file_suffix + ext;

	my_outf.reset();
	std::string gz_str = ".gz";
	if (orig_file_path.find(gz_str) != std::string::npos) {
		my_outf.push(boost_io::gzip_compressor());
	}
	my_outf.push(boost_io::file_sink(ofile));
	return ofile;
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
