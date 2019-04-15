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
#include <chrono>
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


template <typename EigenMat>
EigenMat reduce_mat_to_complete_cases( EigenMat& M,
                                       bool& matrix_reduced,
                                       int n_cols,
                                       std::map< int, bool > incomplete_cases ) {
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



/* Inline functions */

inline size_t numRows(const Eigen::MatrixXd &A) {
	return A.rows();
}

inline size_t numCols(const Eigen::MatrixXd &A) {
	return A.cols();
}

inline void setCol(Eigen::MatrixXd &A, const Eigen::VectorXd &v, size_t col) {
	assert(numRows(v) == numRows(A));
	A.col(col) = v;
}

inline Eigen::VectorXd getCol(const Eigen::MatrixXd &A, size_t col) {
	return A.col(col);
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
				std::string s;
				ss >> s;
				try{
					tmp_d = stod(s);
				} catch (const std::invalid_argument &exc) {
					std::cout << s << " on line " << i << std::endl;
					throw;
				}

				M(i, k) = tmp_d;
			}
			i++;                                                                                     // loop should end at i == n_grid
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
	getline(fg, line);                             // skip header
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
			}
		}
		i++;                                                         // loop should end at i == n_samples
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
	incomplete_cases.insert(missing_phenos.begin(), missing_phenos.end());
	incomplete_cases.insert(missing_envs.begin(), missing_envs.end());
	if(params.pheno_file != "NULL") {
		Y = reduce_mat_to_complete_cases( Y, Y_reduced, n_pheno, incomplete_cases );
	}
	if(params.covar_file != "NULL" || n_covar > 0) {
		W = reduce_mat_to_complete_cases( W, W_reduced, n_covar, incomplete_cases );
	}
	if(params.env_file != "NULL" || n_env > 0) {
		E = reduce_mat_to_complete_cases( E, E_reduced, n_env, incomplete_cases );
	}
	n_samples -= incomplete_cases.size();
	missing_phenos.clear();
	missing_covars.clear();
	missing_envs.clear();

	std::cout << "Reduced to " << n_samples << " samples with complete data";
	std::cout << " across covariates and phenotype." << std::endl;
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
	Zg = EigenDataArrayX::Zero(n_samples);
	Xb2 = EigenDataArrayX::Zero(n_samples);
	Zg2 = EigenDataArrayX::Zero(n_samples);
	int ch = 0;
	long n_matched = 0;
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
					continue;
				}
			} else if(match_snpids) {
				auto it = B_SNPIDS_map.find(SNPID[kk]);
				if (it != B_SNPIDS_map.end()) {
					coeff_index = it->second;
					n_matched++;
				} else {
					continue;
				}
			} else {
				coeff_index = kk + n_total_var;
			}

			Xb += G.col(kk).array() * B(coeff_index, 0);
			if(n_env > 0) {
				Zg += G.col(kk).array() * E.array() * B(coeff_index, 1);
			}

			if(params.coeffs2_file != "NULL") {
				assert(!(match_snpkeys || match_snpids));
				// Not yet implemented
				Xb2 += G.col(kk).array() * B2(coeff_index, 0);
				if(n_env > 0) {
					Zg2 += G.col(kk).array() * E.array() * B2(coeff_index, 1);
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
		std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
	}
	if(match_snpkeys || match_snpids) {
		std::cout << n_matched << " SNPs found matching those given in --coeffs" << std::endl;
	}
	Y = Xb + Zg + Xb2 + Zg2;
}

void Data::gen_pheno() {

	// Get predicted effects
	pred_pheno();

	// Generate random effects
	// http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
	// Random seed
	if(params.random_seed == -1) {
		std::random_device rd;
		params.random_seed = rd();
	}
	std::cout << "Initialising random sample generator with seed " << params.random_seed << std::endl;
	std::mt19937 generator{params.random_seed};

	// target variances
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

	// additive noise
	noise.resize(n_samples);
	sim_gaussian_noise(noise, 1, generator);
	noise    *= std::sqrt(params.sigma / var(noise));

	// additive covar effects
	if(params.hc > 0) {
		EigenDataVector tau(n_covar);
		sim_gaussian_noise(tau, 1, generator);
		Wtau = W * tau;
		double sf = params.sigma * sigma_c / var(Wtau);
		Wtau     *= std::sqrt(sf);
	} else {
		Wtau = EigenDataVector::Zero(n_samples);
	}

	// additive env effects
	if(params.he > 0) {
		EigenDataVector alpha(n_env);
		sim_gaussian_noise(alpha, 1, generator);
		Ealpha = E * alpha;
		double sf = params.sigma * sigma_e / var(Ealpha);
		Ealpha    *= std::sqrt(sf);
	} else {
		Ealpha = EigenDataVector::Zero(n_samples);
	}

	// rescaling for correct heritability
	scalarData var_xb = var(Xb);
	Xb       *= std::sqrt(params.sigma * sigma_b / var_xb);
	B.col(0) *= std::sqrt(params.sigma * sigma_b / var_xb);

	if(n_env > 0) {
		scalarData var_zg = var(Zg);
		Zg       *= std::sqrt(params.sigma * sigma_g / var_zg);
		B.col(1) *= std::sqrt(params.sigma * sigma_g / var_zg);
	}

	if(params.coeffs2_file != "NULL") {
		scalarData var_xb = var(Xb2);
		Xb2       *= std::sqrt(params.sigma * sigma_b2 / var_xb);
		B2.col(0) *= std::sqrt(params.sigma * sigma_b2 / var_xb);

		if(n_env > 0) {
			scalarData var_zg = var(Zg2);
			Zg2       *= std::sqrt(params.sigma * sigma_g2 / var_zg);
			B2.col(1) *= std::sqrt(params.sigma * sigma_g2 / var_zg);
		}
	}

	Y = Wtau + Ealpha + Xb + Zg + Xb2 + Zg2 + noise;
	std::cout << "Heritability partition (ignoring covariance):" << std::endl;
	std::cout << "hb = " << var(Xb) / var(Y) << std::endl;
	std::cout << "hg = " << var(Zg) / var(Y) << std::endl;
	std::cout << "hb2 = " << var(Xb2) / var(Y) << std::endl;
	std::cout << "hg2 = " << var(Zg2) / var(Y) << std::endl;
	std::cout << "he = " << var(Ealpha) / var(Y) << std::endl;
	std::cout << "hc = " << var(Wtau) / var(Y) << std::endl;
	std::cout << "Noise has variance: " << params.sigma << std::endl;
}

void Data::gen2_pheno() {
	// Controlling heritability of generated phenotype with variance of
	// noise subsequently added, a la BSLMM.

	if(params.covar_file != "NULL") {
		read_covar();
	}
	if(params.env_file != "NULL") {
		read_environment();
	} else {
		E = W;
		n_env = n_covar;
		env_names = covar_names;
	}
	if(n_env > 1) {
		std::cout << "WARNING: " << n_env << " cols detected in " << params.env_file << ". Just using first" << std::endl;
		E = E.col(0);
		n_env = 1;
		env_names.resize(1);
	}
	read_coeffs();

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
	Zg = EigenDataArrayX::Zero(n_samples);
	int ch = 0;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		if(ch % 10 == 0) {
			std::cout << "Chunk " << ch+1 << " read (size " << n_var;
			std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
			std::cout << " variants parsed)" << std::endl;
		}
		assert(n_var == G.cols());

		// Add effects to Y
		for (int kk = 0; kk < n_var; kk++) {
			Xb += G.col(kk).array() * B(kk + n_total_var, 0);
			Zg += G.col(kk).array() * E.array() * B(kk + n_total_var, 1);
		}

		n_total_var += n_var;
		ch++;
	}
	if(n_total_var != B.rows()) {
		std::cout << "ERROR: n var read in = " << n_total_var << std::endl;
		std::cout << "ERROR: n coeffs read in = " << B.rows() << std::endl;
		assert(n_total_var == B.rows());
	}
	if(n_constant_variance > 0) {
		std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
	}

	if(params.sim_w_noise) {
		// Generate noise
		double var_xb = (Xb - Xb.mean()).square().sum() / ((double) n_samples - 1.0);
		double var_zg = (Zg - Zg.mean()).square().sum() / ((double) n_samples - 1.0);
		double sigma = (var_xb + var_zg) * (1 - params.hb - params.hg);
		sigma /= (params.hb + params.hg);
		std::cout << "Checking empirical variances" << std::endl;
		std::cout << "Var(Xb) = " << var_xb << std::endl;
		std::cout << "Var(Zg) = " << var_zg << std::endl;
		std::cout << "sigma = " << sigma << std::endl;

		// http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
		std::random_device rd;
		std::mt19937 generator{rd()};
		std::normal_distribution<scalarData> standard_normal(0.0, std::sqrt(sigma));
		std::cout << "Adding white noise with variance: " << sigma << std::endl;
		noise.resize(n_samples);
		for (std::size_t ii = 0; ii < n_samples; ii++) {
			noise(ii) = standard_normal(generator);
		}

		double var_noise = (noise - noise.mean()).square().sum() / ((double) n_samples - 1.0);
		noise *= std::sqrt(sigma / var_noise);

		if(params.covar_file == "NULL") {
			Y = Xb + Zg + noise;
		} else {
			Y = W.rowwise().sum().array() + Xb + Zg + noise;
		}
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
	reduce_to_complete_cases();                             // From read_sids!!
	// std::cout << "num samples: " << n_samples << std::endl;
	while (read_bgen_chunk()) {
		// Raw dosage read in to G
		std::cout << "Chunk " << ch+1 << " read (size " << n_var;
		std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
		std::cout << " variants parsed)" << std::endl;

		for (std::size_t ii = 0; ii < n_var; ii++) {
			outf << chromosome[ii] << " " << rsid[ii] << " " << position[ii];
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
