// File of Data class for use with src/bgen_prog.cpp
#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstddef>     // for ptrdiff_t class
#include <chrono>      // start/end time info
#include <ctime>       // start/end time info
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <string>
#include <stdexcept>
#include "class.h"
#include "utils.hpp"
#include "tools/eigen3.3/Dense"
#include "tools/eigen3.3/Sparse"
#include "tools/eigen3.3/Eigenvalues"

#include "bgen_parser.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace boost_io = boost::iostreams;

inline Eigen::MatrixXd getCols(const Eigen::MatrixXd &X, const std::vector<size_t> &cols);
inline void setCols(Eigen::MatrixXd &X, const std::vector<size_t> &cols, const Eigen::MatrixXd &values);
inline size_t numRows(const Eigen::MatrixXd &A);
inline size_t numCols(const Eigen::MatrixXd &A);
inline void setCol(Eigen::MatrixXd &A, const Eigen::VectorXd &v, size_t col);
inline Eigen::VectorXd getCol(const Eigen::MatrixXd &A, size_t col);
inline Eigen::MatrixXd solve(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b);


class Data
{
	public :
	parameters params;

	// Per chunk
	std::vector< std::string > chromosome, rsid, SNPID;
	std::vector< uint32_t > position;
	std::vector< std::vector< std::string > > alleles;
	std::vector< std::string > keys;

	// cumulative
	std::vector< std::string > chromosome_cum, rsid_cum, SNPID_cum;
	std::vector< uint32_t > position_cum;
	std::vector< std::vector< std::string > > alleles_cum;
	std::vector< std::string > keys_cum;

	int n_pheno;              // number of phenotypes
	int n_covar;              // number of covariates
	int n_samples;            // number of samples
	int n_env;
	std::size_t n_total_var;
	bool bgen_pass;
	int n_var;
	std::size_t n_var_parsed; // Track progress through IndexQuery
	long int n_constant_variance;

	bool Y_reduced;   // Variables to track whether we have already
	bool W_reduced;   // reduced to complete cases or not.
	bool E_reduced;   // reduced to complete cases or not.

	std::vector< double > info;
	std::vector< double > maf;
	std::vector< double > maf_cum;
	std::vector< std::string > incl_rsid_list;
	std::vector< std::string > excl_rsid_list;

	std::map<int, bool> missing_covars; // set of subjects missing >= 1 covariate
	std::map<int, bool> missing_phenos; // set of subjects missing >= phenotype
	std::map<int, bool> missing_envs; // set of subjects missing >= phenotype
	std::map< int, bool > incomplete_cases; // union of samples missing data

	std::vector< std::string > pheno_names;
	std::vector< std::string > covar_names;
	std::vector< std::string > env_names;

	EigenDataMatrix G; // probabilistic genotype matrix
	EigenDataArrayX Y; // phenotype matrix
	EigenDataMatrix W; // covariate matrix
	EigenDataMatrix E;
	EigenDataMatrix B; // matrix of coefficients (beta, gamma)
	EigenDataArrayX Xb, Zg, Ealpha, Wtau, noise;
	genfile::bgen::View::UniquePtr bgenView;
	std::vector< double > beta, tau, neglogP, neglogP_2dof;
	std::vector< std::vector< double > > gamma;

	boost_io::filtering_ostream outf, outf_pred, outf_coeffs;

	std::chrono::system_clock::time_point start;

	// grid things for vbayes
	std::vector< std::string > hyps_names, imprt_names;
	Eigen::MatrixXd hyps_grid, imprt_grid;
	Eigen::VectorXd alpha_init, mu_init;

	// Sum of sample variances
	double s_x, s_z;


	// constructors/destructors
	// data() : bgenView( "NULL" ) {
	// 	bgen_pass = false; // No bgen file set; read_bgen_chunk won't run.
	// }

	Data( const parameters& my_params ) : params(my_params){

		bgenView = genfile::bgen::View::create(my_params.bgen_file);
		bgen_pass = true;
		n_samples = bgenView->number_of_samples();
		n_var_parsed = 0;
		n_total_var = 0;
		n_constant_variance = 0;
		n_covar = 0;
		n_env = 0;

		// system time at start
		start = std::chrono::system_clock::now();
		std::time_t start_time = std::chrono::system_clock::to_time_t(start);
		std::cout << "Starting analysis at " << std::ctime(&start_time) << std::endl;
		std::cout << "Compiled from git branch: master" << std::endl;
	}

	~Data() {
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

	std::string fstream_init(boost_io::filtering_ostream& my_outf,
                             const std::string& orig_file_path,
                             const std::string& file_suffix){

		std::string filepath   = orig_file_path;
		std::string dir        = filepath.substr(0, filepath.rfind("/")+1);
		std::string stem_w_dir = filepath.substr(0, filepath.find("."));
		std::string stem       = stem_w_dir.substr(stem_w_dir.rfind("/")+1, stem_w_dir.size());
		std::string ext        = filepath.substr(filepath.find("."), filepath.size());

		std::string ofile      = dir + stem + file_suffix + ext;

		my_outf.reset();
		std::string gz_str = ".gz";
		if (orig_file_path.find(gz_str) != std::string::npos) {
			my_outf.push(boost_io::gzip_compressor());
		}
		my_outf.push(boost_io::file_sink(ofile.c_str()));
		return ofile;
	}

	void output_init() {
		std::string ofile      = fstream_init(outf, params.out_file, "");


		// Output header for ssv file
		if(params.mode_ssv){
			std::cout << "Writing ssv stats to " << ofile << std::endl;
			outf << "s_x\ts_z\tn\tp" << std::endl;
		} else if(params.mode_gen_pheno || params.mode_gen2_pheno){
			std::string ofile_pred   = fstream_init(outf_pred, params.out_file, "_predicted_effects");

			std::cout << "Writing predicted effects to " << ofile_pred << std::endl;
			outf_pred << "Wtau\tEalpha\tXbeta\tZgamma\tnoise" << std::endl;

			if(params.sim_w_noise){
				std::cout << "Writing simulated phenotype to " << ofile << std::endl;
				outf << "y" << std::endl;
			}

			if(params.rescale_coeffs){
				std::string ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_true_rescaled_coeffs");
				std::cout << "Writing rescaled coeffs to " << ofile_coeffs << std::endl;
				outf_coeffs << "chr rsid pos a0 a1 beta gamma" << std::endl;
			}
		} else if (params.mode_compute_correlations){
			std::cout << "Writing snp-environment correlations to " << ofile << std::endl;
		} else if (params.mode_print_keys){
			std::cout << "Writing snp-keys to " << ofile << std::endl;
			outf << "chr rsid pos a0 a1 maf info" << std::endl;
		}
	}

	void output_results() {
		if(params.mode_ssv){
			outf << s_x << "\t" << s_z << "\t" << n_samples << "\t";
			outf << n_total_var << std::endl;
		} else if(params.mode_gen_pheno || params.mode_gen2_pheno){
			for (std::size_t ii = 0; ii < n_samples; ii++){
				outf_pred << Wtau(ii) <<"\t" << Ealpha(ii) <<"\t" << Xb(ii) << "\t" << Zg(ii) << "\t" << noise(ii) << std::endl;
			}

			if(params.sim_w_noise){
				for (std::size_t ii = 0; ii < n_samples; ii++){
					outf << Y(ii) << std::endl;
				}
			}

			if(params.rescale_coeffs){
				for (std::size_t ii = 0; ii < n_total_var; ii++){
					outf_coeffs << chromosome_cum[ii] << " " << rsid_cum[ii] << " " << position_cum[ii];
					outf_coeffs << " " << alleles_cum[ii][0] << " " << alleles_cum[ii][1] << " " << B(ii, 0) << " " << B(ii, 1) << std::endl;
				}
			}
		} else if (params.mode_print_keys){
			// for (std::size_t ii = 0; ii < n_total_var; ii++){
			// 	outf << chromosome_cum[ii] << " " << rsid_cum[ii] << " " << position_cum[ii];
			// 	outf << " " << alleles_cum[ii][0] << " " << alleles_cum[ii][1];
 			// 	outf << " " << maf_cum[ii] << " " << info_cum[ii] << std::endl;
			// }
		}
	}

	void output_correlations(EigenRefDataArrayXX vec){
		outf << std::scientific << std::setprecision(8);
		for (long jj = 0; jj < n_var; jj++) {
			outf << SNPID[jj] << " " << chromosome[jj] << " " << rsid[jj] << " " << position[jj];
			outf << " " << alleles[jj][0] << " " << alleles[jj][1];
			for (int ll = 0; ll < n_env * n_env; ll++) {
				outf << " " << vec(jj, ll);
			}
			outf << std::endl;
		}
	}

	bool read_bgen_chunk() {
		// Wrapper around BgenView to read in a 'chunk' of data. Remembers
		// if last call hit the EOF, and returns false if so.
		// Assumed that:
		// - commandline args parsed and passed to params
		// - bgenView initialised with correct filename
		// - scale + centering happening internally

		// Exit function if last call hit EOF.
		if (!bgen_pass) return false;

		// Temporary variables to store info from read_variant()
		std::string chr_j ;
		uint32_t pos_j ;
		std::string rsid_j ;
		std::vector< std::string > alleles_j ;
		std::string SNPID_j ; // read but ignored
		std::vector< std::vector< double > > probs ;
		ProbSetter setter( &probs );

		double x, dosage, check, info_j, f1, chunk_missingness;
		double missing_calls = 0.0;
		int n_var_incomplete = 0;

		// Wipe variant context from last chunk
		maf.clear();
		info.clear();
		rsid.clear();
		chromosome.clear();
		position.clear();
		alleles.clear();
		keys.clear();
		SNPID.clear();

		// Resize genotype matrix
		G.resize(n_samples, params.chunk_size);

		std::size_t valid_count, jj = 0;
		while ( jj < params.chunk_size && bgen_pass ) {
			bgen_pass = bgenView->read_variant( &SNPID_j, &rsid_j, &chr_j, &pos_j, &alleles_j );
			if (!bgen_pass) break;
			n_var_parsed++;
			assert( alleles_j.size() > 0 );

			// Read probs + check maf filter
			bgenView->read_genotype_data_block( setter );

			// maf + info filters; computed on valid sample_ids & variants whose alleles
			// sum to 1
			std::map<int, bool> missing_genos;
			EigenDataArrayX dosage_j(n_samples);
			double f2 = 0.0;
			double valid_count = 0;
			std::uint32_t ii_obs = 0;
			for( std::size_t ii = 0; ii < probs.size(); ++ii ) {
				if (incomplete_cases.count(ii) == 0) {
					f1 = dosage = check = 0.0;
					for( std::size_t kk = 0; kk < probs[ii].size(); kk++ ) {
						x = probs[ii][kk];
						dosage += x * kk;
						f1 += x * kk * kk;
						check += x;
					}
					if(check > 0.9999 && check < 1.0001){
						dosage_j[ii_obs] = dosage;
						f2 += (f1 - dosage * dosage);
						valid_count++;
					} else {
						missing_genos[ii_obs] = 1;
						dosage_j[ii_obs] = 0.0;
					}
					ii_obs++;
				}
			}
			assert(ii_obs == n_samples);

			double d1    = dosage_j.sum();
			double maf_j = d1 / (2.0 * valid_count);

			// Flip dosage vector if maf > 0.5
			if(maf_j > 0.5){
				dosage_j = (2.0 - dosage_j);
				for (std::uint32_t ii = 0; ii < n_samples; ii++){
					if (missing_genos.count(ii) > 0){
						dosage_j[ii] = 0.0;
					}
				}

				f2       = 4.0 * valid_count - 4.0 * d1 + f2;
				d1       = dosage_j.sum();
				maf_j    = d1 / (2.0 * valid_count);
			}

			double mu    = d1 / valid_count;
			info_j       = 1.0;
			if(maf_j > 1e-10 && maf_j < 0.9999999999){
				info_j -= f2 / (2.0 * valid_count * maf_j * (1.0 - maf_j));
			}

			// Compute sd
			double sigma = (dosage_j - mu).square().sum();
			sigma = std::sqrt(sigma/(valid_count - 1.0));

			// Filters
			if (params.maf_lim && (maf_j < params.min_maf || maf_j > 1 - params.min_maf)) {
				continue;
			}
			if (params.info_lim && info_j < params.min_info) {
				continue;
			}
			if(!params.keep_constant_variants && d1 < 5.0){
				n_constant_variance++;
				continue;
			}
			if(!params.keep_constant_variants && sigma <= 1e-12){
				n_constant_variance++;
				continue;
			}

			// filters passed; write contextual info
			maf.push_back(maf_j);
			info.push_back(info_j);
			rsid.push_back(rsid_j);
			chromosome.push_back(chr_j);
			position.push_back(pos_j);
			alleles.push_back(alleles_j);
			SNPID.push_back(SNPID_j);

			std::string key_j = chr_j + "~" + std::to_string(pos_j) + "~" + alleles_j[0] + "~" + alleles_j[1];
			keys.push_back(key_j);
			if(params.mode_print_keys){
				jj++;
				continue;
			}

			// Scale + center dosage, set missing to mean
			if(missing_genos.size() > 0){
				n_var_incomplete += 1;
				missing_calls += (double) missing_genos.size();
				for (std::uint32_t ii = 0; ii < n_samples; ii++){
					if (missing_genos.count(ii) > 0){
						dosage_j[ii] = mu;
					}
				}
			}

			if(params.mode_low_mem){
				double L = 2.0;
				double intervalWidth = L / 256.0;
				double invIntervalWidth = 256.0 / L;

				// Compress
				for (std::uint32_t ii = 0; ii < n_samples; ii++){
					dosage_j[ii] = std::floor(std::min(dosage_j[ii], (scalarData) (L - 1e-6)) * invIntervalWidth);
				}

				// Decompress
				dosage_j = (dosage_j + 0.5) * intervalWidth;

				mu    = dosage_j.mean();
				sigma = (dosage_j - mu).square().sum();
				sigma = std::sqrt(sigma/(valid_count - 1.0));
			}

			dosage_j -= mu;
			dosage_j /= sigma;
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

		chunk_missingness = missing_calls / (double) (n_var * n_samples);
		if(chunk_missingness > 0.0){
			std::cout << "Chunk missingness " << chunk_missingness << "(";
 			std::cout << n_var_incomplete << "/" << n_var;
			std::cout << " variants incomplete)" << std::endl;
		}

		if(jj == 0){
			// Immediate EOF
			return false;
		} else {
			return true;
		}
	}

	void read_incl_rsids(){
		read_rsids(params.incl_rsids_file, incl_rsid_list);
	}

	void read_excl_rsids(){
		read_rsids(params.excl_rsids_file, excl_rsid_list);
	}

	void read_rsids(const std::string& filename,
                    std::vector< std::string >& rsid_list){
		rsid_list.clear();
		boost_io::filtering_istream fg;
		fg.push(boost_io::file_source(params.incl_rsids_file.c_str()));
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

	void read_incl_sids(){
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
			[&]( std::string const& id ) { bgen_ids.push_back(id); }
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
		for (long ii = 0; ii < n_samples; ii++){
			it = user_sample_ids.find(bgen_ids[ii]);
			if (it == user_sample_ids.end()){
				incomplete_cases[ii] = true;
			}
		}

		std::cout << "Subsetted down to " << user_sample_ids.size() << " ids from --incl_sample_ids";
		std::cout << std::endl;
	}

	void read_txt_file( std::string filename,
						EigenDataMatrix& M,
						int& n_cols,
						std::vector< std::string >& col_names,
						std::map< int, bool >& incomplete_row ){
		// pass top line of txt file filename to col_names, and body to M.
		// TODO: Implement how to deal with missing values.

		boost_io::filtering_istream fg;
		fg.push(boost_io::file_source(filename.c_str()));
		if (!fg) {
			std::cout << "ERROR: " << filename << " not opened." << std::endl;
			std::exit(EXIT_FAILURE);
		}

		// Reading column names
		std::string line;
		if (!getline(fg, line)) {
			std::cout << "ERROR: " << filename << " not read." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		std::stringstream ss;
		std::string s;
		n_cols = 0;
		ss.clear();
		ss.str(line);
		while (ss >> s) {
			++n_cols;
			col_names.push_back(s);
		}
		std::cout << " Reading matrix of size " << n_samples << " x " << n_cols << " from " << filename << std::endl;

		// Write remainder of file to Eigen matrix M
		incomplete_row.clear();
		M.resize(n_samples, n_cols);
		int i = 0;
		double tmp_d;
		try {
			while (getline(fg, line)) {
				if (i >= n_samples) {
					throw std::runtime_error("ERROR: could not convert txt file (too many lines).");
				}
				ss.clear();
				ss.str(line);
				for (int k = 0; k < n_cols; k++) {
					std::string s;
					ss >> s;
					/// NA
					if (s == "NA" || s == "NAN" || s == "NaN" || s == "nan") {
						tmp_d = params.missing_code;
					} else {
						try{
							tmp_d = stod(s);
						} catch (const std::invalid_argument &exc){
							std::cout << s << " on line " << i << std::endl;
							throw;
						}
					}

					if(tmp_d != params.missing_code) {
						M(i, k) = tmp_d;
					} else {
						M(i, k) = params.missing_code;
						incomplete_row[i] = 1;
					}
				}
				i++; // loop should end at i == n_samples
			}
			if (i < n_samples) {
				throw std::runtime_error("ERROR: could not convert txt file (too few lines).");
			}
		} catch (const std::exception &exc) {
			throw;
		}
	}

	void read_grid_file( const std::string& filename,
						 EigenDataMatrix & M,
						 std::vector< std::string >& col_names){
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
					} catch (const std::invalid_argument &exc){
						std::cout << s << " on line " << i << std::endl;
						throw;
					}

					M(i, k) = tmp_d;
				}
				i++; // loop should end at i == n_grid
			}
			if (i < n_grid) {
				throw std::runtime_error("ERROR: could not convert txt file (too few lines).");
			}
		} catch (const std::exception &exc) {
			throw;
		}
	}

	void center_matrix( EigenDataMatrix & M,
						int& n_cols ){
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

	void scale_matrix( EigenDataMatrix & M,
						int& n_cols,
 						std::vector< std::string >& col_names){
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
			// 	std::cout << reject_names[kk] << std::endl;
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

	void scale_matrix_conserved( EigenDataMatrix& M,
						int& n_cols){
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

	void read_covar( ){
		// Read covariates to Eigen matrix W
		if ( params.covar_file != "NULL" ) {
			read_txt_file( params.covar_file, W, n_covar, covar_names, missing_covars );
		} else {
			throw std::logic_error( "Tried to read NULL covar file." );
		}
		W_reduced = false;
	}

	void read_coeffs( ){
		// Read coefficients to eigen matrix B
		std::vector<std::string> coeff_names;
		if ( params.coeffs_file != "NULL" ) {
			read_grid_file( params.coeffs_file, B, coeff_names );
		} else {
			throw std::logic_error( "Tried to read NULL covar file." );
		}
		std::vector<std::string> true_names = {"beta", "gamma"};
		assert(B.cols() == 2);
		assert(true_names == coeff_names);
		std::cout << B.rows() << " coefficients read in from ";
		std::cout << params.coeffs_file << std::endl;
	}

	void read_environment( ){
		// Read covariates to Eigen matrix W
		if ( params.env_file != "NULL" ) {
			read_txt_file( params.env_file, E, n_env, env_names, missing_envs );
		} else {
			throw std::logic_error( "Tried to read NULL env file." );
		}
		E_reduced = false;
	}

	template <typename EigenMat>
	EigenMat reduce_mat_to_complete_cases( EigenMat& M,
								   bool& matrix_reduced,
								   int n_cols,
								   std::map< int, bool > incomplete_cases ) {
		// Remove rows contained in incomplete_cases
		int n_incomplete;
		EigenMat M_tmp;
		if (matrix_reduced) {
			throw std::runtime_error("ERROR: Trying to remove incomplete cases twice...");
		}

		// Create temporary matrix of complete cases
		n_incomplete = incomplete_cases.size();
		M_tmp.resize(n_samples - n_incomplete, n_cols);

		// Fill M_tmp with non-missing entries of M
		int ii_tmp = 0;
		for (int ii = 0; ii < n_samples; ii++) {
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

	void reduce_to_complete_cases() {
		// Remove any samples with incomplete covariates or phenotypes from
		// Y and W.
		// Note; other functions (eg. read_incl_sids) may add to incomplete_cases
		// Note; during unit testing sometimes only phenos or covars present.

		incomplete_cases.insert(missing_covars.begin(), missing_covars.end());
		incomplete_cases.insert(missing_phenos.begin(), missing_phenos.end());
		incomplete_cases.insert(missing_envs.begin(), missing_envs.end());
		if(params.pheno_file != "NULL"){
			Y = reduce_mat_to_complete_cases( Y, Y_reduced, n_pheno, incomplete_cases );
		}
		if(params.covar_file != "NULL" || n_covar > 0){
			W = reduce_mat_to_complete_cases( W, W_reduced, n_covar, incomplete_cases );
		}
		if(params.env_file != "NULL" || n_env > 0){
			E = reduce_mat_to_complete_cases( E, E_reduced, n_env, incomplete_cases );
		}
		n_samples -= incomplete_cases.size();
		missing_phenos.clear();
		missing_covars.clear();
		missing_envs.clear();

		std::cout << "Reduced to " << n_samples << " samples with complete data";
 		std::cout << " across covariates and phenotype." << std::endl;
	}

	void calc_ssv(){
		// Returns sum of sample variances from columns of X and Z

		// Step 1; Read in raw covariates
		// - also makes a record of missing values
		if(params.covar_file != "NULL"){
			read_covar();
		}

		// Step 2; Reduce raw covariates and phenotypes to complete cases
		// - may change value of n_samples
		// - will also skip these cases when reading bgen later
		reduce_to_complete_cases();

		// Step 3; Normalise covars
		if(n_covar > 0){
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
			for (int kk = 0; kk < n_var; kk++){
				tmp_x = G.col(kk).dot(G.col(kk));
				s_x += tmp_x / ((double) (n_samples - 1.0));

				if(n_covar > 0){
					mean_z = G.col(kk).dot(W.col(0)) / (double) n_samples;
					tmp_z = 0.0;
					for (std::size_t ii = 0; ii < n_samples; ii++){
						tmp_z += (W(ii,0)*G(ii,kk) - mean_z)*(W(ii,0)*G(ii,kk) - mean_z);
					}
					s_z += tmp_z / ((double) (n_samples - 1.0));
				}
			}
			ch++;
		}
		if(n_constant_variance > 0){
			std::cout << " Removed " << n_constant_variance  << " column(s) variance < 1e-12 or mac < 5:" << std::endl;
		}
	}

	void gen_pheno(){
		// Gen vector of fitted effects Y_hat = age + X beta + Z gamma
		// Gaussian noise added when writing to file

		// Step 1; Read in raw covariates
		// - also makes a record of missing values
		if(params.covar_file != "NULL"){
			read_covar();
		}
		if(params.env_file != "NULL"){
			read_environment();
		}
		read_coeffs();

		if(n_env > 1){
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
		if(n_covar > 0 && !params.use_raw_covars){
			center_matrix( W, n_covar );
			scale_matrix( W, n_covar, covar_names );
		}
		if(n_env > 0 && !params.use_raw_env){
			center_matrix( E, n_env );
			scale_matrix( E, n_env, env_names );
		}

		// Generate random effects
		// http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
		// Random seed
		if(params.random_seed == -1){
			std::random_device rd;
			params.random_seed = rd();
		}
		std::cout << "Initialising random sample generator with seed " << params.random_seed << std::endl;
		std::mt19937 generator{params.random_seed};
		EigenDataVector tau(n_covar), alpha(n_env);
		noise.resize(n_samples);
		Ealpha.resize(n_samples);
		Wtau.resize(n_samples);

		std::normal_distribution<scalarData> noise_normal(0.0, std::sqrt(params.sigma));
		std::cout << "Adding white noise with variance: " << params.sigma << std::endl;
		for (std::size_t ii = 0; ii < n_samples; ii++){
			noise(ii) = noise_normal(generator);
		}

		// additive covar effects
		double sigma_c = params.sigma * params.hc / (1.0 - params.hc - params.he - params.hb - params.hg) / (double) n_covar;
		if(params.hc > 0){
			std::normal_distribution<scalarData> covar_normal(0.0, std::sqrt(sigma_c));
			for (int ee = 0; ee < n_covar; ee++){
				tau(ee) = covar_normal(generator);
			}
			Wtau = W * tau;
		} else {
			Wtau = EigenDataVector::Zero(n_samples);
		}

		// additive env effects
		double sigma_e = params.sigma * params.he / (1.0 - params.hc - params.he - params.hb - params.hg) / (double) n_env;
		if(params.he > 0){
			std::normal_distribution<scalarData> env_normal(0.0, std::sqrt(sigma_e));
			for (int ee = 0; ee < n_env; ee++){
				alpha(ee) = env_normal(generator);
			}
			Ealpha = E * alpha;
		} else {
			Ealpha = EigenDataVector::Zero(n_samples);
		}

		// S2; genetic or interaction effects
		Xb = EigenDataArrayX::Zero(n_samples);
		Zg = EigenDataArrayX::Zero(n_samples);
		int ch = 0;
		while (read_bgen_chunk()) {
			// Raw dosage read in to G
			if(ch % 10 == 0){
				std::cout << "Chunk " << ch+1 << " read (size " << n_var;
				std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
				std::cout << " variants parsed)" << std::endl;
			}

			maf_cum.insert(maf_cum.end(),               maf.begin(), maf.end());
			rsid_cum.insert(rsid_cum.end(),             rsid.begin(), rsid.end());
			chromosome_cum.insert(chromosome_cum.end(), chromosome.begin(), chromosome.end());
			position_cum.insert(position_cum.end(),     position.begin(), position.end());
			alleles_cum.insert(alleles_cum.end(),       alleles.begin(), alleles.end());

			assert(G.cols() == n_var);
			if(n_env > 0) assert(E.rows() == n_samples);
			assert(G.rows() == n_samples);
			assert(B.rows() >= n_total_var + n_var);

			// Add effects to Y
			for (int kk = 0; kk < n_var; kk++){
				Xb += G.col(kk).array() * B(kk + n_total_var, 0);
				if(n_env > 0){
					Zg += G.col(kk).array() * E.array() * B(kk + n_total_var, 1);
				}
			}

			n_total_var += n_var;
			ch++;
		}
		if(n_total_var != B.rows()){
			std::cout << "ERROR: n var read in = " << n_total_var << std::endl;
			std::cout << "ERROR: n coeffs read in = " << B.rows() << std::endl;
			assert(n_total_var == B.rows());
		}
		if(n_constant_variance > 0){
			std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
		}

		if(params.rescale_coeffs){
			// // compute variance
			double var_xb = (Xb - Xb.mean()).square().sum() / ((double) n_samples - 1.0);
			double var_noise  = (noise - noise.mean()).square().sum() / ((double) n_samples - 1.0);
			//
			// // rescaling for correct heritability
			double sf_b = std::sqrt(params.sigma * params.hb / (1.0 - params.hc - params.he - params.hb - params.hg) / var_xb);
			std::cout << "Rescaling components for correct trait variance" << std::endl;
			std::cout << "beta scale factor: " << sf_b << std::endl;
			Xb    *= sf_b;
			noise *= std::sqrt(params.sigma / var_noise);
			B.col(0) = B.col(0) * sf_b;

			if(n_env > 0){
				double sf_g;
				if(params.hg > 0){
					double var_zg = (Zg - Zg.mean()).square().sum() / ((double) n_samples - 1.0);
					sf_g = std::sqrt(params.sigma * params.hg / (1.0 - params.hc - params.he - params.hb - params.hg) / var_zg);
				} else {
					sf_g = 0.0;
				}
				std::cout << "gamma scale factor: " << sf_g << std::endl;
				Zg    *= sf_g;
				B.col(1) = B.col(1) * sf_g;
			}

			if(params.he > 0){
				double var_Ea = (Ealpha - Ealpha.mean()).square().sum() / ((double) n_samples - 1.0);
				double sf_e = std::sqrt(sigma_e * (double) n_env / var_Ea);
				std::cout << "Ealpha scale factor: " << sf_e << std::endl;
				Ealpha    *= sf_e;
			}

			if(params.hc > 0){
				double var_Wt = (Wtau - Wtau.mean()).square().sum() / ((double) n_samples - 1.0);
				double sf_c = std::sqrt(sigma_c * (double) n_covar / var_Wt);
				std::cout << "Wtau scale factor: " << sf_c << std::endl;
				Wtau    *= sf_c;

				// double var_W0 = (W.col(0) - W.col(0).mean()).square().sum() / ((double) n_samples - 1.0);
				// std::cout << "var Wt: " << var_Wt << std::endl;
				// std::cout << "var W0: " << var_W0 << std::endl;
				// std::cout << "sigma_c: " << sigma_c << std::endl;
				// std::cout << "n_covar: " << n_covar << std::endl;
			}
		}

		Y = Wtau + Ealpha + Xb + Zg + noise;
	}

	void gen2_pheno(){
		// Controlling heritability of generated phenotype with variance of
		// noise subsequently added, a la BSLMM.

		if(params.covar_file != "NULL"){
			read_covar();
		}
		if(params.env_file != "NULL"){
			read_environment();
		} else {
			E = W;
			n_env = n_covar;
			env_names = covar_names;
		}
		if(n_env > 1){
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
		if(n_covar > 0 && !params.use_raw_covars){
			center_matrix( W, n_covar );
			scale_matrix( W, n_covar, covar_names );
		}
		if(n_env > 0 && !params.use_raw_env){
			center_matrix( E, n_env );
			scale_matrix( E, n_env, env_names );
		}

		// S2; genetic or interaction effects
		Xb = EigenDataArrayX::Zero(n_samples);
		Zg = EigenDataArrayX::Zero(n_samples);
		int ch = 0;
		while (read_bgen_chunk()) {
			// Raw dosage read in to G
			if(ch % 10 == 0){
				std::cout << "Chunk " << ch+1 << " read (size " << n_var;
				std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
				std::cout << " variants parsed)" << std::endl;
			}
			assert(n_var == G.cols());

			// Add effects to Y
			for (int kk = 0; kk < n_var; kk++){
				Xb += G.col(kk).array() * B(kk + n_total_var, 0);
				Zg += G.col(kk).array() * E.array() * B(kk + n_total_var, 1);
			}

			n_total_var += n_var;
			ch++;
		}
		if(n_total_var != B.rows()){
			std::cout << "ERROR: n var read in = " << n_total_var << std::endl;
			std::cout << "ERROR: n coeffs read in = " << B.rows() << std::endl;
			assert(n_total_var == B.rows());
		}
		if(n_constant_variance > 0){
			std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
		}

		if(params.sim_w_noise){
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
			for (std::size_t ii = 0; ii < n_samples; ii++){
				noise(ii) = standard_normal(generator);
			}

			double var_noise = (noise - noise.mean()).square().sum() / ((double) n_samples - 1.0);
			noise *= std::sqrt(sigma / var_noise);

			if(params.covar_file == "NULL"){
				Y = Xb + Zg + noise;
			} else {
				Y = W.rowwise().sum().array() + Xb + Zg + noise;
			}
		}
	}

	void compute_correlations(){
		// Gen vector of fitted effects Y_hat = age + X beta + Z gamma
		// Gaussian noise added when writing to file

		// Step 1; Read in raw covariates
		// - also makes a record of missing values
		read_environment();

		// Step 2; Reduce raw covariates and phenotypes to complete cases
		// - may change value of n_samples
		// - will also skip these cases when reading bgen later
		reduce_to_complete_cases();

		// Step 3; Normalise covars
		center_matrix( E, n_env );
		scale_matrix( E, n_env, env_names );

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

		if(n_constant_variance > 0){
			std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
		}
	}

	void compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk){
		assert(dXtEEX_chunk.rows() == n_var);
		assert(dXtEEX_chunk.cols() == n_env * n_env);
		EigenDataArrayX cl_j;
		for (int jj = 0; jj < n_var; jj++){
			cl_j = G.col(jj);
			for (int ll = 0; ll < n_env; ll++){
				for (int mm = 0; mm <= ll; mm++){
					double x = (cl_j * E.array().col(ll) * E.array().col(mm) * cl_j).sum();
					dXtEEX_chunk(jj, ll*n_env + mm) = x;
					dXtEEX_chunk(jj, mm*n_env + ll) = x;
				}
			}
		}
	}

	void print_keys(){
		// Step 4; compute correlations
		int ch = 0;
		reduce_to_complete_cases(); // From read_sids!!
		// std::cout << "num samples: " << n_samples << std::endl;
		while (read_bgen_chunk()) {
			// Raw dosage read in to G
			std::cout << "Chunk " << ch+1 << " read (size " << n_var;
			std::cout << ", " << n_var_parsed-1 << "/" << bgenView->number_of_variants();
			std::cout << " variants parsed)" << std::endl;

			maf_cum.insert(maf_cum.end(),               maf.begin(), maf.end());
			rsid_cum.insert(rsid_cum.end(),             rsid.begin(), rsid.end());
			chromosome_cum.insert(chromosome_cum.end(), chromosome.begin(), chromosome.end());
			position_cum.insert(position_cum.end(),     position.begin(), position.end());
			alleles_cum.insert(alleles_cum.end(),       alleles.begin(), alleles.end());
			keys_cum.insert(keys_cum.end(),             keys.begin(), keys.end());

			for (std::size_t ii = 0; ii < n_var; ii++){
				outf << chromosome[ii] << " " << rsid[ii] << " " << position[ii];
				outf << " " << alleles[ii][0] << " " << alleles[ii][1];
 				outf << " " << maf[ii] << " " << info[ii] << std::endl;
			}

			n_total_var += n_var;
			ch++;
		}

		if(n_constant_variance > 0){
			std::cout << " Removed " << n_constant_variance  << " column(s) with zero variance:" << std::endl;
		}
	}
};


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

inline void setCols(Eigen::MatrixXd &X, const std::vector<size_t> &cols, const Eigen::MatrixXd &values) {
	assert(cols.size() == numCols(values));
	assert(numRows(X) == numRows(values));

	if (cols.size() == 0) {
		return;
	}

	for (size_t i = 0; i < cols.size(); i++) {
		setCol(X, getCol(values, i), cols[i]);
	}
}

inline Eigen::MatrixXd getCols(const Eigen::MatrixXd &X, const std::vector<size_t> &cols) {
	Eigen::MatrixXd result(numRows(X), cols.size());
	assert(cols.size() == numCols(result));
	assert(numRows(X) == numRows(result));

	if (cols.size() == 0) {
		return result;
	}

	for (size_t i = 0; i < cols.size(); i++) {
		setCol(result, getCol(X, cols[i]), i);
	}

	return result;
}

inline Eigen::MatrixXd solve(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b) {
	Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
	if (fabs((double)((A * x - b).norm()/b.norm())) > 1e-8) {
		std::cout << "ERROR: could not solve covariate scatter matrix." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return x;
}

#endif
