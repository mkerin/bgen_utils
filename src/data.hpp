// File of Data class for use with src/bgen_prog.cpp
#ifndef DATA_H
#define DATA_H

#include "parameters.hpp"
#include "utils.hpp"

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
#include <cstddef>
#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <set>
#include <stdexcept>

namespace boost_io = boost::iostreams;

inline Eigen::MatrixXd solve(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b);
inline scalarData var(const EigenDataArrayX& vv);
void sim_gaussian_noise(EigenRefDataArrayX vv, const double& sigma_sq, std::mt19937& generator);
template <typename EigenMat>
EigenMat reduce_mat_to_complete_cases( EigenMat& M,
                                       bool& matrix_reduced,
                                       int n_cols,
                                       std::map< int, bool > incomplete_cases );

class Data
{
public:
	parameters params;

// Per chunk
	std::vector< std::string > chromosome, rsid, SNPID;
	std::vector< uint32_t > position;
	std::vector< std::vector< std::string > > alleles;
	std::vector< std::string > SNPKEYS;

// cumulative
	std::vector< std::string > chromosome_cum, rsid_cum, SNPID_cum;
	std::vector< uint32_t > position_cum;
	std::vector< std::vector< std::string > > alleles_cum;

	int n_pheno;
	int n_covar;
	int n_samples;
	int n_env;
	long n_total_var;
	bool bgen_pass;
	long n_var;
	long n_var_parsed;
	long int n_constant_variance;

	bool Y_reduced;
	bool W_reduced;
	bool E_reduced;
	std::unordered_map<long, bool> sample_is_invalid;

	std::vector< double > info;
	std::vector< double > maf;
	std::vector< double > maf_cum;
	std::vector< std::string > incl_rsid_list;
	std::vector< std::string > excl_rsid_list;

	std::map<long, bool> missing_covars;
	std::map<long, bool> missing_phenos;
	std::map<long, bool> missing_envs;
	std::map<long, bool > incomplete_cases;

	std::vector< std::string > covar_names;
	std::vector< std::string > env_names;

	EigenDataMatrix G;
	EigenDataArrayX Y;
	EigenDataMatrix W;
	EigenDataMatrix E;
	EigenDataMatrix B, B2;
	EigenDataArrayX Xb, Xg, Zg, Xb2, Xg2, Zg2, Ealpha, Wtau, noise;
	EigenDataMatrix env_profile;

	std::unordered_map<std::string, long> B_SNPKEYS_map;
	std::vector<std::string> B_SNPKEYS;
	bool match_snpkeys;
	std::unordered_map<std::string, long> B_SNPIDS_map;
	std::vector<std::string> B_SNPIDS;
	bool match_snpids;

	genfile::bgen::View::UniquePtr bgenView;

	boost_io::filtering_ostream outf, outf_pred, outf_coeffs;

	std::chrono::system_clock::time_point start;

	double s_x, s_z;


	Data( const parameters& my_params );

	~Data();

	void calc_ssv();

	void pred_pheno();

	void sim_pheno();

	void compute_correlations();

	void print_keys();


	void output_init();

	void output_results();

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

	bool read_bgen_chunk();

	void read_incl_rsids(){
		read_rsids(params.incl_rsids_file, incl_rsid_list);
	}

	void read_excl_rsids(){
		read_rsids(params.excl_rsids_file, excl_rsid_list);
	}

	void read_rsids(const std::string& filename,
	                std::vector< std::string >& rsid_list);

	void read_covar( ){
		// Read covariates to Eigen matrix W
		if ( params.covar_file != "NULL" ) {
			read_txt_file( params.covar_file, W, n_covar, covar_names, missing_covars );
		} else {
			throw std::logic_error( "Tried to read NULL covar file." );
		}
		W_reduced = false;
	}

	void read_coeffs2(){
		read_coeffs_file(params.coeffs2_file, B2);
	}

	void read_coeffs(){
		read_coeffs_file(params.coeffs_file, B);
	}

	void read_coeffs_file(std::string filename, EigenDataMatrix& coeffs_mat){
		// Read coefficients to eigen matrix B
		std::vector<std::string> case1 = {"beta", "gamma"};
		std::vector<std::string> case1b = {"beta"};

		std::vector<std::string> case2 = {"SNPKEY", "beta", "gamma"};
		std::vector<std::string> case2b = {"SNPKEY", "beta"};

		std::vector<std::string> case3 = {"SNPID", "beta", "gamma"};
		std::vector<std::string> case3b = {"SNPID", "beta"};

		std::vector<std::string> coeff_names;
		read_file_header(filename, coeff_names);

		if(coeff_names == case1 || coeff_names == case1b) {
			read_grid_file(filename, coeffs_mat, coeff_names);

			if (coeff_names == case1) assert(coeffs_mat.cols() == 2);
			if (coeff_names == case1b) assert(coeffs_mat.cols() == 1);
		} else if (coeff_names == case2 || coeff_names == case2b) {
			match_snpkeys = true;
			read_txt_file_w_context(filename, 1, coeffs_mat, B_SNPKEYS, coeff_names);
			B_SNPKEYS_map.clear();
			for (long jj = 0; jj < B_SNPKEYS.size(); jj++) {
				B_SNPKEYS_map[B_SNPKEYS[jj]] = jj;
			}

			if (coeff_names == case2) assert(coeffs_mat.cols() == 2);
			if (coeff_names == case2b) assert(coeffs_mat.cols() == 1);
		} else if (coeff_names == case3 || coeff_names == case3b) {
			match_snpids = true;
			read_txt_file_w_context(filename, 1, coeffs_mat, B_SNPIDS, coeff_names);
			B_SNPIDS_map.clear();
			for (long jj = 0; jj < B_SNPIDS.size(); jj++) {
				B_SNPIDS_map[B_SNPIDS[jj]] = jj;
			}

			if (coeff_names == case3) assert(coeffs_mat.cols() == 2);
			if (coeff_names == case3b) assert(coeffs_mat.cols() == 1);
		} else {
			throw std::logic_error("Unexpected header in " + filename);
		}
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

	void read_env_profile(){
		// Read covariates to Eigen matrix W
		std::vector<std::string> placeholder;
		read_grid_file(params.env_profile_file, env_profile, placeholder);
		assert(env_profile.cols() == 1);
		// assert(env_profile.rows() == n_env);
	}

	void reduce_to_complete_cases();


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

	void read_txt_file( const std::string& filename,
	                    EigenDataMatrix& M,
	                    int& n_cols,
	                    std::vector< std::string >& col_names,
	                    std::map<long, bool >& incomplete_row ){
		// pass top line of txt file filename to col_names, and body to M.
		// TODO: Implement how to deal with missing values.

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
						} catch (const std::invalid_argument &exc) {
							std::cout << s << " on line " << i << std::endl;
							throw;
						}
					}

					if(tmp_d != params.missing_code) {
						M(i, k) = tmp_d;
					} else {
						M(i, k) = params.missing_code;
						incomplete_row[i] = true;
					}
				}
				i++;
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
	                     std::vector< std::string >& col_names);


	void read_txt_file_w_context( const std::string& filename,
	                              const int& col_offset,
	                              EigenDataMatrix & M,
	                              std::vector<std::string>& M_snpids,
	                              std::vector<std::string>& col_names);

	void center_matrix( EigenDataMatrix & M,
	                    int& n_cols );

	void scale_matrix( EigenDataMatrix & M,
	                   int& n_cols,
	                   std::vector< std::string >& col_names);

	void scale_matrix_conserved( EigenDataMatrix& M,
	                             int& n_cols);

	void regress_first_mat_from_second(const EigenDataMatrix& A,
	                                   const std::string& Astring,
	                                   const std::vector<std::string>& A_names,
	                                   EigenDataMatrix& yy,
	                                   const std::string& yy_string);

	void compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk);

	std::string fstream_init(boost_io::filtering_ostream& my_outf,
	                         const std::string& orig_file_path,
	                         const std::string& file_suffix);
};



#endif

