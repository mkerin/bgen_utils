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
                                       int n_cols,
                                       std::map< int, bool > incomplete_cases );

class Data
{
public:
	parameters params;

	// Per chunk
	long n_var_chunk;
	std::vector< std::string > chromosome_chunk, rsid_chunk, SNPID_chunk;
	std::vector< uint32_t > position_chunk;
	std::vector< std::vector< std::string > > alleles_chunk;
	std::vector< std::string > SNPKEY_chunk;
	std::vector<long> valid_var_index_chunk;
	std::vector< double > info_chunk;
	std::vector< double > maf_chunk;

	// cumulative
	long n_var;
	std::vector< std::string > chromosome, rsid, SNPID;
	std::vector< uint32_t > position;
	std::vector< std::vector< std::string > > alleles;
	std::vector<long> valid_var_index;
	std::vector< double > maf;

	long n_covar;
	long n_samples;
	long n_env;
	long n_var_parsed;
	long n_constant_variance;
	long n_gxe_components, n_gxe_components2;
	bool bgen_pass;

	std::unordered_map<long, bool> sample_is_invalid;

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
	EigenDataMatrix E, eta;
	EigenDataMatrix B, B2;
	EigenDataArrayX Xb, Zg, Xb2, Zg2, Ealpha, Wtau, noise;
	EigenDataArrayXX Xg, Xg2;
	EigenDataMatrix env_profile;

	genfile::bgen::View::UniquePtr bgenView;

	boost_io::filtering_ostream outf;

	std::chrono::system_clock::time_point start;

	double s_x, s_z;

	Data( const parameters& my_params );

	~Data();

	void sim_pheno();
	void pred_pheno();
	void calc_ssv();
	void print_keys();

	void output_init();
	void output_results();
	bool read_bgen_chunk();
	void gen_genetic_effects();

	/*** Compute dXtEEX ***/
	void compute_correlations();
	void compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk);
	void output_correlations(EigenRefDataArrayXX vec);

	/*** Text file input ***/
	void read_incl_sids();
	void read_covar();
	void read_environment();
	void read_env_profile();
	EigenDataMatrix sort_coeffs(EigenRefDataMatrix orig,
	                            const std::unordered_map<std::string, long>& key_map,
	                            std::string key_type);

	void read_incl_rsids(){
		read_rsids(params.incl_rsids_file, incl_rsid_list);
	}

	void read_excl_rsids(){
		read_rsids(params.excl_rsids_file, excl_rsid_list);
	}

	void read_rsids(const std::string& filename,
	                std::vector< std::string >& rsid_list);

	void read_coeffs2(){
		read_coeffs_file(params.coeffs2_file, B2, n_gxe_components2);
	}

	void read_coeffs(){
		read_coeffs_file(params.coeffs_file, B, n_gxe_components);
	}

	void read_coeffs_file(std::string filename, EigenDataMatrix& coeffs_mat,
	                      long&n_gxe_components);

	void reduce_to_complete_cases();

	void read_txt_file_w_context( const std::string& filename,
	                              const int& col_offset,
	                              EigenDataMatrix & M,
	                              std::vector<std::string>& M_snpids,
	                              std::vector<std::string>& col_names);

	void regress_first_mat_from_second(const EigenDataMatrix& A,
	                                   const std::string& Astring,
	                                   const std::vector<std::string>& A_names,
	                                   EigenDataMatrix& yy,
	                                   const std::string& yy_string);

	void dump_coeffs(const std::string& filename,
	                 const std::vector<std::string>& header,
	                 Eigen::Ref<Eigen::MatrixXd> coeffs,
	                 Eigen::Ref<Eigen::MatrixXd> coeffs2);
};



#endif
