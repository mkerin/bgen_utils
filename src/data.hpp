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
#include <cmath>
#include <cstddef>     // for ptrdiff_t class
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

inline size_t numRows(const Eigen::MatrixXd &A);
inline size_t numCols(const Eigen::MatrixXd &A);
inline void setCol(Eigen::MatrixXd &A, const Eigen::VectorXd &v, size_t col);
inline Eigen::VectorXd getCol(const Eigen::MatrixXd &A, size_t col);
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

int n_pheno;                      // number of phenotypes
int n_covar;                      // number of covariates
int n_samples;                    // number of samples
int n_env;
long n_total_var;
bool bgen_pass;
long n_var;
long n_var_parsed;         // Track progress through IndexQuery
long int n_constant_variance;

bool Y_reduced;           // Variables to track whether we have already
bool W_reduced;           // reduced to complete cases or not.
bool E_reduced;           // reduced to complete cases or not.
	std::unordered_map<long, bool> sample_is_invalid;

std::vector< double > info;
std::vector< double > maf;
std::vector< double > maf_cum;
std::vector< std::string > incl_rsid_list;
std::vector< std::string > excl_rsid_list;

std::map<long, bool> missing_covars;         // set of subjects missing >= 1 covariate
std::map<long, bool> missing_phenos;         // set of subjects missing >= phenotype
std::map<long, bool> missing_envs;         // set of subjects missing >= phenotype
std::map<long, bool > incomplete_cases;         // union of samples missing data

	std::vector< std::string > covar_names;
std::vector< std::string > env_names;

EigenDataMatrix G;         // probabilistic genotype matrix
EigenDataArrayX Y;         // phenotype matrix
EigenDataMatrix W;         // covariate matrix
EigenDataMatrix E;
EigenDataMatrix B, B2;         // matrix of coefficients (beta, gamma)
EigenDataArrayX Xb, Zg, Xb2, Zg2, Ealpha, Wtau, noise;
EigenDataMatrix env_profile;

// Matching snps via the SNPKEY
std::unordered_map<std::string, long> B_SNPKEYS_map;
std::vector<std::string> B_SNPKEYS;
bool match_snpkeys;
std::unordered_map<std::string, long> B_SNPIDS_map;
std::vector<std::string> B_SNPIDS;
bool match_snpids;

genfile::bgen::View::UniquePtr bgenView;

boost_io::filtering_ostream outf, outf_pred, outf_coeffs;

std::chrono::system_clock::time_point start;

// Sum of sample variances
double s_x, s_z;


// constructors/destructors
// data() : bgenView( "NULL" ) {
//  bgen_pass = false; // No bgen file set; read_bgen_chunk won't run.
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

	match_snpkeys = false;                  // Used when reconstructing yhat from coeffs

	// system time at start
	start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout << "Starting analysis at " << std::ctime(&start_time) << std::endl;
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

void output_init() {
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
		outf_pred << "Xbeta\tZgamma" << std::endl;
	} else if(params.mode_gen_pheno || params.mode_gen2_pheno) {
		std::string ofile_pred   = fstream_init(outf_pred, params.out_file, "_predicted_effects");

		std::cout << "Writing predicted effects to " << ofile_pred << std::endl;
		outf_pred << "Wtau\tEalpha\tXbeta\tZgamma\tnoise" << std::endl;

		if(params.sim_w_noise) {
			std::cout << "Writing simulated phenotype to " << ofile << std::endl;
			outf << "y" << std::endl;
		}

		if(params.rescale_coeffs) {
			std::string ofile_coeffs = fstream_init(outf_coeffs, params.out_file, "_true_rescaled_coeffs");
			std::cout << "Writing rescaled coeffs to " << ofile_coeffs << std::endl;
			outf_coeffs << "chr rsid pos a0 a1 beta gamma" << std::endl;
		}
	} else if (params.mode_compute_correlations) {
		std::cout << "Writing snp-environment correlations to " << ofile << std::endl;
	} else if (params.mode_print_keys) {
		std::cout << "Writing snp-keys to " << ofile << std::endl;
		outf << "SNPID chr rsid pos a0 a1 maf info" << std::endl;
	}
}

void output_results() {
	if(params.mode_ssv) {
		outf << s_x << "\t" << s_z << "\t" << n_samples << "\t";
		outf << n_total_var << std::endl;
	} else if (params.mode_pred_pheno) {
		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf_pred << Xb(ii) << "\t" << Zg(ii) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf << Y(ii) << std::endl;
		}
	} else if(params.mode_gen_pheno || params.mode_gen2_pheno) {
		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf_pred << Wtau(ii) <<"\t" << Ealpha(ii) <<"\t" << Xb(ii) << "\t" << Zg(ii) << "\t" << noise(ii) << std::endl;
		}

		for (std::size_t ii = 0; ii < n_samples; ii++) {
			outf << Y(ii) << std::endl;
		}

		if(params.rescale_coeffs) {
			for (std::size_t ii = 0; ii < n_total_var; ii++) {
				outf_coeffs << chromosome_cum[ii] << " " << rsid_cum[ii] << " " << position_cum[ii];
				outf_coeffs << " " << alleles_cum[ii][0] << " " << alleles_cum[ii][1] << " " << B(ii, 0) << " " << B(ii, 1) << std::endl;
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
	std::string chr_j;
	uint32_t pos_j;
	std::string rsid_j;
	std::vector< std::string > alleles_j;
	std::string SNPID_j;                  // read but ignored

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
		if(!std::isfinite(mu) || !std::isfinite(sigma)){
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
			i++;                                 // loop should end at i == n_samples
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

void read_covar( ){
	// Read covariates to Eigen matrix W
	if ( params.covar_file != "NULL" ) {
		read_txt_file( params.covar_file, W, n_covar, covar_names, missing_covars );
	} else {
		throw std::logic_error( "Tried to read NULL covar file." );
	}
	W_reduced = false;
}

void read_coeffs(){
	// Read coefficients to eigen matrix B
	std::vector<std::string> case1 = {"beta", "gamma"};
	std::vector<std::string> case1b = {"beta"};

	std::vector<std::string> case2 = {"SNPKEY", "beta", "gamma"};
	std::vector<std::string> case2b = {"SNPKEY", "beta"};

	std::vector<std::string> case3 = {"SNPID", "beta", "gamma"};
	std::vector<std::string> case3b = {"SNPID", "beta"};

	std::vector<std::string> coeff_names;
	read_file_header(params.coeffs_file, coeff_names);

	if(coeff_names == case1 || coeff_names == case1b) {
		read_grid_file(params.coeffs_file, B, coeff_names);

		if (coeff_names == case1) assert(B.cols() == 2);
		if (coeff_names == case1b) assert(B.cols() == 1);
	} else if (coeff_names == case2 || coeff_names == case2b) {
		match_snpkeys = true;
		read_txt_file_w_context(params.coeffs_file, 1, B, B_SNPKEYS,
		                        coeff_names);
		for (long jj = 0; jj < B_SNPKEYS.size(); jj++) {
			B_SNPKEYS_map[B_SNPKEYS[jj]] = jj;
		}

		if (coeff_names == case2) assert(B.cols() == 2);
		if (coeff_names == case2b) assert(B.cols() == 1);
	} else if (coeff_names == case3 || coeff_names == case3b) {
		match_snpids = true;
		read_txt_file_w_context(params.coeffs_file, 1, B, B_SNPIDS,
								coeff_names);
		for (long jj = 0; jj < B_SNPIDS.size(); jj++) {
			B_SNPIDS_map[B_SNPIDS[jj]] = jj;
		}

		if (coeff_names == case3) assert(B.cols() == 2);
		if (coeff_names == case3b) assert(B.cols() == 1);
	} else {
		throw std::logic_error("Unexpected header in --coeffs");
	}
}

void read_coeffs2( ){
	std::vector<std::string> coeff_names;
	read_grid_file( params.coeffs2_file, B2, coeff_names );
	std::vector<std::string> true_names = {"beta", "gamma"};
	assert(B2.cols() == 2);
	assert(true_names == coeff_names);
	std::cout << B2.rows() << " coefficients read in from ";
	std::cout << params.coeffs2_file << std::endl;
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

void calc_ssv();

void pred_pheno();

void sim_pheno();

//void gen2_pheno();

void compute_correlations();

void regress_first_mat_from_second(const EigenDataMatrix& A,
                                   const std::string& Astring,
                                   const std::vector<std::string>& A_names,
                                   EigenDataMatrix& yy,
                                   const std::string& yy_string);

void compute_correlations_chunk(EigenRefDataArrayXX dXtEEX_chunk);

void print_keys();
};



#endif
