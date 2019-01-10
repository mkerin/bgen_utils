// parse_arguments
#ifndef PARSE_ARGUMENTS_HPP
#define PARSE_ARGUMENTS_HPP

#include <iostream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include "class.h"
#include "version.h"
#include <regex>
#include <stdexcept>

void check_counts(std::string in_str, int i, int num, int argc);
void parse_arguments(parameters &p, int argc, char *argv[]);
void check_file_exists(const std::string& filename);

void check_counts(std::string in_str, int i, int num, int argc) {
	// Stop overflow from argv
	if (i + num >= argc) {
		if (num == 1) {
			std::cout << "ERROR: flag " << in_str << " requres an argument. ";
			std::cout << "Please refer to the manual for usage instructions." << std::endl;
		} else {
			std::cout << "ERROR: flag " << in_str << " seems to require ";
 			std::cout << std::to_string(num) + " arguments. No arguments of ";
			std::cout << "this type should be implemented yet.." << std::endl;
		}
		std::exit(EXIT_FAILURE);
	}
}

void check_file_exists(const std::string& filename){
	// Throw error if given file does not exist.
	// NB: Doesn't check if file is empty etc.
	struct stat buf;
	if(stat(filename.c_str(), &buf) != 0){
		std::cout << "File " << filename << " does not exist" << std::endl;
		throw std::runtime_error("ERROR: file does not exist");
	}
}

void parse_arguments(parameters &p, int argc, char *argv[]) {
	char *in_str;
	int i;
	std::set<std::string> option_list {
		"--bgen",
		"--covar",
		"--environment",
		"--chunk",
		"--range",
		"--maf",
		"--info",
		"--out",
		"--incl_sample_ids",
		"--incl_rsids",
		"--excl_rsids",
		"--rsid",
		"--no_geno_check",
		"--mode_ssv",
		"--mode_gen_pheno",
		"--mode_gen2_pheno",
		"--coeffs",
		"--true_sigma",
		"--true_hb",
		"--true_hg",
		"--true_hc",
		"--true_he",
		"--keep_constant_variants",
		"--print-keys",
		"--mode_low_mem",
		"--use_raw_covars",
		"--use_raw_env",
		"--random_seed",
		"--compute-env-snp-correlations"
	};

	std::set<std::string>::iterator set_it;
	// for --version (& no other args), print splash screen and exit.
	if (argc == 2 && strcmp(argv[1], "--version") == 0) {
		printf("%s\n\n", splash);
		std::exit(EXIT_SUCCESS);
	}

	if (argc == 1) {
		std::cout << "======-----"<< std::endl;
		std::cout << "Matt's BGEN PROG" << std::endl;
		std::cout << "======-----" << std::endl << std::endl;
		std::exit(EXIT_SUCCESS);
	}

	// read in and check option flags
	for (i = 0; i < argc; i++) {
		in_str = argv[i];
		if (strcmp(in_str, "--version") == 0 || strcmp(in_str, "--help") == 0) {
			std::cout << "ERROR: flag '" << in_str << "' cannot be used with any other flags." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	// Ensure some arguments only appear once
	bool check_out = 0;

	for(i = 0; i < argc; i++){
		if(*argv[i] == '-'){
			in_str = argv[i];
			set_it = option_list.find(in_str);

			if(set_it == option_list.end()) {
				std::cout << "ERROR: flag '" << in_str <<
					"' not valid. Please refer to the manual for usage instructions." <<
					std::endl;

				exit(EXIT_FAILURE);
			}

			// flags with parameters should eat their arguments
			// & also make sure they don't go over argc


			// Data inputs
			if(strcmp(in_str, "--bgen") == 0) {
				check_counts(in_str, i, 1, argc);
				p.bgen_file = argv[i + 1]; // bgen file
				check_file_exists(p.bgen_file);
				p.bgi_file = p.bgen_file + ".bgi";
				i += 1;
			}

			if(strcmp(in_str, "--covar") == 0) {
				check_counts(in_str, i, 1, argc);
				p.covar_file = argv[i + 1]; // covar file
				check_file_exists(p.covar_file);
				i += 1;
			}

			if(strcmp(in_str, "--environment") == 0) {
				check_counts(in_str, i, 1, argc);
				p.env_file = argv[i + 1]; // covar file
				check_file_exists(p.env_file);
				i += 1;
			}

			if(strcmp(in_str, "--random_seed") == 0) {
					p.random_seed = std::stoul(argv[i + 1]);
					i += 1;
			}

			if(strcmp(in_str, "--true_sigma") == 0) {
				check_counts(in_str, i, 1, argc);
				p.sigma = std::stod(argv[i + 1]);
				p.sim_w_noise = true;
				i += 1;
			}

			if(strcmp(in_str, "--true_hb") == 0) {
				check_counts(in_str, i, 1, argc);
				p.hb = std::stod(argv[i + 1]);
				p.rescale_coeffs = true;
				i += 1;
			}

			if(strcmp(in_str, "--true_hg") == 0) {
				check_counts(in_str, i, 1, argc);
				p.hg = std::stod(argv[i + 1]);
				p.rescale_coeffs = true;
				i += 1;
			}

			if(strcmp(in_str, "--true_hc") == 0) {
				check_counts(in_str, i, 1, argc);
				p.hc = std::stod(argv[i + 1]);
				i += 1;
			}

			if(strcmp(in_str, "--true_he") == 0) {
				check_counts(in_str, i, 1, argc);
				p.he = std::stod(argv[i + 1]);
				i += 1;
			}

			if(strcmp(in_str, "--coeffs") == 0) {
				check_counts(in_str, i, 1, argc);
				p.coeffs_file = argv[i + 1]; // covar file
				check_file_exists(p.coeffs_file);
				i += 1;
			}

			if(strcmp(in_str, "--out") == 0) {
				if (check_out == 1) {
					std::cout << "ERROR: flag '" << in_str << "' can only be provided once." << std::endl;
					exit(EXIT_FAILURE);
				}
				check_out = 1;
				check_counts(in_str, i, 1, argc);
				p.out_file = argv[i + 1];
				i += 1;
			}

			// Mode
			if(strcmp(in_str, "--mode_ssv") == 0) {
				p.mode_ssv = true;
				i += 0;
			}

			if(strcmp(in_str, "--print-keys") == 0) {
				p.mode_print_keys = true;
				i += 0;
			}

			if(strcmp(in_str, "--use_raw_covars") == 0) {
				p.use_raw_covars = true;
				i += 0;
			}

			if(strcmp(in_str, "--use_raw_env") == 0) {
				p.use_raw_env = true;
				i += 0;
			}

			if(strcmp(in_str, "--compute-env-snp-correlations") == 0) {
				p.mode_compute_correlations = true;
				i += 0;
			}

			if(strcmp(in_str, "--mode_gen_pheno") == 0) {
				p.mode_gen_pheno = true;
				i += 0;
			}

			if(strcmp(in_str, "--mode_gen2_pheno") == 0) {
				p.mode_gen2_pheno = true;
				p.sim_w_noise = true;
				i += 0;
			}

			if(strcmp(in_str, "--mode_low_mem") == 0) {
				std::cout << "WARNING: Simulating effects of low-mem mode" << std::endl;
				p.mode_low_mem = true;
				i += 0;
			}

			// Filters
			if(strcmp(in_str, "--keep_constant_variants") == 0) {
				p.keep_constant_variants = true;
				i += 0;
			}

			if(strcmp(in_str, "--maf") == 0) {
				check_counts(in_str, i, 1, argc);
				p.maf_lim = true;
				p.min_maf = std::stod(argv[i + 1]); // bgen file
				i += 1;
			}

			if(strcmp(in_str, "--info") == 0) {
				check_counts(in_str, i, 1, argc);
				p.info_lim = true;
				p.min_info = std::stod(argv[i + 1]); // bgen file
				i += 1;
			}

			if(strcmp(in_str, "--range") == 0) {
				static bool check = 0;
				if (check == 1) {
					std::cout << "ERROR: flag '" << in_str << "' can only be provided once." << std::endl;
					exit(EXIT_FAILURE);
				}
				check = 1;
				check_counts(in_str, i, 3, argc);
				p.range = true;
				p.chr = argv[i + 1];
				p.start = atoi(argv[i + 2]);
				p.end = atoi(argv[i + 3]);
				i += 3;
			}

			if(strcmp(in_str, "--incl_sample_ids") == 0) {
				check_counts(in_str, i, 1, argc);
				p.incl_sids_file = argv[i + 1]; // include sample ids file
				check_file_exists(p.incl_sids_file);
				i += 1;
			}

			if(strcmp(in_str, "--incl_rsids") == 0) {
				check_counts(in_str, i, 1, argc);
				p.incl_snps = true;
				p.incl_rsids_file = argv[i + 1]; // include variant ids file
				check_file_exists(p.incl_rsids_file);
				i += 1;
			}

			if(strcmp(in_str, "--excl_rsids") == 0) {
				check_counts(in_str, i, 1, argc);
				p.excl_snps = true;
				p.excl_rsids_file = argv[i + 1]; // include variant ids file
				check_file_exists(p.excl_rsids_file);
				i += 1;
			}

			if(strcmp(in_str, "--rsid") == 0) {
				check_counts(in_str, i, 1, argc);
				p.select_rsid = true;
				int jj = i+1;
				while(jj < argc){
					std::string arg_str(argv[jj]);
					if (arg_str.find("--") != std::string::npos) break;
					p.rsid.push_back(argv[jj]);
					jj++;
				}
				i += 1;
			}

			// Other options
			if(strcmp(in_str, "--no_geno_check") == 0) {
				p.geno_check = false;
				i += 0;
			}

			if(strcmp(in_str, "--chunk") == 0) {
				check_counts(in_str, i, 1, argc);
				p.chunk_size = std::stoi(argv[i + 1]); // bgen file
				i += 1;
			}

		}
	}

	// Sanity checks here
	if(p.mode_gen_pheno){
		bool has_bgen = p.bgen_file != "NULL";
		bool has_out = p.out_file != "NULL";
		bool has_all = (has_out && has_bgen);
		if(!has_all){
			std::cout << "ERROR: bgen and out files should be ";
			std::cout << "provided in conjunction with --mode_gen_pheno" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(p.coeffs_file == "NULL"){
			throw std::runtime_error("File of coefficients must be provided in gen pheno mode");
		}
	}
	if(p.mode_gen2_pheno){
		bool has_bgen = p.bgen_file != "NULL";
		bool has_out = p.out_file != "NULL";
		bool has_all = (has_out && has_bgen);
		if(!has_all){
			std::cout << "ERROR: bgen and out files should all be ";
			std::cout << "provided in conjunction with --mode_gen2_pheno" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(p.coeffs_file == "NULL"){
			throw std::runtime_error("File of coefficients must be provided in gen pheno mode");
		}
		if(p.rescale_coeffs){
			std::cout << "--rescale_coeffs not supported by --mode_gen2_pheno" << std::endl;
			p.rescale_coeffs = false;
		}
	}
	if(p.range || p.incl_snps || p.excl_snps){
		struct stat buf;
		p.bgi_file = p.bgen_file + ".bgi";
		if(stat(p.bgi_file.c_str(), &buf) != 0){
			std::cout << "If using --range the BGEN index file " << p.bgi_file << " must exist" << std::endl;
			throw std::runtime_error("ERROR: file does not exist");
		}
	}
}

#endif