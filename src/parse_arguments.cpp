//
// Created by kerin on 2019-05-16.
//

#include "parameters.hpp"
#include "cxxopts.hpp"

#include <cassert>
#include <iostream>
#include <set>
#include <cstring>
#include <sys/stat.h>
#include <regex>
#include <stdexcept>

void check_file_exists(const std::string &filename) {
	// Throw error if given file does not exist.
	// NB: Doesn't check if file is empty etc.
	struct stat buf;
	if(stat(filename.c_str(), &buf) != 0) {
		std::cout << "File " << filename << " does not exist" << std::endl;
		throw std::runtime_error("ERROR: file does not exist");
	}
}

void parse_arguments(parameters &p, int argc, char **argv) {

	cxxopts::Options options("bgen utils", "Utility programme for interacting with bgen files. "
	                         "Includes functionality such as predicting phenotypes "
	                         "given a file of coefficients and simulating phenotypes (ie with added noise).");

	options.add_options("General")
	    ("sim_pheno", "Simulate phenotype of form Y = X beta + epsilon", cxxopts::value<bool>(p.mode_gen_pheno))
	    ("pred_pheno", "Predict phenotype Y = X beta", cxxopts::value<bool>(p.mode_pred_pheno))
	    ("bgen", "BGEN file", cxxopts::value<std::string>(p.bgen_file))
	    ("environment", "Path to file of environments", cxxopts::value<std::string>(p.env_file))
	    ("environment_weights", "Path to file of interaction weights (default: 1 if using a single environment)", cxxopts::value<std::string>(p.env_profile_file))
	    ("coeffs", "Path to file of SNP effect sizes. Header required.", cxxopts::value<std::string>(p.coeffs_file))
	    ("out", "Filepath to output", cxxopts::value<std::string>(p.out_file))
	    ("incl_sample_ids", "Text file of sample ids to include (no header, 1 per line)",
	    cxxopts::value<std::string>(p.incl_sids_file))
	    ("incl_rsids", "Text file of rsids to include (no header, 1 per line)",
	    cxxopts::value<std::string>(p.incl_rsids_file))
	    ("excl_rsids", "Text file of rsids to exclude (no header, 1 per line)",
	    cxxopts::value<std::string>(p.excl_rsids_file))
	    ("range", "Genomic range in format chr:start-end",
	    cxxopts::value<std::string>())
	    ("random_seed", "Seed used when simulating a phenotype (default: random)",
	    cxxopts::value<unsigned int>(p.random_seed))
	    ("true_sigma", "Variance of gaussian noise added to simulated phenotype (default: 1). Must be greater than zero.",
	    cxxopts::value<double>(p.sigma))
	    ("use_raw_dosage", "Use expected dosage instead of normalising columns of the dosage matrix to mean zero and variance one.")
	    ("h, help", "")
	;

	options.add_options("Internal")
	    ("use_raw_env", "")
	    ("coeffs2", "", cxxopts::value<std::string>(p.coeffs2_file))
	    ("environment_weights2", "Path to file of interaction weights (default: 1 if using a single environment)", cxxopts::value<std::string>(p.env_profile_file))
	    ("print_keys", "", cxxopts::value<bool>(p.mode_print_keys))
	    ("maf", "", cxxopts::value<double>())
	    ("compute-env-snp-correlations", "", cxxopts::value<bool>(p.mode_compute_correlations))
	    ("covar", "File of covariables", cxxopts::value<std::string>(p.covar_file))
	;

	try{
		auto opts = options.parse(argc, argv);
		auto args = opts.arguments();

		if (opts.count("help") || args.empty()) {
			std::cout << options.help({"General"}) << std::endl;
			std::exit(0);
		}

		if(p.bgen_file != "NULL") {
			p.bgi_file = p.bgen_file + ".bgi";
			check_file_exists(p.bgen_file);
			check_file_exists(p.bgi_file);
		}
		if(p.env_file != "NULL") {
			check_file_exists(p.env_file);
		}
		if(p.env_profile_file != "NULL") {
			check_file_exists(p.env_profile_file);
		}

		if(p.coeffs_file != "NULL") {
			check_file_exists(p.coeffs_file);
		}

		if(p.incl_sids_file != "NULL") {
			check_file_exists(p.incl_sids_file);
		}

		if(p.incl_rsids_file != "NULL") {
			check_file_exists(p.incl_rsids_file);
		}

		if(p.excl_rsids_file != "NULL") {
			check_file_exists(p.excl_rsids_file);
		}

		if(opts.count("maf")) {
			p.maf_lim = true;
			p.min_maf = opts["maf"].as<double>();
		}

		if(opts.count("info")) {
			p.info_lim = true;
			p.min_info = opts["info"].as<double>();
		}

		if(opts.count("range")) {
			auto ss = opts["range"].as<std::string>();
			p.range = true;
			p.chr = ss.substr(0, ss.find(':'));
			p.start = std::atoi(ss.substr(ss.find(':')+1, ss.find('-')).c_str());
			p.end = std::atoi(ss.substr(ss.find('-')+1, ss.size()).c_str());
		}

		if(opts.count("true_sigma")) {
			assert(p.sigma > 0);
		}

		if(opts.count("use_raw_dosage")) {
			p.normalise_genotypes = false;
		}

		if(opts.count("use_raw_env")) {
			p.use_raw_env = true;
		}

	} catch (const cxxopts::OptionException& e) {
		std::cout << "error parsing options: " << e.what() << std::endl;
		std::exit(1);
	}

	// Sanity checks here
	if(p.mode_gen_pheno) {
		if(p.out_file == "NULL") {
			throw std::runtime_error("Output file must be provided in gen pheno mode");
		}
		if(p.bgen_file == "NULL") {
			throw std::runtime_error("bgen file of coefficients must be provided in gen pheno mode");
		}
		if(p.coeffs_file == "NULL") {
			throw std::runtime_error("File of coefficients must be provided in gen pheno mode");
		}
	}

	if(p.mode_pred_pheno) {
		if(p.out_file == "NULL") {
			throw std::runtime_error("Output file must be provided in pred pheno mode");
		}
		if(p.bgen_file == "NULL") {
			throw std::runtime_error("bgen file of coefficients must be provided in pred pheno mode");
		}
		if(p.coeffs_file == "NULL") {
			throw std::runtime_error("File of coefficients must be provided in pred pheno mode");
		}
	}

	/*** Obsolete code after converting to CXXOPTS ***/
	/*** Commented out functionality deemed unnecessary for simplified branch ***/
//	char *in_str;
//	int i;
//	std::set<std::string> option_list {
//		"--covar",
//		"--environment",
//		"--env_profile",
//		"--chunk",
//		"--rsid",
//		"--no_geno_check",
//		"--mode_ssv",
//		"--mode_gen2_pheno",
//		"--coeffs2",
//		"--true_hb",
//		"--true_hg",
//		"--true_hb2",
//		"--true_hg2",
//		"--true_hc",
//		"--true_he",
//		"--keep_constant_variants",
//		"--print-keys",
//		"--mode_low_mem",
//		"--use_raw_covars",
//		"--use_raw_env",
//		"--compute-env-snp-correlations"
//	};
//
//	std::set<std::string>::iterator set_it;
//	// for --version (& no other args), print splash screen and exit.
//	if (argc == 2 && strcmp(argv[1], "--version") == 0) {
//		printf("%s\n\n", splash);
//		std::exit(EXIT_SUCCESS);
//	}
//
//	// read in and check option flags
//	for (i = 0; i < argc; i++) {
//		std::string in_str1(argv[i]);
//		in_str = argv[i];
//		if (in_str1 == "--version" || in_str1 == "--help") {
//			std::cout << "ERROR: flag '" << in_str << "' cannot be used with any other flags." << std::endl;
//			std::exit(EXIT_FAILURE);
//		}
//	}
//
//	// Ensure some arguments only appear once
//	bool check_out = 0;
//
//	for(i = 0; i < argc; i++){
//		if(*argv[i] == '-'){
//			in_str = argv[i];
//			std::string in_str1(argv[i]);
//			set_it = option_list.find(in_str);
//
//			if(set_it == option_list.end()) {
//				std::cout << "ERROR: flag '" << in_str <<
//					"' not valid. Please refer to the manual for usage instructions." <<
//					std::endl;
//
//				exit(EXIT_FAILURE);
//			}
//
//			// flags with parameters should eat their arguments
//			// & also make sure they don't go over argc
//
//
//			// Data inputs
//
//			if(in_str1 == "--covar") {
//				check_counts(in_str, i, 1, argc);
//				p.covar_file = argv[i + 1]; // covar file
//				check_file_exists(p.covar_file);
//				i += 1;
//			}
//
//			if(in_str1 == "--environment") {
//				check_counts(in_str, i, 1, argc);
//				p.env_file = argv[i + 1]; // covar file
//				check_file_exists(p.env_file);
//				i += 1;
//			}
//
//
//			if(in_str1 == "--true_hb") {
//				check_counts(in_str, i, 1, argc);
//				p.hb = std::stod(argv[i + 1]);
//				p.rescale_coeffs = true;
//				i += 1;
//			}
//
//			if(in_str1 == "--true_hg") {
//				check_counts(in_str, i, 1, argc);
//				p.hg = std::stod(argv[i + 1]);
//				p.rescale_coeffs = true;
//				i += 1;
//			}
//
//			if(in_str1 == "--true_hb2") {
//				check_counts(in_str, i, 1, argc);
//				p.hb2 = std::stod(argv[i + 1]);
//				p.rescale_coeffs = true;
//				i += 1;
//			}
//
//			if(in_str1 == "--true_hg2") {
//				check_counts(in_str, i, 1, argc);
//				p.hg2 = std::stod(argv[i + 1]);
//				p.rescale_coeffs = true;
//				i += 1;
//			}
//
//			if(in_str1 == "--true_hc") {
//				check_counts(in_str, i, 1, argc);
//				p.hc = std::stod(argv[i + 1]);
//				i += 1;
//			}
//
//			if(in_str1 == "--true_he") {
//				check_counts(in_str, i, 1, argc);
//				p.he = std::stod(argv[i + 1]);
//				i += 1;
//			}
//
//			if(in_str1 == "--coeffs2") {
//				check_counts(in_str, i, 1, argc);
//				p.coeffs2_file = argv[i + 1];
//				check_file_exists(p.coeffs2_file);
//				i += 1;
//			}
//
//			// Mode
//			if(in_str1 == "--mode_ssv") {
//				p.mode_ssv = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--print-keys") {
//				p.mode_print_keys = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--use_raw_covars") {
//				p.use_raw_covars = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--use_raw_env") {
//				p.use_raw_env = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--compute-env-snp-correlations") {
//				p.mode_compute_correlations = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--mode_gen2_pheno") {
//				p.mode_gen2_pheno = true;
//				p.sim_w_noise = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--mode_low_mem") {
//				std::cout << "WARNING: Simulating effects of low-mem mode" << std::endl;
//				p.mode_low_mem = true;
//				i += 0;
//			}
//
//			// Filters
//			if(in_str1 == "--keep_constant_variants") {
//				p.keep_constant_variants = true;
//				i += 0;
//			}
//
//			if(in_str1 == "--rsid") {
//				check_counts(in_str, i, 1, argc);
//				p.select_rsid = true;
//				int jj = i+1;
//				while(jj < argc){
//					std::string arg_str(argv[jj]);
//					if (arg_str.find("--") != std::string::npos) break;
//					p.rsid.push_back(argv[jj]);
//					jj++;
//				}
//				i += 1;
//			}
//
//			// Other options
//			if(in_str1 == "--no_geno_check") {
//				p.geno_check = false;
//				i += 0;
//			}
//
//			if(in_str1 == "--chunk") {
//				check_counts(in_str, i, 1, argc);
//				p.chunk_size = std::stoi(argv[i + 1]); // bgen file
//				i += 1;
//			}
//
//			if(in_str1 == "--print_causal_rsids") {
//				p.print_causal_rsids = true;
//				i += 0;
//			}
//
//		}
//	}
//
//	if(p.mode_gen2_pheno){
//		bool has_bgen = p.bgen_file != "NULL";
//		bool has_out = p.out_file != "NULL";
//		bool has_all = (has_out && has_bgen);
//		if(!has_all){
//			std::cout << "ERROR: bgen and out files should all be ";
//			std::cout << "provided in conjunction with --mode_gen2_pheno" << std::endl;
//			std::exit(EXIT_FAILURE);
//		}
//		if(p.coeffs_file == "NULL"){
//			throw std::runtime_error("File of coefficients must be provided in gen pheno mode");
//		}
//		if(p.rescale_coeffs){
//			std::cout << "--rescale_coeffs not supported by --mode_gen2_pheno" << std::endl;
//			p.rescale_coeffs = false;
//		}
//	}
}
