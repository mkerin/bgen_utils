// Adapted from example/bgen_to_vcf.cpp by Gavin Band
// From project at http://bitbucket.org/gavinband/bgen/get/master.tar.gz

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <memory>
#include "parse_arguments.hpp"
#include "class.h"
#include "data.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"
#include "version.h"

void read_directory(const std::string& name, std::vector<std::string>& v);

// TODO: Sensible restructuring of interaction code
// TODO: Use high precision double for pval
// TODO: implement tests for info filter
// TODO: tests for read_pheno, read_covar? Clarify if args for these are compulsory.
// TODO: copy argument_sanity()

// Efficiency changes:
// 1) Use --range to edit query before reading bgen file
// 2) Is there an option to skip sample ids?
//    If so we could fill n_samples at start
//    - read_covar() & read_pheno()
//    - then edit query with incomplete cases before reading bgen

// Notes:
// - Replace missing data in BGEN files with mean

// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
int main( int argc, char** argv ) {
	parameters p;

	try {
		parse_arguments(p, argc, argv);
		data Data( p.bgen_file );
		Data.params = p;

		// filter - incl sample ids
		if(p.incl_sids_file != "NULL"){
			Data.read_incl_sids();
		}
		genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);

		// filter - range
		if (p.range){
			std::cout << "Selecting range..." << std::endl;
			if (p.chr.length() == 1) p.chr = "0" + p.chr;
			std::cout << p.chr << ":" << p.start << "-" << p.end << std::endl;
			// genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);
			genfile::bgen::IndexQuery::GenomicRange rr1(p.chr, p.start, p.end);
			query->include_range( rr1 );
			// query->include_range( rr1 ).initialise();
			// Data.bgenView->set_query( query );
		}

		// filter - incl rsids
		if(p.incl_snps){
			Data.read_incl_rsids();
			std::cout << "Filtering SNPs by rsid, using bgi file: " << p.bgi_file << std::endl;
			// genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);
			query->include_rsids( Data.incl_rsid_list );
			// query->include_rsids( Data.rsid_list ).initialise();
			// Data.bgenView->set_query( query );
		}

		if(p.excl_snps){
			Data.read_excl_rsids();
			std::cout << "Filtering SNPs by rsid, using bgi file: " << p.bgi_file << std::endl;
			// genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);
			query->exclude_rsids( Data.excl_rsid_list );
			// query->include_rsids( Data.rsid_list ).initialise();
			// Data.bgenView->set_query( query );
		}

		// filter - select single rsid
		if(p.select_rsid){
			std::sort(p.rsid.begin(), p.rsid.end());
			std::cout << "Filtering to rsids:" << std::endl;
			for (int kk = 0; kk < p.rsid.size(); kk++) std::cout << p.rsid[kk]<< std::endl;
			// genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);
			query->include_rsids( p.rsid );
			// query->include_rsids( p.rsid ).initialise();
			// Data.bgenView->set_query( query );
		}

		// Summary info
		query->initialise();
		Data.bgenView->set_query( query );
		Data.bgenView->summarise(std::cout);

		std::cout << "Computing sum of column variances" << std::endl;
		Data.output_init();
		if(p.mode_ssv){
			Data.calc_ssv();
		} else if(p.mode_gen_pheno) {
			Data.gen_pheno();
		} else if(p.mode_gen2_pheno) {
			Data.gen2_pheno();
		} else if(p.mode_compute_correlations){
			Data.compute_correlations();
		} else if(p.mode_print_keys){
			Data.print_keys();
		}
		Data.output_results();

		return 0 ;
	}
	catch( genfile::bgen::BGenError const& e ) {
		std::cerr << "!! Uh-oh, error parsing bgen file.\n" ;
		return -1 ;
	}
}


void read_directory(const std::string& name, std::vector<std::string>& v) {
	DIR* dirp = opendir(name.c_str());
	struct dirent * dp;
	while ((dp = readdir(dirp)) != NULL) {
		v.push_back(dp->d_name);
	}
	closedir(dirp);
}
