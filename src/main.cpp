// Adapted from example/bgen_to_vcf.cpp by Gavin Band
// From project at http://bitbucket.org/gavinband/bgen/get/master.tar.gz

#include "parse_arguments.hpp"
#include "parameters.hpp"
#include "data.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"
#include "version.h"

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

int main( int argc, char** argv ) {
	parameters p;

	std::cout << "=================="<< std::endl;
	std::cout << "LEMMA utils v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << std::endl;
	std::cout << "==================" << std::endl << std::endl;

	try {
		parse_arguments(p, argc, argv);
		Data data(p);

		// filter - incl sample ids
		if(p.incl_sids_file != "NULL"){
			data.read_incl_sids();
		}
		genfile::bgen::IndexQuery::UniquePtr query = genfile::bgen::IndexQuery::create(p.bgi_file);

		// filter - range
		if (p.range){
			std::cout << "Selecting range..." << std::endl;
			if (p.chr.length() == 1) p.chr = "0" + p.chr;
			std::cout << p.chr << ":" << p.start << "-" << p.end << std::endl;
			genfile::bgen::IndexQuery::GenomicRange rr1(p.chr, p.start, p.end);
			query->include_range( rr1 );
		}

		// filter - incl rsids
		if(p.incl_snps){
			data.read_incl_rsids();
			std::cout << "Filtering SNPs by rsid, using bgi file: " << p.bgi_file << std::endl;
			query->include_rsids( data.incl_rsid_list );
		}

		if(p.excl_snps){
			data.read_excl_rsids();
			std::cout << "Filtering SNPs by rsid, using bgi file: " << p.bgi_file << std::endl;
			query->exclude_rsids( data.excl_rsid_list );
		}

		// filter - select single rsid
		if(p.select_rsid){
			std::sort(p.rsid.begin(), p.rsid.end());
			std::cout << "Filtering to rsids:" << std::endl;
			for (int kk = 0; kk < p.rsid.size(); kk++) std::cout << p.rsid[kk]<< std::endl;
			query->include_rsids( p.rsid );
		}

		// Summary info
		query->initialise();
		data.bgenView->set_query( query );
		data.bgenView->summarise(std::cout);

		data.output_init();
		if(p.mode_ssv){
			std::cout << "Computing sum of column variances" << std::endl;
			data.calc_ssv();
		} else if(p.mode_pred_pheno){
			data.pred_pheno();
		} else if(p.mode_gen_pheno) {
			data.gen_pheno();
		} else if(p.mode_gen2_pheno) {
			data.gen2_pheno();
		} else if(p.mode_compute_correlations){
			data.compute_correlations();
		} else if(p.mode_print_keys){
			data.print_keys();
		}
		data.output_results();

		return 0 ;
	}
	catch( genfile::bgen::BGenError const& e ) {
		std::cerr << "!! Uh-oh, error parsing bgen file.\n" ;
		return -1 ;
	}
}


