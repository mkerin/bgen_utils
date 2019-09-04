//
// Created by kerin on 2019-05-16.
//

// tests-main.cpp
#include "catch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <sys/stat.h>
#include "../src/tools/eigen3.3/Dense"
#include "../src/parse_arguments.hpp"
#include "../src/data.hpp"


TEST_CASE("parse commandline"){
	parameters p;
	char* argv[] = { (char*) "bin/bgen_utils",
		             (char*) "--sim_pheno",
		             (char*) "--incl_sample_ids", (char*) "unit/data/sample_ids.txt",
		             (char*) "--coeffs", (char*) "unit/data/sample_ids.txt",
		             (char*) "--range", (char*) "12:13-14",
		             (char*) "--out", (char*) "unit/data/tmp.txt",
		             (char*) "--bgen", (char*) "unit/data/n50_p100.bgen"};
	int argc = sizeof(argv)/sizeof(argv[0]);
	parse_arguments(p, argc, argv);

	CHECK(p.incl_sids_file == "unit/data/sample_ids.txt");
	CHECK(p.bgen_file == "unit/data/n50_p100.bgen");
	CHECK(p.covar_file == "NULL");
	CHECK(p.chr == "12");
	CHECK(p.start == 13);
	CHECK(p.end == 14);
}


TEST_CASE("sim_phenotype (SNPKEY)"){
	parameters p;
	p.bgen_file = "unit/data/n50_p100.bgen";
	p.bgi_file = "unit/data/n50_p100.bgen.bgi";
	p.coeffs_file = "unit/data/coeffs_w_snpkey.txt";

	p.mode_gen_pheno = true;
	p.random_seed = 1;
	p.min_maf = 0.01;

	Data data(p);
	data.sim_pheno();
	CHECK(data.Y(0) == Approx(-1.508216515));
	CHECK(data.Y(1) == Approx(-0.6373672634));
	CHECK(data.Y(2) == Approx(-1.0795088689));
	CHECK(data.Y(3) == Approx(1.5720520434));
	CHECK(data.Y(4) == Approx(-2.1022111011));
}

TEST_CASE("sim_phenotype (SNPIDS)"){
	parameters p;
	p.bgen_file = "unit/data/n50_p100.bgen";
	p.bgi_file = "unit/data/n50_p100.bgen.bgi";
	p.coeffs_file = "unit/data/coeffs_w_snpid.txt";

	p.mode_gen_pheno = true;
	p.random_seed = 1;
	p.min_maf = 0.01;

	Data data(p);
	data.sim_pheno();
	CHECK(data.Y(0) == Approx(-1.508216515));
	CHECK(data.Y(1) == Approx(-0.6373672634));
	CHECK(data.Y(2) == Approx(-1.0795088689));
	CHECK(data.Y(3) == Approx(1.5720520434));
	CHECK(data.Y(4) == Approx(-2.1022111011));
}
