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
#include "../src/data.hpp"


TEST_CASE("sim_phenotype"){
	parameters p;
	p.bgen_file = "data/io_test/n50_p100.bgen";
	p.bgi_file = "data/io_test/n50_p100.bgen.bgi";
	p.coeffs_file = "data/io_test/coeffs_w_snpkey.txt";

	p.mode_gen_pheno = true;
	p.random_seed = 1;
	p.min_maf = 0.01;

	Data data(p);
	data.gen_pheno();
	CHECK(data.Y(0) == Approx(-1.4454847077));
	CHECK(data.Y(1) == Approx(-0.5664446173));
	CHECK(data.Y(2) == Approx(-1.0768929511));
	CHECK(data.Y(3) == Approx(1.6308280995));
	CHECK(data.Y(4) == Approx(-2.0189108167));
}

