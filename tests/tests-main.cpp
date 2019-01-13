// tests-main.cpp
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <sys/stat.h>
#include "../src/tools/eigen3.3/Dense"
#include "../src/parse_arguments.hpp"
#include "../src/data.hpp"


TEST_CASE("Data") {
	parameters p;
	p.env_file = "data/io_test/n50_p100_env.txt";

	SECTION("dXtEEX computed correctly (low mem)") {
		p.bgen_file = "data/io_test/n50_p100.bgen";
		p.bgi_file = "data/io_test/n50_p100.bgen.bgi";
		p.mode_low_mem = true;
		Data data(p);
		data.read_environment();
		data.reduce_to_complete_cases();

		CHECK(data.E(0, 0) == Approx(0.785198212));
		data.center_matrix(data.E, data.n_env);
		data.scale_matrix(data.E, data.n_env, data.env_names);
		CHECK(data.E(0, 0) == Approx(0.8957059881));

		data.read_bgen_chunk();
		EigenDataArrayXX dXtEEX_chunk(data.n_var, data.n_env * data.n_env);
		data.compute_correlations_chunk(dXtEEX_chunk);

		CHECK(data.n_env == 4);
		CHECK(data.n_samples == 50);
		CHECK(dXtEEX_chunk(0, 0) == Approx(38.9610805993));
		CHECK(dXtEEX_chunk(1, 0) == Approx(38.2995451744));
		CHECK(dXtEEX_chunk(2, 0) == Approx(33.7077899144));
		CHECK(dXtEEX_chunk(3, 0) == Approx(35.7391671158));

		CHECK(dXtEEX_chunk(0, 4) == Approx(-2.6239467101));
		CHECK(dXtEEX_chunk(1, 4) == Approx(-13.0001255314));
		CHECK(dXtEEX_chunk(2, 4) == Approx(-11.6635557299));
		CHECK(dXtEEX_chunk(3, 4) == Approx(-7.2154836264));
	}

	SECTION("dXtEEX computed correctly (low mem) w covars") {
		p.covar_file = "data/io_test/age.txt";
		p.bgen_file = "data/io_test/n50_p100.bgen";
		p.bgi_file = "data/io_test/n50_p100.bgen.bgi";
		p.mode_low_mem = true;
		Data data(p);
		data.read_environment();
		data.read_covar();
		data.reduce_to_complete_cases();

		CHECK(data.E(0, 0) == Approx(0.785198212));
		data.center_matrix(data.E, data.n_env);
		data.scale_matrix(data.E, data.n_env, data.env_names);
		CHECK(data.E(0, 0) == Approx(0.8957059881));

		data.center_matrix(data.W, data.n_covar);
		data.scale_matrix(data.W, data.n_covar, data.covar_names);
		data.regress_first_mat_from_second(data.W, "covars", data.covar_names, data.E, "env");
		CHECK(data.E(0, 0) == Approx(0.9959851422));

		data.read_bgen_chunk();
		EigenDataArrayXX dXtEEX_chunk(data.n_var, data.n_env * data.n_env);
		data.compute_correlations_chunk(dXtEEX_chunk);

		CHECK(data.n_env == 4);
		CHECK(data.n_samples == 50);
		CHECK(dXtEEX_chunk(0, 0) == Approx(42.2994405499));
		CHECK(dXtEEX_chunk(1, 0) == Approx(43.2979303929));
		CHECK(dXtEEX_chunk(2, 0) == Approx(37.6440444004));
		CHECK(dXtEEX_chunk(3, 0) == Approx(40.9258647207));

		CHECK(dXtEEX_chunk(0, 4) == Approx(-4.0453940676));
		CHECK(dXtEEX_chunk(1, 4) == Approx(-15.6140263169));
		CHECK(dXtEEX_chunk(2, 4) == Approx(-13.2508795732));
		CHECK(dXtEEX_chunk(3, 4) == Approx(-9.8081456731));
	}

	SECTION("dXtEEX computed correctly (low mem), covars, incl sids") {
		p.covar_file = "data/io_test/age.txt";
		p.bgen_file = "data/io_test/n50_p100.bgen";
		p.bgi_file = "data/io_test/n50_p100.bgen.bgi";
		p.incl_sids_file = "data/io_test/sample_ids.txt";
		p.mode_low_mem = true;
		Data data(p);
		data.read_environment();
		data.read_covar();
		data.read_incl_sids();
		data.reduce_to_complete_cases();

		CHECK(data.E(0, 0) == Approx(0.785198212));
		data.center_matrix(data.E, data.n_env);
		data.scale_matrix(data.E, data.n_env, data.env_names);
		CHECK(data.E(0, 0) == Approx(0.7617238014));

		data.center_matrix(data.W, data.n_covar);
		data.scale_matrix(data.W, data.n_covar, data.covar_names);
		data.regress_first_mat_from_second(data.W, "covars", data.covar_names, data.E, "env");
		CHECK(data.E(0, 0) == Approx(0.8123860763));

		data.read_bgen_chunk();
		EigenDataArrayXX dXtEEX_chunk(data.n_var, data.n_env * data.n_env);
		data.compute_correlations_chunk(dXtEEX_chunk);

		CHECK(data.n_env == 4);
		CHECK(data.n_samples == 28);
		CHECK(dXtEEX_chunk(0, 0) == Approx(23.2334219303));
		CHECK(dXtEEX_chunk(1, 0) == Approx(27.9920667408));
		CHECK(dXtEEX_chunk(2, 0) == Approx(24.7041225993));
		CHECK(dXtEEX_chunk(3, 0) == Approx(24.2423580715));

		CHECK(dXtEEX_chunk(0, 4) == Approx(-1.056112897));
		CHECK(dXtEEX_chunk(1, 4) == Approx(-8.526431457));
		CHECK(dXtEEX_chunk(2, 4) == Approx(-6.5950206611));
		CHECK(dXtEEX_chunk(3, 4) == Approx(-3.6842212598));
	}

	SECTION("dXtEEX computed correctly (high mem)") {
		p.bgen_file = "data/io_test/n50_p100.bgen";
		p.bgi_file = "data/io_test/n50_p100.bgen.bgi";
		Data data(p);
		data.read_environment();
		data.reduce_to_complete_cases();
		data.center_matrix(data.E, data.n_env);
		data.scale_matrix(data.E, data.n_env, data.env_names);

		data.read_bgen_chunk();
		EigenDataArrayXX dXtEEX_chunk(data.n_var, data.n_env * data.n_env);
		data.compute_correlations_chunk(dXtEEX_chunk);

		CHECK(data.n_env == 4);
		CHECK(data.n_samples == 50);
		CHECK(dXtEEX_chunk(0, 0) == Approx(38.9390135703));
		CHECK(dXtEEX_chunk(1, 0) == Approx(38.34695));
		CHECK(dXtEEX_chunk(2, 0) == Approx(33.7626));
		CHECK(dXtEEX_chunk(3, 0) == Approx(35.71962));

		CHECK(dXtEEX_chunk(0, 4) == Approx(-2.58481));
		CHECK(dXtEEX_chunk(1, 4) == Approx(-13.04073));
		CHECK(dXtEEX_chunk(2, 4) == Approx(-11.69077));
		CHECK(dXtEEX_chunk(3, 4) == Approx(-7.17068));
	}

//	SECTION("n50_p100_chr2.bgen") {
//		p.bgen_file = "data/io_test/n50_p100_chr2.bgen";
//		p.bgi_file = "data/io_test/n50_p100_chr2.bgen.bgi";
//		Data data(p);
//
//		data.read_non_genetic_data();
//		data.standardise_non_genetic_data();
//		data.read_full_bgen();
//		SECTION("Ex1. bgen read in & standardised correctly") {
//			CHECK(data.G.low_mem);
//			CHECK(data.params.low_mem);
//			CHECK(data.params.flip_high_maf_variants);
//			CHECK(data.G(0, 0) == Approx(-0.7105269065));
//			CHECK(data.G(0, 1) == Approx(-0.6480740698));
//			CHECK(data.G(0, 2) == Approx(-0.7105104917));
//			CHECK(data.G(0, 3) == Approx(-0.586791551));
//			CHECK(data.G(0, 60) == Approx(1.4862052498));
//			CHECK(data.G(0, 61) == Approx(-0.3299831646));
//			CHECK(data.G(0, 62) == Approx(-1.0968694989));
//			CHECK(data.G(0, 63) == Approx(-0.5227553607));
//			CHECK(data.G.compressed_dosage_means(60) == Approx(0.9821875));
//			CHECK(data.G.compressed_dosage_means(61) == Approx(0.10390625));
//			CHECK(data.G.compressed_dosage_means(62) == Approx(0.68328125));
//			CHECK(data.G.compressed_dosage_means(63) == Approx(0.28359375));
//			CHECK(data.n_var == 73);
//		}
//
//		SECTION("Ex1. Confirm calc_dxteex() reorders properly") {
//			data.params.dxteex_file = "data/io_test/inputs/dxteex_mixed.txt";
//			data.read_external_dxteex();
//			data.calc_dxteex();
//			CHECK(data.dXtEEX(0, 0) == Approx(54.8912155253));
//			CHECK(data.n_dxteex_computed == 73);
//		}
//	}
//
//	SECTION("n50_p100_chr2.bgen w/ 2 chunks") {
//		p.bgen_file = "data/io_test/n50_p100_chr2.bgen";
//		p.bgi_file = "data/io_test/n50_p100_chr2.bgen.bgi";
//		p.chunk_size = 72;
//		p.n_bgen_thread = 2;
//		Data data(p);
//
//		data.read_non_genetic_data();
//		data.standardise_non_genetic_data();
//		data.read_full_bgen();
//		SECTION("Ex1. bgen read in & standardised correctly") {
//			CHECK(data.G.low_mem);
//			CHECK(data.params.low_mem);
//			CHECK(data.params.flip_high_maf_variants);
//			CHECK(data.G(0, 0) == Approx(-0.7105269065));
//			CHECK(data.G(0, 1) == Approx(-0.6480740698));
//			CHECK(data.G(0, 2) == Approx(-0.7105104917));
//			CHECK(data.G(0, 3) == Approx(-0.586791551));
//			CHECK(data.G(0, 60) == Approx(1.4862052498));
//			CHECK(data.G(0, 61) == Approx(-0.3299831646));
//			CHECK(data.G(0, 62) == Approx(-1.0968694989));
//			CHECK(data.G(0, 63) == Approx(-0.5227553607));
//			CHECK(data.G.compressed_dosage_means(60) == Approx(0.9821875));
//			CHECK(data.G.compressed_dosage_means(61) == Approx(0.10390625));
//			CHECK(data.G.compressed_dosage_means(62) == Approx(0.68328125));
//			CHECK(data.G.compressed_dosage_means(63) == Approx(0.28359375));
//			CHECK(data.n_var == 73);
//		}
//	}
//
//	SECTION("Check mult_vector_by_chr"){
//		p.bgen_file = "data/io_test/n50_p100_chr2.bgen";
//		p.bgi_file = "data/io_test/n50_p100_chr2.bgen.bgi";
//		Data data(p);
//
//		data.read_non_genetic_data();
//		data.read_full_bgen();
//
//		Eigen::VectorXd vv = Eigen::VectorXd::Ones(data.G.pp);
//		Eigen::VectorXd v1 = data.G.mult_vector_by_chr(1, vv);
//		Eigen::VectorXd v2 = data.G.mult_vector_by_chr(22, vv);
//
//		CHECK(v1(0) == Approx(-9.6711528276));
//		CHECK(v1(1) == Approx(-0.4207388213));
//		CHECK(v1(2) == Approx(-3.0495872499));
//		CHECK(v1(3) == Approx(-9.1478619829));
//
//		CHECK(v2(0) == Approx(-15.6533077013));
//		CHECK(v2(1) == Approx(6.8078348334));
//		CHECK(v2(2) == Approx(-4.4887853578));
//		CHECK(v2(3) == Approx(8.9980192447));
//	}
}
