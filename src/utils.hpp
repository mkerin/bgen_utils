// class for implementation of variational bayes algorithm
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "tools/eigen3.3/Dense"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace boost_io = boost::iostreams;

inline double sigmoid(double x){
	return 1.0 / (1.0 + std::exp(-x));
}

inline std::string gen_snpkey(std::string chr, uint32_t pos,
                              std::vector< std::string > alleles){
	std::string key = chr;
	key += "~" + std::to_string(pos);
	key += "~" + alleles[0];
	key += "~" + alleles[1];
	return key;
}

inline void read_file_header(const std::string& filename,
                             std::vector<std::string>& col_names){

	// Get colnames from file
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
		std::cout << "ERROR: " << filename << " contains zero lines" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::stringstream ss;
	std::string s;
	col_names.clear();
	ss.clear();
	ss.str(line);
	while (ss >> s) {
		col_names.push_back(s);
	}
}

/***************** Typedefs *****************/
#ifdef DATA_AS_FLOAT
using scalarData          = float;
using EigenDataMatrix     = Eigen::MatrixXf;
using EigenDataVector     = Eigen::VectorXf;
using EigenDataArrayXX    = Eigen::ArrayXXf;
using EigenDataArrayX     = Eigen::ArrayXf;
using EigenRefDataMatrix  = Eigen::Ref<Eigen::MatrixXf>;
using EigenRefDataVector  = Eigen::Ref<Eigen::VectorXf>;
using EigenRefDataArrayXX = Eigen::Ref<Eigen::ArrayXXf>;
using EigenRefDataArrayX  = Eigen::Ref<Eigen::ArrayXf>;
#else
using scalarData          = double;
using EigenDataMatrix     = Eigen::MatrixXd;
using EigenDataVector     = Eigen::VectorXd;
using EigenDataArrayXX    = Eigen::ArrayXXd;
using EigenDataArrayX     = Eigen::ArrayXd;
using EigenRefDataMatrix  = Eigen::Ref<Eigen::MatrixXd>;
using EigenRefDataVector  = Eigen::Ref<Eigen::VectorXd>;
using EigenRefDataArrayXX = Eigen::Ref<Eigen::ArrayXXd>;
using EigenRefDataArrayX  = Eigen::Ref<Eigen::ArrayXd>;
#endif

#endif
