// class for implementation of variational bayes algorithm
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include "tools/eigen3.3/Dense"

inline double sigmoid(double x){
	return 1.0 / (1.0 + std::exp(-x));
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
