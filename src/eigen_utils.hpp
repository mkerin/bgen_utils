//
// Created by kerin on 2019-10-20.
//

#ifndef BGEN_UTILS_EIGEN_UTILS_HPP
#define BGEN_UTILS_EIGEN_UTILS_HPP

#include "utils.hpp"
#include "tools/eigen3.3/Dense"

#include <vector>
#include <string>

namespace EigenUtils {

	void center_matrix(EigenDataMatrix & M,
						long& n_cols);

	void scale_matrix(EigenDataMatrix & M,
					   long& n_cols,
					   std::vector< std::string >& col_names);


	EigenDataMatrix subset_rows(const EigenDataMatrix &orig, const std::vector<long> &valid_points);
}

#endif //BGEN_UTILS_EIGEN_UTILS_HPP
