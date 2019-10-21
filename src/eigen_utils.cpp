//
// Created by kerin on 2019-10-20.
//

#include "utils.hpp"
#include "eigen_utils.hpp"
#include "tools/eigen3.3/Dense"

#include <vector>
#include <string>

EigenDataMatrix EigenUtils::subset_rows(const EigenDataMatrix &orig, const std::vector<long> &valid_points) {
	long n_cols = orig.cols(), n_rows = valid_points.size();
	EigenDataMatrix subset(n_rows, n_cols);

	for(int kk = 0; kk < n_rows; kk++) {
		subset.row(kk) = orig.row(valid_points[kk]);
	}
	return subset;
}

void EigenUtils::scale_matrix(EigenDataMatrix &M, long &n_cols, std::vector<std::string> &col_names) {
	// Scale eigen matrix passed by reference.
	// Removes columns with zero variance + updates col_names.
	// Only call on matrixes which have been reduced to complete cases,
	// as no check for incomplete rows.
	long n_samples = M.rows();

	std::vector<size_t> keep;
	std::vector<std::string> keep_names;
	std::vector<std::string> reject_names;
	for (int k = 0; k < n_cols; k++) {
		double sigma = 0.0;
		double count = 0;
		for (int i = 0; i < n_samples; i++) {
			double val = M(i, k);
			sigma += val * val;
			count += 1;
		}

		sigma = sqrt(sigma/(count - 1));
		if (sigma > 1e-12) {
			for (int i = 0; i < n_samples; i++) {
				M(i, k) /= sigma;
			}
			keep.push_back(k);
			keep_names.push_back(col_names[k]);
		} else {
			reject_names.push_back(col_names[k]);
		}
	}

	if (keep.size() != n_cols) {
		for (std::size_t i = 0; i < keep.size(); i++) {
			M.col(i) = M.col(keep[i]);
		}
		M.conservativeResize(M.rows(), keep.size());

		n_cols = keep.size();
		col_names = keep_names;
	}

	if (n_cols == 0) {
		std::cout << "WARNING: No columns left with nonzero variance after scale_matrix()" << std::endl;
	}
}

void EigenUtils::center_matrix(EigenDataMatrix &M, long &n_cols) {
	long n_samples = M.rows();

	for (int k = 0; k < n_cols; k++) {
		double mu = 0.0;
		double count = 0;
		for (int i = 0; i < n_samples; i++) {
			mu += M(i, k);
			count += 1;
		}

		mu = mu / count;
		for (int i = 0; i < n_samples; i++) {
			M(i, k) -= mu;
		}
		// std::cout << "Mean centered matrix:" << std::endl << M << std::endl;
	}
}