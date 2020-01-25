//
// Created by kerin on 2019-10-20.
//

#ifndef BGEN_UTILS_FILE_UTILS_HPP
#define BGEN_UTILS_FILE_UTILS_HPP

#include "tools/eigen3.3/Dense"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <string>
#include <vector>
#include <map>

namespace boost_io = boost::iostreams;

namespace fileUtils {
std::string fstream_init(boost_io::filtering_ostream &my_outf,
                         const std::string &file,
                         const std::string &file_prefix = "",
                         const std::string &file_suffix = "");

std::string filepath_format(const std::string& orig,
                            const std::string& file_prefix,
                            const std::string& file_suffix);

template <typename EigenMat>
void read_matrix(const std::string &filename,
                 const long &n_rows,
                 EigenMat &M,
                 std::vector<std::string> &col_names,
                 std::map<long, bool> &incomplete_row);

void read_matrix(const std::string& filename,
                 Eigen::MatrixXd& M,
                 std::vector< std::string >& col_names);

void read_matrix(const std::string& filename,
                 Eigen::MatrixXd& M);

void write_matrix(const Eigen::MatrixXd& mat,
                  const std::string& filename,
                  const std::vector<std::string>& header);

void write_matrix(const Eigen::ArrayXd& mat,
				  const std::string& filename,
				  const std::vector<std::string>& header);
}

#endif //BGEN_UTILS_FILE_UTILS_HPP
