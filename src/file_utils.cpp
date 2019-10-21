//
// Created by kerin on 2019-10-20.
//

#include "file_utils.hpp"
#include "tools/eigen3.3/Dense"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <map>

namespace boost_io = boost::iostreams;

std::string fileUtils::filepath_format(const std::string& orig,
                                       const std::string& file_prefix,
                                       const std::string& file_suffix){
	std::string filepath   = orig;
	std::string dir        = filepath.substr(0, filepath.rfind('/')+1);
	std::string stem_w_dir = filepath.substr(0, filepath.find('.'));
	std::string stem       = stem_w_dir.substr(stem_w_dir.rfind('/')+1, stem_w_dir.size());
	std::string ext        = filepath.substr(filepath.find('.'), filepath.size());

	std::string ofile      = dir + file_prefix + stem + file_suffix + ext;
	return ofile;
}

void fileUtils::write_matrix(Eigen::Ref<Eigen::MatrixXd> mat,
                             const std::string& filename,
                             const std::vector<std::string>& header){
	assert(header.size() == mat.cols());
	long n_cols = mat.cols();
	long n_rows = mat.rows();

	boost_io::filtering_ostream outf;
	std::string gz_str = ".gz";
	if (filename.find(gz_str) != std::string::npos) {
		outf.push(boost_io::gzip_compressor());
	}
	outf.push(boost_io::file_sink(filename));

	for (long cc = 0; cc < n_cols; cc++) {
		outf << header[cc];
		outf << (cc != n_cols - 1 ? " " : "");
	}
	outf << std::endl;
	for (long ii = 0; ii < n_rows; ii++) {
		for (long cc = 0; cc < n_cols; cc++) {
			outf << mat(ii, cc);
			outf << (cc != n_cols - 1 ? " " : "");
		}
		outf << std::endl;
	}

	outf.pop();
	boost_io::close(outf);
}

std::string fileUtils::fstream_init(boost_io::filtering_ostream &my_outf,
                                    const std::string &file,
                                    const std::string &file_prefix,
                                    const std::string &file_suffix) {

	std::string filepath   = file;
	std::string dir        = filepath.substr(0, filepath.rfind('/')+1);
	std::string stem_w_dir = filepath.substr(0, filepath.find('.'));
	std::string stem       = stem_w_dir.substr(stem_w_dir.rfind('/')+1, stem_w_dir.size());
	std::string ext        = filepath.substr(filepath.find('.'), filepath.size());

	std::string ofile      = dir + file_prefix + stem + file_suffix + ext;

	// Allows prefix to contain subfolders
	boost::filesystem::path bfilepath(ofile);
	if(!boost::filesystem::exists(bfilepath.parent_path())) {
		boost::filesystem::create_directories(bfilepath.parent_path());
	}

	my_outf.reset();
	std::string gz_str = ".gz";
	if (file.find(gz_str) != std::string::npos) {
		my_outf.push(boost_io::gzip_compressor());
	}
	my_outf.push(boost_io::file_sink(ofile));
	return ofile;
}


template <typename EigenMat>
void fileUtils::read_matrix(const std::string &filename,
                            const long &n_rows,
                            EigenMat &M,
                            std::vector<std::string> &col_names,
                            std::map<long, bool> &incomplete_row) {
	/* Assumptions:
	   - n_rows constant (number of samples constant across files)dd
	 */

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
		std::cout << "ERROR: " << filename << " contains zero lines." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::stringstream ss;
	std::string s;
	int n_cols = 0;
	ss.clear();
	ss.str(line);
	while (ss >> s) {
		++n_cols;
		col_names.push_back(s);
	}
	std::cout << "Reading matrix of size " << n_rows << " x " << n_cols << " from " << filename << std::endl;

	// Write remainder of file to Eigen matrix M
	incomplete_row.clear();
	M.resize(n_rows, n_cols);
	int i = 0;
	double tmp_d;
	while (getline(fg, line)) {
		if (i >= n_rows) {
			throw std::runtime_error("ERROR: could not convert txt file (too many lines).");
		}
		ss.clear();
		ss.str(line);
		for (int k = 0; k < n_cols; k++) {
			std::string sss;
			ss >> sss;
			/// NA
			if (sss == "NA" || sss == "NAN" || sss == "NaN" || sss == "nan") {
				tmp_d = 0;
				incomplete_row[i] = true;
			} else {
				try{
					tmp_d = stod(sss);
				} catch (const std::invalid_argument &exc) {
					std::cout << sss << " on line " << i << std::endl;
					throw;
				}
			}
			M(i, k) = tmp_d;
		}
		i++;
	}
	if (i < n_rows) {
		throw std::runtime_error("ERROR: could not convert txt file (too few lines).");
	}
}

void fileUtils::read_matrix(const std::string &filename,
                            Eigen::MatrixXd &M,
                            std::vector <std::string> &col_names) {
	/* Assumptions:
	   - dimensions unknown
	   - assume no missing values
	 */

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

	// Read file twice to acertain number of lines
	std::string line;
	int n_rows = 0;
	getline(fg, line);
	while (getline(fg, line)) {
		n_rows++;
	}
	fg.reset();
	if (filename.find(gz_str) != std::string::npos) {
		fg.push(boost_io::gzip_decompressor());
	}
	fg.push(boost_io::file_source(filename));

	// Reading column names
	if (!getline(fg, line)) {
		std::cout << "ERROR: " << filename << " contains zero lines." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::stringstream ss;
	std::string s1;
	int n_cols = 0;
	ss.clear();
	ss.str(line);
	while (ss >> s1) {
		++n_cols;
		col_names.push_back(s1);
	}
	std::cout << "Reading matrix of size " << n_rows << " x " << n_cols << " from " << filename << std::endl;

	// Write remainder of file to Eigen matrix M
	M.resize(n_rows, n_cols);
	int i = 0;
	double tmp_d;
	while (getline(fg, line)) {
		if (i >= n_rows) {
			throw std::runtime_error("ERROR: could not convert txt file (too many lines).");
		}
		ss.clear();
		ss.str(line);
		for (int k = 0; k < n_cols; k++) {
			std::string s;
			ss >> s;
			try{
				tmp_d = stod(s);
			} catch (const std::invalid_argument &exc) {
				std::cout << s << " on line " << i << std::endl;
				throw;
			}

			M(i, k) = tmp_d;
		}
		i++;
	}
	if (i < n_rows) {
		throw std::runtime_error("ERROR: could not convert txt file (too few lines).");
	}
}

void fileUtils::read_matrix(const std::string& filename,
                            Eigen::MatrixXd& M){
	std::vector<std::string> placeholder;
	fileUtils::read_matrix(filename, M, placeholder);
}

// Explicit instantiation
// https://stackoverflow.com/questions/2152002/how-do-i-force-a-particular-instance-of-a-c-template-to-instantiate

template void fileUtils::read_matrix(const std::string&, const long&, Eigen::MatrixXf&,
                                     std::vector<std::string>&, std::map<long, bool>&);
template void fileUtils::read_matrix(const std::string&, const long&, Eigen::MatrixXd&,
                                     std::vector<std::string>&, std::map<long, bool>&);
