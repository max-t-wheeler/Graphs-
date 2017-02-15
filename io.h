#include <iostream>
#include <fstream>
#include <vector>

bool logging = 0;

//set logging to 1
bool log() {
	logging = 1;
}

//input a value
template <class S>
void input(S entry) {

	std::cin >> entry;

}

//print single value
template <class S>
void output(const S entry) {

	std::cout << entry << '\n';

}

//print large single value w/o sci notation
template <class S>
void output_fixed(const S entry) {

	std::cout << std::fixed << entry << '\n';

}

//print array components
template <class S>
void output_array(const S entry, const int n, const char sep) {

	for (int i = 0; i < n; ++i) {
		std::cout << entry[i] << sep;
	}

	std::cout << '\n';

}

//print matrix by components
template <class S>
void output_matrix(const std::vector<S> M, const int m, const int n) {

	for (int i = 0; i < m; ++i) {

		for (int j = 0; j < n; ++j) {
			std::cout << M[i*n+j] << ' ';
		}

		std::cout << '\n';

	}

}

//print matrix by vector components
template <class S>
void output_matrix(const std::vector<std::vector<S> > M, const int m, const int n) {

	for (int i = 0; i < m; ++i) {

		for (int j = 0; j < n; ++j) {
			std::cout << M[i][j] << ' ';
		}

		std::cout << '\n';

	}

}

//print vector components
template <class S>
void output_vector(const std::vector<S> v) {

	for (int i = 0; i < v.size(); ++i) {
		std::cout << v[i] << ' ';
	}
	
	std::cout << '\n';

}