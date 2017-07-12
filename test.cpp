#include <iostream>
#include <typeinfo>
#include <string>

int main(int argc, char* argv[]) {

	if (argc < 2) {

		std::cout << "Error. Please provide more input." << '\n';

		return 0;
	}

	std::cout << argc << '\n';

	for (int i = 1; i < argc; ++i) {
		std::cout << (int)*argv[i] - 48 << '\n';
	}

	std::string h = "hello";

	std::cout << h << '\n';

	return 0;

}