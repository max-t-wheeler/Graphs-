#include <cstdlib>
#include <cmath>
#include <armadillo>

#include "io.h"
//#include "abstract/numbers_operations_package.h"

//#define pi  3.141592653589793238462643383279502884L

//return index of value in regspace space
template <class S>
int get_ID(S input, const double val, const double precision) {

	if (input.is_empty()) {

		if (logging) {
			output("Input entered in get_ID has not been initialized");
		}

		return 0;

	}

	int n;
	unsigned int count = 0;
	
	for (int i = 0; i < input.size(); ++i) {

		if (std::abs(val - input[i]) < precision) {

			count++;
			n = i;

		}

	}

	if(count == 0) {

		if (logging) {
			output("Value entered int get_ID not in input set");
		}

		return -1;

	}
	else {
		return n;
	}

}

//vector contains zero
template <class S>
bool contains_zero(const S input, const double tol) {

	if (input.is_empty()) {

		if (logging) {
			output("Input entered in contains has not been initialized");
		}

		return 0;
		
	}

	unsigned int count = 0;

	for (int i = 0; i <  input.size(); ++i) {

		if (std::abs(input[i]) < tol) {
			count++;
		}

	}

	if(count > 0) {
		return 1;
	}
	else {
		return 0;
	}

}

//determine whether a number is prime
template <class S>
bool is_prime(const S n) {

	if (n < 0) {

		if (logging) {
			output("Value entered in is_prime must be a positive integer");
		}

		return 0;

	}

	long count = 0;

	for (int i = 2; i < n; ++i) {

		if (n%i == 0) {
			count++;
		}

	}
	
	if (count > 0) {
		return false;
	}
	else {
		return true;
	}

}

//find the factorial of an integer
long long factorial(const int n) {

	if (n < 0) {

		if (logging) {
			output("Value entered in factorial must be a positive integer");
		}

		return 0;

	}

    long long ans = 1;

    for (int i = 1; i <= n; ++i) {
        ans *= i;
    }

    return ans;

}

//find the binomial coefficient of two integes (http://www.geeksforgeeks.org/dynamic-programming-set-9-binomial-coefficient/)
template <class S>
long long binomial_coefficient(const S n, const S k) {

	if (k > n) {

		if (logging) {
			output("Criterion not met in binomial_coefficient");
		}

		return 0;

	}

    S C[k+1];

    memset(C, 0, sizeof(C));
 
    C[0] = 1;
 
    for (int i = 1; i <= n; i++) {
        for (int j = std::min(i, k); j > 0; j--) {
            C[j] = C[j] + C[j-1];
        }
    }

    return C[k];
}

//return a random value between to bounds
double Rand(const double min, const double max) {

    double x = (double)rand() / RAND_MAX;
    
    return min + x * (max - min);

}