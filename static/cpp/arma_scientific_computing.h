#include "math_package.h"

//machine epsilon (Holmes Table 1.2)
double eps(const std::string precision) {

	if (precision == "single") {
		return 1/(pow(2, 23));
	}
	else if (precision == "double") {
		return 1/(pow(2, 52));
	}
	else if (precision == "quadruple") {
		return 1/(pow(2, 112));
	}
	else {

		if (logging) {
			output("String entered in eps must be single, double or quadruple");
		}

		return 0;

	}

}

//return maximum exponent of specified precision (Holmes Chapter 1.2.1)
double E_M(const std::string precision) {

	if (precision == "single") {
		return pow(2, 7) - 1;
	}
	else if (precision == "double") {
		return pow(2, 10) - 1;
	}
	else if (precision == "quadruple") {
		return pow(2, 14) - 1;
	}
	else {

		if (logging) {
			output("String entered in E_M must be single, double or quadruple");
		}

		return 0;
		
	}

}

//return maximum exponent of specified precision (Holmes Chapter 1.2.1)
double E_m(const std::string precision) {

	if (precision != "single" && precision != "double" && precision != "quadruple") {
		
		if (logging) {
			output("String entered in E_m must be single, double or quadruple");
		}

		return 0;
	
	}
	
	return 1-E_M(precision);

}

//return maximum machine value (Holmes Chapter 1.2.2)
double x_M(const std::string precision) {

	if (precision != "single" && precision != "double" && precision != "quadruple") {
		
		if (logging) {
			output("String entered in x_M must be single, double or quadruple");
		}

		return 0;
	
	}

	return pow(2, E_M(precision))*(1-(eps(precision)/2))*2;

}

//return minimum machine value (Holmes Chapter 1.2.2)
double x_m(const std::string precision) {

	if (precision != "single" && precision != "double" && precision != "quadruple") {
		
		if (logging) {
			output("String entered in x_m must be single, double or quadruple");
		}

		return 0;
	
	}

	return pow(2, E_m(precision));

}

//return absolute error (Holmes Chapter 1.4)
//NOTE: exact solution must be known
template <class S>
S abs_err(const S exact, const S computed) {

	return std::abs(exact-computed);

}

//return relative error (Holmes Chapter 1.4)
//NOTE: exact solution must be known
template <class S>
S rel_err(const S exact, const S computed) {

	S tol = eps("double");

	if (std::abs(exact) <= tol) {
		std::cout << "Error" << '\n';
	}
	else {
		return std::abs(exact - computed)/std::abs(exact);
	}
	
} 

//return error vector (Holmes Chapter 3.6)
//NOTE: exact solution must be known
template <class S>
S error_vec(const S exact, const S computed) {

	return exact - computed;

}

//return residual vector (Holmes Chapter 3.6)
//NOTE: exact solution must be known
template <class S, class T>
S residual_vec(const S b, const T A, const S computed) {

	return b - (A * computed);

}