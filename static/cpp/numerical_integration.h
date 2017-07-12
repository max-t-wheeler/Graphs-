#include "arma_algos_package.h"

////////////////////methods based on interpolation////////////////////////

//numerical integration using the composite midpoint rule (Holmes Chapter 6.2.1)
double composite_midpoint_integration(const double a, const double b, const int n, arma::vec (*f)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in composite_midpoint_integration; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in composite_midpoint_integration must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);

	double ans = 0;

	arma::vec x = arma::regspace(a, h, b);
	x += h/2;

	arma::vec y = f(x);

	for (int i = 0; i < x.size()-1; ++i) {
		ans += y(i);
	}

	ans *= h;

	return ans;

}

//numerical integration using the composite midpoint rule (Holmes Chapter 6.3.4)
//NOTE: requires knowledge of first derivative of ys
double corrected_composite_midpoint_integration(const double a, const double b, const int n, arma::vec (*f)(arma::vec), arma::vec (*f_p)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in corrected_composite_midpoint_integration; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in corrected_composite_midpoint_integration must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);

	double ans = 0;

	arma::vec x = arma::regspace(a, h, b);

	x += h/2;

	arma::vec y = f(x);
	arma::vec y_p = f_p(x);

	for (int i = 0; i < x.size()-1; ++i) {
		ans += y(i);
	}

	ans *= h;

	ans -= pow(h, 2)*(y_p(0) - y_p(y_p.size()-1))/24;

	return ans;

}

//numerical integration using the composite trapezoidal rule (Holmes Chapter 6.3.1)
double composite_trapezoidal_integration(const double a, const double b, const int n, arma::vec (*f)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in composite_trapezoidal_integration; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in composite_trapezoidal_integration must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);

	arma::vec x = arma::regspace(a, h, b);
	arma::vec y = f(x);

	double ans = y(0)/2;

	for (int i = 1; i < y.size()-1; ++i) {
		ans += y(i);
	}

	ans += y(y.size()-1)/2;

	ans *= h;

	return ans;

}

//numerical integration using the composite simpson's rule (Holmes Chapter 6.3.2)
double composite_simpson_integration(const double a, const double b, const int n, arma::vec (*f)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in composite_simpson_integration; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in composite_simpson_integration must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);

	arma::vec x = arma::regspace(a, h, b);
	arma::vec y = f(x);

	double ans = y(0);

	for (int i = 1; i < y.size()-1; ++i) {

		if (i%2 == 0) {
			ans += 2*y(i);
		}
		else {
			ans += 4*y(i);
		}

	}

	ans += y(y.size()-1);

	ans *= h/3;

	return ans;

}

//numerical integration using the composite hermite rule (Holmes Chapter 6.3.3)
//NOTE: requires knowledge of first derivative of y
double composite_hermite_integration(const double a, const double b, const int n, arma::vec (*f)(arma::vec), arma::vec (*f_p)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in composite_hermite_integration; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in composite_hermite_integration must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);

	arma::vec x = arma::regspace(a, h, b);
	arma::vec y = f(x);
	arma::vec y_p = f_p(x);

	//composite trapezoidal rule
	double ans = y(0)/2;

	for (int i = 1; i < y.size()-1; ++i) {
		ans += y(i);
	}

	ans += y(y.size()-1)/2;

	ans *= h;

	//adjustment
	ans += pow(h, 2)*(y_p(0) - y_p(y_p.size()-1))/12;

	return ans;

}

/////////////////////methods based on precision//////////////////////

//2-point gaussian quadrature (Holmes Chapter 6.4.2)
double two_point_guass_quad(const double a, const double b, const double n, arma::vec (*f)(arma::vec)) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in two_point_guass_quad; adjust bounds");
		}

		return 0;

	}

	if (n <= 0) {

		if (logging) {
			output("Third argument in two_point_guass_quad must be greater than zero");
		}

		return 0;

	}

	double h = (b-a)/(n-1);
	double ans = 0;

	arma::vec x = arma::regspace(a, h, b);
	arma::vec z_plus(x.size()-1);
	arma::vec z_minus(x.size()-1);

	for (int i = 0; i < x.size()-1; ++i) {

		z_plus(i) = ((h/sqrt(3)) + x(i) + x(i+1))/2;
		z_minus(i) = (-(h/sqrt(3)) + x(i) + x(i+1))/2;

	}

	arma::vec f_plus = f(z_plus);
	arma::vec f_minus = f(z_minus);

	for (int i = 0; i < x.size()-1; ++i) {
		ans += f_plus(i) + f_minus(i);
	}

	ans *= h/2;

	return ans;

}

//iterative integration algorithm (Holmes Table 6.9)
double iterative_integration(const double a, const double b, arma::vec (*f)(arma::vec), const double tol, const std::string type) {

	if (a >= b) {

		if (logging) {
			output("Second argument must be less than first argument in iterative_integration; adjust bounds");
		}

		return 0;

	}

	if (tol <= 0) {

		if (logging) {
			output("Precision argument in iterative_integration must be greater than zero");
		}

		return 0;

	}

	double A, B, R;

	unsigned int k = 1;

	if (type == "simpson") {

		A = 0;
		R = 1;

		while(std::abs(R-A) >= tol) {

			A = composite_simpson_integration(a, b, pow(2, k), f);
			B = composite_simpson_integration(a, b, pow(2, k+1),  f);
			R = ((16*B) - A)/15;

			k++;

		}

		return R;

	}
	else if (type == "gaussquad2") {
		
		A = 0;
		R = 1;
		
		while(std::abs(R-A) >= tol) {

			A = two_point_guass_quad(a, b, pow(2, k), f);
			B = two_point_guass_quad(a, b, pow(2, k+1), f);
			R = ((16*B) - A)/15;

			k++;

		}

		return R;

	}
	else {

		if (logging) {
			output("Integration type in iterative_integration must be simpson or gaussquad2");
		}

		return 0;

	}

}