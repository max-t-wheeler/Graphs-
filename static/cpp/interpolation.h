#include "arma_lin_alg_package.h"

//interpolation using a global interpolation polynomial (Holmes Chapter 5.2)
//Note: Requires unique x(i) values
arma::vec global_interpolation_poly(const arma::vec x, const arma::vec y, arma::vec& D, const double range, const double precision) {

	if (x.empty() || y.empty()) {

		if(logging) {
			output("One or more vectors entered in global_interpolation_poly have not been initialized");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (range < 0 || precision < 0) {

		if(logging) {
			output("The range and precision arguments entered in global_interpolation_poly must be greater than or equal to zero");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (x.size() != y.size()) {

		if(logging) {
			output("Error: vectors in global_interpolation_poly are not the same length.");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;
	
	}

	D = arma::regspace(-range, precision, range);

	arma::vec l = arma::ones(D.size());
	arma::vec p = arma::zeros(D.size());

	for (int i = 0; i < x.size(); ++i) {

		l.ones();

		for (int j = 0; j < x.size(); ++j) {			
			if (i != j) {
				l %= (D - x(j))/(x(i) - x(j));
			}	
		}

		p += (y(i) * l);
	
	}

	return p;

}

//chebyshev interpolation (Holmes Chapter 5.5.3)
//Note: Requires unique x(i) values such that x(i) = (a + b + ((b-a)*z))/2 where z = cos((2*i+1)*M_PI/(2*n))
//Note: Suggested that range = (b-a)/2
arma::vec chebyshev_interpolation(const arma::vec x, const arma::vec y, arma::vec& D, const unsigned int n) {

	if (x.empty() || y.empty()) {

		if(logging) {
			output("One or more vectors entered in chebyshev_interpolation have not been initialized");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (x.size() != y.size()) {
		
		if (logging) {
			output("Error: vectors in chebyshev_interpolation are not the same length.");
		}

		D.zeros(1);
		
		return arma::zeros(1);
	}

	D.set_size(n+1);

	for (int i = 0; i < D.size(); ++i) {
		D(i) = cos((2*i-1)*M_PI/(2*n));
	}

	arma::vec l = arma::ones(D.size());
	arma::vec p = arma::zeros(D.size());

	for (int i = 0; i < x.size(); ++i) {

		l.ones();
		
		for (int j = 0; j < x.size(); ++j) {			
			if (i != j) {
				l %= (D - x(j))/(x(i) - x(j));
			}	
		}

		p += (y(i) * l);
	
	}

	return p;

}

//linear piecewise interpolation (Holmes Chapter 5.3)
//Note: Requires unique x(i) values
//NOTE: Necessary to include a point at each end to have proper definition of G_1, G_n+1
arma::vec linear_piecewise_interpolation(const arma::vec x, const arma::vec y, arma::vec& D, const double precision) {

	if (x.empty() || y.empty()) {

		if(logging) {
			output("One or more vectors entered in linear_piecewise_interpolation have not been initialized");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (precision < 0) {

		if(logging) {
			output("The precision argument entered in linear_piecewise_interpolation must be greater than or equal to zero");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (x.size() != y.size()) {
		
		if (logging) {
			output("Error: vectors in linear_piecewise_interpolation are not the same length.");
		}

		D.zeros(1);
		
		return arma::zeros(1);
	}

	double h = x(1) - x(0);
	arma::vec X = arma::regspace(x(0) - h, h, x(x.size()-1) + h);

	D = arma::regspace(X(0), precision, X(X.size()-1));
	
	arma::vec g = arma::zeros(D.size());
	arma::mat G = arma::zeros(D.size(), X.size()); 

	for (int j = 1; j < G.n_cols-1; ++j) {

		for (int i = 0; i < G.n_rows; ++i) {

			if (X(j-1) <= D(i) && D(i) <= X(j)) {
				G(i, j) = (D(i) - X(j-1))/(X(j) - X(j-1));
			}

			if(X(j) <= D(i) && D(i) <= X(j+1)) {
				G(i, j) = (D(i) - X(j+1))/(X(j) - X(j+1));
			}

		}

		g += (y(j-1) * G.col(j));

	}

	return g;

}

//cubic piecewise interpolation using natural cubic B splines (Holmes Chapter 5.4)
//Note: Requires unique x(i) values
//NOTE: Necessary to include two points at each end to have proper definition of B_0, B_n+2
arma::vec cubic_piecewise_interpolation(const arma::vec x, const arma::vec y, arma::vec& D, const double precision) {

	if (x.empty() || y.empty()) {

		if(logging) {
			output("One or more vectors entered in cubic_piecewise_interpolation have not been initialized");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (precision < 0) {

		if(logging) {
			output("The precision argument entered in cubic_piecewise_interpolation must be greater than or equal to zero");
		}

		D.zeros(1);

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (x.size() != y.size()) {
		
		if (logging) {
			output("Error: vectors in cubic_piecewise_interpolation are not the same length.");
		}

		D.zeros(1);
		
		return arma::zeros(1);
	}

	if (y.size() < 5) {

		if (logging){
			output("Warning. Data set entered in cubic_pieceswise_interpolation is has less than 5 datapoints.");
		}

	}

	//define a
	arma::mat A = tri_diag_mat(4, 1, 1, y.size()-2, y.size()-2);
	arma::vec a = arma::ones(A.n_rows+4);
	arma::vec z(A.n_rows);		

	for (int i = 0; i < z.size(); ++i) {
		z(i) = 6*y(i+1);
	}

	z(0) -= y(0);
	z(z.size()-1) -= y(y.size()-1);

	arma::vec middle = thomas(A, z);

	for (int i = 0; i < middle.size(); ++i) {
		a(i+2) = middle(i);
	}

	a(1) = y(0);
	a(0) = 2*a(1) - a(2);
	a(a.size()-2) = y(y.size()-1);
	a(a.size()-1) = 2*a(a.size()-2) - a(a.size()-3);

	//define B
	double h = x(1) - x(0);
	arma::vec X = arma::regspace(x(0) - h, h, x(x.size()-1) + h);
	D = arma::regspace(X(0), precision, X(X.size()-1));
	arma::vec s = arma::zeros(D.size());
	arma::mat B = arma::zeros(D.size(), X.size()); 

	for (int j = 0; j < B.n_cols; ++j) {

		arma::vec v = (D - X(j))/h;

		for (int i = 0; i < B.n_rows; ++i) {

			if (std::abs(v(i)) <= 1) {
				B(i, j) = 2.0/3.0 - (pow(v(i), 2)*(1-(std::abs(v(i))/2)));
			}
			else if (1 <= std::abs(v(i)) && std::abs(v(i)) <= 2) {
				B(i, j) = pow((2-std::abs(v(i))), 3)/6;
			}
			else {
				B(i, j) = 0;
			}

		}

		s += (a(j) * B.col(j));
	
	}

	return s;

}