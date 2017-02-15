#include "interpolation.h"

//summation reducing error (Holmes Chapter 1.3)
template <class S>
double compensated_summation(const S x) {

	if (x.empty()) {

		if(logging) {
			output("Vector entered in compensated_summation has not been initialized");
		}

		return 0;

	}

	double sum = 0;
	double err = 0;
	double z, q;

	for (int i = 0; i < x.size(); ++i) {

		z = x[i] + err;
		q = sum;
		sum = q + z;
		err = (q-sum) + z;

	}

	return sum;
}

//Horner's method for polynomial evaluation (Holmes Chapter 1.3)
template <class S>
double horner(const S a, const double x) {

	if (a.empty()) {

		if(logging) {
			output("Vector entered in horner has not been initialized");
		}

		return 0;

	}

	int n = a.size();
	double p = a[n-1];

	for (int i = 2; i <= n; ++i) {
		p = a[n-i] + p*x;
	}

	return p;
}

////root finding methods////

//bisection method (Holmes Chapter 2.3)
template <class S>
double bisection_method(const S x, const S y, const double lb, const double ub, const double precision) {

	if (x.empty() || y.empty()) {

		if (logging) {
			output("One or more vectors entered in bisection_method have not been initialized");
		}

		return 0;

	}

	if (lb < x(0) || ub > x(x.size()-1)) {

		if (logging) {
			output("Arguments out of bounds");
		}

		return 0;

	}

	double tol = 1e-7;      //don't drop below this value, runtime increases

	int a = get_ID(x, lb, precision);
	int b = get_ID(x, ub, precision);

	if ((a>b) || (y[a]*y[b] > 0)) {

		if (logging) {
			output("Pick different points, criterion not met in bisection_method");
		}

		return 0;

	}

	int c;
	double err = (x[b]-x[a])/2;
	
	while (err > tol) {

		c = (a+b)/2;

		if (std::abs(y[c]) < tol) {

			return x[c];
			break;
		
		}
		else if (y[a]*y[c] < 0) {

			b = c;
			a = a;
		
		}
		else {
		
			a= c;
			b = b;
		
		}
		
		err = (x[b] - x[a])/2;
	
	}

}

//newton's method (Holmes Chapter 2.4)
template <class S>
double newtons_method(const S x, const S y, const S dy, double guess, const double precision) {

	if (x.empty() || y.empty()) {

		if (logging) {
			output("One or more vectors entered in newtons_method have not been initialized");
		}

		return 0;

	}

	if (precision <= 0) {

		if (logging) {
			output("Precision argument in newtons_method must be greater than zero");
		}

		return 0;

	}

	double tol = 1e-6;      //don't drop below this value, runtime increases

	int a = 0;
	int b = x.size()-1;

	if ((x[a] > guess) || (guess > x[b]) || (contains_zero(dy, tol) == 1)) {

		if (logging) {
			output("Criterion not met in newtons_method; pick different points");
		}
		
		return 0;

	}

	int g;
	double err = 1;
	int I = 100;
	int i = 0;
	double z;

	while (err > tol) {

		g = get_ID(x, guess, precision);
		z = y[g]/dy[g];
		err = std::abs(z);
		guess = guess - z;
		i = i + 1;

		if ((guess < x[a]) || (x[b] < guess) || (I < i)) {
			break;
		}

	}

	return guess;

}

//secant method (Holmes Chapter 2.5)
template <class S>
double secant_method(const S x, const S y, double guess1, double guess2, const double precision) {

	if (x.empty() || y.empty()) {

		if (logging) {
			output("One or more vectors entered in secant_method have not been initialized");
		}

		return 0;

	}

	if (precision <= 0) {

		if (logging) {
			output("Precision argument in secant_method must be greater than zero");
		}

		return 0;

	}

	int a = 0;
	int b = x.size()-1;

	if ((x(a) > guess1) || (guess1 > x(b)) || (x(a) > guess2) || (guess2 > x(b))) {
		
		if (logging) {
			output("Criterion not met in secant_method; pick different points");
		}
		
		return 0;

	}

	int g1, g2;
	int i = 0;	
	int I = 200;
	double err = 1;
	double tol = 1e-7;      //don't drop below 1e-7, runtime increases incredibly, or shit breaks down, unsure
	double z;

	while (err > tol) {

		g1 = get_ID(x, guess1, precision);
		g2 = get_ID(x, guess2, precision);
		z = y(g1)*((guess1-guess2)/(y(g1)-y(g2)));

		if (std::abs(y(g1)-y(g2)) < tol) {
			break;
		}

		err = std::abs(z);
		guess2 = guess1;
		guess1 = guess1 - z;
		i = i + 1;

		if ((guess1 < x(a)) || (x(b) < guess1) || (I < i)) {
			break;
		}

	}


	return guess2;

}

////sorting methods////

//insertion sort
template <class S>
void insertion_sort(S A[], const int n, const std::string order) {

	int i, j;
	double val;

	if (order == "forward" || order == "F") {

		for (j = 1; j < n; ++j) {

			val = A[j];
			i = j - 1;

			while (i >= 0 && A[i] > val) {

				A[i+1] = A[i];
				i = i - 1;

			}

			A[i+1]=val;

		}

	}
	else if (order == "backward" || order == "B") {

		for (j = 1; j < n; ++j) {

			val = A[j];
			i = j - 1;

			while (i >= 0 && A[i] < val) {

				A[i+1] = A[i];
				i = i - 1;

			}

			A[i+1]=val;

		}

	}

}

//insertion sort
void insertion_sort(arma::umat& A, const std::string order, const bool reduced = 0) {

	if (A.is_empty()) {

		if (logging) {
			output("Matrix entered in insertion_sort has not been initialized");
		}

		return;

	}

	unsigned int count = 0;
	bool swap;
	int i;
	double val;
	arma::urowvec temp;

	if (order == "forward" || order == "F") {

		for (int j = 1; j < A.n_rows; ++j) {

			temp = A.row(j);
			i = j - 1;

			while (i >= 0) {

				for (int k = 0; k < A.n_cols; ++k) {

					val = temp(k);

					if (A(i,k) < val) {

						swap = 0;

						break;

					}
					else if (A(i, k) > val) {

						A.row(i+1) = A.row(i);
						A.row(i) = temp;
						swap = 1;

						break;

					}

				}

				i = i - 1;

			}

			if (swap) {
				A.row(i+1) = temp;
			}

		}

	}
	else if (order == "backward" || order == "B") {

		for (int j = 1; j < A.n_rows; ++j) {

			temp = A.row(j);
			i = j - 1;

			while (i >= 0) {

				for (int k = 0; k < A.n_cols; ++k) {

					val = temp(k);

					if (A(i,k) > val) {

						swap = 0;

						break;

					}
					else if (A(i, k) < val) {

						A.row(i+1) = A.row(i);
						A.row(i) = temp;
						swap = 1;

						break;

					}

				}

				i = i - 1;

			}

			if (swap) {
				A.row(i+1) = temp;
			}

		}

	}

	if (reduced) {

		count = 0;

		for (int j = 1; j < A.n_rows; ++j) {

			if (arma::sum(A.row(j-1) == A.row(j)) == A.n_cols) {

				A.shed_row(j-1);
				j -= 1;

			}

		}

	}

}

////number theory algorithms////

//Euclid's algorithm (Hall Chapter 4.5)
template <class S>
S Euclid(const S a, const S b){

	if (b <= 0){

		if (logging) {
			output("Second argument in Euclid must be greater than zero");
		}

		return 0;
	}

	int a_n, b_n, q_n, r_n, h;

	a_n = a;
	b_n = b;
	q_n = a_n/b_n;
	r_n = a - q_n*b;

	if (r_n < 0){           //ensure positivity
		r_n = 0 - r_n;
	}

	if (r_n == 0){          //determine whether a is a multiple of b
		h = b_n;
	}

	while (r_n >= 1) {

		h = r_n;
		a_n = b_n;
		b_n = r_n;
		q_n = a_n/b_n;
		r_n = a_n - q_n*b_n;

		if (r_n <= 0){
			break;
		}

	}

	if (h < 0){          //ensure positivity of h
		h = 0-h;
	}

	return h;

}

//prime sieve
template <class S>
void sieve(const S n) {

	if (n <= 0){

		if (logging) {
			output("Value entered in sieve must be greater than zero");
		}

		return;
	}

	bool A[n-1];

	for (int i = 0; i < n-1; ++i) {
		A[i] = 1;
	}

	for (int i = 0; i <= floor(sqrt(n))-1; ++i) {

		if(A[i] == 1) {

			for (int j = (i+2)*(i+2)-2; j <= n-1; j += (i+2)) {
				A[j] = 0;
			}

		}

	}

	for (int i = 0; i < n; ++i) {

		if(A[i] == 1) {
			output(i+2);
		}

	}

}