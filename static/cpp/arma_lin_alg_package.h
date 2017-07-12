#include "arma_scientific_computing.h"

//matrices

//determine if matrix can be factored using LU decomposition (Holmes Chapter 3.2.1)
template <class S>
bool is_LU(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_LU has not been initialized");
		}

		return 0;

	}

	if (A.n_rows == A.n_cols) {
		return 1;
	}
	else {
		return 0;
	}

}

//determine if matrix is tridiagonal
template <class S>
bool is_tridiagonal(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_tridiagonal has not been initialized");
		}

		return 0;

	}

	unsigned int count = 0;

	for (int i = 0; i < A.n_rows-2; ++i) {
		for (int j = i+2; j < A.n_cols; ++j) {
			if (A(i,j) != 0) {
				count++;
			}
		}
	}

	for (int i = 2; i < A.n_rows; ++i) {
		for (int j = i-2; j < A.n_cols-2; ++j) {
			if (A(i,j) != 0) {
				count++;
			}
		}
	}

	if(count > 0) {
		return 0;
	}

	else {
		return 1;
	}

}

//determine if matrix is symmetric
template <class S>
bool is_symmetric(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_symmetric has not been initialized");
		}

		return 0;

	}

	if (A.n_rows != A.n_cols) {

		if (logging) {
			output("The matrix entered into is_symmetric is not nxn");
		}
		
		return 0;
	}
	else {

		unsigned int count = 0;

		for (int i = 0; i < A.n_rows-1; ++i) {
			for (int j = i+1; j < A.n_cols; ++j) {
				if (A(i,j) != A(j,i)) {
					count++;
				}
			}
		}

		if(count == 0) {
			return 1;
		}
		else {
			return 0;
		}

	}

}

//determine if matrix is strictly diagonal dominant (Holmes Chapter 3.7.1)
template <class S>
bool is_strictly_diagonal_dominant(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_strictly_diagonal_dominant has not been initialized");
		}

		return 0;

	}

	unsigned int count = 0;

	for (int i = 0; i < A.n_rows; ++i) {

		unsigned int sum = 0;
		
		for (int j = 0; j < A.n_cols; ++j) {
			if (i != j) {
				sum += std::abs(A(i,j));
			}
		}
		if (std::abs(A(i,i)) > sum) {
			count++;
		}
	}

	if (count == A.n_rows) {
		return 1;
	}
	else {
		return 0;
	}

}

//determine if matrix is strictly positive definite
template <class S>
bool is_strictly_positive_definite(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_strictly_positive_definite has not been initialized");
		}

		return 0;

	}

	if (!is_symmetric(A)) {

		if (logging) {
			output ("Matrix entered into is_strictly_positive_definite is not symmetric");
		}
		
		return 0;
	}
	else if (!is_strictly_diagonal_dominant(A)) {

		if (logging) {
			output ("Matrix entered into is_strictly_positive_definite is not strictly diagonal dominant");
		}
		
		return 0;
	}
	else {
		return 1;
	}

}

//determine if matrix is positive definite (Holmes Chapter 3.7.1)
template <class S>
bool is_positive_definite(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in is_positive_definite has not been initialized");
		}

		return 0;

	}

	unsigned int count = 0;

	for (int i = 0; i < A.n_rows; ++i) {
		if (A(i,i) <= 0) {
			return 0;
		}
	}

	if(is_symmetric(A) && is_strictly_positive_definite(A)) {
		return 1;
	}	
	else {

		if (logging) {
			output("Unable to determine if matrix entered in is_positive_definite is positive definite.");
		}

		return 0;
	
	}

}

//determine if matrix might be positive definite (Holmes Chapter 3.7.1)
template <class S>
bool might_be_positive_definite(S A) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in might_be_positive_definite has not been initialized");
		}

		return 0;

	}

	if (!is_symmetric(A)) {

		if (logging) { 
			output ("Matrix entered into might_be_positive_definite is not symmetric");
		}

		return 0;

	}
	else {
		unsigned int count = 0;
		unsigned int tally = 0;

		for (int i = 0; i < A.n_rows; ++i) {
			if (A(i, i) <= 0) {
				tally++;
			}
		}

		if (tally == 0) {
			count++;
		}

		for (int i = 0; i < A.n_rows-1; ++i) {
			for (int j = 1; j < A.n_rows; ++j) {
				if (A(i,j) >=  arma::max(arma::max(A))) {
					count++;
				}
			}
		}

		if (count == 1) {
			return 1;
		}
		else {
			return 0;
		}

	}

}

//create vandermonde matrix
arma::mat vandermonde_matrix(arma::vec x, unsigned int m, unsigned int n) {

	if (x.empty()) {

		if (logging) {
			output("Vector entered in vandermonde_matrix has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);

		return exit;

	}

	arma::mat A(m, n);

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			A(i, j) = pow(x(i), j);
		}
	}

	return A;

}

//create tridiagonal matrix
arma::mat tri_diag_mat(double a, double b, double c, unsigned int m, unsigned int n) {

	arma::mat A = arma::zeros(m, n);

	for (int i = 0; i < A.n_rows-1; ++i) {

		A(i, i) = a;
		A(i+1, i) = b;
		A(i, i+1) = c;

	}

	A(A.n_rows-1, A.n_cols-1) = a;

	return A;

}

//Gram-Schmidt procedure - used to transform a list of linearly independent vectors into an orthonormal basis (Holmes Chapter 4.3.2)
//NOTE: columns in C must be linearly independent
arma::mat gram_schmidt(arma::mat C) {

	if (C.empty()) {

		if (logging) {
			output("Matrix entered in gram_schmidt has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);

		return exit;

	}

	arma::mat E = arma::zeros(C.n_rows, C.n_cols);

	E.col(0) = C.col(0);

	for (int j = 1; j < C.n_cols; ++j) {

		E.col(j) = C.col(j);

		for (int k = 0; k < j; ++k) {

			double r = arma::dot(C.col(j), E.col(k))/arma::dot(E.col(k), E.col(k));

			E.col(j) -= (r * E.col(k));

		}

	}

	for (int j = 0; j < C.n_cols; ++j) {
		E.col(j) = E.col(j)/arma::norm(E.col(j));
	}

	return E;

}

//modified Gram-Schmidt procedure - used to transform a list of linearly independent vectors into an orthonormal basis (Holmes Chapter 4.3.2)
//NOTE: columns in C must be linearly independent
arma::mat modified_gram_schmidt(arma::mat C) {

	if (C.empty()) {

		if (logging) {
			output("Matrix entered in modified_gram_schmidt has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);

		return exit;

	}

	arma::mat E = arma::zeros(C.n_rows, C.n_cols);
	
	E.col(0) = C.col(0);

	for (int j = 1; j < C.n_cols; ++j) {

		E.col(j) = C.col(j);

		for (int k = 0; k < j; ++k) {

			double r = arma::dot(E.col(j), E.col(k))/arma::dot(E.col(k), E.col(k));

			E.col(j) -= (r * E.col(k));

		}

	}

	for (int j = 0; j < C.n_cols; ++j) {
		E.col(j) = E.col(j)/arma::norm(E.col(j));
	}

	return E;

}

//Thomas algorithm (Holmes Table 3.7)                    
//NOTE: A must be tridiagonal
arma::vec thomas(arma::mat A, arma::vec z) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in thomas has not been initialized");
		}

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (z.empty()) {

		if (logging) {
			output("Vector entered in thomas has not been initialized");
		}

		arma::vec exit = arma::zeros(1);

		return exit;

	}

	if (A.n_rows != A.n_cols) {
		
		if(logging) {
			output("The matrix entered into thomas is not nxn.");
		}

		arma::vec exit = arma::zeros(1);

		return exit;
	}
	else {

		if (A(0,0) == 0) {
			
			if (logging) {
				output("In thomas, the first entry in A is nonzero. The algorithm will likey error out.");
			}

		}

		arma::vec x(A.n_rows, arma::fill::zeros);
		arma::vec v(A.n_rows);

		double w = A(0, 0);

		x(0) = z(0)/w;

		for (int i = 1; i < A.n_rows; ++i) {

			v(i) = A(i-1, i)/w;
			w = A(i, i) - (A(i, i-1)*v(i));
			x(i) = (z(i) - (A(i, i-1)*x(i-1)))/w;

		}

		for (int j = A.n_rows-2; j >= 0; j--) {
			x(j) = x(j) - (v(j+1)*x(j+1));
		}

		return x;

	}

}

////Eigen-decomposition and factorization methods////

//power method - used to solve for the greatest eigenvalue (by absolute value) (Holmes Table 4.1)
double power_method(arma::mat A, arma::vec y, double tol, unsigned int K) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in power_method has not been initialized");
		}

		return 0;

	}

	if (y.empty()) {

		if (logging) {
			output("Vector entered in power_method has not been initialized");
		}

		return 0;

	}

	if(!is_symmetric(A)) {
		
		if (logging) {
			output("No guarantee eigenvalues are real in power_method");
		}

		return 0;

	}
	else if (A.n_rows != y.size()) {

		if(logging) {
			output("Matrix and guess vector not the same length in power_method");
		}

		return 0;
		
	}

	arma::vec z = A * y;

	double v0 = arma::dot(y, z)/arma::dot(y, y);
	double vk, vk1;

	for(int k = 0; k < K; ++k) {

		vk = v0;
		y = z/arma::norm(z);
		z = A * y;
		vk1 = arma::dot(y, z);

		if (abs_err(vk1,vk)/std::abs(vk1) < tol) {
			break;
		}

	}

	return vk1;

}

//inverse power method - used to solve for the eigenvalue closest to zero (Holmes Table 4.4)
double inv_power_method(arma::mat A, arma::vec y, double tol, unsigned int K) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in inv_power_method has not been initialized");
		}

		return 0;

	}

	if (y.empty()) {

		if (logging) {
			output("Vector entered in inv_power_method has not been initialized");
		}

		return 0;

	}

	if(!is_symmetric(A)) {
		
		if (logging) {
			output("No guarantee eigenvalues are real in inv_power_method");
		}

		return 0;

	}
	else if (A.n_rows != y.size()) {

		if(logging) {
			output("Matrix and guess vector not the same length in inv_power_method");
		}

		return 0;

	}

	arma::vec z = solve(A, y);
	double v0 = arma::dot(y, z)/arma::dot(y, y);
	double vk, vk1;

	for(int k = 0; k < K; ++k) {

		vk = v0;
		y = z/arma::norm(z);
		z = solve(A, y);
		vk1 = arma::dot(y, z);

		if (abs_err(vk1,vk)/std::abs(vk1) < tol) {
			break;
		}

	}

	return vk1;

}

//orthogonal iteration - solve for a matrix's n greatest eigenvalues (evalorevec=0) or eigenvectors (evalorevec=1) (Holmes Table 4.5)
//NOTE: if the columns of B are not orthogonal, this algorithm will yield erroneous results
arma::mat orthogonal_iteration(arma::mat A, unsigned int m, unsigned int n, unsigned int K, bool evalorevec) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in orthogonal_iteration has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);
		
		return exit;

	}

	arma::mat B = arma::randu(m, n);
	arma::mat Q;

	Q = modified_gram_schmidt(B);

	for (int k = 0; k < K; ++k) {

		B = A * Q;

		Q = modified_gram_schmidt(B);

	}

	if (!evalorevec) {

		arma::vec evals(n);

		for (int i = 0; i < n; ++i) {
			evals(i) = arma::dot(Q.col(i), A*Q.col(i))/arma::dot(Q.col(i), Q.col(i));
		}

		return evals;

	}
	else {
		return Q;
	}

}

//inverse orthogonal iteration - solve for a matrix's n least eigenvalues (evalorevec=0) or eigenvectors (evalorevec=1) (Holmes Table 4.10)
//NOTE: if the columns of B are not orthogonal, this algorithm will yield erroneous results
arma::mat inv_orthogonal_iteration(arma::mat A, int m, int n, int K, bool evalorevec) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in inv_orthogonal_iteration has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);
		
		return exit;

	}

	arma::mat B = arma::randu(m, n);
	arma::mat L, U, Q;

	Q = modified_gram_schmidt(B);

	arma::lu(L, U, A);

	for (int k = 0; k < K; ++k) {
		B = solve(L*U, Q);
		Q = modified_gram_schmidt(B);
	}

	if (!evalorevec) {

		arma::vec evals(n);

		for (int i = 0; i < n; ++i) {
			evals(i) = arma::dot(Q.col(i), A*Q.col(i))/arma::dot(Q.col(i), Q.col(i));
		}

		return evals;
	}

	else {
		return Q;
	}
}

//solve for a matrix's eigenvalues (Holmes Chapter 4.3.4)
arma::mat QR_method(arma::mat A, unsigned int K) {

	if (A.empty()) {

		if (logging) {
			output("Matrix entered in QR_method has not been initialized");
		}

		arma::mat exit = arma::zeros(1, 1);
		
		return exit;

	}

	arma::mat B(A.n_rows, A.n_cols);
	arma::mat R(A.n_rows, A.n_cols);
	arma::mat Q(A.n_rows, A.n_cols);

	B = A;

	for (int k = 0; k < K; ++k) {
		Q = modified_gram_schmidt(B);
		R = Q.t() * B;
		B = R * Q;
	}

	return B;

}