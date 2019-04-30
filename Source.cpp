#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "mvector.h"

// function implementing the steepest descent algorithm, returning the number of 
// iterations needed for x to a solution (i.e. number of calculations needed for 
// the residual to become sufficiently small)
int SDLS(const MMatrix& A, const MVector& b, MVector& x0, 
	int maxIterations = 1000, double tol = 1.0e-6, bool output = false)
{
	if (output) // if output==true, write x values for each iteration to a file 
	{
		std::ofstream sdlsoutput;
		sdlsoutput.open("sdls.csv");
		MMatrix At = A.transpose(); // transpose of A
		MVector r = At * (b - A * x0); // residual at x
		double alpha; // step length
		for (int iter = 0; iter <= maxIterations; iter++)
		{
			if (r.dot(r) < tol*tol) return iter; // if ||r|| sufficiently small
												 // return no. of iterations
			else // calculate values for next iteration
			{
				alpha = r.dot(r) / (A*r).dot(A*r);
				x0 = x0 + r.scale(alpha);
				sdlsoutput << x0 << "\n";
				r = At*(b - A*x0);
			}
		}
		return -1; // if it fails to find a soln within the no. of maxIterations,
				   // return -1
	}
	else
	{
		MMatrix At = A.transpose(); // transpose of A
		MVector r = At * (b - A * x0); // residual at x
		double alpha; // step length
		for (int iter = 0; iter <= maxIterations; iter++)
		{
			if (r.dot(r) < tol*tol) return iter; // if ||r|| sufficiently small
												 // return no. of iterations
			else // calculate values for next iteration
			{
				alpha = r.dot(r) / (A*r).dot(A*r);
				x0 = x0 + r.scale(alpha);
				r = At * (b - A * x0);
			}
		}
		return -1; // if it fails to find a soln within the no. of maxIterations,
				   // return -1
	}
}

// function implementing the normalised iterative hard thresholding algortihm,
// returning the number of iterations needed for x to converge to a solution (i.e.
// no. of calculations needed for residual to become sufficiently small)
int NIHT(const MMatrix& A, const MVector& b, MVector&x0, int k,
	int maxIterations=1000, double tol=1.0e-6)
{

	// initialise starting vector
	x0 = A.transpose()*b;
	x0.threshold(k);
	MMatrix At = A.transpose();
	MVector r = At * (b - A * x0), r_ip1(x0.size());
	double alpha, change_in_r;
	for (int iter = 0; iter <= maxIterations; iter++)
	{
		if (r.dot(r) < tol*tol) { return iter; }
		else
		{
			alpha = r.dot(r) / (A*r).dot(A*r);
			x0 = x0 + At.scale(alpha)*(b - A * x0);
			x0.threshold(k);
			r_ip1 = At * (b - A * x0);
			change_in_r = abs(r.dot(r) - r_ip1.dot(r_ip1)); // change in residuals
			if (change_in_r < 1e-16) return -1; // if residuals don't sufficiently
												// change, nonconvergence occurs
			else r = r_ip1; // continue with algorithm
		}
	}
	return -1; // return -1 if it does not converge
}

// function that tests whether the residual is changing in value enough between 
// iterations of the NIHT algorithm
int niht_residual_test(const MMatrix& A, const MVector& b, MVector&x0, int k,
	int maxIterations=1000, double tol=1.0e-6)
{
	// initialise starting vector
	x0 = A.transpose()*b;
	x0.threshold(k);
	MMatrix At = A.transpose();
	MVector r = At * (b - A * x0), r_ip1(x0.size());
	double alpha, close;
	for (int iter = 0; iter <= maxIterations; iter++)
	{
		alpha = r.dot(r) / (A*r).dot(A*r);
		x0 = x0 + At.scale(alpha)*(b - A * x0);
		x0.threshold(k);
		r_ip1 = At * (b - A * x0);
		close = abs(r.dot(r) - r_ip1.dot(r_ip1));
		if (close < 1e-16) { return iter; }
		else r = r_ip1;
	}
	return -1; // return -1 if it does not converge
}


// test to check that NIHT is working as expected
/*
int main()
{
	std::srand(std::time(NULL));
	int m = 150, n = 200, k = 10, maxIter = 1000;
	double tol = 1e-6;
	MMatrix A(m, n);
	A.initialise_normal();
	MVector x(n);
	x.initialise_normal(k);
	MVector b = A * x;
	MVector xout(n);
	std::cout << "x = \n" << x << "\n";
	std::cout << NIHT(A,b,xout,k,maxIter,tol) << "\n";
	std::cout << "xout = \n" << xout << "\n";
	x = x - xout;
	if (x.has_converged()) std::cout << "NIHT converged to correct x";
	return 0;
}
*/


// residuals not changing test
int main()
{
	std::srand(std::time(NULL));
	int m, n, k, niht_it, resfail_it;
	for (int it = 1; it <= 200; it++)
	{
		n = rand() % 91 + 10; // random int in [10,100]
		m = rand() % (n-9) + 9; // random int in [9,n-1]
		k = rand() % ((m / 2) - 4) + 4; // random in [4,m/2 - 1]
		std::cout << "(m,n,k) = " << m << "," << n << "," << k << "\n";
		MMatrix A(m, n);
		A.initialise_normal();
		MVector x(n), x0(n), x00(n);
		x.initialise_normal(k);
		MVector b = A * x;
		niht_it = NIHT(A, b, x0, k);
		// if it converges to wrong x, set niht_it to -1
		if (niht_it != -1 && !(x - x0).has_converged()) niht_it = -1;
		resfail_it = niht_residual_test(A, b, x00, k);
		std::cout << "no. of iter to converge = " << niht_it <<
			" | | r stops changing at = " << resfail_it << "\n";
	}
	return 0;
}


// write no. of NIHT iterations to a file
/*
int main()
{
	std::ofstream NIHTtest;
	NIHTtest.open("NIHTtest.csv");
	NIHTtest << "m,n,k,iterations" << "\n";
	double tol = 1e-6, nihtnum;
	int maxIter = 500;
	for (int n = 5; n < 200; n++)
		for (int m = 4; m < n; m++)
			for (int k = 5; 2 * k < m; k++)
			{
				MMatrix A(m, n);
				A.initialise_normal();
				MVector x(n), xout(n);
				x.initialise_normal(k);
				MVector b = A * x;
				nihtnum = NIHT(A, b, xout, k, maxIter, tol);
				if (nihtnum != -1 && !(x - xout).zero())
				{
					std::cout << "||converged to wrong x,   ";
					nihtnum = -1;
				}
				std::cout << "m = " << m << ", n = " << n << ", k = " << k <<
					", nihtnum = " << nihtnum << "\n";
				NIHTtest << m << "," << n << "," << k << "," << nihtnum << "\n";
			}
	return 0;
}
*/

// task for SDLS algortihm to find convergence of x to a minimum for the system
// initialised below
/*
int main()
{
	std::srand(std::time(NULL)); // avoid same seq of random values each time the
								 // the program is run
	MMatrix A(3, 2);
	MVector b(3);
	MVector x0(2);
	// initialise A as in task
	A.set(0, 0, 1.0);
	A.set(1, 0, 2.0);
	A.set(2, 0, 1.8);
	A.set(0, 1, 2.0);
	A.set(1, 1, 1.0);
	A.set(2, 1, -2.0);
	// initialise b as in task
	b[0] = 10;
	b[1] = -1;
	b[2] = 0;

	std::cout << SDLS(A, b, x0, 1000, 1.0e-06, true); // also writes x trajectory
	return 0;										  // to a file
}
*/

// write p(m) values to a file for 4<=m<=199 (steps of 5)
// screen prints are for making sure program is running
int main()
{
	std::srand(std::time(NULL)); // avoid same seq of random values
	std::ofstream niht_successrate; // open file for writing values
	niht_successrate.open("niht_successrate.csv");
	niht_successrate << "m,n,k,p(m)" << "\n";
	int n = 200, k = 20, T = 100; // T the no. of times NIHT is run
	double pm; // pm the successful recovery rate
	for (int m = 4; m < 200; m += 5)
	{
		pm = 0.0;
		for (int t = 0; t < T; t++)
		{
			MVector x(n), x0(n); // x random, x0 for NIHT
			x.initialise_normal(k);
			MMatrix A(m, n);
			A.initialise_normal();
			MVector b = A * x;
			int it = NIHT(A, b, x0, k);
			// compare x and x0 to see if they are really the same
			if (it != -1 && (x-x0).has_converged()) {
				std::cout << "no. of iterations = " << it << "\n";
				pm += 1.0;
			}
			else std::cout << "does not converge \n";
		}
		std::cout << "m = " << m << " is done\n";
		pm /= T;
		niht_successrate << m << "," << n << "," << k << "," << pm << "\n";
	}
	return 0;
}


// simple NIHT test
/*
int main()
{
	int m = 19, n = 23, k = 6, maxIt = 1000;
	double tol = 1e-6;
	MMatrix A(m, n);
	MVector x(n);
	A.initialise_normal();
	x.initialise_normal(k);
	MVector b = A * x;
	std::cout << NIHT(A, b, x, k, maxIt, tol) << "\n";
	return 0;
}
*/