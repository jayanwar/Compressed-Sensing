#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>

// returns a standard normal value (i.e. mean 0 and variance 1)
double rand_normal()
{
	static const double pi = 3.141592653589793238;
	double u = 0, v;
	while (u == 0) // loop to ensure u non-zero for log
		u = rand() / static_cast<double>(RAND_MAX);
	v = rand() / static_cast<double>(RAND_MAX);
	return sqrt(-2.0*log(u))*cos(2.0*pi*v);
}


// Class that represents a mathematical vector
class MVector
{
public:
	/* constructors */
	MVector() : v(0) {}

	// initialise empty vector of size n 
	explicit MVector(int n) : v(n) {}

	// initialise vector of size n with values x
	MVector(int n, double x) : v(n, x) {}

	/* sorting functions */

	// returns no. of elements in vector
	unsigned size() const { return v.size(); }

	// swaps two elements in a vector
	void swap(int i, int j)
	{
		double temp = v[i];
		v[i] = v[j];
		v[j] = temp;
	}

	// compares order of two elements in a vector
	bool cmp(int i, int j)
	{
		if (v[i] < v[j]) return true;
		else return false;
	}


	/* initialise random vector */

	// for an n-vector, give each place a random value in the range [xmin,xmax]
	void initialise_random(double xmin, double xmax)
	{
		size_t s = v.size();
		for (size_t i = 0; i < s; i++)
		{
			v[i] = xmin + (xmax - xmin)*rand() / static_cast<double>(RAND_MAX);
		}
	}
	
	// initialise random k-sparse vector with (standard) normally distributed 
	// values
	void initialise_normal(int k)
	{
		int i = 0;
		while (i < k)
		{
			int j = rand() % v.size(); // produce random index for vector v
			if (v[j] == 0.0)
			{
				v[j] = rand_normal(); // if index empty (i.e. = 0) place random
				i++;				  // value here
			}
		}
	}
	
	// returns the dot product of vectors v and w
	double dot(MVector w) const
	{
		double dp = 0.0;
		for (unsigned i = 0; i < v.size(); i++) dp += v[i] * w[i];
		return dp;
	}

	// implements median of three rule, places the median of the first, last
	// middle elements of the subsequence of the vector at the start
	void MedOfThree(int start, int end)
	{
		int mid = start + (end - start) / 2;
		if (this->cmp(end, start)) this->swap(end, start);
		if (this->cmp(mid, start)) this->swap(mid, start);
		if (this->cmp(end, mid)) this->swap(end, mid);
		this->swap(start, mid);
	}

	// implements the quick sorting step - places the pivot in the correct places,
	// all elements less than pivot before it and all elements more than pivot 
	// above it
	int CorrectPivotPlacement(int start, int end)
	{
		this->MedOfThree(start, end); // choose pivot and place at start
		int p_goes_here = start + 1; // index of pivot
		for (int i = start + 1; i <= end; i++)
		{
			// if v[i] < pivot, increase index of pivot by 1
			if (this->cmp(i, start))
			{
				this->swap(p_goes_here, i);
				p_goes_here++;
			}
		}
		p_goes_here--;
		this->swap(start, p_goes_here); // place pivot in correct position
		return p_goes_here; // return index of pivot
	}

	// recursive call for quick sort
	void quick_recursive(int start, int end)
	{
		if (start == end) { ; } // do nothing if vector is empty or 1 dim.
		else
		{
			int pivot_index = this->CorrectPivotPlacement(start, end);
			// do quicksort on either side of pivot if possible
			if (pivot_index != start) this->quick_recursive(start, pivot_index - 1);
			if (pivot_index != end) this->quick_recursive(pivot_index + 1, end);
		}
	}

	// wrapper for quick recursive
	void quick() { this->quick_recursive(0, v.size() - 1); }

	// makes all values in vector absolute
	void absall() { for (int i = 0; i < v.size(); i++) v[i] = abs(v[i]); }

	// sets all but k absolutely largest elements to 0 in the vector
	void threshold(int k)
	{
		// copy vector v to w and quicksort to find kth largest element
		MVector w(v.size());
		w = *this;
		w.absall();
		w.quick();
		// set all elements absolutely smaller than kth largest element of w to 0
		double kthmax = w[v.size() - k];
		for (int i = 0; i < v.size(); i++) { if (abs(v[i]) < kthmax) v[i] = 0.0; }
	}

	// returns the Euclidean norm of a vector v
	double norm()
	{
		double modv = 0.0;
		for (unsigned i = 0; i < v.size(); i++) modv += v[i] * v[i];
		modv = std::sqrt(modv);
		return modv;
	}

	// scales the vector v with respect to scalar a
	MVector scale(double a)
	{
		MVector x(v.size());
		for (int i = 0; i < v.size(); i++) x[i] = a * v[i];
		return x;
	}

	// checks for convergence to the vector in the SDLS and NIHT algorithms
	// (within an error of 1.0e-5 for each index)
	bool has_converged()
	{
		for (int i = 0; i < v.size(); i++) if (abs(v[i]) > 1.0e-5) return false;
		return true;
	}
	
	/* overloaded operator functions */
	// access element (lvalue)
	double &operator[](int index)
	{
		assert(index >= 0 && index < v.size());
		return v[index];
	}

	// access element (rvalue)
	double operator[](int index) const
	{
		assert(index >= 0 && index < v.size());
		return v[index];
	}

	// set one vector to take the values of another vector
	void operator=(MVector const w)
	{
		assert(v.size() == w.size());
		for (int i = 0; i < v.size(); i++) v[i] = w[i];
	}

	// coordinate-wise vector addition
	MVector operator+(MVector w) const
	{
		assert(v.size() == w.size());
		MVector x(v.size());
		for (int i = 0; i < v.size(); i++) x[i] = v[i] + w[i];
		return x;
	}

	// coordinate-wise vector subtraction
	MVector operator-(MVector w) const
	{
		assert(v.size() == w.size());
		MVector x(v.size());
		for (int i = 0; i < v.size(); i++) x[i] = v[i] - w[i];
		return x;
	}

	// print n-vector in the format v0,v1,....,vn,
	friend std::ostream& operator<<(std::ostream& os, const MVector& v);

private:
	std::vector<double> v;
};

// Class that represents a mathematical matrix
class MMatrix
{
public:
	// Constructors and destructors (see also vector class)
	explicit MMatrix() : mRows(0), nCols(0) {}
	MMatrix(int n, int m) : mRows(n), nCols(m), mtrx(n, std::vector<double>(m)) {}
	// Include other methods here

	// initialise a random matrix of values in the range [xmin,xmax]
	void initialise_random(double xmin, double xmax)
	{
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < nCols; j++)
				mtrx[i][j] = xmin + (xmax - xmin)*rand()
				/ static_cast<double>(RAND_MAX);
	}

	// intialise a random matrix of normally distributed values with mean zero and
	// variance 1/m (for an m x n matrix)
	void initialise_normal()
	{
		double sd = 1.0 / sqrt(mRows);
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < nCols; j++)
				mtrx[i][j] = rand_normal() * sd; // place random value here
	}


	// matrix-matrix multiplication
	MMatrix operator*(MMatrix B)
	{
		assert(nCols == B.mRows); // check if dimensions are sufficient for 
								  // matrix- matrix multiplication
		MMatrix C(mRows, B.nCols);
		for (unsigned i = 0; i < mRows; i++)
			for (unsigned j = 0; j < B.nCols; j++)
			{
				double Cij = 0.0;
				for (unsigned k = 0; k < nCols; k++)
					Cij += mtrx[i][k] * B.mtrx[k][j];
				C.mtrx[i][j] = Cij;
			}
		return C;
	}
	
	// matrix-vector multiplication
	MVector operator*(MVector v) const
	{
		assert(v.size() == nCols); // check if dimensions are sufficient for 
								   // matrix-vector multiplication
		MVector w(mRows);
		for (unsigned i = 0; i < mRows; i++)
		{
			double wi = 0.0;
			for (unsigned k = 0; k < v.size(); k++)
			{
				wi += mtrx[i][k] * v[k];
			}
			w[i] = wi;
		}
		return w;
	}

	// return the transpose of a matrix
	MMatrix transpose() const
	{
		MMatrix B(nCols, mRows);
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < nCols; j++)
				B.mtrx[j][i] = mtrx[i][j];
		return B;
	}

	// scale a matrix with respect to a scalar a
	MMatrix scale(double a) const
	{
		MMatrix B(mRows, nCols);
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < nCols; j++)
				B.mtrx[i][j] = a * mtrx[i][j];
		return B;
	}

	// set the value in the ith row and jth column to the double a
	void set(int i, int j, double a)
	{
		assert(i >= 0 && j >= 0 && i < mRows && j < nCols);
		mtrx[i][j] = a;
	}

	// print matrix to file
	void print()
	{
		std::ofstream matrix;
		matrix.open("randmatrix.csv");
		for (int i = 0; i < mRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				matrix << mtrx[i][j] << ",";
			}
			matrix << "\n";
		}
	}

	friend std::ostream& operator<<(std::ostream& os, const MMatrix& A);

private:
	unsigned int mRows, nCols; // dimensions of matrix
	std::vector<std::vector<double>> mtrx; // matrix a vector of double vectors
};

// print MVector in the format v0 ,v1 ,.... , vn ,
std::ostream & operator <<(std::ostream & os, const MVector & v)
{
	for (unsigned i = 0; i < v.size(); i++) { os << v[i] << ","; }
	return os;
}

// overload << operator to print matrix
std::ostream & operator <<(std::ostream & os, const MMatrix & A)
{
	for (unsigned i = 0; i < A.mtrx.size(); i++)
	{
		os << "|";
		for (unsigned j = 0; j < A.mtrx[i].size(); j++) 
			os << A.mtrx[i][j] << ",";
		os << "|\n";
	}
	return os;
}


#endif