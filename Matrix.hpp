/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* template class Matrix
*
* Comments: optimized from original vector of vectors
*
* $Id$
*
***********************************************************************/

#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include <iostream>
#include <vector>
#include "RegisterSVNId.hpp"

namespace {
RegisterSVNId Matrix_hpp_id("$Id$");
}

template <class T> class Matrix
{
    private:

	// Data Methods
	T* mdata;
	size_t mrows; 
	size_t mcols;
	size_t msize;

    public:

	// Constructors
	Matrix() :
	    mdata(NULL), mrows(0), mcols(0), msize(0)
	{ }

	Matrix(const Matrix& src) :
	    mdata(new T[src.mrows*src.mcols]),
	    mrows(src.mrows), mcols(src.mcols), msize(src.msize)
	{
	    std::copy(src.mdata, src.mdata + src.msize, mdata);
	}

	Matrix(size_t m, size_t n) :
	    mrows(m), mcols(n), msize(mrows*mcols)
	{
	    mdata = new T[mrows*mcols];
	    std::fill(mdata, mdata + msize, T(0));
	}

	// Destructor
	~Matrix()
	{
	    delete[] mdata;
	}
	// Methods
	Matrix& operator=(const Matrix& src)
	{
	    if (this == &src)	return *this;
	    if (msize != src.msize)
	    {
		delete[] mdata;
		mdata = new T[src.msize];
	    }
	    std::copy(src.mdata, src.mdata + src.msize, mdata);
	    mrows = src.mrows;
	    mcols = src.mcols;
	    msize = src.msize;
	    return *this;
	}

	Matrix& operator=(T value)
	{
	    std::fill_n(mdata, msize, value);
	    return *this;
	}

	void resize(size_t m, size_t n, T value=T(0))
	{
	    if (m == mrows && n == mcols)   return;
	    T* resized = new T[m*n];
	    std::fill(resized, resized + m*n, value);
	    for (   size_t i = 0, offset = 0;
		    i != std::min(m, mrows); ++i, offset += n )
	    {
		for (size_t j = 0; j != std::min(n, mcols); ++j)
		{
		    resized[offset+j] = (*this)(i, j);
		}
	    }
	    delete[] mdata;
	    mdata = resized;
	    mrows = m;
	    mcols = n;
	    msize = m*n;
	}

	void clear()
	{
	    delete[] mdata;
	    mdata = NULL;
	    mrows = mcols = msize = 0;
	}

	std::vector<T> rowVector(size_t i)
	{
	    return std::vector<T>(mdata + mcols*i, mdata + mcols*(i+1));
	}

	std::vector<T> columnVector(size_t j)
	{
	    std::vector<T> column(mrows);
	    typename std::vector<T>::iterator vii;
	    vii = column.begin();
	    T* pt = mdata + j;
	    for (; vii != column.end(); ++vii, pt += mcols)   *vii = *pt;
	    return column;
	}

	inline size_t rows()
	{
	    return mrows;
	}

	inline size_t columns()
	{
	    return mcols;
	}

	inline T* operator[](size_t i)
	{
	    return mdata + mcols*i;
	}

	inline T& operator()(size_t i, size_t j)
	{
	    return *(mdata + i*mcols + j);
	}

	Matrix<T> transposed() 
	{ 
	    Matrix<T> mxt(mcols, mrows);
	    T* pt = mxt.mdata;
	    for (size_t j = 0; j != mcols; ++j)
	    {
		for (size_t i = 0; i != mrows; ++i, ++pt)
		{
		    *pt = (*this)(i, j);
		}
	    }
	    return mxt;
	}

	Matrix<T> operator*(Matrix<T>& B)
	{
	    Matrix<T>& A = *this;
	    Matrix<T> C(A.mrows, B.mcols);
	    if (A.mcols != B.mrows)
	    {
		const char* emsg = "Inconsistent matrix multiplication";
		throw emsg;
	    }
	    for (size_t i = 0; i != A.mrows; ++i)
	    {
		for (size_t j = 0; j != B.mcols; ++j)
		{
		    T& Cij = C(i,j);
		    Cij = 0;
		    for (size_t k = 0; k != A.mcols; ++k)
		    {
			Cij += A(i,k)*B(k,j);
		    }
		}
	    }
	    return C;
	}

	template <class T1>
	    friend std::ostream& operator<<(std::ostream &out, Matrix<T1> &mx);

	void print() 
	{
	    using namespace std;
	    for (size_t i = 0; i != mrows; ++i)
	    {
		cout << i << ": ";
		for (size_t j = 0; j != mcols; ++j)
		{
		    cout << (*this)(i,j) << " ";
		}
		cout << endl;
	    }
	}
};

template <class T>
std::ostream& operator<<(std::ostream &out, Matrix<T> &mx)
{
    using namespace std;
    out << "( ";
    for (size_t i = 0; i != mx.mrows; ++i)
    {
	out << "(";
	for (size_t j = 0; j != mx.mcols; ++j)
	{
	    cout << mx(i,j);
	    if (j != (mx.mcols-1))  cout << ", "; 
	    else		    cout << ")";
	}
	if (i != (mx.mrows-1))	    cout << ", ";
	else			    cout << " )\n";
    }
    return out;
}

#endif	// MATRIX_HPP_INCLUDED
