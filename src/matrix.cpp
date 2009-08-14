/******************************************************************************
 *  _ _ _      ____                                                           *
 * | (_) |__  / ___|_ __ __ _ ___ ___ _ __ ___   __ _ _ __  _ __              *
 * | | | '_ \| |  _| '__/ _` / __/ __| '_ ` _ \ / _` | '_ \| '_ \             *
 * | | | |_) | |_| | | | (_| \__ \__ \ | | | | | (_| | | | | | | |            *
 * |_|_|_.__/ \____|_|  \__,_|___/___/_| |_| |_|\__,_|_| |_|_| |_|            *
 *                                                                            *
 *                                                                            *
 * The files in this directory and elsewhere which refer to this LICENCE      *
 * file are part of Grassmann, the library for the high-level management of   *
 * linear algebra structs.                                                    *
 *                                                                            *
 * Copyright (C) 2009 by BlackLight, <blacklight@autistici.org>               *
 * Web: http://0x00.ath.cx                                                    *
 *                                                                            *
 * Grassmann is free software; you can redistribute it and/or modify it under *
 * the terms of the GNU General Public License as published by the Free       *
 * Software Foundation; either version 3 or (at your option) any later        *
 * version.                                                                   *
 *                                                                            *
 * Grassmann is distributed in the hope that it will be useful, but WITHOUT   *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License      *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with Grassmann; if not, write to the Free Software Foundation, Inc.,       *
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.                     *
 *                                                                            *
 * As a special exception, if other files instantiate templates or use        *
 * macros or inline functions from these files, or you compile these          *
 * files and link them with other works to produce a work based on these      *
 * files, these files do not by themselves cause the resulting work to be     *
 * covered by the GNU General Public License. However the source code for     *
 * these files must still be made available in accordance with section (3)    *
 * of the GNU General Public License.                                         *
 *                                                                            *
 * This exception does not invalidate any other reasons why a work based on   *
 * this file might be covered by the GNU General Public License.              *
 *******************************************************************************/

#include <algorithm>
#include <sstream>
#include <cstdarg>
#include <cmath>

#include "grassmann.hpp"
#include "grassmann_exception.hpp"

using std::vector;
using std::string;
using std::stringstream;
using std::exception;

namespace grassmann  {
Matrix::Matrix()  {}

Matrix::Matrix(const Matrix &m)  {
	matrix.assign(m.matrix.begin(), m.matrix.end());
}

Matrix::Matrix(double **m, size_t r, size_t c)  {
	matrix = vector<Vector>(r);

	for (size_t i=0; i<r; i++)
		for (size_t j=0; j<c; j++)
			matrix[i].push_back(m[i][j]);
}

Matrix::Matrix(size_t r, size_t c)  {
	matrix = vector<Vector>(r);

	for (size_t i=0; i<r; i++)
		matrix[i] = vector<double>(c, 0.0);
}

size_t Matrix::rows()  {
	return matrix.size();
}

size_t Matrix::cols()  {
	if (matrix.size() == 0)
		return 0;

	return matrix[0].size();
}

string Matrix::toString()  {
	stringstream ss (stringstream::in | stringstream::out);
	vector<size_t> maxlen(cols(), 0);
	
	for (size_t i=0; i < rows(); i++)  {
		for (size_t j=0; j < cols(); j++)  {
			stringstream tmps (stringstream::in | stringstream::out);

			tmps << matrix[i][j];

			if (tmps.str().length() > maxlen[j])
				maxlen[j] = tmps.str().length();
		}
	}

	for (size_t i=0; i < rows(); i++)  {
		for (size_t j=0; j < cols(); j++)  {
			if (j == 0)
				ss << "(";

			ss << matrix[i][j];

			stringstream tmps (stringstream::in | stringstream::out);
			tmps << matrix[i][j];

			for (size_t k=0; k < maxlen[j] - tmps.str().length(); k++)
				ss << " ";

			if (j == cols()-1)  {
				ss << ")";

				if (i < rows()-1)
					ss << std::endl;
			} else
				ss << "  ";
		}
	}

	return ss.str();
}

Matrix& operator+ (Matrix a, Matrix b) throw(UnequalMatrixSizeException)  {
	if ( (a.rows() != b.rows()) || (a.cols() != b.cols()) )
		throw UnequalMatrixSizeException();
	
	Matrix *c = new Matrix(a.rows(), a.cols());

	for (size_t i=0; i < c->rows(); i++)
		for (size_t j=0; j < c->cols(); j++)
			c->matrix[i][j] = a(i,j) + b(i,j);

	return *c;
}

Matrix& operator- (Matrix a, Matrix b) throw(UnequalMatrixSizeException)  {
	if ( (a.rows() != b.rows()) || (a.cols() != b.cols()) )
		throw UnequalMatrixSizeException();
	
	Matrix *c = new Matrix(a.rows(), a.cols());

	for (size_t i=0; i < c->rows(); i++)
		for (size_t j=0; j < c->cols(); j++)
			c->matrix[i][j] = a(i,j) - b(i,j);

	return *c;
}

Matrix& operator* (Matrix a, Matrix b) throw(UnequalMatrixSizeException)  {
	Matrix *c = NULL;

	if (a.cols() != b.rows())  {
		throw UnequalMatrixSizeException();
		return *c;
	}

	c = new Matrix(a.rows(), b.cols());
	
	for (size_t i=0; i < a.rows(); i++)
		for (size_t j=0; j < b.cols(); j++)
			for (size_t k=0; k < a.cols(); k++)
				(*c).matrix[i][j] += a.matrix[i][k] * b.matrix[k][j];

	return *c;
}

Vector operator* (Matrix a, Vector v) throw(UnequalMatrixSizeException)  {
	if (a.cols() != v.size())
		throw UnequalMatrixSizeException();

	Vector res(a.rows(), 0.0);

	for (size_t i=0; i < a.rows(); i++)
		for (size_t j=0; j<v.size(); j++)
			res[i] = res[i] + (a.matrix[i][j]*v[j]);

	return res;
}

Matrix& operator* (Matrix a, double l)  {
	Matrix *c = new Matrix(a.rows(), a.cols());

	for (size_t i=0; i < a.rows(); i++)
		for (size_t j=0; j < a.cols(); j++)
			c->matrix[i][j] = l*a.matrix[i][j];

	return *c;
}

Matrix& operator* (double l, Matrix a)  {
	return a*l;
}

bool operator== (Matrix a, Matrix b)  {
	if ( (a.rows() != b.cols()) || (a.cols() != b.cols()) )
		return false;

	for (size_t i=0; i < a.rows(); i++)
		for (size_t j=0; j < a.cols(); j++)
			if (a(i,j) != b(i,j))
				return false;

	return true;
}

bool operator!= (Matrix a, Matrix b)  {
	return !(a == b);
}

double& Matrix::operator() (size_t i, size_t j) throw(InvalidMatrixIndexException)  {
	if (i >= rows())
		throw InvalidMatrixIndexException();
	
	if (j >= cols())
		throw InvalidMatrixIndexException();

	return matrix[i][j];
}

Vector& Matrix::operator[] (size_t i) throw (InvalidMatrixIndexException)  {
	if (i >= rows())
		throw InvalidMatrixIndexException();

	return matrix[i];
}

void Matrix::insertColumn (Vector v, size_t index) throw(InvalidMatrixIndexException)  {
	if (index >= cols())
		throw InvalidMatrixIndexException();

	if (v.size() != rows())
		throw InvalidMatrixIndexException();

	for (size_t i=0; i < rows(); i++)
		matrix[i][index] = (double) v[i];
}

void Matrix::insertColumn (size_t index, ...) throw(InvalidMatrixIndexException)  {
	if ((index<0) || (index>=cols()))
		throw InvalidMatrixIndexException();

	va_list col;
	va_start(col, index);

	for (size_t i=0; i < rows(); i++)  {
		double val = va_arg(col,double);
		matrix[i][index] = val;
	}

	va_end(col);
}

void Matrix::insertRow (Vector v, size_t index) throw(InvalidMatrixIndexException)  {
	if ((index<0) || (index>=rows()))
		throw InvalidMatrixIndexException();

	if (v.size() != cols())
		throw InvalidMatrixIndexException();

	for (size_t i=0; i < cols(); i++)
		matrix[index][i] = (double) v[i];
}

void Matrix::insertRow (size_t index, ...) throw(InvalidMatrixIndexException)  {
	if ((index<0) || (index>=rows()))
		throw InvalidMatrixIndexException();

	va_list row;
	va_start(row, index);

	for (size_t i=0; i < cols(); i++)  {
		double val = va_arg(row,double);
		matrix[index][i] = val;
	}

	va_end(row);
}

Matrix Matrix::validate() throw(NullMatrixDiagException) {
	size_t steps=0;
	size_t min = (rows() < cols()) ? rows() : cols();
	bool valid = true;
	Matrix m;

	if (min==1)
		return m;

	m = Matrix(*this);
	m = m.removeNullRows();
	m = m.removeNullCols();

	do  {
		valid=true;

		for (size_t i=0; i<min; i++)  {
			if (m.matrix[i][i] == 0)  {
				if (i==0)
					m.swapRows(i,i+1);
				else if (i==min-1)
					m.swapRows(i,i-1);
				else  {
					if (matrix[i+1][i] != 0)
						m.swapRows(i,i+1);
					else
						m.swapRows(i,i-1);
				}
			
				steps++;
				valid=false;
			}
		}
	} while ((!valid) && (steps < rows()*(rows()-1)));

	if (!valid)
		throw NullMatrixDiagException();

	return m;
}

Matrix Matrix::validate(size_t& steps) throw(NullMatrixDiagException)  {
	steps=0;
	size_t min = (rows() < cols()) ? rows() : cols();
	bool valid=true;
	Matrix m;

	if (min==1)
		return m;

	m = Matrix(*this);
	
	do  {
		valid=true;

		for (size_t i=0; i<min; i++)  {
			if (m.matrix[i][i] == 0)  {
				if (i==0)
					m.swapRows(i,i+1);
				else if (i==min-1)
					m.swapRows(i,i-1);
				else  {
					if (matrix[i+1][i] != 0)
						m.swapRows(i,i+1);
					else
						m.swapRows(i,i-1);
				}
			
				steps++;
				valid=false;
			}
		}
	} while ((!valid) && (steps < rows()*(rows()-1)));

	if (!valid)
		throw NullMatrixDiagException();

	return m;
}

Matrix Matrix::validate(Vector& v) throw(NullMatrixDiagException)  {
	size_t steps=0;
	size_t min = (rows() < cols()) ? rows() : cols();
	bool valid=true;
	Matrix m;

	if (min==1)
		return m;

	m = Matrix(*this);

	do  {
		valid=true;

		for (size_t i=0; i < min; i++)  {
			if (m.matrix[i][i] == 0)  {
				if (i==0)  {
					m.swapRows(i,i+1);
					
					double tmp = v[i];
					v[i] = v[i+1];
					v[i+1] = tmp;
				} else if (i==min-1) {
					m.swapRows(i,i-1);
					
					double tmp = v[i];
					v[i] = v[i-1];
					v[i-1] = tmp;
				} else {
					if (matrix[i+1][i] != 0)  {
						m.swapRows(i,i+1);
						
						double tmp = v[i];
						v[i] = v[i+1];
						v[i+1] = tmp;
					} else {
						m.swapRows(i,i-1);
						
						double tmp = v[i];
						v[i] = v[i+1];
						v[i+1] = tmp;
					}
				}
			
				steps++;
				valid=false;
			}
		}
	} while ((!valid) && (steps < rows()*(rows()-1)));

	if (!valid)
		throw NullMatrixDiagException();

	return m;
}

Matrix triang(Matrix m, size_t& steps) throw(NullMatrixDiagException)  {
	Vector v;
	Matrix a(m);

	for (size_t n=0; n < a.cols(); n++)  {
		Matrix b(a.rows(), a.cols());
		v.clear();

		if (a.matrix[n][n] == 0)  {
			a = a.validate(steps);

			if (a.cols() == 0)
				throw NullMatrixDiagException();
		}

		double pivot = a.matrix[n][n];

		for (size_t i=0; i < a.rows(); i++)  {
			if (i<n)
				v.push_back(0);
			else if (i==n)
				v.push_back(1);
			else  { 
				double tmp = -(a.matrix[i][n]);
				tmp /= pivot;
				v.push_back(tmp);
			}
		}

		b.insertColumn(v,n);
		
		for (size_t i=0; i < a.cols(); i++)  {
			if (i!=n)  {
				v.clear();

				for (size_t j=0; j < a.rows(); j++)  {
					if (i==j)
						v.push_back(1);
					else
						v.push_back(0);
				}

				b.insertColumn(v,i);
			}
		}

		Matrix c = b*a;
		
		for (size_t i=0; i < c.cols(); i++)
			for (size_t j=0; j < c.rows(); j++)
				a.matrix[i][j] = c.matrix[i][j];
	}

	return a;
}


Matrix triang(Matrix m) throw(NullMatrixDiagException)  {
	Vector v;
	Matrix a(m);

	for (size_t n=0; n < a.cols(); n++)  {
		Matrix b(a.rows(), a.cols());
		v.clear();
		double pivot = a.matrix[n][n];

		if (pivot==0)
			throw NullMatrixDiagException();

		for (size_t i=0; i < a.rows(); i++)  {
			if (i<n)
				v.push_back(0);
			else if (i==n)
				v.push_back(1);
			else  { 
				double tmp = -(a.matrix[i][n]);
				tmp /= pivot;
				v.push_back(tmp);
			}
		}

		b.insertColumn(v,n);

		for (size_t i=0; i < a.cols(); i++)  {
			if (i!=n)  {
				v.clear();

				for (size_t j=0; j < a.rows(); j++)  {
					if (i==j)
						v.push_back(1);
					else
						v.push_back(0);
				}

				b.insertColumn(v,i);
			}
		}

		Matrix c = b*a;
		
		for (size_t i=0; i < c.cols(); i++)
			for (size_t j=0; j < c.rows(); j++)
				a.matrix[i][j] = c.matrix[i][j];
	}

	return a;
}

double Matrix::det() throw(NotSquareMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();

	if (rows() == 1)
		return matrix[0][0];

	for (size_t i=0; i < rows(); i++)
		if (rowAt(i).isNull())
			return 0;
	
	for (size_t i=0; i < cols(); i++)
		if (colAt(i).isNull())
			return 0;

	Matrix a(*this);

	size_t steps=0;
	Matrix c(a.rows(), a.cols());

	try  {
		c = triang(a, steps);
	}

	catch (exception e)  {
		return 0;
	}

	if (c.rows() == 0)
		return 0;

	double d=1;

	for (size_t i=0; i < c.rows(); i++)
		d *= c.matrix[i][i];

	return ((steps%2) ? -d : d);
}

Vector solve (Matrix a, Vector coeff) throw(NotSquareMatrixException, UndefinedSystemException)  {
	Vector x,v;

	if (a.cols() != a.rows())
		throw NotSquareMatrixException();

	if (a.rank() != a.cols())
		throw UndefinedSystemException();

	x = vector<double>(coeff.size(), 0.0);

	if (a.cols() != coeff.size())
		return x;

	for (size_t n=0; n < a.cols(); n++)  {
		Matrix b(a.rows(), a.cols());
		v.clear();

		if (a(n,n) == 0)
			a = a.validate(coeff);

		double pivot = a(n,n);

		for (size_t i=0; i < a.rows(); i++)  {
			if (i<n)
				v.push_back(0);
			else if (i==n)
				v.push_back(1);
			else  {
				double tmp = -(a(i,n));
				tmp /= pivot;
				v.push_back(tmp);
			}
		}

		b.insertColumn(v,n);

		for (size_t i=0; i < a.cols(); i++)  {
			if (i!=n)  {
				v.clear();

				for (size_t j=0; j < a.rows(); j++)  {
					if (i==j)
						v.push_back(1);
					else
						v.push_back(0);
				}

				b.insertColumn(v,i);
			}
		}
		
		a = b*a;
		coeff = b*coeff;
	}

	for (size_t i = coeff.size()-1; i >= 0; i--)  {
		double tmp=0.0;
		
		for (size_t j=i+1; j<x.size(); j++)
			tmp += (a(i,j)*x[j]);

		x[i] = (coeff[i] - tmp) / a(i,i);
	}

	return x;
}

Vector Matrix::rowAt (size_t index) throw(IndexOutOfBoundsException)  {
	if (index > rows())
		throw IndexOutOfBoundsException();

	return Vector(matrix[index]);
}

Vector Matrix::colAt (size_t index) throw(IndexOutOfBoundsException)  {
	if (index > cols())
		throw IndexOutOfBoundsException();
	
	Vector r;

	for (size_t i=0; i < rows(); i++)
		r.push_back(matrix[i][index]);

	return r;
}

inline int min (int a, int b)  {
	return ((a<b) ? a : b);
}

void Matrix::swapRows (size_t a, size_t b) throw(InvalidMatrixIndexException)  {
	if ((a > rows()-1) || (b > rows()-1))
		throw InvalidMatrixIndexException();

	Vector row1 = rowAt(a);
	Vector row2 = rowAt(b);

	for (size_t i=0; i<cols(); i++)
		matrix[b][i] = row1[i];
	for (size_t i=0; i<cols(); i++)
		matrix[a][i] = row2[i];
}

void Matrix::swapCols (size_t a, size_t b) throw(InvalidMatrixIndexException)  {
	if ((a > cols()-1) || (b > cols()-1))
		throw InvalidMatrixIndexException();

	Vector row1 = colAt(a);
	Vector row2 = colAt(b);

	for (size_t i=0; i < rows(); i++)
		matrix[i][b] = row1[i];
	for (size_t i=0; i < rows(); i++)
		matrix[i][a] = row2[i];
}

Matrix Matrix::removeMultipleRows()  {
	vector<int> dup;

	for (size_t i=0; i < rows(); i++)  {
		Vector row1 = rowAt(i);

		for (size_t j=i; j < rows(); j++)  {
			if (i!=j)  {
				Vector row2 = rowAt(j);

				if (row1.multiple(row2))
					if (!binary_search(dup.begin(), dup.end(), j))
						dup.push_back(j);
			}
		}
	}

	Matrix b(rows() - dup.size(), cols());
	int pos = 0;

	for (size_t i=0; i < rows(); i++)  {
		if (!binary_search(dup.begin(), dup.end(), i))  {
			Vector row = rowAt(i);
			b.insertRow(row,pos++);
		}
	}

	return b;
}

Matrix Matrix::removeMultipleCols()  {
	vector<int> dup;

	for (size_t i=0; i < cols(); i++)  {
		Vector col1 = colAt(i);

		for (size_t j=i; j<cols(); j++)  {
			if (i!=j)  {
				Vector col2 = colAt(j);

				if (col1.multiple(col2))
					if (!binary_search(dup.begin(), dup.end(), j))
						dup.push_back(j);
			}
		}
	}

	Matrix b(rows(), cols() - dup.size());
	int pos = 0;

	for (size_t i=0; i<cols(); i++)  {
		if (!binary_search(dup.begin(), dup.end(), i))  {
			Vector col = colAt(i);
			b.insertColumn(col,pos++);
		}
	}

	return b;
}

Matrix Matrix::removeNullRows()  {
	vector<int> nulls;

	for (size_t i=0; i<rows(); i++)
		if (rowAt(i).isNull())
			nulls.push_back(i);

	Matrix b( rows()-nulls.size(), cols() );
	int pos=0;

	for (size_t i=0; i<rows(); i++)  {
		if (!binary_search(nulls.begin(), nulls.end(), i))  {
			Vector tmp = rowAt(i);
			b.insertRow(tmp,pos++);
		}
	}

	return b;
}

Matrix Matrix::removeNullCols()  {
	vector<int> nulls;

	for (size_t i=0; i<cols(); i++)
		if (colAt(i).isNull())
			nulls.push_back(i);

	Matrix b( rows(), cols()-nulls.size() );
	int pos=0;

	for (size_t i=0; i<cols(); i++)  {
		if (!binary_search(nulls.begin(), nulls.end(), i))  {
			Vector tmp = colAt(i);
			b.insertColumn(tmp,pos++);
		}
	}

	return b;
}

size_t Matrix::rank()  {
	Vector v;
	size_t r = (rows() < cols()) ? rows() : cols();

	if (r == 1)
		return 1;

	Matrix a(rows(), cols());

	for (size_t i=0; i<r; i++)
		a.insertColumn(colAt(i),i);

	try  {
		if (a.det())
			return r;
	}

	catch (NullMatrixDiagException e)  {}

	a = a.removeNullRows();
	a = a.removeNullCols();
	a = a.removeMultipleRows();
	a = a.removeMultipleCols();

	r = min(a.rows(), a.cols());

	while (r>0)  {
		Matrix b(r,r);

		try  {
			a = a.validate();
		}

		catch (exception e)  {
			r--;
			continue;
		}

		for (size_t i=0; i<r; i++)
			for (size_t j=0; j<r; j++)
				b.matrix[i][j] = a.matrix[i][j];

		try  {
			if (b.det())
				return r;
		}

		catch (exception e)  {}
	}

	return 0;
}

double Matrix::norm_inf()  {
	double t,norm=0;

	for (size_t i=0; i<rows(); i++)  {
		t=0;
		
		for (size_t j=0; j<cols(); j++)
			t += abs(matrix[i][j]);

		if (t>norm)
			norm=t;
	}

	return norm;
}

double Matrix::norm1()  {
	double t,norm=0;

	for (size_t j=0; j<cols(); j++)  {
		t=0;

		for (size_t i=0; i<rows(); i++)
			t += abs(matrix[i][j]);

		if (t>norm)
			norm=t;
	}

	return norm;
}

Matrix identityMatrix (size_t n)  {
	Matrix a(n,n);

	for (size_t i=0; i < n; i++)
		for (size_t j=0; j < n; j++)  {
			if (i == j)
				a(i,i) = 1;
			else
				a(i,j) = 0;
		}

	return a;
}

Matrix Matrix::transpose()  {
	Matrix a(cols(), rows());

	for (size_t i=0; i < cols(); i++)
		for (size_t j=0; j < rows(); j++)
			a(i,j) = matrix[j][i];

	return a;
}

Vector Matrix::eigenValues() throw(NotSquareMatrixException, SingularMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();

	if (rank() < rows())
		throw SingularMatrixException();

	if (isDiagonal())  {
		Vector v;

		for (size_t i=0; i < rows(); i++)
			v.push_back(matrix[i][i]);

		return v;
	}

	Matrix A(*this);

	for (size_t n=0; n < 100; n++)  {
		Vector x = A.colAt(0);
		double alpha = (x[0] >= 0) ? x.norm2() : - (x.norm2());

		Vector u = x;
		u[0] -= alpha;

		Vector v = u;

		for (size_t i=0; i < v.size(); i++)
			v[i] /= u.norm2();

		Matrix V(rows(), cols());

		for (size_t i=0; i < V.rows(); i++)
			for (size_t j=0; j < V.cols(); j++)
				V(i,j) = v[i]*v[j];

		Matrix Q = identityMatrix(rows()) - 2*V;
		Matrix R = Q.transpose() * A;
		A = R*Q;
	}

	Vector eigen;

	for (size_t i=0; i < rows(); i++)
		eigen.push_back(A(i,i));

	return eigen;
}

Matrix Matrix::subMatrix (size_t x1, size_t y1, size_t x2, size_t y2) throw(IndexOutOfBoundsException)  {
	if (x1 > rows() || x2 > rows() || y1 > cols() || y2 > cols())
		throw IndexOutOfBoundsException();

	Matrix a(y2-y1+1, x2-x1+1);

	for (size_t i=x1; i <= x2; i++)
		for (size_t j=y1; j <= y2; j++)
			a(i-x1, j-y1) = matrix[i][j];

	return a;
}

bool Matrix::isSingular() throw(NotSquareMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();

	return (rank() < rows()) ? true : false;
}

bool Matrix::isDiagonal() throw(NotSquareMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();

	for (size_t i=0; i < rows(); i++)
		for (size_t j=0; j < cols(); j++)
			if ( (i != j) && (matrix[i][j] != 0) )
				return false;

	return true;
}

bool Matrix::isUpperTriangular() throw(NotSquareMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();
	
	for (size_t i=0; i < rows(); i++)
		for (size_t j=0; j < cols(); j++)
			if ( (i > j) && (matrix[i][j] != 0) )
				return false;

	return true;
}

bool Matrix::isLowerTriangular() throw(NotSquareMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();
	
	for (size_t i=0; i < rows(); i++)
		for (size_t j=0; j < cols(); j++)
			if ( (i < j) && (matrix[i][j] != 0) )
				return false;

	return true;
}

Matrix Matrix::inverse() throw(NotSquareMatrixException, SingularMatrixException)  {
	if (rows() != cols())
		throw NotSquareMatrixException();

	if (isSingular())
		throw SingularMatrixException();

	Matrix A(rows(), 2*cols());

	for (size_t i=0; i < rows(); i++)
		for (size_t j=0; j < 2*cols(); j++)  {
			if (j < cols())
				A(i,j) = matrix[i][j];
			else  {
				if (i != j-cols())
					A(i,j) = 0.0;
				else
					A(i,j) = 1.0;
			}
		}

	for (size_t i=0; i < A.rows(); i++)  {
		if (A(i,i) == 0.0)  {
			for (size_t k=0; k < A.rows(); k++)  {
				if ( (i != k) && (A(k,i) != 0.0) )
					A.swapRows(i,k);
			}
		}

		double pivot = A(i,i);

		for (size_t j=0; j < A.cols(); j++)
			A(i,j) /= pivot;

		for (size_t k=0; k < A.rows(); k++)  {
			double tmp = A(k,i);

			if (i != k)  {
				for (size_t j=0; j < A.cols(); j++)
					A(k,j) -= (tmp*A(i,j));
			}
		}
	}

	size_t r = rows();
	size_t c = cols();
	return A.subMatrix(0, c, r-1, 2*c - 1);
}
}

