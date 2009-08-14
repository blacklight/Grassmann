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

/**
 * @file grassmann.hpp
 * @version 0.2
 * @author BlackLight < blacklight@autistici.org >
 */

#ifndef __cplusplus
#error  "This is a C++ library, you know, so you would need a C++ compiler to use it"
#else

#ifndef _GRASSMANN_H
#define _GRASSMANN_H

#include <string>
#include <vector>

#include "grassmann_exception.hpp"

/**
 * @namespace grassmann
 * Main namespace of the library
 */
namespace grassmann {

/**
 * @class Vector
 * @brief Class for operations on vectors
 * @author BlackLight
 */

	class Vector {
		std::vector <double> v;

public:
	/**
	 * @brief Default constructor. It simply does nothing.
	 */
		Vector();

	/**
	 * @brief Constructor, given a vector<double> object
	 * @param v Input vector<double>
	 */
		 Vector(std::vector <double> v);

	/**
	 * @brief Constructor, given a double[] array
	 * @param v double[] array
	 * @param size Size
	 */
		 Vector(double *v, size_t size);

	/**
	 * @brief Constructor, given a list of double values
	 * @param size Number of elements to be included
	 */
		 Vector(size_t size, ...);

	/**
	 * @brief Constructor, returns a vector from a list of
	 * terms given as string, with a delimitator between
	 * each of them (default delimitator: ','). i.e.,
	 * "2,3,4" with delim = ',' will return a vector containing
	 * terms (2 3 4)
	 * @param str String containing the terms
	 * @param delim Delimitator between each term
	 */
		 Vector(std::string, char delim = ',');

	/**
	 * @brief Gets the 2-norm of a double vector of size n
	 * @return 2-norm (Euclidean distance)
	 */
		double norm2();

	/**
	 * @brief Synonim for norm2() function (modulus, 2-norm, Euclidean distance)
	 */
		double modulus();

	/**
	 * @brief Sum of two vectors
	 * @param a First vector
	 * @param b Second vector
	 * @return The vector that represents the sum between a and b
	 * @exception UnequalVectorSizeException If the size of a and b is not the same
	 */
		friend Vector operator+(Vector& a, Vector& b) throw(UnequalVectorSizeException);

	/**
	 * @brief Difference between two vectors
	 * @param a First vector
	 * @param b Second vector
	 * @return The vector that represents the difference between a and b
	 * @exception UnequalVectorSizeException If the size of a and b is not the same
	 */
		friend Vector operator-(Vector& a, Vector& b) throw(UnequalVectorSizeException);

	/**
	 * @brief Scalar product between two vectors
	 * @param a First vector
	 * @param b Second vector
	 * @return Scalar product
	 * @exception UnequalVectorSizeException If the size of a and b is not the same
	 */
		friend double operator*(Vector& a, Vector& b) throw(UnequalVectorSizeException);

	/**
	 * @brief Product between a vector and a scalar
	 * @param l The scalar
	 * @param v The vector
	 * @return The product l*v
	 */
		friend Vector operator*(double l, Vector v);
		friend Vector operator*(Vector v, double l);

	/**
	 * @brief Operator the vectorial product z = v x w between two vectors, v and w
	 * @param v First vector
	 * @param w Second vector
	 * @return A vector representing the vectorial product v x w
	 * @throw InvalidVectorSizeException When the size of the vectors is != 3
	 */
		friend Vector& operator% (Vector& v, Vector& w) throw(UnequalVectorSizeException, InvalidVectorSizeException);

	/**
	 * @brief Check if two vectors are equal, i.e. if their size is the same and all of their elements are equal
	 * @param a First vector
	 * @param b Second vector
	 * @return true if they are equal, false otherwise
	 */
		friend bool operator==(Vector a, Vector b);

	/**
	 * @brief Check if two vectors are different, i.e. if their size is not the same or their elements are different
	 * @param a First vector
	 * @param b Second vector
	 * @return true if they are different, false otherwise
	 */
		friend bool operator!=(Vector a, Vector b);

	/**
	 * @brief Construct the vector from a string
	 * i.e.: Vector v(3); v = "2,3,4";
	 * @param s String containing the elements of the vector (warning: you must use ',' as delimiter
	 * between the terms)
	 */
		void operator= (std::string s);

	/**
	 * @brief Gets the i-th member of a vector
	 * @param i Index to get
	 * @return The value to be read, if exists
	 * @exception IndexOutOfBoundsException If trying to access an index outside of the vector
	 */
		double& operator[] (size_t i) throw(IndexOutOfBoundsException);

		Vector& operator<< (double element);

	/**
	 * @brief Gets the size of the vector
	 */
		size_t size();

	/**
	 * @brief Gets the maximum element of a vector
	 */
		double max();

	/**
	 * @brief Push a new element into the array
	 * @param el Element to be pushed
	 */
		void push_back(double el);

	/**
	 * @brief Checks if the vector contains only null elements
	 * @return true if is empty, false otherwise
	 */
		bool isNull();

	/**
	 * @brief Checks if a vector is multiple of another (e.g. linear dependent)
	 * @param v1 Vector to be compared
	 * @return true if the vector are linear dependent, false otherwise
	 */
		bool multiple(Vector v1);

	/**
	 * @brief Clears the vector
	 */
		void clear();

	/**
	 * @brief Return the vector's content as a string
	 */
		 std::string toString();
	};

/**
 * @class Matrix
 * @brief Class for operations on matrices
 * @author BlackLight
 */

	class Matrix {
		std::vector <Vector> matrix;

	/**
	 * @brief Swaps two rows of the matrix
	 * @param a Index of the first row
	 * @param b Index of the second row
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void swapRows(size_t a, size_t b) throw(InvalidMatrixIndexException);

	/**
	 * @brief Swaps two columns of the matrix
	 * @param a Index of the first column
	 * @param b Index of the second column
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void swapCols(size_t a, size_t b) throw(InvalidMatrixIndexException);

	/**
	 * @brief Removes the null rows from a matrix
	 * @return Reference to the matrix so generated
	 */
		Matrix removeNullRows();

	/** 
	 * @brief Removes the null columns from a matrix
	 * @return Reference to the matrix so generated
	 */
		Matrix removeNullCols();

	/**
	 * @brief Tries to remove the rows in the matrix that are linear dependent from other rows
	 * @return Reference to the matrix so generated
	 */
		Matrix removeMultipleRows();

	/**
	 * @brief Tries to remove the columns in the matrix that are linear dependent from other rows
	 * @return Reference to the matrix so generated
	 */
		Matrix removeMultipleCols();

public:
	/**
	 * @brief Empty constructor for the class
	 */
		Matrix();

	/**
	 * @brief A constructor for the class
	 * @param m A Matrix object
	 */
		 Matrix(const Matrix & m);

	/**
	 * @brief A constructor for the class
	 * @param m A double double-array
	 * @param r Number of rows
	 * @param c Number of cols
	 */
		 Matrix(double **m, size_t r, size_t c);

	/**
	 * @brief A constructor for the class
	 * @param r Number of rows
	 * @param c Number of cols
	 */
		 Matrix(size_t r, size_t c);

	/**
	 * @brief This gets the number of rows of a matrix
	 */
		size_t rows();

	/**
	 * @brief This gets the number of cols of a matrix
	 */
		size_t cols();

	/**
	 * @brief Get the content of the row at index 'index'
	 * @param index Index to be accessed
	 * @return A Vector object representing that row
	 * @exception IndexOutOfBoundsException If trying to access an index outside of the matrix
	 */
		Vector rowAt(size_t index) throw(IndexOutOfBoundsException);

	/**
	 * @brief Get the content of the column at index 'index'
	 * @param index Index to be accessed
	 * @return A Vector object representing that column
	 * @exception IndexOutOfBoundsException If trying to access an index outside of the matrix
	 */
		Vector colAt(size_t index) throw(IndexOutOfBoundsException);

	/**
	 * @brief Insert a new vector in the matrix as column
	 * @param v Vector to be inserted
	 * @param i Index in which the vector should be inserted
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void insertColumn(Vector v, size_t i) throw(InvalidMatrixIndexException);

	/**
	 * @brief Insert a new column in the matrix
	 * @param index Index of the column to be inserted
	 * @param ... Values to be copied
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void insertColumn(size_t index, ...) throw(InvalidMatrixIndexException);

	/**
	 * @brief Insert a new vector in the matrix as row
	 * @param v Vector to be inserted
	 * @param i Index in which the vector should be inserted
	 * @return InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void insertRow(Vector v, size_t i) throw(InvalidMatrixIndexException);

	/**
	 * @brief Insert a new row in the matrix
	 * @param index Index of the column to be inserted
	 * @param ... Values to be copied
	 * @return InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		void insertRow(size_t index, ...) throw(InvalidMatrixIndexException);

	/**
	 * @brief Return the content of a matrix as a string
	 */
		 std::string toString();

	/**
	 * @brief This calculates the determinant of a square matrix
	 * @return The determinant of the matrix
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 */
		double det() throw(NotSquareMatrixException);

	/**
	 * @brief This calculates the rank of a matrix
	 */
		size_t rank();

	/**
	 * @brief Get a submatrix from the matrix, from (x1,y1) to (x2,y2)
	 * @exception IndexOutOfBoundsException If trying to access an index outside of the matrix
	 */
		Matrix subMatrix(size_t x1, size_t y1, size_t x2,
				 size_t y2) throw(IndexOutOfBoundsException);

	/**
	 * @brief Compute an approximation of the eigenValues of a square matrix, if they exist, using the iterative QR algorithm
	 * @return A vector containing an approximation of the eigenValues
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 * @exception SingularMatrixException If trying to compute the eigenvalues of a singular matrix
	 */
		Vector eigenValues() throw(NotSquareMatrixException, SingularMatrixException);

	/**
	 * @brief Compute the inverse matrix of the matrix
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 * @exception SingularMatrixException If trying to compute the eigenvalues of a singular matrix
	 */
		Matrix inverse() throw(NotSquareMatrixException, SingularMatrixException);

	/**
	 * @brief Return the transposed matrix
	 */
		Matrix transpose();

	/**
	 * @brief Re-definition of '+' operator for the sum of two matrixes
	 * @param a First matrix
	 * @param b Second matrix
	 * @exception UnequalMatrixSizeException If trying to operate on two matrices whose size is not the same
	 */
		friend Matrix & operator+(Matrix a, Matrix b) throw(UnequalMatrixSizeException);

	/**
	 * @brief Re-definition of '-' operator for subtrating two matrixes
	 * @param a First matrix
	 * @param b Second matrix
	 * @exception UnequalMatrixSizeException If trying to operate on two matrices whose size is not the same
	 */
		friend Matrix & operator-(Matrix a, Matrix b) throw(UnequalMatrixSizeException);

	/**
	 * @brief Re-definition of '*' operator for the product between two matrixes
	 * @param a First matrix
	 * @param b Second matrix
	 * @return A reference to a matrix representing the product between a and b
	 * @exception UnequalMatrixSizeException If trying to operate on two matrices whose size is not the same
	 */
		friend Matrix & operator*(Matrix a, Matrix b) throw(UnequalMatrixSizeException);

	/**
	 * @brief Re-definition of '*' operator for the product between a matrix and a vector
	 * @param a Matrix
	 * @param v Vector
	 * @return A vector representing the product between a and v
	 * @exception UnequalMatrixSizeException If trying to operate on two matrices whose size is not the same
	 */
		friend Vector operator*(Matrix a, Vector v) throw(UnequalMatrixSizeException);

	/**
	 * @brief Re-definition of '*' operator for the product between a matrix and a scalar
	 * @param a Matrix
	 * @param l Scalar
	 * @return A matrix representing the product between a and l
	 */
		friend Matrix & operator*(Matrix a, double l);
		friend Matrix & operator*(double l, Matrix a);

	/**
	 * @brief Check if two matrices are equal (i.e. the size is the same and all the elements are equal)
	 * @param a First matrix
	 * @param b Second matrix
	 * @return true if they are equal, false otherwise
	 */
		friend bool operator==(Matrix a, Matrix b);

	/**
	 * @brief Check if two matrices are not equal (i.e. the size is not the same or the elements are not equal)
	 * @param a First matrix
	 * @param b Second matrix
	 * @return true if they are not equal, false otherwise
	 */
		friend bool operator!=(Matrix a, Matrix b);

	/**
	 * @brief Re-definition of '()' operator to access the elements of a matrix. To access
	 *  the elements of a matrix you will use parenthesis. e.g.:<pre>
	 *  Matrix A(2,2);
	 *  A(0,0) = 1;
	 *  A(0,1) = 0;
	 *  A(1,0) = 0;
	 *  A(1,1) = 1;
	 *
	 *  for (size_t i=0; i < A.rows(); i++)  {
	 *    for (size_t j=0; j < A.cols(); j++)
	 *      cout << A(i,j) << " ";
	 *    cout << endl;
	 *  }
	 *
	 * @param i Index of the row
	 * @param j Index of the column
	 * @return A reference to the element at pos (i,j)
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		double& operator() (size_t i, size_t j) throw(InvalidMatrixIndexException);

	/**
	 * @brief Return the i-th row of the matrix
	 * @param i Index to be accessed
	 * @exception InvalidMatrixIndexException If trying to access an index outside of the matrix
	 */
		Vector& operator[] (size_t i) throw(InvalidMatrixIndexException);

	/**
	 * @brief Returns a triangular matrix associated to our matrix through Gauss' method. To be done if and only if
	 *        the matrix is not singular, else an exception will be raised
	 * @param m Matrix to triangularize
	 * @return Reference to the triangularized matrix
	 * @exception NullMatrixDiagException If the matrix has null elements along the diagonal and they can't be reduced
	 */
		friend Matrix triang(Matrix m) throw(NullMatrixDiagException);
		friend Matrix triang(Matrix m, size_t &steps) throw(NullMatrixDiagException);

	/**
	 * @brief Tries to validate a matrix (i.e. removing eventual null rows or column or linear dependent rows or columns)
	 * @return Reference to the validated matrix
	 * @exception NullMatrixDiagException If the matrix has null elements along the diagonal and they can't be reduced
	 */
		Matrix validate() throw(NullMatrixDiagException);

	/**
	 * @brief Tries to validate a matrix (i.e. removing eventual null rows or column or linear dependent rows or columns)
	 * @param steps Reference to an integer that will contain the number of iterations made to get a valid matrix
	 * @return Reference to the validated matrix
	 * @exception NullMatrixDiagException If the matrix has null elements along the diagonal and they can't be reduced
	 */
		Matrix validate(size_t &steps) throw(NullMatrixDiagException);

	/**
	 * @brief Tries to validate a matrix (i.e. removing eventual null rows or column or linear dependent rows or columns)
	 *        to be used to solve a linear system, so all the swaps applied to the matrix are also applied to the vector
	 *        of coefficients of the system
	 * @param v Reference to the vector of system's coefficients
	 * @return Reference to the validated matrix
	 * @exception NullMatrixDiagException If the matrix has null elements along the diagonal and they can't be reduced
	 */
		Matrix validate(Vector & v) throw(NullMatrixDiagException);

	/**
	 * @brief Gets norm one of a matrix
	 */
		double norm1();

	/**
	 * @brief Gets norm infinity of a matrix
	 */
		double norm_inf();

	/**
	 * @brief Check if the matrix is singular
	 * @return true if it is singular, false otherwise
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 */
		bool isSingular() throw(NotSquareMatrixException);

	/**
	 * @brief Check if the matrix is diagonal
	 * @return true if it is diagonal, false otherwise
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 */
		bool isDiagonal() throw(NotSquareMatrixException);

	/**
	 * @brief Check if the matrix is upper triangular
	 * @return true if it is upper triangular, false otherwise
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 */
		bool isUpperTriangular() throw(NotSquareMatrixException);

	/**
	 * @brief Check if the matrix is lower triangular
	 * @return true if it is lower triangular, false otherwise
	 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
	 */
		bool isLowerTriangular() throw(NotSquareMatrixException);
	};

/**
 * @brief Solves a linear system, for which is given the matrix a of variables' coefficients and the vector coeff
 *        of numeric coefficients. To be done if and only if the matrix is not singular and the vector does not contain
 *        null values. For example, let's consider the linear system
 *
 *<pre>		x1 + x2 + x3 = 6
 * 		2*x1 + x2 - x3 = 1
 * 		x1 + x2 + 2*x3 = 9</pre>
 *
 * 		Then, A will be
 * 		
 *<pre> 		(1  1  1)
 * 		(2  1 -1)
 * 		(1  1  2)</pre>
 *
 * 		and coeff will be (6,1,9)
 * 		A vector containing the solution (1,2,3) will be returned
 *
 * @param coeff Vector of numeric coefficients
 * @param a Matrix containing the coefficients of the variables
 * @return A vector containing, if the system is possible and the solution is unique, the solution for the system
 * @exception NotSquareMatrixException If trying to compute the determinant of a non-square matrix
 * @exception UndefinedSystemException If the system is undefined
 */
	Vector solve(Matrix a, Vector coeff) throw(NotSquareMatrixException, UndefinedSystemException);

/**
 * @brief Return the identity matrix (i.e. with all '1' along the main diagonal) by size n
 * @param n Size of the matrix
 * @return Identity matrix by size n
 */
	Matrix identityMatrix(size_t n);

/**
 * @brief Return the null vector by size n
 * @param n Size of the vector
 * @return Null vector by size n
 */
	Vector nullVector(size_t n);

/**
 * @brief From n points on a cartesian plan, given as an array of Vector(s), compute the
 * polynomial that interpolates (i.e. 'touches') all the points, and its
 * value in 'x'. e.g. the polynomial that interpolates (-2,4), (0,0) and (1,1) is p(x) = x^2,
 * so if I consider x=3 I get p(3) = 9. To do this, I'll construct a Vector[] array containing
 * the points first:
 *
 * <pre>vector<Vector> v(3);
 * v[0] = "-2,4";
 * v[1] = "0,0";
 * v[2] = "1,1";
 *
 * double x = 3.0;
 * double px = polynomialInterpolation(x,v);
 * // x = 9</pre>
 *
 * @param x Value in which I'm going to compute the interpolation polynomial
 * @param points vector of Vector(s) containing the coordinates of the points to be interpolated
 * @return The value of the interpolation polynomial in the point
 * @throw CoincidentPointsException Exception raised when two or more points in 'points' vector
 * have the same abscissa
 */
	double polynomialInterpolation (double x, std::vector<Vector> points) throw(CoincidentPointsException);

/**
 * @brief From n points on a cartesian plan, given as an array of Vector(s), compute the
 * value in 'x' variable from the linear interpolation of the points
 * @param x Value in which I'm going to compute the interpolation polynomial
 * @param points vector of Vector(s) containing the coordinates of the points to be interpolated
 * @return The value of the linear interpolation in the point
 * @throw ValueOutOfRangeException Exception raised when trying to get the linear interpolation in a
 * 'x' value smaller than the smallest value in the range of points or greater than the greatest
 * value in the range
 * @throw CoincidentPointsException Exception raised when two or more points in 'points' vector
 * have the same abscissa
 */
	double linearInterpolation (double x, std::vector<Vector> points)
		throw(ValueOutOfRangeException, CoincidentPointsException);
}

#endif
#endif

