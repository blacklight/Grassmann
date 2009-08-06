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
 * @file grassmann_exception.h
 * @author BlackLight < blacklight@autistici.org >
 */

#ifndef __cplusplus
#error  "This is a C++ library, you know, so you would need a C++ compiler to use it"
#else

#ifndef _GRASSMANN_EXCEPTION_H
#define _GRASSMANN_EXCEPTION_H

#include <exception>

namespace grassmann {

/**
 * @class UnequalMatrixSizeException
 * @brief Exception raised when trying to sum or doing operations between matrix having incompatible sizes
 */
	class UnequalMatrixSizeException : public std::exception {
		public:
			UnequalMatrixSizeException() {}

			const char* what() const throw()  {
				return "Invalid operation between matrixes with incompatible sizes";
			}
	};

/**
 * @class InvalidMatrixIndexException
 * @brief Exception raised when trying to access to a member of a matrix on an invalid position
 */
	class InvalidMatrixIndexException : public std::exception {
		public:
			InvalidMatrixIndexException() {}
			
			const char* what() const throw() {
				return "Trying access to an invalid index for the matrix";
			}
	};

/**
 * @class NullMatrixDiagException
 * @brief Exception raised when trying to triangularize or compute the determinant of a matrix that may be singular
 */
	class NullMatrixDiagException : public std::exception {
		public:
			NullMatrixDiagException() {}
			
			const char* what() const throw() {
				return "Singular matrix";
			}
	};

/**
 * @class NotSquareMatrixException
 * @brief Exception raised when trying to compute the determinant or solving a linear system with a matrix that is not square
 */
	class NotSquareMatrixException : public std::exception {
		public:
			NotSquareMatrixException() {}
			
			const char* what() const throw() {
				return "The given matrix is not square";
			}
	};

/**
 * @class UndefinedSystemException
 * @brief Exception raised when trying to solve an undefined linear system
 */
	class UndefinedSystemException : public std::exception {
		public:
			UndefinedSystemException() {}
			const char* what() const throw() {
				return "The given linear system is undefined";
			}
	};

/**
 * @class UnequalVectorSizeException
 * @brief Exception raised when the size of two vector is not the same
 */

	class UnequalVectorSizeException : public std::exception {
		public:
			UnequalVectorSizeException() {}
			
			const char* what() const throw() {
				return "The size of the two vectors is not the same";
			}
	};

/**
 * @class IndexOutOfBoundsException
 * @brief Exception raised when trying to access an element out of the vector's bounds
 */

	class IndexOutOfBoundsException : public std::exception {
		public:
			IndexOutOfBoundsException() {}
			
			const char* what() const throw() {
				return "Index out of bounds";
			}
	};

/**
 * @class SingularMatrixException
 * @brief Exception raised when a square matrix is singular
 */

	class SingularMatrixException:public std::exception {
		public:
			SingularMatrixException() {}
			
			const char* what() const throw() {
				return "The matrix is singular";
			}
	};
}

#endif
#endif

