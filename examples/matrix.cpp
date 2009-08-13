/**
 * Little program to get some information from a 2x2 arbitrary matrix
 * After the installation of Grassmann, compile it using
 * g++ -o matrix matrix.cpp -lgrassmann
 */

#include <iostream>
#include <grassmann.hpp>

using namespace std;
using namespace grassmann;

int main()  {
	Matrix A(2,2);

	for (size_t i=0; i < A.rows(); i++)
		for (size_t j=0; j < A.cols(); j++)  {
			cout << "Element [" << (i+1) << "][" << (j+1) << "]: ";
			cin >> A(i,j);
		}

	/**
	 * If you, instead, would like to define a matrix through preset values:
	 *
	 * Matrix A(2,2);
	 * A[0] = "1,2";
	 * A[1] = "3,2";
	 */

	cout << "\nMatrix:\n\n"
		<< A.toString() << endl << endl
		<< "rank: "
		<< A.rank() << endl << endl
		<< "product A*A:\n\n"
		<< (A*A).toString() << endl << endl;
	
	try  {
		cout << "determinant: "
			<< A.det() << endl << endl
			<< "inverse matrix:\n\n"
			<< A.inverse().toString() << endl << endl
			<< "eigenvalues (probably an approximation):\n"
			<< A.eigenValues().toString() << endl;
	}

	catch (SingularMatrixException e)  {
		cerr << "Ooops...there was an error while working with your matrix:\n"
			<< "\t-> " << e.what() << endl;
	}

	return 0;
}

