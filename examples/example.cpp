/**
 * Little program to get some information from a 2x2 arbitrary matrix
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

	cout << "\nMatrix:\n\n"
		<< A.toString() << endl
		<< "rank: "
		<< A.rank() << endl << endl
		<< "product A*A:\n\n"
		<< (A*A).toString() << endl;
	
	try  {
		cout << "determinant: "
			<< A.det() << endl << endl
			<< "inverse matrix:\n\n"
			<< A.inverse().toString() << endl
			<< "eigenvalues (probably an approximation):\n\n"
			<< A.eigenValues().toString();
	}

	catch (SingularMatrixException e)  {
		cerr << "Ooops...there was an error while working with your matrix:\n"
			<< "\t-> " << e.what() << endl;
	}

	return 0;
}

