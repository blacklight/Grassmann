/**
 * Simple program to show how to manage vectors using Grassmann
 * To compile it:
 * g++ -o vector vector.cpp -lgrassmann
 */

#include <iostream>
#include <grassmann.hpp>

using namespace std;
using namespace grassmann;

main()  {
	// Two vectors, each with 3 elements
	Vector u(3), v(3);

	// First vector -> you can define its constant values as a simple string.
	// Each value is separated by another by a comma ","
	u = "1,2,0";

	// You can, for example, get the second vector via stdin
	for (size_t i=0; i < 3; i++)  {
		cout << "Element #" << (i+1) << ": ";
		cin >> v[i];
	}
	
	cout << endl << "u: " << u.toString() << endl
		<< "v: " << v.toString() << endl << endl
		<< "v's modulus: " << v.modulus() << endl << endl
		<< "u+v: " << (u+v).toString() << endl << endl
		<< "u-v: " << (u-v).toString() << endl << endl
		<< "scalar product u*v: " << u*v << endl << endl
		<< "vectorial product u x v: " << (u%v).toString() << endl;
}

