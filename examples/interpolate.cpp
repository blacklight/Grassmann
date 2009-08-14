/**
 * Little program to show how to use Grassmann to do some basic
 * linear or polynomial interpolation operations.
 * After the installation of Grassmann, compile it using
 * g++ -o interpolate interpolate.cpp -lgrassmann
 */

#include <iostream>
#include <grassmann.hpp>

using namespace std;
using namespace grassmann;

main()  {
	double x = 0.5;
	vector<Vector> points(3);

	points[0] = "-2,4";
	points[1] = "0,0";
	points[2] = "1,1";

	try  {
		cout << "Polynomial interpolation in x = " << x << ": "
			<< polynomialInterpolation(x,points) << endl
			<< "Linear interpolation in x = " << x << ": "
			<< linearInterpolation(x,points) << endl;
	}

	catch (exception e)  {
		cerr << "Oops...an error occurred while computing the interpolation:" << endl
			<< e.what() << endl;
	}
}

