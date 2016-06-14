#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
#include <complex>
#include <vector>
#include <math.h>
using namespace std;


//Function prototypes
vector<complex<long double>> FFT(vector<complex<long double>>, bool, int);
void computeRootsofUnity(int);
void randomPoly(vector<complex<long double>> Poly, int);
complex<long double> operator*(complex<long double>& lhs, complex<long double>& rhs);

//Global consts and variables
const int SIZE = 1024;
vector<complex<long double>> RootsofUnity;
vector<complex<long double>> InverseRootsofUnity;

int main(void) {

	//Polynomials and Variables
	cout << "Allocating first array" << endl;
	vector<complex<long double>> Poly1;
	Poly1.resize(SIZE);
	cout << "Allocating second array" << endl;
	vector<complex<long double>> Poly2;
	Poly2.resize(SIZE);
	vector<complex<long double>> PQ;
	PQ.resize(SIZE);
	double totalTime = 0;
	int AVERAGE = 10;


	//Find the roots of unity that will be used for FFT()
	computeRootsofUnity(SIZE);

	//Randomize the coefficients in both polynomials
	randomPoly(Poly1, SIZE/2);
	randomPoly(Poly2, SIZE/2);

	//Pad Poly1 and Poly2 with 0s for FFT Algorithm
	for (int i = SIZE / 2 + 1; i < SIZE; i++) {
		Poly1[i].real(0);
		Poly2[i].real(0);
	}


	//Multiply the two polynomials
	//AVERAGE over the total number of runs
	for (int i = 0; i < AVERAGE; i++) {
		cout << "i: " << i << std::endl;
		clock_t begin = clock();	//start collecting time

		vector<complex<long double>> Sol_P;
		Sol_P.resize(SIZE);
		Sol_P = FFT(Poly1, false, SIZE);
		
		
		vector<complex<long double>> Sol_Q;
		Sol_Q.resize(SIZE);
		Sol_Q = FFT(Poly2, false, SIZE);
		
		
		vector<complex<long double>> Sol_PQ;
		Sol_PQ.resize(SIZE);

		for (int i = 0; i < SIZE; i++) {
			Sol_PQ[i] = Sol_P[i] * Sol_Q[i];
		}

		//Interpolation using Fast Fourier Transform
		//true is sent so that InverseRootsofUnity will be used
		PQ = FFT(Sol_PQ, true, SIZE);

		//Divide each coefficient by n to get the product
		for (int i = 0; i < SIZE; i++) {
			PQ[i].real( PQ[i].real() / SIZE);
		}

		clock_t end = clock();	//end collecting time
		double elapsed_secs = static_cast<double>(end - begin) / CLOCKS_PER_SEC;

		if (elapsed_secs == 0) {
			i--;
		}
		else {
			totalTime += elapsed_secs;
		}
	}

	totalTime /= AVERAGE;
	cout << "Average Seconds: " << totalTime << endl;

	//No longer used for large polynomial numbers for timed experiments
	//
	//std::cout << "The first polynomial P = ";
	//printPoly(Poly1, n);
	//std::cout << "The second polynomial Q = ";
	//printPoly(Poly2, n);
	//std::cout << "PQ = ";
	//printPoly(PQ, SIZE);

	std::cout << "Program ending..." << std::endl;
	std::cin.get();
	return 0;

}


//Fast Fourier Transform algorithm to multiply two polynomials. O(n log n)
vector<complex<long double>> FFT(vector<complex<long double>> Poly, bool useInverse, int n) {

	if (n == 1) {
		vector<complex<long double>> temp;
		temp.push_back(Poly[0]);
		return temp;
	}
	else {
		vector<complex<long double>> Poly_even;
		Poly_even.resize(n / 2);
		vector<complex<long double>> Poly_odd;
		Poly_odd.resize(n / 2);


		//initialize arrays into the even and odd sections
		for (int i = 0; i < n / 2; i++) {
			Poly_even[i].real(Poly[2*i].real());
			Poly_even[i].imag(Poly[2*i].imag());
			Poly_odd[i].real(Poly[2 * i + 1].real());
			Poly_odd[i].imag(Poly[2 * i + 1].imag());
		}

		//Solutions of the odd and even sides of the arrays.
		vector<complex<long double>> Sol_even;
		Sol_even.resize(n);
		vector<complex<long double>> Sol_odd;
		Sol_odd.resize(n);


		//initialize new arrays
		for (int i = 0; i < n; i++) {
			Sol_even[i].real(0);
			Sol_even[i].imag(0);
			Sol_odd[i].real(0);
			Sol_odd[i].imag(0);
		}

		//Recursively solve the solutions to the odd and even sides of the array
		Sol_even = FFT(Poly_even, useInverse, n / 2);
		Sol_odd = FFT(Poly_odd, useInverse, n / 2);

		//Allocate new solution vector to calculate values
		//initialize solution to all 0.
		vector<complex<long double>> Solution;
		Solution.resize(n);
		for (int i = 0; i < n; i++) {
			Solution[i].real(0);
			Solution[i].imag(0);
		}	

		int shift = RootsofUnity.size() / n;

		//Calculate the solution using the even and odd sides of the array
		//The odd must be multiplied by the Roots of Unity squared.
		//Shift up the array based on the size of Poly[]
		//Example:
		//[w0, w1, w2, w3, w4, w5, w6, w7]
		//[w0, w2, w4, w6]
		//[w0, w4]
		if (useInverse == true) {
		
			for (int i = 0; i < n/2; i++) {
				Solution[i] = Sol_even[i] + InverseRootsofUnity[i*shift] * Sol_odd[i];
				Solution[i + n / 2] = Sol_even[i] - InverseRootsofUnity[i*shift] * Sol_odd[i];
			}
		} 
		else {

			for (int i = 0; i < n/2; i++) {
				Solution[i] = Sol_even[i] + RootsofUnity[i*shift] * Sol_odd[i];
				Solution[i + n / 2] = Sol_even[i] - RootsofUnity[i*shift] * Sol_odd[i];
			}
		}


		return Solution;
	}
}

//Computes the roots of unity angles based on the number "num" of Coefficients
void computeRootsofUnity(int num) {

	//Resize the array to hold all angles of Omega "W"
	RootsofUnity.resize(num);
	InverseRootsofUnity.resize(num);

	for (int i = 0; i < num; i++) {
		RootsofUnity[i].real(cos((2.0 * M_PI * static_cast<long double>(i)) / num));
		InverseRootsofUnity[i].real(RootsofUnity[i].real());
		RootsofUnity[i].imag(sin((2.0 * M_PI * static_cast<long double>(i)) / num));
		InverseRootsofUnity[i].imag(RootsofUnity[i].imag() * -1.0); //Invert the imaginary portion
	}
}

//Overloaded * operator for multiplying complex numbers.
complex<long double> operator*(complex<long double>& lhs,complex<long double>& rhs) {
	/*
	//Multiplcation, 4 times technique
	long double a, b;
	complex<long double> temp;
	a = (lhs.real() * rhs.real()) - (lhs.imag() * rhs.imag());
	b = (lhs.real() * rhs.imag()) + (rhs.real() * lhs.imag());
	temp.real(a);
	temp.imag(b);
	*/
	
	//Multiplication, 3 times technique
	long double real, imaginary, combine;
	complex<long double> temp;

	real = lhs.real() * rhs.real();
	imaginary = lhs.imag() * rhs.imag();
	combine = (lhs.real() + lhs.imag()) * (rhs.real() + rhs.imag());

	temp.real(real - imaginary);
	temp.imag(combine - (real + imaginary));
	
	return temp;
}

//Generate random polynomial coefficients
void randomPoly(vector<complex<long double>> Poly, int n) {

	int maxSize = 10;
	int minSize = -10;

	srand(time(NULL));

	//Random value range from minSize to maxSize.
	for (int i = 0; i < n; i++) {
		Poly[i] = minSize + (std::rand() % (maxSize - minSize + 1));
		Poly[i] /= 10; //range between -1.0 and 1.0
	}
}
