#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

using namespace std;

//returns conditional value based on given matrix
double findCond(MatDoub mat, int dim) {

	//the condition number is the maximum value of the sum at a specific row
	double conditional = 0.0;

	//walks through each row of array finding max at each, then sums it
	for(int i = 0; i < dim; i++){
		
		//initializes row value to be zero
		double rowVal = 0.0;

		for(int h = 0; h < dim; h++){

			//adds absolute value of item h in row i to total value of row
			rowVal = rowVal + abs(mat[i][h]);

		}

		//once total value of that row is calculated, compare to current condition
		//number value
		if(rowVal > conditional){
			conditional = rowVal;
		}

	}

	return conditional;

}

int main(){

	//so that way I can't mess up dimensions
	const int dim = 2;

	//initializes a matrix
	MatDoub a(dim, dim);

	//can set values of a since they are constant
	a[0][0] = 0.789;
	a[0][1] = 0.56301;
	a[1][0] = 1.182711;
	a[1][1] = 0.843952;

	//sets up a so it's been de-compositioned (is that a verb?)
    LUdcmp alu(a);


    //this component of code dependent on parts B or C
    string part;
    cout << "which part of the problem are you working on?\n";
    cin >> part;


    if(part == "b"){
    	//creates a matrix to find inverse of A
		MatDoub a_inv(dim, dim);

		//calculate inverse
		alu.inverse(a_inv);

		//I recognize this could have been done more elegantly with nested for loops
		//but it would have taken up as much space
		cout << a_inv[0][0] << "   ";
		cout << a_inv[0][1] << "\n";
		cout << a_inv[1][0] << "   ";
		cout << a_inv[1][1] << "\n";

		double conditionalA = findCond(a, dim);
		double conditional_Ainv = findCond(a_inv, dim);

		cout << "condition no. of A: " << conditionalA << endl;
		cout << "condition no. of A inverse: " << conditional_Ainv << endl;

		cout << "actual condition number is:  " << conditionalA * conditional_Ainv << endl;


    }
    else if(part == "c"){
    	//initalizes b and x vectors to actually solve for x
		VecDoub b(dim), x(dim);

		b[0] = 0.22599;
		b[1] = 0.338759;

		alu.solve(b,x);

		//outputs x values to compare precisions
		cout << setprecision(15) << x[0] << endl;
		cout << setprecision(15) << x[1] << endl;
    }

    else{
    	cout << "you typed in something wrong... answer b or c\n";
    }

	return 0;
}