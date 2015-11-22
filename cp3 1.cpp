#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

using namespace std;

//need to solve Ax=b for x, where A is linear equations, x is current values, b
//voltage (divided by ohms so units work out)

MatDoub setA(int dim, vector<double> resist){

        MatDoub a(dim, dim);

        //for loop will go down array by line, setting each linear equation for A
        for(int i = 0; i < dim; i++){


                for(int x = 0; x < dim; x++){

                        //accounts for Kirkoff's Current Law
                        if( i == 0){
                                a[i][x] = -1.0;

                                //makes sure I1 is positive
                                if(x == 0){
                                        a[i][x] = 1.0;
                                }
                        }

                        else{
                             	//makes sure first item in array is equal to R1
                                if(x == 0){
                                        a[i][x] = resist[0];
                                }
                                //if on the same row collum that corresponds with the ith resistance
                                //put resistance value into matrix
                                else if( x == i){
                                        a[i][x] = resist[i];
                                }else{
                                //otherwise, other resistance does not need to be accounted for
                                        a[i][x] = 0.0;
                                }
                        }


                }

        }
        return a;

}

VecDoub setB(int dim, vector<double> volts){

        VecDoub b(dim);

        //set first value to 0 because of KCL
        b[0] = 0.0;

        //sets rest as combination of V1 and voltage for specific linear combo
        for(int i = 1; i < dim; i++){
                b[i] = volts[0] + volts[i];
        }

        return b;
}


int main() {
 //need to set value of n (number of eq'ns and dimensions of square matrix)
        double n_in;
        cout << "please input n (max number of equations):\n";
        cin >> n_in;

        const double n = n_in;
        const int nN = (int) n;


        //this component of code dependent on parts C, D, and E
        string part;
        cout << "which part of the problem are you working on?\n";
        cin >> part;

        if(part == "c"){

                //allocate a square matrix
                MatDoub a(nN,nN);

                //need to allocate vectors
                //b is right hand side (solutions to linear equations)
                //x is coefficients (in this case currents) that need to be solved for
                VecDoub b(nN), x(nN);


                //need vector that holds resistances
                vector<double> resistances(nN);

                //need vector that holds voltages
               vector<double> voltage(nN);

                //sets values for test of code
               for(int t = 0; t < n; t++){
                    cout << "\nplease input the " << t + 1 << " voltage value;\n";
                    cin >> voltage[t];
                    cout << "\nplease input the " << t + 1 << " resistance value:\n";
                    cin >> resistances[t];
               }

                //sets A based on given values
                a = setA(nN, resistances);

                //sets b based on voltages given
                b = setB(nN, voltage);
                 //solves using LU decomposition
                LUdcmp alu(a);
                alu.solve(b, x);

                //prints out current values for error checking
                for(int c = 0; c < n; c++){
                    cout << "current val " << c + 1 << ": " << x[c] << endl;
                }


        }
    else if(part == "d"){

            //writes part d to a .dat file to graph later
            ofstream myfilee;
            myfilee.open("partd.dat");
            myfilee << "#n    I1    predicted I1\n";


                for(double var = 1.0; var <= n; var++){
                        int n_var = (int) var;

                        //need to make array
                        MatDoub a(n_var, n_var);

                        //need to allocate vectors
                        //b is right hand side (solutions to linear equations)
                        //x is coefficients (in this case currents) that need to be solved for
                        VecDoub b(n_var), x(n_var);

                       //need vector that holds resistances
                        vector<double> resistances(n_var);

                        //need vector that holds voltages
                       vector<double> voltage(n_var);

                       //sets voltages and resistances based on premise of problem
                       for(int i = 0; i < n_var; i++){
                            voltage[i] = 0.0;
                            resistances[i] = 1.0;
                       }

                        //voltage for v1 is actually 1.0
                       voltage[0] = 1.0;

                        a = setA(n_var, resistances);
                        b = setB(n_var, voltage);

                        LUdcmp alu(a);
                        alu.solve(b,x);

                        double predictedVal = 1.0 - 1.0/var;

                        myfilee << var << " " << x[0] << " " << predictedVal << endl;
                }
            myfilee.close();
        }
 else if(part == "e"){
            ofstream myfile;
            myfile.open("data.dat");

          for(double var = 1.0; var <= n; var++){
                    //to avoid complications, have as both double and int
                    int n_var = (int) var;

                    //need to make array
                    MatDoub a(n_var, n_var);

                    //need to allocate vectors
                    //b is right hand side (solutions to linear equations)
                    //x is coefficients (in this case currents) that need to be solved for
                    VecDoub b(n_var), x(n_var);

                   //need vector that holds resistances
                    vector<double> resistances(n_var);

                    //need vector that holds voltages
                   vector<double> voltage(n_var);

                    //assigns voltages of 1 and resistances equivalent to ith rank in code
                   for(int i = 0; i < n_var; i++){
                        //converting loop into a double
                        double iD = (double) i;

                        voltage[i] = 1.0;
                        //need to add 1 since c++ counts from 0, since resistance
                        //in this problem corresponds to given index
                        resistances[i] = iD + 1.0;
                   }


                    a = setA(n_var, resistances);
                    b = setB(n_var, voltage);

                    LUdcmp alu(a);
                    alu.solve(b,x);

                    //writes data to file
                    myfile << var << " " << x[0] << endl;
            }

            myfile.close();
        }

 return 0;


}
