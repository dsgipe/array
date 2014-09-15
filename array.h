//********************************************
// Specification file for Arr class
//********************************************
#include "lpkinterface.h"//Used for class Arr
#include <iostream>
#include "lpkinterface.h"
#include "cstddef"//for null
#include <cmath>
#include <iomanip>
using namespace std;
class Arr{
public:
    //constructors and destructors
    Arr();
    Arr(const Arr& obj);//copy constructor
    Arr(double* valInit, int m,int n);
    //Initializers
    void Init(double* valInit, int m,int n);
    void Init(int m,int n);
    //Operators
    Arr& operator=(const Arr& obj);
    Arr operator/(const Arr& obj);//matrix divide
    Arr operator*(const Arr& obj);//matrix multiply
    Arr operator%(const Arr& obj);//divide each value
    Arr operator%(const double intval);//divide each value
    Arr operator,(const Arr& obj);  //multiple each value
    Arr operator,(const int intval);//multiple each value
    Arr operator+(const Arr& obj);  //multiple each value
    Arr operator+(const double intval);//multiple each value
    ~Arr();
    void print();
    void print(const char * message);
    double * val;
    int M;
    int N;
    Arr transpose();
    Arr cholesky();
    double & element(int i, int j);
};
Arr times(const double intval, const Arr& obj);
struct Array{
    double * val;
    int size;
};

Arr mu_num_den(Arr& UPsiX,Arr& Y,Arr& oneN);
void buildPsi(Arr& x, double* theta, Arr& CKPsixRtn);
void Cholesky_arr(const Arr& S,Arr& D);


