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
    Arr(double valInit, int m,int n);
    Arr(int m,int n);
    //Initializers
    void Init(double* valInit, int m,int n);
    void Init(double valInit, int m,int n);
    void Init(int m,int n);
    //Operators
    Arr& operator=(const Arr& obj);
    Arr operator/(const Arr& obj);//matrix divide
    Arr operator*(const Arr& obj);//matrix multiply
    Arr operator%(const Arr& obj);//divide each value
    Arr operator%(const double intval);//divide each value
    Arr operator,(const Arr& obj);  //multiply each value
    Arr operator,(const int intval);//multiply each value
    Arr operator+(const Arr& obj);  //add each value
    Arr operator+(const double intval);//add each value
    Arr operator-(const Arr& obj);  //subtract each value
    ~Arr();
    void print();
    void print(const char * message);
    double * val;
    int M;
    int N;
    Arr transpose();
    Arr cholesky();
    double & element(int i, int j);
    void push(double value, int i,int j); 
};
Arr times(const double intval, const Arr& obj);
struct Array{
    double * val;
    int size;
};

Arr concatinate(Arr& obj1, Arr& obj2, int dim);


