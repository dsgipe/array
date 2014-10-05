//------------------------------------------------
//       Arr class function descriptions
//-----------------------------------------------
#include "array.h"
Arr::Arr(){
    //constructor
    M=1;
    N=1;
    val = NULL;

 }
//************************************************
Arr::Arr(double* valInit,int m, int n){
    M = m;
    N = n;
    val = new double[m*n];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = valInit[ii];}
}
//************************************************
Arr::Arr(double valInit,int m, int n){
    M = m;
    N = n;
    val = new double[m*n];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = valInit;}
}
//************************************************
Arr::Arr(const Arr& obj){//copy constuctor
    M = obj.M;
    N = obj.N;
    val = new double[M*N];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = obj.val[ii];}
}
//************************************************
void Arr::Init(double* valInit,int m, int n){
    //constructor
    M = m;
    N = n;
    val = new double[m*n];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = valInit[ii];}
}
void Arr::Init(double valInit,int m, int n){
    //constructor
    M = m;
    N = n;
    val = new double[m*n];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = valInit;}
}

void Arr::Init(int m, int n){
    //constructor
    M = m;
    N = n;
    val = new double[m*n];
    for(int ii = 0; ii < M*N;ii++){ val[ii] = 0;}
}

//************************************************
Arr::~Arr(){
    if (val !=NULL){
       delete [] val;
       val = NULL;
    }
}
//************************************************
void Arr::print(const char * message){
    cout << message << endl;
    Arr::print();
}
void Arr::print(){
    //---------------------------------------------//
    //prints out what is a matrix in row major format
    //inputs:
    //     Write A, in columns of m, and rows of n
    //     outputs A to the  screen 
    //---------------------------------------------//
    int counter = 0;
    int b = 0;
    for(int ii = 0; ii < M;ii++){ 
        if(ii ==0){  cout << setw(9);}
        for(int jj = 0; jj < N;jj++){ 
            cout  << setiosflags(ios::fixed) << setprecision(2)<<  val[b+counter*N] ;
            cout << setw(9);
            //cout << setiosflags(ios::fixed) << setprecision(4) << A[b+counter*n]<< "\t";
            counter ++;
            if(counter == M){
               counter =0;
               b++;
               cout << endl;
            }
        }
    }
}
//************************************************
Arr Arr::transpose(){
    //---------------------------------------------//
    // Transpose a row ordered array, arr,
    // Inputs:
    //     arr: 1d array in column ordered format
    //     nc:  number of columns
    //     nr:  number of rows
    // Returns:
    //     transpose of arr
    //---------------------------------------------//
    Arr RtnArray;
    RtnArray.Init(N,M);
    for (int ii = 0; ii < M; ii++){
        for (int jj = 0; jj <N; jj++){
            RtnArray.push(element(jj,ii),ii,jj);
        }
    }
    return RtnArray;

}
Arr& Arr::operator=(const Arr& obj){
    if (val ==NULL){
        val = new double[M*N];
    }else{
        delete[] val;
        val = new double[M*N];
    }
    M = obj.M;
    N = obj.N;
    for(int ii = 0; ii < M*N;ii++){ val[ii] = obj.val[ii];}
}
//************************************************
Arr Arr::cholesky(){
    //---------------------------------------------//
    // Solve cholesky of a square array
    // Inputs: 
    //    S: input array
    //    D: output array. Cholesky of input array. 
    //---------------------------------------------//
    Arr D;
    D.Init(M,N);
    int d = M;
    D.M = M;D.N=N;
    for(int k=0;k<d;++k){
        double sum=0.;
        for(int p=0;p<k;++p)sum+=D.val[k*d+p]*D.val[k*d+p];
        D.val[k*d+k]=sqrt(val[k*d+k]-sum);
        for(int i=k+1;i<d;++i){
           double sum=0.;
           for(int p=0;p<k;++p)sum+=D.val[i*d+p]*D.val[k*d+p];
           D.val[i*d+k]=(val[i*d+k]-sum)/D.val[k*d+k];
        }
    }
    return D;
}
//************************************************
Arr Arr::operator/(const Arr& obj){
    //---------------------------------------------//
    //         Matrix solution to Y=M*X
    //---------------------------------------------//
    int NRHS = obj.N;
    int NLHS = obj.M;
    int *IPIV = new int[N+1];
    int INFO;
    //Create temporary variables so the fortran code doesn't change
    //input variables
    double A_tmp[M*N];for(int ii =0;ii<M*N;ii++){A_tmp[ii]=val[ii];}
    double B_rtn[obj.N*obj.M];for(int ii =0;ii<obj.N*obj.M;ii++){B_rtn[ii]=obj.val[ii];}
    dgesv_(&NLHS,&NRHS,A_tmp,&NLHS,IPIV,B_rtn,&N,&INFO);
    
    delete [] IPIV;
    return Arr (B_rtn,obj.M,obj.N);
}
//************************************************
Arr Arr::operator%(const Arr& obj){
    //---------------------------------------------//
    // divide every element
    //---------------------------------------------//
    double rtnArr[M*N]; 
    if (M*N != obj.M*obj.N)
        cout << "Size not compatible, results are likely wrong!\n";
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]/obj.val[ii];
    }
   return Arr (rtnArr,M,N); 
}
//************************************************
Arr Arr::operator%(const double intval){
    double rtnArr[M*N]; 
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]/intval;
    }
   return Arr (rtnArr,M,N); 
}
//************************************************
Arr Arr::operator+(const Arr& obj){
    //---------------------------------------------//
    // divide every element
    //---------------------------------------------//
    double rtnArr[M*N]; 
    if (M*N != obj.M*obj.N)
        cout << "Size not compatible, results are likely wrong!\n";
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]+obj.val[ii];
    }
   return Arr (rtnArr,M,N); 
}
//************************************************
Arr Arr::operator+(const double intval){
    double rtnArr[M*N]; 
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]+intval;
    }
   return Arr (rtnArr,M,N); 
}
Arr Arr::operator-(const Arr& obj){
    //---------------------------------------------//
    // divide every element
    //---------------------------------------------//
    double rtnArr[M*N]; 
    if (M*N != obj.M*obj.N)
        cout << "Size not compatible, results are likely wrong!\n";
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]-obj.val[ii];
    }
   return Arr (rtnArr,M,N); 
}
//************************************************

//************************************************
Arr times(const double intval, const Arr& obj){
    double rtnArr[obj.M*obj.N]; 
    for (int ii = 0; ii < obj.M*obj.N;ii++){
        rtnArr[ii]=obj.val[ii]*intval;
    }
   return Arr (rtnArr,obj.M,obj.N); 
}
//************************************************
Arr Arr::operator,(const int intval){
    double rtnArr[M*N]; 
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]*intval;
    }
   return Arr (rtnArr,M,N); 
}
//************************************************
Arr Arr::operator,(const Arr& obj){
    //---------------------------------------------//
    // muliply every element
    //---------------------------------------------//
    double rtnArr[M*N]; 
    if (M*N != obj.M*obj.N)
        cout << "Size not compatible, results are likely wrong!\n";
    for (int ii = 0; ii < M*N;ii++){
        rtnArr[ii]=val[ii]*obj.val[ii];
    }
   return Arr (rtnArr,M,N); 
}
//************************************************
Arr Arr::operator*(const Arr& obj){
    //---------------------------------------------//
    // Multiple A by B
    // Inputs: 
    //    A: MxK matrix
    //    B: K+N matrix
    //    M: size of A
    //    N: size of B
    //    K: size shared by A and B
    //  Returns:
    //    C: MxN matrix
    //---------------------------------------------//
    int K = obj.M;
    int N_input = obj.N;
    int M_input = N;
    if (obj.M!=M){
        cout << "Likely problem with matrix, please check results\n";
    }
    double A_tmp[K*M_input];for(int ii =0;ii<K*M_input;ii++){A_tmp[ii]=val[ii];}
    double B_tmp[K*N_input];for(int ii =0;ii<K*N_input;ii++){B_tmp[ii]=obj.val[ii];}
    double C_rtn[M_input*N_input];for(int ii =0;ii<<M_input*N_input;ii++){C_rtn[ii]=0;}
    char transa = 'n';
    char transb = 't';
    double alpha =1;double beta = 0;
    dgemm_(&transa, &transb, &M_input, &N_input, &K,&alpha,A_tmp,&M_input,B_tmp,&N_input,&beta,C_rtn, &K );
 
    return Arr (C_rtn,M_input,N_input);
}
//************************************************
double& Arr::element(int i, int j){
    //---------------------------------------------//
    // inputs:
    //    i: Row index
    //    j: Column index
    // Convert single array indexing to easy find 
    // row and column info
    //---------------------------------------------//
    return val[j*N+i];
}
void Arr::push(double value, int i, int j){
    //---------------------------------------------//
    // inputs:
    //    i: Row index
    //    j: Column index
    // Convert single array indexing to easy find 
    // row and column info
    //---------------------------------------------//
    val[j*N+i] = value;
}
Arr concatinate(Arr& obj1, Arr& obj2,int dim){
    Arr ArrRtn;
    int M,N;
    //If the concatinate rows
    if (dim == 1){
        M = obj1.M;
        N = obj1.N+obj2.N;
    }else{//concatinate columns
        M = obj1.M+obj2.M;
        N = obj1.N;
    }
    ArrRtn.Init(0.0,M,N);
    //fill using first array
    for (int ii = 0; ii < obj1.N; ii++){
        for (int jj = 0; jj < obj1.M; jj++){
            ArrRtn.push(obj1.element(ii,jj),ii,jj);
        }
    } 
    //---------------------------------------------//
    //       fill using second array
    //---------------------------------------------//
    //concatinate rows
    if (dim == 1){
        for (int ii = obj1.N; ii < N; ii++){
            for (int jj = 0; jj < M; jj++){
                ArrRtn.push(obj2.element(ii-obj1.N,jj),ii,jj);
            }
        } 
    }else{//concatinate columns
        for (int ii = 0; ii < N; ii++){
            for (int jj = obj1.M; jj < M; jj++){
                ArrRtn.push(obj2.element(ii,jj-obj1.M),ii,jj);
            }
        } 
    }
    return ArrRtn;
}

