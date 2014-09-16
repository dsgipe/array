array
=====

Array class that performs MATLAB like operations to ease code translations from MATLAB to C/C++
dependencies: lapack

uses: One of the big obstacles holding back companies and engineers from moving to C/C++ is a dependency on MATLAB. Many engineering students have become dependent on MATLAB, but are frustrated with speed and memory obsticles. It is the author's goal, to create an easy method of converting MATLAB code to C/C++ to make use of the significant efficiency improvement.

To do: Create an Init overloaded operator to fill array with a constant, to mimick ones and zeros behavior in MATLAB
   create a concatinate function
   Determine whether or not [] can be used to concatinate Arr similar to MATLAB
