# Required
- gcc compiler
# Compile
- `bash compile1.sh` - hw1
- `bash compile2.sh` - hw2
# Problems - hw1
1. Doesn't have any associated code
1. ./hw1 2 \<system to solve (A augmented with b)\> \<solution\>
1. ./hw1 3 \<matrix to decompose\> -> if there is a NaN in the result its not decomposeable
1. ./hw1 4 \<system to solve (A augmented with b)\> \<solution\>
1. ./hw1 5 \<system to solve (A augmented with b)\>

# Problems - hw2
1. ./hw2 1 \<method (Jacobi, Gauss-Seidel or SOR)\> \<max number of iterations\> \<system to solve\> \<omega (only required for SOR method)\> | ./hw2 2 \<system to solve\> (second command is for generating data for the plot)
3. ./hw2 3 \<matrix to find dominant eigenvalue for\> \<inital vector (optional)\>

display1.py and display2.py were used to generate the plots for the first and third problems

# Notes
- include/matrix.h contains comments for what the different matrix functions do
- src/hw1.c contains the implementations of each problem
- There are not a lot of comments in src/hw1.c because the function and variable names are as explicit as they sound, and the code "chunked" according to sections of logic (one chunk might be parsing the input, then the next is reducing the matrices, then the final chunk prints the result)
- Input format: "r1c1, r1c2, r1c3; r2c1, r2c2, r2c3; r3c1, r3c2, r3c3"
- **All input matrices must be surrounded with parenthesis**
