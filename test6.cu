#include <iostream>     // std::cout

#include <cusp/csr_matrix.h>
#include <cusp/precond/diagonal.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/gmres.h>
#include <cusp/print.h>
#include <typeinfo>

// where to perform the computation
typedef cusp::device_memory MemorySpace;

// which floating point type to use
typedef double ValueType;

using namespace std;

extern "C"
{
    void fortran_solve_csr_(int n,int m,int nnz,int *rowoffset,int *col,double *val,double *rhs,double *x)
    {
        // create an empty sparse matrix structure (CSR format)
        cusp::csr_matrix<int, ValueType, MemorySpace> A(n,m,nnz);

        cout << "Row offset\n";
        for(int i=0;i<=n;i++)
        {
            A.row_offsets[i] = rowoffset[i];
            //cout << rowoffset[i] << " ";
        }

        cout << "cloumn and val\n";
        for(int i=0;i<nnz;i++)
        {
            A.column_indices[i] = col[i]; A.values[i] = val[i];

            //cout << col[i] << " " << val[i] << "\n";
        }
        
        //cusp::print(A);

        // allocate storage for solution (x) and right hand side (b)
        cusp::array1d<ValueType, MemorySpace> X(A.num_rows, 0);
        cusp::array1d<ValueType, MemorySpace> B(A.num_rows, 0);

        for(int i=0;i<n;i++)
        {
            B[i] = rhs[i];
        }
        //cusp::print(B);

         //std::cout << typeid(A.row_offsets).name() << '\n';
        // set stopping criteria:
        //  iteration_limit    = 20000
        //  relative_tolerance = 1e-15
        //  absolute_tolerance = 1e-10
        cusp::monitor<ValueType> monitor(B, 20000, 1e-15, 1e-10, false);

        // setup preconditioner
        /*cusp::precond::diagonal<ValueType, MemorySpace> M(A);
        //cusp::identity_operator<ValueType, MemorySpace> M(A.num_rows, A.num_rows);

        // solve the linear system A * x = b with the BiConjugate Gradient Stabilized method
        cusp::krylov::bicgstab(A, X, B, monitor, M);*/

        cusp::identity_operator<ValueType, MemorySpace> M(A.num_rows, A.num_rows);

        cusp::krylov::gmres(A, X, B, 50, monitor, M);


        for(int i=0;i<n;i++)
        {
            x[i] = X[i];
        }
    }
}