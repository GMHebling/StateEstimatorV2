#ifndef ENGINE_H
#define ENGINE_H

#include <cholmod.h>
#include "SuiteSparseQR_C.h"

#include <iostream>

enum solver_method {
    LU,
    CHOLESKY,
    SPARSEQR,
    GRADIENT
};

class engine {
    public:
        int solver_method = SPARSEQR;
        int order_method = 0;
        double tol = 1e-4;
        int max_it = 20;

        std::vector<double> engine_output;

        void engine_run(std::vector<std::vector<double> > A, std::vector<double> b)
        {
            //std::cout << "| START ENGINE RUN |\n";
            if (solver_method == LU)
            {
                //std::cout << "| USING LU FACTORIZATION |\n";
                engine_output = pivot_lu(A, b);
                //std::cout << "| RESULT |\n";
                
            }
            else if (solver_method == SPARSEQR)
            {
                //std::cout << "| USING SPARSE QR FACTORIZATION |\n";
                engine_output = sparse_QR(A, b, order_method);
                //std::cout << "| RESULT |\n";
            }
        }

        //LU
        //Choleksy
        //SparseQR
        //SteepestDescent
        //GradientDescent

    
    private:
        std::vector<double> pivot_lu(std::vector<std::vector<double> > A, std::vector<double> b)
        {
            int i,j;
            int n, m;
            int row, col;
            int col_count;
            int pivot_row;

            double max_element;
            double pivot;

            bool pivot_flag;

            std::vector<double> result;
            std::vector<double> temp_x;
            std::vector<double> y;
            std::vector<double> swap_row;
            

            std::vector<std::vector<double> > L;
            std::vector<std::vector<double> > U;
            std::vector<std::vector<double> > TA;

            n = A.size();
            m = A.size();

            std::vector<double> vector_row(n,0);
            for (i=0;i<m;i++)
            {
                result.push_back(0.0);
                temp_x.push_back(0.0);
                y.push_back(0.0);
                swap_row.push_back(0.0);
                L.push_back(vector_row);
                U.push_back(vector_row);
                TA.push_back(vector_row);
            }

            for (i=0;i<m;i++)
            {
                L[i][i]=1;
            }
            for (i=0;i<m;i++)
            {
                for (j=0;j<m;j++)
                {
                    TA[i][j] = A[i][j];
                }
            }

            for(col=0;col<m;col++){
                max_element = TA[col][col];
                pivot_row=0;
                pivot_flag=false;
                for(row=col;row<m;row++)
                {
                    if (TA[row][col] > max_element)
                    {
                        max_element = TA[row][col];
                        pivot_row = row;
                        pivot_flag = true;
                    }
                }
            

                if (pivot_flag)
                {
                    for (col_count=0;col_count<m;col_count++)
                    {
                        swap_row[col_count] = TA[pivot_row][col_count];
                        TA[pivot_row][col_count] = TA[col][col_count];
                        TA[col][col_count] = swap_row[col_count];
                    }
                    double aux_b=0.0;
                    aux_b = b[pivot_row];
                    b[pivot_row] = b[col];
                    b[col] = aux_b;
                }
                for (row=col+1;row<m;row++)
                {
                    pivot = -TA[row][col]/TA[col][col];
                    for (col_count=col;col_count<m;col_count++)
                    {
                        TA[row][col_count] = TA[row][col_count] + pivot*TA[col][col_count];
                    }
                    TA[row][col] = -pivot;
                }
            }

            //std::cout << "| FACTORIZATION DONE |\n";
            
            //obtain L and U matrices
            for(j=0;j<m;j++)
            {
                for(i=j+1;i<m;i++)
                {
                    L[i][j] = TA[i][j];
                }
            }

            for(j=0;j<m;j++)
            {
                for(i=0;i<j+1;i++)
                {
                    U[i][j] = TA[i][j];
                }
            }

            y = forwards_substitution(L, b);
            //std::cout << "| FORWARDS DONE |\n";
            temp_x = backwards_substitution(U, y);
            //std::cout << "| BACKWARDS DONE |\n";
            return (temp_x);

        }
        std::vector<double> forwards_substitution(std::vector<std::vector<double> > L, std::vector<double> b)
        {
            int m;
            int i,j;
            m = L.size();
            std::vector<double> y(m,0);
            double sum = 0.0;
            for (i=0;i<m;i++)
            {
                sum = b[i];
                for (j=0;j<i;j++)
                {
                    sum = sum - L[i][j]*y[j];
                }
                y[i] = sum/L[i][i];
            }
            return(y);

        }
        std::vector<double> backwards_substitution(std::vector<std::vector<double> > U, std::vector<double> y)
        {
            int m;
            int i,j;
            m = U.size();
            std::vector<double> x(m,0);
            double sum = 0.0;
            for (i=m-1;i>=0;i--)
            {
                sum = y[i];
                for (j=i+1;j<m;j++)
                {
                    sum = sum - U[i][j]*x[j];
                }
                x[i] = sum/U[i][i];
            }
            return(x);
        }
        std::vector<double> sparse_QR(std::vector<std::vector<double> > A, std::vector<double> b, int order_method)
        {
            int m = A.size();

            int rows = A.size();
            int cols = A[0].size();

            int i;
            double *Dx = NULL;
            cholmod_sparse *A_SS = NULL;
            cholmod_dense *b_SS = NULL;
            cholmod_dense *X_SS = NULL;
            cholmod_triplet *T_SS = NULL;
            cholmod_common Common, *c;
            c = &Common;
            cholmod_l_start(c);
            T_SS = cholmod_l_allocate_triplet(rows, cols, rows * cols, 0, CHOLMOD_REAL, c);
            A_SS = cholmod_l_allocate_sparse(rows, cols, rows * cols, 0, 0, 0, CHOLMOD_REAL, c);
            b_SS = cholmod_l_allocate_dense(rows, 1, rows, CHOLMOD_REAL, c);
            X_SS = cholmod_l_allocate_dense(cols, 1, cols, CHOLMOD_REAL, c);
            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int r = 0; r < cols; r++)
                {
                    if (A[i][r] != 0)
                    {
                        ((long int *)T_SS->i)[index] = i;
                        ((long int *)T_SS->j)[index] = r;
                        ((double *)T_SS->x)[index] = A[i][r];

                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
            }
            for (int i = 0; i < rows; i++)
            {
                ((double *)b_SS->x)[i] = b[i];
            }
            A_SS = cholmod_l_triplet_to_sparse(T_SS, rows * cols, c);
            X_SS = SuiteSparseQR_C_backslash(order_method, SPQR_DEFAULT_TOL, A_SS, b_SS, c);
            Dx = (double *)X_SS->x;
            cholmod_l_finish(c);
            std::vector<double> result(cols);
            for (i=0;i<cols;i++)
            {
                result[i] = Dx[i];
            }
            return (result);
            
        }
        
};
#endif