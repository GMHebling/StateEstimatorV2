#include "data_structures.h"
#include "mat_utils.h"
#include <complex>
#include <iostream>

int main()
{
    mat_utils mu;

    //Test - Complex Matrix allocation
    // std::complex<double> **T;
    // int n = 3;
    // T = mu.c_matAloca(n);
    // int i, j;
    // for (i = 0; i < n; i ++)
    // {
    //     for (j = 0; j < n; j++)
    //     {
    //         std::cout << T[i][j].real() << "\n";
    //     }
    // }

    std::complex<double> **T2, **result, **T3, **T4, **T5, **T6;
    int n = 2;
    T2 = mu.c_matAloca(n);
    T3 = mu.c_matAloca(n);
    T4 = mu.c_matAloca(n);
    result = mu.c_matAloca(n);

    std::complex<double> c(1.0, 1.0);
    int i, j;
    for (i = 0; i < n; i ++)
    {
        for (j = 0; j < n; j++)
        {
            T2[i][j] = std::complex<double> (1.0,1.0);
            T3[i][j] = std::complex<double> (1.0,1.0);
        }
    }
    T2 = mu.c_matMultEsc(T2, c, n);
    // std::cout << T2[0][0] << "\n";
    
    T2 = mu.c_matConj(T2, n);
    // std::cout << T2[0][0] << "\n";

    T4 = mu.c_matIgual(T2, n);
    // std::cout << T2[0][0] << "\n";

    T5 = mu.c_matAloca(3);
    T6 = mu.c_matAloca(3);

    
    T5[0][0] = 1;
    T5[1][1] = 3;
    T5[2][2] = 7;

    T6 = mu.c_matInversaZ3(T5, 3);
    std::cout << T6[0][0] << "\n";
    std::cout << T6[1][1] << "\n";
    std::cout << T6[2][2] << "\n";
}