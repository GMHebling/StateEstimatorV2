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

    std::complex<double> **T2, **result, **T3, **T4;
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
    std::cout << T2[0][0] << "\n";
    
    T2 = mu.c_matConj(T2, n);
    std::cout << T2[0][0] << "\n";

    T4 = mu.c_matIgual(T2, n);
    std::cout << T2[0][0] << "\n";

}