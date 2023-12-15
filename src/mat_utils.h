#ifndef MAT_UTILS_H
#define MAT_UTILS_H
#include <complex>
#include <cstdlib>
class mat_utils {
    public:
        std::complex<double>** c_matAloca(int n)
        {
            int i, j;
            
            std::complex<double>** output = new std::complex<double> *[n];
            // std::complex<double>* row = new std::complex<double>[n];
            for (i = 0; i < n; i++)
            {
                output[i] = new std::complex<double>[n];
                for (j = 0; j < n; j++)
                {
                    output[i][j] = std::complex<double>(0,0);
                }
            }
            return output;
        }

        std::complex<double>* c_vetAloca(int n)
        {
            int j;
            std::complex<double>* row = new std::complex<double>[n];
            for (j = 0; j < n; j++){ 
                row[j] = std::complex<double> (0,0); 
            }

            return row; 
        }

        std::complex<double>** c_matMultEsc(std::complex<double> **A, std::complex<double> b, int n)
        {
            int i,j;
            std::complex<double> **output;
            output = c_matAloca(n);
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    output[i][j] = b*A[i][j];
                }
            }   
            return (output);
        }

        std::complex<double>** c_matConj(std::complex<double> **A, int n)
        {
            int i,j;
            std::complex<double> **output;
            output = c_matAloca(n);
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    output[i][j] = std::conj(A[i][j]);
                }
            }   
            return (output);
        }
        std::complex<double>** c_matIgual(std::complex<double> **B, int n)
        {
            int i,j;
            std::complex<double> **output;
            output = c_matAloca(n);
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    output[i][j] = B[i][j];
                }
            } 
            return (output);  
        }

        std::complex<double>** c_matInversaZ(std::complex<double> **A, int n)
        {
            
            int i,j,k;
            std::complex<double> **Zaux, piv;
                
            //Aloca Matriz
            Zaux = c_matAloca(3);
            Zaux[0][0] = 1;
            Zaux[1][1] = 1;
            Zaux[2][2] = 1;
            
            //Etapa Forward
            for (k=0; k<3-1; k++)
            {
                for (i=(k+1); i<3; i++)
                {
                    if (A[k][k] != std::complex<double>(0,0) && A[i][k] != std::complex<double>(0,0)){
                        piv=-((A[i][k])/(A[k][k]));
                        //printf("%lf %lf\n", real(piv), imag(piv));
                        //getchar();
                        for (j=0; j<3; j++)
                        {                      
                            A[i][j]=(A[i][j]+(piv*A[k][j]));
                            Zaux[i][j]=(Zaux[i][j]+(piv*Zaux[k][j]));
                        }
                    }
                }
            }
            //Etapa Diagonalização
            for (i=0; i<3; i++)
            {
                if (A[i][i]!= std::complex<double>(0,0)){
                    piv=(std::complex<double>(1,0)/(A[i][i]));
                    // printf("%lf %lf\n", real(piv), imag(piv));
                    for (j=0; j<3; j++)
                    {
                        A[i][j]=(piv*A[i][j]);
                        Zaux[i][j]=(piv*Zaux[i][j]);
                    }
                }
                else{
                    Zaux[i][i]=std::complex<double>(0,0);
                }
            }
            //Etapa Backward
            for (k=0; k<3-1; k++)
            {
                for (i=(k+1); i<3; i++)
                {
                    piv=-(A[k][i]);
                    //printf("%lf %lf\n", real(piv), imag(piv));
                    for (j=0; j<3; j++)
                    {
                            A[k][j]=(A[k][j]+(piv*A[i][j]));
                            Zaux[k][j]=(Zaux[k][j]+(piv*Zaux[i][j]));
                    }

                }
            }
            return (Zaux);
        }

        std::complex<double>** c_matTransp(std::complex<double>  **A, int n)
        {
            int i,j;
            std::complex<double> **At;
            At = c_matAloca(n);
            
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    At[i][j] = A[j][i];
                }
            }
            return(At);

        }

        std::complex<double>** c_matInversaZ3(std::complex<double> **A,  int n)
        {
            std::complex<double> **Inv, det;
            std::complex<double> a, b, c, d, e, f, g, h, i;
            Inv = c_matAloca(n);
            a = A[0][0];
            b = A[1][0];
            c = A[2][0];
            d = A[0][1];
            e = A[1][1];
            f = A[2][1];
            g = A[0][2];
            h = A[1][2];
            i = A[2][2];

            if (a == 0.0)
            {
                a = 1.0;
            }

            if (e == 0.0)
            {
                e = 1.0;
            }

            if (i == 0.0)
            {
                i = 1.0;
            }

            det = a*(e*i - f*h) + b*(f*g - d*i) + c*(d*h - g*e);

            if (det == std::complex<double>(0,0)) {
                exit(EXIT_FAILURE);
            }
            Inv[0][0] = (e*i - f*h)/det;
            Inv[1][0] = (b*i - c*h)/det;
            Inv[2][0] = (c*e - b*f)/det;
            Inv[0][1] = (d*i - f*g)/det;
            Inv[1][1] = (a*i - c*g)/det;
            Inv[2][1] = (a*f - c*d)/det;
            Inv[0][2] = (g*e - d*h)/det;
            Inv[1][2] = (a*h - b*g)/det;
            Inv[2][2] = (a*e - b*d)/det;

            return(Inv);
        }

        double norma_inf(std::vector<double> a)
        {
            int i;
            double c;
            int n = a.size();

            c = std::abs(a[0]);
            for (i=0;i<n;i++)
            {
                if (std::abs(a[i]) > c) c = std::abs(a[i]);
            }

            return(c);
        }

    private:
};
#endif