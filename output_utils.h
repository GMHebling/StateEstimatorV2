#ifndef OUTPUT_UTILS_H
#define OUTPUT_UTILS_H

#include "data_structures.h"
#include "topology_utils.h"

#include <complex>
#include <iostream>
#include <fstream>
#include <string>

class output_utils{
    public:

    void writeMatrixToFile(std::vector<std::vector<double> > Matrix)
    {
        int rows = Matrix.size();
        int cols = Matrix[0].size();

        std::ofstream outputFile ("Matrix.txt");
        
        int i,j;

        for (i = 0; i < rows; i++)
        {
            for (j = 0; j < cols; j++)
            {
                if (Matrix[i][j] != 0.0)
                {
                    outputFile << i << ',' << j << ',' << Matrix[i][j] << std::endl;
                }
                
            }
        }
        outputFile.close();
    }

    void writeVectorToFile(std::vector<double> Vector)
    {
        int rows = Vector.size();
        int i;
        std::ofstream outputFile ("Vector.txt");
        

        for (i = 0; i < rows; i++)
        {
            outputFile << i << ',' << Vector[i] << std::endl;
        }

        outputFile.close();
    }


    private:
};

#endif