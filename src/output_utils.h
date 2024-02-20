#ifndef OUTPUT_UTILS_H
#define OUTPUT_UTILS_H

#include "data_structures.h"
#include "topology_utils.h"

#include <complex>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

class output_utils{
    public:

    bool itHeader = false;
    bool resHeader = false;
    bool stHeader = false;
    std::string outputFolder = "../OutputFiles/";

    int idSIM;

    void writeMatrixToFile(std::vector<std::vector<double> > Matrix)
    {
        int rows = Matrix.size();
        int cols = Matrix[0].size();

        std::ofstream outputFile (outputFolder + "Matrix.txt");
        
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
        std::ofstream outputFile (outputFolder + "Vector.txt");
        

        for (i = 0; i < rows; i++)
        {
            outputFile << i << ',' << Vector[i] << std::endl;
        }

        outputFile.close();
    }

    void appendIterationInfo(ESTIMATOR_TYPE esType, int itCount, double duration, double normaX, int order_method)
    {
        int simID = 1;
        std::ofstream outfile;
          
        outfile.open(outputFolder + "iterationInfo_" + std::to_string(idSIM) + ".txt", std::ios_base::app);
        if (itHeader==false)
        {
            outfile << "METHOD;" << "ITERATION" << ";" << "TIME" << ";" << "NORMX" << ";" << "ORDERING" << std::endl;
            itHeader = true;
        }
        if (esType == WLS)
        {
            outfile << "WLS;" << itCount << ";" << duration << ";" << normaX << ";" << order_method << std::endl;
        }
    }

    void writeDMEDResult(DMED *medidas, int nmed)
    {
        int i;
        std::ofstream outfile;
          
        outfile.open(outputFolder + "residueInfo_" + std::to_string(idSIM) + ".txt", std::ios_base::app);
        if (resHeader == false){
            outfile << "ID;" << "ZMED" << ";" << "HX" << ";" << "RESIDUE"  << std::endl;
            resHeader = true;
        }
        for (i = 0; i<nmed; i++)
        {
            int idMed = medidas[i].id;
            double zmed = medidas[i].zmed;
            double hx = medidas[i].h;
            double residue = zmed - hx;
            outfile << idMed <<";"<< zmed << ";" << hx << ";" << residue << std::endl;

        }
    }

    void writeEstimatedState(GRAFO *grafo, int numeroBarras)
    {
        int i;
        std::ofstream outfile;
          
        outfile.open(outputFolder + "stateInfo_" + std::to_string(idSIM) + ".txt", std::ios_base::app);
        if (stHeader == false)
        {
            outfile << "IDNODE;" << "FASE" << ";" << "REAL" << ";" << "IMAG" << std::endl;
            stHeader = true;
        }

        for (i = 0; i<numeroBarras; i++)
        {
            int id = grafo[i].idNo;
            int f = 0;
            for (f = 0; f<3; f++)
            {
                double rV = grafo[i].V[f].real();
                double iV = grafo[i].V[f].imag();
                outfile << id << ";" << f << ";" << rV << ";" << iV << std::endl;
            }
            
        }
    }

    void writeReferenceFile(GRAFO *grafo, long int numeroBarras, double Sbase){
    int i,j,k;
    
    FILE *arquivo;
    
    arquivo = fopen("refPowerflow.txt","wt");
    
    //----------------------------------------------------------------------
    //Fluxo de potência ativa em kW
    for(i=0;i<numeroBarras;i++){    
        //Percorre os ramos adjacentes
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            
            switch (grafo[i].adjacentes[k].ramo->fases){
                case 1:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).real() * Sbase);
                    break;
                case 2:
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).real() * Sbase);
                    break;
                case 3:
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).real() * Sbase);
                    break;
                case 4:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).real() * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).real() * Sbase);
                    break;
                case 5:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).real() * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).real() * Sbase);
                    break;
                case 6:
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).real() * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).real() * Sbase);
                    break;
                case 7:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).real() * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).real() * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).real() * Sbase);
                    break;
            }
        }        
    }
    //----------------------------------------------------------------------
    //Fluxo de potência reativa em kVAr
    for(i=0;i<numeroBarras;i++){    
        //Percorre os ramos adjacentes
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            
            switch (grafo[i].adjacentes[k].ramo->fases){
                case 1:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).imag() * Sbase);
                    break;
                case 2:
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).imag() * Sbase);
                    break;
                case 3:
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).imag() * Sbase);
                    break;
                case 4:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).imag() * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).imag() * Sbase);
                    break;
                case 5:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).imag() * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).imag() * Sbase);
                    break;
                case 6:
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).imag() * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).imag() * Sbase);
                    break;
                case 7:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[0]).imag() * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[1]).imag() * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(grafo[i].adjacentes[k].S[2]).imag() * Sbase);
                    break;
            }
        }        
    }
    
    
    
    //----------------------------------------------------------------------
    //Injeção de potência ativa em kW
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).real() * Sbase);
                break;
            case 2:
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).real() * Sbase);
                break;
            case 3:
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).real() * Sbase);
                break;
            case 4:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).real() * Sbase);
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).real() * Sbase);
                break;
            case 5:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).real() * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).real() * Sbase);
                break;
            case 6:
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).real() * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).real() * Sbase);
                break;
            case 7:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).real() * Sbase);
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).real() * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).real() * Sbase);
                break;
        }
    }
    //----------------------------------------------------------------------
    //Injeção de potência ativa em kVAr
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).imag() * Sbase);
                break;
            case 2:
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).imag() * Sbase);
                break;
            case 3:
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).imag() * Sbase);
                break;
            case 4:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).imag() * Sbase);
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).imag() * Sbase);
                break;
            case 5:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).imag() * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).imag() * Sbase);
                break;
            case 6:
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).imag() * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).imag() * Sbase);
                break;
            case 7:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[0]).imag() * Sbase);
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[1]).imag() * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(grafo[i].S[2]).imag() * Sbase);
                break;
        }
    }
    //----------------------------------------------------------------------
    //Magnitudes de tensão em kV
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[1]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[1]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*std::abs(grafo[i].V[2]));
                break;
        }
    }
//    //----------------------------------------------------------------------
//    //Ângulo de tensão em graus
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[1]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[1]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*std::arg(grafo[i].V[2]));
                break;
        }
    }
    fclose(arquivo);
}


    private:
};

#endif