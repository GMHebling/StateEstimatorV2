#ifndef WLS_UTILS_H
#define WLS_UTILS_H

#include "SparseSystem.h"
#include "engine.h"
#include "data_structures.h"
#include "electrical_utils.h"
#include "input_data.h"
#include "mat_utils.h"
#include "output_utils.h"

#include <vector>
class wls_utils {
public:
  std::vector<double> z;  // vetor medidas
  std::vector<std::vector<double> > H; // Matriz Jacobiana
  std::vector<double> x;  // vetor de estado
  std::vector<double> Dx;
  std::vector<double> regua;

  std::vector<double> b;
  std::vector<std::vector<double> > A;

  int nvar, nmed;
  int ref_1, ref_2;

  std::vector<double> last_state;
  std::vector<double> curr_state;

  double max_result_value = 0.0;

  void initializeStructures(input_data networkData) 
  {
    nvar = calcula_nvar(networkData.grafo, networkData.numeroBarras);
    nmed = calcula_nmed(networkData.numeroMedidas);

    z.resize(nmed);
    b.resize(nmed);
    x.resize(nvar);
    regua.resize(nvar);
    H.resize(nmed, std::vector<double>(nvar));
    
    
    
    monta_z(z, networkData.numeroMedidas, networkData.medidas_pu);
    monta_regua(networkData.grafo, networkData.numeroBarras, nvar, regua);
    tratamento_referencia(&ref_1, &ref_2, networkData.alimentador->noRaiz, nvar);
    tira_refs_regua(nvar, ref_1, ref_2, regua); 
    nvar = nvar - (ref_2 - ref_1 +1); 
    incializa_vetor_x(networkData.grafo, networkData.numeroBarras, networkData.alimentador, networkData.numeroAlimentadores,x,nvar);
  }

  void updateStructures(input_data networkData) 
  {
    atualiza_Rede(networkData.grafo, networkData.numeroBarras);
    atualiza_Modelo(networkData.grafo, networkData.numeroBarras, nmed, networkData.medidas_pu);
    atualizaMatrizW(networkData.medidas_pu, networkData.numeroMedidas);
    atualiza_H(networkData.grafo, networkData.numeroBarras, networkData.allRamos, networkData.medidas_pu, nmed);
  }

  void buildNormalEquation(input_data networkData) 
  {
    montaMatrizA(networkData.medidas_pu, nmed, nvar);
    montaVetorB(networkData.medidas_pu);


  }

  void solveNormalEquation(int order_method) {
    Solver.order_method = order_method;
    Solver.solve(H, b);
    Dx = Solver.output;
    for (int i = 0; i<nvar; i++){
        x[i] += Solver.output[i];
    }
  }

  void updateSolution(input_data networkData)
  {
    atualiza_estado(networkData.grafo, x, nvar);
  }

  CONV_STATUS verifyConvergence(std::vector<double> result, double tol) 
  {
    int i;
    for (i = 0; i < nvar; i++)
        {
            if (std::abs(result[i]) >= tol)
            {
                return NOT_CONVERGENCE;
            }
        }
    return CONVERGENCE;
  }

  void writeSystemToFile()
    {
        ou.writeMatrixToFile(H);
        ou.writeVectorToFile(b);
    }

  

  void monta_W(double **W, long int nmed, DMED *medidas);
  void monta_W_Ident(double **W, long int nmed, DMED *medidas);
  
  void incializa_vetor_x_leitura(GRAFO *grafo, long int numeroBarras,
                                 ALIMENTADOR *alimentadores,
                                 long int numeroAlimentadores, double *x,
                                 double *regua, long int nVariaveis);

  
  void atualiza_h(GRAFO *grafo, long int numeroBarras, long int nmed,
                  DMED *medidas);
  
  void atualiza_H_ret(GRAFO *grafo, long int numeroBarras, DRAM *ramos,
                      DMED *medidas, long int numeroMedidas);



private:
  sparseSystem Solver;
  electrical_utils eu;
  mat_utils matUtils;
  output_utils ou;

  void monta_z(std::vector<double>& z, int nmed, DMED* medida)
  {
    int i;
    
    for(i=0;i<nmed;i++){
        z[i] = medida[i].zmed;
    }  
  }

    void monta_regua(GRAFO *grafo, int numeroBarras, int nvar, std::vector<double>& regua)
    {
        int i, j;
        double aux;
        j=0;
        for(i=0;i<numeroBarras;i++){
            aux = (double) grafo[i].idNo;
            aux += 0.01;
            switch (grafo[i].fases){
                case 1:
                    regua[j] = aux;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 2:
                    regua[j] = aux + 0.1;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 3:
                    regua[j] = aux + 0.2;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 4:
                    regua[j] = aux;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    regua[j] = aux + 0.1;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 5:
                    regua[j] = aux;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    regua[j] = aux + 0.2;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 6:
                    regua[j] = aux+0.1;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    regua[j] = aux + 0.2;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;
                case 7:
                    regua[j] = aux;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    regua[j] = aux + 0.1;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    regua[j] = aux + 0.2;
                    regua[j + (int) nvar/2] = -regua[j];
                    j++;
                    break;    
            }
        }
    }

    int calcula_nvar(GRAFO *grafo, int numeroBarras)
    {
        int _nvar = 0;
        int i;
        for (i = 0; i < numeroBarras; i++){
            switch (grafo[i].fases){
                case 1:
                    _nvar +=2;
                    break;
                case 2:
                    _nvar +=2;
                    break;
                case 3:
                    _nvar +=2;
                    break;
                case 4:
                    _nvar +=4;
                    break;    
                case 5:
                    _nvar +=4;
                    break;    
                case 6:
                    _nvar +=4;
                    break;    
                case 7:
                    _nvar +=6;
                    break;    
            }
        }
        return (_nvar);   

    }

    int calcula_nmed(int numeroMedidas)
    {
        int _nmed = 0;
        _nmed = numeroMedidas;
        return (_nmed);
    }
    void tratamento_referencia(int *ref_1, int *ref_2, int noRaiz, int nVariaveis)
    { 
        int j, k, n_ref = 3;
        for(j=0; j<nVariaveis; j++){
            if (regua[j] < 0){
                k = (int)regua[j];
                k = -k;
                if (k == noRaiz){
                    ref_1[0] = j;
                    ref_2[0] = j + n_ref - 1;
                    break;
                }
            }
        }        
    }
    void tira_refs_regua(int m, int col1, int col2, std::vector<double>& regua)
    {
    int j,n_cols;
    n_cols = col2 - col1 + 1;
    for (j=0;j<m;j++)
        {
            if (j<col1){
                regua[j] = regua[j];
            }
            else if (j>col2){
                regua[(j-n_cols)] = regua[j];
            }
            
        }
    
    }

    void incializa_vetor_x(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, std::vector<double>& x, long int nVariaveis)
    { 
        int i, k, fase;
        int *visitado;
        visitado = new int[numeroBarras];
        std::complex<double> V0[3];
        
        //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
        for(i=0; i<numeroBarras; i++){ 
            visitado[i] = FALSE;
        }
        for(i=0; i<numeroAlimentadores; i++)
        {
            //Tensão Inicial da subestação
            V0[0] = grafo[alimentadores[i].noRaiz].barra->Vinicial[0];//1.0*(cos(0) + I*sin(0));
            V0[1] = grafo[alimentadores[i].noRaiz].barra->Vinicial[1];//1.0*(cos(-120*PI/180) + I*sin(-120*PI/180));
            V0[2] = grafo[alimentadores[i].noRaiz].barra->Vinicial[2];//1.0*(cos(120*PI/180) + I*sin(120*PI/180));
            
            FILABARRAS *barraAtual = &alimentadores[i].rnp[0];
            
            int de = barraAtual->idNo;
            grafo[de].V[0] = V0[0];
            grafo[de].V[1] = V0[1];
            grafo[de].V[2] = V0[2];

            while(barraAtual != NULL)
            {
                de = barraAtual->idNo;
                int n_adj = grafo[de].numeroAdjacentes;
                for(k=0;k< n_adj;k++){
                    int para = grafo[de].adjacentes[k].idNo;
                    if (visitado[para] == FALSE){ 
                        if (grafo[de].adjacentes[k].tipo == 1){ //Atualiza o V0 para trafo visto a ligação e tap
                            grafo[para].V[0] = grafo[de].V[0];
                            grafo[para].V[1] = grafo[de].V[1];
                            grafo[para].V[2] = grafo[de].V[2];
                            
                            if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2)){                                
                                grafo[para].V[0] = std::abs(grafo[de].V[0])*std::complex<double>(cos(-30.0*PI/180), sin(-30.0*PI/180));
                                grafo[para].V[1] = std::abs(grafo[de].V[1])*std::complex<double>(cos(-150.0*PI/180), sin(-150.0*PI/180));
                                grafo[para].V[2] = std::abs(grafo[de].V[2])*std::complex<double>(cos(90.0*PI/180), sin(90.0*PI/180));                                    
                            }
                            else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2)){                                
                                grafo[para].V[0] = std::abs(grafo[de].V[0])*std::complex<double>(cos(-30.0*PI/180), sin(-30.0*PI/180));
                                grafo[para].V[1] = std::abs(grafo[de].V[1])*std::complex<double>(cos(-150.0*PI/180), sin(-150.0*PI/180));
                                grafo[para].V[2] = std::abs(grafo[de].V[2])*std::complex<double>(cos(90.0*PI/180), sin(90.0*PI/180));                                   
                            }
                            else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1)){
                                if (grafo[de].adjacentes[k].ramo->k == de){
                                    grafo[para].V[0] = std::abs(grafo[de].V[0])*std::complex<double>(cos(-30.0*PI/180), sin(-30.0*PI/180));
                                    grafo[para].V[1] = std::abs(grafo[de].V[1])*std::complex<double>(cos(-150.0*PI/180), sin(-150.0*PI/180));
                                    grafo[para].V[2] = std::abs(grafo[de].V[2])*std::complex<double>(cos(90.0*PI/180), sin(90.0*PI/180)); 
                                }
                                else{
                                    grafo[para].V[0] = std::abs(grafo[de].V[0])*std::complex<double>(cos(0.0*PI/180), sin(0.0*PI/180));
                                    grafo[para].V[1] = std::abs(grafo[de].V[1])*std::complex<double>(cos(-120.0*PI/180), sin(-120.0*PI/180));
                                    grafo[para].V[2] = std::abs(grafo[de].V[2])*std::complex<double>(cos(120.0*PI/180), sin(120.0*PI/180));
                                }
                            }            
                        }
                        else if (grafo[de].adjacentes[k].tipo == 2){ //Para o caso de regulador de tensão
                            grafo[para].V[0] = grafo[de].V[0]*grafo[de].adjacentes[k].ramo->tap_pri[0]*grafo[de].adjacentes[k].ramo->tap_sec[0];
                            grafo[para].V[1] = grafo[de].V[1]*grafo[de].adjacentes[k].ramo->tap_pri[1]*grafo[de].adjacentes[k].ramo->tap_sec[1];
                            grafo[para].V[2] = grafo[de].V[2]*grafo[de].adjacentes[k].ramo->tap_pri[2]*grafo[de].adjacentes[k].ramo->tap_sec[2];
                        }
                        else{
                            grafo[para].V[0] = grafo[de].V[0];
                            grafo[para].V[1] = grafo[de].V[1];
                            grafo[para].V[2] = grafo[de].V[2];
                        }
                    }
                }
                visitado[de] = TRUE;
                barraAtual = barraAtual->prox;
            }        
        }
        //Montagem do vetor x e da régua
        for(i=0; i<nVariaveis; i++)
        {
            k = (int)regua[i];
            fase = (int) std::abs((regua[i] - k)*10);
            if (regua[i] > 0){ //magnitude de tensão
                x[i] = std::abs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude            
            }
            else{ //ângulo de tensão
                // fase = fase;
                k = -k;
                x[i] = std::arg(grafo[k].V[fase]);
            }
        }
    }
    void atualiza_H(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas)
    {
        int i,j,k, de, para,ramo,opt,fase;
        std::complex<double> *dSaux, *dIaux;
        
        dSaux = matUtils.c_vetAloca(3);
        dIaux = matUtils.c_vetAloca(3);
        //Atualiza a matriz H de acordo com a régua
        for(i=0;i<numeroMedidas;i++){
            switch (medidas[i].tipo){
                case 0: //0: Fluxo de Potência Ativa - kW
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        de = medidas[i].k;
                        para = medidas[i].m;
                        ramo = medidas[i].ramo;
                        
                        if (ramos[ramo].k == de){
                            if(medidas[i].reguaH[j]>=0){
                                if (de == k) opt = 0;
                                else opt = 2;
                            }
                            else{
                                if (de == -k) opt = 1;
                                else opt = 3;
                            }
                            eu.dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                            medidas[i].H[j] = dSaux[medidas[i].fases - 1].real();
                        }
                        else{
                            if(medidas[i].reguaH[j]>=0){
                                if (de == k) opt = 2;
                                else opt = 0;
                            }
                            else{
                                if (de == -k) opt = 3;
                                else opt = 1;
                            }
                            eu.dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                            medidas[i].H[j] = dSaux[medidas[i].fases - 1].real();
                        }
                    }
                    break;
                case 1: //1: Fluxo de Potência Reativa - kVAr
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        de = medidas[i].k;
                        para = medidas[i].m;
                        ramo = medidas[i].ramo;
                        
                        if (ramos[ramo].k == de){
                            if(medidas[i].reguaH[j]>=0){
                                if (de == k) opt = 0;
                                else opt = 2;
                            }
                            else{
                                if (de == -k) opt = 1;
                                else opt = 3;
                            }
                            eu.dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                            medidas[i].H[j] = dSaux[medidas[i].fases - 1].imag();
                        }
                        else{
                            if(medidas[i].reguaH[j]>=0){
                                if (de == k) opt = 2;
                                else opt = 0;
                            }
                            else{
                                if (de == -k) opt = 3;
                                else opt = 1;
                            }
                            eu.dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                            medidas[i].H[j] = dSaux[medidas[i].fases - 1].imag();
                        }
                    }
                    break;    
                case 2: //2: Injeção de Potência Ativa - kW
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        de = medidas[i].k;
                        
                        if(medidas[i].reguaH[j]>=0){
                            opt = 0;
                        }
                        else{
                            opt = 1;
                        }
                        eu.dSk(grafo, de, dSaux, opt, std::abs(k), fase);
                    
                        medidas[i].H[j] = dSaux[medidas[i].fases - 1].real();
                    }
                    break;     
                case 3: //3: Injeção de Potência Reativa - kVAr
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        de = medidas[i].k;
                        
                        if(medidas[i].reguaH[j]>=0){
                            opt = 0;
                        }
                        else{
                            opt = 1;
                        }
                        eu.dSk(grafo, de, dSaux, opt, std::abs(k), fase);
                    
                        medidas[i].H[j] = dSaux[medidas[i].fases - 1].imag();
                    }break;
                case 4: //4: Magnitude de tensão
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]>0))
                            medidas[i].H[j] = 1;
                    }
                    break;
                case 5: //5: Ângulo de tensão
                    for(j=0;j<medidas[i].nvar;j++){
                        k = (int) medidas[i].reguaH[j];
                        fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                        if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]<0))
                            medidas[i].H[j] = 1;
                    }
                    break;
                case 7: //7: Magnitude de Corrente
                    // for(j=0;j<medidas[i].nvar;j++){
                    //     k = (int) medidas[i].reguaH[j];
                    //     fase = (int) std::abs((medidas[i].reguaH[j] - k)*10);
                        
                    //     de = medidas[i].k;
                    //     para = medidas[i].m;
                    //     ramo = medidas[i].ramo;
                        
                    //     if (ramos[ramo].k == de){
                    //         if(medidas[i].reguaH[j]>=0){
                    //             if (de == k) opt = 0;
                    //             else opt = 2;
                    //         }
                    //         else{
                    //             if (de == -k) opt = 1;
                    //             else opt = 3;
                    //         }
                    //         if (it == 0) eu.dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    //         else eu.dIkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    //         medidas[i].H[j] = dSaux[medidas[i].fases - 1];
                    //     }
                    //     else{
                    //         if(medidas[i].reguaH[j]>=0){
                    //             if (de == k) opt = 2;
                    //             else opt = 0;
                    //         }
                    //         else{
                    //             if (de == -k) opt = 3;
                    //             else opt = 1;
                    //         }
                    //         eu.dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    //         medidas[i].H[j] = dSaux[medidas[i].fases - 1].real();
                    //     }
                    // }
                    // break;      
                    break;     
            }        
        }
        free(dSaux);free(dIaux);
    }

    void atualiza_Modelo(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas)
    {
        int i,k,idMed,fase;
        
        for(idMed=0;idMed<nmed;idMed++){
            if (medidas[idMed].PARA == -1){ //Medidor instalado em uma barra
                i = medidas[idMed].k;
                fase = medidas[idMed].fases - 1;//Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo

                switch (medidas[idMed].tipo){
                    case 2: //2: Injeção de Potência Ativa - kW
                        medidas[idMed].h = grafo[i].S[fase].real();
                        break;
                    case 3: //3: Injeção de Potência Reativa - kVAr
                        medidas[idMed].h = grafo[i].S[fase].imag();
                        break;    
                    case 4: //4: Magnitude de Tensão - kV
                        medidas[idMed].h = std::abs(grafo[i].V[fase]);
                        break;     
                    case 5: //5: Ângulo de Tensão - graus
                        medidas[idMed].h = std::arg(grafo[i].V[fase]);
                        break;     
                    case 8: //8: Injeção Magnitude de Corrente - A  //REVISAR A MEDIDA DE INJEÇÂO DE PMU NO ESTIMADOR HIBRIDO
                        medidas[idMed].h = std::abs(grafo[i].Cur[fase]);
                        break;     
                    case 9: //9: Injeção Ângulo de Corrente) - graus
                        medidas[idMed].h = std::arg(grafo[i].Cur[fase]);
                        break; 
                }
            }
            else{ //Medidor instalado em um ramo
                i = medidas[idMed].k;
                fase = medidas[idMed].fases - 1;//Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo
                for(k=0;k<grafo[i].numeroAdjacentes;k++){
                    if (grafo[i].adjacentes[k].idNo == medidas[idMed].m){
                        switch (medidas[idMed].tipo){
                            case 0: //0: Fluxo de Potência Ativa - kW
                                medidas[idMed].h = grafo[i].adjacentes[k].S[fase].real();
                                break;
                            case 1: //1: Fluxo de Potência Reativa - kVAr
                                medidas[idMed].h = grafo[i].adjacentes[k].S[fase].imag();
                                break;    
                            case 6: //6: Fluxo Magnitude de Corrente - A
                                medidas[idMed].h = std::abs(grafo[i].adjacentes[k].Cur[fase]);
                                break;     
                            case 7: //7: Fluxo Ângulo de Corrente) - graus
                                medidas[idMed].h = std::arg(grafo[i].adjacentes[k].Cur[fase]);
                                break; 
                        }
                    }
                }
            }
        }
    }

    void atualiza_Rede(GRAFO *grafo, long int numeroBarras)
    {
        int i,k;
        std::complex<double> *Saux, *Iaux;

        Saux = matUtils.c_vetAloca(3);

        Iaux = matUtils.c_vetAloca(3);;    
        //Percorre o grafo atualizando o cálculo de h(x))
        for(i=0;i<numeroBarras;i++){        
            eu.Sk(grafo, i, Saux);
            grafo[i].S[0] = Saux[0];
            grafo[i].S[1] = Saux[1];
            grafo[i].S[2] = Saux[2];

            grafo[i].Cur[0] = conj(Saux[0]/grafo[i].V[0]);
            grafo[i].Cur[1] = conj(Saux[1]/grafo[i].V[1]);
            grafo[i].Cur[2] = conj(Saux[2]/grafo[i].V[2]);
            
            //Percorre os ramos adjacentes
            for(k=0;k<grafo[i].numeroAdjacentes;k++){
                if (i == grafo[i].adjacentes[k].ramo->k){
                    eu.Skm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Saux);
                    eu.Ikm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Iaux);
                }
                else{
                    eu.Smk(&grafo[grafo[i].adjacentes[k].idNo],&grafo[i], grafo[i].adjacentes[k].ramo, Saux);
                    eu.Imk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Iaux);
                }
                grafo[i].adjacentes[k].S[0] = Saux[0];
                grafo[i].adjacentes[k].S[1] = Saux[1];
                grafo[i].adjacentes[k].S[2] = Saux[2];

                grafo[i].adjacentes[k].Cur[0] = Iaux[0];
                grafo[i].adjacentes[k].Cur[1] = Iaux[1];
                grafo[i].adjacentes[k].Cur[2] = Iaux[2];
            }        
        }
        // free(Saux);free(Iaux);
    }

    void montaMatrizA(DMED *medidas, int nmed, int nvar)
    {
        int i, j, r;
        for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            medidas[i].H[j] = medidas[i].H[j] / medidas[i].sigma;
            for(r = 0;r<nvar;r++){
                if (std::abs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    H[i][r] = medidas[i].H[j];
                    break;
                }
            }
        }
     }  
    }

    void montaVetorB(DMED *medidas)
    {
        int i;
        for (i = 0; i < nmed; i++)
        {
            b[i] = (medidas[i].zmed - medidas[i].h) / medidas[i].sigma;
        }
    }

    void atualiza_estado(GRAFO *grafo, std::vector<double> x, int nVariaveis)
    {
        int i,k,fase;
        for(i=0; i<nVariaveis; i++)
        {
            k = (int)regua[i];
            fase = (int) std::abs((regua[i] - k)*10);
            if (regua[i] > 0){ //magnitude de tensão
                grafo[k].V[fase] = x[i]*grafo[k].V[fase]/std::abs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude            
            }
            else{ //ângulo de tensão
                fase = fase;
                k = -k;
                grafo[k].V[fase] = std::abs(grafo[k].V[fase])*std::complex<double>(cos(x[i]), sin(x[i])); //mantém a magnitude e altera o ângulo
            }     
        }
    }

    void atualizaMatrizW(DMED *medidas, int nmed)
    {
        int i,j;
        double prec, menorSigma = 1000000;
        double fundoEscala, auxFE;
        
        //Matriz W diagonal - inverso da variância
        for(i=0;i<nmed;i++){
            auxFE = round(medidas[i].zmed * 1.25 * 1000);
            fundoEscala = auxFE/1000;
            switch (medidas[i].tipo){
                case 0:
                case 1:
    //                medidas[i].sigma = 0.002;
    //                break;
                case 2:    
                case 3:
                    //W[i][i] = 40000;
    //                medidas[i].sigma = 0.010;
                    prec = medidas[i].prec; //0.02; //5% para SCADA de potência
                    medidas[i].sigma = 0.33333*prec*std::abs(medidas[i].zmed);
    //                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
                    break;
                case 4:
                case 6:
                    //W[i][i] = 40000;
    //                medidas[i].sigma = 0.001;
                    prec =  medidas[i].prec; //0.01; //1% para magnitude de tensão ou corrente
                    medidas[i].sigma = 0.33333*prec*std::abs(medidas[i].zmed);
    //                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
                    break;
                case 5:
                case 7:
                    //W[i][i] = 1000000;
                    medidas[i].sigma = 0.0001;
                    prec = medidas[i].prec; //0.001; //0.1% para PMU de ângulo de tensão ou corrente
                    break;
                case 8:
                case 9:
                case 10:
                case 11:    
                case 12:
                case 13:
                    //W[i][i] = 1000000;
    //                medidas[i].sigma = 0.0001;
                    prec = medidas[i].prec; //0.001; //0.1% para PMUs retangulares
                    medidas[i].sigma = 0.33333*prec*std::abs(medidas[i].zmed);
                    break;
            }
            //medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed); //Ponderação de acordo com o valor medido (Fórmula B. Pal)
            if (std::abs(medidas[i].zmed) > 0.00001){
                if (medidas[i].sigma < menorSigma) menorSigma = medidas[i].sigma;
            }
            medidas[i].sigma = medidas[i].sigma; 
    //        W[i][i] = 1/(pow(medidas[i].sigma,2)); 
        }
        for(i=0;i<nmed;i++){ //Tratamento da medida virtual e medidas proximo de zero
            if (std::abs(medidas[i].zmed) < 0.00001){
                medidas[i].sigma = menorSigma; 
                //medidas[i].sigma = 0.01*menorSigma; 
    //            W[i][i] = 1/(pow(medidas[i].sigma,2));
            }
        }
    }
                       

};
#endif