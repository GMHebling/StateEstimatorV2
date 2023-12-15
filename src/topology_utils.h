#ifndef TOPOLOGY_UTILS_H
#define TOPOLOGY_UTILS_H

#include "data_structures.h"
#include "mat_utils.h"
#include <math.h>


class topology_utils {
    public:
        ALIMENTADOR* feeder;

        void initializeTopology(GRAFO *grafo, int numeroBarras, int numeroAlimentadores){
            feeder = buscaProfundidadeAlimentadores(grafo, numeroBarras, numeroAlimentadores);
        }

        void buildNetworkModel(GRAFO *grafo, int numeroBarras, DRAM *allRamos, int numeroRamos, double Sbase)
        {
            calculaPU(grafo, numeroBarras, allRamos,  numeroRamos, Sbase);
            atualizaTaps(allRamos, numeroRamos);
        }


    private:
        mat_utils matUtils;

        ALIMENTADOR* buscaProfundidadeAlimentadores(GRAFO *grafo, int numeroBarras, int numeroAlimentadores) 
        {
            int i, idAlim = 0;
            FILABARRAS *barraAtual = NULL;
            
            int visitado[numeroBarras];
            ALIMENTADOR* m_feeder = new ALIMENTADOR[numeroAlimentadores];
            
            // if (((*alimentadores)= (ALIMENTADOR *)malloc( numeroAlimentadores * sizeof(ALIMENTADOR)))==NULL)
            // {
            //     printf("Erro -- Nao foi possivel alocar espaco de memoria para alimentadores !!!!");
            //     exit(1); 
            // }
            for(i=0; i<numeroBarras; i++){ 
                visitado[i] = 0;
                if(grafo[i].tipo == 2){
                    (m_feeder)[idAlim].idAlim = idAlim;
                    (m_feeder)[idAlim].noRaiz = i;
                    (m_feeder)[idAlim].numeroNos = 1;
                    (m_feeder)[idAlim].idRaiz = grafo[i].barra->ID;
                    (m_feeder)[idAlim].rnp[0].idNo = i;
                    (m_feeder)[idAlim].rnp[0].profundidade = 0;
                    (m_feeder)[idAlim].rnp[0].prox = NULL;
                    idAlim++;
                }
            }
            
            for(i=0; i<numeroAlimentadores; i++)
            {
                //buscaLargura(grafo, (m_feeder), i, (m_feeder)[i].noRaiz, visitado);
                barraAtual = &(m_feeder)[i].rnp[0];
                buscaProfundidade(barraAtual,(m_feeder)[i].noRaiz,0,visitado,grafo,i);
                //printf("\n alimentador %d Raiz: %d   Nos: %d",i,(m_feeder)[i].idRaiz,(*alimentadores)[i].numeroNos);
            }
            return (m_feeder);
            
        }
        void buscaProfundidade(FILABARRAS *barraAtual, long int idNo, int profundidade,  int *visitado, GRAFO * grafo, long int idAlim)
        {
            //Depth-Search Algorithm - busca no e a sua profundidade (gera RNP))
            long int barraAdj,i = 0;
            
            visitado[idNo] = 1;
            barraAtual->profundidade = profundidade;
            GRAFO * no = &grafo[idNo];
            grafo[idNo].idAlim = idAlim;
            grafo[idNo].profundidade = profundidade;
            //printf("\nidNo: %d  -  %d", grafo[idNo].barra->ID, profundidade);
            profundidade++;
            for(i = 0; i < no->numeroAdjacentes; i++)
            {   
                barraAdj = no->adjacentes[i].idNo;
                if (visitado[barraAdj]== 0)
                    {
                        idNo= barraAdj;
                        //adicionaNo(&barraAtual, idNo);
                        adicionaNoNaFila(&barraAtual, idNo);
                        apontaProxNoNaFila(&barraAtual);
                        buscaProfundidade(barraAtual, idNo, profundidade, visitado, grafo, idAlim);
                    } 
            }
        }
        void adicionaNoNaFila(FILABARRAS ** fila, long int idNo) 
        {
            FILABARRAS *novoVertice = NULL;
            FILABARRAS *aux = NULL;
        
            novoVertice = (FILABARRAS *)malloc(sizeof(FILABARRAS));
        
            if(novoVertice == NULL) {
                printf("erro insere_fila\n");
                exit(EXIT_FAILURE);
            }
            
            novoVertice->idNo = idNo;
            novoVertice->prox = NULL;
            
            if(*fila == NULL)
                *fila = novoVertice;
            else {
                aux = *fila;
                while(aux->prox !=NULL) aux = aux->prox;
                aux->prox = novoVertice;
            }
        }
        void apontaProxNoNaFila(FILABARRAS ** fila) 
        {
            FILABARRAS *aux = NULL;
            
            aux = *fila;
            while(aux->prox !=NULL) aux = aux->prox;
            *fila = aux;
            
        }
        void atualizaTaps(DRAM *allRamos, int numeroRamos)
        {
            int i;        
            for(i=0;i<numeroRamos;i++){
                if (allRamos[i].tipo == 2){
                    allRamos[i].tap_pri[0] = (1 + allRamos[i].regulador.tap[0]*allRamos[i].regulador.regulacao/allRamos[i].regulador.ntaps);
                    allRamos[i].tap_pri[1] = (1 + allRamos[i].regulador.tap[1]*allRamos[i].regulador.regulacao/allRamos[i].regulador.ntaps);
                    allRamos[i].tap_pri[2] = (1 + allRamos[i].regulador.tap[2]*allRamos[i].regulador.regulacao/allRamos[i].regulador.ntaps);
                    
                    allRamos[i].tap_sec[0] = 1;
                    allRamos[i].tap_sec[1] = 1;
                    allRamos[i].tap_sec[2] = 1;
                    
                    allRamos[i].Ypp[0][0] = pow(allRamos[i].tap_pri[0],2)*allRamos[i].Ypp[0][0];
                    allRamos[i].Ypp[1][1] = pow(allRamos[i].tap_pri[1],2)*allRamos[i].Ypp[1][1];
                    allRamos[i].Ypp[2][2] = pow(allRamos[i].tap_pri[2],2)*allRamos[i].Ypp[2][2];
                    
                    allRamos[i].Yps[0][0] = allRamos[i].tap_sec[0]*allRamos[i].tap_pri[0]*allRamos[i].Yps[0][0];
                    allRamos[i].Yps[1][1] = allRamos[i].tap_sec[1]*allRamos[i].tap_pri[1]*allRamos[i].Yps[1][1];
                    allRamos[i].Yps[2][2] = allRamos[i].tap_sec[2]*allRamos[i].tap_pri[2]*allRamos[i].Yps[2][2];
                    
                    allRamos[i].Ysp[0][0] = allRamos[i].tap_sec[0]*allRamos[i].tap_pri[0]*allRamos[i].Ysp[0][0];
                    allRamos[i].Ysp[1][1] = allRamos[i].tap_sec[1]*allRamos[i].tap_pri[1]*allRamos[i].Ysp[1][1];
                    allRamos[i].Ysp[2][2] = allRamos[i].tap_sec[2]*allRamos[i].tap_pri[2]*allRamos[i].Ysp[2][2];
                    
                    allRamos[i].Yss[0][0] = pow(allRamos[i].tap_sec[0],2)*allRamos[i].Yss[0][0];
                    allRamos[i].Yss[1][1] = pow(allRamos[i].tap_sec[1],2)*allRamos[i].Yss[1][1];
                    allRamos[i].Yss[2][2] = pow(allRamos[i].tap_sec[2],2)*allRamos[i].Yss[2][2];
                }
                else if (allRamos[i].tipo == 1){ //REVISAR PARA TAP FORA DA NOMINAL EM TRAFOS
                    allRamos[i].tap_pri[0] = 1;
                    allRamos[i].tap_pri[1] = 1;
                    allRamos[i].tap_pri[2] = 1;
                    
                    allRamos[i].tap_sec[0] = 1;
                    allRamos[i].tap_sec[1] = 1;
                    allRamos[i].tap_sec[2] = 1;
                }
            }
        }
        void calculaPU(GRAFO *grafo, int numeroBarras, DRAM *ramos, int numeroRamos, double Sbase) 
        {
            int i, idNo, idRam;
            double Vbase;
            
            //Transforma em PU informações dos nós do grafo
            for (idNo=0;idNo<numeroBarras;idNo++){
                //Transforma em PU os shunts da rede elétrica
                for(i=0;i<grafo[idNo].barra->nshunts;i++){
                    grafo[idNo].barra->shunts[i].Qnom[0] = 1000*grafo[idNo].barra->shunts[i].Qnom[0]/Sbase; 
                    grafo[idNo].barra->shunts[i].Qnom[1] = 1000*grafo[idNo].barra->shunts[i].Qnom[1]/Sbase;
                    grafo[idNo].barra->shunts[i].Qnom[2] = 1000*grafo[idNo].barra->shunts[i].Qnom[2]/Sbase;
                    montaQuadripoloShunt(&grafo[idNo],&grafo[idNo].barra->shunts[i]);
                }
            }
            
            //Transforma em PU informações dos ramos do grafo
            for (idRam=0;idRam<numeroRamos;idRam++){
                Vbase = grafo[ramos[idRam].m].Vbase;
                //Transforma as impedâncias em pu
                switch(ramos[idRam].tipo){
                    case 0:
                        ramos[idRam].linha.Zaa = ramos[idRam].linha.Zaa/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Zab = ramos[idRam].linha.Zab/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Zac = ramos[idRam].linha.Zac/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Zbb = ramos[idRam].linha.Zbb/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Zbc = ramos[idRam].linha.Zbc/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Zcc = ramos[idRam].linha.Zcc/((pow(Vbase,2))/Sbase);
                        
                        ramos[idRam].linha.Baa = ramos[idRam].linha.Baa/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Bab = ramos[idRam].linha.Bab/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Bac = ramos[idRam].linha.Bac/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Bbb = ramos[idRam].linha.Bbb/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Bbc = ramos[idRam].linha.Bbc/((pow(Vbase,2))/Sbase);
                        ramos[idRam].linha.Bcc = ramos[idRam].linha.Bcc/((pow(Vbase,2))/Sbase);
                        
                        
                        montaQuadripoloLinha(&ramos[idRam], &ramos[idRam].linha);
                        break;
                    case 1:
                        ramos[idRam].trafo.R = 3.0*ramos[idRam].trafo.R/((pow(Vbase,2))/Sbase);
                        ramos[idRam].trafo.X = 3.0*ramos[idRam].trafo.X/((pow(Vbase,2))/Sbase);
                        
                        montaQuadripoloTrafo(&ramos[idRam], &ramos[idRam].trafo);
                        break;
                    case 2:
                        ramos[idRam].regulador.R = 3.0*ramos[idRam].regulador.R/((pow(Vbase,2))/Sbase);
                        ramos[idRam].regulador.X = 3.0*ramos[idRam].regulador.X/((pow(Vbase,2))/Sbase);
                        
                        montaQuadripoloRegulador(&ramos[idRam], &ramos[idRam].regulador);
                        break;    
                }
            }   
        }
        void montaQuadripoloShunt(GRAFO *no, DSHNT *shunt)
        {
            int i, j;
            std::complex<double> ya,yb,yc, **Yi, **Yii, **Yiii, **Ysh;
            
            //Aloca Matrizes de Quadripolos
            Ysh = matUtils.c_matAloca(3);
            
            
            //Aloca Matrizes de Ligação D-Y-YN
            Yi = matUtils.c_matAloca(3);
            Yii = matUtils.c_matAloca(3);
            Yiii = matUtils.c_matAloca(3);
            
            //Admitância do trafo
            
            ya = std::complex<double> (0, shunt->Qnom[0]/(pow(shunt->Vbase/sqrt(3.0),2)/pow(no->Vbase,2)));
            yb = std::complex<double> (0, shunt->Qnom[1]/(pow(shunt->Vbase/sqrt(3.0),2)/pow(no->Vbase,2)));
            yc = std::complex<double> (0, shunt->Qnom[2]/(pow(shunt->Vbase/sqrt(3.0),2)/pow(no->Vbase,2)));
            
            //Yi
            Yi[0][0] = std::complex<double> (1, 0)*ya;
            Yi[0][1] = 0;
            Yi[0][2] = 0;
            Yi[1][0] = 0;
            Yi[1][1] = std::complex<double> (1, 0)*yb;
            Yi[1][2] = 0;
            Yi[2][0] = 0;
            Yi[2][1] = 0;
            Yi[2][2] = std::complex<double> (1, 0)*yc;
            //Yii
            Yii[0][0] = std::complex<double> (2, 0)*ya;
            Yii[0][1] = std::complex<double> (-1, 0)*yb;
            Yii[0][2] = std::complex<double> (-1, 0)*yc;
            Yii[1][0] = std::complex<double> (-1, 0)*ya;
            Yii[1][1] = std::complex<double> (2, 0)*yb;
            Yii[1][2] = std::complex<double> (-1, 0)*yc;
            Yii[2][0] = std::complex<double> (-1, 0)*ya;
            Yii[2][1] = std::complex<double> (-1, 0)*yb;
            Yii[2][2] = std::complex<double> (2, 0)*yc;
            //Yiii
            Yiii[0][0] = std::complex<double> (-1/(pow(3.0,0.5)), 0)*ya;
            Yiii[0][1] = std::complex<double> (1/(pow(3.0,0.5)), 0)*yb;
            Yiii[0][2] = 0;
            Yiii[1][0] = 0;
            Yiii[1][1] = std::complex<double> (-1/(pow(3.0,0.5)), 0)*yb;
            Yiii[1][2] = std::complex<double> (1/(pow(3.0,0.5)), 0)*yc;
            Yiii[2][0] = std::complex<double> (1/(pow(3.0,0.5)), 0)*ya;
            Yiii[2][1] = 0;
            Yiii[2][2] = std::complex<double> (-1/(pow(3.0,0.5)), 0)*yc;
            
            if(shunt->lig == 1){ //Ligação YN
                Ysh = matUtils.c_matIgual(Yi, 3);
            }
            else if(shunt->lig == 2){ //Ligação D trifásico
                Ysh = matUtils.c_matIgual(Yiii, 3);
                Ysh[2][0] = Ysh[2][0] + std::complex<double> (1.0, 0); //Sequência zero para manter o quadripólo com posto completo
                Ysh[2][1] = Ysh[2][1] + std::complex<double> (1.0, 0);
                Ysh[2][2] = Ysh[2][2] + std::complex<double> (1.0, 0);
            }
            else if(shunt->lig == 3){ //Ligação Y
                Ysh = matUtils.c_matIgual(Yii, 3);
            }
            for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                    no->Ysh[i][j] += Ysh[i][j];
                }
            }
        }
        void montaQuadripoloLinha(DRAM *ramo, DLIN *linha)
        {
            std::complex<double> **Zl,**B;
            
            //Aloca Matrizes de Quadripolos
            ramo->Ypp = matUtils.c_matAloca(3);
            ramo->Yps = matUtils.c_matAloca(3);
            ramo->Ysp = matUtils.c_matAloca(3);
            ramo->Yss = matUtils.c_matAloca(3);
            
            //Aloca Matrizes de Impedância e Admitância
            Zl = matUtils.c_matAloca(3);
            B = matUtils.c_matAloca(3);
            
            //Matriz Impedância da linha
            Zl[0][0] = linha->Zaa;
            Zl[0][1] = linha->Zab;
            Zl[0][2] = linha->Zac;
            Zl[1][0] = linha->Zab;
            Zl[1][1] = linha->Zbb;
            Zl[1][2] = linha->Zbc;
            Zl[2][0] = linha->Zac;
            Zl[2][1] = linha->Zbc;
            Zl[2][2] = linha->Zcc;
            
            //Matriz Susceptãncia Shunt da linha
            B[0][0] = std::complex<double>(0,1) * linha->Baa/std::complex<double>(2,0);
            B[0][1] = std::complex<double>(0,1) * linha->Bab/std::complex<double>(2,0);
            B[0][2] = std::complex<double>(0,1) * linha->Bac/std::complex<double>(2,0);
            B[1][0] = std::complex<double>(0,1) * linha->Bab/std::complex<double>(2,0);
            B[1][1] = std::complex<double>(0,1) * linha->Bbb/std::complex<double>(2,0);
            B[1][2] = std::complex<double>(0,1) * linha->Bbc/std::complex<double>(2,0);
            B[2][0] = std::complex<double>(0,1) * linha->Bac/std::complex<double>(2,0);
            B[2][1] = std::complex<double>(0,1) * linha->Bbc/std::complex<double>(2,0);
            B[2][2] = std::complex<double>(0,1) * linha->Bcc/std::complex<double>(2,0);
            
            //Inversa de Z - salva na variável Zl
            Zl = matUtils.c_matInversaZ(Zl, 3);
            
            ramo->Ypp = matUtils.c_matIgual(Zl, 3);
            ramo->Yss = matUtils.c_matIgual(Zl, 3);
            Zl = matUtils.c_matMultEsc(Zl, -1, 3);
            ramo->Yps = matUtils.c_matIgual(Zl, 3);
            ramo->Ysp = matUtils.c_matIgual(Zl, 3);
            
            
            ramo->Ypp[0][0] = ramo->Ypp[0][0] + B[0][0];
            ramo->Ypp[0][1] = ramo->Ypp[0][1] + B[0][1];
            ramo->Ypp[0][2] = ramo->Ypp[0][2] + B[0][2];
            ramo->Ypp[1][0] = ramo->Ypp[1][0] + B[1][0];
            ramo->Ypp[1][1] = ramo->Ypp[1][1] + B[1][1];
            ramo->Ypp[1][2] = ramo->Ypp[1][2] + B[1][2];
            ramo->Ypp[2][0] = ramo->Ypp[2][0] + B[2][0];
            ramo->Ypp[2][1] = ramo->Ypp[2][1] + B[2][1];
            ramo->Ypp[2][2] = ramo->Ypp[2][2] + B[2][2];
            
            ramo->Yss[0][0] = ramo->Yss[0][0] + B[0][0];
            ramo->Yss[0][1] = ramo->Yss[0][1] + B[0][1];
            ramo->Yss[0][2] = ramo->Yss[0][2] + B[0][2];
            ramo->Yss[1][0] = ramo->Yss[1][0] + B[1][0];
            ramo->Yss[1][1] = ramo->Yss[1][1] + B[1][1];
            ramo->Yss[1][2] = ramo->Yss[1][2] + B[1][2];
            ramo->Yss[2][0] = ramo->Yss[2][0] + B[2][0];
            ramo->Yss[2][1] = ramo->Yss[2][1] + B[2][1];
            ramo->Yss[2][2] = ramo->Yss[2][2] + B[2][2];
        }
        void montaQuadripoloTrafo(DRAM *ramo, DTRF *trafo)
        {
            std::complex<double> y, **Yi, **Yii,**Yiii;
            
            //Aloca Matrizes de Quadripolos
            ramo->Ypp = matUtils.c_matAloca(3);
            ramo->Yps = matUtils.c_matAloca(3);
            ramo->Ysp = matUtils.c_matAloca(3);
            ramo->Yss = matUtils.c_matAloca(3);
            
            //Aloca Matrizes de Ligação D-Y-YN
            Yi = matUtils.c_matAloca(3);
            Yii = matUtils.c_matAloca(3);
            Yiii = matUtils.c_matAloca(3);
            
            //Admitância do trafo
            y = std::complex<double>(1,0)/std::complex<double>(trafo->R, trafo->X);
            
            //Yi
            Yi[0][0] = 1.0;
            Yi[0][1] = 0;
            Yi[0][2] = 0;
            Yi[1][0] = 0;
            Yi[1][1] = 1.0;
            Yi[1][2] = 0;
            Yi[2][0] = 0;
            Yi[2][1] = 0;
            Yi[2][2] = 1.0;
            //Yii
            Yii[0][0] = 2.0;
            Yii[0][1] = -1.0;
            Yii[0][2] = -1.0;
            Yii[1][0] = -1.0;
            Yii[1][1] = 2.0;
            Yii[1][2] = -1.0;
            Yii[2][0] = -1.0;
            Yii[2][1] = -1.0;
            Yii[2][2] = 2.0;
            //Yiii
            Yiii[0][0] = -1.0/(pow(3.0,0.5));
            Yiii[0][1] = 1.0/(pow(3.0,0.5));
            Yiii[0][2] = 0;
            Yiii[1][0] = 0;
            Yiii[1][1] = -1.0/(pow(3.0,0.5));
            Yiii[1][2] = 1.0/(pow(3.0,0.5));
            Yiii[2][0] = 1.0/(pow(3.0,0.5));
            Yiii[2][1] = 0;
            Yiii[2][2] = -1.0/(pow(3.0,0.5));
            
            Yi = matUtils.c_matMultEsc(Yi, y, 3);
            Yii = matUtils.c_matMultEsc(Yii, y/std::complex<double>(3,0), 3);
            Yiii = matUtils.c_matMultEsc(Yiii, y, 3);
                
            if((trafo->lig_pri == 1)&& (trafo->lig_sec == 1)){ //Ligação YN-YN
                ramo->Ypp = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1, 3);
                ramo->Ysp = matUtils.c_matIgual(Yi, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1, 3);
                ramo->Yss = matUtils.c_matIgual(Yi, 3);
            }
            else if((trafo->lig_pri == 1)&& (trafo->lig_sec == 3)){ //Ligação YN-Y
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*std::abs(y);
                ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 1)&& (trafo->lig_sec == 2)){ //Ligação YN-D
                ramo->Ypp = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matTransp(ramo->Ysp, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Yss[2][0] = ramo->Yss[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yss[2][1] = ramo->Yss[2][1] + 100*std::abs(y);
                ramo->Yss[2][2] = ramo->Yss[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 1)){ //Ligação Y-YN
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Yss[2][0] = ramo->Yss[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yss[2][1] = ramo->Yss[2][1] + 100*std::abs(y);
                ramo->Yss[2][2] = ramo->Yss[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 3)){ //Ligação Y-Y
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Yss[2][0] = ramo->Yss[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yss[2][1] = ramo->Yss[2][1] + 100*std::abs(y); //Pensar se trata-se de simplificação, pois pode existir tensão de neutro neste caso, mas não irá te corrente de neutro.
                ramo->Yss[2][2] = ramo->Yss[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 2)){ //Ligação Y-D
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matTransp(ramo->Ysp, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Yss[2][0] = ramo->Yss[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yss[2][1] = ramo->Yss[2][1] + 100*std::abs(y);
                ramo->Yss[2][2] = ramo->Yss[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 1)){ //Ligação D-YN
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matTransp(ramo->Ysp, 3);
                ramo->Yss = matUtils.c_matIgual(Yi, 3);
                ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*std::abs(y);
                ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*std::abs(y);
                ramo->Yps[2][0] = ramo->Yps[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yps[2][1] = ramo->Yps[2][1] + 100*std::abs(y);
                ramo->Yps[2][2] = ramo->Yps[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 3)){ //Ligação D-Y
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yiii, 3);
                ramo->Yps = matUtils.c_matTransp(ramo->Yps, 3);
                ramo->Ysp = matUtils.c_matIgual(Yiii, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*std::abs(y);
                ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*std::abs(y);
                
                // Aproximação para o caso do Y não aterrado - tensão de sequência é igual a zero no secundário - fica igual DYn
                ramo->Yps[2][0] = ramo->Yps[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yps[2][1] = ramo->Yps[2][1] + 100*std::abs(y);
                ramo->Yps[2][2] = ramo->Yps[2][2] + 100*std::abs(y);
            }
            else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 2)){ //Ligação D-D
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual( Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
                ramo->Yss[2][0] = ramo->Yss[2][0] + 100*std::abs(y); //Sequência zero para manter o quadripólo com posto completo
                ramo->Yss[2][1] = ramo->Yss[2][1] + 100*std::abs(y);
                ramo->Yss[2][2] = ramo->Yss[2][2] + 100*std::abs(y);
            }
        }
        void montaQuadripoloRegulador(DRAM *ramo, DREG *reg)
        {
            std::complex<double> y, **Yi, **Yii,**Yiii;
            
            //Aloca Matrizes de Quadripolos
            ramo->Ypp = matUtils.c_matAloca(3);
            ramo->Yps = matUtils.c_matAloca(3);
            ramo->Ysp = matUtils.c_matAloca(3);
            ramo->Yss = matUtils.c_matAloca(3);
            
            //Aloca Matrizes de Ligação D-Y-YN
            Yi = matUtils.c_matAloca(3);
            Yii = matUtils.c_matAloca(3);
            Yiii = matUtils.c_matAloca(3);
            
            //Admitância do trafo
            y = std::complex<double>(1, 0)/std::complex<double>(reg->R, reg->X);
            
            //Yi
            Yi[0][0] = 1.0;
            Yi[0][1] = 0;
            Yi[0][2] = 0;
            Yi[1][0] = 0;
            Yi[1][1] = 1.0;
            Yi[1][2] = 0;
            Yi[2][0] = 0;
            Yi[2][1] = 0;
            Yi[2][2] = 1.0;
            //Yii
            Yii[0][0] = 2.0;
            Yii[0][1] = -1.0;
            Yii[0][2] = -1.0;
            Yii[1][0] = -1.0;
            Yii[1][1] = 2.0;
            Yii[1][2] = -1.0;
            Yii[2][0] = -1.0;
            Yii[2][1] = -1.0;
            Yii[2][2] = 2.0;
            //Yiii
            Yiii[0][0] = -1.0/(pow(3.0,0.5));
            Yiii[0][1] = 1.0/(pow(3.0,0.5));
            Yiii[0][2] = 0;
            Yiii[1][0] = 0;
            Yiii[1][1] = -1.0/(pow(3.0,0.5));
            Yiii[1][2] = 1.0/(pow(3.0,0.5));
            Yiii[2][0] = 1.0/(pow(3.0,0.5));
            Yiii[2][1] = 0;
            Yiii[2][2] = -1.0/(pow(3.0,0.5));
            
            
            //Feito somente para Yn
            switch (ramo->fases){
                case 1:
                    Yi[1][1] = 0;
                    Yi[2][2] = 0;
                    break;    
                case 2:
                    Yi[0][0] = 0;
                    Yi[2][2] = 0;
                    break;
                case 3:
                    Yi[1][1] = 0;
                    Yi[0][0] = 0;
                    break;
                case 4:
                    Yi[2][2] = 0;
                    break;    
                case 5:
                    Yi[1][1] = 0;
                    break;
                case 6:
                    Yi[0][0] = 0;
                    break; 
            }   
            
            Yi = matUtils.c_matMultEsc(Yi, y, 3);
            Yii = matUtils.c_matMultEsc(Yii, y/std::complex<double>(3,0), 3);
            Yiii = matUtils.c_matMultEsc(Yiii, y, 3);
            
            //Feito Somente para YN!!!!!!!!!!!!!!!!!!!!!!!!!
            if(reg->lig == 1){ //Ligação YN
                ramo->Ypp = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1.0, 3);
                ramo->Ysp = matUtils.c_matIgual(Yi, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1.0, 3);
                ramo->Yss = matUtils.c_matIgual(Yi, 3);
            }
            else if(reg->lig == 2){ //Ligação D
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1.0, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1.0, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
            }
            else if(reg->lig == 3){ //Ligação Y
                ramo->Ypp = matUtils.c_matIgual(Yi, 3);
                ramo->Yps = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matIgual(Yiii, 3);
                ramo->Ysp = matUtils.c_matTransp(ramo->Ysp, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
            }
            else if(reg->lig == 4){ //Ligação OYn
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1.0, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1.0, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
            }
            else if(reg->lig == 5){ //Ligação OD
                ramo->Ypp = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matIgual(Yii, 3);
                ramo->Yps = matUtils.c_matMultEsc(ramo->Yps, -1.0, 3);
                ramo->Ysp = matUtils.c_matIgual(Yii, 3);
                ramo->Ysp = matUtils.c_matMultEsc(ramo->Ysp, -1.0, 3);
                ramo->Yss = matUtils.c_matIgual(Yii, 3);
            }
        }
};
#endif