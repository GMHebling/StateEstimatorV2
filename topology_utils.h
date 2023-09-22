#ifndef TOPOLOGY_UTILS_H
#define TOPOLOGY_UTILS_H

class topology_utils {
    public:
        ALIMENTADOR* feeder;

        void initializeTopology(GRAFO *grafo, int numeroBarras, int numeroAlimentadores){
            feeder = buscaProfundidadeAlimentadores(grafo, numeroBarras, numeroAlimentadores);
        }

        void buildNetworkModel(DRAM *allRamos, int numeroRamos)
        {
            //atualizaTaps(allRamos, numeroRamos);
        }


    private:

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
        void adicionaNoNaFila(FILABARRAS ** fila, long int idNo) {
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
        void apontaProxNoNaFila(FILABARRAS ** fila) {
            FILABARRAS *aux = NULL;
            
            aux = *fila;
            while(aux->prox !=NULL) aux = aux->prox;
            *fila = aux;
            
        }

        void atualizaTaps(DRAM *allRamos, int numeroRamos){
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
};
#endif