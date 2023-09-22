#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "data_structures.h"
#include "topology_utils.h"

#include <string>
#include <iostream>
#include <complex>

class input_data {
    public:
        //Dados mestres
        DBAR *barras;
        DLIN *linhas;

        DRAM *ramos;
        DRAM *trafos;
        DRAM *reguladores;
        DRAM *chaves;

        //Juncao de todos os dados do tipo RAMO
        //linhas, trafos, reguladores e chaves
        DRAM *allRamos;

        DLOAD *load;
        DSHNT *shunt;
        DGD *dgd;
        // DTRF *trafo;
        DREG *reg;
        ALIMENTADOR *alimentador;
        //Estrutura de medidas
        DMED *medidas;
        //Estrutura de medidas apos processamento - conversao para PU
        DMED *medidas_pu;
        GRAFO *grafo;

        int numeroBarras;
        int numeroRamos;
        int numeroTrafos;
        int numeroRegs;
        int numeroChaves;
        int numeroAlimentadores;

        int numeroMedidas;

        int ierror;

        int i,j;

        double Sbase = 1000000; //VA

        std::string dataDir;
        std::string dataFolder;

        void initialize(){
            readConfigFile(dataDir, dataFolder);
            std::cout << dataDir << "\n";
            std::cout << dataFolder << "\n";

            barras = leituraDBAR(numeroBarras, numeroAlimentadores);
            ierror = leituraDSHNT(barras, numeroBarras);
            if (ierror) {
                std::cout << "DSHNT not found\n";
            }
            ierror = leituraDGD(barras, numeroBarras);
            if (ierror) {
                std::cout << "DGD not found\n";
            }
            ramos = leituraDLIN(numeroRamos, barras, numeroBarras);
            trafos = leituraDTRF(numeroTrafos, barras, numeroBarras);
            reguladores = leituraDREG(numeroRegs, barras, numeroBarras);
            chaves = leituraDSWTC(numeroChaves, barras, numeroBarras);

            allRamos = joinDRAM(ramos, trafos, reguladores, chaves);
            numeroRamos = numeroRamos + numeroTrafos + numeroRegs + numeroChaves;

            ierror = leituraVinicial(barras, numeroBarras);
            if (ierror) {
                std::cout << "Vinicial not found\n";
            }
            grafo = geraGrafo(barras, numeroBarras, allRamos, numeroRamos);

            medidas = leituraMedidas(allRamos, numeroRamos, barras, numeroBarras, grafo, Sbase, numeroMedidas);
            medidas_pu = processaMedidas(medidas, numeroMedidas, grafo, numeroBarras);

            topology.initializeTopology(grafo, numeroBarras, numeroAlimentadores);
            alimentador = topology.feeder;
            topology.buildNetworkModel(allRamos, numeroRamos);
            
        }


    private:
        topology_utils topology;

        void readConfigFile(std::string &dataDir, std::string &dataFolder)
        {
            char linha[1000], *pasta, *folder;
            std::cout << "Leitura de dados da rede eletrica...\n";
            
            
            //Recebe o nome da pasta com os dados a serem lidos - arquivo config.txt
            FILE *config = NULL;
            config = fopen("config.txt", "r");
            if(config == NULL){
                    std::cout << "Erro ao abrir arquivo config.txt !!!\n";
                    exit(1);
            };

            fgets(linha, 1000, config);
            folder = getField(linha,1);
            dataDir = folder;
            fgets(linha, 1000, config);
            pasta = getField(linha,1);
            dataFolder = pasta;
            fclose(config);
        }
        char* getField(char* lin, int num){
            char* tok;
            char *line;
            line = strdup(lin);
            //printf("\nteste -  %s",line);
            tok = strtok(line, ",\t\n\r");
            while (tok != NULL){
                if (!--num){
                    return tok;
                }
                tok = strtok (NULL, ",\t\n\r");
            }
            return NULL;
        }
        DBAR* leituraDBAR(int &numeroBarras, int &numeroAlimentadores)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DBAR.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "Erro no arquivo DBAR.csv\n";
                exit(1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int contador =0, i, aux, k; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            double PA,PB,PC,QA,QB,QC;
            dados = new char[100];
            //Aloca na memória espaço para as barras
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            DBAR* m_barras = new DBAR[numLinhas];
            //m_barras = (DBAR *)malloc((numLinhas+1)*sizeof(DBAR));
            rewind(arquivo);
            
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                
                dados = blocoLeitura;
                //printf("%s\n", blocoLeitura);
                
                //Verifica se a barra já foi criada
                aux = -1;
                for (i=0;i<contador;i++){
                    if (m_barras[i].ID == atoi(getField(dados,1))){
                        aux = i;
                    }
                }
                
                if (aux == -1){ //Criando novo DBAR
                    m_barras[contador].ID = atoi(getField(dados,1));
                    m_barras[contador].i = contador;
                    m_barras[contador].ligacao = static_cast<LIGACAO>(atoi(getField(dados,2)));
                    m_barras[contador].fases = static_cast<FASES>(atoi(getField(dados,3)));
                    m_barras[contador].Vbase = atof(getField(dados,4))/pow(3,0.5);
                    m_barras[contador].tipo = 0;
                    
                    //m_barras[contador].loads = (DLOAD *)malloc( 1 * sizeof(DLOAD));
                    m_barras[contador].nloads = 0;
                    //m_barras[contador].shunts = (DSHNT *)malloc( 1 * sizeof(DSHNT));
                    m_barras[contador].nshunts = 0;
                    //m_barras[contador].gds = (DGD *)malloc( 1 * sizeof(DGD));
                    m_barras[contador].ngds = 0;
                    
                    m_barras[contador].Vinicial[0] = std::complex<double>(1*(cos(0), sin(0)));
                    m_barras[contador].Vinicial[1] = std::complex<double>(1*(cos(-120*PI/180), sin(-120*PI/180)));
                    m_barras[contador].Vinicial[2] = std::complex<double>(1*(cos(120*PI/180), sin(120*PI/180)));
                    
                    
                    //Leitura das Cargas
                    PA = atof(getField(dados,5));
                    PB = atof(getField(dados,6));
                    PC = atof(getField(dados,7));
                    QA = atof(getField(dados,8));
                    QB = atof(getField(dados,9));
                    QC = atof(getField(dados,10));
                    
                    if ((PA != 0) || (PB != 0) || (PC != 0) || (QA != 0) || (QB != 0) || (QC != 0)){
                        m_barras[contador].nloads++;

                        k = m_barras[contador].nloads - 1;
                        
                        m_barras[contador].loads[k].ID = m_barras[contador].ID;
                        m_barras[contador].loads[k].Vbase = m_barras[contador].Vbase;
                        m_barras[contador].loads[k].fases = m_barras[contador].fases;
                        m_barras[contador].loads[k].lig = m_barras[contador].ligacao;
                        m_barras[contador].loads[k].ZIP = atof(getField(dados,11));
                        
                        m_barras[contador].loads[k].Pnom[0] = PA;
                        m_barras[contador].loads[k].Pnom[1] = PB;
                        m_barras[contador].loads[k].Pnom[2] = PC;
                        m_barras[contador].loads[k].Qnom[0] = QA;
                        m_barras[contador].loads[k].Qnom[1] = QB;
                        m_barras[contador].loads[k].Qnom[2] = QC;
                        
                    }
                    
                    //Leitura da Barra de Referência
                    if (getField(dados,13) != NULL){
                        m_barras[contador].tipo = 2;
                        numeroAlimentadores++;
                        
                        double VA = atof(getField(dados,12));
                        double VB = atof(getField(dados,13));
                        double VC = atof(getField(dados,14));
                        double TA = atof(getField(dados,15));
                        double TB = atof(getField(dados,16));
                        double TC = atof(getField(dados,17));
                        
                        m_barras[contador].Vref[0] = std::complex<double>(VA*cos(TA*M_PI/180), VA*sin(TA*M_PI/180));
                        m_barras[contador].Vref[1] = std::complex<double>(VB*cos(TB*M_PI/180), VB*sin(TB*M_PI/180));
                        m_barras[contador].Vref[2] = std::complex<double>(VC*cos(TC*M_PI/180), VC*sin(TC*M_PI/180));
                        // __real__ m_barras[contador].Vref[0] = VA*cos(TA*M_PI/180);
                        // __imag__ m_barras[contador].Vref[0] = VA*sin(TA*M_PI/180);
                        // __real__ m_barras[contador].Vref[1] = VB*cos(TB*M_PI/180);
                        // __imag__ m_barras[contador].Vref[1] = VB*sin(TB*M_PI/180);
                        // __real__ m_barras[contador].Vref[2] = VC*cos(TC*M_PI/180);
                        // __imag__ m_barras[contador].Vref[2] = VC*sin(TC*M_PI/180);
                        
                        m_barras[contador].Vinicial[0] = std::complex<double>(VA*(cos(TA*M_PI/180), sin(TA*M_PI/180)));
                        m_barras[contador].Vinicial[1] = std::complex<double>(VB*(cos(TB*M_PI/180), sin(TB*M_PI/180)));
                        m_barras[contador].Vinicial[2] = std::complex<double>(VC*(cos(TC*M_PI/180), sin(TC*M_PI/180)));
                        
                    }
                    contador++;
                    
                }
                else{ // Inserindo nova carga em DBAR já existente
                    //Leitura das Cargas
                    PA = atof(getField(dados,5));
                    PB = atof(getField(dados,6));
                    PC = atof(getField(dados,7));
                    QA = atof(getField(dados,8));
                    QB = atof(getField(dados,9));
                    QC = atof(getField(dados,10));
                            
                    if ((PA != 0) || (PB != 0) || (PC != 0) || (QA != 0) || (QB != 0) || (QC != 0)){
                        //Inseri um valor de load
                        m_barras[aux].nloads++;

                        k = m_barras[aux].nloads - 1;
                        m_barras[aux].loads[k].ID = m_barras[aux].ID;
                        m_barras[aux].loads[k].lig = static_cast<LIGACAO>(atoi(getField(dados,2)));
                        m_barras[aux].loads[k].fases = static_cast<FASES>(atoi(getField(dados,3)));
                        m_barras[aux].loads[k].Vbase = atof(getField(dados,4));
                        m_barras[aux].loads[k].ZIP = atof(getField(dados,11));
                        
                        m_barras[aux].loads[k].Pnom[0] = PA;
                        m_barras[aux].loads[k].Pnom[1] = PB;
                        m_barras[aux].loads[k].Pnom[2] = PC;
                        m_barras[aux].loads[k].Qnom[0] = QA;
                        m_barras[aux].loads[k].Qnom[1] = QB;
                        m_barras[aux].loads[k].Qnom[2] = QC;
                    }
                }
            }
            numeroBarras = contador;
            return (m_barras);
            
        }
        int leituraDSHNT(DBAR *barras, int numeroBarras){
            std::string dbarFolder = dataDir + dataFolder + "/DSHNT.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL) {
                return (1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int i, aux, k; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            dados = new char[100];
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                //Verifica se a barra já foi criada
                aux = -1;
                for (i=0;i<numeroBarras;i++){
                    if ((barras)[i].ID == atoi(getField(dados,1))){
                        aux = i;
                    }
                }
                if(aux != -1){
                    (barras)[aux].nshunts++;
                    barras[aux].shunts = new DSHNT[barras[aux].nshunts];
                    // if (((barras)[aux].shunts = (DSHNT *)realloc((barras)[aux].shunts, (barras)[aux].nshunts  * sizeof(DSHNT)))==NULL)
                    // {
                    //     printf("Erro -- Nao foi possivel alocar espaco de memoria para shunts !!!!");
                    //     exit(1); 
                    // }
                    k = (barras)[aux].nshunts - 1;

                    (barras)[aux].shunts[k].ID = (barras)[aux].ID;
                    (barras)[aux].shunts[k].lig = static_cast<LIGACAO>(atoi(getField(dados,2)));
                    (barras)[aux].shunts[k].fases = static_cast<FASES>(atoi(getField(dados,3)));
                    (barras)[aux].shunts[k].Vbase = atof(getField(dados,4));

                    (barras)[aux].shunts[k].Qnom[0] = atof(getField(dados,5));
                    (barras)[aux].shunts[k].Qnom[1] = atof(getField(dados,6));
                    (barras)[aux].shunts[k].Qnom[2] = atof(getField(dados,7));
                    (barras)[aux].shunts[k].controle = atoi(getField(dados,8));

                    if (getField(dados,10) != NULL){
                        (barras)[aux].shunts[k].num = atoi(getField(dados,9));
                        (barras)[aux].shunts[k].DV = atof(getField(dados,10));
                        (barras)[aux].shunts[k].Vset[0] = atof(getField(dados,11));
                        (barras)[aux].shunts[k].Vset[1] = atof(getField(dados,12));
                        (barras)[aux].shunts[k].Vset[2] = atof(getField(dados,13));
                    }
                }
            }
            fclose(arquivo);
            return (0);
        }
        int leituraDGD(DBAR *barras, int numeroBarras){
            std::string dbarFolder = dataDir + dataFolder + "/DGD.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL) {
                return (1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int i, aux, k; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            dados = new char[100];
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                //Verifica se a barra já foi criada
                aux = -1;
                for (i=0;i<numeroBarras;i++){
                    if ((barras)[i].ID == atoi(getField(dados,1))){
                        aux = i;
                    }
                }
                if(aux != -1){
                    (barras)[aux].ngds++;
                    barras[aux].gds = new DGD[barras[aux].ngds];
                    // if (((barras)[aux].gds = (DGD *)realloc((barras)[aux].gds, (barras)[aux].ngds * sizeof(DGD)))==NULL)
                    // {
                    //     printf("Erro -- Nao foi possivel alocar espaco de memoria para gds !!!!");
                    //     exit(1); 
                    // }
                    k = (barras)[aux].ngds - 1;

                    (barras)[aux].gds[k].ID = (barras)[aux].ID;
                    (barras)[aux].gds[k].lig = static_cast<LIGACAO>(atoi(getField(dados,2)));
                    (barras)[aux].gds[k].fases = static_cast<FASES>(atoi(getField(dados,3)));
                    (barras)[aux].gds[k].Vbase = atof(getField(dados,4));
                    (barras)[aux].gds[k].Snominal = atof(getField(dados,5));

                    (barras)[aux].gds[k].Pnom[0] = atof(getField(dados,6));
                    (barras)[aux].gds[k].Pnom[1] = atof(getField(dados,7));
                    (barras)[aux].gds[k].Pnom[2] = atof(getField(dados,8));
                    (barras)[aux].gds[k].Qnom[0] = atof(getField(dados,9));
                    (barras)[aux].gds[k].Qnom[1] = atof(getField(dados,10));
                    (barras)[aux].gds[k].Qnom[2] = atof(getField(dados,11));
                    (barras)[aux].gds[k].controle = atoi(getField(dados,12));

                    if (getField(dados,14) != NULL){
                        (barras)[aux].gds[k].Qmin = atof(getField(dados,13));
                        (barras)[aux].gds[k].Qmax = atof(getField(dados,14));
                        (barras)[aux].gds[k].Vset[0] = atof(getField(dados,15));
                        (barras)[aux].gds[k].Vset[1] = atof(getField(dados,16));
                        (barras)[aux].gds[k].Vset[2] = atof(getField(dados,17));
                        (barras)[aux].gds[k].controlePV = atoi(getField(dados,18));
                    }
                }
            }
            fclose(arquivo);
            return (0);
        }
        DRAM* leituraDLIN(int &numeroRamos, DBAR *barras, int numeroBarras)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DLIN.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "Erro no arquivo DLIN.csv\n";
                exit(1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            long int contador =0, i; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            dados = new char[100];
            
            //Aloca na memória espaço para as linhas
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            rewind(arquivo);
            DRAM *m_ramos = new DRAM[numLinhas];
            // if (((m_ramos) = (DRAM *)malloc( (numLinhas) * sizeof(DRAM)))==NULL)
            // {
            //     printf("Erro -- Nao foi possivel alocar espaco de memoria para as linhas !!!!");
            //     exit(1); 
            // }
            
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                (m_ramos)[contador].DE = atoi(getField(dados,1));
                (m_ramos)[contador].PARA = atoi(getField(dados,2));
                (m_ramos)[contador].tipo = static_cast<RAMO>(0);
                (m_ramos)[contador].estado = static_cast<ESTADO>(1);
                (m_ramos)[contador].fases = static_cast<FASES>(atoi(getField(dados,3)));
                        
                for(i=0;i<numeroBarras;i++){
                    if((m_ramos)[contador].DE == (barras)[i].ID ){
                        (m_ramos)[contador].k = (barras)[i].i;
                    }
                    if((m_ramos)[contador].PARA == (barras)[i].ID ){
                        (m_ramos)[contador].m = (barras)[i].i;
                    }
                }
                
                //Preenche o (m_ramos)[contador].linha
                (m_ramos)[contador].linha.fases = static_cast<FASES>(atoi(getField(dados,3)));
                (m_ramos)[contador].linha.comprimento = atof(getField(dados,4));
                double comp = atof(getField(dados,4));
                (m_ramos)[contador].linha.Zaa.real(comp*atof(getField(dados,5)));
                (m_ramos)[contador].linha.Zaa.imag(comp*atof(getField(dados,6)));
                (m_ramos)[contador].linha.Zab.real(comp*atof(getField(dados,7)));
                (m_ramos)[contador].linha.Zab.imag(comp*atof(getField(dados,8)));
                (m_ramos)[contador].linha.Zac.real(comp*atof(getField(dados,9)));
                (m_ramos)[contador].linha.Zac.imag(comp*atof(getField(dados,10)));
                (m_ramos)[contador].linha.Zbb.real(comp*atof(getField(dados,11)));
                (m_ramos)[contador].linha.Zbb.imag(comp*atof(getField(dados,12)));
                (m_ramos)[contador].linha.Zbc.real(comp*atof(getField(dados,13)));
                (m_ramos)[contador].linha.Zbc.imag(comp*atof(getField(dados,14)));
                (m_ramos)[contador].linha.Zcc.real(comp*atof(getField(dados,15)));
                (m_ramos)[contador].linha.Zcc.imag(comp*atof(getField(dados,16)));
                //__real__ (m_ramos)[contador].linha.Zaa = comp*atof(getField(dados,5));
                // __imag__ (m_ramos)[contador].linha.Zaa = comp*atof(getField(dados,6));
                // __real__ (m_ramos)[contador].linha.Zab = comp*atof(getField(dados,7));
                // __imag__ (m_ramos)[contador].linha.Zab = comp*atof(getField(dados,8));
                // __real__ (m_ramos)[contador].linha.Zac = comp*atof(getField(dados,9));
                // __imag__ (m_ramos)[contador].linha.Zac = comp*atof(getField(dados,10));
                // __real__ (m_ramos)[contador].linha.Zbb = comp*atof(getField(dados,11));
                // __imag__ (m_ramos)[contador].linha.Zbb = comp*atof(getField(dados,12));
                // __real__ (m_ramos)[contador].linha.Zbc = comp*atof(getField(dados,13));
                // __imag__ (m_ramos)[contador].linha.Zbc = comp*atof(getField(dados,14));
                // __real__ (m_ramos)[contador].linha.Zcc = comp*atof(getField(dados,15));
                // __imag__ (m_ramos)[contador].linha.Zcc = comp*atof(getField(dados,16));
                (m_ramos)[contador].linha.Baa = comp*atof(getField(dados,17))/1000000;
                (m_ramos)[contador].linha.Bab = comp*atof(getField(dados,18))/1000000;
                (m_ramos)[contador].linha.Bac = comp*atof(getField(dados,19))/1000000;
                (m_ramos)[contador].linha.Bbb = comp*atof(getField(dados,20))/1000000;
                (m_ramos)[contador].linha.Bbc = comp*atof(getField(dados,21))/1000000;
                (m_ramos)[contador].linha.Bcc = comp*atof(getField(dados,22))/1000000;
                
                contador++;     
                
            }
            numeroRamos = contador;
            return (m_ramos);
        }
        DRAM* leituraDTRF(int &numeroTrafos, DBAR *barras, int numeroBarras)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DTRF.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "DTRF not found.csv\n";
                numeroTrafos = 0;
                return(NULL);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int contador =0, i;/* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            dados = new char[100];
            
            //Aloca na memória espaço para os trafos
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            rewind(arquivo);
            // DRAM *m_ramos;
            DRAM *m_ramos = new DRAM[numLinhas];
            // if (((ramos) = (DRAM *)realloc((ramos), (numeroRamos + numLinhas) * sizeof(DRAM)))==NULL)
            // {
            //     printf("Erro -- Nao foi possivel alocar espaco de memoria para as linhas !!!!");
            //     exit(1); 
            // }
                
            contador = 0;
            // Le o arquivo até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                (m_ramos)[contador].DE = atoi(getField(dados,1));
                (m_ramos)[contador].PARA = atoi(getField(dados,2));
                (m_ramos)[contador].tipo = static_cast<RAMO>(1);
                (m_ramos)[contador].estado = static_cast<ESTADO>(1);
                (m_ramos)[contador].fases = static_cast<FASES>(atoi(getField(dados,3)));
                
                for(i=0;i<numeroBarras;i++){
                    if((m_ramos)[contador].DE == (barras)[i].ID ){
                        (m_ramos)[contador].k = (barras)[i].i;
                    }
                    if((m_ramos)[contador].PARA == (barras)[i].ID ){
                        (m_ramos)[contador].m = (barras)[i].i;
                    }
                }
                
                //Preencher o (m_ramos)[contador].trafo
                (m_ramos)[contador].trafo.fases = static_cast<FASES>(atoi(getField(dados,3)));
                (m_ramos)[contador].trafo.Vpri = atof(getField(dados,4));
                (m_ramos)[contador].trafo.Vsec = atof(getField(dados,5));
                (m_ramos)[contador].trafo.Snominal = atof(getField(dados,6));
                (m_ramos)[contador].trafo.R = atof(getField(dados,7)) * pow((m_ramos)[contador].trafo.Vsec,2)/((m_ramos)[contador].trafo.Snominal*1000)/3; //R*(Vpri^2/(S*1000))
                (m_ramos)[contador].trafo.X = atof(getField(dados,8)) * pow((m_ramos)[contador].trafo.Vsec,2)/((m_ramos)[contador].trafo.Snominal*1000)/3;
                (m_ramos)[contador].trafo.lig_pri = static_cast<LIGACAO>(atoi(getField(dados,9)));
                (m_ramos)[contador].trafo.lig_sec = static_cast<LIGACAO>(atoi(getField(dados,10)));
                (m_ramos)[contador].trafo.defasamento = atoi(getField(dados,11));
                (m_ramos)[contador].trafo.tap_pri = atof(getField(dados,12));
                (m_ramos)[contador].trafo.tap_sec = atof(getField(dados,13));
                    
                contador++;
            }
            numeroTrafos = numLinhas;
            return (m_ramos);
        }
        DRAM* leituraDREG(int &numeroRegs, DBAR *barras, int numeroBarras)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DREG.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "DREG not found\n";
                numeroRegs = 0;
                return(NULL);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int contador =0, i; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            dados = new char[100];
            
            //Aloca na memória espaço para os reguladores
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            rewind(arquivo);
            DRAM *m_ramos = new DRAM[numLinhas];
            // if (((ramos) = (DRAM *)realloc((ramos), (numeroRamos + numLinhas) * sizeof(DRAM)))==NULL)
            // {
            //     printf("Erro -- Nao foi possivel alocar espaco de memoria para os reguladores de tensão !!!!");
            //     exit(1); 
            // }
            contador = 0;
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                (m_ramos)[contador].DE = atoi(getField(dados,1));
                (m_ramos)[contador].PARA = atoi(getField(dados,2));
                (m_ramos)[contador].tipo = static_cast<RAMO>(2);
                (m_ramos)[contador].estado = static_cast<ESTADO>(1);
                (m_ramos)[contador].fases = static_cast<FASES>(atoi(getField(dados,3)));
                
                for(i=0;i<numeroBarras;i++){
                    if((m_ramos)[contador].DE == (barras)[i].ID ){
                        (m_ramos)[contador].k = (barras)[i].i;
                    }
                    if((m_ramos)[contador].PARA == (barras)[i].ID ){
                        (m_ramos)[contador].m = (barras)[i].i;
                    }             
                }
                
                //Preencher o (*ramos)[contador].regulador
                (m_ramos)[contador].regulador.fases = static_cast<FASES>(atoi(getField(dados,3)));
                (m_ramos)[contador].regulador.Vnom = atof(getField(dados,4));
                (m_ramos)[contador].regulador.regulacao = atof(getField(dados,5));
                (m_ramos)[contador].regulador.ntaps = atoi(getField(dados,6));
                (m_ramos)[contador].regulador.Snominal = atof(getField(dados,7));
                (m_ramos)[contador].regulador.R = atof(getField(dados,8)) * pow((m_ramos)[contador].regulador.Vnom,2)/((m_ramos)[contador].regulador.Snominal*1000)/3;
                (m_ramos)[contador].regulador.X = atof(getField(dados,9)) * pow((m_ramos)[contador].regulador.Vnom,2)/((m_ramos)[contador].regulador.Snominal*1000)/3;
                (m_ramos)[contador].regulador.lig = static_cast<LIGACAO>(atoi(getField(dados,10)));
                (m_ramos)[contador].regulador.TP = atof(getField(dados,11));
                (m_ramos)[contador].regulador.TC = atof(getField(dados,12));
                (m_ramos)[contador].regulador.deltaV = atof(getField(dados,13));
                (m_ramos)[contador].regulador.R1 = atof(getField(dados,14));
                (m_ramos)[contador].regulador.X1 = atof(getField(dados,15));
                (m_ramos)[contador].regulador.R2 = atof(getField(dados,16));
                (m_ramos)[contador].regulador.X2 = atof(getField(dados,17));
                (m_ramos)[contador].regulador.R3 = atof(getField(dados,18));
                (m_ramos)[contador].regulador.X3 = atof(getField(dados,19));
                (m_ramos)[contador].regulador.V1 = atof(getField(dados,20));
                (m_ramos)[contador].regulador.V2 = atof(getField(dados,21));
                (m_ramos)[contador].regulador.V3 = atof(getField(dados,22));
                (m_ramos)[contador].regulador.controle = atoi(getField(dados,23));
                (m_ramos)[contador].regulador.tap[0] = atof(getField(dados,24));
                (m_ramos)[contador].regulador.tap[1] = atof(getField(dados,25));
                (m_ramos)[contador].regulador.tap[2] = atof(getField(dados,26));
                
                //Se tiver parametros de controle reverso
                /*
                (*ramos)[contador].regulador.deltaVr = atof(getField(dados,27));
                (*ramos)[contador].regulador.R1r = atof(getField(dados,28));
                (*ramos)[contador].regulador.X1r = atof(getField(dados,29));
                (*ramos)[contador].regulador.R2r = atof(getField(dados,30));
                (*ramos)[contador].regulador.X2r = atof(getField(dados,31));
                (*ramos)[contador].regulador.R3r = atof(getField(dados,32));
                (*ramos)[contador].regulador.X3r = atof(getField(dados,33));
                (*ramos)[contador].regulador.V1r = atof(getField(dados,34));
                (*ramos)[contador].regulador.V2r = atof(getField(dados,35));
                (*ramos)[contador].regulador.V3r = atof(getField(dados,36));
                */
                contador++;
            }
            numeroRegs = numLinhas;
            return(m_ramos);
        }
        DRAM* leituraDSWTC(int &numeroChaves, DBAR *barras, int numeroBarras)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DSWTC.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "DSWTC not found\n";
                numeroChaves = 0;
                return(NULL);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int contador =0, i; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            dados = new char[100];
            
            //Aloca na memória espaço para os trafos
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            rewind(arquivo);
            DRAM *m_ramos = new DRAM[numLinhas];
            // if (((*ramos) = (DRAM *)realloc((*ramos), (numeroRamos[0] + numLinhas) * sizeof(DRAM)))==NULL)
            // {
            //     printf("Erro -- Nao foi possivel alocar espaco de memoria para as chaves !!!!");
            //     exit(1); 
            // }
                
            contador = 0;
            // Le o arquivo até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                (m_ramos)[contador].DE = atoi(getField(dados,1));
                (m_ramos)[contador].PARA = atoi(getField(dados,2));
                (m_ramos)[contador].tipo = static_cast<RAMO>(3);
                
                for(i=0;i<numeroBarras;i++){
                    if((m_ramos)[contador].DE == (barras)[i].ID ){
                        (m_ramos)[contador].k = (barras)[i].i;
                    }
                    if((m_ramos)[contador].PARA == (barras)[i].ID ){
                        (m_ramos)[contador].m = (barras)[i].i;
                    }
                }
                (m_ramos)[contador].fases = static_cast<FASES>(atoi(getField(dados,3)));
                (m_ramos)[contador].estado = static_cast<ESTADO>(atoi(getField(dados,4)));
                
                contador++;
            }
            numeroChaves =  numLinhas;
            return(m_ramos);
        }
        DRAM* joinDRAM(DRAM *ramos, DRAM *trafos, DRAM *reguladores, DRAM *chaves){
            DRAM *m_ramos = new DRAM[numeroRamos+numeroTrafos+numeroRegs+numeroChaves];
            for (i = 0; i<numeroRamos; i++)
            {
                m_ramos[i].DE = ramos[i].DE;
                m_ramos[i].PARA = ramos[i].PARA;
                m_ramos[i].k = ramos[i].k;
                m_ramos[i].m = ramos[i].m;
                m_ramos[i].fases = ramos[i].fases;
                m_ramos[i].tipo = ramos[i].tipo;
                m_ramos[i].estado = ramos[i].estado;

                m_ramos[i].linha = ramos[i].linha;
                m_ramos[i].trafo = ramos[i].trafo;
                m_ramos[i].regulador = ramos[i].regulador;

                m_ramos[i].Ypp = ramos[i].Ypp;
                m_ramos[i].Yps = ramos[i].Yps;
                m_ramos[i].Ysp = ramos[i].Ysp;
                m_ramos[i].Yss = ramos[i].Yss;
                
                for (j = 0; j<3; j++){
                    m_ramos[i].tap_pri[j] = ramos[i].tap_pri[j];
                    m_ramos[i].tap_sec[j] = ramos[i].tap_sec[j];
                }
                
            }
            if (trafos != NULL){
                for (i = 0; i<numeroTrafos; i++)
                {
                    m_ramos[i+numeroRamos].DE = trafos[i].DE;
                    m_ramos[i+numeroRamos].PARA = trafos[i].PARA;
                    m_ramos[i+numeroRamos].k = trafos[i].k;
                    m_ramos[i+numeroRamos].m = trafos[i].m;
                    m_ramos[i+numeroRamos].fases = trafos[i].fases;
                    m_ramos[i+numeroRamos].tipo = trafos[i].tipo;
                    m_ramos[i+numeroRamos].estado = trafos[i].estado;

                    m_ramos[i+numeroRamos].linha = trafos[i].linha;
                    m_ramos[i+numeroRamos].trafo = trafos[i].trafo;
                    m_ramos[i+numeroRamos].regulador = trafos[i].regulador;

                    m_ramos[i+numeroRamos].Ypp = trafos[i].Ypp;
                    m_ramos[i+numeroRamos].Yps = trafos[i].Yps;
                    m_ramos[i+numeroRamos].Ysp = trafos[i].Ysp;
                    m_ramos[i+numeroRamos].Yss = trafos[i].Yss;

                    for (j = 0; j<3; j++){
                        m_ramos[i+numeroRamos].tap_pri[j] = trafos[i].tap_pri[j];
                        m_ramos[i+numeroRamos].tap_sec[j] = trafos[i].tap_sec[j];
                    }
                }
            }
            if (reguladores != NULL){
                    for (i = 0; i<numeroTrafos; i++)
                    {
                        m_ramos[i+numeroRamos+numeroTrafos].DE = reguladores[i].DE;
                        m_ramos[i+numeroRamos+numeroTrafos].PARA = reguladores[i].PARA;
                        m_ramos[i+numeroRamos+numeroTrafos].k = reguladores[i].k;
                        m_ramos[i+numeroRamos+numeroTrafos].m = reguladores[i].m;
                        m_ramos[i+numeroRamos+numeroTrafos].fases = reguladores[i].fases;
                        m_ramos[i+numeroRamos+numeroTrafos].tipo = reguladores[i].tipo;
                        m_ramos[i+numeroRamos+numeroTrafos].estado = reguladores[i].estado;

                        m_ramos[i+numeroRamos+numeroTrafos].linha = reguladores[i].linha;
                        m_ramos[i+numeroRamos+numeroTrafos].trafo = reguladores[i].trafo;
                        m_ramos[i+numeroRamos+numeroTrafos].regulador = reguladores[i].regulador;

                        m_ramos[i+numeroRamos+numeroTrafos].Ypp = reguladores[i].Ypp;
                        m_ramos[i+numeroRamos+numeroTrafos].Yps = reguladores[i].Yps;
                        m_ramos[i+numeroRamos+numeroTrafos].Ysp = reguladores[i].Ysp;
                        m_ramos[i+numeroRamos+numeroTrafos].Yss = reguladores[i].Yss;

                        for (j = 0; j<3; j++){
                            m_ramos[i+numeroRamos+numeroTrafos].tap_pri[j] = reguladores[i].tap_pri[j];
                            m_ramos[i+numeroRamos+numeroTrafos].tap_sec[j] = reguladores[i].tap_sec[j];
                        }
                    }
            }
            

            if (chaves != NULL){
                for (i = 0; i<numeroChaves; i++)
                    {
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].DE = chaves[i].DE;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].PARA = chaves[i].PARA;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].k = chaves[i].k;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].m = chaves[i].m;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].fases = chaves[i].fases;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].tipo = chaves[i].tipo;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].estado = chaves[i].estado;

                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].linha = chaves[i].linha;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].trafo = chaves[i].trafo;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].regulador = chaves[i].regulador;

                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].Ypp = chaves[i].Ypp;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].Yps = chaves[i].Yps;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].Ysp = chaves[i].Ysp;
                        m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].Yss = chaves[i].Yss;

                        for (j = 0; j<3; j++){
                            m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].tap_pri[j] = chaves[i].tap_pri[j];
                            m_ramos[i+numeroRamos+numeroTrafos+numeroChaves].tap_sec[j] = chaves[i].tap_sec[j];
                        }
                    }
                
            }
            return (m_ramos);
        }
        int leituraVinicial(DBAR *barras, int numeroBarras){
            std::string dbarFolder = dataDir + dataFolder + "/Vinicial.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                return (1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            int i, aux; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            double Va,Vb,Vc,Ta,Tb,Tc;
            dados = new char[100];
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                //Verifica se a barra já foi criada
                aux = -1;
                for (i=0;i<numeroBarras;i++){
                    if ((barras)[i].ID == atoi(getField(dados,1))){
                        aux = i;
                    }
                }
                if(aux != -1){
                    Va = atof(getField(dados,2));
                    Vb = atof(getField(dados,3));
                    Vc = atof(getField(dados,4));
                    Ta = atof(getField(dados,5))*PI/180;
                    Tb = atof(getField(dados,6))*PI/180;
                    Tc = atof(getField(dados,7))*PI/180;
                    
                    if (Va != Va) Va = 1;
                    if (Vb != Vb) Vb = 1;
                    if (Vc != Vc) Vc = 1;
                    if (Ta != Ta) Ta = 0*PI/180;
                    if (Tb != Tb) Tb = -120*PI/180;
                    if (Tc != Tc) Tc = 120*PI/180;
                    
                    (barras)[aux].Vinicial[0] = std::complex<double>(Va*cos(Ta), Va*sin(Ta));
                    (barras)[aux].Vinicial[1] = std::complex<double>(Vb*cos(Tb), Vb*sin(Tb));
                    (barras)[aux].Vinicial[2] = std::complex<double>(Vc*cos(Tc), Vc*sin(Tc));
                    
                }
            }
            return (0);
        }
        GRAFO* geraGrafo(DBAR *barras, int numeroBarras, DRAM *ramos, int numeroRamos){
            long int i,j,k;
            long int barraDe, barraPara,nadj;
            GRAFO* m_grafo = new GRAFO[numeroBarras];
            
            for(i=0;i<numeroBarras;i++){
                (m_grafo)[i].idNo = i;//barras[i].ID;
                (m_grafo)[i].tipo = barras[i].tipo;
                (m_grafo)[i].fases = barras[i].fases;
                (m_grafo)[i].Vbase = barras[i].Vbase;
                (m_grafo)[i].barra = &barras[i];
                (m_grafo)[i].nmed = 0;
                (m_grafo)[i].numeroAdjacentes = 0;
                
                //(*grafo)[i].medidores = (DMED **)malloc(30 * sizeof(DMED*));
                (m_grafo)[i].medidores = (DMED **)malloc(30 * sizeof(DMED*));
                for(j=0;j<3;j++){
                    for(k=0;k<3;k++){
                        (m_grafo)[i].Ysh[j][k] = 0;
                    }
                }
                
            }
            
            for(i=0;i<numeroRamos;i++){
                barraDe = ramos[i].k;
                barraPara = ramos[i].m;
                nadj = (m_grafo)[barraDe].numeroAdjacentes;
                
                (m_grafo)[barraDe].adjacentes[nadj].ramo = &ramos[i];
                (m_grafo)[barraDe].adjacentes[nadj].idNo = barraPara;
                (m_grafo)[barraDe].adjacentes[nadj].tipo = ramos[i].tipo;
                (m_grafo)[barraDe].adjacentes[nadj].estado = ramos[i].estado;
                (m_grafo)[barraDe].adjacentes[nadj].nmed = 0;
                (m_grafo)[barraDe].adjacentes[nadj].idram = i;
                (m_grafo)[barraDe].numeroAdjacentes++;
                if(ramos[i].tipo == 1){
                    (m_grafo)[barraDe].adjacentes[nadj].relacao = ramos[i].trafo.Vsec/ramos[i].trafo.Vpri;
                }
                (m_grafo)[barraDe].adjacentes[nadj].medidores = (DMED **)malloc(30 * sizeof(DMED*));
                
                nadj = (m_grafo)[barraPara].numeroAdjacentes;
                (m_grafo)[barraPara].adjacentes[nadj].ramo = &ramos[i];
                (m_grafo)[barraPara].adjacentes[nadj].idNo = barraDe;
                (m_grafo)[barraPara].adjacentes[nadj].tipo = ramos[i].tipo;
                (m_grafo)[barraPara].adjacentes[nadj].estado = ramos[i].estado;
                (m_grafo)[barraPara].adjacentes[nadj].nmed = 0;
                (m_grafo)[barraPara].adjacentes[nadj].idram = i;
                (m_grafo)[barraPara].numeroAdjacentes++;
                if(ramos[i].tipo == 1){
                    (m_grafo)[barraPara].adjacentes[nadj].relacao = ramos[i].trafo.Vpri/ramos[i].trafo.Vsec;
                }
                (m_grafo)[barraPara].adjacentes[nadj].medidores = (DMED **)malloc(30 * sizeof(DMED*)); 
            } 
            return (m_grafo);
        }

        DMED* leituraMedidas(DRAM *ramos, int numeroRamos, DBAR *barras, int numeroBarras, GRAFO *grafo, double Sbase, int &numeroMeds)
        {
            std::string dbarFolder = dataDir + dataFolder + "/DMED.csv";
            const char *barFolder = dbarFolder.c_str();
            FILE *arquivo = fopen(barFolder, "r");
            if (arquivo == NULL){
                std::cout << "Erro no arquivo DMED.csv\n";
                exit(1);
            }
            char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
            char *dados; /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
            long int contador =0, i,j; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
            int carac,numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
            long int **numeroMedidas, nmed;
            dados = new char[100];
            
            numeroMedidas = (long int**)malloc(14 * sizeof(long int*)); 
            for (i = 0; i < 14; i++){ 
                numeroMedidas[i] = (long int*) malloc(8 * sizeof(long int));
                for (j = 0; j < 8; j++){
                    numeroMedidas[i][j] = 0;
                }
            }
            
            //Aloca na memória espaço para as linhas
            while ((carac = fgetc(arquivo)) != EOF) {
            if (carac == '\n')
                numLinhas++;
            }
            rewind(arquivo);
            DMED* m_medidas = new DMED[numLinhas];
            
            // Le o arquivo de curva de cargas até o fim
            while( (fgets(blocoLeitura, 2000, arquivo))!= NULL ){
                dados = blocoLeitura;
                
                (m_medidas)[contador].ligado = static_cast<ESTADO>(atoi(getField(dados,1)));
                (m_medidas)[contador].tipo = atoi(getField(dados,2));
                (m_medidas)[contador].DE = atoi(getField(dados,3));
                (m_medidas)[contador].PARA = atoi(getField(dados,4));
                (m_medidas)[contador].fases = static_cast<FASES>(atoi(getField(dados,5)));
                (m_medidas)[contador].id = contador;
                (m_medidas)[contador].par = -1;
                
                (m_medidas)[contador].h = 0;
                (m_medidas)[contador].zmed = atof(getField(dados,6));
                (m_medidas)[contador].sigma = 1;//atof(getField(dados,7));
                (m_medidas)[contador].prec = atof(getField(dados,7));
                
                numeroMedidas[(m_medidas)[contador].tipo][(m_medidas)[contador].fases-1]++;
                switch((m_medidas)[contador].tipo){
                    case 0: //Medidas nos ramos
                    case 1:
                    case 6:
                    case 7:
                    case 12: //PMU de corrente retangular real
                    case 13: //PMU de corrente retangular imaginário
                        for(i=0;i<numeroRamos;i++){
                            if(((m_medidas)[contador].DE == ramos[i].DE ) && ((m_medidas)[contador].PARA == ramos[i].PARA )){
                                (m_medidas)[contador].k = ramos[i].k;
                                (m_medidas)[contador].m = ramos[i].m;
                                (m_medidas)[contador].ramo = i;
                            }
                            if(((m_medidas)[contador].DE == ramos[i].PARA ) && ((m_medidas)[contador].PARA == ramos[i].DE )){
                                (m_medidas)[contador].k = ramos[i].m;
                                (m_medidas)[contador].m = ramos[i].k;
                                (m_medidas)[contador].ramo = i;
                            }
                        }
                        break;
                    case 2: //Medidas nas barras
                    case 3:
                    case 4: 
                    case 5:
                    case 8: //PMU de tensão retangular real
                    case 9: //PMU de tensão retangular imaginário
                    case 10: //PMU de injeção de corrente retangular real
                    case 11: //PMU de injeção de corrente retangular imaginário
                        for(i=0;i<numeroBarras;i++){
                            if((m_medidas)[contador].DE == barras[i].ID ){
                                (m_medidas)[contador].k = barras[i].i;
                            }
                            (m_medidas)[contador].m = -1;
                            (m_medidas)[contador].ramo = -1;
                        }
                        break;
                }
                contador++;             
            }
            fclose(arquivo); 
            //Associa as medidas ao grafo e transforma em pu os dados medidos
            nmed = 0;
            for (i = 0; i < 14; i++){ 
                for (j = 0; j < 8; j++){
                    nmed = nmed + numeroMedidas[i][j];
                }
            }
            numeroMeds = nmed;
            return (m_medidas);
        }
        DMED* processaMedidas(DMED *medidas, int numeroMedidas, GRAFO *grafo, int numeroBarras)
        {
            int i, k, m, l;
            int ind, adj;
            double regua;
            double Vbase = 1;
            DMED* m_medidas = new DMED[numeroMedidas];
            for (i = 0; i < numeroMedidas; i++){ 
                k = (m_medidas)[i].k;
                m = (m_medidas)[i].m;
                
                (m_medidas)[i].idAlim = grafo[k].idAlim;
                
                //Associa a medida ao grafo e transforma em pu o valor medido e sigma
                switch ((m_medidas)[i].tipo) {
                    case 0: //Medida de Fluxo de Potência Ativa em kW
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Sbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Sbase/1000);
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 1: //Medida de Fluxo de Potência Reativa em kVAr
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Sbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Sbase/1000);
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 2: //Medida de Injeção de Potência Ativa em kW
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Sbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Sbase/1000);
                                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 6+6*grafo[k].numeroAdjacentes;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        for(l=0;l<grafo[k].numeroAdjacentes;l++){
                            regua = grafo[k].adjacentes[l].idNo;
                            regua += 0.01;
                            for(j=0;j<3;j++){
                                (m_medidas)[i].reguaH[6+6*l+j] = regua;
                                (m_medidas)[i].reguaH[6+6*l+j+3] = -regua;
                                regua += 0.1;
                            }
                        }
                        break;    
                    case 3: //Medida de Injeção de Potência Reativa em kVAr
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Sbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Sbase/1000);
                                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 6+6*grafo[k].numeroAdjacentes;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        for(l=0;l<grafo[k].numeroAdjacentes;l++){
                            regua = grafo[k].adjacentes[l].idNo;
                            regua += 0.01;
                            for(j=0;j<3;j++){
                                (m_medidas)[i].reguaH[6+6*l+j] = regua;
                                (m_medidas)[i].reguaH[6+6*l+j+3] = -regua;
                                regua += 0.1;
                            }
                        }
                        break;
                    case 4: //Medida de Magnitude de Tensão - kV
                        switch ((m_medidas)[i].fases){
                            case 1:
                            case 2:
                            case 3:
                                Vbase = grafo[k].Vbase;
                                break;
                            case 4:
                            case 5:
                            case 6:
                                Vbase = grafo[k].Vbase*(pow(3,0.5));
                                break;
                        }                
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Vbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Vbase/1000);
                                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 6;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            regua += 0.1;
                        }
                        regua = (double)k;
                        regua += 0.01;
                        for(j=3;j<6;j++){
                            (m_medidas)[i].reguaH[j] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 5: //Medida de Ângulo de tensão - graus
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed *PI/180;
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma *PI/180;;
                                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 3;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 6: //Medida de Magnitude de Corrente em A - Fluxo
                        Vbase = grafo[k].Vbase;
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / ((Sbase/1000)/(pow(3,0.5)*(Vbase/1000)));
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / ((Sbase/1000)/(pow(3,0.5)*(Vbase/1000)));
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 7: //Medida de Ângulo de Corrente em graus
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed *PI/180;
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma *PI/180;
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 8: //Medida PMU de tensão retangular Real
                        switch ((m_medidas)[i].fases){
                            case 1:
                            case 2:
                            case 3:
                                Vbase = grafo[k].Vbase;
                                break;
                            case 4:
                            case 5:
                            case 6:
                                Vbase = grafo[k].Vbase*(pow(3,0.5));
                                break;
                        }                
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Vbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Vbase/1000);
                                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 3;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            regua += 0.1;
                        }
                        break;
                    case 9: //Medida PMU de tensão retangular Imaginário
                        switch ((m_medidas)[i].fases){
                            case 1:
                            case 2:
                            case 3:
                                Vbase = grafo[k].Vbase;
                                break;
                            case 4:
                            case 5:
                            case 6:
                                Vbase = grafo[k].Vbase*(pow(3,0.5));
                                break;
                        }                
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / (Vbase/1000);
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / (Vbase/1000);
                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 3;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 10: //Medida PMU de injeção de corrente retangular real
                        Vbase = grafo[k].Vbase;
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / ((Sbase/1000)/((Vbase/1000)));
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / ((Sbase/1000)/((Vbase/1000)));
                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 6+6*grafo[k].numeroAdjacentes;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        for(l=0;l<grafo[k].numeroAdjacentes;l++){
                            regua = grafo[k].adjacentes[l].idNo;
                            regua += 0.01;
                            for(j=0;j<3;j++){
                                (m_medidas)[i].reguaH[6+6*l+j] = regua;
                                (m_medidas)[i].reguaH[6+6*l+j+3] = -regua;
                                regua += 0.1;
                            }
                        }
                        break;
                    case 11: //Medida PMU de injeção de corrente retangular imaginário
                        Vbase = grafo[k].Vbase;
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / ((Sbase/1000)/((Vbase/1000)));
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / ((Sbase/1000)/((Vbase/1000)));
                        
                        ind = grafo[k].nmed;
                        grafo[k].medidores[ind] = &(m_medidas)[i];
                        grafo[k].nmed++;
                        
                        (m_medidas)[i].nvar = 6+6*grafo[k].numeroAdjacentes;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        for(l=0;l<grafo[k].numeroAdjacentes;l++){
                            regua = grafo[k].adjacentes[l].idNo;
                            regua += 0.01;
                            for(j=0;j<3;j++){
                                (m_medidas)[i].reguaH[6+6*l+j] = regua;
                                (m_medidas)[i].reguaH[6+6*l+j+3] = -regua;
                                regua += 0.1;
                            }
                        }
                        
                        break;
                    case 12: //Medida PMU de corrente retangular real
                        Vbase = grafo[k].Vbase;
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / ((Sbase/1000)/((Vbase/1000)));
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / ((Sbase/1000)/((Vbase/1000)));
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                    case 13: //Medida PMU de corrente retangular Imaginário
                        Vbase = grafo[k].Vbase;
                        (m_medidas)[i].zmed = (m_medidas)[i].zmed / ((Sbase/1000)/((Vbase/1000)));
                        (m_medidas)[i].sigma = (m_medidas)[i].sigma / ((Sbase/1000)/((Vbase/1000)));
                        
                        for(j=0;j<grafo[k].numeroAdjacentes;j++){
                            if (grafo[k].adjacentes[j].idNo == m){
                                adj = j;
                            }
                        }
                        ind = grafo[k].adjacentes[adj].nmed;
                        grafo[k].adjacentes[adj].medidores[ind] = &(m_medidas)[i];
                        grafo[k].adjacentes[adj].nmed++;
                        
                        (m_medidas)[i].nvar = 12;
                        (m_medidas)[i].reguaH = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        (m_medidas)[i].H = (double*) malloc ((m_medidas)[i].nvar * sizeof(double));
                        for(j=0;j<(m_medidas)[i].nvar;j++){
                            (m_medidas)[i].H[j] = 0;
                        }
                        
                        regua = (double)k;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j] = regua;
                            (m_medidas)[i].reguaH[j+3] = -regua;
                            regua += 0.1;
                        }
                        regua = (double)m;
                        regua += 0.01;
                        for(j=0;j<3;j++){
                            (m_medidas)[i].reguaH[j+6] = regua;
                            (m_medidas)[i].reguaH[j+9] = -regua;
                            regua += 0.1;
                        }
                        break;
                }        
            }
            //Associa os medidores ao grafo e transforma medidas em pu (equivalente de corrente para o AMB e Baran precisa do par de medida)
            for (i = 0; i < numeroMedidas; i++){ 
                k = (m_medidas)[i].k;
                m = (m_medidas)[i].m;
                
                //Associa a medida ao grafo e transforma em pu o valor medido e sigma
                switch ((m_medidas)[i].tipo) {
                    case 0: //Medida de Fluxo de Potência Ativa em kW
                        for (j = 0; j < numeroMedidas; j++){ 
                            if (((m_medidas)[j].tipo == 1) && ((m_medidas)[j].k == k) && ((m_medidas)[j].m == m)){
                                (m_medidas)[i].par = j;
                                (m_medidas)[j].par = i;
                                j=numeroMedidas;
                            }
                        }
                        break;
                    case 2: //Medida de Injeção de Potência Ativa em kW
                        for (j = 0; j < numeroMedidas; j++){ 
                            if (((m_medidas)[j].tipo == 3) && ((m_medidas)[j].k == k) && ((m_medidas)[j].m == m)){
                                (m_medidas)[i].par = j;
                                (m_medidas)[j].par = i;
                                j=numeroMedidas;
                            }
                        }
                        break;    
                    case 4: //Medida de Magnitude de Tensão - kV
                        for (j = 0; j < numeroMedidas; j++){ 
                            if (((m_medidas)[j].tipo == 5) && ((m_medidas)[j].k == k) && ((m_medidas)[j].m == m)){
                                (m_medidas)[i].par = j;
                                (m_medidas)[j].par = i;
                                j=numeroMedidas;
                            }
                        }
                        break;
                    case 6: //Medida de Magnitude de Corrente em A
                        for (j = 0; j < numeroMedidas; i++){ 
                            if (((m_medidas)[j].tipo == 7) && ((m_medidas)[j].k == k) && ((m_medidas)[j].m == m)){
                                (m_medidas)[i].par = j;
                                (m_medidas)[j].par = i;
                                j=numeroMedidas;
                            }
                        }
                        break;   
                }        
            }

            return (m_medidas);
        }
};
#endif