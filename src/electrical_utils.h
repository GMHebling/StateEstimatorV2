#ifndef ELECTRICAL_UTILS_H
#define ELECTRICAL_UTILS_H

#include "data_structures.h"
#include "input_data.h"
#include "mat_utils.h"
#include <math.h>
class electrical_utils {
public:
  // Funções das grandezas elétricas
  void Skm(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *S) {
    int i, k;
    std::complex<double> Vp[3], Vs[3];

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    for (k = 0; k < 3; k++) {
      S[k] = std::complex<double>(0, 0);
      for (i = 0; i < 3; i++) {
        S[k] = S[k] + Vp[k] * conj(ramo->Ypp[k][i]) * conj(Vp[i]) +
               Vp[k] * conj(ramo->Yps[k][i]) * conj(Vs[i]);
      }
    }
    switch (ramo->fases) {
    case 1:
      S[1] = std::complex<double>(0, 0);
      ;
      S[2] = std::complex<double>(0, 0);
      ;
      break;
    case 2:
      S[0] = std::complex<double>(0, 0);
      ;
      S[2] = std::complex<double>(0, 0);
      ;
      break;
    case 3:
      S[1] = std::complex<double>(0, 0);
      ;
      S[0] = std::complex<double>(0, 0);
      ;
      break;
    case 4:
      S[2] = std::complex<double>(0, 0);
      ;
      break;
    case 5:
      S[1] = std::complex<double>(0, 0);
      ;
      break;
    case 6:
      S[0] = std::complex<double>(0, 0);
      ;
      break;
    }
  }
  void Smk(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *S) {
    int i, k;
    std::complex<double> Vp[3], Vs[3];

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    for (k = 0; k < 3; k++) {
      S[k] = std::complex<double>(0, 0);
      ;
      for (i = 0; i < 3; i++) {
        S[k] = S[k] + Vs[k] * conj(ramo->Ysp[k][i]) * conj(Vp[i]) +
               Vs[k] * conj(ramo->Yss[k][i]) * conj(Vs[i]);
      }
    }
    switch (ramo->fases) {
    case 1:
      S[1] = 0;
      S[2] = 0;
      break;
    case 2:
      S[0] = 0;
      S[2] = 0;
      break;
    case 3:
      S[1] = 0;
      S[0] = 0;
      break;
    case 4:
      S[2] = 0;
      break;
    case 5:
      S[1] = 0;
      break;
    case 6:
      S[0] = 0;
      break;
    }
  }
  void Sk(GRAFO *grafo, long int k, std::complex<double> *S) {
    int i, j;
    std::complex<double> *Saux;
    Saux = matUtils.c_vetAloca(3);

    for (j = 0; j < 3; j++) {
      S[j] = std::complex<double>(0, 0);
      ;
    }
    for (i = 0; i < grafo[k].numeroAdjacentes; i++) {
      if (k == grafo[k].adjacentes[i].ramo->k)
        Skm(&grafo[k], &grafo[grafo[k].adjacentes[i].idNo],
            grafo[k].adjacentes[i].ramo, Saux);
      else
        Smk(&grafo[grafo[k].adjacentes[i].idNo], &grafo[k],
            grafo[k].adjacentes[i].ramo, Saux);

      for (j = 0; j < 3; j++) {
        S[j] = S[j] - Saux[j];
      }
    }
    // Potência nos shunts
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) {
        S[j] = S[j] +
               grafo[k].V[j] * conj(grafo[k].Ysh[j][i]) * conj(grafo[k].V[i]);
      }
    }

    //free(Saux);
  }
  void Ikm(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *Ikm) {
    int i, k;
    std::complex<double> Vp[3], Vs[3];

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    for (k = 0; k < 3; k++) {
      Ikm[k] = std::complex<double>(0, 0);
      ;
      for (i = 0; i < 3; i++) {
        Ikm[k] = Ikm[k] + ramo->Ypp[k][i] * Vp[i] + ramo->Yps[k][i] * Vs[i];
      }
    }
  }
  void Imk(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *Imk) {
    int i, k;
    std::complex<double> Vp[3], Vs[3];

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    for (k = 0; k < 3; k++) {
      Imk[k] = std::complex<double>(0, 0);
      for (i = 0; i < 3; i++) {
        Imk[k] = Imk[k] + ramo->Ysp[k][i] * Vp[i] + ramo->Yss[k][i] * Vs[i];
      }
    }
  }

  //------------------------------------------------------------------------------
  // Funções das derivadas
  void dSkm(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *dS,
            long int opt, long int i) {
    int j, k;
    std::complex<double> *Vp, *Vs, **J = NULL;

    J = matUtils.c_matAloca(3);
    J[i][i] = std::complex<double>(1, 0);

    Vp = matUtils.c_vetAloca(3);
    Vs = matUtils.c_vetAloca(3);

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    switch (opt) {
    case 0: // Derivada em relação a tensão do noP - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = std::complex<double>(0, 0);
        ;
        for (j = 0; j < 3; j++) {
          dS[k] += J[k][k] * std::complex<double>(1, 0) / std::abs(Vp[k]) *
                       Vp[k] * std::conj(ramo->Ypp[k][j]) * std::conj(Vp[j]) +
                   J[k][k] * std::complex<double>(1, 0) / std::abs(Vp[k]) *
                       Vp[k] * std::conj(ramo->Yps[k][j]) * std::conj(Vs[j]);
        }
        dS[k] += Vp[k] * std::conj(ramo->Ypp[k][i]) *
                 std::complex<double>(1, 0) / std::abs(Vp[i]) *
                 std::conj(Vp[i]);
      }
      break;
    case 1: // Derivada em relação ao ângulo do noP - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = std::complex<double>(0, 0);
        ;
        for (j = 0; j < 3; j++) {
          dS[k] += std::complex<double>(0, 1) * Vp[k] * J[k][k] *
                       std::conj(ramo->Ypp[k][j]) * std::conj(Vp[j]) +
                   std::complex<double>(0, 1) * Vp[k] * J[k][k] *
                       std::conj(ramo->Yps[k][j]) * std::conj(Vs[j]);
        }
        dS[k] += Vp[k] * std::conj(ramo->Ypp[k][i]) *
                 (std::complex<double>(-1, 0)) * std::complex<double>(0, 1) *
                 std::conj(Vp[i]);
      }
      break;
    case 2: // Derivada em relação a tensão do noS - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = Vp[k] * std::conj(ramo->Yps[k][i]) *
                std::complex<double>(1, 0) / std::abs(Vs[i]) * std::conj(Vs[i]);
      }
      break;
    case 3: // Derivada em relação ao ângulo do noS - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = Vp[k] * std::conj(ramo->Yps[k][i]) *
                (std::complex<double>(-1, 0)) * std::complex<double>(0, 1) *
                std::conj(Vs[i]);
      }
      break;
    }
    switch (ramo->fases) {
    case 1:
      dS[1] = 0;
      dS[2] = 0;
      break;
    case 2:
      dS[0] = 0;
      dS[2] = 0;
      break;
    case 3:
      dS[1] = 0;
      dS[0] = 0;
      break;
    case 4:
      dS[2] = 0;
      break;
    case 5:
      dS[1] = 0;
      break;
    case 6:
      dS[0] = 0;
      break;
    }
    // free(Vp);
    // free(Vs);
    // for (i = 0; i < 3; i++)
    //   //free(J[i]);
    // //free(J);
  }

  void dSmk(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *dS,
            long int opt, long int i) {
    int j, k;
    std::complex<double> *Vp, *Vs, **J = NULL;

    J = matUtils.c_matAloca(3);
    J[i][i] = std::complex<double>(1, 0);

    Vp = matUtils.c_vetAloca(3);
    Vs = matUtils.c_vetAloca(3);

    Vp[0] = noP->V[0];
    Vp[1] = noP->V[1];
    Vp[2] = noP->V[2];
    Vs[0] = noS->V[0];
    Vs[1] = noS->V[1];
    Vs[2] = noS->V[2];

    switch (opt) {
    case 0: // Derivada em relação a tensão do noP - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = Vs[k] * std::conj(ramo->Ysp[k][i]) *
                std::complex<double>(1, 0) / std::abs(Vp[i]) * std::conj(Vp[i]);
      }
      break;
    case 1: // Derivada em relação ao ângulo do noP - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = Vs[k] * std::conj(ramo->Ysp[k][i]) *
                (std::complex<double>(-1, 0)) * std::complex<double>(0, 1) *
                std::conj(Vp[i]);
      }
      break;
    case 2: // Derivada em relação a tensão do noS - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = 0;
        for (j = 0; j < 3; j++) {
          dS[k] += J[k][k] * std::complex<double>(1, 0) / std::abs(Vs[k]) *
                       Vs[k] * std::conj(ramo->Ysp[k][j]) * std::conj(Vp[j]) +
                   J[k][k] * std::complex<double>(1, 0) / std::abs(Vs[k]) *
                       Vs[k] * conj(ramo->Yss[k][j]) * conj(Vs[j]);
        }
        dS[k] += Vs[k] * std::conj(ramo->Yss[k][i]) *
                 std::complex<double>(1, 0) / std::abs(Vs[i]) *
                 std::conj(Vs[i]);
      }
      break;
    case 3: // Derivada em relação ao ângulo do noS - DSp / DVpi
      for (k = 0; k < 3; k++) {
        dS[k] = 0;
        for (j = 0; j < 3; j++) {
          dS[k] += std::complex<double>(0, 1) * Vs[k] * J[k][k] *
                       conj(ramo->Ysp[k][j]) * conj(Vp[j]) +
                   std::complex<double>(0, 1) * Vs[k] * J[k][k] *
                       conj(ramo->Yss[k][j]) * conj(Vs[j]);
        }
        dS[k] += Vs[k] * std::conj(ramo->Yss[k][i]) *
                 (std::complex<double>(-1, 0)) * std::complex<double>(0, 1) *
                 std::conj(Vs[i]);
      }
      break;
    }
    switch (ramo->fases) {
    case 1:
      dS[1] = 0;
      dS[2] = 0;
      break;
    case 2:
      dS[0] = 0;
      dS[2] = 0;
      break;
    case 3:
      dS[1] = 0;
      dS[0] = 0;
      break;
    case 4:
      dS[2] = 0;
      break;
    case 5:
      dS[1] = 0;
      break;
    case 6:
      dS[0] = 0;
      break;
    }
    //free(Vp);
    //free(Vs);
    // for (i = 0; i < 3; i++)
      //free(J[i]);
    //free(J);
  }
  void dSk(GRAFO *grafo, long int k, std::complex<double> *dS, long int opt,
           long int barra, long int fase) {
    int i, j, opt1;
    std::complex<double> **J = NULL;

    J = matUtils.c_matAloca(3);
    J[fase][fase] = std::complex<double>(1, 0);

    std::complex<double> *dSaux;
    dSaux = matUtils.c_vetAloca(3);

    for (j = 0; j < 3; j++) {
      dS[j] = std::complex<double>(0, 0);
    }
    if (barra == k) {
      for (i = 0; i < grafo[k].numeroAdjacentes; i++) {
        if (k == grafo[k].adjacentes[i].ramo->k) {
          if (opt == 0)
            opt1 = 0;
          else
            opt1 = 1;
          dSkm(&grafo[k], &grafo[grafo[k].adjacentes[i].idNo],
               grafo[k].adjacentes[i].ramo, dSaux, opt1, fase);
        } else {
          if (opt == 0)
            opt1 = 2;
          else
            opt1 = 3;
          dSmk(&grafo[grafo[k].adjacentes[i].idNo], &grafo[k],
               grafo[k].adjacentes[i].ramo, dSaux, opt1, fase);
        }
        for (j = 0; j < 3; j++) {
          dS[j] = dS[j] - dSaux[j];
        }
      }

      // Derivada da potência nos shunts
      if (opt == 0) { // Derivada em relação a tensão do noP - DSp / DVpi
        for (j = 0; j < 3; j++) {
          for (i = 0; i < 3; i++) {
            dS[j] += J[j][j] * std::complex<double>(1, 0) /
                     std::abs(grafo[k].V[j]) * grafo[k].V[j] *
                     std::conj(grafo[k].Ysh[j][i]) * std::conj(grafo[k].V[i]);
          }
          dS[j] += grafo[k].V[j] * std::conj(grafo[k].Ysh[j][fase]) *
                   std::complex<double>(1, 0) / std::abs(grafo[k].V[fase]) *
                   std::conj(grafo[k].V[fase]);
        }
      } else { // opt==1 >> Derivada em relação ao angulo do noP - DSp / DVpi
        for (j = 0; j < 3; j++) {
          for (i = 0; i < 3; i++) {
            dS[j] += std::complex<double>(0, 1) * grafo[k].V[j] * J[j][j] *
                     std::conj(grafo[k].Ysh[j][i]) * std::conj(grafo[k].V[i]);
          }
          dS[j] += grafo[k].V[j] * std::conj(grafo[k].Ysh[j][fase]) *
                   (std::complex<double>(-1, 0)) * std::complex<double>(0, 1) *
                   std::conj(grafo[k].V[fase]);
        }
      }
    } else {
      for (i = 0; i < grafo[k].numeroAdjacentes; i++) {
        if ((k == grafo[k].adjacentes[i].ramo->k) &&
            (barra == grafo[k].adjacentes[i].ramo->m)) {
          if (opt == 0)
            opt1 = 2;
          else
            opt1 = 3;
          dSkm(&grafo[k], &grafo[grafo[k].adjacentes[i].idNo],
               grafo[k].adjacentes[i].ramo, dSaux, opt1, fase);
          break;
        } else if ((barra == grafo[k].adjacentes[i].ramo->k) &&
                   (k == grafo[k].adjacentes[i].ramo->m)) {
          if (opt == 0)
            opt1 = 0;
          else
            opt1 = 1;
          dSmk(&grafo[grafo[k].adjacentes[i].idNo], &grafo[k],
               grafo[k].adjacentes[i].ramo, dSaux, opt1, fase);
          break;
        }
      }
      for (j = 0; j < 3; j++) {
        dS[j] = dS[j] - dSaux[j];
      }
    }

    //free(dSaux);
//     for (i = 0; i < 3; i++)
//       //free(J[i]);
//     //free(J);
  }

  // //coordenadas retangulares
  // void dSkm_ret(GRAFO *noP, GRAFO *noS, DRAM *ramo, std::complex<double> *dS,
  // long int opt, long int i); void dSmk_ret(GRAFO *noP, GRAFO *noS, DRAM
  // *ramo, std::complex<double> *dS, long int opt, long int i); void
  // dSk_ret(GRAFO *grafo, long int k, std::complex<double> *dS, long int opt,
  // long int barra, long int fase);

private:
  mat_utils matUtils;
};
#endif