#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "data_structures.h"
#include "input_data.h"
#include "mat_utils.h"
#include "wls_utils.h"
#include <chrono>

class estimator {
public:
  void init() {
    std::cout <<"|-------Starting State Estimator------|\n";
    networkData.initialize();
    getEstimatorType();
    verifyVirtualMeasurements();
    std::cout <<"|-------Initialization Completed------|\n";
  }

  void run() {
    switch (esType) {
    case WLS:
      WLSEstimator();
      break;
    case HATCHEL:
      HATCHELEstimator();
      break;
    case BRANCHCURRENT:
      BCEstimator();
      break;
    case AMB:
      AMBEstimator();
      break;
    }
  }

  ESTIMATOR_TYPE esType;
  double tol = 0.000001;
  int MAX_IT = 30;
  int itCount;
  CONV_STATUS conv = NOT_CONVERGENCE;

private:
  input_data networkData;

  void getEstimatorType() {
    this->esType = static_cast<ESTIMATOR_TYPE>(networkData.estimatorType);
  }

  void verifyVirtualMeasurements() {
    if (this->esType == HATCHEL) {
      networkData.processVirtualMeasurements();
    }
  }

  void WLSEstimator() {
    std::cout << "|---------WLS State Estimator---------|\n";
    wls_utils wlsUtils;
    mat_utils mu;
    auto start_est_process = std::chrono::high_resolution_clock::now();
    wlsUtils.initializeStructures(networkData);
    for (itCount = 0; itCount < MAX_IT; itCount++) {
      auto start_iteration = std::chrono::high_resolution_clock::now();
      wlsUtils.updateStructures(networkData);
      wlsUtils.buildNormalEquation(networkData);
      if (itCount == 0)
      {
        wlsUtils.writeSystemToFile();      
      }
      wlsUtils.solveNormalEquation();
      if (wlsUtils.verifyConvergence(wlsUtils.x, tol) == CONVERGENCE) {
        auto end_est_process = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_est_process - start_est_process);
        std::cout << "| Estimation process completed in " << duration.count() << " ms |" << std::endl;
        break;
      }
      wlsUtils.updateSolution(networkData); 
      auto end_iteration = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_iteration - start_iteration);
      std::cout << "| Iteration: " << itCount << " completed in " << duration.count() << " ms -- ||x|| = " << mu.norma_inf(wlsUtils.Dx) << " |" << std::endl;
    }
  }
  void HATCHELEstimator() { std::cout << "HATCHEL Estimator\n"; }
  void BCEstimator() { std::cout << "Branch Current Estimator\n"; }
  void AMBEstimator() { std::cout << "AMB Estimator\n"; }
};
#endif