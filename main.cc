#include "estimator.h"
#include "input_data.h"
#include "output.h"
#include "electrical_utils.h"
#include "data_structures.h"
#include "mat_utils.h"
#include "optimization_utils.h"
#include "topology_utils.h"
#include "wls_utils.h"
#include "output_utils.h"

#include <complex.h>
#include <iostream>


int main(int argc, char* argv[]){
    estimator EESEP;
    EESEP.init();
    EESEP.run();
}