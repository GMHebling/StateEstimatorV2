#ifndef SPARSESYSTEM_H
#define SPARSESYSTEM_H

#include "engine.h"

#include <iostream>

class sparseSystem {
    public: 
        int solver_method;
        int order_method = 0;
        int PrintFlag = 0;
        std::vector<double> output;

        
        int solve(std::vector<std::vector<double> > A, std::vector<double> b)
        {
            solver_engine.order_method = this->order_method;
            solver_engine.engine_run(A, b);
            output = solver_engine.engine_output;
        }


    private:
        engine solver_engine;
        
        


};
#endif