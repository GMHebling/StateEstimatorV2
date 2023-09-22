#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "input_data.h"

class estimator {
    public:
        void init() {
            networkData.initialize();
        }

    private:
        input_data networkData;
};
#endif