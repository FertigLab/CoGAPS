#ifndef _TRANS_H
#define _TRANS_H

#include <vector>

std::vector<double> logit(std::vector<double> data) {
    std::vector<double> transformation(data.size());

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}

#endif
