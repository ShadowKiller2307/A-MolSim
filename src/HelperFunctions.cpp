//
// Created by andreass on 01/11/2023.
//

#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include <array>

class HelperFunctions {
public:
    static double euclideanNorm(const std::array<double, 3> &arr) {
        double sum = 0.0;
        for (const auto &element: arr) {
            sum += element * element;
        }
        return std::sqrt(sum);
    }

    static void scalarOperations(std::array<double, 3> &array, double scalar, bool isDivision) {
        if (isDivision) {
            for (double &i: array) {
                i /= scalar;
            }
        } else {
            for (double &i: array) {
                i *= scalar;
            }
        }
    }
};