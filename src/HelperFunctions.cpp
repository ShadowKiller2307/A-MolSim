#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include <array>

/**
 * @brief Helper function to quickly calculate the euclidean norm of an array
 *
 * @param arr The array in question to calculate the euclidean norm
 * @return The euclidean norm of the array as a double
 */
double HelperFunctions::euclideanNorm(const std::array<double, 3> &arr) {
    double sum = 0.0;
    for (const auto &element: arr) {
        sum += element * element;
    }
    return std::sqrt(sum);
}

/**
 * This function is used in the velocity calculation.
 * @param array The array to divide by the scalar. The values are changed inplace.
 * @param scalar The double to divide by.
 * @param isDivision The mode of the operation, where false equals multiplication and true equals division.
 */
void HelperFunctions::scalarOperations(std::array<double, 3> &array, double scalar, bool isDivision) {
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
