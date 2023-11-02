/**
 * @brief This class will provide two static methods, on for calculating the euclidean Norm and
 * one for Scalar-Vector multiplication and division
 */

#pragma once

class HelperFunctions {
public:
    static double euclideanNorm(const std::array<double, 3> &arr);
    static void scalarOperations(std::array<double, 3> &array, double scalar, bool isDivision);
};
