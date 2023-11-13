/**
 * @brief This class will provide two static methods, on for calculating the euclidean Norm and
 * one for Scalar-Vector multiplication and division
 */

#pragma once

class HelperFunctions
{
public:
    /**
     * @brief Helper function to quickly calculate the euclidean norm of an array
     * @param arr The array in question to calculate the euclidean norm
     * @return The euclidean norm of the array as a double
     */
    static double euclideanNorm(const std::array<double, 3> &arr);

    /**
     * @brief This function is used to perform vector scalar multiplication and division.
     * @param array The array to divide by the scalar. The values are changed inplace.
     * @param scalar The double to divide by.
     * @param isDivision The mode of the operation, where false equals multiplication and true equals division.
     * @return void
     */
    static void scalarOperations(std::array<double, 3> &array, double scalar, bool isDivision);
};
