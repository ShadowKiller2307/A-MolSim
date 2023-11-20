/**
 * @brief This class will provide two static methods, on for calculating the euclidean Norm and
 * one for Scalar-Vector multiplication and division
 */

#pragma once

#include <string>
#include <iostream>
#include <sstream>


class HelperFunctions {
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

    /***
     * Inspired from https://www.techiedelight.com/convert-int-array-string-cpp/
     * @tparam T The type of the std::array
     * @tparam N The size of the std::array
     * @param array The array that is transformed to a string
     * @return Returns the string representation of the array
     */
    template<typename T, size_t N>
    static std::string arrayToString(std::array<T, N> &array) {
        std::ostringstream buffer;
        std::string comma = ",";
        buffer << "[";
        for (size_t i = 0; i < N; ++i) {
            buffer << array[i];
            if (i != N - 1) {
                buffer << comma;
            }

        }
        buffer << "]";
        return buffer.str();
    }

    /***
     * Inspired from https://www.techiedelight.com/convert-int-array-string-cpp/
     * @tparam T The type of the std::vector
     * @param vector The vector that is transformed to a string
     * @return Returns the string representation of the vector
     */
    template<typename T>
    static std::string vectorToString(std::vector<T> &vector) {
        std::ostringstream buffer;
        std::string comma = ",";
        buffer << "[";
        for (int i = 0; i < vector.size(); ++i) {
            buffer << vector.at(i);
            if (i != vector.size() - 1) {
                buffer << comma;
            }

        }
        buffer << "]";
        return buffer.str();
    }
};
