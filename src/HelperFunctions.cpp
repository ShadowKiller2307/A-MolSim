#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include <array>

double HelperFunctions::euclideanNorm(const std::array<double, 3> &arr)
{
    double sum = 0.0;
    for (const auto &element : arr)
    {
        sum += element * element;
    }
    return std::sqrt(sum);
}

void HelperFunctions::scalarOperations(std::array<double, 3> &array, double scalar, bool isDivision)
{
    if (isDivision)
    {
        for (double &i : array)
        {
            i /= scalar;
        }
    }
    else
    {
        for (double &i : array)
        {
            i *= scalar;
        }
    }
}


