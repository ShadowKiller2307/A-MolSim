#include <vector>
#include <complex>
#include <iostream>
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include "ForceCalculator.h"
#include "HelperFunctions.h"
#include "ForceV1.h"

void ForceV1::calculateForces(std::vector<Particle> &particles)
{
    for (auto &p : particles)
    {
        auto oldForce = p.getF();
        std::array<double, 3> zero = {0.0, 0.0, 0.0};
        p.setF(zero);
        p.setOldF(oldForce);
    }

    for (size_t i = 0; i < particles.size() - 1; ++i)
    {
        Particle &pi = particles.at(i);
        for (size_t j = i + 1; j < particles.size(); ++j)
        {
            Particle &pj = particles.at(j);
            double scalar =
                pi.getM() * pj.getM() / std::pow(HelperFunctions::euclideanNorm(pi.getX() - pj.getX()), 3);
            std::array<double, 3> force = scalar * (pj.getX() - pi.getX());
            std::array<double, 3> resultingForce = pi.getF() + force;
            pi.setF(resultingForce);

            HelperFunctions::scalarOperations(force, -1.0, false);
            resultingForce = pj.getF() + force;
            pj.setF(resultingForce);
        }
    }
}
