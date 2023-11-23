#include "LennardJonesForce.h"
#include <vector>
#include "HelperFunctions.h"
#include "utils/ArrayUtils.h"
#include <iostream>

LennardJonesForce::LennardJonesForce(double epsilon, double sigma)
{
    this->epsilon = epsilon;
    this->sigma = sigma;
}

// TODO: abstract parts of the force calculation to make the code more beautiful
void LennardJonesForce::calculateForces(std::vector<Particle> &particles)
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
            std::array<double, 3> xDiff = pi.getX() - pj.getX();
            double norm = HelperFunctions::euclideanNorm(xDiff);
            double scalar = -(24 * epsilon) / (std::pow(norm, 2));
            double scalar2 = std::pow(sigma / norm, 6);
            double scalar3 = 2 * std::pow(sigma / norm, 12);
            HelperFunctions::scalarOperations(xDiff, scalar * (scalar2 - scalar3), false);
            // neue Force ist jetzt in xDiff
            std::array<double, 3> resultingForce = pi.getF() + xDiff;
            pi.setF(resultingForce);
            HelperFunctions::scalarOperations(xDiff, -1.0, false);
            resultingForce = pj.getF() + xDiff;
            pj.setF(resultingForce);
        }
    }
}

void LennardJonesForce::calculateForcesWithLambda(ParticleContainer &container)
{
    /*double epsilonCapture = this->epsilon;
    double sigmaCapture = this->sigma;
    auto forceLambda = [epsilonCapture, sigmaCapture](Particle &pi, Particle &pj) {
        std::array<double, 3> xDiff = pi.getX() - pj.getX();
        double norm = HelperFunctions::euclideanNorm(xDiff);
        double scalar = -(24 * epsilonCapture) / (std::pow(norm, 2));
        double scalar2 = std::pow(sigmaCapture / norm, 6);
        double scalar3 = 2 * std::pow(sigmaCapture / norm, 12);
        HelperFunctions::scalarOperations(xDiff, scalar * (scalar2 - scalar3), false);
        // neue Force ist jetzt in xDiff
        std::array<double, 3> resultingForce = pi.getF() + xDiff;
        pi.setF(resultingForce);
        HelperFunctions::scalarOperations(xDiff, -1.0, false);
        resultingForce = pj.getF() + xDiff;
        pj.setF(resultingForce);
    };
    container.iterOverPairs(forceLambda);*/
}