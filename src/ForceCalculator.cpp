#include "ForceCalculator.h"
#include <iostream>

/*
void calculateF() {
    for (size_t i = 0; i < particles.size(); ++i)
    {
        Particle &p = particles.at(i);
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
            double scalar = pi.getM() * pj.getM() / std::pow(euclideanNorm(pi.getX() - :q
            getX()), 3);
            std::array<double, 3> force = scalar * (pj.getX() - pi.getX());
            std::array<double, 3> resultingforce = pi.getF() + force;
            pi.setF(resultingforce);

            scalarOperations(force, -1.0, false);
            resultingforce = pj.getF() + force;
            pj.setF(resultingforce);
        }
    }*/
}
