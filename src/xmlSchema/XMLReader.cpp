
#include <iostream>
#include "XMLReader.h"


XMLReader::XMLReader(std::string &path) {
    try {
        this->simulation = Configuration(path, xml_schema::flags::dont_validate);
        simulationConstructor = SimulationConstructor();
    }
    catch (xml_schema::exception &e) {
        throw std::invalid_argument(std::string(e.what()));
    }
}


void XMLReader::extractSimulationParameters() {
    double t = simulation->t_end();
    int level = simulation->logLevel();
    double delta = simulation->delta_t();

    std::string name = simulation->baseName();
    int frequency = simulation->writeFrequency();
    std::array<double, 3> domainSize = createDoubleArray(simulation->domainSize());
    std::string containerType = simulation->containerType();
    std::string boundaries = simulation->boundaries();
    double cutOffRadius = simulation->cutOffRadius();



    simulationConstructor.setAllSimulationParameters(t, delta, level,
                                                     frequency, domainSize,
                                                     containerType, name,boundaries,cutOffRadius);

}

void XMLReader::extractCuboid() {
    auto cuboids = simulation->Cuboid();
    for (auto &cuboid: cuboids) {
        auto llfc = createDoubleArray(cuboid.llfc());
        auto particlesPerDimension =
                createUnsignedIntArray(cuboid.particlePerDimension());
        auto velocity = createDoubleArray(cuboid.particleVelocity());
        double h = cuboid.h();
        double mass = cuboid.mass();
        int type = cuboid.generateNumber();

        CuboidConstructor cuboidConstructor(llfc, particlesPerDimension, velocity, h, mass, type);
        cuboidConstructors.push_back(cuboidConstructor);

    }

}

void XMLReader::extractSphere() {
    auto spheres = simulation->Sphere();
    for (auto &sphere: spheres) {
        std::array<double, 3> cCoord = createDoubleArray(sphere.centerCoordinates());
        std::array<double, 3> iVel = createDoubleArray(sphere.initialVelocity());
        int radius = sphere.radius();
        double distance = sphere.distance();
        double mass = sphere.mass();

        SphereConstructor sphereConstructor(cCoord, iVel, radius, distance, mass);
        sphereConstructors.push_back(sphereConstructor);
    }
}


std::array<double, 3> XMLReader::createDoubleArray(arrayOfThreeDoubles arrayOfThreeDoubles) {
    std::array<double, 3> createdArray =
            {arrayOfThreeDoubles.at(0), arrayOfThreeDoubles.at(1), arrayOfThreeDoubles.at(2)};
    return createdArray;

}


std::array<unsigned int, 3> XMLReader::createUnsignedIntArray(arrayOfThreeUnsignedInts arrayOfThreeUnsignedInts) {
    std::array<unsigned int, 3> createdArray =
            {arrayOfThreeUnsignedInts.at(0), arrayOfThreeUnsignedInts.at(1), arrayOfThreeUnsignedInts.at(2)};
    return createdArray;
}


SimulationConstructor XMLReader::getSimulationConstructor() {
    return this->simulationConstructor;
}

std::vector<CuboidConstructor> XMLReader::getCuboidConstructors() {
    return this->cuboidConstructors;
}

std::vector<SphereConstructor> XMLReader::getSphereConstructors() {
    return this->sphereConstructors;
}




