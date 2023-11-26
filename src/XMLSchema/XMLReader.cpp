
#include <iostream>
#include "XMLReader.h"


XMLReader::XMLReader(const std::string& path) {
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
    bool inputG = simulation->inputGenerator();
    bool inputP = simulation->inputPicture();
    bool inputX = simulation->inputXML();
    std::string name = simulation->baseName();
    int frequency = simulation->writeFrequency();
    simulationConstructor.setAllSimulationParameters(t,delta,level,frequency,inputG,
                                                     inputP,inputX,name);

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





SimulationConstructor XMLReader::getSimulationConstructor(){
    return this->simulationConstructor;
}

std::vector<CuboidConstructor> XMLReader::getCuboidConstructors() {
    return this->cuboidConstructors;
}


