
#include <iostream>
#include "XMLReader.h"


XMLReader::XMLReader(std::string &path) {
    try {
        this->simulation = Configuration(path, xml_schema::flags::dont_validate);
    }
    catch (xml_schema::exception &e) {
        throw std::invalid_argument(std::string(e.what()));
    }
}

/***
 * @brief Extracts each parameter for the simulation from the XML file.
 * It also extracts each Cuboid from the file and uses CuboidConstructor to create them based on the file.
 */

void XMLReader::extractSimulationParameters() {
    this->t_end = simulation->t_end();
    this->logLevel = simulation->logLevel();
    this->delta_t = simulation->delta_t();
    this->inputGenerator = simulation->inputGenerator();
    this->inputPicture = simulation->inputPicture();
    this->inputXML = simulation->inputXML();
    this->baseName = simulation->baseName();
    this->writeFrequency = simulation->writeFrequency();


    auto cuboids = simulation->Cuboid();


    //iterates over every cuboid and creates a cuboidConstructor for each one.
    //The cuboidConstructor is then added to the vector

    for (auto it = cuboids.begin(); it != cuboids.end(); ++it) {
        auto llfc = createDoubleArray(it->llfc());
        auto particlesPerDimension =
                createUnsignedIntArray(it->particlePerDimension());
        auto velocity = createDoubleArray(it->particleVelocity());
        double h = it->h();
        double mass = it->mass();
        int type = it->generateNumber();

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