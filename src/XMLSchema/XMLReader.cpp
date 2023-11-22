
#include <iostream>
#include "XMLReader.h"



XMLReader::XMLReader(std::string &path) {
    try {
        this->simulation = Configuration(path, xml_schema::flags::dont_validate);
        extractSimulationParameters();
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

}

std::vector<CuboidConstructor> XMLReader::extractCuboid(){
    auto cuboids = simulation->Cuboid();
    for (auto & cuboid : cuboids) {
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
    return cuboidConstructors;
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

// Getter functions

double XMLReader::getT_end(){
    return this->t_end;
}
double XMLReader::getDelta_t(){
    return this->delta_t;
}
int XMLReader::getLogLevel(){
    return this->logLevel;
}

bool XMLReader::isInputGenerator(){
    return this->inputGenerator;
}
bool XMLReader::isInputPicture(){
    return this->inputPicture;
}
bool XMLReader::isInputXML(){
    return this->inputXML;
}

std::string XMLReader::getBaseName(){
    return this->baseName;
}
int XMLReader::getWriteFrequency(){
    return this->writeFrequency;
}
