#pragma once
#include "schema.h"
#include <iostream>
#include <array>


class CuboidConstructor{
private:
    std::array<double,3> llfc;
    std::array<unsigned int,3> particlesPerDimension;
    std::array<double,3> particleVelocity;
    double h;
    double mass;
    double type;

public:
    CuboidConstructor(std::array<double,3>& llfc,std::array<unsigned int,3> particlesPerDimension,
                      std::array<double,3> particleVelocity, double h, double mass, double type){
        this->h = h;
        this->llfc = llfc;
        this->mass = mass;
        this->particleVelocity = particleVelocity;
        this->particlesPerDimension = particlesPerDimension;
        this->type = type;
    }
};

class XMLReader {

private:
    std::unique_ptr<simulationConfig> simulation;
    double t_end;
    double delta_t;
    int logLevel;
    bool inputGenerator, inputPicture,inputXML;
    std::string baseName;
    double writeFrequency;
    std::array<double,3> cuboidLLFC;
    std::array<unsigned int,3> cuboidParticlesPerDimension;
    std::array<double,3> particleVelocity;
    double cuboidH;
    double cuboidMass;
    double cuboidType;
    std::vector<CuboidConstructor> cuboidConstructors;


public:
    explicit XMLReader(std::string &path);

    void extractSimulationParameters();

    double getT_end(){
        return this->t_end;
    }
    double getDelta_t(){
        return this->delta_t;
    }
    int getLogLevel(){
        return this->logLevel;
    }

    bool isInputGenerator(){
        return this->inputGenerator;
    }
    bool isInputPicture(){
        return this->inputPicture;
    }
    bool isInputXML(){
        return this->inputXML;
    }

    std::string getBaseName(){
        return this->baseName;
    }
    double getWriteFrequency(){
        return this->writeFrequency;
    }

    /***
 * @brief Constructs an std::array<double,3> from the one in the XML file
 * @param arrayOfThreeDoubles The complex type from the XML file
 * @return the newly constructed std::array<double,3>
 */
    std::array<double,3>  createDoubleArray(arrayOfThreeDoubles arrayOfThreeDoubles);

    /***
 * @brief Constructs an std::array<unsigned int,3> from the one in the XML file
 * @param arrayOfThreeUnsignedInts The complex type from the XML file
 * @return the newly constructed std::array<unsigned int,3>
 */
    std::array<unsigned int,3> createUnsignedIntArray(arrayOfThreeUnsignedInts arrayOfThreeUnsignedInts);

};





