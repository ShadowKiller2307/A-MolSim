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
    int writeFrequency;
    std::array<double,3> cuboidLLFC;
    std::array<unsigned int,3> cuboidParticlesPerDimension;
    std::array<double,3> particleVelocity;
    double cuboidH;
    double cuboidMass;
    double cuboidType;
    std::vector<CuboidConstructor> cuboidConstructors;


public:
    explicit XMLReader(std::string &path);

    /// @brief Extracts the simulation parameters from the XML file
    void extractSimulationParameters();
    /// @brief Extracts the cuboids from the XML file
    std::vector<CuboidConstructor> extractCuboid();



    /***
 * @brief Constructs an std::array<double,3> from the one in the XML file
 * @param arrayOfThreeDoubles The complex type from the XML file
 * @return the newly constructed std::array<double,3>
 */
    static std::array<double,3>  createDoubleArray(arrayOfThreeDoubles arrayOfThreeDoubles);

    /***
 * @brief Constructs an std::array<unsigned int,3> from the one in the XML file
 * @param arrayOfThreeUnsignedInts The complex type from the XML file
 * @return the newly constructed std::array<unsigned int,3>
 */
    static std::array<unsigned int,3> createUnsignedIntArray(arrayOfThreeUnsignedInts arrayOfThreeUnsignedInts);


    double getT_end();
    double getDelta_t();
    int getLogLevel();
    bool isInputGenerator();
    bool isInputPicture();
    bool isInputXML();
    std::string getBaseName();
    int getWriteFrequency();


};





