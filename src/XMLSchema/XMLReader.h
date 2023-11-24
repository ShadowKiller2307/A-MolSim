#pragma once

#include "schema.h"
#include <iostream>
#include <array>
#include "ParticleGenerator.h"
#include "XMLSchema/Constructors/CuboidConstructor.h"
#include "XMLSchema/Constructors/SimulationConstructor.h"


class XMLReader {

private:
    std::unique_ptr<simulationConfig> simulation;

    ///The vector containing the extracted cuboids
    std::vector<CuboidConstructor> cuboidConstructors;

    ///The object containing all the extracted simulation parameters
    SimulationConstructor simulationConstructor;


public:
    /***
     * The use of explicit helps avoiding unexpected conversions
     * @param path The path to the XMl file we want to read
     */
    explicit XMLReader(std::string &path);

    /// @brief Extracts the simulation parameters from the XML file and sets the parameters of the simulationConstructor
    void extractSimulationParameters();

    /// @brief Extracts the cuboids from the XML file
    void extractCuboid();


    /***
 * @brief Constructs an std::array<double,3> from the one in the XML file
 * @param arrayOfThreeDoubles The complex type from the XML file
 * @return the newly constructed std::array<double,3>
 */
    static std::array<double, 3> createDoubleArray(arrayOfThreeDoubles arrayOfThreeDoubles);

    /***
 * @brief Constructs an std::array<unsigned int,3> from the one in the XML file
 * @param arrayOfThreeUnsignedInts The complex type from the XML file
 * @return the newly constructed std::array<unsigned int,3>
 */
    static std::array<unsigned int, 3> createUnsignedIntArray(arrayOfThreeUnsignedInts arrayOfThreeUnsignedInts);



    ///Getter for the simulationConstructor
    SimulationConstructor getSimulationConstructor();
    ///Getter for the vector of cuboidConstructors
    std::vector<CuboidConstructor> getCuboidConstructors();

};





