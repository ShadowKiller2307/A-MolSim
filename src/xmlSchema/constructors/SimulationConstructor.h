#pragma once
#include <iostream>
#include <array>
class SimulationConstructor {

private:
    //TODO set default values
    double t_end;
    double delta_t;
    int logLevel;
    int writeFrequency;
    std::array<double,3> domainSize;
    std::string containerType;
    std::string baseName;


public:
    SimulationConstructor() = default;

    void setAllSimulationParameters(double t_end, double delta_t, int logLevel, int writeFrequency,
                                    std::array<double,3>& domainSize,std::string& containerType,
                                    std::string& baseName);


    double getT_end() const;
    double getDelta_t() const;
    int getLogLevel() const;
    int getWriteFrequency() const;
    std::string getBaseName();
    std::string getContainerType() const;
    std::array<double,3> getDomainSize() const;

};




