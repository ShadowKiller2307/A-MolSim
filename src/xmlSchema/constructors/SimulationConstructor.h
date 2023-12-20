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
    std::string boundaries;
    double cutOffRadius;

    double gGrav = 10;
    bool useThermostat = false;


    double initialTemp;
    unsigned int nThermostat;
    double tempTarget;
    double maxDiff;

public:
    SimulationConstructor() = default;

    void setAllSimulationParameters(double t_end, double delta_t, int logLevel, int writeFrequency,
                                    std::array<double,3>& domainSize,std::string& containerType,
                                    std::string& baseName,std::string& bound,double cutOffR,
                                    double g_Grav, bool use_Thermostat, double initial_Temp, unsigned int
                                    n_Thermostat, double temp_Target, double max_Diff);


    double getT_end() const;
    double getDelta_t() const;
    int getLogLevel() const;
    int getWriteFrequency() const;
    std::string getBaseName();
    std::string getContainerType() const;
    std::array<double,3> getDomainSize() const;
    std::string getBoundaries() const;
    double getCutOffRadius() const;
    double getGrav() const;
    bool isUseThermostat() const;
    double getInitialTemp() const;

};




