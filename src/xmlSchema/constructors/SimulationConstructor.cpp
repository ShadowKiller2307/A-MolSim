
#include "SimulationConstructor.h"




void SimulationConstructor::setAllSimulationParameters(double t, double delta, int level, int frequency,
                                                       std::array<double,3>& dSize,std::string& cType,
                                                       std::string& name,std::string& bound,double cutOffR,
                                                       double g_Grav, bool use_Thermostat, double initial_Temp, unsigned int
                                                       n_Thermostat, double temp_Target, double max_Diff) {
    this->t_end = t;
    this->delta_t = delta;
    this->logLevel = level;
    this->writeFrequency = frequency;
    this->domainSize = dSize;
    this->containerType = cType;
    this->baseName = name;
    this->boundaries = bound;
    this->cutOffRadius = cutOffR;
    this->gGrav = g_Grav;
    this->useThermostat = use_Thermostat;
    this->initialTemp = initial_Temp;
    this->nThermostat = n_Thermostat;
    this->tempTarget = temp_Target;
    this->maxDiff = max_Diff;
}






double SimulationConstructor::getT_end() const{
    return t_end;
}

double SimulationConstructor::getDelta_t() const{
    return delta_t;
}

int SimulationConstructor::getLogLevel() const{
    return logLevel;
}

int SimulationConstructor::getWriteFrequency() const{
    return writeFrequency;
}

std::string SimulationConstructor::getContainerType() const{
    return containerType;
}


std::string SimulationConstructor::getBaseName(){
    return baseName;
}

std::array<double, 3> SimulationConstructor::getDomainSize() const {
    return domainSize;
}

std::string SimulationConstructor::getBoundaries() const {
    return boundaries;
}

double SimulationConstructor::getCutOffRadius() const {
    return cutOffRadius;
}




