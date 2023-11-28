
#include "SimulationConstructor.h"




void SimulationConstructor::setAllSimulationParameters(double t, double delta, int level, int frequency,
                                                       std::array<double,3>& dSize,std::string& cType,
                                                       std::string& name) {
    this->t_end = t;
    this->delta_t = delta;
    this->logLevel = level;
    this->writeFrequency = frequency;
    this->domainSize = dSize;
    this->containerType = cType;
    this->baseName = name;
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




