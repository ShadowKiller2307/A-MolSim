
#include "SimulationConstructor.h"




void SimulationConstructor::setAllSimulationParameters(double t, double delta, int level, int frequency,
                                                       bool inputG, bool inputP, bool inputX,
                                                       std::string& name) {
    this->t_end = t;
    this->delta_t = delta;
    this->logLevel = level;
    this->writeFrequency = frequency;
    this->inputGenerator = inputG;
    this->inputPicture = inputP;
    this->inputXML = inputX;
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

bool SimulationConstructor::isInputGenerator() const{
    return inputGenerator;
}
bool SimulationConstructor::isInputPicture() const{
    return inputPicture;
}

bool SimulationConstructor::isInputXML() const{
    return inputXML;
}

std::string SimulationConstructor::getBaseName(){
    return baseName;
}




