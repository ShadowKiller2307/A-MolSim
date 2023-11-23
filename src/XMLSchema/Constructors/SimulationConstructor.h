#pragma once
#include <iostream>
class SimulationConstructor {

private:
    //TODO set default values
    double t_end;
    double delta_t;
    int logLevel;
    int writeFrequency;

    bool inputGenerator, inputPicture, inputXML;
    std::string baseName;


public:
    SimulationConstructor() = default;

    void setAllSimulationParameters(double t_end, double delta_t, int logLevel, int writeFrequency,
                                    bool inputGenerator, bool inputPicture, bool inputXML, std::string& baseName);


    double getT_end() const;
    double getDelta_t() const;
    int getLogLevel() const;
    int getWriteFrequency() const;
    bool isInputGenerator() const;
    bool isInputPicture() const;
    bool isInputXML() const;
    std::string getBaseName();

};




