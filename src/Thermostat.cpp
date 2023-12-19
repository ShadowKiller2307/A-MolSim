#include <cmath>
#include "Thermostat.h"
#include "utils/ArrayUtils.h"
#include "../src/logOutputManager/LogManager.h"


Thermostat::Thermostat(double initT,double tempTarget, double maxDiff) {

   /* initialTemperature = initT;
    maxTemperatureDifference = maxDiff;
    temperatureTarget = tempTarget;*/
}
/*

void Thermostat::scaleVelocity(std::vector<cell> &cells, double newTemp) {
    if(fabs(currentTemperature-newTemp)>maxTemperatureDifference){
        LogManager::warnLog("Temperature difference exceeds the max allowed: {}\n",
                            maxTemperatureDifference);
        return;
    }
    double div = currentTemperature / newTemp;
    double beta = std::pow(div, 1 / 2.0);
    for (auto &cell: cells) {
        for (auto &p: cell) {
            auto scaledVelocity = beta * p.getV();
            p.setV(scaledVelocity);
        }
    }
}

void Thermostat::gradualScaleVelocity(std::vector<cell> &cells, double tempToAdd){
    double newTemp = currentTemperature + tempToAdd;
    scaleVelocity(cells,newTemp);
}

void Thermostat::adjustTemperature() {

}
double Thermostat::calculateKinEnergy(std::vector<cell> &cells){
    double sum = 0.0;
    for(auto &cell:cells){
        for(auto &p:cell){
            double massP = p.getM();
            double innerProduct = std::inner_product(p.getV().begin(),
                                                     p.getV().end(),p.getV().begin(),0.0);
            double div = (massP*innerProduct)/2;
            sum+=div;
        }
    }
    return sum;
}
double Thermostat::calculateCurrentTemperature(int nrParticles, double kineticEnergy) {
    double div = 2.0/(NR_DIMENSIONS * nrParticles);
    currentTemperature = div*kineticEnergy;
    return div*kineticEnergy;

}

void Thermostat::initializeTemperature(std::vector<cell> &cells) {
    // if a system doesn't have an initial velocity
    // also wir haben einen cuboid an partikeln
    // deren initial velocity 0 ist
    // nun muss deren velocity mittels brownian motion so initialisiert werden
    // dass danach die temperatur des systems Tinit entspricht




}


void Thermostat::regulateTemperature(std::vector<cell> &cells, int nrParticles) {
    if(currentTemperature>temperatureTarget){
        double diff = currentTemperature - temperatureTarget;
        gradualScaleVelocity(cells,-diff);



    }
    else if(currentTemperature<temperatureTarget){
        double diff = temperatureTarget-currentTemperature;
        gradualScaleVelocity(cells,diff);
    }
    double kineticEnergy = calculateKinEnergy(cells);
    calculateCurrentTemperature(nrParticles,kineticEnergy);
}
*/
