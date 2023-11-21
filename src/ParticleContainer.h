#pragma once
#include "LogManager.h"
#include "spdlog/spdlog.h"
#include "Particle.h"

class ParticleContainer {

protected:
    template<typename T, typename ...Args>
static void debugLog(T& message,Args ...args){
    std::stringstream stream;
    stream << message;


    LogManager::getInstance().getLogger()->log(isDebug(),message);
}

    template<typename T>
    static void infoLog(T& message){
        LogManager::getInstance().getLogger()->log(isInfo(),message);
    }

static spdlog::level::level_enum isDebug();

static spdlog::level::level_enum isInfo();

public:
    /**
    * @brief iterate over every particle pair in the container and apply the lambda function
    * @param f the lambda function
    */
    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);

};
