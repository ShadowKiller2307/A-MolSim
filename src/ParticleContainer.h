#pragma once

#include "LogManager.h"
#include "spdlog/spdlog.h"
#include "Particle.h"

class ParticleContainer {

protected:
    template<typename T, typename ...Args>
    static void debugLog(T &message, Args &&...args) {
        LogManager::getInstance().getLogger()->log(isDebug(), message, std::forward<Args>(args)...);
    }


    static spdlog::level::level_enum isDebug();


public:
    /**
    * @brief iterate over every particle pair in the container and apply the lambda function
    * @param f the lambda function
    */
    virtual void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);
    /**
    * @brief method to add an Particle to an LC container, thus adding it to the according cell
    */
    virtual void add(Particle a);

};
