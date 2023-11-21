#pragma once

#include "LogManager.h"
#include "spdlog/spdlog.h"
#include "Particle.h"

class ParticleContainer
{

public:
    /***
     * @brief Logs only if the level of the LogManager is debug
     * @tparam T The type of message to log
     * @tparam Args The types of the extra arguments
     * @param message The actual message that needs to be logged
     * @param args The arguments used for formatting the message
     */
    template <typename T, typename... Args>
    static void debugLog(T &message, Args &&...args)
    {
        if(LogManager::getInstance().getLevel()==spdlog::level::debug) {
            LogManager::getInstance().getLogger()->log(spdlog::level::debug, message, std::forward<Args>(args)...);
        }
    }




    /**
     * @brief iterate over every particle pair in the container and apply the lambda function
     * @param f the lambda function
     */
    virtual void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);
    /**
     * @brief method to add an Particle to an LC container, thus adding it to the according cell
     */
    // virtual void add(Particle a);
};
