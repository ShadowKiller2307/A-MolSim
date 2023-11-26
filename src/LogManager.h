#pragma once

#include "spdlog/spdlog.h"
#include <spdlog/sinks/rotating_file_sink.h>

class LogManager {
public:

    ///Returns the Singleton instance of the LogManager
    static LogManager &getInstance();

    ///Returns the logger of the LogManager
    std::shared_ptr<spdlog::logger> getLogger();

    ///Sets the log level of the LogManager
    void setLogLevel(spdlog::level::level_enum level);

    ///Returns the log level of the LogManager
    spdlog::level::level_enum getLevel();


    /***
    * Debug logging only if log level is: debug
    * @tparam T The type of the message
    * @tparam Args The type of the arguments for formatting
    * @param message The message that needs to be logged
    * @param args The parameters used for formatting
    */
    template<typename T, typename... Args>
    static void debugLog(T &message, Args &&...args) {
        if (LogManager::getInstance().getLevel() == spdlog::level::debug) {
            LogManager::getInstance().getLogger()->log(spdlog::level::debug, message, std::forward<Args>(args)...);
        }
    }

    /***
    * Info logging only if log level is: info,debug
    * @tparam T The type of the message
    * @tparam Args The type of the arguments for formatting
    * @param message The message that needs to be logged
    * @param args The parameters used for formatting
    */
    template<typename T, typename... Args>
    static void infoLog(T &message, Args &&...args) {
        spdlog::level::level_enum level = LogManager::getInstance().getLevel();
        if (level == spdlog::level::info || level == spdlog::level::debug) {
            LogManager::getInstance().getLogger()->log(spdlog::level::info, message, std::forward<Args>(args)...);
        }
    }

    /***
    * Warning logging only if log level is: warn,info,debug
    * @tparam T The type of the message
    * @tparam Args The type of the arguments for formatting
    * @param message The message that needs to be logged
    * @param args The parameters used for formatting
    */
    template<typename T, typename... Args>
    static void warnLog(T &message, Args &&...args) {
        spdlog::level::level_enum level = LogManager::getInstance().getLevel();
        if (level == spdlog::level::info || level == spdlog::level::debug || level == spdlog::level::warn) {
            LogManager::getInstance().getLogger()->log(spdlog::level::warn, message, std::forward<Args>(args)...);
        }
    }

    /***
     * Error logging only if log level is: err,warn,info,debug
     * @tparam T The type of the message
     * @tparam Args The type of the arguments for formatting
     * @param message The message that needs to be logged
     * @param args The parameters used for formatting
     */
    template<typename T, typename... Args>
    static void errorLog(T &message, Args &&...args) {
        spdlog::level::level_enum level = LogManager::getInstance().getLevel();
        if (level == spdlog::level::info || level == spdlog::level::debug || level == spdlog::level::warn
            || level == spdlog::level::err) {
            LogManager::getInstance().getLogger()->log(spdlog::level::err, message, std::forward<Args>(args)...);
        }
    }


    ///Deletes the copy constructor of the class
    LogManager(const LogManager &) = delete;

    ///Deletes the assignment operator
    LogManager &operator=(const LogManager &) = delete;

private:
    std::shared_ptr<spdlog::logger> logger;
    const int maxSize = 5 * 1024 * 1024;
    const int maxFiles = 6;
    spdlog::level::level_enum logLevel;


    LogManager();

    ~LogManager() = default;


};



