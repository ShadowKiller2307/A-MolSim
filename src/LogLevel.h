#pragma once

/// @brief how much output the program should generate to stdout and (log-)files
enum LogLevel
{
    error,
    warn,
    debug,           // debug, print all, slowest
    standard,        // default, no debug
    off,
    noFiles,
    noCout
};