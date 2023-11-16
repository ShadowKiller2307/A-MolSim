#pragma once

enum LogLevel
{
    debug,           // debug, print all, slowest
    standard,        // default, no debug
    noCOut,          // disable std::cout and logging but still write Files, minor performance improvement
    noFiles,         // don't write output to files but still std::cout and logging, big performance improvement but no data :(
    onlyCalculations // combination of noCOut and noFiles, no output whatsoever, fastest
};