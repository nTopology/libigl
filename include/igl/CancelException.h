#pragma once

#include <stdexcept>

class CancelException : public std::runtime_error{
  public:
    explicit CancelException(const std::string& what) : std::runtime_error(what) {}
    explicit CancelException(const char* what) : std::runtime_error(what) {}
};
