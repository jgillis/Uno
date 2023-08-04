// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Logger.hpp"

#ifdef UNO_SHARED
Level Logger::level = INFO;

void Logger::set_logger(const std::string& logger_level) {
   if (logger_level == "ERROR") {
      Logger::level = ERROR;
   }
   else if (logger_level == "WARNING") {
      Logger::level = WARNING;
   }
   else if (logger_level == "INFO") {
      Logger::level = INFO;
   }
   else if (logger_level == "DEBUG") {
      Logger::level = DEBUG;
   }
   else if (logger_level == "DEBUG2") {
      Logger::level = DEBUG2;
   }
   else {
      throw std::out_of_range("The logger level " + logger_level + " was not found");
   }
}
#endif