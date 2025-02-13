// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_STATISTICS_H
#define UNO_STATISTICS_H

#include <string>
#include <map>
#include "tools/Options.hpp"

class Statistics {
public:
   explicit Statistics(const Options& options);

   static int int_width;
   static int double_width;
   static int char_width;

   bool serialize_iterations = false;

   void add_column(std::string name, int width, int order);
   void add_statistic(std::string name, std::string value);
   void add_statistic(std::string name, int value);
   void add_statistic(std::string name, size_t value);
   void add_statistic(std::string name, double value);
   void print_header(bool first_occurrence);
   void print_current_line();
   void print_footer();
   void new_line();

   void add_iteration();
   void serialize();
   std::ostream& write_line(std::ostream& os, std::map<std::string, std::string>& line);

private:
   size_t iteration{0};
   std::map<int, std::string> columns{};
   std::map<std::string, int> widths{};
   std::map<std::string, std::string> current_line{};
   std::map<std::string, std::map<std::string, std::string>> iteration_info;

   size_t print_header_every_iterations{};
   static const std::string& symbol(const std::string& value);
};

#endif // UNO_STATISTICS_H