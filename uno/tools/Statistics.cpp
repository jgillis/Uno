// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <iomanip>
#include <fstream>
#include "Statistics.hpp"

// TODO move this to the option file
int Statistics::int_width = 7;
int Statistics::double_width = 17;
int Statistics::char_width = 7;

Statistics::Statistics(const Options& options): print_header_every_iterations(options.get_unsigned_int("statistics_print_header_every_iterations")) {
}

void Statistics::add_column(std::string name, int width, int order) {
   this->columns[order] = name;
   this->widths[std::move(name)] = width;
}

void Statistics::add_statistic(std::string name, std::string value) {
   this->current_line[std::move(name)] = std::move(value);
}

void Statistics::add_statistic(std::string name, int value) {
   add_statistic(std::move(name), std::to_string(value));
}

void Statistics::add_statistic(std::string name, size_t value) {
   add_statistic(std::move(name), std::to_string(value));
}

void Statistics::add_statistic(std::string name, double value) {
   std::ostringstream stream;
   stream << std::defaultfloat << std::setprecision(7) << value;
   add_statistic(std::move(name), stream.str());
}

void Statistics::print_header(bool first_occurrence) {
   /* line above */
   std::cout << (first_occurrence ? Statistics::symbol("top-left") : Statistics::symbol("left-mid"));
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << (first_occurrence ? Statistics::symbol("top-mid") : Statistics::symbol("mid-mid"));
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("top");
      }
      k++;
   }
   std::cout << (first_occurrence ? Statistics::symbol("top-right") : Statistics::symbol("right-mid")) << '\n';
   /* headers */
   std::cout << Statistics::symbol("left");
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("middle");
      }
      std::string header = element.second;
      std::cout << " " << header;
      for (int j = 0; j < this->widths[header] - static_cast<int>(header.size()) - 1; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol("right") << '\n';
}

void Statistics::print_current_line() {
   if (this->iteration % this->print_header_every_iterations == 0) {
      this->print_header(this->iteration == 0);
   }
   std::cout << Statistics::symbol("left-mid");
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("mid-mid");
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("bottom");
      }
      k++;
   }
   std::cout << Statistics::symbol("right-mid") << '\n';
   /* headers */
   std::cout << Statistics::symbol("left");
   k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("middle");
      }
      const std::string& header = element.second;
      int size;
      try {
         std::string value = this->current_line.at(header);
         std::cout << " " << value;
         size = 1 + static_cast<int>(value.size());
      }
      catch (const std::out_of_range&) {
         std::cout << " -";
         size = 2;
      }
      int number_spaces = (size <= this->widths[header]) ? this->widths[header] - size : 0;
      for (int j = 0; j < number_spaces; j++) {
         std::cout << " ";
      }
      k++;
   }
   std::cout << Statistics::symbol("right") << '\n';
   this->iteration++;
}

void Statistics::print_footer() {
   std::cout << Statistics::symbol("bottom-left");
   int k = 0;
   for (const auto& element: this->columns) {
      if (0 < k) {
         std::cout << Statistics::symbol("bottom-mid");
      }
      std::string header = element.second;
      for (int j = 0; j < this->widths[header]; j++) {
         std::cout << Statistics::symbol("bottom");
      }
      k++;
   }
   std::cout << Statistics::symbol("bottom-right") << '\n';
}

void Statistics::new_line() {
   this->current_line.clear();
}

// Add the stats of one iteration to the map of all iterations
void Statistics::add_iteration() {
   std::string iterations_number = std::to_string(this->iteration);
   this->iteration_info[std::move(iterations_number)] = std::move(this->current_line);
}

std::ostream& Statistics::write_line(std::ostream& os, std::map<std::string, std::string>& line)
{
   os << "    {\n";

   std::map<int, std::string>::iterator it;
   for (it = this->columns.begin(); it != this->columns.end(); it++) {
      
      if (it != this->columns.begin()){
         os << ",\n"; //Write the comma
      }
      const std::string& header = it->second;
      try {
         std::string value = line.at(header);
         os << "        \"" << header << "\": " << value;
      }
      catch (const std::out_of_range&) {
         os << "        \"" << header << "\": -";
      }
   }

    os << "    }";
    return os;
}


void Statistics::serialize() {

   std::ofstream output_file; 
   output_file.open("uno_statistics.json");
   output_file << "{\n";
   std::map<std::string, std::map<std::string, std::string>>::iterator it;
   for (it = this->iteration_info.begin(); it != this->iteration_info.end(); it++)
   {
      if (it != this->iteration_info.begin()){
         output_file << ",\n"; //Write the comma
      }
      output_file << "    \"" << it->first << "\":\n";
      this->write_line(output_file, it->second);
   }

   output_file << "}";
   output_file.close();
}


const std::string& Statistics::symbol(const std::string& value) {
   static std::map<std::string, std::string> symbols = {
         {"top", "─"},
         {"top-mid", "┬"},
         {"top-left", "┌"},
         {"top-right", "┐"},
         {"bottom", "─"},
         {"bottom-mid", "┴"},
         {"bottom-left", "└"},
         {"bottom-right", "┘"},
         {"left", "│"},
         {"left-mid", "├"},
         {"mid", "─"},
         {"mid-mid", "┼"},
         {"right", "│"},
         {"right-mid", "┤"},
         {"middle", "│"}
   };
   return symbols[value];
}