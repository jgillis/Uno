// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include <cmath>
#include "tools/Logger.hpp"
#include "tools/range.hpp"

enum Norm {
   L1_NORM = 1,
   L2_NORM = 2,
   L2_SQUARED_NORM,
   INF_NORM
};

inline Norm norm_from_string(const std::string& norm_string) {
   if (norm_string == "L1") {
      return L1_NORM;
   }
   else if (norm_string == "L2") {
      return L2_NORM;
   }
   else if (norm_string == "INF") {
      return INF_NORM;
   }
   throw std::invalid_argument("The norm " + norm_string + " is not known");
}

// result <- x + scaling_factor * y
template <typename T>
void add_vectors(const std::vector<T>& x, const std::vector<T>& y, T scaling_factor, std::vector<T>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t i: range(x.size())) {
      result[i] = x[i] + scaling_factor * y[i];
   }
}

template <typename T>
void initialize_vector(std::vector<T>& x, T value) {
   for (T& xi: x) {
      xi = T(value);
   }
}

/*
template <typename T>
void scale(std::vector<T>& x, T scaling_factor) {
   for (T& xi: x) {
      xi *= scaling_factor;
   }
}
*/

/*
template <typename T>
T dot(const std::vector<T>& x, const std::vector<T>& y) {
   assert(x.size() == y.size() && "The vectors do not have the same size.");

   T dot_product = 0.;
   for (size_t i: Range(x.size())) {
      dot_product += x[i]*y[i];
   }
   return dot_product;
}
*/

template <typename T>
void copy_from(std::vector<T>& destination, const std::vector<T>& source, size_t length = std::numeric_limits<size_t>::max()) {
   length = std::min(length, std::min(source.size(), destination.size()));
   const auto source_start_position = std::cbegin(source);
   const auto source_end_position = source_start_position + length;
   const auto destination_position = std::begin(destination);
   std::copy(source_start_position, source_end_position, destination_position);
}

// norms of any array with elements of any type

// compute l1 norm = sum |x|_i
template <typename ARRAY, typename T = typename ARRAY::value_type>
T norm_1(const ARRAY& x) {
   T norm{0};
   for (size_t i = 0; i < x.size(); i++) {
      norm += std::abs(x[i]);
   }
   return norm;
}

// l1 norm of several arrays
template<typename ARRAY, typename... ARRAYS, typename T = typename ARRAY::value_type>
T norm_1(const ARRAY& x, ARRAYS... other_arrays) {
   return norm_1(x) + norm_1(other_arrays...);
}

// compute l2 squared norm = sum x_i^2
template <typename ARRAY, typename T = typename ARRAY::value_type>
T norm_2_squared(const ARRAY& x) {
   T norm_squared{0};
   for (size_t i = 0; i < x.size(); i++) {
      const T xi = x[i];
      norm_squared += xi * xi;
   }
   return norm_squared;
}

// l2 squared norm of several arrays
template<typename ARRAY, typename... ARRAYS, typename T = typename ARRAY::value_type>
T norm_2_squared(const ARRAY& x, ARRAYS... other_arrays) {
   return norm_2_squared(x) + norm_2_squared(other_arrays...);
}

// compute ||x||_2
template <typename ARRAY, typename T = typename ARRAY::value_type>
T norm_2(const ARRAY& x) {
   return std::sqrt(norm_2_squared(x));
}

// l2 norm of several arrays
template<typename ARRAY, typename... ARRAYS, typename T = typename ARRAY::value_type>
T norm_2(const ARRAY& x, ARRAYS... other_arrays) {
   return std::sqrt(norm_2_squared(x) + norm_2_squared(other_arrays...));
}

// compute ||x||_inf
template <typename ARRAY, typename T = typename ARRAY::value_type>
T norm_inf(const ARRAY& x) {
   T norm{0};
   for (size_t i = 0; i < x.size(); i++) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

// inf norm of several arrays
template<typename ARRAY, typename... ARRAYS, typename T = typename ARRAY::value_type>
T norm_inf(const ARRAY& x, ARRAYS... other_arrays) {
   return std::max(norm_inf(x), norm_inf(other_arrays...));
}

template <typename ARRAY, typename T = typename ARRAY::value_type>
T norm(const ARRAY& x, Norm norm) {
   // choose the right norm
   if (norm == L1_NORM) {
      return norm_1(x);
   }
   else if (norm == L2_NORM) {
      return norm_2(x);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared(x);
   }
   else if (norm == INF_NORM) {
      return norm_inf(x);
   }
   throw std::invalid_argument("The norm is not known");
}

// these methods take:
// - a callback as argument whose parameter is the current index. This avoids forming the vector explicitly
// - an array of arbitrary type (can be Range, std::vector, std::array, etc)
template <typename T, typename ARRAY>
T norm_1(const std::function<T(size_t i)>& ith_component, const ARRAY& array) {
   T norm = T(0);
   for (size_t i: array) {
      norm += std::abs(ith_component(i));
   }
   return norm;
}

template <typename T, typename ARRAY>
T norm_inf(const std::vector<T>& x, const ARRAY& array) {
   T norm = T(0);
   for (size_t i: array) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

template <typename T, typename ARRAY>
T norm_inf(const std::function<T(size_t i)>& ith_component, const ARRAY& array) {
   T norm = T(0);
   for (size_t i: array) {
      norm = std::max(norm, std::abs(ith_component(i)));
   }
   return norm;
}

template <typename T, typename ARRAY>
T norm_2_squared(const std::function<T(size_t i)>& ith_component, const ARRAY& array) {
   T norm = T(0);
   for (size_t i: array) {
      const T x_i = ith_component(i);
      norm += x_i * x_i;
   }
   return norm;
}

template <typename T, typename ARRAY>
T norm_2(const std::function<T(size_t /*i*/)>& ith_component, const ARRAY& array) {
   return std::sqrt(norm_2_squared(ith_component, array));
}

template <typename T, typename ARRAY>
T norm(const std::function<T(size_t /*i*/)>& ith_component, const ARRAY& array, Norm norm) {
   // choose the right norm
   if (norm == L1_NORM) {
      return norm_1(ith_component, array);
   }
   else if (norm == L2_NORM) {
      return norm_2(ith_component, array);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared(ith_component, array);
   }
   else if (norm == INF_NORM) {
      return norm_inf(ith_component, array);
   }
   throw std::invalid_argument("The norm is not known");
}

// use && to allow temporaries (such as std::cout or logger DEBUG, WARNING, etc)
template <typename ARRAY, typename STREAM>
void print_vector(STREAM&& stream, const ARRAY& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max()) {
   for (size_t i: range(start, std::min(start + length, x.size()))) {
      stream << x[i] << " ";
   }
   stream << '\n';
}

// check that an array of integers is in increasing order (x[i] <= x[i+1])
template <typename ARRAY>
bool in_increasing_order(const ARRAY& array, size_t length) {
   size_t i = 0;
   while (i < length-1) {
      if (array[i] > array[i + 1]) {
         return false;
      }
      i++;
   }
   return true;
}

#endif // UNO_VECTOR_H
