// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HSLLoader.hpp"
#include <cstdlib>
#include "tools/Logger.hpp"
#include "fortran_interface.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// stringize the Fortran-mangled symbol name (e.g. FC_GLOBAL(ma57id, MA57ID) -> "ma57id_")
#define UNO_STRINGIZE_(x) #x
#define UNO_STRINGIZE(x) UNO_STRINGIZE_(x)

#if defined(_WIN32)
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.dll"
#elif defined(__APPLE__)
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.dylib"
#else
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.so"
#endif

namespace uno {
   ma57id_fp hsl_ma57id = nullptr;
   ma57ad_fp hsl_ma57ad = nullptr;
   ma57bd_fp hsl_ma57bd = nullptr;
   ma57cd_fp hsl_ma57cd = nullptr;
   ma57dd_fp hsl_ma57dd = nullptr;
   ma57ed_fp hsl_ma57ed = nullptr;
   ma27id_fp hsl_ma27id = nullptr;
   ma27ad_fp hsl_ma27ad = nullptr;
   ma27bd_fp hsl_ma27bd = nullptr;
   ma27cd_fp hsl_ma27cd = nullptr;

   namespace {
#ifdef _WIN32
      using LibraryHandle = HMODULE;
      LibraryHandle open_library(const char* name) { return LoadLibraryA(name); }
      void* load_symbol(LibraryHandle handle, const char* symbol) {
         return reinterpret_cast<void*>(GetProcAddress(handle, symbol));
      }
#else
      using LibraryHandle = void*;
      LibraryHandle open_library(const char* name) { return dlopen(name, RTLD_LAZY | RTLD_GLOBAL); }
      void* load_symbol(LibraryHandle handle, const char* symbol) { return dlsym(handle, symbol); }
#endif

      bool load_attempted = false;
      LibraryHandle hsl_handle = nullptr;

      template <typename FunctionPointer>
      void resolve(LibraryHandle handle, FunctionPointer& function_pointer, const char* symbol) {
         function_pointer = reinterpret_cast<FunctionPointer>(load_symbol(handle, symbol));
      }
   } // anonymous namespace

   bool load_hsl_library(const std::string& library_name) {
      if (load_attempted) {
         return hsl_handle != nullptr;
      }
      load_attempted = true;

      std::string name = library_name;
      if (name.empty()) {
         if (const char* env = std::getenv("UNO_HSL_LIBRARY")) {
            name = env;
         }
      }
      if (name.empty()) {
         name = UNO_HSL_DEFAULT_LIBRARY;
      }

      hsl_handle = open_library(name.c_str());
      if (hsl_handle == nullptr) {
         DEBUG << "Uno: could not load the HSL library '" << name << "' at runtime\n";
         return false;
      }
      DEBUG << "Uno: loaded the HSL library '" << name << "' at runtime\n";

      resolve(hsl_handle, hsl_ma57id, UNO_STRINGIZE(FC_GLOBAL(ma57id, MA57ID)));
      resolve(hsl_handle, hsl_ma57ad, UNO_STRINGIZE(FC_GLOBAL(ma57ad, MA57AD)));
      resolve(hsl_handle, hsl_ma57bd, UNO_STRINGIZE(FC_GLOBAL(ma57bd, MA57BD)));
      resolve(hsl_handle, hsl_ma57cd, UNO_STRINGIZE(FC_GLOBAL(ma57cd, MA57CD)));
      resolve(hsl_handle, hsl_ma57dd, UNO_STRINGIZE(FC_GLOBAL(ma57dd, MA57DD)));
      resolve(hsl_handle, hsl_ma57ed, UNO_STRINGIZE(FC_GLOBAL(ma57ed, MA57ED)));
      resolve(hsl_handle, hsl_ma27id, UNO_STRINGIZE(FC_GLOBAL(ma27id, MA27ID)));
      resolve(hsl_handle, hsl_ma27ad, UNO_STRINGIZE(FC_GLOBAL(ma27ad, MA27AD)));
      resolve(hsl_handle, hsl_ma27bd, UNO_STRINGIZE(FC_GLOBAL(ma27bd, MA27BD)));
      resolve(hsl_handle, hsl_ma27cd, UNO_STRINGIZE(FC_GLOBAL(ma27cd, MA27CD)));
      return true;
   }

   bool ma57_symbols_available() {
      load_hsl_library();
      return hsl_ma57id && hsl_ma57ad && hsl_ma57bd && hsl_ma57cd && hsl_ma57dd && hsl_ma57ed;
   }

   bool ma27_symbols_available() {
      load_hsl_library();
      return hsl_ma27id && hsl_ma27ad && hsl_ma27bd && hsl_ma27cd;
   }
} // namespace

// Mirrors the symbol the linked coinhsl meta-library exports; the
// SymmetricIndefiniteLinearSolverFactory uses it to gate MA27/MA57.
extern "C" bool LIBHSL_isfunctional() {
   return uno::ma57_symbols_available() || uno::ma27_symbols_available();
}
