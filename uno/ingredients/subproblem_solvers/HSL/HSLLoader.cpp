// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HSLLoader.hpp"
#include <cctype>
#include <cstdlib>
#include <string>
#include "tools/Logger.hpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// default library name: libhsl.<platform shared-lib extension> (matches IPOPT's hsllib)
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
      void* raw_symbol(LibraryHandle handle, const char* symbol) {
         return reinterpret_cast<void*>(GetProcAddress(handle, symbol));
      }
#else
      using LibraryHandle = void*;
      LibraryHandle open_library(const char* name) {
         // match upstream IPOPT: resolve now, do not export the HSL symbols globally
         int flags = RTLD_NOW;
#if defined(UNO_HSL_DEEPBIND) && defined(RTLD_DEEPBIND)
         // opt-in (HSL_RUNTIME_DEEPBIND): mirror the jgillis/Ipopt-1 .mod patch that
         // ORs in RTLD_DEEPBIND so libhsl prefers its own symbols. glibc-only.
         flags |= RTLD_DEEPBIND;
#endif
         return dlopen(name, flags);
      }
      void* raw_symbol(LibraryHandle handle, const char* symbol) { return dlsym(handle, symbol); }
#endif

      bool load_attempted = false;
      LibraryHandle hsl_handle = nullptr;

      // Resolve a Fortran symbol trying the manglings IPOPT tries, so the runtime
      // libhsl can have been built by any compiler regardless of how Uno was:
      // base, base_, lower_, lower, UPPER_, UPPER.
      void* resolve_symbol(LibraryHandle handle, const std::string& base) {
         std::string lower = base, upper = base;
         for (char& c: lower) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
         for (char& c: upper) { c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
         const std::string candidates[] = {base, base + "_", lower + "_", lower, upper + "_", upper};
         for (const std::string& candidate: candidates) {
            if (void* symbol = raw_symbol(handle, candidate.c_str())) {
               return symbol;
            }
         }
         return nullptr;
      }

      template <typename FunctionPointer>
      void resolve(LibraryHandle handle, FunctionPointer& function_pointer, const std::string& base) {
         function_pointer = reinterpret_cast<FunctionPointer>(resolve_symbol(handle, base));
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

      resolve(hsl_handle, hsl_ma57id, "ma57id");
      resolve(hsl_handle, hsl_ma57ad, "ma57ad");
      resolve(hsl_handle, hsl_ma57bd, "ma57bd");
      resolve(hsl_handle, hsl_ma57cd, "ma57cd");
      resolve(hsl_handle, hsl_ma57dd, "ma57dd");
      resolve(hsl_handle, hsl_ma57ed, "ma57ed");
      resolve(hsl_handle, hsl_ma27id, "ma27id");
      resolve(hsl_handle, hsl_ma27ad, "ma27ad");
      resolve(hsl_handle, hsl_ma27bd, "ma27bd");
      resolve(hsl_handle, hsl_ma27cd, "ma27cd");
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
