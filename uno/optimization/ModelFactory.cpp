// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifdef WITH_AMPL

#include "ModelFactory.hpp"
#include "EqualityConstrainedModel.hpp"
#include "ScaledModel.hpp"
<<<<<<< HEAD
#include "interfaces/AMPL/AMPLModel.hpp"
=======
#include "BoundRelaxedModel.hpp"
>>>>>>> a561f0d13e470766366b6ecfc9a83501ad3044a1
#include "preprocessing/Scaling.hpp"

// note: transfer ownership of the pointer
std::unique_ptr<Model> ModelFactory::reformulate(std::unique_ptr<Model> model, Iterate& first_iterate, const Options& options) {
   // optional: scale the problem using the evaluations at the first iterate
   if (options.get_string("scale_functions") == "yes") {
      model = std::make_unique<ScaledModel>(std::move(model), first_iterate, options);
   }

   // if an equality-constrained problem is required (e.g. barrier or AL), reformulate the model with slacks
   if (options.get_string("subproblem") == "barrier") {
<<<<<<< HEAD
      // transfer ownership of the pointer
=======
      // introduce slacks to obtain equality constraints
>>>>>>> a561f0d13e470766366b6ecfc9a83501ad3044a1
      model = std::make_unique<EqualityConstrainedModel>(std::move(model));
      first_iterate.set_number_variables(model->number_variables);
   }
   return model;
}

#endif // WITH_AMPL