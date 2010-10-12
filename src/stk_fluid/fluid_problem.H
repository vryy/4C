#ifndef FLUID_PROBLEM_H
#define FLUID_PROBLEM_H

#ifdef STKADAPTIVE

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace STK
{
  class Mesh;
  class MetaMesh;
  class Discretization;

  namespace FLD
  {
    class Fluid;

    /// central fluid problem definition
    class Problem
    {
    public:

      Problem();

      /// setup all mesh part definitions
      void Setup();

      /// execute the whole thing
      void Execute();

    private:

      /// mesh definitions
      Teuchos::RCP<MetaMesh> meta_;

      /// the whole mesh
      Teuchos::RCP<Mesh> mesh_;

      /// the one and only discretization
      Teuchos::RCP<STK::Discretization> dis_;

      /// the fluid algorithm
      Teuchos::RCP<STK::FLD::Fluid> fluid_;
    };
  }
}

#endif
#endif
