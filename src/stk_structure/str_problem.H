#ifndef STR_PROBLEM_H
#define STR_PROBLEM_H

#ifdef STKADAPTIVE

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace STK
{
  class Mesh;
  class MetaMesh;
  class Discretization;

  namespace STR
  {
    class Structure;

    class Problem
    {
    public:

      Problem();

      /// setup all mesh part definitions
      void Setup();

      void Execute();

    private:

      /// mesh definitions
      Teuchos::RCP<MetaMesh> meta_;

      /// the whole mesh
      Teuchos::RCP<Mesh> mesh_;

      /// the one and only discretization
      Teuchos::RCP<STK::Discretization> dis_;

      /// the fluid algorithm
      Teuchos::RCP<STK::STR::Structure> structure_;
    };
  }
}

#endif
#endif
