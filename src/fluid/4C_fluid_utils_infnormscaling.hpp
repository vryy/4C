// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_UTILS_INFNORMSCALING_HPP
#define FOUR_C_FLUID_UTILS_INFNORMSCALING_HPP


#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class MapExtractor;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace FLD
{
  namespace Utils
  {
    class FluidInfNormScaling
    {
     public:
      //! constructor
      FluidInfNormScaling(Core::LinAlg::MapExtractor& mapextractor);

      //! destructor
      virtual ~FluidInfNormScaling() = default;

      //! perform infnorm-scaling of linear system
      void scale_system(
          Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Core::LinAlg::Vector<double>& b);

      //! perform un-scaling of solution (and the system, just to be on the safe side)
      void unscale_solution(Teuchos::RCP<Core::LinAlg::SparseOperator> matrix,
          Core::LinAlg::Vector<double>& x, Core::LinAlg::Vector<double>& b);

     private:
      //! processor id
      const int myrank_;

      //! Extractor for splitting of velocity and pressure dofs
      Core::LinAlg::MapExtractor& velpressplitter_;

      Teuchos::RCP<Core::LinAlg::Vector<double>> srowsum_;
      Teuchos::RCP<Core::LinAlg::Vector<double>> scolsum_;
      Teuchos::RCP<Core::LinAlg::Vector<double>> prowsum_;
      Teuchos::RCP<Core::LinAlg::Vector<double>> pcolsum_;

      // flags
      const bool leftscale_momentum_;
      const bool leftscale_continuity_;

    };  // class FluidInfNormScaling

  }  // namespace Utils
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
