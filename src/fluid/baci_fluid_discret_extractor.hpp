/*-----------------------------------------------------------*/
/*! \file

\brief creates a second discretization as part of the complete discretization for inflow
generation


\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_DISCRET_EXTRACTOR_HPP
#define FOUR_C_FLUID_DISCRET_EXTRACTOR_HPP

#include "baci_config.hpp"

#include "baci_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class FluidDiscretExtractor
  {
   public:
    /*!
   \brief Constructor

   */
    FluidDiscretExtractor(Teuchos::RCP<DRT::Discretization> actdis,  //! parent discretization
        const std::string& condition,  //! condition for separation of domain
        bool yescondition);  //! (unused) bool to distinguish between all nodes having the condition
                             //! and all nodes not having it

    /*!
   \brief Destructor

   */
    virtual ~FluidDiscretExtractor() = default;

    //! get child discretization
    Teuchos::RCP<DRT::Discretization> GetChildDiscretization() { return childdiscret_; }
    //! get node to node coupling in case of periodic boundary conditions (column and row version)
    Teuchos::RCP<std::map<int, std::vector<int>>> GetCoupledColNodesChildDiscretization()
    {
      return col_pbcmapmastertoslave_;
    }
    Teuchos::RCP<std::map<int, std::vector<int>>> GetCoupledRowNodesChildDiscretization()
    {
      return row_pbcmapmastertoslave_;
    }

   private:
    //! the parent discretization
    Teuchos::RCP<DRT::Discretization> parentdiscret_;
    //! the child discretization
    Teuchos::RCP<DRT::Discretization> childdiscret_;
    //! periodic boundary condition: node to node coupling (column and row version)
    Teuchos::RCP<std::map<int, std::vector<int>>> col_pbcmapmastertoslave_;
    Teuchos::RCP<std::map<int, std::vector<int>>> row_pbcmapmastertoslave_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
