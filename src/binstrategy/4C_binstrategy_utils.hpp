/*----------------------------------------------------------------------*/
/*! \file

\brief utils class for use of binning strategy

\level 2

*----------------------------------------------------------------------*/


#ifndef FOUR_C_BINSTRATEGY_UTILS_HPP
#define FOUR_C_BINSTRATEGY_UTILS_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::Nodes
{
  class Node;
}

namespace BINSTRATEGY::UTILS
{
  /*!
   * \brief Type of elements assigned to a bin
   */
  enum class BinContentType : int
  {
    Scatra,       ///< scatra element
    Fluid,        ///< fluid element
    BELE3,        ///< bele3 element
    Beam,         ///< beam element
    RigidSphere,  ///< rigid sphere element
    Solid  ///< solid element (all elements derived from So_base, if more distinction is needed,
           ///< split this type)
  };

  /*!
   * \brief Extend ghosting of discretization according to extended element col map
   *
   * @param[in] discret discretization
   * @param[in] extendedelecolmap extended element column map
   * @param[in] assigndegreesoffreedom assign degrees of freedom, enters fill_complete call
   * @param[in] initelements init elements, enters fill_complete call
   * @param[in] doboundaryconditions do boundary conditions, enters fill_complete call
   */
  void ExtendDiscretizationGhosting(Teuchos::RCP<Core::FE::Discretization> discret,
      Teuchos::RCP<Epetra_Map> const& extendedelecolmap, bool assigndegreesoffreedom,
      bool initelements, bool doboundaryconditions);

  /*!
   * \brief convert element to bin content type
   *
   * @param[in] eleptr
   *
   * @return bin content type
   */
  BinContentType ConvertElementToBinContentType(Core::Elements::Element const* const eleptr);

  /*!
   * \brief communicate elements that get a new owner
   *
   * @param[in] discret discretization
   * @param[in] toranktosendeles key: new owner std::vector: elements that are sended to new owner
   */
  void CommunicateElements(Teuchos::RCP<Core::FE::Discretization>& discret,
      std::map<int, std::vector<Core::Elements::Element*>> const& toranktosendeles);

  /*!
   * \brief communicate distribution of transferred elements to bins
   *
   * @param[in] discret discretization
   * @param[in] toranktosendbinids sent bin ids (std::vector) to rank (key)
   * @param[out] bintorowelemap
   */
  void CommunicateDistributionOfTransferredElementsToBins(
      Teuchos::RCP<Core::FE::Discretization>& discret,
      std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosendbinids,
      std::map<int, std::set<int>>& bintorowelemap);

  /*!
   * \brief get current position of node
   *
   * @param[in] discret discretization
   * @param[in] node node
   * @param[in] disnp current displacement state
   * @param[out] currpos current position of node
   */
  void GetCurrentNodePos(Teuchos::RCP<const Core::FE::Discretization> const discret,
      Core::Nodes::Node const* node, Teuchos::RCP<const Epetra_Vector> const disnp,
      double* currpos);

}  // namespace BINSTRATEGY::UTILS


FOUR_C_NAMESPACE_CLOSE

#endif
