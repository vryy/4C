/*----------------------------------------------------------------------*/
/*! \file

\brief utils class for use of binning strategy

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#ifndef BINNING_STRATEGY_UTILS_H
#define BINNING_STRATEGY_UTILS_H

#include <set>
#include <vector>

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Vector;
class Epetra_Map;
namespace DRT
{
  class Discretization;
  class Element;
  class Node;
}  // namespace DRT

namespace BINSTRATEGY
{
  namespace UTILS
  {
    /*!
     * \brief Type of elements assigned to a bin
     */
    enum BinContentType
    {
      Scatra,       ///< scatra element
      Fluid,        ///< fluid element
      BELE3,        ///< bele3 element
      Beam,         ///< beam element
      RigidSphere,  ///< rigid sphere element
      Solid,  ///< solid element (all elements derived from So_base, if more distinction is needed,
              ///< split this type)

      enumsize  ///< NOTE: This entry gives the size of the enum and is used somewhere else in baci
                ///< therefore put it ALWAYS in last place
    };

    /*!
     * \brief Extend ghosting of discretization according to extended element col map
     *
     * @param[in] discret discretization
     * @param[in] extendedelecolmap extended element column map
     * @param[in] assigndegreesoffreedom assign degrees of freedom, enters FillComplete call
     * @param[in] initelements init elements, enters FillComplete call
     * @param[in] doboundaryconditions do boundary conditions, enters FillComplete call
     */
    void ExtendDiscretizationGhosting(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Map> const& extendedelecolmap, bool assigndegreesoffreedom,
        bool initelements, bool doboundaryconditions);

    /*!
     * \brief convert element to bin content type
     *
     * @param[in] eleptr
     *
     * @return bin content type
     */
    BinContentType ConvertElementToBinContentType(DRT::Element const* const eleptr);

    /*!
     * \brief communicate elements that get a new owner
     *
     * @param[in] discret discretization
     * @param[in] toranktosendeles key: new owner std::vector: elements that are sended to new owner
     */
    void CommunicateElements(Teuchos::RCP<DRT::Discretization>& discret,
        std::map<int, std::vector<DRT::Element*>> const& toranktosendeles);

    /*!
     * \brief communicate distribution of transferred elements to bins
     *
     * @param[in] discret discretization
     * @param[in] toranktosendbinids sent bin ids (std::vector) to rank (key)
     * @param[out] bintorowelemap
     */
    void CommunicateDistributionOfTransferredElementsToBins(
        Teuchos::RCP<DRT::Discretization>& discret,
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
    void GetCurrentNodePos(Teuchos::RCP<const DRT::Discretization> const discret,
        DRT::Node const* node, Teuchos::RCP<const Epetra_Vector> const disnp, double* currpos);

  }  // namespace UTILS
}  // namespace BINSTRATEGY

#endif
