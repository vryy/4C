/*!----------------------------------------------------------------------
\file drt_utils_parallel.H

\brief A collection of helper methods for namespace DRT

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_PARALLEL_H
#define UTILS_PARALLEL_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "../drt_adapter/adapter_coupling.H"


namespace DRT
{
/// forward declarations
class Discretization;

/// collection of useful and important but less central stuff
namespace UTILS
{
  //! Given a nodal row and coloumn map, this routines creates the corresponding element maps and call FillComplete()
  void RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    );

  /*! \brief Redistribute discretizations in input vector using BinningStrategy
   *
  \param vector_of_discretizations (in) : vector containing RCPs to discretizations */
  void RedistributeDiscretizationsByBinning(
      const std::vector<Teuchos::RCP<DRT::Discretization> >& vector_of_discretizations,
      bool revertextendedghosting = false);

  /*!
    \brief Ghost the discretization handed in to this method on all procs

    \param distobeghosted (in) : RCP to the discretisation which is ghosted after execution of this method

    \return void */
  void GhostDiscretizationOnAllProcs(
      const Teuchos::RCP<DRT::Discretization> distobeghosted );

  /*! \brief Redistribute nodes matching to another discretization.
   *
   * \note The discretization serving as template and the discretization
   *       supposed to be redistributed have to be matching! That means,
   *       nodal positions have to coincide.
   *
  \param dis_template        (in) : discretization with parallel distr., serving as template
  \param dis_to_redistribute (in) : discretization which is redistributed matching to dis_template  */
  void MatchNodalDistributionOfMatchingDiscretizations(
    DRT::Discretization&     dis_template,
    DRT::Discretization&     dis_to_redistribute );

  /*! \brief Redistribute elements matching to another discretization.
   *
   *    The difference to \ref MatchNodalDistributionOfMatchingDiscretizations
   *    is, that a equal distribution of nodes may lead to unequal distribution
   *    of elements, depending on the algorithm, determining which elements are
   *    owned and which elements are ghosted. This method makes sure, that element
   *    ownerships match and nodes are distributed accordingly.
   *
   * \note The discretization serving as template and the discretization
   *       supposed to be redistributed have to be matching! That means,
   *       nodal positions have to coincide.
   *
  \param dis_template        (in) : discretization with parallel distr., serving as template
  \param dis_to_redistribute (in) : discretization which is redistributed matching to dis_template  */
  void MatchElementDistributionOfMatchingDiscretizations(
    DRT::Discretization&     dis_template,
    DRT::Discretization&     dis_to_redistribute );

  /// print information about parallel mesh distribution
  void PrintParallelDistribution(const DRT::Discretization& dis);

  /*! \brief Get a column vector made of a row vector.

  Using this method, a reference to a vector is returned.
  If the vector is supplied in dis.DofColMap() the vector itself will be returned.
  If the vector is NOT supplied in dis.DofColMap(), but in dis.DofRowMap(),
  a new vector with column map is allocated and the supplied vector is exported to it and returned.

  \note The very same functionality is used in Discretization::SetState()!

  \param name (in): discretization
  \param state (in): vector of some data  */
  Teuchos::RCP<const Epetra_Vector> GetColVersionOfRowVector(
      const Teuchos::RCP<const DRT::Discretization> dis,
      const Teuchos::RCP<const Epetra_Vector> state,
      const int nds = 0);


  /// recompute nodecolmap of standard discretization to include all nodes as of subdicretization
  Teuchos::RCP<Epetra_Map> ComputeNodeColMap
  (
      const Teuchos::RCP<DRT::Discretization> sourcedis,  ///< standard discretization we want to redistribute
      const Teuchos::RCP<DRT::Discretization> subdis      ///< subdiscretization prescribing ghosting
  );

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_PARALLEL_H
