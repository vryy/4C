/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef UTILS_PARALLEL_H
#define UTILS_PARALLEL_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace UTILS
  {
    /*! \brief Redistribute discretizations in input vector using BinningStrategy
     *
    \param vector_of_discretizations (in) : vector containing RCPs to discretizations */
    void RedistributeDiscretizationsByBinning(
        const std::vector<Teuchos::RCP<DRT::Discretization>>& vector_of_discretizations,
        bool revertextendedghosting = false);

    /*!
      \brief Ghost the discretization handed in to this method on all procs

      \param distobeghosted (in) : RCP to the discretisation which is ghosted after execution of
      this method

      \return void */
    void GhostDiscretizationOnAllProcs(const Teuchos::RCP<DRT::Discretization> distobeghosted);

    /*! \brief Redistribute nodes matching to another discretization.
     *
     * \note The discretization serving as template and the discretization
     *       supposed to be redistributed have to be matching! That means,
     *       nodal positions have to coincide.
     *
    \param dis_template        (in) : discretization with parallel distr., serving as template
    \param dis_to_redistribute (in) : discretization which is redistributed matching to dis_template
  */
    void MatchNodalDistributionOfMatchingDiscretizations(
        DRT::Discretization& dis_template, DRT::Discretization& dis_to_redistribute);

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
    \param dis_to_redistribute (in) : discretization which is redistributed matching to dis_template
  */
    void MatchElementDistributionOfMatchingDiscretizations(
        DRT::Discretization& dis_template, DRT::Discretization& dis_to_redistribute);

    /*! \brief Redistribute conditioned elements matching other conditioned elements.
     *
     * Redistribute the elements of discretization dis_to_redistribute, wearing the
     * condition condname_redistribute to match the parallel distribution of the
     * elements of discretization dis_template, wearing the condition condname_redistribute.
     *
     * \note The elements serving as template and the elements supposed to be redistributed
     * have to be matching! That means, nodal positions have to coincide.
     *
    \param dis_template        (i) : discretization with parallel distr., serving as template
    \param dis_to_redistribute (i) : discretization which is redistributed matching to dis_template
    \param condname_template   (i) : condition on template dis with template distribution
    \param condname_redistribute (i) : condition on elements to be redistributed */
    void MatchElementDistributionOfMatchingConditionedElements(DRT::Discretization& dis_template,
        DRT::Discretization& dis_to_redistribute, const std::string& condname_template,
        const std::string& condname_redistribute);

    /// print information about parallel mesh distribution
    void PrintParallelDistribution(const DRT::Discretization& dis);

    /*! \brief Get a column vector made of a row vector.

    Using this method, a reference to a vector is returned.
    If the vector is supplied in dis.DofColMap() the vector itself will be returned.
    If the vector is NOT supplied in dis.DofColMap(), but in dis.DofRowMap(),
    a new vector with column map is allocated and the supplied vector is exported to it and
    returned.

    \note The very same functionality is used in Discretization::SetState()!

    \param name (in): discretization
    \param state (in): vector of some data  */
    Teuchos::RCP<const Epetra_Vector> GetColVersionOfRowVector(
        const Teuchos::RCP<const DRT::Discretization> dis,
        const Teuchos::RCP<const Epetra_Vector> state, const int nds = 0);


    /// recompute nodecolmap of standard discretization to include all nodes as of subdicretization
    Teuchos::RCP<Epetra_Map> ComputeNodeColMap(
        const Teuchos::RCP<DRT::Discretization>
            sourcedis,  ///< standard discretization we want to redistribute
        const Teuchos::RCP<DRT::Discretization> subdis  ///< subdiscretization prescribing ghosting
    );

    /*! \brief Fill processor local row and col vectors with element ids fitting the desired
    parallel distribution.

    \author Andreas Rauch
    \date   10/16

    \param dis_template        (i): template for parallel distribution
    \param dis_to_redistribute (i): discretization supposed to be distributed matching template
    \param row_id_vec_to_fill  (o): on exit this vector contains the matched element row gids
    \param col_id_vec_to_fill  (o): on exit this vector contains the matched element col gids  */
    void MatchElementRowColDistribution(const DRT::Discretization& dis_template,
        const DRT::Discretization& dis_to_redistribute, std::vector<int>& row_id_vec_to_fill,
        std::vector<int>& col_id_vec_to_fill);

    /*! \brief Fill processor local row and col vectors with node ids fitting the desired parallel
    distribution.

    \author Andreas Rauch
    \date   10/16

    \param dis_template        (i): template for parallel distribution
    \param dis_to_redistribute (i): discretization supposed to be distributed matching template
    \param row_id_vec_to_fill  (o): on exit this vector contains the matched node row gids
    \param col_id_vec_to_fill  (o): on exit this vector contains the matched node col gids  */
    void MatchNodalRowColDistribution(const DRT::Discretization& dis_template,
        const DRT::Discretization& dis_to_redistribute, std::vector<int>& row_id_vec_to_fill,
        std::vector<int>& col_id_vec_to_fill);

    /// \brief redistribute a map in accordance with a already redistributed reference map
    /** \note This function is supposed to work for row as well as for column maps.
     *
     *  \param(in) ref_red_map : already redistributed reference map. This map
     *                           can be larger than the non-redistributed map but has to
     *                           contain the non-redistributed map as a subset.
     *  \param(in) unred_map   : original non-redistributed map
     *  \param(out) red_map    : redistributed map, following the same parallel
     *                           distribution as the reference map
     *
     *  \author hiermeier \date 01/18 */
    void RedistributeInAccordanceWithReference(const Epetra_Map& ref_red_map,
        const Epetra_Map& unred_map, Teuchos::RCP<Epetra_Map>& red_map);

    /// \brief redistribute a map in accordance with a already redistributed reference map
    /** This is a short function wrapper which first generates a new Epetra_Map
     *  RCP before the magic is happening. See the called function for more info.
     *
     *  \return the redistributed version of the non-redistributed map in accordance
     *  with the already redistributed reference map.
     *
     *  \author hiermeier \date 01/18 */
    Teuchos::RCP<Epetra_Map> RedistributeInAccordanceWithReference(
        const Epetra_Map& ref_red_map, const Epetra_Map& unred_map);

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_PARALLEL_H
