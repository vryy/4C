/*-----------------------------------------------------------*/
/*! \file
\brief calc utils for beam interaction framework


\level 3
*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_CALC_UTILS_HPP
#define FOUR_C_BEAMINTERACTION_CALC_UTILS_HPP

#include "4C_config.hpp"

#include "4C_binstrategy_utils.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_mapextractor.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}  // namespace Core::LinAlg

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

namespace Core::Geo
{
  namespace MeshFree
  {
    class BoundingBox;
  }
}  // namespace Core::Geo

namespace BEAMINTERACTION
{
  class CrosslinkingParams;
  class BeamLink;
  namespace UTILS
  {
    /// specific MultiMapExtractor to handle different types of element during beam interaction
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        beam = 0,
        sphere = 1,
        solid = 2
      };

      MAP_EXTRACTOR_VECTOR_METHODS(Beam, beam)
      MAP_EXTRACTOR_VECTOR_METHODS(Sphere, sphere)
      MAP_EXTRACTOR_VECTOR_METHODS(Solid, solid)
    };

    /// class for comparing Core::Elements::Element* (and Core::Nodes::Node*) in a std::set
    /*! -------------------------------------------------------------------------
     * overwrites standard < for pointers, this is necessary to ensure same order
     * of neighboring elements for crosslink check and therefore for random numbers
     * independent of pointer addresses. Without this,
     * simulation with crosslinker is not wrong, but depends randomly on memory
     * allocation, i.e. pointer addresses. Without random numbers, everything is fine
     * with default compare operator
    *  \author J. Eichinger March 2017
     -------------------------------------------------------------------------*/
    class Less
    {
     public:
      template <typename ELEMENT>
      bool operator()(ELEMENT const* first, ELEMENT const* second) const
      {
        return first->Id() < second->Id();
      }
    };

    /*! -------------------------------------------------------------------------
     * class for comparing std::set< std::pair < int, int > >
     *  \author J. Eichinger March 2017
     -------------------------------------------------------------------------*/
    class StdPairComparatorOrderCounts
    {
     public:
      bool operator()(std::pair<int, int> const& lhs, std::pair<int, int> const& rhs) const
      {
        return (lhs.first == rhs.first) ? (lhs.second < rhs.second) : (lhs.first < rhs.first);
      }
    };


    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsBeamElement(Core::Elements::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsRigidSphereElement(Core::Elements::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsBeamNode(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsRigidSphereNode(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsBeamCenterlineNode(Core::Nodes::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void PeriodicBoundaryConsistentDisVector(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<const Core::Geo::MeshFree::BoundingBox> const& pbb,
        Teuchos::RCP<const Core::FE::Discretization> const& discret);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    inline int CalculateNumberOfBeamElementsFromNumberOfNodesOnFilament(
        int const numnodes, int const numnodesperele)
    {
      // from: nodesperfil = nodesperele + ( numele - 1 ) * ( nodesperele - 1 )
      return ((numnodes - numnodesperele) / (numnodesperele - 1.0)) + 1.0;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    std::vector<int> Permutation(int number);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetCurrentElementDis(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetCurrentUnshiftedElementDis(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        Core::Geo::MeshFree::BoundingBox const& pbb, std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    template <typename T>
    void SetFilamentBindingSpotPositions(
        Teuchos::RCP<Core::FE::Discretization> discret, Teuchos::RCP<T> params);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendGhostingForFilamentBspotSetup(
        std::set<int>& relevantfilaments, Teuchos::RCP<Core::FE::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
        std::set<int>& relevantfilaments, std::set<int>& setofrequirednodes,
        Teuchos::RCP<Core::FE::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ComputeFilamentLengthAndSortItsElements(
        std::vector<Core::Elements::Element*>& sortedfilamenteles, std::vector<int> const* nodeids,
        double& filreflength, Teuchos::RCP<Core::FE::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void SetBindingSpotsPositionsOnFilament(
        std::vector<Core::Elements::Element*>& sortedfilamenteles, double start,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, int numbspot,
        double filamentbspotinterval, double tol);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(Core::Elements::Element* ele,
        Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& pbb,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        Core::LinAlg::Matrix<3, 1>& bspotpos, Core::LinAlg::Matrix<3, 3>& bspottriad,
        std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(Core::FE::Discretization const& discret,
        Core::Elements::Element* ele, Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
        Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& pbb,
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        Core::LinAlg::Matrix<3, 1>& bspotpos, Core::LinAlg::Matrix<3, 3>& bspottriad);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool IsDistanceOutOfRange(Core::LinAlg::Matrix<3, 1> const& pos1,
        Core::LinAlg::Matrix<3, 1> const& pos2, double const lowerbound, double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool IsEnclosedAngleOutOfRange(Core::LinAlg::Matrix<3, 1> const& direction1,
        Core::LinAlg::Matrix<3, 1> const& direction2, double const lowerbound,
        double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool DoBeamElementsShareNodes(
        Core::Elements::Element const* const beam, Core::Elements::Element const* const nbbeam);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void fe_assemble_ele_force_stiff_into_system_vector_matrix(
        const Core::FE::Discretization& discret, std::vector<int> const& elegid,
        std::vector<Core::LinAlg::SerialDenseVector> const& elevec,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector> fe_sysvec,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> fe_sysmat);


    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the number of centerline DOFs for a given beam element.
     * @param ele (in) Pointer to the element.
     */
    unsigned int GetNumberOfElementCenterlineDof(const Core::Elements::Element* elerline_gid);

    /**
     * \brief Get the global indices of the centerline DOFs of a beam element.
     * @param discret (in) Pointer to the discretization.
     * @param ele (in) Pointer to the element.
     * @param ele_centerline_dof_indices (out) Vector with global indices of centerline DOFs in the
     * element.
     */
    template <unsigned int n_centerline_dof>
    void GetElementCenterlineGIDIndices(Core::FE::Discretization const& discret,
        const Core::Elements::Element* ele,
        Core::LinAlg::Matrix<n_centerline_dof, 1, int>& centerline_gid);

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Pointer to the discretization.
     * @param ele (in) Pointer to the element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     */
    void GetElementCenterlineDOFIndices(Core::FE::Discretization const& discret,
        const Core::Elements::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
        unsigned int& num_dof);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void AssembleCenterlineDofForceStiffIntoElementForceStiff(
        Core::FE::Discretization const& discret, std::vector<int> const& elegid,
        std::vector<Core::LinAlg::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<Core::LinAlg::SerialDenseVector>* eleforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>* elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assemble a matrix with columns based on centerline DOFs of an element into a matrix
     * with columns based on all DOFs of the element. Example: Mortar coupling matrices as the rows
     * correspond the Lagrange multipliers and the columns correspond to the centerline DOFs.
     * @param discret (in) Pointer to the discretization.
     * @param element (in) Pointer to the element.
     * @param row_matrix_centerlineDOFs (in) Matrix where the columns correspond to the centerline
     * DOFs.
     * @param row_matrix_elementDOFs (out) Matrix where the columns correspond to all Element DOFs
     * (the rest will be 0).
     */
    void AssembleCenterlineDofColMatrixIntoElementColMatrix(Core::FE::Discretization const& discret,
        const Core::Elements::Element* element,
        Core::LinAlg::SerialDenseMatrix const& row_matrix_centerlineDOFs,
        Core::LinAlg::SerialDenseMatrix& row_matrix_elementDOFs);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void ExtractPosDofVecAbsoluteValues(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_absolutevalues);
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void ExtractPosDofVecValues(Core::FE::Discretization const& discret,
        Core::Elements::Element const* ele, Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_values);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceToParentElements(Core::FE::Discretization const& discret,
        Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<Core::LinAlg::SerialDenseVector> const& bspotforce,
        std::vector<Core::LinAlg::SerialDenseVector>& eleforce);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotStiffToParentElements(Core::FE::Discretization const& discret,
        Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& bspotstiff,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceStiffToParentElements(Core::FE::Discretization const& discret,
        Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<Core::LinAlg::SerialDenseVector> const& bspotforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& bspotstiff,
        std::vector<Core::LinAlg::SerialDenseVector>& eleforce,
        std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void SetupEleTypeMapExtractor(Teuchos::RCP<const Core::FE::Discretization> const& discret,
        Teuchos::RCP<Core::LinAlg::MultiMapExtractor>& eletypeextractor);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void UpdateDofMapOfVector(Teuchos::RCP<Core::FE::Discretization> discret,
        Teuchos::RCP<Epetra_Vector>& dofmapvec, Teuchos::RCP<Epetra_Vector> old = Teuchos::null);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    long long CantorPairing(std::pair<int, int> const& pair);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    std::pair<int, int> CantorDePairing(long long z);

    //! convert element @p ele to bin content type
    Core::Binstrategy::Utils::BinContentType ConvertElementToBinContentType(
        const Core::Elements::Element* ele);

  }  // namespace UTILS
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
