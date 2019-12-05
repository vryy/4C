/*-----------------------------------------------------------*/
/*! \file
\brief calc utils for beam interaction framework

\maintainer Jonas Eichinger

\level 3
*/
/*-----------------------------------------------------------*/

#ifndef BEAMINTERACTION_CALC_UTILS_H_
#define BEAMINTERACTION_CALC_UTILS_H_

#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_beaminteraction.H"


// forward declaration
class Epetra_FEVector;
class Epetra_Vector;

namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}  // namespace LINALG

namespace DRT
{
  class Node;
  class Element;
  class Discretization;
}  // namespace DRT

namespace GEO
{
  namespace MESHFREE
  {
    class BoundingBox;
  }
}  // namespace GEO

namespace BEAMINTERACTION
{
  class CrosslinkingParams;
  class BeamLink;
  namespace UTILS
  {
    /// specific MultiMapExtractor to handle different types of element during beam interaction
    class MapExtractor : public LINALG::MultiMapExtractor
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

    /// class for comparing DRT::Element* (and DRT::Node*) in a std::set
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
    bool IsBeamElement(DRT::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsRigidSphereElement(DRT::Element const& element);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsBeamNode(DRT::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsRigidSphereNode(DRT::Node const& node);

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsBeamCenterlineNode(DRT::Node const& node);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void PeriodicBoundaryConsistentDisVector(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<const GEO::MESHFREE::BoundingBox> const& pbb,
        Teuchos::RCP<const DRT::Discretization> const& discret);

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
    void GetCurrentElementDis(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp, std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetCurrentUnshiftedElementDis(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp, GEO::MESHFREE::BoundingBox const& pbb,
        std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    template <typename T>
    void SetFilamentBindingSpotPositions(
        Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<T> params);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendGhostingForFilamentBspotSetup(
        std::set<int>& relevantfilaments, Teuchos::RCP<DRT::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
        std::set<int>& relevantfilaments, std::set<int>& setofrequirednodes,
        Teuchos::RCP<DRT::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ComputeFilamentLengthAndSortItsElements(std::vector<DRT::Element*>& sortedfilamenteles,
        std::vector<int> const* nodeids, double& filreflength,
        Teuchos::RCP<DRT::Discretization> discret);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void SetBindingSpotsPositionsOnFilament(std::vector<DRT::Element*>& sortedfilamenteles,
        double start, INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int numbspot,
        double filamentbspotinterval, double tol);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(DRT::Element* ele,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        LINALG::Matrix<3, 1>& bspotpos, LINALG::Matrix<3, 3>& bspottriad,
        std::vector<double>& eledisp);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(DRT::Discretization const& discret, DRT::Element* ele,
        Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        LINALG::Matrix<3, 1>& bspotpos, LINALG::Matrix<3, 3>& bspottriad);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool IsDistanceOutOfRange(LINALG::Matrix<3, 1> const& pos1, LINALG::Matrix<3, 1> const& pos2,
        double const lowerbound, double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool IsEnclosedAngleOutOfRange(LINALG::Matrix<3, 1> const& direction1,
        LINALG::Matrix<3, 1> const& direction2, double const lowerbound, double const upperbound);

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    bool DoBeamElementsShareNodes(DRT::Element const* const beam, DRT::Element const* const nbbeam);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void FEAssembleEleForceStiffIntoSystemVectorMatrix(const DRT::Discretization& discret,
        std::vector<int> const& elegid, std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assembles element coupling contributions into global coupling matrices and force
     * vectors needed for partitioned algorithms
     *
     * The function expects a vector object of length 2 containing discretizations
     *
     */
    void FEAssembleEleForceStiffIntoSystemVectorMatrices(
        const Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>>& discretizations,
        std::vector<int> const& elegid, std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector>& ff, Teuchos::RCP<Epetra_FEVector>& fb,
        Teuchos::RCP<LINALG::SparseMatrix>& cff, Teuchos::RCP<LINALG::SparseMatrix>& cbb,
        Teuchos::RCP<LINALG::SparseMatrix>& cfb, Teuchos::RCP<LINALG::SparseMatrix>& cbf);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Pointer to the discretization.
     * @param ele (in) Pointer to the element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     */
    void GetElementCenterlineDOFIndices(DRT::Discretization const& discret, const DRT::Element* ele,
        std::vector<unsigned int>& ele_centerline_dof_indices, unsigned int& num_dof);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Get the local indices of the centerline DOFs of an element.
     * @param discret (in) Pointer to the fluid or structure discretization.
     * @param ele (in) Pointer to the fluid or beam element.
     * @param ele_centerline_dof_indices (out) Vector with local indices of centerline DOFs in the
     * element.
     * @param num_dof (out) Number total DOFs on the element.
     *
     */
    void GetFBIElementCenterlineDOFIndices(DRT::Discretization const& discret,
        const DRT::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
        unsigned int& num_dof);


    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void AssembleCenterlineDofForceStiffIntoElementForceStiff(DRT::Discretization const& discret,
        std::vector<int> const& elegid,
        std::vector<LINALG::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<LINALG::SerialDenseVector>* eleforce,
        std::vector<std::vector<LINALG::SerialDenseMatrix>>* elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     * \brief Assembles centerline DOF coupling contributions into element coupling matrices and
     * force vectors
     *
     * The function expects a vector object of length 2 containing discretizations
     *
     * \param[in] discretizations vector to the structure and then fluid discretization
     * \param[in] elegid global id for the beam and fluid element
     * \param[in] eleforce_centerlineDOFs force vectors containing contributions from the centerline
     * dofs in the case of beams, otherwise just the classical element DOFs
     * \param[in]
     * elestiffcenterlineDOFs stiffness matrices containing contributions from the centerline dofs
     * in the case of beams, otherwise just the classical element DOFs
     * \param[inout] eleforce
     * element force vector
     * \param[inout] elestiff element stiffness matrix
     *
     */

    void AssembleCenterlineDofForceStiffIntoFBIElementForceStiff(
        const Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>>& discretizations,
        std::vector<int> const& elegid,
        std::vector<LINALG::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<LINALG::SerialDenseVector>* eleforce,
        std::vector<std::vector<LINALG::SerialDenseMatrix>>* elestiff);

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
    void AssembleCenterlineDofColMatrixIntoElementColMatrix(DRT::Discretization const& discret,
        const DRT::Element* element, LINALG::SerialDenseMatrix const& row_matrix_centerlineDOFs,
        LINALG::SerialDenseMatrix& row_matrix_elementDOFs);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void ExtractPosDofVecAbsoluteValues(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_absolutevalues);
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void ExtractPosDofVecValues(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_absolutevalues);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<LINALG::SerialDenseVector> const& bspotforce,
        std::vector<LINALG::SerialDenseVector>& eleforce);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotStiffToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& bspotstiff,
        std::vector<std::vector<LINALG::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceStiffToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<LINALG::SerialDenseVector> const& bspotforce,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& bspotstiff,
        std::vector<LINALG::SerialDenseVector>& eleforce,
        std::vector<std::vector<LINALG::SerialDenseMatrix>>& elestiff);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void SetupEleTypeMapExtractor(Teuchos::RCP<const DRT::Discretization> const& discret,
        Teuchos::RCP<LINALG::MultiMapExtractor>& eletypeextractor);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void UpdateDofMapOfVector(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Vector>& dofmapvec, Teuchos::RCP<Epetra_Vector> old = Teuchos::null);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    long long CantorPairing(std::pair<int, int> const& pair);

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    std::pair<int, int> CantorDePairing(long long z);


  }  // namespace UTILS
}  // namespace BEAMINTERACTION

#endif
