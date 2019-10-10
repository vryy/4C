/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for line-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_LINEBASED_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_LINEBASED_H_

#include "poromultiphase_scatra_artery_coupling_base.H"
#include "../drt_inpar/inpar_bio.H"

#include <Epetra_Time.h>

// forward declaration
class Epetra_FEVector;
namespace DRT
{
  class Element;
  class PoromultiphaseScatraArteryCouplingExchangeManager;
}  // namespace DRT
namespace LINALG
{
  class SerialDenseVector;
}
namespace POROMULTIPHASESCATRA
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplLineBased : public PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplLineBased(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraArtCouplLineBased() {}

    // check if initial fields on coupled DOFs are equal
    void CheckInitialFields(
        Teuchos::RCP<const Epetra_Vector> vec_cont, Teuchos::RCP<const Epetra_Vector> vec_art);

    //! access artery (1D) dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const;

    //! access full dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    //! set-up of global system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
        Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art);

    /// setup global vector
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    void SetupVector(Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art);

    /// extract single field vectors
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    void ExtractSingleFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont, Teuchos::RCP<const Epetra_Vector>& vec_art);

    //! set solution of fields
    void SetSolutionVectors(Teuchos::RCP<const Epetra_Vector> phinp_cont,
        Teuchos::RCP<const Epetra_Vector> phin_cont, Teuchos::RCP<const Epetra_Vector> phinp_art);

    //! set the element pairs that are close as found by search algorithm
    void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs);

    /// init the strategy
    void Init();

    /// setup the strategy
    void Setup();

    /// apply mesh movement (on artery elements)
    void ApplyMeshMovement();

    /// access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction();

   private:
    //! create interaction pairs
    void CreateCouplingPairs();

    //! pre-evaluate the pairs
    void PreEvaluateCouplingPairs();

    //! fill the length not changed through deformation
    void FillUnaffectedArteryLength();

    //! calculate the volume fraction occupied by blood vessels
    void CalculateBloodVesselVolumeFraction();

    //! create gid to seg
    void CreateGIDToSegmentVector();

    //! fill the gid to segment vector
    void FillGIDToSegmentVector(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        std::map<int, std::vector<double>>& gid_to_seglength);

    //! evaluate the pairs
    void EvaluateCouplingPairs(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art);

    //! FE-assemble into force and stiffness matrix
    void FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid, const int& ele2gid,
        std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art);

    //! FE-assemble into D, M and kappa (MP case)
    void FEAssembleDMKappa(const int& ele1gid, const int& ele2gid,
        const LINALG::SerialDenseMatrix& D_ele, const LINALG::SerialDenseMatrix& M_ele,
        const LINALG::SerialDenseVector& Kappa_ele);

    //! sum D and M into global force and stiffness matrix
    void SumDMIntoGlobalForceStiff(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art);

    //! invert kappa
    void InvertKappa();

    //! get the segment lengths of element 'artelegid'
    std::vector<double> GetEleSegmentLengths(const int artelegid);

    //! check for duplicates
    bool IsDuplicateSegment(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
            possible_duplicate);

    //! check for duplicates
    bool IsIdenticalSegment(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        const int& ele1gid, const double& etaA, const double& etaB, int& elepairID);

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
    CreateNewArteryCouplingPair(std::vector<DRT::Element const*> const& ele_ptrs);

    //! print output of mesh tying pairs
    void OutputSummary() const;

    //! fill the vectors
    void FillFunctionAndScaleVectors();

    //! set the right-hand side factor for time integration
    void SetTimeFacRhs();

    //! print out the coupling method
    void PrintOutCouplingMethod() const;

    //! the parameters
    const Teuchos::ParameterList& couplingparams_;

    //! have the managers been set?
    bool porofluidmanagersset_;

    //! has Setup() been called
    bool issetup_;

    //! is it a pure fluid problem
    bool porofluidprob_;

    //! right hand side factor for artery time integration
    double timefacrhs_art_;
    //! right hand side factor for time integration of 2D/3D discretization
    double timefacrhs_cont_;

    //! result of brute force search
    std::map<int, std::set<int>> nearbyelepairs_;

    //! interacting pairs of artery and continuous-discretization elements
    std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>
        coupl_elepairs_;

    //! phinp for continuous dis
    Teuchos::RCP<const Epetra_Vector> phinp_cont_;

    //! phin for continuous dis
    Teuchos::RCP<const Epetra_Vector> phin_cont_;

    //! phinp for artery dis
    Teuchos::RCP<const Epetra_Vector> phinp_art_;

    //! zeros for continuous dis
    Teuchos::RCP<const Epetra_Vector> zeros_cont_;

    //! zeros for artery dis
    Teuchos::RCP<const Epetra_Vector> zeros_art_;

    //! scale and function-vector
    std::vector<std::vector<int>> scale_vec_;
    std::vector<std::vector<int>> funct_vec_;

    //! type of coupling method
    INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod coupling_method_;

    //! mortar coupling matrices
    Teuchos::RCP<LINALG::SparseMatrix> D_;
    Teuchos::RCP<LINALG::SparseMatrix> M_;
    Teuchos::RCP<Epetra_FEVector> kappaInv_;

    //! penalty parameter
    double pp_;

    //! maximum number of segments per artery element
    int maxnumsegperartele_;

    //! length of artery elements unaffected by deformation
    Teuchos::RCP<Epetra_FEVector> unaffected_seg_lengths_artery_;

    //! length of artery elements in current configuration
    Teuchos::RCP<Epetra_FEVector> current_seg_lengths_artery_;

    //! volume fraction of blood vessels (for output)
    Teuchos::RCP<Epetra_Vector> bloodvesselvolfrac_;

    //! gid to segment: stores [GID; [eta_a eta_b]_1, [eta_a eta_b]_2, ..., [eta_a eta_b]_n]
    //  of artery elements in column format, i.e. fully overlapping
    std::map<int, std::vector<double>> gid_to_segment_;

    //! gid to segment length: stores [GID; seglength_1, seglength_2, ..., seglength_n]
    //  of artery elements in column format, i.e. fully overlapping (only used for
    //  porofluid-problems)
    std::map<int, std::vector<double>> gid_to_seglength_;

    //! timer for measurement of duration of brute force search
    Epetra_Time timersearch_;

    //! coupling matrix (FE)
    Teuchos::RCP<LINALG::SparseMatrix> FEmat_;

    //! coupling rhs-vector (FE)
    Teuchos::RCP<Epetra_FEVector> FErhs_;
  };

}  // namespace POROMULTIPHASESCATRA

#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_LINEBASED_H_ */
