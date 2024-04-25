/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for non-conforming coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NONCONFORMING_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NONCONFORMING_HPP

#include "4C_config.hpp"

#include "4C_inpar_bio.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_base.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}  // namespace DRT
namespace CORE::LINALG
{
  class SerialDenseVector;
}
namespace POROMULTIPHASESCATRA
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNonConforming : public PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplNonConforming(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

   protected:
    //! Evaluate the 1D-3D coupling
    void Evaluate(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override;

    //! set-up of linear system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const Epetra_Map> dbcmap_art,
        Teuchos::RCP<const Epetra_Map> dbcmap_art_with_collapsed);

    //! setup the strategy
    void Setup() override;

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    virtual void EvaluateAdditionalLinearizationofIntegratedDiam() = 0;

    //! FE-assemble into global force and stiffness matrix
    virtual void FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid,
        const int& ele2gid, const double& integrated_diam,
        std::vector<CORE::LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs);

    //! set flag if varying diameter has to be calculated
    virtual void SetVaryingDiamFlag();

    //! print out the coupling method
    void PrintOutCouplingMethod() const override;

    //! the parameters
    const Teuchos::ParameterList& couplingparams_;

    //! name of the condition
    const std::string condname_;

    //! have the managers been set?
    bool porofluidmanagersset_;

    //! has Setup() been called
    bool issetup_;

    //! is it a pure fluid problem
    bool porofluidprob_;

    //! does it have a varying (by-function)-diameter
    bool has_varying_diam_;

    //! flag to determine if 'free-hanging' elements should be deleted
    bool delete_free_hanging_eles_;

    //! small connected components whose size is smaller than this fraction of the overall network
    //! size are additionally deleted (even if they have a Dirichlet boundary condition)
    double delete_free_hanging_eles_threshold_;

    //! interacting pairs of artery and continuous-discretization elements
    std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>
        coupl_elepairs_;

    //! vector with 1D coupling nodes for ntp coupling
    std::vector<int> couplingnodes_ntp_;

    //! type of coupling method
    INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod coupling_method_;

    //! phinp for artery dis
    Teuchos::RCP<const Epetra_Vector> phinp_art_;

    //! coupling matrix (FE)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> FEmat_;

   private:
    //! check if initial fields on coupled DOFs are equal
    //!  \note not performed here since penalty approach will force solution to be
    //!        equal anyway
    void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override{};

    //! access artery (1D) dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    //! access full dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const override;

    //! Setup global vector
    void SetupVector(Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
        Teuchos::RCP<const Epetra_Vector> vec_art) override;

    void ExtractSingleFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont,
        Teuchos::RCP<const Epetra_Vector>& vec_art) override;

    //! set solution vector of single fields
    void SetSolutionVectors(Teuchos::RCP<const Epetra_Vector> phinp_cont,
        Teuchos::RCP<const Epetra_Vector> phin_cont,
        Teuchos::RCP<const Epetra_Vector> phinp_art) override;

    //! init the strategy
    void Init() override;

    //! set the element pairs that are close as found by search algorithm
    void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs) override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction() override;

    //! create interaction pairs for line- or surfacebased coupling (GPTS and MP)
    void CreateCouplingPairsLineSurfBased();

    //! create interaction pairs for ntp coupling
    void CreateCouplingPairsNTP();

    //! calculate the volume fraction occupied by blood vessels
    void CalculateBloodVesselVolumeFraction();

    //! get the Id from 1D nodes for coupling
    void GetCouplingIdsfromInput();

    //! evaluate the pairs
    void EvaluateCouplingPairs(
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs);

    //! FE-assemble into global D, M and kappa (MP case)
    void FEAssembleDMKappa(const int& ele1gid, const int& ele2gid,
        const CORE::LINALG::SerialDenseMatrix& D_ele, const CORE::LINALG::SerialDenseMatrix& M_ele,
        const CORE::LINALG::SerialDenseVector& Kappa_ele);

    //! sum D and M into global force and stiffness matrix
    void SumDMIntoGlobalForceStiff(
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs);

    //! invert kappa vector
    void InvertKappa();

    //! return appropriate internal implementation class (acts as a simple factory to create single
    //! pairs)
    static Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
    CreateNewArteryCouplingPair(std::vector<DRT::Element const*> const& ele_ptrs);

    //! set the artery diameter in material to be able to use it on 1D discretization
    virtual void SetArteryDiamInMaterial() = 0;

    //! get the segment lengths of element 'artelegid'
    virtual std::vector<double> GetEleSegmentLengths(const int artelegid) = 0;

    //! reset the integrated diameter vector to zero
    virtual void ResetIntegratedDiamToZero() = 0;

    //! fill Function and Scale vectors as read from input file
    void FillFunctionAndScaleVectors();

    //! set the right-hand side factor for time integration
    void SetTimeFacRhs();

    //! right hand side factor for artery time integration
    double timefacrhs_art_;
    //! right hand side factor for time integration of 2D/3D discretization
    double timefacrhs_cont_;

    //! result of brute force search
    std::map<int, std::set<int>> nearbyelepairs_;

    //! phinp for continuous dis
    Teuchos::RCP<const Epetra_Vector> phinp_cont_;

    //! phin for continuous dis
    Teuchos::RCP<const Epetra_Vector> phin_cont_;

    //! zeros for continuous dis
    Teuchos::RCP<const Epetra_Vector> zeros_cont_;

    //! zeros for artery dis
    Teuchos::RCP<const Epetra_Vector> zeros_art_;

    //! scale and function-vector
    std::vector<std::vector<int>> scale_vec_;
    std::vector<std::vector<int>> funct_vec_;

    //! mortar coupling matrices
    Teuchos::RCP<CORE::LINALG::SparseMatrix> d_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> m_;
    Teuchos::RCP<Epetra_FEVector> kappa_inv_;

    //! penalty parameter
    double pp_;

    //! coupling rhs-vector (FE)
    Teuchos::RCP<Epetra_FEVector> fe_rhs_;
  };
}  // namespace POROMULTIPHASESCATRA


FOUR_C_NAMESPACE_CLOSE

#endif
