/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_BASE_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_BASE_H_

#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace POROMULTIPHASESCATRA
{
  //! base class for coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplBase
  {
   public:
    // constructor
    PoroMultiPhaseScaTraArtCouplBase(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& meshtyingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraArtCouplBase(){};

    /// access to full map
    const Teuchos::RCP<const Epetra_Map>& FullMap() const;

    /// global extractor
    const Teuchos::RCP<LINALG::MultiMapExtractor>& GlobalExtractor() const;

    // check if initial fields on coupled DOFs are equal
    virtual void CheckInitialFields(
        Teuchos::RCP<const Epetra_Vector> vec_cont, Teuchos::RCP<const Epetra_Vector> vec_art) = 0;

    //! access artery (1D) dof row map
    virtual Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const = 0;

    //! access full dof row map
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const = 0;

    //! print out the coupling method
    virtual void PrintOutCouplingMethod() const = 0;

    //! set-up of global system of equations of coupled problem
    virtual void SetupSystem(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
        Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art) = 0;

    //! set solution of fields
    virtual void SetSolutionVectors(Teuchos::RCP<const Epetra_Vector> phinp_cont,
        Teuchos::RCP<const Epetra_Vector> phin_cont, Teuchos::RCP<const Epetra_Vector> phinp_art);

    //! set the element pairs that are close as found by search algorithm
    virtual void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs);

    /// setup global vector
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    virtual void SetupVector(Teuchos::RCP<Epetra_Vector> vec,
        Teuchos::RCP<const Epetra_Vector> vec_cont, Teuchos::RCP<const Epetra_Vector> vec_art) = 0;

    /// extract single field vectors
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    virtual void ExtractSingleFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
        Teuchos::RCP<const Epetra_Vector>& vec_cont,
        Teuchos::RCP<const Epetra_Vector>& vec_art) = 0;

    /// init the strategy
    virtual void Init() = 0;

    /// setup the strategy
    virtual void Setup() = 0;

    /// apply mesh movement (on artery elements)
    virtual void ApplyMeshMovement() = 0;

    //! return blood vessel volume fraction inside each 2D/3D element
    virtual Teuchos::RCP<const Epetra_Vector> BloodVesselVolumeFraction() = 0;

   protected:
    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    //! artery (1D) discretization
    Teuchos::RCP<DRT::Discretization> arterydis_;
    //! continous field (2D, 3D) discretization
    Teuchos::RCP<DRT::Discretization> contdis_;

    //! coupled dofs of the single fields
    std::vector<int> coupleddofs_art_;
    std::vector<int> coupleddofs_cont_;
    int num_coupled_dofs_;

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! global extractor
    Teuchos::RCP<LINALG::MultiMapExtractor> globalex_;

    //! myrank
    const int myrank_;

    //! so far, it is assumed that artery elements always follow the deformation of the underlying
    //! porous medium
    //  hence, we actually have to evalute them in current configuration
    //  if this flag is set to true, artery elements will not move and are evaluated in reference
    //  configuration
    bool evaluate_in_ref_config_;

   private:
    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;
  };

}  // namespace POROMULTIPHASESCATRA



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_BASE_H_ */
