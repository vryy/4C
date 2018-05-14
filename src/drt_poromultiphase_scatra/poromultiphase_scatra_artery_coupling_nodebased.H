/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_artery_coupling_nodebased.H

 \brief base algorithm for node-based coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_H_

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include "../linalg/linalg_utils.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace POROMULTIPHASESCATRA
{

  //! Node based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplNodeBased
  {

  public:

    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplNodeBased(
        Teuchos::RCP<DRT::Discretization>  arterydis,
        Teuchos::RCP<DRT::Discretization>  contdis,
        const Teuchos::ParameterList&      meshtyingparams,
        const std::string&                 condname,
        const std::string&                 artcoupleddofname,
        const std::string&                 contcoupleddofname
   );

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraArtCouplNodeBased(){}

    /// access to full map
    const Teuchos::RCP<const Epetra_Map>& FullMap() const;

    /// global extractor
    const Teuchos::RCP<LINALG::MultiMapExtractor>& GlobalExtractor() const;

    //! set-up of global system of equations of coupled problem
    void SetupSystem(
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector>                 rhs,
        Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art,
        Teuchos::RCP<const Epetra_Vector>           rhs_cont,
        Teuchos::RCP<const Epetra_Vector>           rhs_art
        );

    /// setup global vector
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    void SetupVector(
        Teuchos::RCP<Epetra_Vector>                 vec,
        Teuchos::RCP<const Epetra_Vector>           vec_cont,
        Teuchos::RCP<const Epetra_Vector>           vec_art
        );

    /// extract single field vectors
    /// (o)    vec:      combined vector containing both artery and continuous field quantities
    /// (i)    vec_cont: vector containing quantities from continuous field
    /// (i)    vec_art:  vector containing quantities from artery field
    void ExtractSingleFieldVectors(
      Teuchos::RCP<const Epetra_Vector>       globalvec,
      Teuchos::RCP<const Epetra_Vector>&      vec_cont,
      Teuchos::RCP<const Epetra_Vector>&      vec_art
      );

    // check if initial fields on coupled DOFs are equal
    void CheckInitialFields(
        Teuchos::RCP<const Epetra_Vector>      vec_cont,
        Teuchos::RCP<const Epetra_Vector>      vec_art
        );

    //! access artery (1D) dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const;

    //! access full dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() const;


  private:

    //! set-up of global rhs vector of coupled problem
    void SetupRHS(
        Teuchos::RCP<Epetra_Vector>                 rhs,
        Teuchos::RCP<const Epetra_Vector>           rhs_cont,
        Teuchos::RCP<const Epetra_Vector>           rhs_art
        );

    //! set-up of global matrix of coupled problem
    void SetupMatrix(
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art
        );

    //! setup map extractor
    //  full map -> all DOFs
    //  maps(0)  -> coupled DOFs
    //  maps(1)  -> uncoupled DOFs
    void SetupMapExtractor(
        Teuchos::RCP<LINALG::MultiMapExtractor>    mapextractor,       //!< the map extractor to setup
        Teuchos::RCP<DRT::Discretization>          dis,                //!< discretization
        const std::vector<int>&                    coupleddofs         //!< vector with DOFs to couple
        );

    //! check if dirichlet BC defined on coupled dofs, which is not possible
    void CheckDbcOnCoupledDofs(
        Teuchos::RCP<DRT::Discretization>          dis,                //!< discretization
        const Teuchos::RCP<const Epetra_Map>&      coupleddofmap       //!< map with coupled DOFs
        );

    //! artery (1D) discretization
    Teuchos::RCP<DRT::Discretization>  arterydis_;
    //! continous field (2D, 3D) discretization
    Teuchos::RCP<DRT::Discretization>  contdis_;

    //! name of the condition
    const std::string condname_;

    //! coupled dofs of the single fields
    std::vector<int> coupleddofs_art_;
    std::vector<int> coupleddofs_cont_;
    int num_coupled_dofs_;

    //! dof row map (not splitted)
    Teuchos::RCP<Epetra_Map> fullmap_;

    //! dof row map splitted in (field) blocks
    Teuchos::RCP<LINALG::MultiMapExtractor> blockrowdofmap_;

    //! extractors for continous field and artery field, maps(0) -> Coupled Dofs, maps(1) uncoupled Dofs
    Teuchos::RCP<LINALG::MultiMapExtractor> contfieldex_;
    Teuchos::RCP<LINALG::MultiMapExtractor> artex_;

    //! global extractor
    Teuchos::RCP<LINALG::MultiMapExtractor>   globalex_;

    //! coupling adapter
    Teuchos::RCP<ADAPTER::Coupling> artcontfieldcoup_;

    /// needed for matrix transforms
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> sbbtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> sbitransform_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> sibtransform_;

  };

}  // namespace POROMULTIPHASESCATRA


#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_NODEBASED_H_ */
