/*----------------------------------------------------------------------*/
/*! \file

\brief Special version of block matrix that includes the AMGnxn preconditioner

\maintainer Matthias Mayr

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_OVERLAPPREC_AMGNXN_H
#define FSI_OVERLAPPREC_AMGNXN_H

#ifdef HAVE_MueLu

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"

#include "fsi_overlapprec.H"

namespace FSI
{
  class OverlappingBlockMatrixAMGnxn : public OverlappingBlockMatrix
  {
   public:
    // All these arguments are required to construct the base class OverlappingBlockMatrix
    OverlappingBlockMatrixAMGnxn(const LINALG::MultiMapExtractor& maps,
        ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
        bool structuresplit, std::string xml_file, FILE* err = NULL,
        std::string ProblemType = "FSI");

    // This function is used to apply the preconditioner
    virtual void SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      P_->ApplyInverse(X, Y);
      return;
    }

    // This function is used to create the preconditioner
    virtual void SetupPreconditioner();

   private:
    // This is the preconditioner
    Teuchos::RCP<Epetra_Operator> P_;        // The underlying preconditioner object
    Teuchos::RCP<BlockSparseMatrixBase> A_;  // A own copy of the system matrix
    std::string xml_file_;                   // XMLfile containing the amg nxn prec
    std::string ProblemType_;
  };

  class AMGnxnInterfaceFSI
  {
   public:
    AMGnxnInterfaceFSI(Teuchos::RCP<LINALG::Preconditioner> structuresolver,
        Teuchos::RCP<LINALG::Preconditioner> fluidsolver,
        Teuchos::RCP<LINALG::Preconditioner> alesolver, std::string amgnxn_xml);

    std::vector<int> GetNumPdes() { return num_pdes_; }
    std::vector<int> GetNullSpacesDim() { return null_spaces_dim_; }
    std::vector<Teuchos::RCP<std::vector<double>>> GetNullSpacesData() { return null_spaces_data_; }
    Teuchos::ParameterList GetPreconditionerParams() { return prec_params_; }
    Teuchos::ParameterList GetSmoothersParams() { return smoo_params_; }
    std::string GetPreconditionerType() { return prec_type_; }
    virtual void ParseXML();
    void Check();

   protected:
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList prec_params_;
    Teuchos::ParameterList smoo_params_;
    std::string prec_type_;
    Teuchos::RCP<LINALG::Preconditioner> structuresolver_;
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver_;
    Teuchos::RCP<LINALG::Preconditioner> alesolver_;
    std::string amgnxn_xml_;

   private:
    void Setup();
    void Default_AMG_BGS(Teuchos::ParameterList& params);
  };

  class AMGnxnInterfaceLungFSI : public AMGnxnInterfaceFSI
  {
   public:
    AMGnxnInterfaceLungFSI(Teuchos::RCP<LINALG::Preconditioner> structuresolver,
        Teuchos::RCP<LINALG::Preconditioner> fluidsolver,
        Teuchos::RCP<LINALG::Preconditioner> alesolver, std::string amgnxn_xml);
    virtual void ParseXML();

   private:
    void Default_SCHUR_ILU_KLU(Teuchos::ParameterList& params);
  };



}  // namespace FSI

#endif  // HAVE_MueLu
#endif  // FSI_OVERLAPPREC_AMGNXN_H
