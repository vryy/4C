/*-----------------------------------------------------------------------*/
/*! \file
\brief Meshtying of porous media using Lagrange multipliers

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Matthias Mayr
// Originates from contact_poro_lagrange_strategy


\level 3
*/
/*-----------------------------------------------------------------------*/

#include "4C_contact_meshtying_poro_lagrange_strategy.hpp"

#include "4C_inpar_contact.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_SerialComm.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                      h.Willmann    2015|
 *----------------------------------------------------------------------*/
CONTACT::PoroMtLagrangeStrategy::PoroMtLagrangeStrategy(const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<MORTAR::Interface>> interface, int dim, Teuchos::RCP<Epetra_Comm> comm,
    double alphaf, int maxdof)
    : MtLagrangeStrategy(DofRowMap, NodeRowMap, params, interface, dim, comm, alphaf, maxdof)
{
  return;
}


/*----------------------------------------------------------------------*
 | Poro Meshtying initialization calculations         h.Willmann    2015|
 *----------------------------------------------------------------------*/
void CONTACT::PoroMtLagrangeStrategy::InitializePoroMt(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteffoffdiag)

{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kteffmatrix =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(kteffoffdiag);

  fvelrow_ = Teuchos::rcp(new Epetra_Map(kteffmatrix->OperatorDomainMap()));

  return;
}


/*----------------------------------------------------------------------*
 | Poro Meshtying method regarding coupling terms      h.Willmann   2015|
 *----------------------------------------------------------------------*/
void CONTACT::PoroMtLagrangeStrategy::EvaluateMeshtyingPoroOffDiag(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteffoffdiag)
{
  // system type
  INPAR::CONTACT::SystemType systype =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::SystemType>(Params(), "SYSTEM");

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed ||
      systype == INPAR::CONTACT::system_condensed_lagmult)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      FOUR_C_THROW("Condensation only for dual LM");

    // h.Willmann actual method

    // complete stiffness matrix
    // (this is a prerequisite for the Split2x2 methods to be called later)
    kteffoffdiag->Complete();

    /**********************************************************************/
    /* Split kteffoffdiag into 3 block matrix rows                        */
    /**********************************************************************/
    // we want to split k into 3 rows n, m and s
    Teuchos::RCP<CORE::LINALG::SparseMatrix> cn, cm, cs;

    // temporarily we need the block row csm
    // (FIXME: because a direct SplitMatrix3x1 is missing here!)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> csm;


    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap1;
    Teuchos::RCP<Epetra_Map> tempmap2;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> tempmtx3;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> tempmtx4;

    Teuchos::RCP<CORE::LINALG::SparseMatrix> kteffmatrix =
        Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(kteffoffdiag);

    //    std::cout<< " kteffmatrix " << std::endl;
    //    kteffmatrix->DomainMap().Print(std::cout);

    if (ParRedist())  // asdf
    {
      FOUR_C_THROW(
          "no parallel redistribution of poro meshtying implemented - feel free to implement");
    }

    // first split into slave/master block row + remaining part
    CORE::LINALG::SplitMatrix2x2(
        kteffmatrix, gsmdofrowmap_, gndofrowmap_, fvelrow_, tempmap1, csm, tempmtx1, cn, tempmtx2);

    //    std::cout<< " tempmap1 " << std::endl;
    //    tempmap1->Print(std::cout);

    // second split slave/master block row
    CORE::LINALG::SplitMatrix2x2(
        csm, gsdofrowmap_, gmdofrowmap_, fvelrow_, tempmap2, cs, tempmtx3, cm, tempmtx4);

    // store some stuff for the recovery of the lagrange multiplier
    cs_ = cs;


    /**********************************************************************/
    /* Build the final matrix block row                                   */
    /**********************************************************************/
    // cn: nothing to do

    // cm: add T(mbar)*cs
    Teuchos::RCP<CORE::LINALG::SparseMatrix> cmmod =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(*gmdofrowmap_, 100));
    cmmod->Add(*cm, false, 1.0, 1.0);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> cmadd =
        CORE::LINALG::MLMultiply(*GetMHat(), true, *cs, false, false, false, true);
    cmmod->Add(*cmadd, false, 1.0, 1.0);
    cmmod->Complete(cm->DomainMap(), cm->RowMap());

    // cs: nothing to do as it remains zero

    /**********************************************************************/
    /* Global setup of kteffoffdiagnew,  (including meshtying)            */
    /**********************************************************************/
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kteffoffdiagnew =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(
            *ProblemDofs(), 81, true, false, kteffmatrix->GetMatrixtype()));

    // add n matrix row
    kteffoffdiagnew->Add(*cn, false, 1.0, 1.0);

    // add m matrix row
    kteffoffdiagnew->Add(*cmmod, false, 1.0, 1.0);

    // s matrix row remains zero (thats what it was all about)

    kteffoffdiagnew->Complete(kteffmatrix->DomainMap(), kteffmatrix->RangeMap());

    kteffoffdiag = kteffoffdiagnew;
  }
  else
  {
    FOUR_C_THROW("Trying to use not condensed PoroMeshtying --- Feel Free to implement!");
  }
  return;
}


/*----------------------------------------------------------------------*
 | Poro Recovery method for structural displacement LM  h.Willmann  2015|
 *----------------------------------------------------------------------*/
void CONTACT::PoroMtLagrangeStrategy::RecoverCouplingMatrixPartofLMP(
    Teuchos::RCP<Epetra_Vector> veli)
{
  Teuchos::RCP<Epetra_Vector> zfluid = Teuchos::rcp(new Epetra_Vector(z_->Map(), true));

  Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  cs_->Multiply(false, *veli, *mod);
  zfluid->Update(-1.0, *mod, 1.0);
  Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*zfluid));
  GetDInverse()->Multiply(true, *zcopy, *zfluid);
  zfluid->Scale(1 / (1 - alphaf_));

  z_->Update(1.0, *zfluid, 1.0);  // Add FluidCoupling Contribution to LM!
  return;
}

FOUR_C_NAMESPACE_CLOSE
