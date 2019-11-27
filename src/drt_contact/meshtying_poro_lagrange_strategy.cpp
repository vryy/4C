/*-----------------------------------------------------------------------*/
/*! \file
\brief Meshtying of porous media using Lagrange multipliers

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Matthias Mayr
// Originates from contact_poro_lagrange_strategy

\maintainer Matthias Mayr

\level 3
*/
/*-----------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "meshtying_poro_lagrange_strategy.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                      h.Willmann    2015|
 *----------------------------------------------------------------------*/
CONTACT::PoroMtLagrangeStrategy::PoroMtLagrangeStrategy(const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface>> interface, int dim,
    Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
    : MtLagrangeStrategy(DofRowMap, NodeRowMap, params, interface, dim, comm, alphaf, maxdof)
{
  return;
}


/*----------------------------------------------------------------------*
 | Poro Meshtying initialization calculations         h.Willmann    2015|
 *----------------------------------------------------------------------*/
void CONTACT::PoroMtLagrangeStrategy::InitializePoroMt(
    Teuchos::RCP<LINALG::SparseMatrix>& kteffoffdiag)

{
  Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteffoffdiag);

  fvelrow_ = Teuchos::rcp(new Epetra_Map(kteffmatrix->OperatorDomainMap()));

  return;
}


/*----------------------------------------------------------------------*
 | Poro Meshtying method regarding coupling terms      h.Willmann   2015|
 *----------------------------------------------------------------------*/
void CONTACT::PoroMtLagrangeStrategy::EvaluateMeshtyingPoroOffDiag(
    Teuchos::RCP<LINALG::SparseMatrix>& kteffoffdiag)
{
  // system type
  INPAR::CONTACT::SystemType systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(), "SYSTEM");

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

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
      dserror("Condensation only for dual LM");

    // h.Willmann actual method

    // complete stiffness matrix
    // (this is a prerequisite for the Split2x2 methods to be called later)
    kteffoffdiag->Complete();

    /**********************************************************************/
    /* Split kteffoffdiag into 3 block matrix rows                        */
    /**********************************************************************/
    // we want to split k into 3 rows n, m and s
    Teuchos::RCP<LINALG::SparseMatrix> cn, cm, cs;

    // temporarily we need the block row csm
    // (FIXME: because a direct SplitMatrix3x1 is missing here!)
    Teuchos::RCP<LINALG::SparseMatrix> csm;


    // some temporary Teuchos::RCPs
    Teuchos::RCP<Epetra_Map> tempmap1;
    Teuchos::RCP<Epetra_Map> tempmap2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx4;

    Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix =
        Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteffoffdiag);

    //    std::cout<< " kteffmatrix " << std::endl;
    //    kteffmatrix->DomainMap().Print(std::cout);

    if (ParRedist())  // asdf
    {
      dserror("no parallel redistribution of poro meshtying implemented - feel free to implement");
    }

    // first split into slave/master block row + remaining part
    LINALG::SplitMatrix2x2(
        kteffmatrix, gsmdofrowmap_, gndofrowmap_, fvelrow_, tempmap1, csm, tempmtx1, cn, tempmtx2);

    //    std::cout<< " tempmap1 " << std::endl;
    //    tempmap1->Print(std::cout);

    // second split slave/master block row
    LINALG::SplitMatrix2x2(
        csm, gsdofrowmap_, gmdofrowmap_, fvelrow_, tempmap2, cs, tempmtx3, cm, tempmtx4);

    // store some stuff for the recovery of the lagrange multiplier
    cs_ = cs;


    /**********************************************************************/
    /* Build the final matrix block row                                   */
    /**********************************************************************/
    // cn: nothing to do

    // cm: add T(mbar)*cs
    Teuchos::RCP<LINALG::SparseMatrix> cmmod =
        Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_, 100));
    cmmod->Add(*cm, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> cmadd =
        LINALG::MLMultiply(*GetMHat(), true, *cs, false, false, false, true);
    cmmod->Add(*cmadd, false, 1.0, 1.0);
    cmmod->Complete(cm->DomainMap(), cm->RowMap());

    // cs: nothing to do as it remains zero

    /**********************************************************************/
    /* Global setup of kteffoffdiagnew,  (including meshtying)            */
    /**********************************************************************/
    Teuchos::RCP<LINALG::SparseMatrix> kteffoffdiagnew = Teuchos::rcp(
        new LINALG::SparseMatrix(*ProblemDofs(), 81, true, false, kteffmatrix->GetMatrixtype()));

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
    dserror("Trying to use not condensed PoroMeshtying --- Feel Free to implement!");
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
