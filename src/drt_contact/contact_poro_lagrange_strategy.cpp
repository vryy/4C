/*---------------------------------------------------------------------*/
/*! \file
\brief Contact Strategy handling the porous no penetraction condition on the active contact
interface

\level 3


*/
/*---------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "contact_poro_lagrange_strategy.H"
#include "contact_interface.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../drt_lib/drt_utils.H"

// this is not nice, but since there is no chance at them moment to set the
// full displacement during a restart from outside before the contact tries to
// evaluate I'm force to do this atm...
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              ager 08/14|
 *----------------------------------------------------------------------*/
CONTACT::PoroLagrangeStrategy::PoroLagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr, const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
    Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof, bool poroslave, bool poromaster)
    : MonoCoupledLagrangeStrategy(
          data_ptr, DofRowMap, NodeRowMap, params, interface, dim, comm, alphaf, maxdof),
      no_penetration_(params.get<bool>("CONTACTNOPEN")),
      nopenalpha_(0.0),
      poroslave_(poroslave),
      poromaster_(poromaster)
{
  if (!poroslave_ and !poromaster_)
    dserror(
        "you called a poroelastic meshtying method without participating poroelastic domains on "
        "your interface");
  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      ager 12/16|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::DoReadRestart(IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  Teuchos::RCP<DRT::Discretization> discret = DRT::Problem::Instance()->GetDis("structure");
  if (discret == Teuchos::null) dserror("didn't get my discretization");

  Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*discret->DofColMap(), true));
  // it's clear that we get some zeros here ... but poroelast monolithic fixes this a little bit
  // later by doing the same thing with correct displacements again :-)
  LINALG::Export(*dis, *global);
  SetParentState("displacement", global, discret);

  // Call (nearly absolute)Base Class
  CONTACT::CoAbstractStrategy::DoReadRestart(reader, dis, cparams_ptr);
  return;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                ager 12/16 |
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::Setup(bool redistributed, bool init)
{
  // Call Base Class
  CONTACT::CoAbstractStrategy::Setup(redistributed, init);

  if (no_penetration_) SetupNoPenetrationCondition();
}

/*----------------------------------------------------------------------*
 | Activate No-Penetration for the contact surface (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SetupNoPenetrationCondition()
{
  lambda_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  lambdaold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));
  if (!poroslave_ and poromaster_)
    dserror("poroelastic meshtying/contact method needs the slave side to be poroelastic");
  /*
   *  The first error may also occur when there is no single element on the interface (it is empty)
   * but POROELASTICITY DYNAMIC coupling algorithmus (COUPALGO) is chosen as
   * poro_monolithicmeshtying that means that the method creates an interface but has nothing to
   * fill in
   *
   *  the second error shows that there are no elements with PoroCoupling on the slave side of the
   * interface and the master side is fully poroelastic
   *
   *  having the condition with structure slave and poro master the problem is that there is no
   * diagonal D Matrix, from the porofluid lagrange multiplier equality condition, that can be
   * inverted easily to condense out the respective lagrange multipliers
   *
   *  there are two alternatives:
   *  1. do not condense out the lagrange multiplier for the porofluid meshtying condition in
   *  this constellation (master/slave) and solve the saddlepoint system
   *
   *  2. Invert the M Matrix to condense out the Lagrange multiplier
   *  - which is only applicable with adequate costs for small problems as the M Matrix is not at
   * all diagonal, except for matching meshes - for matching meshes it doesn't hurt to chose the
   * poro side as the slave side. Being applicable only for small problems this alternative is too
   * restricted in its application overall to be implemented.
   */
}

/*----------------------------------------------------------------------*
 | initialize global poro contact variables                             |
 |                            for next Newton step (public)   ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroInitialize(
    ADAPTER::Coupling& coupfs, Teuchos::RCP<const Epetra_Map> fluiddofs, bool fullinit)
{
  if (fullinit)  // fullinit is true by default, but needed when this method is called for
                 // meshtying, as the maps and matrix mapping stay the same for meshtying. would
                 // work without this but would do things repeatedly unnecessarily
  {
    if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
    {
      //  (1)                                                          //
      //      Get required fluid maps from structural maps             //
      //                                                               //
      fgsdofrowmap_ = coupfs.MasterToSlaveMap(gsdofrowmap_);
      fgmdofrowmap_ = coupfs.MasterToSlaveMap(gmdofrowmap_);
      fgsmdofrowmap_ = coupfs.MasterToSlaveMap(gsmdofrowmap_);
      fgndofrowmap_ = LINALG::SplitMap(*fluiddofs,
          *fgsmdofrowmap_);  // Not equal to transforming gndofrowmap_ (pressure dofs missing!)
      fgactivedofs_ = coupfs.MasterToSlaveMap(gactivedofs_);
      falldofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(*fluiddofs));
      fgactiven_ = coupfs.MasterToSlaveMap(gactiven_);
      fgactivet_ = coupfs.MasterToSlaveMap(gactivet_);
    }
  }
  //  (2)                                                          //
  //      Initialize Matrices                                      //
  //                                                               //
  if (fullinit)
  {
    if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
    {
      // (re)setup global nCoup Vector
      NCoup_ = Teuchos::rcp(new Epetra_Vector(*gactiven_));

      // (re)setup global linearisation matrices of nCoup
      NCoup_lindisp_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_, 10));
      NCoup_linvel_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_, 10));

      // (re)setup global tangential and lin(tangential)*lambda matrices
      Tangential_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 10));
      linTangentiallambda_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 10));

      // (re)setup global lin of D & M * lambda - Matrix
      porolindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *gsdofrowmap_, 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
      porolinmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *gmdofrowmap_, 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
    }
  }
  else
  {
    // In the case of meshtying resetting the matrices is sufficient as they retain their size
    NCoup_->PutScalar(0.0);

    // (re)setup global linearisation matrices of nCoup
    NCoup_lindisp_->Zero();
    NCoup_linvel_->Zero();

    // (re)setup global tangential and lin(tangential)*lambda matrices
    Tangential_->Zero();
    linTangentiallambda_->Zero();

    // (re)setup global lin of D & M * lambda - Matrix
    porolindmatrix_->Zero();
    porolinmmatrix_->Zero();
  }
  //  (3)                                                          //
  //      Assemble Matrices                                        //
  //                                                               //
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
    {
      interface_[i]->AssembleNCoup(*NCoup_);

      interface_[i]->AssembleNCoupLin(*NCoup_lindisp_, coupfs);
      interface_[i]->AssembleNCoupLin(*NCoup_linvel_, coupfs, true);

      interface_[i]->AssembleTN(Tangential_, Teuchos::null);
      interface_[i]->AssembleTNderiv(linTangentiallambda_, Teuchos::null,
          true);  // use lambda(n +1) for tangential condition!!!

      interface_[i]->AssembleLinDM(*porolindmatrix_, *porolinmmatrix_, true);
    }
  }

  //  (4)                                                          //
  //      Complete Matrices                                        //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  {
    NCoup_lindisp_->Complete(*ProblemDofs(), *gactiven_);
    NCoup_linvel_->Complete(*fluiddofs, *gactiven_);

    Tangential_->Complete(*gactivedofs_, *gactivet_);
    linTangentiallambda_->Complete(*gsmdofrowmap_, *gactivet_);

    porolindmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);
    porolinmmatrix_->Complete(*gsmdofrowmap_, *gmdofrowmap_);
  }

  //  (5)                                                          //
  //      Reset Matrix Transform Objects                           //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  {
    linncoupveltransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    linncoupdisptransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    tanginvtransform_ = Teuchos::rcp(new LINALG::MatrixRowColTransform);
    lintangentlambdatransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    porolindmatrixtransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    porolinmmatrixtransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    mhataamtransform_ = Teuchos::rcp(new LINALG::MatrixRowColTransform);  // h.Willmann
    dhattransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    doldtransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    moldtransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
    invDatransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
  }

  //  (6)                                                          //
  //      Transform Matrices from structural dofs to fluid dofs    //
  //                                                               //
  if (no_penetration_ && (IsInContact() || WasInContact() || WasInContactLastTimeStep()))
  // eventually a coupling object just on the mortar interface would make sense!!! ChrAg
  {
    // transform matrices coming from contact to fluid maps, as they are all in structure maps!
    //
    // A generell problem here is that we would need to update coupling objects in everey newton
    // step if the active set changes. To avoid this, a 'bigger' coupling object is used - but
    // therefore now the Row-Maps of the created Sparse Matrixes are to big! --- Leads to problems
    // for Matrix - Multiplications where the Row - Map is used!
    //
    // At the moment this is solved by using SplitMatrix2x2 to cut out just the relevant part of the
    // matrix, but this shouldn't be the final solution!
    //
    //************************************************************************************************
    //
    Teuchos::RCP<Epetra_Vector> tmpfullncoup =
        Teuchos::rcp(new Epetra_Vector(*coupfs.MasterDofMap()));
    LINALG::Export(*NCoup_, *tmpfullncoup);
    tmpfullncoup = coupfs.MasterToSlave(tmpfullncoup);
    fNCoup_ = Teuchos::rcp(new Epetra_Vector(*fgactiven_));
    LINALG::Export(*tmpfullncoup, *fNCoup_);
    //
    //************************************************************************************************
    //
    fdoldtransp_ = Teuchos::rcp<LINALG::SparseMatrix>(
        new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*doldtransform_)(
        *dold_->Transpose(), 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fdoldtransp_, false);
    fdoldtransp_->Complete(dold_->RowMap(), *fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    if (poromaster_)
    {
      fmoldtransp_ = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*moldtransform_)(
          *mold_->Transpose(), 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fmoldtransp_, false);
      fmoldtransp_->Complete(mold_->RowMap(), *fgmdofrowmap_);
    }
    //
    //************************************************************************************************
    //
    fporolindmatrix_ = Teuchos::rcp<LINALG::SparseMatrix>(
        new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*porolindmatrixtransform_)(
        *porolindmatrix_, 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fporolindmatrix_, false);
    fporolindmatrix_->Complete(porolindmatrix_->DomainMap(), *fgsdofrowmap_);
    //
    //************************************************************************************************
    //
    if (poromaster_)
    {
      fporolinmmatrix_ = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*porolinmmatrixtransform_)(*porolinmmatrix_, 1.0, ADAPTER::CouplingMasterConverter(coupfs),
          *fporolinmmatrix_, false);
      fporolinmmatrix_->Complete(porolinmmatrix_->DomainMap(), *fgmdofrowmap_);
    }
    // porolinmmatrixtransform_ is no longer missing as twosided poro meshtying is considered!
    //
    //************************************************************************************************
    if (poromaster_)
    //
    {
      Teuchos::RCP<LINALG::SparseMatrix> tmpfmhataam = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));

      fmhataam_ = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*fgmdofrowmap_, 1, true, false));
      (*mhataamtransform_)(*mhataam_, 1.0, ADAPTER::CouplingMasterConverter(coupfs),
          ADAPTER::CouplingMasterConverter(coupfs), *tmpfmhataam, false, false);

      tmpfmhataam->Complete(*fgmdofrowmap_, *falldofrowmap_);

      // better solution to get maps as wanted? -- for this matrix map as important as there will be
      // a matrix-matrix multiplication

      Teuchos::RCP<Epetra_Map> restfgmdofrowmap, restfgactivedofs;
      Teuchos::RCP<LINALG::SparseMatrix> tmpm1, tmpm2, tmpm3;

      // This should just be a temporary solution to change the row map of the matrix ...
      LINALG::SplitMatrix2x2(tmpfmhataam, fgactivedofs_, restfgactivedofs, fgmdofrowmap_,
          restfgmdofrowmap, fmhataam_, tmpm1, tmpm2, tmpm3);
      fmhataam_->Complete(*fgmdofrowmap_, *fgactivedofs_);
    }
    //
    //************************************************************************************************
    //
    fdhat_ = Teuchos::rcp<LINALG::SparseMatrix>(
        new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
    (*dhattransform_)(*dhat_, 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fdhat_, false);
    fdhat_->Complete(dhat_->DomainMap(), *fgactivedofs_);
    // fdhat is expected to be zero in this method but still used for condensation
    //
    //************************************************************************************************
    //
    if (gactivedofs_->NumGlobalElements())
    {
      Teuchos::RCP<LINALG::SparseMatrix> tanginvD =
          LINALG::MLMultiply(*Tangential_, false, *invda_, true, false, false, true);
      Teuchos::RCP<LINALG::SparseMatrix> tmpftanginvD = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*tanginvtransform_)(*tanginvD, 1.0, ADAPTER::CouplingMasterConverter(coupfs),
          ADAPTER::CouplingMasterConverter(coupfs), *tmpftanginvD, false, false);
      tmpftanginvD->Complete(*fgactivedofs_, *falldofrowmap_);
      // better solution to get maps as wanted? -- for this matrix map as important as there will be
      // a matrix-matrix multiplication

      Teuchos::RCP<Epetra_Map> restfgactivet, restfgactivedofs;
      Teuchos::RCP<LINALG::SparseMatrix> tmpm1, tmpm2, tmpm3;

      // This should just be a temporary solution to change the row map of the matrix ...
      LINALG::SplitMatrix2x2(tmpftanginvD, fgactivet_, restfgactivet, fgactivedofs_,
          restfgactivedofs, ftanginvD_, tmpm1, tmpm2, tmpm3);

#if (0)
      // Some solutions that do not work!

      //    ftanginvD_ = Teuchos::rcp<LINALG::SparseMatrix>(new
      //    LINALG::SparseMatrix(*fgactivet_,1,true,false));
      //  ftanginvD_->Assign(View,*tmpftanginvD);

      //   ftanginvD_->Add(*tmpftanginvD,false,1.0,1.0);
      //   ftanginvD_->Complete(*fgactivedofs_,*fgactivet_);

      //  int err = ftanginvD_->EpetraMatrix()->ReplaceRowMap(*fgactivet_);
      // if (err) dserror("RRM STOPPED .... with err! %d",err);
#endif
      //
      //************************************************************************************************
      //
      fNCoup_linvel_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_, 108, false));

      (*linncoupveltransform_)(
          *NCoup_linvel_, 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fNCoup_linvel_, false);
      fNCoup_linvel_->Complete(*falldofrowmap_, *fgactiven_);
      //
      //************************************************************************************************
      //
      fNCoup_lindisp_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_, 81, false));
      (*linncoupdisptransform_)(
          *NCoup_lindisp_, 1.0, ADAPTER::CouplingMasterConverter(coupfs), *fNCoup_lindisp_, false);
      fNCoup_lindisp_->Complete(*ProblemDofs(), *fgactiven_);
      //
      //************************************************************************************************
      //
      flinTangentiallambda_ = Teuchos::rcp(new LINALG::SparseMatrix(*falldofrowmap_, 81, false));
      (*lintangentlambdatransform_)(*linTangentiallambda_, 1.0,
          ADAPTER::CouplingMasterConverter(coupfs), *flinTangentiallambda_, false);
      flinTangentiallambda_->Complete(*gsdofrowmap_, *fgactivet_);

#if (0)  // just in case
         // only for parallel redistribution case
         //     if (parredist) //care about that if its in contact ... ChrAg Todo
         //     {
         //       NCoup_lindisp_    = MORTAR::MatrixRowColTransform(NCoup_lindisp_,
         //       interface_[i]->NormalDofs(), masterdofrowmap_); NCoup_linvel_     =
         //       MORTAR::MatrixRowColTransform(NCoup_linvel_, interface_[i]->NormalDofs(),
         //       slavedofrowmap_);
         //     }

      // only for parallel redistribution case
      //     if (parredist) //care about that if its in contact ... ChrAg Todo
      //     {
      //       Tangential_    = MORTAR::MatrixRowColTransform(Tangential_,
      //       interface_[i]->TangentialDofs(), slavedofrowmap_); linTangentiallambda_     =
      //       MORTAR::MatrixRowColTransform(linTangentiallambda_, interface_[i]->TangentialDofs(),
      //       slavedofrowmap_);
      //     }


      //  Assemble lin of D & M - Matrix --
      //  Teuchos::RCP<LINALG::SparseMatrix> lindglobal = Teuchos::rcp(new
      //  LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      //  //ChrAg Teuchos::RCP<LINALG::SparseMatrix> linmglobal = Teuchos::rcp(new
      //  LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      //  //ChrAg



      // porolindmatrix_ = lindglobal;
      // porolinmmatrix_ = linmglobal;
#endif
      //
      //************************************************************************************************
      //
      finvda_ = Teuchos::rcp<LINALG::SparseMatrix>(
          new LINALG::SparseMatrix(*falldofrowmap_, 1, true, false));
      (*invDatransform_)(*invda_, 1.0, ADAPTER::CouplingMasterConverter(coupfs), *finvda_, false);
      finvda_->Complete(invda_->DomainMap(), *fgactivedofs_);
    }
  }
  return;
}  // CONTACT::PoroLagrangeStrategy::PoroInitialize

/*-------------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration             |
 |  condition on contact surface (pure porous problem)(public)   ager 07/15|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluatePoroNoPenContact(
    Teuchos::RCP<LINALG::SparseMatrix>& k_fseff, Teuchos::RCP<LINALG::SparseMatrix>& Feff,
    Teuchos::RCP<Epetra_Vector>& feff)
{
  EvaluateMatPoroNoPen(k_fseff, feff);

  EvaluateOtherMatPoroNoPen(Feff, 0);
}

/*-------------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration             |
 |  condition on contact surface(public)                         ager 07/15|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluatePoroNoPenContact(
    Teuchos::RCP<LINALG::SparseMatrix>& k_fseff,
    std::map<int, Teuchos::RCP<LINALG::SparseMatrix>*>& Feff, Teuchos::RCP<Epetra_Vector>& feff)
{
  EvaluateMatPoroNoPen(k_fseff, feff);

  // Take care of the alternative condensation of the off-diagonal blocks!!!
  std::map<int, Teuchos::RCP<LINALG::SparseMatrix>*>::iterator matiter;
  for (matiter = Feff.begin(); matiter != Feff.end(); ++matiter)
  {
    EvaluateOtherMatPoroNoPen(*(matiter->second), matiter->first);
  }
}

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration          |
 |  condition on contact surface (public)                     ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluateMatPoroNoPen(
    Teuchos::RCP<LINALG::SparseMatrix>& k_fseff, Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ || (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()))
    return;
  // h.Willmann this method should be renamed as it handles twosided poro meshtying now aswell and
  // is able to handle fluid coupling for twosided contact

  nopenalpha_ = alphaf_;  // to use different alpha for nopen condition (not used at the moment)

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************

  // double-check if this is a dual LM system
  if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    dserror("Condensation only for dual LM");

  /**********************************************************************/
  /* (2) Add contact stiffness terms to kteff                           */
  /**********************************************************************/

  // transform if necessary
  if (ParRedist())
  {
    dserror("CHECK ME!!!");
    lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_, pgsdofrowmap_);
    linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_, pgmdofrowmap_);
  }

  k_fseff->UnComplete();
  k_fseff->Add(*fporolindmatrix_, false, (1.0 - nopenalpha_) * 1.0, 1.0);

  if (poromaster_)
    k_fseff->Add(*fporolinmmatrix_, false, (1.0 - nopenalpha_) * 1.0,
        1.0);  // is needed only for twosided poro contact or meshtying

  k_fseff->Complete(*ProblemDofs(), *falldofrowmap_);  // gets bigger because of linearisation
                                                       // w.r.t. to pure structural displacements!
  /**********************************************************************/
  /* (3) Split k_fseff and Feff into 3x3 matrix blocks                             */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_ss, k_fs_sm, k_fs_sn, k_fs_ms, k_fs_mm, k_fs_mn, k_fs_ns,
      k_fs_nm, k_fs_nn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_smsm, k_fs_smn, k_fs_nsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<Epetra_Map> ftempmap1, ftempmap2, ftempmap3;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;

  // split into slave/master part + structure part
  Teuchos::RCP<LINALG::SparseMatrix> k_fseffmatrix =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fseff);
  if (ParRedist())
  {
    dserror("CHECK ME!");
    // split and transform to redistributed maps
    //        LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    //        ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    //        ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    //        knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitMatrix2x2(k_fseffmatrix, fgsmdofrowmap_, fgndofrowmap_, gsmdofrowmap_,
        gndofrowmap_, k_fs_smsm, k_fs_smn, k_fs_nsm, k_fs_nn);
  }

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(k_fs_smsm, fgsdofrowmap_, fgmdofrowmap_, gsdofrowmap_, gmdofrowmap_,
      k_fs_ss, k_fs_sm, k_fs_ms, k_fs_mm);
  LINALG::SplitMatrix2x2(k_fs_smn, fgsdofrowmap_, fgmdofrowmap_, gndofrowmap_, tempmap, k_fs_sn,
      tempmtx1, k_fs_mn, tempmtx2);
  LINALG::SplitMatrix2x2(k_fs_nsm, fgndofrowmap_, ftempmap1, gsdofrowmap_, gmdofrowmap_, k_fs_ns,
      k_fs_nm, tempmtx1, tempmtx2);

  /**********************************************************************/
  /* (4) Split feff into 3 subvectors                                   */
  /**********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (ParRedist())
  {
    dserror("CHECK ME!");
    // split and transform to redistributed maps
    LINALG::SplitVector(*ProblemDofs(), *feff, pgsmdofrowmap_, fsm, gndofrowmap_, fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    LINALG::Export(*fsm, *fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitVector(*falldofrowmap_, *feff, fgsmdofrowmap_, fsm, fgndofrowmap_, fn);
  }

  // abbreviations for slave  and master set
  int sset = fgsdofrowmap_->NumGlobalElements();
  int mset = fgmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));

  // do the vector splitting sm -> s+m
  LINALG::SplitVector(*fgsmdofrowmap_, *fsm, fgsdofrowmap_, fs, fgmdofrowmap_, fm);

  // store some stuff for static condensation of poro no pen. LM

  ffs_ = fs;

  cfssn_ = k_fs_sn;
  cfssm_ = k_fs_sm;
  cfsss_ = k_fs_ss;

  /**********************************************************************/
  /* (5) Split slave quantities into active / inactive                  */
  /**********************************************************************/

  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_aa, k_fs_ai, k_fs_ia, k_fs_ii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_an, k_fs_in, k_fs_am, k_fs_im, k_fs_ma, k_fs_mi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;
  Teuchos::RCP<Epetra_Map> fgidofs;

  // some more temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap1, tempmap2;
  Teuchos::RCP<Epetra_Map> ftempmap4, ftempmap5, ftempmap6, ftempmap7;

  // do the splitting
  LINALG::SplitMatrix2x2(
      k_fs_ss, fgactivedofs_, fgidofs, gactivedofs_, gidofs, k_fs_aa, k_fs_ai, k_fs_ia, k_fs_ii);
  LINALG::SplitMatrix2x2(k_fs_sn, fgactivedofs_, fgidofs, gndofrowmap_, tempmap1, k_fs_an, tempmtx1,
      k_fs_in, tempmtx2);
  LINALG::SplitMatrix2x2(k_fs_sm, fgactivedofs_, fgidofs, gmdofrowmap_, tempmap2, k_fs_am, tempmtx1,
      k_fs_im, tempmtx2);
  LINALG::SplitMatrix2x2(k_fs_ms, fgmdofrowmap_, ftempmap4, gactivedofs_, gidofs, k_fs_ma, k_fs_mi,
      tempmtx1, tempmtx2);

  // abbreviations for active and inactive set
  int aset = fgactivedofs_->NumGlobalElements();
  int iset = fgidofs->NumGlobalElements();

  // we want to split fsmod into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*fgidofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*fgsdofrowmap_, *fs, fgactivedofs_, fa, fgidofs, fi);

  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /* where K stands for k_fs and F!!!                                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do
  //---------------------------------------------------------- SECOND LINE --- Will just exist when
  // starting with two sided poro contact!!!
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_mnmod;

  Teuchos::RCP<LINALG::SparseMatrix> k_fs_mmmod;

  Teuchos::RCP<LINALG::SparseMatrix> k_fs_mimod;

  Teuchos::RCP<LINALG::SparseMatrix> k_fs_mamod;
  if (mset)
  {
    // kmn: add T(mhataam)*kan
    k_fs_mnmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_, 100));
    k_fs_mnmod->Add(*k_fs_mn, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> k_fs_mnadd =
        LINALG::MLMultiply(*fmhataam_, true, *k_fs_an, false, false, false, true);
    k_fs_mnmod->Add(*k_fs_mnadd, false, 1.0, 1.0);
    k_fs_mnmod->Complete(k_fs_mn->DomainMap(), k_fs_mn->RowMap());

    // kmm: add T(mhataam)*kam
    k_fs_mmmod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_, 100));
    k_fs_mmmod->Add(*k_fs_mm, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> k_fs_mmadd =
        LINALG::MLMultiply(*fmhataam_, true, *k_fs_am, false, false, false, true);
    k_fs_mmmod->Add(*k_fs_mmadd, false, 1.0, 1.0);
    k_fs_mmmod->Complete(k_fs_mm->DomainMap(), k_fs_mm->RowMap());

    // kmi: add T(mhataam)*kai
    if (iset)
    {
      k_fs_mimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_, 100));
      k_fs_mimod->Add(*k_fs_mi, false, 1.0, 1.0);
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_miadd =
          LINALG::MLMultiply(*fmhataam_, true, *k_fs_ai, false, false, false, true);
      k_fs_mimod->Add(*k_fs_miadd, false, 1.0, 1.0);
      k_fs_mimod->Complete(k_fs_mi->DomainMap(), k_fs_mi->RowMap());
    }

    // kma: add T(mhataam)*kaa
    if (aset)
    {
      k_fs_mamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_, 100));
      k_fs_mamod->Add(*k_fs_ma, false, 1.0, 1.0);
      Teuchos::RCP<LINALG::SparseMatrix> k_fs_maadd =
          LINALG::MLMultiply(*fmhataam_, true, *k_fs_aa, false, false, false, true);
      k_fs_mamod->Add(*k_fs_maadd, false, 1.0, 1.0);
      k_fs_mamod->Complete(k_fs_ma->DomainMap(), k_fs_ma->RowMap());
    }
  }

  //----------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------
  // fdhat is expected to be zero here but still used for condensation
  // kin: subtract T(dhat)*kan --
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_inmod =
      Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs, 100));
  k_fs_inmod->Add(*k_fs_in, false, 1.0, 1.0);
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_inadd =
      LINALG::MLMultiply(*fdhat_, true, *k_fs_an, false, false, false, true);
  k_fs_inmod->Add(*k_fs_inadd, false, -1.0, 1.0);
  k_fs_inmod->Complete(k_fs_in->DomainMap(), k_fs_in->RowMap());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_immod =
      Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs, 100));
  k_fs_immod->Add(*k_fs_im, false, 1.0, 1.0);
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_imadd =
      LINALG::MLMultiply(*fdhat_, true, *k_fs_am, false, false, false, true);
  k_fs_immod->Add(*k_fs_imadd, false, -1.0, 1.0);
  k_fs_immod->Complete(k_fs_im->DomainMap(), k_fs_im->RowMap());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_iimod;
  if (iset)
  {
    k_fs_iimod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs, 100));
    k_fs_iimod->Add(*k_fs_ii, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> k_fs_iiadd =
        LINALG::MLMultiply(*fdhat_, true, *k_fs_ai, false, false, false, true);
    k_fs_iimod->Add(*k_fs_iiadd, false, -1.0, 1.0);
    k_fs_iimod->Complete(k_fs_ii->DomainMap(), k_fs_ii->RowMap());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_iamod;
  if (iset && aset)
  {
    k_fs_iamod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs, 100));
    k_fs_iamod->Add(*k_fs_ia, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> k_fs_iaadd =
        LINALG::MLMultiply(*fdhat_, true, *k_fs_aa, false, false, false, true);
    k_fs_iamod->Add(*k_fs_iaadd, false, -1.0, 1.0);
    k_fs_iamod->Complete(k_fs_ia->DomainMap(), k_fs_ia->RowMap());
  }

  //---------------------------------------------------------- FOURTH LINE
  // nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_anmod;
  if (aset)
    k_fs_anmod = LINALG::MLMultiply(*ftanginvD_, false, *k_fs_an, false, false, false, true);

  // kam: multiply tmatrix with invda and kam
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_ammod;
  if (aset)
    k_fs_ammod = LINALG::MLMultiply(*ftanginvD_, false, *k_fs_am, false, false, false, true);

  // kai: multiply tmatrix with invda and kai
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_aimod;
  if (aset && iset)
    k_fs_aimod = LINALG::MLMultiply(*ftanginvD_, false, *k_fs_ai, false, false, false, true);

  // kaa: multiply tmatrix with invda and kaa
  Teuchos::RCP<LINALG::SparseMatrix> k_fs_aamod;
  if (aset)
    k_fs_aamod = LINALG::MLMultiply(*ftanginvD_, false, *k_fs_aa, false, false, false, true);

  /**********************************************************************/
  /* (8) Build the final f blocks                                       */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // fm: add alphaf * old contact forces (t_n)
  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (IsSelfContact())
  {
    dserror("CHECK ME!");
    //        Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    //        Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new
    //        Epetra_Vector(mold_->DomainMap())); Teuchos::RCP<Epetra_Vector> zoldexp  =
    //        Teuchos::rcp(new Epetra_Vector(mold_->RowMap())); if
    //        (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    //        mold_->Multiply(true,*zoldexp,*tempvecm2);
    //        if (mset) LINALG::Export(*tempvecm2,*tempvecm);
    //        fm->Update(alphaf_,*tempvecm,1.0);
  }
  // if there is no self contact everything is ok
  else if (poromaster_)
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    fmoldtransp_->Multiply(false, *lambdaold_, *tempvecm);
    fm->Update(nopenalpha_, *tempvecm, 1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (IsSelfContact())
  {
    dserror("CHECK ME!");
    //        Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new
    //        Epetra_Vector(dold_->DomainMap())); Teuchos::RCP<Epetra_Vector> zoldexp  =
    //        Teuchos::rcp(new Epetra_Vector(dold_->RowMap())); if
    //        (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    //        dold_->Multiply(true,*zoldexp,*tempvec);
    //        if (sset) LINALG::Export(*tempvec,*fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    fdoldtransp_->Multiply(false, *lambdaold_, *fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*fgactivedofs_));
    LINALG::Export(*fsadd, *faadd);

    fa->Update(-nopenalpha_, *faadd, 1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod;
  if (mset)
  {
    fmmod = Teuchos::rcp(new Epetra_Vector(*fgmdofrowmap_));
    if (aset) fmhataam_->Multiply(true, *fa, *fmmod);
    fmmod->Update(1.0, *fm, 1.0);
  }

  //----------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*fgidofs));
    LINALG::Export(*fsadd, *fiadd);
    fi->Update(-nopenalpha_, *fiadd, 1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*fgidofs));
  if (aset) fdhat_->Multiply(true, *fa, *fimod);
  fimod->Update(1.0, *fi, -1.0);

  //---------------------------------------------------------- FOURTH LINE
  // gactive: nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // fa: multiply tmatrix with invda and fa
  Teuchos::RCP<Epetra_Vector> famod;
  if (aset)
  {
    famod = Teuchos::rcp(new Epetra_Vector(*fgactivet_));
    ftanginvD_->Multiply(false, *fa, *famod);
  }
  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (ParRedist())
  {
    dserror("CHECK ME!");
    //        //----------------------------------------------------------- FIRST LINE
    //        // nothing to do (ndof-map independent of redistribution)
    //
    //        //---------------------------------------------------------- SECOND LINE
    //        kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    //        kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    //        if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
    //        if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);
    //
    //        //----------------------------------------------------------- THIRD LINE
    //        if (iset)
    //        {
    //          kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
    //          kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
    //          kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
    //          if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
    //        }
    //
    //        //---------------------------------------------------------- FOURTH LINE
    //        if (aset) smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);
    //
    //        //----------------------------------------------------------- FIFTH LINE
    //        if (aset)
    //        {
    //          kanmod = MORTAR::MatrixRowTransform(kanmod,pgsdofrowmap_);
    //          kammod = MORTAR::MatrixRowTransform(kammod,pgsdofrowmap_);
    //          kaamod = MORTAR::MatrixRowTransform(kaamod,pgsdofrowmap_);
    //          if (iset) kaimod = MORTAR::MatrixRowTransform(kaimod,pgsdofrowmap_);
    //          tderivmatrix_ = MORTAR::MatrixRowTransform(tderivmatrix_,pgsdofrowmap_);
    //        }
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> k_fs_effnew = Teuchos::rcp(
      new LINALG::SparseMatrix(*falldofrowmap_, 81, true, false, k_fseff->GetMatrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*falldofrowmap_);

  //----------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  k_fs_effnew->Add(*k_fs_nn, false, 1.0, 1.0);
  k_fs_effnew->Add(*k_fs_nm, false, 1.0, 1.0);
  if (sset) k_fs_effnew->Add(*k_fs_ns, false, 1.0, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  if (mset)
  {
    k_fs_effnew->Add(*k_fs_mnmod, false, 1.0, 1.0);
    k_fs_effnew->Add(*k_fs_mmmod, false, 1.0, 1.0);
    if (iset) k_fs_effnew->Add(*k_fs_mimod, false, 1.0, 1.0);
    if (aset) k_fs_effnew->Add(*k_fs_mamod, false, 1.0, 1.0);
  }

  //----------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) k_fs_effnew->Add(*k_fs_inmod, false, 1.0, 1.0);
  if (iset) k_fs_effnew->Add(*k_fs_immod, false, 1.0, 1.0);
  if (iset) k_fs_effnew->Add(*k_fs_iimod, false, 1.0, 1.0);
  if (iset && aset) k_fs_effnew->Add(*k_fs_iamod, false, 1.0, 1.0);

  //---------------------------------------------------------- FOURTH LINE
  // add a submatrices to kteffnew
  if (aset) k_fs_effnew->Add(*fNCoup_lindisp_, false, 1.0, 1.0);

  //----------------------------------------------------------- FIFTH LINE
  // add a submatrices to kteffnew
  if (aset) k_fs_effnew->Add(*k_fs_anmod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->Add(*k_fs_ammod, false, -1.0, 1.0);
  if (aset && iset) k_fs_effnew->Add(*k_fs_aimod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->Add(*k_fs_aamod, false, -1.0, 1.0);
  if (aset) k_fs_effnew->Add(*flinTangentiallambda_, false, 1.0, 1.0);

  // FillComplete kteffnew (square)
  k_fs_effnew->Complete(*ProblemDofs(), *falldofrowmap_);
  /**********************************************************************/
  /* (11) Global setup of feffnew (including contact)                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));

  LINALG::Export(*fn, *fnexp);

  feffnew->Update(1.0, *fnexp, 1.0);
  //---------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  if (mset)
  {
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    LINALG::Export(*fmmod, *fmmodexp);
    feffnew->Update(1.0, *fmmodexp, 1.0);
  }
  //----------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    LINALG::Export(*fimod, *fimodexp);
    feffnew->Update(1.0, *fimodexp, 1.0);
  }

  //---------------------------------------------------------- FOURTH LINE
  // add weighted nCoup vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> nCoupexp;
  if (aset)
  {
    nCoupexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    LINALG::Export(*fNCoup_, *nCoupexp);
    feffnew->Update(-1.0, *nCoupexp, 1.0);
  }

  //----------------------------------------------------------- FIFTH LINE
  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> famodexp;
  if (aset)
  {
    famodexp = Teuchos::rcp(new Epetra_Vector(*falldofrowmap_));
    LINALG::Export(*famod, *famodexp);
    feffnew->Update(-1.0, *famodexp, 1.0);
  }

  // finally do the replacement
  k_fseff = k_fs_effnew;
  feff = feffnew;
  return;
}  // CONTACT::PoroLagrangeStrategy::EvaluatePoroNoPen

/*----------------------------------------------------------------------*
 |  evaluate poro coupling contact matrices for no penetration          |
 |  condition on contact surface (public)                     ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::EvaluateOtherMatPoroNoPen(
    Teuchos::RCP<LINALG::SparseMatrix>& Feff, int Column_Block_Id)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ || (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()))
    return;
  // this method should be renamed as it handles twosided poro meshtying now aswell and is
  // able to handle fluid coupling for twosided contact

  nopenalpha_ = alphaf_;  // to use different alpha for nopen condition (not used at the moment)

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  Feff->Complete();

  Teuchos::RCP<Epetra_Map> domainmap = Teuchos::rcp(new Epetra_Map(Feff->DomainMap()));

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************

  // double-check if this is a dual LM system
  if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    dserror("Condensation only for dual LM");

  /**********************************************************************/
  /* (3) Split k_fseff and Feff into 3x3 matrix blocks                             */
  /**********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> F_s, F_m, F_n;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> F_sm, F_sm0, F_n0, F_m0, F_s0;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap0;
  Teuchos::RCP<Epetra_Map> tempmap1;
  Teuchos::RCP<Epetra_Map> ftempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<LINALG::SparseMatrix> Feffmatrix =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(Feff);
  if (ParRedist())
  {
    dserror("CHECK ME!");
    // split and transform to redistributed maps
    //        LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    //        ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    //        ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    //        knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitMatrix2x2(
        Feffmatrix, fgsmdofrowmap_, fgndofrowmap_, domainmap, tempmap0, F_sm, F_sm0, F_n, F_n0);
  }

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(
      F_sm, fgsdofrowmap_, fgmdofrowmap_, domainmap, tempmap0, F_s, F_s0, F_m, F_m0);

  // store some stuff for static condensation of LM
  cfx_s_.insert(std::pair<int, Teuchos::RCP<LINALG::SparseMatrix>>(Column_Block_Id, F_s));

  /**********************************************************************/
  /* (5) Split slave quantities into active / inactive                  */
  /**********************************************************************/

  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> F_a, F_a0, F_i, F_i0;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> F_an, F_in, F_am, F_im, F_ma, F_mi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> fgidofs;

  LINALG::SplitMatrix2x2(F_s, fgactivedofs_, fgidofs, domainmap, tempmap1, F_a, F_a0, F_i, F_i0);

  // abbreviations for active and inactive set
  int aset = fgactivedofs_->NumGlobalElements();
  int iset = fgidofs->NumGlobalElements();

  // abbreviations for slave  and master set
  // int sset = fgsdofrowmap_->NumGlobalElements(); //usally slave should anyway exist!
  int mset = fgmdofrowmap_->NumGlobalElements();

  /**********************************************************************/
  /* (7) Build the final K blocks                                       */
  /* where K stands for k_fs and F!!!                                   */
  /**********************************************************************/

  //----------------------------------------------------------- FIRST LINE
  // kn: nothing to do
  //---------------------------------------------------------- SECOND LINE --- Will just exist when
  // starting with two sided poro contact!!!
  // km: add T(mhataam)*kan
  Teuchos::RCP<LINALG::SparseMatrix> F_mmod =
      Teuchos::rcp(new LINALG::SparseMatrix(*fgmdofrowmap_, 100));
  F_mmod->Add(*F_m, false, 1.0, 1.0);
  if (aset && mset)
  {
    Teuchos::RCP<LINALG::SparseMatrix> F_madd =
        LINALG::MLMultiply(*fmhataam_, true, *F_a, false, false, false, true);
    F_mmod->Add(*F_madd, false, 1.0, 1.0);
  }
  F_mmod->Complete(F_m->DomainMap(), F_m->RowMap());

  //----------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------

  //--- For using non diagonal D-Matrix, it should be checked if this assumtion isn't anywhere
  // else!!!

  // kin: subtract T(dhat)*kan --
  Teuchos::RCP<LINALG::SparseMatrix> F_imod = Teuchos::rcp(new LINALG::SparseMatrix(*fgidofs, 100));
  F_imod->Add(*F_i, false, 1.0, 1.0);
  if (aset)
  {
    Teuchos::RCP<LINALG::SparseMatrix> F_iadd =
        LINALG::MLMultiply(*fdhat_, true, *F_a, false, false, false, true);
    F_imod->Add(*F_iadd, false, -1.0, 1.0);
  }
  F_imod->Complete(F_i->DomainMap(), F_i->RowMap());

  //---------------------------------------------------------- FOURTH LINE
  // nothing to do
  //----------------------------------------------------------- FIFTH LINE
  // kan: multiply tmatrix with invda and kan
  Teuchos::RCP<LINALG::SparseMatrix> F_amod;
  if (aset)
  {
    F_amod = LINALG::MLMultiply(*ftanginvD_, false, *F_a, false, false, false, true);
  }

  /********************************************************************/
  /* (9) Transform the final K blocks                                 */
  /********************************************************************/
  // The row maps of all individual matrix blocks are transformed to
  // the parallel layout of the underlying problem discretization.
  // Of course, this is only necessary in the parallel redistribution
  // case, where the contact interfaces have been redistributed
  // independently of the underlying problem discretization.

  if (ParRedist())
  {
    dserror("CHECK ME!");
    //        //----------------------------------------------------------- FIRST LINE
    //        // nothing to do (ndof-map independent of redistribution)
    //
    //        //---------------------------------------------------------- SECOND LINE
    //        kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    //        kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    //        if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
    //        if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);
    //
    //        //----------------------------------------------------------- THIRD LINE
    //        if (iset)
    //        {
    //          kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
    //          kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
    //          kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
    //          if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
    //        }
    //
    //        //---------------------------------------------------------- FOURTH LINE
    //        if (aset) smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);
    //
    //        //----------------------------------------------------------- FIFTH LINE
    //        if (aset)
    //        {
    //          kanmod = MORTAR::MatrixRowTransform(kanmod,pgsdofrowmap_);
    //          kammod = MORTAR::MatrixRowTransform(kammod,pgsdofrowmap_);
    //          kaamod = MORTAR::MatrixRowTransform(kaamod,pgsdofrowmap_);
    //          if (iset) kaimod = MORTAR::MatrixRowTransform(kaimod,pgsdofrowmap_);
    //          pmatrix_ = MORTAR::MatrixRowTransform(tderivmatrix_,pgsdofrowmap_);
    //        }
  }
  /**********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                  */
  /**********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> F_effnew = Teuchos::rcp(
      new LINALG::SparseMatrix(*falldofrowmap_, 108, true, false, Feff->GetMatrixtype()));

  //----------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  F_effnew->Add(*F_n, false, 1.0, 1.0);

  //---------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  if (mset) F_effnew->Add(*F_mmod, false, 1.0, 1.0);
  //----------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) F_effnew->Add(*F_imod, false, 1.0, 1.0);
  //---------------------------------------------------------- FOURTH LINE
  // add a submatrices to kteffnew //assume that Column_Block_Id==0 is the porofluid coupling
  // block!!!
  if (Column_Block_Id == 0 && aset) F_effnew->Add(*fNCoup_linvel_, false, 1.0, 1.0);
  //----------------------------------------------------------- FIFTH LINE
  // add a submatrices to kteffnew
  if (aset) F_effnew->Add(*F_amod, false, -1.0, 1.0);

  // FillComplete kteffnew (square)
  F_effnew->Complete(*domainmap, *falldofrowmap_);

  // finally do the replacement
  Feff = F_effnew;
  return;
}  // CONTACT::PoroLagrangeStrategy::EvaluateOtherMatPoroNoPen

/*----------------------------------------------------------------------------*
 | Poro Recovery method for no penetration LM (pure porous problem) ager 08/14|
 *---------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::RecoverPoroNoPen(
    Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> inc)
{
  std::map<int, Teuchos::RCP<Epetra_Vector>> incm;
  incm.insert(std::pair<int, Teuchos::RCP<Epetra_Vector>>(0, inc));

  RecoverPoroNoPen(disi, incm);
  return;
}

/*----------------------------------------------------------------------*
 | Poro Recovery method for no penetration LM                 ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::RecoverPoroNoPen(
    Teuchos::RCP<Epetra_Vector> disi, std::map<int, Teuchos::RCP<Epetra_Vector>> inc)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!no_penetration_ || (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()))
    return;

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

  //**********************************************************************
  //**********************************************************************
  // CASE: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> finvda;
    Teuchos::RCP<Epetra_Map> tempmap1, tempmap2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(finvda_, fgactivedofs_, tempmap1, gactivedofs_, tempmap2, finvda,
        tempmtx1, tempmtx2, tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> finvdmod =
        Teuchos::rcp(new LINALG::SparseMatrix(*fgsdofrowmap_, 10));
    finvdmod->Add(*finvda, false, 1.0, 1.0);
    finvdmod->Complete(*gsdofrowmap_, *fgsdofrowmap_);

    /**********************************************************************/
    /* Update Lagrange multipliers lambda_n+1                                  */
    /**********************************************************************/
    {
      Teuchos::RCP<Epetra_Vector> flambda = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_, true));

      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*fgsdofrowmap_));

      cfssn_->Multiply(false, *disin, *mod);
      flambda->Update(-1.0, *mod, 0.0);
      cfssm_->Multiply(false, *disim, *mod);
      flambda->Update(-1.0, *mod, 1.0);
      cfsss_->Multiply(false, *disis, *mod);
      flambda->Update(-1.0, *mod, 1.0);

      // loop over all offdiag blocks!!!
      std::map<int, Teuchos::RCP<LINALG::SparseOperator>>::iterator matiter;
      std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator inciter;
      for (matiter = cfx_s_.begin(); matiter != cfx_s_.end(); ++matiter)
      {
        inciter = inc.find(matiter->first);
        if (inciter == inc.end())
          dserror(
              "CONTACT::PoroLagrangeStrategy::RecoverPoroNoPen: Couldn't find increment block %d "
              "for recovery of the lagrange multiplier!",
              matiter->first);

        matiter->second->Multiply(false, *inciter->second, *mod);
        flambda->Update(-1.0, *mod, 1.0);
      }

      flambda->Update(1.0, *ffs_, 1.0);



      fdoldtransp_->Multiply(false, *lambdaold_, *mod);
      flambda->Update(-nopenalpha_, *mod, 1.0);

      Teuchos::RCP<Epetra_Vector> lambdacopy = Teuchos::rcp(new Epetra_Vector(*flambda));

      finvdmod->Multiply(true, *lambdacopy, *lambda_);  // should be lambda_ at the end!!!

      lambda_->Scale(
          (1 - alphaf_) / (1 - nopenalpha_));  //-- is already scaled by this factor by scaling
                                               // invda_!!! --- scale it back to with nopenalpha_...
    }
  }
  // store updated LM into nodes
  SetState(MORTAR::state_lagrange_multiplier, *lambda_);

  return;
}

/*------------------------------------------------------------------------------*
 |  Additional update and output poro contact at end of time step     ager 08/14|
 *-----------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::UpdatePoroContact()
{
  if (no_penetration_)
  {
    // std::cout << "print lambda: " << *lambda_ << std::endl;
    lambdaold_->Update(1.0, *lambda_, 0.0);
  }
  return;
}

/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                          ager 08/14|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SetState(
    const enum MORTAR::StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case MORTAR::state_fvelocity:
    case MORTAR::state_svelocity:
    case MORTAR::state_lagrange_multiplier:
    case MORTAR::state_fpressure:
    {
      // set state on interfaces
      for (int i = 0; i < (int)interface_.size(); ++i)
      {
        // interface_[i]->SetState(statename, vec);
        DRT::Discretization& idiscret_ = interface_[i]->Discret();

        switch (statetype)
        {
          case MORTAR::state_fvelocity:
          case MORTAR::state_svelocity:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.DofColMap(), true));

            LINALG::Export(vec, *global);

            // loop over all nodes to set current velocity
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.NumMyColNodes(); ++i)
            {
              CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_.lColNode(i));
              const int numdof = node->NumDof();
              std::vector<double> myvel(numdof);
              std::vector<int> lm(numdof);

              for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];
              DRT::UTILS::ExtractMyValues(*global, myvel, lm);

              // add myvel[2]=0 for 2D problems
              if (myvel.size() < 3) myvel.resize(3);
              // set current configuration
              for (int j = 0; j < 3; ++j)
              {
                if (statetype == MORTAR::state_fvelocity)
                  node->CoPoroData().fvel()[j] = myvel[j];
                else
                  node->CoPoroData().svel()[j] = myvel[j];
              }
            }
            break;
          }
          case MORTAR::state_lagrange_multiplier:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.DofColMap(), true));
            LINALG::Export(vec, *global);

            // loop over all nodes to set current velocity
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.NumMyColNodes(); ++i)
            {
              CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_.lColNode(i));

              const int numdof = node->NumDof();
              std::vector<double> mylm(numdof);
              std::vector<int> lm(numdof);

              for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

              DRT::UTILS::ExtractMyValues(*global, mylm, lm);

              // add myvel[2]=0 for 2D problems
              if (mylm.size() < 3) mylm.resize(3);
              // set current configuration
              if (node->IsSlave())
              {
                for (int j = 0; j < 3; ++j)
                {
                  node->CoPoroData().poroLM()[j] = mylm[j];
                }
              }
            }
            break;
          }
          case MORTAR::state_fpressure:
          {
            // alternative method to get vec to full overlap
            Teuchos::RCP<Epetra_Vector> global =
                Teuchos::rcp(new Epetra_Vector(*idiscret_.DofColMap(), true));
            LINALG::Export(vec, *global);

            // loop over all nodes to set current pressure
            // (use fully overlapping column map)
            for (int i = 0; i < idiscret_.NumMyColNodes(); ++i)
            {
              CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_.lColNode(i));

              double myfpres;
              int fpres;

              fpres = node->Dofs()[0];  // here get ids of first component of node

              myfpres = global->Values()[global->Map().LID(fpres)];

              *node->CoPoroData().fpres() = myfpres;
            }
            break;
          }
          default:
          {
            dserror("Shouldn't happen!");
            break;
          }
        }  // end inner switch statement
      }    // end loop over all interfaces
      break;
    }
    default:
    {
      CONTACT::CoAbstractStrategy::SetState(statetype, vec);
      break;
    }
  }  // end outer switch statement
  return;
}


// this should add the displacement of the parent element to the contact element into the
// datacontainer...
//...and it should additionally add the dof IDs to its Data container
/*------------------------------------------------------------------------*
 | Assign generell poro contact state!                          ager 10/14|
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::SetParentState(const std::string& statename,
    const Teuchos::RCP<Epetra_Vector> vec, const Teuchos::RCP<DRT::Discretization> dis)
{
  if (statename == "displacement")
  {
    Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap(), true));
    LINALG::Export(*vec, *global);

    // set state on interfaces
    for (int i = 0; i < (int)interface_.size(); ++i)
    {
      DRT::Discretization& idiscret_ = interface_[i]->Discret();

      if (poroslave_)
      {
        for (int j = 0; j < interface_[i]->SlaveColElements()->NumMyElements();
             ++j)  // will just work for onesided poro contact as the porosity is just on slave
                   // side!!!
        {
          int gid = interface_[i]->SlaveColElements()->GID(j);

          MORTAR::MortarElement* ele =
              dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          // this gets values in local order
          ele->ParentElement()->LocationVector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*global, myval, lm);

          ele->MoData().ParentDisp() = myval;
          ele->MoData().ParentDof() = lm;
        }
      }
      if (poromaster_)  // add master parent element displacements
      {
        for (int j = 0; j < interface_[i]->MasterColElements()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->MasterColElements()->GID(j);

          MORTAR::MortarElement* mele =
              dynamic_cast<MORTAR::MortarElement*>(idiscret_.gElement(gid));

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;

          // this gets values in local order
          mele->ParentElement()->LocationVector(*dis, lm, lmowner, lmstride);

          std::vector<double> myval;
          DRT::UTILS::ExtractMyValues(*global, myval, lm);

          mele->MoData().ParentDisp() = myval;
          mele->MoData().ParentDof() = lm;
        }
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | Initialize poro meshtying matrices with corresponding maps  h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroMtInitialize()
{
  // (re)setup global lin of D & M * lambda - Matrix
  porolindmatrix_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*gsdofrowmap_, 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
  porolinmmatrix_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*gmdofrowmap_, 100, true, false, LINALG::SparseMatrix::FE_MATRIX));

  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 10));
  mmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_, 100));

  // this is needed because on the mesh tying way to build this strategy it is not done elsewhere
  SetupNoPenetrationCondition();

  isincontact_ = true;      // simply set true for meshtying
  wasincontact_ = true;     // as meshtying interfaces stay the same and are fully in contact
  wasincontactlts_ = true;  // this is necessary for other methods in this strategy
  return;
}  // CONTACT::PoroLagrangeStrategy::PoroLinkDM

/*------------------------------------------------------------------------*
 | assemble porofluid meshtying matrices                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroMtPrepareFluidCoupling()
{
  // reset D and M matrices for new Newton iteration

  dmatrix_->Zero();
  mmatrix_->Zero();

  // Assemble D and M matrices
  for (int i = 0; i < (int)interface_.size(); ++i)
  {
    interface_[i]->AssembleDM(*dmatrix_, *mmatrix_);
  }

  // complete D and M matrices
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);

  // as mhataam-, dhat_ and invda_ are not computed in poro - meshtying before this point it is
  // necessary here
  PoroMtSetCouplingMatrices();

  return;
}  // CONTACT::PoroLagrangeStrategy::PoroAssembleFluidCoupling()

/*------------------------------------------------------------------------*
 | assemble porofluid meshtying matrices                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroMtSetCouplingMatrices()
{
  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

  Teuchos::RCP<Epetra_Map> gidofs;

  int aset = gactivedofs_->NumGlobalElements();

  // invert dmatrix to invd
  invd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_, true);
  int err = 0;

  // extract diagonal of invd into diag
  invd_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  //  for (int i = 0; i < diag->MyLength(); ++i)
  //    if ((*diag)[i] == 0.0)
  //      (*diag)[i] = 1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err > 0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd_->ReplaceDiagonalValues(*diag);
  // inversion end

  // active part of invd
  Teuchos::RCP<LINALG::SparseMatrix> invda;
  LINALG::SplitMatrix2x2(
      invd_, gactivedofs_, gidofs, gactivedofs_, gidofs, invda, tempmtx1, tempmtx2, tempmtx3);

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)(not considered for poro
  // mt/contact yet-> ==0)
  Teuchos::RCP<LINALG::SparseMatrix> dai;
  LINALG::SplitMatrix2x2(
      dmatrix_, gactivedofs_, gidofs, gactivedofs_, gidofs, tempmtx1, dai, tempmtx2, tempmtx3);

  int iset = gidofs->NumGlobalElements();
  // do the multiplication dhat = invda * dai
  Teuchos::RCP<LINALG::SparseMatrix> dhat =
      Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_, 10));
  if (aset && iset) dhat = LINALG::MLMultiply(*invda, false, *dai, false, false, false, true);
  dhat->Complete(*gidofs, *gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixa;
  LINALG::SplitMatrix2x2(mmatrix_, gactivedofs_, gidofs, gmdofrowmap_, tempmap, mmatrixa, tempmtx1,
      tempmtx2, tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<LINALG::SparseMatrix> mhataam =
      Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_, 10));
  if (aset) mhataam = LINALG::MLMultiply(*invda, false, *mmatrixa, false, false, false, true);
  mhataam->Complete(*gmdofrowmap_, *gactivedofs_);

  // scaling of invd and dai
  invda->Scale(1 / (1 - alphaf_));

  SaveCouplingMatrices(dhat, mhataam, invda);

  return;
}  // CONTACT::PoroLagrangeStrategy::PoroMtPrepareFluidCoupling

/*------------------------------------------------------------------------*
 | update old meshtying matrices and LMP                       h.Willmann |
 *------------------------------------------------------------------------*/
void CONTACT::PoroLagrangeStrategy::PoroMtUpdate()
{
  // set dold_ mold_ and lambdaold_ after every time step
  dold_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*dmatrix_));
  mold_ = Teuchos::rcp<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*mmatrix_));
  dold_->Complete(dmatrix_->DomainMap(), dmatrix_->RangeMap());
  mold_->Complete(mmatrix_->DomainMap(), mmatrix_->RangeMap());
  UpdatePoroContact();
  return;
}  // CONTACT::PoroLagrangeStrategy::PoroMtUpdate()
