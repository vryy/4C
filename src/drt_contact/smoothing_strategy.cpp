/*---------------------------------------------------------------------*/
/*!
\file smoothing_strategy.cpp

\brief strategy for interface smoothing

\level 3

\maintainer Alexander Popp

*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 01/15|
 *----------------------------------------------------------------------*/
#include "smoothing_strategy.H"
#include "contact_interface.H"
#include "contact_defines.H"
#include "friction_node.H"

#include "../drt_mortar/mortar_interface.H"

#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             farah 01/15|
 *----------------------------------------------------------------------*/
CONTACT::SmoothingStrategy::SmoothingStrategy(
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > cinterface,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface> > mtinterface,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoAbstractStrategy(Teuchos::rcp(new CONTACT::AbstractStratDataContainer()),
        DofRowMap, NodeRowMap, params, dim, comm, alphaf, maxdof),
      activesetssconv_(false),
      activesetconv_(false),
      activesetsteps_(1),
      interface_(cinterface),
      minterface_(mtinterface)
{
  // get number of contact lm dofs as offset
  globaloffset_ = glmdofrowmap_->NumGlobalElements();

  // setup maps
  SetupMT(false);

  // evaluate mortar coupling
  MortarCoupling(Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::SmoothingStrategy::SmoothingStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > cinterface,
    std::vector<Teuchos::RCP<MORTAR::MortarInterface> > mtinterface,
    int dim,
    Teuchos::RCP<Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoAbstractStrategy(data_ptr, DofRowMap, NodeRowMap, params,
        dim, comm, alphaf, maxdof),
      activesetssconv_(false),
      activesetconv_(false),
      activesetsteps_(1),
      interface_(cinterface),
      minterface_(mtinterface)
{
  // get number of contact lm dofs as offset
  globaloffset_ = glmdofrowmap_->NumGlobalElements();

  // setup maps
  SetupMT(false);

  // evaluate mortar coupling
  MortarCoupling(Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | setup this strategy object                               farah 01/15 |
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::SetupMT(bool redistributed)
{
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  MT_gsdofrowmap_ = Teuchos::null;
  MT_gmdofrowmap_ = Teuchos::null;
  MT_gsmdofrowmap_ = Teuchos::null;
  MT_glmdofrowmap_ = Teuchos::null;
  MT_gsnoderowmap_ = Teuchos::null;
  MT_gmnoderowmap_ = Teuchos::null;
  gdisprowmap_ = Teuchos::null;

  // reset global neutral map
  if (!redistributed)
    gndofrowmap_ = Teuchos::null;

  // merge interface maps to global maps
  for (int i = 0; i < (int) minterface_.size(); ++i) {
    // build Lagrange multiplier dof map
    minterface_[i]->UpdateLagMultSets(globaloffset_);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    MT_glmdofrowmap_ = LINALG::MergeMap(MT_glmdofrowmap_,
        minterface_[i]->LagMultDofs());
    globaloffset_ += MT_glmdofrowmap_->NumGlobalElements();
    if (globaloffset_ < 0)
      globaloffset_ = 0;

    // merge interface master, slave maps to global master, slave map
    MT_gsdofrowmap_ = LINALG::MergeMap(MT_gsdofrowmap_,
        minterface_[i]->SlaveRowDofs());
    MT_gmdofrowmap_ = LINALG::MergeMap(MT_gmdofrowmap_,
        minterface_[i]->MasterRowDofs());
    MT_gsnoderowmap_ = LINALG::MergeMap(MT_gsnoderowmap_,
        minterface_[i]->SlaveRowNodes());
    MT_gmnoderowmap_ = LINALG::MergeMap(MT_gmnoderowmap_,
        minterface_[i]->MasterRowNodes());
  }

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the dicretization dof map)
  // (no need to rebuild this map after redistribution)

  if (!redistributed) {
    gndofrowmap_ = LINALG::SplitMap(*ProblemDofs(), *MT_gsdofrowmap_);
    gndofrowmap_ = LINALG::SplitMap(*gndofrowmap_, *MT_gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  MT_gsmdofrowmap_ = LINALG::MergeMap(*MT_gsdofrowmap_, *MT_gmdofrowmap_,
      false);
  gdisprowmap_ = LINALG::MergeMap(*gndofrowmap_, *MT_gsmdofrowmap_, false);

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // setup Lagrange multiplier vectors
  MT_z_ = Teuchos::rcp(new Epetra_Vector(*MT_gsdofrowmap_));
  MT_zincr_ = Teuchos::rcp(new Epetra_Vector(*MT_gsdofrowmap_));
  MT_zold_ = Teuchos::rcp(new Epetra_Vector(*MT_gsdofrowmap_));

  // setup constraint rhs vector
  //  constrrhs_ = Teuchos::null; // only for saddle point problem formulation

  return;
}


/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step  farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::Initialize()
{
  if (friction_)
    dserror("ERROR: Frictional contact not implemented for interface smoothing!");
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    dserror("ERROR: constraint xyz orientation not implemented for interface smoothing!");

  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*gsdofrowmap_, 100, true, false,
          LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*gmdofrowmap_, 100, true, false,
          LINALG::SparseMatrix::FE_MATRIX));

  if (stype_ == INPAR::CONTACT::solution_lagmult)
  {
    // (re)setup global tangent matrix
    tmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 3));

    // (re)setup global matrix containing gap derivatives
    smatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_, 3));

    // inactive rhs for the saddle point problem
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,
        *gactivedofs_);
    inactiverhs_ = LINALG::CreateVector(*gidofs, true);

    // tangential rhs
    tangrhs_ = LINALG::CreateVector(*gactivet_, true);
    tderivmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_, 3));
  }

  if (stype_ == INPAR::CONTACT::solution_penalty)
  {
    // (re)setup global vector containing lagrange multipliers
    z_ = LINALG::CreateVector(*gsdofrowmap_, true);

    // (re)setup global matrix containing lagrange multiplier derivatives
    linzmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate frictional contact (public)                      farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::EvaluateFriction(
    Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff)
{
  dserror("ERROR: Frictional contact not implemented for interface smoothing!");
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate contact (public)                                farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::EvaluateContact(
    Teuchos::RCP<LINALG::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up

  // system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::SystemType>(Params(), "SYSTEM");

  // shape function
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");

  if (stype_ == INPAR::CONTACT::solution_penalty)
  {
    bool isincontact = false;
    bool activesetchange = false;

    for (int i=0; i<(int)interface_.size(); ++i)
    {
      bool localisincontact = false;
      bool localactivesetchange = false;

      // evaluate lagrange multipliers (regularized forces) in normal direction
      // and nodal derivz matrix values, store them in nodes
      interface_[i]->AssembleRegNormalForces(localisincontact, localactivesetchange);

      isincontact = isincontact || localisincontact;
      activesetchange = activesetchange || localactivesetchange;
    }

    // broadcast contact status & active set change
    int globalcontact, globalchange = 0;
    int localcontact = isincontact;
    int localchange = activesetchange;

    Comm().SumAll(&localcontact, &globalcontact, 1);
    Comm().SumAll(&localchange, &globalchange, 1);

    if (globalcontact>=1)
    {
      isincontact_=true;
      wasincontact_=true;
    }
    else
      isincontact_=false;

    if( (Comm().MyPID()==0) && (globalchange>=1) )
      std::cout << "ACTIVE CONTACT SET HAS CHANGED..." << std::endl;

    // (re)setup active global Epetra_Maps
    // the map of global active nodes is needed for the penalty case, too.
    // this is due to the fact that we want to monitor the constraint norm
    // of the active nodes
    gactivenodes_ = Teuchos::null;
    gslipnodes_ = Teuchos::null;
    gactivedofs_=Teuchos::null;

    // update active sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    for (int i=0;i<(int)interface_.size();++i)
    {
      interface_[i]->BuildActiveSet();
      gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
      gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
      gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
    }

    // check if contact contributions are present,
    // if not we can skip this routine to speed things up
    if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

    // since we will modify the graph of kteff by adding additional
    // meshtyong stiffness entries, we have to uncomplete it
    kteff->UnComplete();

    // assemble contact quantities on all interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      // assemble global lagrangian multiplier vector
      interface_[i]->AssembleLM(*z_);
      // assemble global derivatives of lagrangian multipliers
      interface_[i]->AssembleLinZ(*linzmatrix_);
      // assemble global derivatives of mortar D and M matrices
      interface_[i]->AssembleLinDM(*lindmatrix_, *linmmatrix_);
    }

    // FillComplete() global matrices LinD, LinM, LinZ
    lindmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);
    linmmatrix_->Complete(*gsmdofrowmap_, *gmdofrowmap_);
    linzmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);

    //----------------------------------------------------------------------
    // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
    //----------------------------------------------------------------------
    // Concretely, we apply the following transformations:
    // LinD      ---->   T^(-T) * LinD
    // D         ---->   D * T^(-1)
    //----------------------------------------------------------------------

    // **********************************************************************
    // Build Contact Stiffness #1
    // **********************************************************************
    // involving contributions of derivatives of D and M:
    //  Kc,1 = delta[ 0 -M(transpose) D] * LM

    // transform if necessary
    if (ParRedist())
    {
      lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
      linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
    }

    // add to kteff
    kteff->Add(*lindmatrix_, false, 1.0-alphaf_, 1.0);
    kteff->Add(*linmmatrix_, false, 1.0-alphaf_, 1.0);

    // **********************************************************************
    // Build Contact Stiffness #2
    // **********************************************************************
    // involving contributions of derivatives of lagrange multipliers:
    //  Kc,2= [ 0 -M(transpose) D] * deltaLM

    // multiply Mortar matrices D and M with LinZ
    Teuchos::RCP<LINALG::SparseMatrix> dtilde = LINALG::MLMultiply(*dmatrix_,true,*linzmatrix_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> mtilde = LINALG::MLMultiply(*mmatrix_,true,*linzmatrix_,false,false,false,true);

    // transform if necessary
    if (ParRedist())
    {
      dtilde = MORTAR::MatrixRowTransform(dtilde,pgsdofrowmap_);
      mtilde = MORTAR::MatrixRowTransform(mtilde,pgmdofrowmap_);
    }

    // add to kteff
    kteff->Add(*dtilde, false, 1.0-alphaf_, 1.0);
    kteff->Add(*mtilde, false, -(1.0-alphaf_), 1.0);

    // **********************************************************************
    // Build RHS
    // **********************************************************************
    // feff += -alphaf * fc,n - (1-alphaf) * fc,n+1,k

    {
      // we initialize fcmdold with dold-rowmap instead of gsdofrowmap
      // (this way, possible self contact is automatically included)

      Teuchos::RCP<Epetra_Vector> fcmdold = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
      dold_->Multiply(true, *zold_, *fcmdold);
      Teuchos::RCP<Epetra_Vector> fcmdoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fcmdold, *fcmdoldtemp);
      feff->Update(-alphaf_, *fcmdoldtemp, 1.0);
    }

    {
      // we initialize fcmmold with mold-domainmap instead of gmdofrowmap
      // (this way, possible self contact is automatically included)

      Teuchos::RCP<Epetra_Vector> fcmmold = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
      mold_->Multiply(true, *zold_, *fcmmold);
      Teuchos::RCP<Epetra_Vector> fcmmoldtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fcmmold, *fcmmoldtemp);
      feff->Update(alphaf_, *fcmmoldtemp, 1.0);
    }

    {
      Teuchos::RCP<Epetra_Vector> fcmd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      dmatrix_->Multiply(true, *z_, *fcmd);
      Teuchos::RCP<Epetra_Vector> fcmdtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fcmd, *fcmdtemp);
      feff->Update(-(1-alphaf_), *fcmdtemp, 1.0);
    }

    {
      Teuchos::RCP<Epetra_Vector> fcmm = LINALG::CreateVector(*gmdofrowmap_, true);
      mmatrix_->Multiply(true, *z_, *fcmm);
      Teuchos::RCP<Epetra_Vector> fcmmtemp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fcmm, *fcmmtemp);
      feff->Update(1-alphaf_, *fcmmtemp, 1.0);
    }

    //==================== contributions due to meshtying (LM-method)==============

    // add current meshtying forces MT
    Teuchos::RCP<Epetra_Vector> MT_fs = Teuchos::rcp(
        new Epetra_Vector(*MT_gsdofrowmap_));
    MT_dmatrix_->Multiply(true, *MT_z_, *MT_fs);
    Teuchos::RCP<Epetra_Vector> MT_fsexp = Teuchos::rcp(
        new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*MT_fs, *MT_fsexp);
    feff->Update(-(1.0 - alphaf_), *MT_fsexp, 1.0);

    Teuchos::RCP<Epetra_Vector> MT_fm = Teuchos::rcp(
        new Epetra_Vector(*MT_gmdofrowmap_));
    MT_mmatrix_->Multiply(true, *MT_z_, *MT_fm);
    Teuchos::RCP<Epetra_Vector> MT_fmexp = Teuchos::rcp(
        new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*MT_fm, *MT_fmexp);
    feff->Update(1.0 - alphaf_, *MT_fmexp, 1.0);

    // add old meshtying forces (t_n) for mt1 and mt2
    Teuchos::RCP<Epetra_Vector> MT_fsold = Teuchos::rcp(
        new Epetra_Vector(*MT_gsdofrowmap_));
    MT_dmatrix_->Multiply(true, *MT_zold_, *MT_fsold);
    Teuchos::RCP<Epetra_Vector> MT_fsoldexp = Teuchos::rcp(
        new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*MT_fsold, *MT_fsoldexp);
    feff->Update(-alphaf_, *MT_fsoldexp, 1.0);

    Teuchos::RCP<Epetra_Vector> MT_fmold = Teuchos::rcp(
        new Epetra_Vector(*MT_gmdofrowmap_));
    MT_mmatrix_->Multiply(true, *MT_zold_, *MT_fmold);
    Teuchos::RCP<Epetra_Vector> MT_fmoldexp = Teuchos::rcp(
        new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*MT_fmold, *MT_fmoldexp);
    feff->Update(alphaf_, *MT_fmoldexp, 1.0);
  }

  if (stype_ == INPAR::CONTACT::solution_lagmult)
  {
    if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
      return;

    // complete stiffness matrix
    // (this is a prerequisite for the Split2x2 methods to be called later)
    kteff->Complete();

    /**********************************************************************/
    /* export weighted gap vector to gactiveN-map                         */
    /**********************************************************************/
    Teuchos::RCP<Epetra_Vector> gact;
    gact = LINALG::CreateVector(*gactivenodes_, true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_, *gact);
      gact->ReplaceMap(*gactiven_);
    }

    /**********************************************************************/
    /* calculate                                                          */
    /**********************************************************************/
    /* build global matrix t with tangent vectors of active nodes         */
    /* and global matrix s with normal derivatives of active nodes        */
    /* and global matrix p with tangent derivatives of active nodes       */
    /* and inactive right-hand side with old lagrange multipliers (incr)  */
    /* and tangential right-hand side (incr)                              */
    /**********************************************************************/

    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      interface_[i]->AssembleTN(tmatrix_);
      interface_[i]->AssembleS(*smatrix_);
      interface_[i]->AssembleTNderiv(tderivmatrix_);
      interface_[i]->AssembleLinDM(*lindmatrix_, *linmmatrix_);

      if (systype != INPAR::CONTACT::system_condensed)
      {
        interface_[i]->AssembleInactiverhs(*inactiverhs_);
        interface_[i]->AssembleTangrhs(*tangrhs_);
      }
    }

    // FillComplete() global matrix T
    tmatrix_->Complete(*gactivedofs_, *gactivet_);

    // FillComplete() global matrix S
    smatrix_->Complete(*gsmdofrowmap_, *gactiven_);

    // FillComplete() global matrix P
    // (actually gsdofrowmap_ is in general sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    tderivmatrix_->Complete(*gsmdofrowmap_, *gactivet_);

    // FillComplete() global matrices LinD, LinM
    // (again for linD gsdofrowmap_ is sufficient as domain map,
    // but in the edge node modification case, master entries occur!)
    lindmatrix_->Complete(*gsmdofrowmap_, *gsdofrowmap_);
    linmmatrix_->Complete(*gsmdofrowmap_, *gmdofrowmap_);

    //**********************************************************************
    // Lagrange CASE A: CONDENSED SYSTEM (DUAL)
    //**********************************************************************
    if (systype == INPAR::CONTACT::system_condensed)
    {
      dserror("ERROR: Condensed system not implemented!");

      // double-check if this is a dual LM system
      if (shapefcn != INPAR::MORTAR::shape_dual
          && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
        dserror("Condensation only for dual LM");
    }

    //************************************************************************
    // Lagrange CASE B: SADDLE POINT SYSTEM
    //************************************************************************
    else
    {
      // add contact stiffness
      kteff->UnComplete();
      kteff->Add(*lindmatrix_, false, 1.0 - alphaf_, 1.0);
      kteff->Add(*linmmatrix_, false, 1.0 - alphaf_, 1.0);
      kteff->Complete();

      // add contact force terms
      Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      dmatrix_->Multiply(true, *z_, *fs);
      Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fs, *fsexp);
      feff->Update(-(1.0 - alphaf_), *fsexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(
          new Epetra_Vector(*gmdofrowmap_));
      mmatrix_->Multiply(true, *z_, *fm);
      Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fm, *fmexp);
      feff->Update(1.0 - alphaf_, *fmexp, 1.0);

      // add old contact forces (t_n)
      Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(
          new Epetra_Vector(*gsdofrowmap_));
      dold_->Multiply(true, *zold_, *fsold);
      Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fsold, *fsoldexp);
      feff->Update(-alphaf_, *fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(
          new Epetra_Vector(*gmdofrowmap_));
      mold_->Multiply(true, *zold_, *fmold);
      Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fmold, *fmoldexp);
      feff->Update(alphaf_, *fmoldexp, 1.0);

      // add current meshtying forces MT
      Teuchos::RCP<Epetra_Vector> MT_fs = Teuchos::rcp(
          new Epetra_Vector(*MT_gsdofrowmap_));
      MT_dmatrix_->Multiply(true, *MT_z_, *MT_fs);
      Teuchos::RCP<Epetra_Vector> MT_fsexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*MT_fs, *MT_fsexp);
      feff->Update(-(1.0 - alphaf_), *MT_fsexp, 1.0);

      Teuchos::RCP<Epetra_Vector> MT_fm = Teuchos::rcp(
          new Epetra_Vector(*MT_gmdofrowmap_));
      MT_mmatrix_->Multiply(true, *MT_z_, *MT_fm);
      Teuchos::RCP<Epetra_Vector> MT_fmexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*MT_fm, *MT_fmexp);
      feff->Update(1.0 - alphaf_, *MT_fmexp, 1.0);

      // add old meshtying forces (t_n) for mt1 and mt2
      Teuchos::RCP<Epetra_Vector> MT_fsold = Teuchos::rcp(
          new Epetra_Vector(*MT_gsdofrowmap_));
      MT_dmatrix_->Multiply(true, *MT_zold_, *MT_fsold);
      Teuchos::RCP<Epetra_Vector> MT_fsoldexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*MT_fsold, *MT_fsoldexp);
      feff->Update(-alphaf_, *MT_fsoldexp, 1.0);

      Teuchos::RCP<Epetra_Vector> MT_fmold = Teuchos::rcp(
          new Epetra_Vector(*MT_gmdofrowmap_));
      MT_mmatrix_->Multiply(true, *MT_zold_, *MT_fmold);
      Teuchos::RCP<Epetra_Vector> MT_fmoldexp = Teuchos::rcp(
          new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*MT_fmold, *MT_fmoldexp);
      feff->Update(alphaf_, *MT_fmoldexp, 1.0);

    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | Setup saddle point system for contact problems           farah 01/15 |
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::BuildSaddlePointSystem(
    Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd,
    Teuchos::RCP<Epetra_Vector> sold,
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
    Teuchos::RCP<Epetra_Operator>& blockMat,
    Teuchos::RCP<Epetra_Vector>& blocksol,
    Teuchos::RCP<Epetra_Vector>& blockrhs)
{
  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->FullMap())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp( new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp, *dirichtoggle);

  // get system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::SystemType>(Params(), "SYSTEM");

  if (stype_ == INPAR::CONTACT::solution_penalty)
  {
    //**********************************************************************
    // FRICTIONAL CONTACT
    //**********************************************************************
    if (friction_)
      dserror("ERROR: Friction and Smoothing not implemented yet!");

    //**********************************************************************
    // FRICTIONLESS CONTACT
    //**********************************************************************

    // the standard stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiffmt = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kdd);

    // initialize merged system (matrix, rhs, sol)
    Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(ProblemDofs(),MT_glmdofrowmap_,false);
    Teuchos::RCP<LINALG::SparseMatrix> mergedmt = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> mergedrhs = LINALG::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedsol = LINALG::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedzeros = LINALG::CreateVector(*mergedmap);

    // initialize constraint matrices for 2x2 system
    //        [ dd  | dm ]     [ dd | dm ]  (displacement)
    //  k  =  [ md  | mm ]  =  [ md | 0  ]  (meshtying)


    Teuchos::RCP<LINALG::SparseMatrix> kdm = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
    Teuchos::RCP<LINALG::SparseMatrix> kmd = Teuchos::rcp(new LINALG::SparseMatrix(*MT_gsdofrowmap_,100,false,true));

    // initialize transformed constraint matrices
    Teuchos::RCP<LINALG::SparseMatrix> trkdm,trkmd;

    //**********************************************************************
    // build matrix and vector blocks
    //**********************************************************************

    // build constraint matrix kdm
    kdm->Add(*MT_dmatrix_, true, (1.0 - alphaf_), 1.0);
    kdm->Add(*MT_mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdm->Complete(*MT_gsdofrowmap_, *gdisprowmap_);

    // build constraint matrix kmd
    kmd->Add(*MT_dmatrix_, false, 1.0, 1.0);
    kmd->Add(*MT_mmatrix_, false, -1.0, 1.0);
    kmd->Complete(*gdisprowmap_, *MT_gsdofrowmap_);

    //**********************************************************************
    // transform matrix blocks to lm-dof-map
    //**********************************************************************
    //TODO Parallel Redistribution (siehe contact_lagrange_str)

    // transform constraint matrix kmd (MatrixRowTransform)
    trkmd = MORTAR::MatrixRowTransformGIDs(kmd, MT_glmdofrowmap_);

    // transform constraint matrix kdm (MatrixColTransform)
    trkdm = MORTAR::MatrixColTransformGIDs(kdm, MT_glmdofrowmap_);

    //**********************************************************************
    // build and solve saddle point system
    //**********************************************************************

    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*dirichtoggle, *dirichtoggleexp);

    int err=0;

    Teuchos::RCP<Epetra_Vector> MT_lmDBC = Teuchos::rcp(new Epetra_Vector(*MT_gsdofrowmap_));
    LINALG::Export(*dirichtoggle, *MT_lmDBC);

    err=MT_lmDBC->ReplaceMap(*MT_glmdofrowmap_);
    if(err!=0)
      dserror("ReplaceMap failed");

    Teuchos::RCP<Epetra_Vector> MT_lmDBCexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*MT_lmDBC, *MT_lmDBCexp);

    err = dirichtoggleexp->Update(1., *MT_lmDBCexp, 1.);
    if(err!=0)
      dserror("failed");

    // apply Dirichlet conditions to (1,0)
    trkmd->ApplyDirichlet(MT_lmDBC, false);

    // apply Dirichlet conditions to (0,1)
    trkdm->ApplyDirichlet(dirichtoggle, false);

    // apply Dirichlet conditions to (0,0) block
    Teuchos::RCP<Epetra_Vector> zeros   = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));
    Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
    LINALG::ApplyDirichlettoSystem(stiffmt, sold, rhscopy, zeros, dirichtoggle);

    // row map (equals domain map) extractor
    LINALG::MapExtractor rowmapext(*mergedmap, MT_glmdofrowmap_, ProblemDofs());
    LINALG::MapExtractor dommapext(*mergedmap, MT_glmdofrowmap_, ProblemDofs());

    // build 2x2 block matrix for SIMPLER
    blockMat = Teuchos::rcp( new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext, rowmapext, 81, false, false));
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat = Teuchos::rcp_dynamic_cast< LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> >( blockMat);
    mat->Assign(0, 0, LINALG::View, *stiffmt);
    mat->Assign(0, 1, LINALG::View, *trkdm);
    mat->Assign(1, 0, LINALG::View, *trkmd);
    mat->Complete();

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(// force_residual_merged_export
        new Epetra_Vector(*mergedmap));
    LINALG::Export(*fd, *fresmexp);


    mergedrhs->Update(1.0, *fresmexp, 1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));

    if(constrrhs_ != Teuchos::null)
      LINALG::Export(*constrrhs_, *constrexp);

    err=mergedrhs->Update(1.0, *constrexp, 1.0);
    if(err!=0)
      dserror("failed");

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    LINALG::ApplyDirichlettoSystem(mergedsol, mergedrhs, mergedzeros,
        dirichtoggleexp);

    blocksol = mergedsol;
    blockrhs = mergedrhs;

    return;
  }


  if (stype_ == INPAR::CONTACT::solution_lagmult)
  {
    //**********************************************************************
    // FRICTIONAL CONTACT
    //**********************************************************************
    if (friction_)
      dserror("ERROR: Friction and Smoothing not implemented yet!");

    //**********************************************************************
    // FRICTIONLESS CONTACT
    //**********************************************************************

    // the standard stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiffmt = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kdd);

    // initialize merged system (matrix, rhs, sol)
    Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(ProblemDofs(),LMDoFRowMapPtr(true),false);
    mergedmap = LINALG::MergeMap(mergedmap,MT_glmdofrowmap_,false);
    Teuchos::RCP<LINALG::SparseMatrix> mergedmt = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> mergedrhs = LINALG::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedsol = LINALG::CreateVector(*mergedmap);
    Teuchos::RCP<Epetra_Vector> mergedzeros = LINALG::CreateVector(*mergedmap);

    // initialize constraint matrices for 3x3 system
    //      [ dd | dz | dm ]    [ dd | dz | dm ]  (displacement)
    //  k = [ zd | zz | zm ]  = [ zd | zz | 0  ]  (contact)
    //      [ md | mz | mm ]    [ md | 0  | 0  ]  (meshtying)

    Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
    Teuchos::RCP<LINALG::SparseMatrix> kdm = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));

    Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,true));
    Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));

    Teuchos::RCP<LINALG::SparseMatrix> kmd = Teuchos::rcp(new LINALG::SparseMatrix(*MT_gsdofrowmap_,100,false,true));

    // initialize transformed constraint matrices
    Teuchos::RCP<LINALG::SparseMatrix> trkdz,
    trkdm,trkzd,trkzz,trkzm,trkmd,trkmz,trkmm  ;

    //**********************************************************************
    // build matrix and vector blocks
    //**********************************************************************

    // build constraint matrix kdz
    kdz->Add(*dmatrix_, true, 1.0 - alphaf_, 1.0);
    kdz->Add(*mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdz->Complete(*gsdofrowmap_, *gdisprowmap_);

    // build constraint matrix kdm
    kdm->Add(*MT_dmatrix_, true, (1.0 - alphaf_), 1.0);
    kdm->Add(*MT_mmatrix_, true, -(1.0 - alphaf_), 1.0);
    kdm->Complete(*MT_gsdofrowmap_, *gdisprowmap_);

    // build constraint matrix kzd
    if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    {
      if (gactivedofs_->NumGlobalElements())
      {
        kzd->Add(*smatrix_, false, 1.0, 1.0);
        kzd->Add(*tderivmatrix_, false, 1.0, 1.0);
      }
    }
    else
    {
      if (gactiven_->NumGlobalElements())
        kzd->Add(*smatrix_, false, 1.0, 1.0);
      if (gactivet_->NumGlobalElements())
        kzd->Add(*tderivmatrix_, false, 1.0, 1.0);
    }

    kzd->Complete(*gdisprowmap_, *gsdofrowmap_);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,
        *gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(
        new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    if (gidofs->NumGlobalElements())
      kzz->Add(*onesdiag, false, 1.0, 1.0);
    if (gactivet_->NumGlobalElements())
      kzz->Add(*tmatrix_, false, 1.0, 1.0);
    kzz->Complete(*gsdofrowmap_, *gsdofrowmap_);

    // build constraint matrix kmd
    kmd->Add(*MT_dmatrix_, false, 1.0, 1.0);
    kmd->Add(*MT_mmatrix_, false, -1.0, 1.0);
    kmd->Complete(*gdisprowmap_, *MT_gsdofrowmap_);

    //**********************************************************************
    // transform matrix blocks to lm-dof-map
    //**********************************************************************
    //TODO Parallel Redistribution (siehe contact_lagrange_str)

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = MORTAR::MatrixRowTransformGIDs(kzd, LMDoFRowMapPtr(true));

    // transform constraint matrix kmd (MatrixRowTransform)
    trkmd = MORTAR::MatrixRowTransformGIDs(kmd, MT_glmdofrowmap_);

    // transform constraint matrix kdz (MatrixColTransform)
    trkdz = MORTAR::MatrixColTransformGIDs(kdz, LMDoFRowMapPtr(true));

    // transform constraint matrix kdm (MatrixColTransform)
    trkdm = MORTAR::MatrixColTransformGIDs(kdm, MT_glmdofrowmap_);

    // transform constraint matrix kzz (MatrixRowColTransform)
    trkzz = MORTAR::MatrixRowColTransformGIDs(kzz, LMDoFRowMapPtr(true), LMDoFRowMapPtr(true));

    //**********************************************************************
    // build and solve saddle point system
    //**********************************************************************
    if (systype == INPAR::CONTACT::system_saddlepoint)
    {
      // like for contact_lagrange_strategy.cpp

      Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      LINALG::Export(*dirichtoggle, *dirichtoggleexp);

      Teuchos::RCP<Epetra_Vector> lmDBC = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      LINALG::Export(*dirichtoggle, *lmDBC);

      lmDBC->ReplaceMap(LMDoFRowMap(true));
      Teuchos::RCP<Epetra_Vector> lmDBCexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      LINALG::Export(*lmDBC, *lmDBCexp);

      int err = dirichtoggleexp->Update(1., *lmDBCexp, 1.);
      if(err!=0)
        dserror("failed");

      // apply Dirichlet conditions to (1,0)
      trkzd->ApplyDirichlet(lmDBC, false);

      //----------
      Teuchos::RCP<Epetra_Vector> MT_lmDBC = Teuchos::rcp(new Epetra_Vector(*MT_gsdofrowmap_));
      LINALG::Export(*dirichtoggle, *MT_lmDBC);

      err=MT_lmDBC->ReplaceMap(*MT_glmdofrowmap_);
      if(err!=0)
        dserror("ReplaceMap failed");

      Teuchos::RCP<Epetra_Vector> MT_lmDBCexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      LINALG::Export(*MT_lmDBC, *MT_lmDBCexp);

      err = dirichtoggleexp->Update(1., *MT_lmDBCexp, 1.);
      if(err!=0)
        dserror("ERROR: Update procedure failed!");

      // apply Dirichlet conditions to (2,0)
      trkmd->ApplyDirichlet(MT_lmDBC, false);

      // apply Dirichlet conditions to (0,2)
      trkdm->ApplyDirichlet(dirichtoggle, false);

      //----------


      // apply Dirichlet conditions to (1,1) block
      trkzz->Complete();
      trkzz->ApplyDirichlet(lmDBC, true);

      // apply Dirichlet conditions to (0,0) block
      Teuchos::RCP<Epetra_Vector> zeros   = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));
      Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
      LINALG::ApplyDirichlettoSystem(stiffmt, sold, rhscopy, zeros, dirichtoggle);

      //apply Dirichlet conditions to (0,1) block
      trkdz->ApplyDirichlet(dirichtoggle, false);


      // set a helper flag for the CheapSIMPLE preconditioner (used to detect, if Teuchos::nullspace has to be set explicitely)
      // do we need this? if we set the Teuchos::nullspace when the solver is constructed?
      //solver.Params().set<bool>("CONTACT",true); // for simpler precond

      // reduce 3x3 system to 2x2 system
      //      [ trdd | trdg ]  (displacement)
      //  k = [ trgd | trgg ]  (contact and meshtying)

      Teuchos::RCP<Epetra_Map> MTco_mergedmap = LINALG::MergeMap(MT_glmdofrowmap_,LMDoFRowMapPtr(true),false);

      Teuchos::RCP<LINALG::SparseMatrix> trkdg = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
      Teuchos::RCP<LINALG::SparseMatrix> trkgd = Teuchos::rcp(new LINALG::SparseMatrix(*MTco_mergedmap,100,false,true));
      Teuchos::RCP<LINALG::SparseMatrix> trkgg = Teuchos::rcp(new LINALG::SparseMatrix(*MTco_mergedmap,100,false,true));

      trkdg->Add(*trkdz, false, 1.0, 1.0);
      trkdg->Add(*trkdm, false, 1.0, 1.0);
      trkgd->Add(*trkzd, false, 1.0, 1.0);
      trkgd->Add(*trkmd, false, 1.0, 1.0);
      trkgg->Add(*trkzz, false, 1.0, 1.0);

      // row map (equals domain map) extractor
      LINALG::MapExtractor rowmapext(*mergedmap, MTco_mergedmap, ProblemDofs());
      LINALG::MapExtractor dommapext(*mergedmap, MTco_mergedmap, ProblemDofs());

      // build 2x2 block matrix for SIMPLER
      blockMat = Teuchos::rcp( new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext, rowmapext, 81, false, false));
      Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat = Teuchos::rcp_dynamic_cast< LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> >( blockMat);
      mat->Assign(0, 0, LINALG::View, *stiffmt);
      mat->Assign(0, 1, LINALG::View, *trkdg);
      mat->Assign(1, 0, LINALG::View, *trkgd);
      mat->Assign(1, 1, LINALG::View, *trkgg);
      mat->Complete();

      // we also need merged rhs here
      Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(// force_residual_merged_export
          new Epetra_Vector(*mergedmap));
      LINALG::Export(*fd, *fresmexp);

      mergedrhs->Update(1.0, *fresmexp, 1.0);
      Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));

      if(constrrhs_ != Teuchos::null)
        LINALG::Export(*constrrhs_, *constrexp);

      err=mergedrhs->Update(1.0, *constrexp, 1.0);
      if(err!=0)
        dserror("ERROR: Update procedure failed!");

      // apply Dirichlet B.C. to mergedrhs and mergedsol
      LINALG::ApplyDirichlettoSystem(mergedsol, mergedrhs, mergedzeros,
          dirichtoggleexp);

      blocksol = mergedsol;
      blockrhs = mergedrhs;

      return;
    }

    //**********************************************************************
    // invalid system types
    //**********************************************************************
    else
    {
      dserror("ERROR: Invalid system type in SaddlePointSolve");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | set current deformation state                            farah 01/15 |
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::SetState(
    const enum MORTAR::StateType& statetype,
    const Epetra_Vector& vec)
{
  switch(statetype)
  {
    case MORTAR::state_new_displacement:
    case MORTAR::state_old_displacement:
    {
      // set state on contact interfaces
      for (std::size_t i = 0; i < interface_.size(); ++i)
        interface_[i]->SetState(statetype, vec);

      // set state on meshtying interfaces
      for (std::size_t i = 0; i < minterface_.size(); ++i)
        minterface_[i]->SetState(statetype, vec);
      break;
    }
    default:
    {
      dserror("Unsupported state type! (state type = %s)",
          MORTAR::StateType2String(statetype).c_str());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  do mortar coupling in reference configuration            farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::MortarCoupling(
    const Teuchos::RCP<Epetra_Vector> dis)
{
  //********************************************************************
  // initialize and evaluate interfaces
  //********************************************************************
  // for all interfaces
  for (int i = 0; i < (int) minterface_.size(); ++i)
  {
    // initialize / reset interfaces
    minterface_[i]->Initialize();

    // evaluate interfaces
    minterface_[i]->Evaluate();
  }

  //********************************************************************
  // restrict mortar treatment to actual meshtying zone
  //********************************************************************
  //  RestrictMeshtyingZone();

  //********************************************************************
  // initialize and evaluate global mortar stuff
  //********************************************************************
  // (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors
  MT_dmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*MT_gsdofrowmap_, 10));
  MT_mmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*MT_gsdofrowmap_, 100));
  MT_g_ = LINALG::CreateVector(*MT_gsdofrowmap_, true);

  // assemble D- and M-matrix on all interfaces
  for (int i = 0; i < (int) minterface_.size(); ++i)
    minterface_[i]->AssembleDM(*MT_dmatrix_, *MT_mmatrix_);

  // FillComplete() global Mortar matrices
  MT_dmatrix_->Complete();
  MT_mmatrix_->Complete(*MT_gmdofrowmap_, *MT_gsdofrowmap_);

  // compute g-vector at global level
  Teuchos::RCP<Epetra_Vector> xs = LINALG::CreateVector(*MT_gsdofrowmap_, true);
  Teuchos::RCP<Epetra_Vector> xm = LINALG::CreateVector(*MT_gmdofrowmap_, true);
  AssembleCoords("slave", true, xs);
  AssembleCoords("master", true, xm);
  Teuchos::RCP<Epetra_Vector> Dxs = Teuchos::rcp(
      new Epetra_Vector(*MT_gsdofrowmap_));
  MT_dmatrix_->Multiply(false, *xs, *Dxs);
  Teuchos::RCP<Epetra_Vector> Mxm = Teuchos::rcp(
      new Epetra_Vector(*MT_gsdofrowmap_));
  MT_mmatrix_->Multiply(false, *xm, *Mxm);
  MT_g_->Update(1.0, *Dxs, 1.0);
  MT_g_->Update(-1.0, *Mxm, 1.0);

  return;
}


/*----------------------------------------------------------------------*
 | create vector with coords                                 farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::AssembleCoords(const std::string& sidename,
    bool ref, Teuchos::RCP<Epetra_Vector> vec)
{
  // NOTE:
  // An alternative way of doing this would be to loop over
  // all interfaces and to assemble the coordinates there.
  // In thast case, one would have to be very careful with
  // edge nodes / crosspoints, which must not be assembled
  // twice. The solution would be to overwrite the corresp.
  // entries in the Epetra_Vector instead of using Assemble().

  // decide which side (slave or master)
  Teuchos::RCP<Epetra_Map> sidemap = Teuchos::null;
  if (sidename == "slave")
    sidemap = MT_gsnoderowmap_;
  else if (sidename == "master")
    sidemap = MT_gmnoderowmap_;
  else
    dserror("ERROR: Unknown sidename");

  // loop over all row nodes of this side (at the global level)
  for (int j = 0; j < sidemap->NumMyElements(); ++j) {
    int gid = sidemap->GID(j);

    // find this node in interface discretizations
    bool found = false;
    DRT::Node* node = NULL;
    for (int k = 0; k < (int) minterface_.size(); ++k) {
      found = minterface_[k]->Discret().HaveGlobalNode(gid);
      if (found) {
        node = minterface_[k]->Discret().gNode(gid);
        break;
      }
    }
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* mtnode = dynamic_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(Dim());
    std::vector<int> lm(Dim());
    std::vector<int> lmowner(Dim());

    for (int k = 0; k < Dim(); ++k) {
      // reference (true) or current (false) configuration
      if (ref)
        val[k] = mtnode->X()[k];
      else
        val[k] = mtnode->xspatial()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*vec, val, lm, lmowner);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Recovery method                                           farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::Recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
    return;

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<
      INPAR::MORTAR::ShapeFcn>(Params(), "LM_SHAPEFCN");
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::SystemType>(Params(), "SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL)
  //**********************************************************************
  //**********************************************************************
  if (systype == INPAR::CONTACT::system_condensed)
  {
    dserror("ERROR: Condensed system not implemented!");

    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual
        && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");
  }

  //**********************************************************************
  //**********************************************************************
  // CASE B: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // do nothing (z_ was part of solution already)
  }

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);

  return;
}


/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)     farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::UpdateActiveSetSemiSmooth()
{
  // get out gof here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all interfaces
    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      // loop over all slave nodes on the current interface
      for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements();++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node)
          dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // The nested active set strategy cannot deal with the case of
        // active nodes that have no integration segments/cells attached,
        // as this leads to zero rows in D and M and thus to singular systems.
        // However, this case might possibly happen when slave nodes slide
        // over the edge of a master body within one fixed active set step.
        // (Remark: Semi-smooth Newton has no problems in this case, as it
        // updates the active set after EACH Newton step, see below, and thus
        // would always set the corresponding nodes to INACTIVE.)
        if (cnode->Active() && !cnode->HasSegment())
          dserror("ERROR: Active node %i without any segment/cell attached",
              cnode->Id());
      }
    }
    return;
  }

  // get input parameter ftype
  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<
      INPAR::CONTACT::FrictionType>(Params(), "FRICTION");

  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weigh the two violations via cn!
  double cn_input = Params().get<double>("SEMI_SMOOTH_CN");
  double ct_input = Params().get<double>("SEMI_SMOOTH_CT");

  // this is the complementarity parameter we use for the decision.
  // it might be scaled with a mesh-size dependent factor
  double cn = 0.;
  double ct = 0.;

  // do we use a mesh-size scaling for cn and ct?
  bool adaptive_cn = DRT::INPUT::IntegralValue<int>(Params(),
      "MESH_ADAPTIVE_CN");
  bool adaptive_ct = DRT::INPUT::IntegralValue<int>(Params(),
      "MESH_ADAPTIVE_CT");

  // do we apply a nodal scaling for better matrix condition
  bool scale = DRT::INPUT::IntegralValue<int>(scontact_, "LM_NODAL_SCALE");

  // assume that active set has converged and check for opposite
  activesetconv_ = true;

  // loop over all interfaces
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j = 0; j < interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);

      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // get scaling factor
      double scalefac = 1.;
      if (scale == true && cnode->MoData().GetScale() != 0.)
        scalefac = cnode->MoData().GetScale();

      // compute weighted gap
      double wgap = cnode->CoData().Getg() / scalefac;

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k = 0; k < 3; ++k) {
        nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];
        nzold += cnode->MoData().n()[k] * cnode->MoData().lmold()[k];
      }

      // calculate mesh-size scaled version of cn
      cn = cn_input;
      ct = ct_input;
      if (adaptive_cn || adaptive_ct)
        if (cnode->MoData().GetD().size() != 0) // only do that if there is a D matrix
        {
          // row sum of D matrix
          double sumd = 0.;
          GEN::pairedvector<int, double>::const_iterator p;
          for (p = cnode->MoData().GetD().begin();
              p != cnode->MoData().GetD().end(); p++)
            sumd += p->second;
          double mesh_h = pow(sumd, 1. / ((double) Dim() - 1.));
          if (adaptive_cn && mesh_h != 0.)
            cn /= mesh_h;
          if (adaptive_ct && mesh_h != 0.)
            ct /= mesh_h;
        }

      // friction
      std::vector<double> tz(Dim() - 1, 0);
      std::vector<double> tjump(Dim() - 1, 0);
      double euclidean = 0.0;

      if (friction_) {
        // static cast
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute tangential parts and of Lagrange multiplier and incremental jumps
        for (int i = 0; i < Dim(); ++i) {
          tz[0] += frinode->CoData().txi()[i] * frinode->MoData().lm()[i];
          if (Dim() == 3)
            tz[1] += frinode->CoData().teta()[i] * frinode->MoData().lm()[i];

          if (DRT::INPUT::IntegralValue<int>(Params(), "GP_SLIP_INCR")
              == false) {
            tjump[0] += frinode->CoData().txi()[i]
                                                * frinode->FriData().jump()[i];
            if (Dim() == 3)
              tjump[1] += frinode->CoData().teta()[i]
                                                   * frinode->FriData().jump()[i];
          }
        }

        if (DRT::INPUT::IntegralValue<int>(Params(), "GP_SLIP_INCR") == true) {
          tjump[0] = frinode->FriData().jump_var()[0];
          if (Dim() == 3)
            tjump[1] = frinode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm |tz+ct.tjump|
        std::vector<double> sum(Dim() - 1, 0);
        sum[0] = tz[0] + ct * tjump[0];
        if (Dim() == 3)
          sum[1] = tz[1] + ct * tjump[1];
        if (Dim() == 2)
          euclidean = abs(sum[0]);
        if (Dim() == 3)
          euclidean = sqrt(sum[0] * sum[0] + sum[1] * sum[1]);
      }

      // adhesion
      double adhbound = 0.0;
      if (DRT::INPUT::IntegralValue<INPAR::CONTACT::AdhesionType>(Params(),
          "ADHESION") == INPAR::CONTACT::adhesion_bound)
        adhbound = interface_[i]->IParams().get<double>("ADHESION_BOUND");

      // check nodes of inactive set *************************************
      if (cnode->Active() == false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << std::endl;

        // check for penetration and/or tensile contact forces
        if (nz - cn * wgap > 0) // no averaging of Lagrange multipliers
          //if ((0.5*nz+0.5*nzold) - cn*wgap > 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = true;
          activesetconv_ = false;

          // friction
          if (friction_) {
            // nodes coming into contact
            dynamic_cast<FriNode*>(cnode)->FriData().Slip() = true;
          }
        }
      }

      // check nodes of active set ***************************************
      else
      {
        // check for fulfillment of contact condition
        //if (abs(wgap) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << std::endl;

        //adhesion modification
        nz += adhbound;

        // check for tensile contact forces and/or penetration
        if (nz - cn * wgap <= 0) // no averaging of Lagrange multipliers
          //if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = false;
          activesetconv_ = false;

          // friction
          if (friction_)
            dynamic_cast<FriNode*>(cnode)->FriData().Slip() = false;
        }

        // only do something for friction
        else
        {
          // friction tresca
          if (ftype == INPAR::CONTACT::friction_tresca) {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);

            // CAREFUL: friction bound is now interface-local (popp 08/2012)
            double frbound = interface_[i]->IParams().get<double>("FRBOUND");

            if (frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if (euclidean - frbound <= 0)
              {
              }
              // do nothing (stick was correct)
              else
              {
                frinode->FriData().Slip() = true;
                activesetconv_ = false;
              }
            }
            else
            {
              // check (euclidean)-frbound > 0
              if (euclidean - frbound > 0)
              {
              }
              // do nothing (slip was correct)
              else
              {
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (fytpe=="tresca")

          // friction coulomb
          if (ftype == INPAR::CONTACT::friction_coulomb)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);

            // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
            double frcoeff = interface_[i]->IParams().get<double>("FRCOEFF");

            if (frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if (euclidean - frcoeff * (nz - cn * wgap) <= 1e-10)
              {
              }
              // do nothing (stick was correct)
              else {
                frinode->FriData().Slip() = true;
                activesetconv_ = false;
              }
            } else {
              // check (euclidean)-frbound > 0
              if (euclidean - frcoeff * (nz - cn * wgap) > -1e-10)
              {
              }
              // do nothing (slip was correct)
              else {
                frinode->FriData().Slip() = false;
                activesetconv_ = false;
              }
            }
          } // if (ftype == INPAR::CONTACT::friction_coulomb)
        } // if (nz - cn*wgap <= 0)
      } // if (cnode->Active()==false)
    } // loop over all slave nodes
  } // loop over all interfaces

  // broadcast convergence status among processors
  int convcheck = 0;
  int localcheck = activesetconv_;
  Comm().SumAll(&localcheck, &convcheck, 1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck != Comm().NumProc()) {
    activesetconv_ = false;
    activesetsteps_ += 1;
  }

  // also update special flag for semi-smooth Newton convergence
  activesetssconv_ = activesetconv_;

  // update zig-zagging history (shift by one)
  if (zigzagtwo_ != Teuchos::null)
    zigzagthree_ = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_ != Teuchos::null)
    zigzagtwo_ = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_ != Teuchos::null)
    zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int) interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,
        interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(),
        false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(),
        false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(),
        false);

    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(),
          false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(),
          false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // CHECK FOR ZIG-ZAGGING / JAMMING OF THE ACTIVE SET
  // *********************************************************************
  // A problem of the active set strategy which sometimes arises is known
  // from optimization literature as jamming or zig-zagging. This means
  // that within a load/time-step the semi-smooth Newton algorithm can get
  // stuck between more than one intermediate solution due to the fact that
  // the active set decision is a discrete decision. Hence the semi-smooth
  // Newton algorithm fails to converge. The non-uniquenesss results either
  // from highly curved contact surfaces or from the FE discretization.
  // *********************************************************************
  // To overcome this problem we monitor the development of the active
  // set scheme in our contact algorithms. We can identify zig-zagging by
  // comparing the current active set with the active set of the second-
  // and third-last iteration. If an identity occurs, we interfere and
  // let the semi-smooth Newton algorithm restart from another active set
  // (e.g. intermediate set between the two problematic candidates), thus
  // leading to some kind of damped / modified semi-smooth Newton method.
  // This very simple approach helps stabilizing the contact algorithm!
  // *********************************************************************
  int zigzagging = 0;
  // FIXGIT: For friction zig-zagging is not eliminated
  if (ftype != INPAR::CONTACT::friction_tresca
      && ftype != INPAR::CONTACT::friction_coulomb) {
    // frictionless contact
    if (ActiveSetSteps() > 2) {
      if (zigzagtwo_ != Teuchos::null) {
        if (zigzagtwo_->SameAs(*gactivenodes_)) {
          // detect zig-zagging
          zigzagging = 1;
        }
      }

      if (zigzagthree_ != Teuchos::null) {
        if (zigzagthree_->SameAs(*gactivenodes_)) {
          // detect zig-zagging
          zigzagging = 2;
        }
      }
    }
  } // if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)

  // output to screen
  if (Comm().MyPID() == 0) {
    if (zigzagging == 1) {
      std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................."
          << std::endl;
    } else if (zigzagging == 2) {
      std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................"
          << std::endl;
    } else {
      // do nothing, no zig-zagging
    }
  }

  // reset zig-zagging history
  if (activesetconv_ == true) {
    zigzagone_ = Teuchos::null;
    zigzagtwo_ = Teuchos::null;
    zigzagthree_ = Teuchos::null;
  }

  // output of active set status to screen
  if (Comm().MyPID() == 0 && activesetconv_ == false)
    std::cout << "ACTIVE CONTACT SET HAS CHANGED... CHANGE No. "
    << ActiveSetSteps() - 1 << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements()) {
    isincontact_ = true;
    wasincontact_ = true;
  } else
    isincontact_ = false;

  return;
}


/*----------------------------------------------------------------------*
 |  UpdateDisplacementsAndLMincrements                       farah 01/15|
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::UpdateDisplacementsAndLMincrements(
    Teuchos::RCP<Epetra_Vector> sold,
    Teuchos::RCP<Epetra_Vector> blocksol)
{
  if (stype_ == INPAR::CONTACT::solution_penalty)
  {
    //**********************************************************************
    // extract results only for displacement and LM-meshtying increments
    //**********************************************************************

    //Teuchos::RCP<Epetra_Map> MTco_mergedmap = LINALG::MergeMap(MT_glmdofrowmap_,glmdofrowmap_,false);
    Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(ProblemDofs(),MT_glmdofrowmap_,false);

    Teuchos::RCP<Epetra_Vector> MT_sollm = Teuchos::rcp(new Epetra_Vector(*MT_glmdofrowmap_));

    LINALG::MapExtractor mapext(*mergedmap,ProblemDofs(),MT_glmdofrowmap_);
    mapext.ExtractCondVector(blocksol,sold);
    mapext.ExtractOtherVector(blocksol,MT_sollm);

    MT_sollm->ReplaceMap(*MT_gsdofrowmap_);

    if (IsSelfContact())
      // for self contact, slave and master sets may have changed,
      // thus we have to reinitialize the LM vector map
    {
      //      zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
      //      LINALG::Export(*z_, *zincr_);                     // change the map of z_
      //      z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
      //      zincr_->Update(1.0, *sollm, 0.0);                 // save sollm in zincr_
      //      z_->Update(1.0, *zincr_, 1.0);                    // update z_
      dserror("ERROR: Selfcontact not yet implemented for smoothing!!!");
    }
    else
    {

      MT_zincr_->Update(1.0, *MT_sollm, 0.0);
      MT_z_->Update(1.0, *MT_zincr_, 1.0);

    }
  }

  if (stype_ == INPAR::CONTACT::solution_lagmult)
  {
    //**********************************************************************
    // extract results for displacement and LM (from meshtying and contact) increments
    //**********************************************************************

    Teuchos::RCP<Epetra_Map> MTco_mergedmap = LINALG::MergeMap(MT_glmdofrowmap_,LMDoFRowMapPtr(true),false);
    Teuchos::RCP<Epetra_Map> mergedmap = LINALG::MergeMap(ProblemDofs(),LMDoFRowMapPtr(true),false);
    mergedmap = LINALG::MergeMap(mergedmap,MT_glmdofrowmap_,false);

    Teuchos::RCP<Epetra_Vector> MTco_sollm = Teuchos::rcp(new Epetra_Vector(*MTco_mergedmap));
    Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*LMDoFRowMapPtr(true)));
    Teuchos::RCP<Epetra_Vector> MT_sollm = Teuchos::rcp(new Epetra_Vector(*MT_glmdofrowmap_));

    LINALG::MapExtractor mapext(*mergedmap,ProblemDofs(),MTco_mergedmap);
    mapext.ExtractCondVector(blocksol,sold);
    mapext.ExtractOtherVector(blocksol,MTco_sollm);

    //**********************************************************************
    // extract results for Contact_LM and Meshtying_lm increments
    //**********************************************************************

    LINALG::MapExtractor MTco_mapext(*MTco_mergedmap,LMDoFRowMapPtr(true),MT_glmdofrowmap_);
    MTco_mapext.ExtractCondVector(MTco_sollm,sollm);
    MTco_mapext.ExtractOtherVector(MTco_sollm,MT_sollm);

    sollm->ReplaceMap(*gsdofrowmap_);
    MT_sollm->ReplaceMap(*MT_gsdofrowmap_);

    if (IsSelfContact())
      // for self contact, slave and master sets may have changed,
      // thus we have to reinitialize the LM vector map
    {
      zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
      LINALG::Export(*z_, *zincr_);                     // change the map of z_
      z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
      zincr_->Update(1.0, *sollm, 0.0);                 // save sollm in zincr_
      z_->Update(1.0, *zincr_, 1.0);                    // update z_
    }
    else
    {
      zincr_->Update(1.0, *sollm, 0.0);
      z_->Update(1.0, *zincr_, 1.0);

      MT_zincr_->Update(1.0, *MT_sollm, 0.0);
      MT_z_->Update(1.0, *MT_zincr_, 1.0);

    }
  }
  //-------------

  return;
}


/*----------------------------------------------------------------------*
 | EvalConstrRHS                                            farah 01/15 |
 *----------------------------------------------------------------------*/
void CONTACT::SmoothingStrategy::EvalConstrRHS()
{
  if (stype_ == INPAR::CONTACT::solution_penalty)
    return;

  // get system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");
  if (systype == INPAR::CONTACT::system_condensed)
    return;

  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    constrrhs_          = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) * z^(i)_(n+1)

  // export weighted gap vector
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_,*gact);
    }
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gactiven_->NumGlobalElements())
    {
      LINALG::Export(*g_,*gact);
      gact->ReplaceMap(*gactiven_);
    }
  }

  Teuchos::RCP<Epetra_Vector> gact_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*gact,*gact_exp);

  constrrhs->Update(-1.0,*gact_exp,1.0);

  // export inactive rhs
  Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*inactiverhs_, *inactiverhsexp);

  // build constraint rhs (1)
  constrrhs->Update(1.0, *inactiverhsexp, 1.0);

  // *** CASE 1: FRICTIONLESS CONTACT *******************************************************
  if (!friction_)
  {
    // export tangential rhs
    Teuchos::RCP<Epetra_Vector> tangrhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*tangrhs_, *tangrhs_exp);

    // build constraint rhs (2)
    constrrhs->Update(1.0, *tangrhs_exp, 1.0);

  }
  // *** CASE 2: FRICTIONAL CONTACT *******************************************************
  else
  {
    dserror("FRICTION!!!");
  }

  constrrhs->ReplaceMap(LMDoFRowMap(true));

  constrrhs_ = constrrhs;                                 // set constraint rhs vector
  return;
}

