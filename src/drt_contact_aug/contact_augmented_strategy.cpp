/*!----------------------------------------------------------------------
\file contact_augmented_strategy.cpp

<pre>
Created on: Apr 7, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_lagrange_strategy.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                        hiermeier 04/14 |
 *----------------------------------------------------------------------*/
CONTACT::AugmentedLagrangeStrategy::AugmentedLagrangeStrategy(DRT::Discretization& probdiscret,
                                                              Teuchos::ParameterList params,
                                                              std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces,
                                                              int dim,
                                                              Teuchos::RCP<Epetra_Comm> comm,
                                                              double alphaf,
                                                              int maxdof) :
CoLagrangeStrategy(probdiscret,params,interfaces,dim,comm,alphaf,maxdof),
gFdCheck_(false)
{
  // cast to augmented interfaces
  for (int i=0; i<(int) interfaces.size(); ++i)
  {
    interface_.push_back(Teuchos::rcp_dynamic_cast<CONTACT::AugmentedInterface>(interfaces[i]));
    if (interface_[i]==Teuchos::null)
      dserror("AugmentedLagrangeStartegy: Interface-cast failed!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for              hiermeier 04/14 |
 | next Newton step                                                     |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::Initialize()
{
  // *** (re)setup global matrices ***
  dGLmSlLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  dGLmMaLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  dGGSlLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  dGGMaLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  augDnMatrix_  = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveNDofs_,100));
  augMnMatrix_  = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveNDofs_,100));

  aWGapLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveNDofs_,100));

  dLmNWGapLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveNDofs_,100));
  dLmTLmTMatrix_     = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveTDofs_,100));
  dLmTLmTLinMatrix_  = Teuchos::rcp(new LINALG::SparseMatrix(*gAugActiveSlaveTDofs_,100));

  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = LINALG::SplitMap(*gsdofrowmap_, *gAugActiveSlaveDofs_);

  augInactiveMatrix_    = Teuchos::rcp(new LINALG::SparseMatrix(*gAugInactiveSlaveDofs,100));
  augInactiveLinMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gAugInactiveSlaveDofs,100));

  // *** (re)setup global augmented Epetra_Vectors ***
  augLm_          = LINALG::CreateVector(*gAugActiveSlaveNDofs_, true);
  aWGapRhs_       = LINALG::CreateVector(*gAugActiveSlaveNDofs_, true);
  dLmTLmTRhs_     = LINALG::CreateVector(*gAugActiveSlaveTDofs_,true);
  augInactiveRhs_ = LINALG::CreateVector(*gAugInactiveSlaveDofs, true);

  dLmNWGapRhs_ = LINALG::CreateVector(*gAugActiveSlaveNDofs_, true);

  augfs_lm_ = LINALG::CreateVector(*gsdofrowmap_,true);
  augfm_lm_ = LINALG::CreateVector(*gmdofrowmap_,true);

  augfs_g_  = LINALG::CreateVector(*gsdofrowmap_,true);
  augfm_g_  = LINALG::CreateVector(*gmdofrowmap_,true);

  if (friction_)
    dserror("AugmentedLagrangeStrategy::Initialize: Frictional case is not yet considered!");

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate contact contributions (frictionless)        hiermeier 04/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvaluateContact(
                                Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                Teuchos::RCP<Epetra_Vector>& feff)
{
  // Call standard version
//  if (!gFdCheck_)
//    CoLagrangeStrategy::EvaluateContact(kteff,feff);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  /**********************************************************************/
  /*             UPDATE THE CONTACT RIGHT HAND SIDE                     */
  /**********************************************************************/
  /*---------------------------------------------------------------------*
   | calculate                                                           |
   *---------------------------------------------------------------------*
   | Build Dn and Mn matrices                                            |
   | and inactive right-hand side with old lagrange multipliers (incr)   |
   | and tangential right-hand side (incr)                               |
   *---------------------------------------------------------------------*/
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // *** augmented Lagrange formulation ********************************
    // *** Force Balance *************************************************
    interface_[i]->AssembleAugDnMnMatrix(*augDnMatrix_,*augMnMatrix_);
    interface_[i]->AssembleAugLmVector(*augLm_);
    interface_[i]->AssembleAWGapRhs(*aWGapRhs_);
    // *** CONSTRAINTS ***************************************************
    // active - normal direction
    interface_[i]->AssembleDLmNWGapRhs(*dLmNWGapRhs_);
    // active - tangential direction
    interface_[i]->AssembleDLmTLmTRhs(*dLmTLmTRhs_);
    // inactive - all directions
    interface_[i]->AssembleAugInactiveRhs(*augInactiveRhs_);
  }

  augDnMatrix_->Complete(*gsdofrowmap_,*gAugActiveSlaveNDofs_);
  augMnMatrix_->Complete(*gmdofrowmap_,*gAugActiveSlaveNDofs_);

  // If we are only interested in the RHS-update, we do not consider any
  // stiffness matrix entries and just return.
  if (iterls_ > -1) return;

  /**********************************************************************/
  /*             UPDATE THE TANGENTIAL STIFFNESS MATRIX                 */
  /**********************************************************************/
  // *** augmented Lagrange formulation *********************************
  for (int i=0; i<(int) interface_.size();++i)
  {
    // *** Force Balance ***************************************************
    // linearization w.r.t. displ.
    interface_[i]->AssembleDGLmLinMatrix(*dGLmSlLinMatrix_,
                                         *dGLmMaLinMatrix_);
    interface_[i]->AssembleDGGLinMatrix(*dGGSlLinMatrix_,
                                        *dGGMaLinMatrix_);
    // *** Constraints *****************************************************
    // linearization w.r.t. LM
//    interface_[i]->AssembleAugTMatrix(*augTMatrix_);
    interface_[i]->AssembleDLmTLmTMatrix(*dLmTLmTMatrix_);
    interface_[i]->AssembleAugInactiveMatrix(*augInactiveMatrix_);
    // linearization w.r.t. displ.
    // active - normal direction
    interface_[i]->AssembleDLmNWGapLinMatrix(*dLmNWGapLinMatrix_);
    // active - tangential direction
//    interface_[i]->AssembleLmTLinMatrix(*lmTLinMatrix_);
    interface_[i]->AssembleDLmTLmTLinMatrix(*dLmTLmTLinMatrix_);
    interface_[i]->AssembleAugInactiveLinMatrix(*augInactiveLinMatrix_);
  }

  // >>> START - FillComplete matrices <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // domainmap: columnmap | rangemap: rowmap
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = LINALG::SplitMap(*gsdofrowmap_, *gAugActiveSlaveDofs_);

  // *** Force Balance ***************************************************
  // linearization w.r.t. displ.
  dGLmSlLinMatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  dGLmMaLinMatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);
  dGGSlLinMatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  dGGMaLinMatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);
  // *** Constraints *****************************************************
  // linearization w.r.t. LM
//  augTMatrix_->Complete(*gsdofrowmap_,*gAugActiveSlaveTDofs_);
  dLmTLmTMatrix_->Complete(*gAugActiveSlaveTDofs_,*gAugActiveSlaveTDofs_);
  augInactiveMatrix_->Complete(*gAugInactiveSlaveDofs,*gAugInactiveSlaveDofs);
  // linearization w.r.t. displ.
  dLmNWGapLinMatrix_->Complete(*gsmdofrowmap_,*gAugActiveSlaveNDofs_);
//  lmTLinMatrix_->Complete(*gsdofrowmap_,*gAugActiveSlaveTDofs_);
  dLmTLmTLinMatrix_->Complete(*gsdofrowmap_,*gAugActiveSlaveTDofs_);
  augInactiveLinMatrix_->Complete(*gsdofrowmap_,*gAugInactiveSlaveDofs);
  // >>> END - FillComplete matrices <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // *** DEBUGGING ***
  // Finite Difference check at Gauss-point level
  AugFDCheckGP();

  // Finite Difference check at global level
  AugFDCheckGlobal();

  return;
}

/*----------------------------------------------------------------------*
 | Update and build local and global active set         hiermeier 05/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::BuildGlobalAugActiveSet(const int it)
{
  // Do the active set decision
  activesetconv_ = interface_[it]->UpdateAugActiveSetSemiSmooth();

  // also update special flag for semi-smooth Newton convergence
  activesetssconv_ = activesetconv_;

  // if we are at the first contact interface, we start with (re)setup of the
  // global Epetra_maps.
  if (it==0)
  {
    gAugActiveSlaveNodes_ = Teuchos::null;
    gAugActiveSlaveDofs_  = Teuchos::null;
    gAugActiveSlaveNDofs_ = Teuchos::null;
    gAugActiveSlaveTDofs_ = Teuchos::null;
  }

  // Update Active set
  gAugActiveSlaveNodes_ = LINALG::MergeMap(gAugActiveSlaveNodes_,interface_[it]->AugActiveSlaveNodes(),false);
  gAugActiveSlaveDofs_  = LINALG::MergeMap(gAugActiveSlaveDofs_,interface_[it]->AugActiveSlaveDofs(),false);
  gAugActiveSlaveNDofs_ = LINALG::MergeMap(gAugActiveSlaveNDofs_,interface_[it]->AugActiveSlaveNDofs(),false);
  gAugActiveSlaveTDofs_ = LINALG::MergeMap(gAugActiveSlaveTDofs_,interface_[it]->AugActiveSlaveTDofs(),false);

  // update flag for global contact status
  if (gAugActiveSlaveNodes_->NumGlobalElements())
  {
    isincontact_=true;
    wasincontact_=true;
  }
  else
    isincontact_=false;


  return;
}

/*----------------------------------------------------------------------*
 | Update structural right hand side                    hiermeier 05/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::UpdateStructuralRHS(Teuchos::RCP<Epetra_Vector>& feff)
{
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // For self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (IsSelfContact())
    dserror("ERROR: Augmented Lagrange Formulation: Self contact is not yet considered!");
  // if there is no self contact everything is ok
  else
  {
    // add contact force terms
    // *** Slave side ***
    Teuchos::RCP<Epetra_Vector> augfs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    augDnMatrix_->Multiply(true,*augLm_,*augfs);
    Teuchos::RCP<Epetra_Vector> augfs_exp = Teuchos::rcp(new Epetra_Vector(feff->Map()));
    LINALG::Export(*augfs,*augfs_exp);
    feff->Update((1.0-alphaf_),*augfs_exp,1.0);

    // Master side
    Teuchos::RCP<Epetra_Vector> augfm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    augMnMatrix_->Multiply(true,*augLm_,*augfm);
    Teuchos::RCP<Epetra_Vector> augfm_exp = Teuchos::rcp(new Epetra_Vector(feff->Map()));
    LINALG::Export(*augfm,*augfm_exp);
    feff->Update((1.0-alphaf_),*augfm_exp,1.0);

    // Check linear and angular momentum conservation
#ifdef CHECKCONSERVATIONLAWS
    CheckConservationLaws(*augfs,*augfm);
#endif
    if (alphaf_!=0.0) dserror("ERROR: No contact dynamics! Fco_old is still missing!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Check conservation of linear and angular momentum    hiermeier 05/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::CheckConservationLaws(const Epetra_Vector& augfs,
                                                               const Epetra_Vector& augfm)
{
  /****************** SLAVE SIDE ********************************************************/
  // standard Lagrange multiplier fraction
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Epetra_Vector> z_exp = Teuchos::rcp(new Epetra_Vector(*gAugActiveSlaveNDofs_));
  LINALG::Export(*z_,*z_exp);
  augDnMatrix_->Multiply(true,*z_exp,*augfs_lm);
  // regularization fraction
  Teuchos::RCP<Epetra_Vector> augfs_g = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  augDnMatrix_->Multiply(true,*aWGapRhs_,*augfs_g);
  // *** START: DEBUGGING OUTPUT  ***
//  // COMPARE two different assembly strategies
//  double cn= Params().get<double>("SEMI_SMOOTH_CN");
//  Teuchos::RCP<Epetra_Vector> augfs_check = Teuchos::rcp(new Epetra_Vector(*augfs_lm));
//  augfs_check->Update(-cn,*augfs_g,1.0);
//
//  std::cout << ">>> augfs <<<" << std::endl;
//  std::cout << augfs << std::endl;
//  std::cout << "======================================" << std::endl;
//  std::cout << ">>> augfs_lm <<<" << std::endl;
//  std::cout << *augfs_lm << std::endl;
//  std::cout << ">>> augfs_g <<<" << std::endl;
//  std::cout << *augfs_g << std::endl;
//  std::cout << "======================================" << std::endl;
//  std::cout << ">>> augfs_check<<<" << std::endl;
//  std::cout << *augfs_check << std::endl;
  // *** END: DEBUGGING OUTPUT  ***
  /****************** MASTER SIDE *******************************************************/
  // standard lagrange multiplier fraction
  Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  augMnMatrix_->Multiply(true,*z_exp,*augfm_lm);
  // regularization fraction
  Teuchos::RCP<Epetra_Vector> augfm_g = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  augMnMatrix_->Multiply(true,*aWGapRhs_,*augfm_g);
  // *** START: DEBUGGING OUTPUT  ***
//  // COMPARE two different assembly strategies
//  Teuchos::RCP<Epetra_Vector> augfm_check = Teuchos::rcp(new Epetra_Vector(*augfm_lm));
//  augfm_check->Update(-cn,*augfm_g,1.0);
//
//    std::cout << ">>> augfm <<<" << std::endl;
//    std::cout << augfm << std::endl;
//    std::cout << "======================================" << std::endl;
//    std::cout << ">>> augfm_lm <<<" << std::endl;
//    std::cout << *augfm_lm << std::endl;
//    std::cout << ">>> augfm_g <<<" << std::endl;
//    std::cout << *augfm_g << std::endl;
//    std::cout << "======================================" << std::endl;
//    std::cout << ">>> augfm_check<<<" << std::endl;
//    std::cout << *augfm_check << std::endl;
  // *** END: DEBUGGING OUTPUT  ***
  /*-------------------------------*
   | LINEAR MOMENTUM CONSERVATION  |
   *-------------------------------*/
  double ssum = 0.0;
  double msum = 0.0;
  double csum = 0.0;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Linear Momentum Conservation      <<" << std::endl;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Standard terms (lm)               <<" << std::endl;
  for (int i=0;i<augfs_lm->GlobalLength();++i) ssum+=(*augfs_lm)[i];
  std::cout << "SLAVE:   " << std::setw(14) << ssum<< std::endl;
  for (int i=0;i<augfm_lm->GlobalLength();++i) msum+=(*augfm_lm)[i];
  std::cout << "MASTER:  " << std::setw(14) << msum << std::endl;
  csum = ssum+msum;
  if (abs(csum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
  std::cout << "Balance: " << std::setw(14) << csum << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
  ssum=0.0;
  for (int i=0;i<augfs_g->GlobalLength();++i) ssum+=(*augfs_g)[i];
  std::cout << "SLAVE:   " << std::setw(14) << ssum<< std::endl;
  msum=0.0;
  for (int i=0;i<augfm_g->GlobalLength();++i) msum+=(*augfm_g)[i];
  std::cout << "MASTER:  " << std::setw(14) << msum << std::endl;
  csum = ssum+msum;
  if (abs(csum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
  std::cout << "Balance: " << std::setw(14) << csum << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << ">>      Complete                          <<" << std::endl;
  ssum=0.0;
  for (int i=0;i<augfs.GlobalLength();++i) ssum+= augfs[i];
  std::cout << "SLAVE:   " << std::setw(14) << ssum<< std::endl;
  msum=0.0;
  for (int i=0;i<augfm.GlobalLength();++i) msum+= augfm[i];
  std::cout << "MASTER:  " << std::setw(14) << msum << std::endl;
  csum = ssum+msum;
  if (abs(csum)>1.0e-11) dserror("Conservation of linear momentum is not fulfilled!");
  std::cout << "Balance: " << std::setw(14) << csum << std::endl;

  /*-------------------------------*
   | ANGULAR MOMENTUM CONSERVATION |
   *-------------------------------*/
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  std::cout << ">>      Angular Momentum Conservation     <<" << std::endl;
  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    std::cout << ">>----- Interface " << std::setw(2) << i;
    std::cout << " ---------------------<<" << std::endl;
    std::cout << ">>      Standard terms (lm)               <<" << std::endl;
    interface_[i]->EvalResultantMoment(*augfs_lm,*augfm_lm);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->EvalResultantMoment(*augfs_g,*augfm_g);

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << ">>      Complete                          <<" << std::endl;
    for (int i=0; i<(int)interface_.size(); ++i)
      interface_[i]->EvalResultantMoment(augfs,augfm);
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 | Evaluate the augmented contact forces (Sl & Ma)      hiermeier 06/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::AugForces(Epetra_Vector& gAugFs_lm, Epetra_Vector& gAugFs_g,
                                                   Epetra_Vector& gAugFm_lm, Epetra_Vector& gAugFm_g)
{
  if (!IsInContact()) return;

  double cn= Params().get<double>("SEMI_SMOOTH_CN");
  /****************** SLAVE SIDE ********************************************************/
  // *** standard Lagrange multiplier fraction ***
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Epetra_Vector> z_exp = Teuchos::rcp(new Epetra_Vector(*gAugActiveSlaveNDofs_));
  LINALG::Export(*z_,*z_exp);
  augDnMatrix_->Multiply(true,*z_exp,*augfs_lm);

  // Export
  LINALG::Export(*augfs_lm,gAugFs_lm);

  // *** regularization fraction ***
  Teuchos::RCP<Epetra_Vector> augfs_g = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  Teuchos::RCP<Epetra_Vector> aWGap = Teuchos::rcp(new Epetra_Vector(*aWGapRhs_));
  aWGap->Scale(-cn);
  augDnMatrix_->Multiply(true,*aWGap,*augfs_g);
  // Export
  LINALG::Export(*augfs_g,gAugFs_g);

  /****************** MASTER SIDE *******************************************************/
  // *** standard lagrange multiplier fraction ***
  Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  augMnMatrix_->Multiply(true,*z_exp,*augfm_lm);
  // Export
  LINALG::Export(*augfm_lm,gAugFm_lm);
  // *** regularization fraction ***
  Teuchos::RCP<Epetra_Vector> augfm_g = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  augMnMatrix_->Multiply(true,*aWGap,*augfm_g);
  // Export
  LINALG::Export(*augfm_g,gAugFm_g);

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses       hiermeier 07/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::OutputStresses()
{
  // reset contact stress class variables
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

      // get nodal normal and tangential directions
      double* nn  = cnode->MoData().n();
      double* nt1 = cnode->CoData().txi();
      double* nt2 = cnode->CoData().teta();
      double lmn  = cnode->MoData().lm()[0];
      double lmt1 = cnode->MoData().lm()[1];
      double lmt2 = cnode->MoData().lm()[2];

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(),"LM_NODAL_SCALE")==false
            || cnode->MoData().GetScale()==0.)
          (*stressnormal_)[locindex[dof]] = -lmn*nn[dof];
        else
          (*stressnormal_)[locindex[dof]] = -lmn*nn[dof]/cnode->MoData().GetScale();
      }

      // tangential stress components
      for (int dof=0;dof<dim;++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->Dofs()[dof]);
        if (DRT::INPUT::IntegralValue<int>(Params(),"LM_NODAL_SCALE")==false
            || cnode->MoData().GetScale()==0.)
          (*stresstangential_)[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof];
        else
          (*stresstangential_)[locindex[dof]] = -lmt1*nt1[dof]-lmt2*nt2[dof]/cnode->MoData().GetScale();
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Update structural stiffness matrix                   hiermeier 06/14 |
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::UpdateStructuralStiff(Teuchos::RCP<LINALG::SparseOperator>& kteff)
{
  // if there was no contact evaluation
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;
  // If we are only interested in the RHS-update, we do not consider any
  // stiffness matrix entries and just return.
  if (iterls_ > -1) return;

  // transform if necessary
  if (ParRedist())
  {
    dGLmSlLinMatrix_ = MORTAR::MatrixRowTransform(dGLmSlLinMatrix_,pgsdofrowmap_);
    dGLmMaLinMatrix_ = MORTAR::MatrixRowTransform(dGLmMaLinMatrix_,pgmdofrowmap_);
    dGGSlLinMatrix_  = MORTAR::MatrixRowTransform(dGGSlLinMatrix_,pgsdofrowmap_);
    dGGMaLinMatrix_ = MORTAR::MatrixRowTransform(dGGMaLinMatrix_,pgmdofrowmap_);
  }

  // build structural matrix kdd
  // NOTE: This step is moved from the EvaluateContact function here to completely separate
  kteff->UnComplete();
  kteff->Add(*dGLmSlLinMatrix_,false,-(1.0-alphaf_),1.0);
  kteff->Add(*dGLmMaLinMatrix_,false,-(1.0-alphaf_),1.0);
  kteff->Add(*dGGSlLinMatrix_,false,-(1.0-alphaf_),1.0);
  kteff->Add(*dGGMaLinMatrix_,false,-(1.0-alphaf_),1.0);
  kteff->Complete();

  return;
}

/*----------------------------------------------------------------------*
 | Calculate augmented constraint RHS entries            hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::EvalConstrRHS()
{
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
  {
    // (re)setup the vector
    augConstrRhs_          = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> augConstrRhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!

  // ToDo Three export vectors seem unnecessary and can be replaced by scaling a general export vector
  // with zero, between the different steps! Check for possible performance gain!
  // Add active constraints in normal direction:
  Teuchos::RCP<Epetra_Vector> dLmNWGapRhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dLmNWGapRhs_,*dLmNWGapRhs_exp);
  augConstrRhs->Update(-1.0,*dLmNWGapRhs_exp,1.0);

  // Add inactive constraints
  Teuchos::RCP<Epetra_Vector> augInactiveRhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*augInactiveRhs_,*augInactiveRhs_exp);
  augConstrRhs->Update(-1.0,*augInactiveRhs_exp,1.0);

  // Add tangential frictionless constraints
  Teuchos::RCP<Epetra_Vector> dLmTLmTRhs_exp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  LINALG::Export(*dLmTLmTRhs_,*dLmTLmTRhs_exp);
  augConstrRhs->Update(-1.0,*dLmTLmTRhs_exp,1.0);

  // replace row map
  augConstrRhs->ReplaceMap(*glmdofrowmap_);

  // set constraint rhs vector
  augConstrRhs_ = augConstrRhs;

  return;
}

/*----------------------------------------------------------------------*
 | Solve linear system of saddle point type              hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedLagrangeStrategy::SaddlePointSolve(LINALG::Solver& solver,
    LINALG::Solver& fallbacksolver,
    Teuchos::RCP<LINALG::SparseOperator> kdd,  Teuchos::RCP<Epetra_Vector> fd,
    Teuchos::RCP<Epetra_Vector>  sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
    int numiter)
{
  // create old style dirichtoggle vector (supposed to go away)
  // the use of a toggle vector is more flexible here. It allows to apply dirichlet
  // conditions on different matrix blocks separately.
  Teuchos::RCP<Epetra_Vector> dirichtoggle = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->FullMap())));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp,*dirichtoggle);

  // get system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");
  if (systype!=INPAR::CONTACT::system_saddlepoint)
    dserror("ERROR: Invalid system type in SaddlePointSolve");

  // check if contact contributions are present,
  // if not we make a standard solver call to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep())
  {
    //std::cout << "##################################################" << std::endl;
    //std::cout << " USE FALLBACK SOLVER (pure structure problem)" << std::endl;
    //std::cout << fallbacksolver.Params() << std::endl;
    //std::cout << "##################################################" << std::endl;

    // standard solver call
    fallbacksolver.Solve(kdd->EpetraOperator(),sold,fd,true,numiter==0);
    return;
  }

  //**********************************************************************
  // prepare saddle point system
  //**********************************************************************
  // the standard stiffness matrix
  Teuchos::RCP<LINALG::SparseMatrix> stiffmt = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kdd);

  // initialize merged system (matrix, rhs, sol)
  Teuchos::RCP<Epetra_Map>           mergedmap   = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
  Teuchos::RCP<Epetra_Vector>        mergedrhs   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedsol   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedzeros = LINALG::CreateVector(*mergedmap);

  // initialize transformed constraint matrices
  Teuchos::RCP<LINALG::SparseMatrix> trkdz, trkzd, trkzz;

  //**********************************************************************
  // build matrix blocks
  //**********************************************************************
  // *** CASE 1: FRICTIONLESS CONTACT ************************************
  if (friction_) dserror("Friction is not yet supported for the augmented contact formulation!");

  // build constraint matrix kdz
  Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
  kdz->Add(*augDnMatrix_,true,-(1.0-alphaf_),1.0);
  kdz->Add(*augMnMatrix_,true,-(1.0-alphaf_),1.0);
  kdz->Complete(*gsdofrowmap_,*gdisprowmap_);

  // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
  trkdz = MORTAR::MatrixColTransformGIDs(kdz,glmdofrowmap_);

  // transform parallel row distribution of constraint matrix kdz
  // (only necessary in the parallel redistribution case)
  if (ParRedist()) trkdz = MORTAR::MatrixRowTransform(trkdz,ProblemDofs());

  // build constraint matrix kzd
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = LINALG::SplitMap(*gsdofrowmap_,*gAugActiveSlaveDofs_);
  Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
  kzd->Add(*dLmNWGapLinMatrix_,false,1.0,1.0);
  kzd->Add(*dLmTLmTLinMatrix_,false,1.0,1.0);
  kzd->Add(*augInactiveLinMatrix_,false,1.0,1.0);
  kzd->Complete(*gdisprowmap_,*gsdofrowmap_);

  // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
  trkzd = MORTAR::MatrixRowTransformGIDs(kzd,glmdofrowmap_);

  // transform parallel column distribution of constraint matrix kzd
  // (only necessary in the parallel redistribution case)
  if (ParRedist()) trkzd = MORTAR::MatrixColTransform(trkzd,ProblemDofs());

  // build constraint matrix kzz
  Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
  kzz->Add(*augInactiveMatrix_,false,1.0,1.0);
  kzz->Add(*dLmTLmTMatrix_,false,1.0,1.0);
  kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);

  // transform constraint matrix kzz to lmdofmap (MatrixRowColTransform)
  trkzz = MORTAR::MatrixRowColTransformGIDs(kzz,glmdofrowmap_,glmdofrowmap_);

  //**********************************************************************
  // build and solve saddle point system
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  LINALG::Export(*dirichtoggle,*dirichtoggleexp);

  // apply Dirichlet conditions to (0,0) and (0,1) blocks
  Teuchos::RCP<Epetra_Vector> zeros   = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));
  Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
  LINALG::ApplyDirichlettoSystem(stiffmt,sold,rhscopy,zeros,dirichtoggle);
  trkdz->ApplyDirichlet(dirichtoggle,false);

  // row map (equals domain map) extractor
  LINALG::MapExtractor rowmapext(*mergedmap,glmdofrowmap_,ProblemDofs());
  LINALG::MapExtractor dommapext(*mergedmap,glmdofrowmap_,ProblemDofs());

  // set a helper flag for the CheapSIMPLE preconditioner (used to detect, if Teuchos::nullspace has to be set explicitely)
  // do we need this? if we set the Teuchos::nullspace when the solver is constructed?
  solver.Params().set<bool>("CONTACT",true); // for simpler precond

  // build block matrix for SIMPLER
  Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));
  mat->Assign(0,0,View,*stiffmt);
  mat->Assign(0,1,View,*trkdz);
  mat->Assign(1,0,View,*trkzd);
  mat->Assign(1,1,View,*trkzz);
  mat->Complete();

  // we also need merged rhs here
  Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  LINALG::Export(*fd,*fresmexp);
  mergedrhs->Update(1.0,*fresmexp,1.0);
  Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
  LINALG::Export(*augConstrRhs_,*constrexp);
  mergedrhs->Update(1.0,*constrexp,1.0);

  // apply Dirichlet B.C. to mergedrhs and mergedsol
  LINALG::ApplyDirichlettoSystem(mergedsol,mergedrhs,mergedzeros,dirichtoggleexp);

  // SIMPLER preconditioning solver call
  solver.Solve(mat->EpetraOperator(),mergedsol,mergedrhs,true,numiter==0);

  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  LINALG::MapExtractor mapext(*mergedmap,ProblemDofs(),glmdofrowmap_);
  mapext.ExtractCondVector(mergedsol,sold);
  mapext.ExtractOtherVector(mergedsol,sollm);
  sollm->ReplaceMap(*gsdofrowmap_);

  if (IsSelfContact())
  // for self contact, slave and master sets may have changed,
  // thus we have to reinitialize the LM vector map
  {
    zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
    LINALG::Export(*z_, *zincr_);                   // change the map of z_
    z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
    zincr_->Update(1.0, *sollm, 0.0);               // save sollm in zincr_
    z_->Update(1.0, *zincr_, 1.0);                  // update z_

  }
  else
  {
    zincr_->Update(1.0, *sollm, 0.0);
    z_->Update(1.0, *zincr_, 1.0);
  }

  return;
}
