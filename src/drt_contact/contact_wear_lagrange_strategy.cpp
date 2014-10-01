/*!----------------------------------------------------------------------
\file contact_wear_lagrange_srategy.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                   farah 09/13 |
 *----------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"

#include "contact_wear_lagrange_strategy.H"
#include "contact_wear_interface.H"

#include "contact_lagrange_strategy.H"
#include "contact_interface.H"
#include "contact_defines.H"
#include "contact_integrator.H"
#include "contact_element.H"

#include "friction_node.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include <Epetra_FEVector.h>

/*----------------------------------------------------------------------*
 | ctor (public)                                             farah 09/13|
 *----------------------------------------------------------------------*/
CONTACT::WearLagrangeStrategy::WearLagrangeStrategy(DRT::Discretization& probdiscret,
                                                Teuchos::ParameterList params,
                                                std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces,
                                                int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof) :
CoLagrangeStrategy(probdiscret,params,interfaces,dim,comm,alphaf,maxdof),
sswear_(DRT::INPUT::IntegralValue<int>(Params(),"SSWEAR"))
{
  // cast to  wearinterfaces
  for (int z=0; z<(int)interfaces.size();++z)
  {
    interface_.push_back(Teuchos::rcp_dynamic_cast<CONTACT::WearInterface>(interfaces[z]));
    if (interface_[z]==Teuchos::null)
      dserror("WearLagrangeStrategy: Interface-cast failed!");
  }

  // set wear contact status
  INPAR::CONTACT::WearType wtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(Params(),"WEARTYPE");
  INPAR::CONTACT::WearSide wside =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(Params(),"BOTH_SIDED_WEAR");
  INPAR::CONTACT::WearTimeScale wtime =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearTimeScale>(Params(),"WEAR_TIMESCALE");

  if (wtype == INPAR::CONTACT::wear_impl)
    wearimpl_ = true;
  else
    wearimpl_=false;

  // set wear contact discretization
  if (wtype == INPAR::CONTACT::wear_discr)
    weardiscr_ = true;
  else
    weardiscr_=false;

  // both sided wear for discrete wear
  if (wside == INPAR::CONTACT::wear_both_discr)
    wearbothdiscr_=true;
  else
    wearbothdiscr_=false;

  // different wear timescales?
  if (wtime == INPAR::CONTACT::wear_time_different)
    weartimescales_=true;
  else
    weartimescales_=false;

  // max dof number -- disp dofs and lm dofs considered
  maxdofwear_ = maxdof_ + glmdofrowmap_->NumGlobalElements();

  SetupWear(false,true);

  return;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                               farah 09/13 |
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::SetupWear(bool redistributed, bool init)
{
  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gminvolvednodes_  = Teuchos::null;  // all involved master nodes
  gminvolveddofs_   = Teuchos::null;  // all involved master dofs
  gwdofrowmap_      = Teuchos::null;
  gwmdofrowmap_     = Teuchos::null;
  gslipn_           = Teuchos::null;  // vector dummy for wear - slave slip dofs
  gsdofnrowmap_     = Teuchos::null;  // vector dummy for wear - slave all dofs
  gwinact_          = Teuchos::null;  // vector dummy for wear - slave inactive dofss

  gmdofnrowmap_     = Teuchos::null;  // vector dummy for wear - master all dofs
  gmslipn_          = Teuchos::null;
  gwminact_         = Teuchos::null;

  galldofnrowmap_   = Teuchos::null;
  gwalldofrowmap_   = Teuchos::null;

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // ****************************************************
    // for wear as own variable
    // ****************************************************
    if(weardiscr_)
    {
      // build wear dof map
      interface_[i]->UpdateWSets(offset_if,maxdofwear_,wearbothdiscr_);

      // merge interface slave wear dof maps to global slave wear dof map
      gwdofrowmap_ = LINALG::MergeMap(gwdofrowmap_, interface_[i]->WDofs());
      offset_if = gwdofrowmap_->NumGlobalElements();
      if (offset_if < 0) offset_if = 0;

      // merge interface master wear dof maps to global slave wear dof map
      if (wearbothdiscr_)
      {
        gwmdofrowmap_ = LINALG::MergeMap(gwmdofrowmap_, interface_[i]->WMDofs());
        offset_if = gwmdofrowmap_->NumGlobalElements();
        if (offset_if < 0) offset_if = 0;
      }

      //slavenode normal part (first entry)
      interface_[i]->SplitSlaveDofs();
      gsdofnrowmap_ = LINALG::MergeMap(gsdofnrowmap_, interface_[i]->SNDofs());

      //masternode normal part (first entry)
      if (wearbothdiscr_)
      {
        interface_[i]->SplitMasterDofs();
        gmdofnrowmap_ = LINALG::MergeMap(gmdofnrowmap_, interface_[i]->MNDofs());

        interface_[i]->BuildActiveSetMaster();
        gmslipn_= LINALG::MergeMap(gmslipn_,interface_[i]->SlipMasterNDofs(),false);
      }

      // initialize nodal wcurr for integrator (mod. gap)
      for (int j=0; j<(int)interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cnode = dynamic_cast<FriNode*>(node);

        cnode->FriDataPlus().wcurr()[0]=0.0;
      }

      if (wearbothdiscr_)
      {
        for (int j=0; j<(int)interface_[i]->MasterColNodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->MasterColNodes()->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cnode = dynamic_cast<FriNode*>(node);

          cnode->FriDataPlus().wcurr()[0]=0.0;
        }
      }
    }
    // ****************************************************
    // both-sided wear specific
    // ****************************************************
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(Params(),"BOTH_SIDED_WEAR")
        == INPAR::CONTACT::wear_both_map)
    {
      gminvolvednodes_ = LINALG::MergeMap(gminvolvednodes_, interface_[i]->InvolvedNodes(), false);
      gminvolveddofs_  = LINALG::MergeMap(gminvolveddofs_, interface_[i]->InvolvedDofs(), false);
    }
  }

  // get the normal dof of slip nodes for wear definition
  if(weardiscr_)
  {
    gslipn_  =  LINALG::SplitMap(*gslipdofs_,*gslipt_);

    if(sswear_)
      gwinact_ =  LINALG::SplitMap(*gsdofnrowmap_,*gactiven_);
    else
      gwinact_ =  LINALG::SplitMap(*gsdofnrowmap_,*gslipn_);

    if(wearbothdiscr_)
    {
      gwminact_ = LINALG::SplitMap(*gmdofnrowmap_,*gmslipn_);

      // complete wear dofs (s+m)
      galldofnrowmap_ = LINALG::MergeMap(*gsdofnrowmap_,*gmdofnrowmap_,false);
      gwalldofrowmap_ = LINALG::MergeMap(*gwdofrowmap_,*gwmdofrowmap_,false);
    }
  }

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------
  if (weardiscr_)
  {
    // initialize vectors and matrices
    if (!redistributed)
    {
      // setup Lagrange multiplier vectors
      w_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      wincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    }

    // In the redistribution case, first check if the vectors and
    // matrices have already been defined, If yes, transform them
    // to the new redistributed maps. If not, initialize them.
    // Moreover, store redistributed quantities into nodes!!!
    else
    {
      // setup Lagrange multiplier vectors
      if (w_ == Teuchos::null) {
        w_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> neww = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
        LINALG::Export(*w_,*neww);
        w_ = neww;
      }

      if (wincr_ == Teuchos::null) {
        wincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      }
      else
      {
        Teuchos::RCP<Epetra_Vector> newwincr = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
        LINALG::Export(*wincr_,*newwincr);
        wincr_ = newwincr;
      }
    }

    // both sided wear
    if(wearbothdiscr_)
    {
      // initialize vectors and matrices
      if (!redistributed)
      {
        // setup Lagrange multiplier vectors
        wm_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        wmincr_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
      }

      // In the redistribution case, first check if the vectors and
      // matrices have already been defined, If yes, transform them
      // to the new redistributed maps. If not, initialize them.
      // Moreover, store redistributed quantities into nodes!!!
      else
      {
        // setup Lagrange multiplier vectors
        if (wm_ == Teuchos::null) {
          wm_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        }
        else
        {
          Teuchos::RCP<Epetra_Vector> neww = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
          LINALG::Export(*wm_,*neww);
          wm_ = neww;
        }

        if (wmincr_ == Teuchos::null) {
          wmincr_ = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
        }
        else
        {
          Teuchos::RCP<Epetra_Vector> newwincr = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
          LINALG::Export(*wmincr_,*newwincr);
          wmincr_ = newwincr;
        }
      }
    }
  }

  // output wear ... this is for the unweighted wear*n vector
  wearoutput_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  wearoutput2_ = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));


  return;
}

/*----------------------------------------------------------------------*
 | initialize wear stuff for next Newton step                farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::InitMortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  if (IsSelfContact()) UpdateMasterSlaveSetsGlobal();

  /**********************************************************************/
  /* initialize Dold and Mold if not done already                       */
  /**********************************************************************/
  if (dold_==Teuchos::null)
  {
    dold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    dold_->Zero();
    dold_->Complete();
  }
  if (mold_==Teuchos::null)
  {
    mold_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));
    mold_->Zero();
    mold_->Complete(*gmdofrowmap_,*gsdofrowmap_);
  }

  /**********************************************************************/
  /* (re)setup global Mortar LINALG::SparseMatrices and Epetra_Vectors  */
  /**********************************************************************/
  dmatrix_  = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  d2matrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  //d2matrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,10));
  mmatrix_  = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100));

  // global gap
  g_        = LINALG::CreateVector(*gsnoderowmap_, true);

  /**********************************************************************/
  /* (re)setup global wear Epetra_Vector (for all wear problems)        */
  /**********************************************************************/
  if (wear_) wearvector_ = LINALG::CreateVector(*gsnoderowmap_, true);

  /**********************************************************************/
  /* (re)setup global matrix A for tsi problems                         */
  /**********************************************************************/
  if (tsi_) amatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));

  /**********************************************************************/
  /* in the case of dual quad 3D, the modified D matrices are setup     */
  /**********************************************************************/
  if (friction_ && Dualquadslave3d())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_==Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
      doldmod_->Zero();
      doldmod_->Complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble wear stuff for next Newton step                  farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::AssembleMortar()
{
  // for all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //************************************************************
    // assemble D-, M-matrix and g-vector, store them globally
    interface_[i]->AssembleDM(*dmatrix_,*mmatrix_);
    interface_[i]->AssembleG(*g_);

    //************************************************************
    // only assemble D2 for both-sided wear --> unweights the
    // weighted wear increment in master side
    // --> based on weak dirichlet bc!
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(Params(),"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
      interface_[i]->AssembleD2(*d2matrix_);

    //************************************************************
    // assemble wear vector
    if (wear_) interface_[i]->AssembleWear(*wearvector_);

    //************************************************************
    // assemble matrix A for tsi with frictional contact
    if (tsi_ and friction_) interface_[i]->AssembleA(*amatrix_);

#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
    interface_[i]->FDCheckNormalDeriv();
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDNORMAL

#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
    dmatrix_->Complete();
    if( dmatrix_->NormOne() )
      interface_[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if( mmatrix_->NormOne() )
        interface_[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARM
  } // end interface loop

  // *********************************************************************************
  // modify gap vector towards wear, only if no structure with ale is applied
  // This additional gap is also required for an implicit ALE-wear algorithm !!!
  // THIS IS THE EXPLICIT WEAR ALGORITHM
  // wearvector_ only updated at the end of a time step --> this newton-step-wise
  // update is not elegant!
  // *********************************************************************************
  if (wear_ and Params().get<int>("PROBTYPE")!=INPAR::CONTACT::structalewear and !wearimpl_)
  {
    g_->Update(1.0,*wearvector_,1.0);
  }
  else
  {
    if (wearimpl_)
    {
      // update the gap function with the current wear increment
      // this is for the implicit wear algorithm
      // we have to update the wear-increment in every newton-step and not just
      // after a time step!
      StoreNodalQuantities(MORTAR::StrategyBase::wear);
      interface_[0]->AssembleWear(*wearvector_);
      g_->Update(1.0,*wearvector_,1.0);

      //update alle gap function entries for slave nodes!
      for (int i=0; i<(int)interface_.size(); ++i)
      {
        for (int j=0; j<(int)interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->SlaveRowNodes()->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cnode = dynamic_cast<FriNode*>(node);

          if(cnode->FriData().Slip()==true)
            cnode->CoData().Getg()+=cnode->FriDataPlus().Wear();
        }
      }
    }
  }

  //********************************************
  // FillComplete() global Mortar matrices     *
  //********************************************
  dmatrix_->Complete();
  mmatrix_->Complete(*gmdofrowmap_,*gsdofrowmap_);

  //********************************************
  // FillComplete() matrix for both-sided wear *
  //********************************************
  d2matrix_->Complete(*gmdofrowmap_,*gmdofrowmap_);

  //********************************************
  // FillComplete() matrix for tsi problems    *
  //********************************************
  if (tsi_) amatrix_->Complete();

  return;
}

/*----------------------------------------------------------------------*
 | initialize global contact variables for next Newton step  farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::Initialize()
{
  // (re)setup global matrices containing fc derivatives
  // must use FE_MATRIX type here, as we will do non-local assembly!
  lindmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  linmmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // (re)setup global tangent matrix
  tmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactivet_,3));

  // (re)setup global matrix containing gap derivatives
  smatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3));

  // inactive rhs for the saddle point problem
  Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_, *gactivedofs_);
  inactiverhs_ = LINALG::CreateVector(*gidofs, true);

  // (re)setup of global friction
  // here the calculation of gstickt is necessary
  Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
  linstickLM_  = Teuchos::rcp(new LINALG::SparseMatrix(*gstickt,3));
  linstickDIS_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickt,3));
  linstickRHS_ = LINALG::CreateVector(*gstickt,true);

  linslipLM_   = Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));
  linslipDIS_  = Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));
  linslipRHS_  = LINALG::CreateVector(*gslipt_,true);

  if (wearimpl_)
  {
    // create matrices for implicite wear: these matrices are due to
    // the gap-change in the compl. fnc.
    // Here are only the lin. w.r.t. the lagr. mult.
    wlinmatrix_   = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3));
    wlinmatrixsl_ = Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));

  #ifdef CONSISTENTSTICK
    wlinmatrixst_ = Teuchos::rcp(new LINALG::SparseMatrix(*gstickt,3));
  #endif
  }

  if(weardiscr_)
  {
    if(sswear_)
    {
      // basic matrices
      twmatrix_        = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,100)); //gsdofnrowmap_
      ematrix_         = Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,100)); //gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdis_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      lintlm_          =  Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      linedis_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
    }
    else
    {
      // basic matrices
      twmatrix_        = Teuchos::rcp(new LINALG::SparseMatrix(*gslipn_,100)); //gsdofnrowmap_
      ematrix_         = Teuchos::rcp(new LINALG::SparseMatrix(*gslipn_,100)); //gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdis_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      lintlm_          =  Teuchos::rcp(new LINALG::SparseMatrix(*gslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      linedis_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
    }

    // linearizations w.r.t w
    smatrixW_        =  Teuchos::rcp(new LINALG::SparseMatrix(*gactiven_,3)); //gactiven_
    linslipW_        =  Teuchos::rcp(new LINALG::SparseMatrix(*gslipt_,3));

    // w - rhs
    inactiveWearRhs_ = LINALG::CreateVector(*gwinact_, true);

    if(sswear_)
      WearCondRhs_     = LINALG::CreateVector(*gactiven_, true);
    else
      WearCondRhs_     = LINALG::CreateVector(*gslipn_, true);

    // both-sided discr wear
    if(wearbothdiscr_)
    {
      // basic matrices
      twmatrixM_        = Teuchos::rcp(new LINALG::SparseMatrix(*gmslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX)); //gsdofnrowmap_
      ematrixM_         = Teuchos::rcp(new LINALG::SparseMatrix(*gmslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX)); //gsdofnrowmap_

      // linearizations w.r.t d and z
      lintdisM_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gmslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      lintlmM_          =  Teuchos::rcp(new LINALG::SparseMatrix(*gmslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
      linedisM_         =  Teuchos::rcp(new LINALG::SparseMatrix(*gmslipn_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

      // w - rhs
      inactiveWearRhsM_ = Teuchos::rcp(new Epetra_FEVector(*gwminact_));
      WearCondRhsM_     = Teuchos::rcp(new Epetra_FEVector(*gmslipn_));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | condense wear and lm. for impl/expl wear algorithm        farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::CondenseWearImplExpl(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                   Teuchos::RCP<Epetra_Vector>& feff,
                                                   Teuchos::RCP<Epetra_Vector>& gact)
{
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");

  // double-check if this is a dual LM system
  if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    dserror("Condensation only for dual LM");

  // get stick map
  Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);

  /********************************************************************/
  /* (1) Multiply Mortar matrices: m^ = inv(d) * m                    */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> invd = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication mhat = inv(D) * M
  mhatmatrix_ = LINALG::MLMultiply(*invd,false,*mmatrix_,false,false,false,true);

  /********************************************************************/
  /* (2) Add contact stiffness terms to kteff                         */
  /********************************************************************/

  // transform if necessary
  if (ParRedist())
  {
    lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
    linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
  }

  kteff->UnComplete();
  kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
  kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
  kteff->Complete();

  /********************************************************************/
  /* (3) Split kteff into 3x3 matrix blocks                           */
  /********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);
  if (ParRedist())
  {
    // split and transform to redistributed maps
    LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  }

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  /********************************************************************/
  /* (4) Split feff into 3 subvectors                                 */
  /********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (ParRedist())
  {
    // split and transform to redistributed maps
    LINALG::SplitVector(*ProblemDofs(),*feff,pgsmdofrowmap_,fsm,gndofrowmap_,fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    LINALG::Export(*fsm,*fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitVector(*ProblemDofs(),*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);
  }

  // abbreviations for slave and master set
  int sset = gsdofrowmap_->NumGlobalElements();
  int mset = gmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;

  //--------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //--------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  //--------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> temp3 = LINALG::MLMultiply(*trafo_,false,*invd_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> temp4 = LINALG::MLMultiply(*trafo_,false,*mhatmatrix_,false,false,false,true);
    dmatrix_    = temp2;
    invd_       = temp3;
    mhatmatrix_ = temp4;
  }

  /********************************************************************/
  /* (5) Split slave quantities into active / inactive, stick / slip  */
  /********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split kaa into 2 groups sl,st = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> temp1map;
  Teuchos::RCP<LINALG::SparseMatrix> temp1mtx4,temp1mtx5;

  // we will get the stick rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gstdofs;

  LINALG::SplitMatrix2x2(kaa,gactivedofs_,gidofs,gstdofs,gslipdofs_,kast,kasl,temp1mtx4,temp1mtx5);

  // abbreviations for active and inactive set, stick and slip set
  int aset = gactivedofs_->NumGlobalElements();
  int iset = gidofs->NumGlobalElements();
  int stickset = gstdofs->NumGlobalElements();
  int slipset = gslipdofs_->NumGlobalElements();

  // we want to split fs into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*gsdofrowmap_,*fs,gactivedofs_,fa,gidofs,fi);

  // we want to split fa into 2 groups sl,st
  Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
  Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

  // do the vector splitting a -> sl+st
  if(aset)
    LINALG::SplitVector(*gactivedofs_,*fa,gslipdofs_,fsl,gstdofs,fst);

  /********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix             */
  /********************************************************************/

  // active, stick and slip part of invd
  Teuchos::RCP<LINALG::SparseMatrix> invda, invdsl, invdst;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gstdofs,gslipdofs_,invdst,tempmtx1,tempmtx2,tempmtx3);

  // for Implicit Wear
  Teuchos::RCP<LINALG::SparseMatrix> wa, wi, wsl, wst;          // split lin.matrix for Cgap
  Teuchos::RCP<LINALG::SparseMatrix> wsla, wsli, wslsl, wslst;  // split lin.matrix for Cslip
  Teuchos::RCP<Epetra_Vector> fw = Teuchos::rcp(new Epetra_Vector(*gactiven_));
  Teuchos::RCP<Epetra_Vector> fwsl = Teuchos::rcp(new Epetra_Vector(*gslipt_));

  if (wearimpl_)
  {
    Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));

    LINALG::SplitVector( *gsdofrowmap_, *z_,gactivedofs_, za, gidofs, zi);

    // split wlinmatrix into wa for active dofs
    LINALG::SplitMatrix2x2(wlinmatrix_,gactiven_,gactiven_,gactivedofs_,gidofs,wa,tempmtx1,tempmtx2,tempmtx3);

    // split wlinmatrixsl into wa for active dofs
    LINALG::SplitMatrix2x2(wlinmatrixsl_,gslipt_,gslipt_,gactivedofs_,gidofs,wsla,tempmtx1,tempmtx2,tempmtx3);

    // WEAR LIN. W.R.T. Z *****************************************************
    // for normal contact
    if (aset)
      wa->Multiply(false,*za,*fw);

    // WEAR LIN. W.R.T. Z *****************************************************
    // for slip contact
    if (slipset)
      wsla->Multiply(false,*za,*fwsl);
  }

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<LINALG::SparseMatrix> dai;
  LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,tempmtx1,dai,tempmtx2,tempmtx3);

   // do the multiplication dhat = invda * dai
  Teuchos::RCP<LINALG::SparseMatrix> dhat = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
  if (aset && iset) dhat = LINALG::MLMultiply(*invda,false,*dai,false,false,false,true);
  dhat->Complete(*gidofs,*gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixa;
  LINALG::SplitMatrix2x2(mmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mmatrixa,tempmtx1,tempmtx2,tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<LINALG::SparseMatrix> mhataam = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
  if (aset) mhataam = LINALG::MLMultiply(*invda,false,*mmatrixa,false,false,false,true);
  mhataam->Complete(*gmdofrowmap_,*gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<LINALG::SparseMatrix> mhata;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);

  // scaling of invd and dai
  invda->Scale(1/(1-alphaf_));
  invdsl->Scale(1/(1-alphaf_));
  invdst->Scale(1/(1-alphaf_));
  dai->Scale(1-alphaf_);

  /********************************************************************/
  /* (7) Build the final K blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //-------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<LINALG::SparseMatrix> kmnmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  kmnmod->Add(*kmn,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhataam,true,*kan,false,false,false,true);
  kmnmod->Add(*kmnadd,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<LINALG::SparseMatrix> kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  kmmmod->Add(*kmm,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhataam,true,*kam,false,false,false,true);
  kmmmod->Add(*kmmadd,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<LINALG::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmimod->Add(*kmi,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhataam,true,*kai,false,false,false,true);
    kmimod->Add(*kmiadd,false,1.0,1.0);
    kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<LINALG::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmamod->Add(*kma,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhataam,true,*kaa,false,false,false,true);
    kmamod->Add(*kmaadd,false,1.0,1.0);
    kmamod->Complete(kma->DomainMap(),kma->RowMap());
  }

  //--------------------------------------------------------- THIRD LINE
  // kin: subtract T(dhat)*kan
  Teuchos::RCP<LINALG::SparseMatrix> kinmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
  kinmod->Add(*kin,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kinadd = LINALG::MLMultiply(*dhat,true,*kan,false,false,false,true);
  kinmod->Add(*kinadd,false,-1.0,1.0);
  kinmod->Complete(kin->DomainMap(),kin->RowMap());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<LINALG::SparseMatrix> kimmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
  kimmod->Add(*kim,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kimadd = LINALG::MLMultiply(*dhat,true,*kam,false,false,false,true);
  kimmod->Add(*kimadd,false,-1.0,1.0);
  kimmod->Complete(kim->DomainMap(),kim->RowMap());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<LINALG::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kiimod->Add(*kii,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kiiadd = LINALG::MLMultiply(*dhat,true,*kai,false,false,false,true);
    kiimod->Add(*kiiadd,false,-1.0,1.0);
    kiimod->Complete(kii->DomainMap(),kii->RowMap());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<LINALG::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kiamod->Add(*kia,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kiaadd = LINALG::MLMultiply(*dhat,true,*kaa,false,false,false,true);
    kiamod->Add(*kiaadd,false,-1.0,1.0);
    kiamod->Complete(kia->DomainMap(),kia->RowMap());
  }

  //-------------------------------------------------------- FOURTH LINE

  // create matrices for implicit wear
  Teuchos::RCP<LINALG::SparseMatrix> kgnmod;
  Teuchos::RCP<LINALG::SparseMatrix> kgmmod;
  Teuchos::RCP<LINALG::SparseMatrix> kgimod;
  Teuchos::RCP<LINALG::SparseMatrix> kgaamod;

  if (wearimpl_)
  {
    if (aset)
    {
      kgnmod = LINALG::MLMultiply(*wa,false,*invda,true,false,false,true);
      kgnmod = LINALG::MLMultiply(*kgnmod,false,*kan,false,false,false,true);
    }

    if (aset)
    {
      kgmmod = LINALG::MLMultiply(*wa,false,*invda,true,false,false,true);
      kgmmod = LINALG::MLMultiply(*kgmmod,false,*kam,false,false,false,true);
    }

    if (aset and iset)
    {
      kgimod = LINALG::MLMultiply(*wa,false,*invda,true,false,false,true);
      kgimod = LINALG::MLMultiply(*kgimod,false,*kai,false,false,false,true);
    }

    if (aset)
    {
      kgaamod = LINALG::MLMultiply(*wa,false,*invda,true,false,false,true);
      kgaamod = LINALG::MLMultiply(*kgaamod,false,*kaa,false,false,false,true);
    }
  }
  //--------------------------------------------------------- FIFTH LINE
  // blocks for complementary conditions (stick nodes)

  // kstn: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstnmod;
  if (stickset)
  {
    kstnmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstnmod = LINALG::MLMultiply(*kstnmod,false,*kan,false,false,false,true);
  }

  // kstm: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstmmod;
  if(stickset)
  {
    kstmmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstmmod = LINALG::MLMultiply(*kstmmod,false,*kam,false,false,false,true);
  }

  // ksti: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstimod;
  if(stickset && iset)
  {
    kstimod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstimod = LINALG::MLMultiply(*kstimod,false,*kai,false,false,false,true);
  }

  // kstsl: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstslmod;
  if(stickset && slipset)
  {
    kstslmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstslmod = LINALG::MLMultiply(*kstslmod,false,*kasl,false,false,false,true);
  }

  // kststmod: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kststmod;
  if (stickset)
  {
    kststmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kststmod = LINALG::MLMultiply(*kststmod,false,*kast,false,false,false,true);
  }

  //--------------------------------------------------------- SIXTH LINE
  // blocks for complementary conditions (slip nodes)

  // ksln: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslnmod;
  if(slipset)
  {
    kslnmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslnmod = LINALG::MLMultiply(*kslnmod,false,*kan,false,false,false,true);
  }

  // kslm: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslmmod;
  if(slipset)
  {
    kslmmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslmmod = LINALG::MLMultiply(*kslmmod,false,*kam,false,false,false,true);
  }

  // ksli: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslimod;
  if (slipset && iset)
  {
    kslimod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslimod = LINALG::MLMultiply(*kslimod,false,*kai,false,false,false,true);
  }

  // kslsl: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslslmod;
  if(slipset)
  {
    kslslmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslslmod = LINALG::MLMultiply(*kslslmod,false,*kasl,false,false,false,true);
  }

  // slstmod: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslstmod;
  if (slipset && stickset)
  {
    kslstmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslstmod = LINALG::MLMultiply(*kslstmod,false,*kast,false,false,false,true);
  }

  // create matrices for implicit wear
  Teuchos::RCP<LINALG::SparseMatrix> kslwnmod;
  Teuchos::RCP<LINALG::SparseMatrix> kslwmmod;
  Teuchos::RCP<LINALG::SparseMatrix> kslwimod;
  Teuchos::RCP<LINALG::SparseMatrix> kslwaamod;

  if (wearimpl_)
  {
    if (slipset)
    {
      kslwnmod = LINALG::MLMultiply(*wsla,false,*invda,true,false,false,true);
      kslwnmod = LINALG::MLMultiply(*kslwnmod,false,*kan,false,false,false,true);
    }

    if (slipset)
    {
      kslwmmod = LINALG::MLMultiply(*wsla,false,*invda,true,false,false,true);
      kslwmmod = LINALG::MLMultiply(*kslwmmod,false,*kam,false,false,false,true);
    }

    if (slipset and iset)
    {
      kslwimod = LINALG::MLMultiply(*wsla,false,*invda,true,false,false,true);
      kslwimod = LINALG::MLMultiply(*kslwimod,false,*kai,false,false,false,true);
    }

    if (slipset)
    {
      kslwaamod = LINALG::MLMultiply(*wsla,false,*invda,true,false,false,true);
      kslwaamod = LINALG::MLMultiply(*kslwaamod,false,*kaa,false,false,false,true);
    }
  }

  /********************************************************************/
  /* (8) Build the final f blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // fm: add alphaf * old contact forces (t_n)
  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(mold_->RowMap()));
    if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    mold_->Multiply(true,*zoldexp,*tempvecm2);
    if (mset) LINALG::Export(*tempvecm2,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->Multiply(true,*zold_,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    dold_->Multiply(true,*zoldexp,*tempvec);
    if (sset) LINALG::Export(*tempvec,*fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    dold_->Multiply(true,*zold_,*fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*fsadd,*faadd);
    fa->Update(-alphaf_,*faadd,1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  //--------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*fsadd,*fiadd);
    fi->Update(-alphaf_,*fiadd,1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (aset) dhat->Multiply(true,*fa,*fimod);
  fimod->Update(1.0,*fi,-1.0);

  //-------------------------------------------------------- FOURTH LINE
  Teuchos::RCP<Epetra_Vector> fgmod;
  if (wearimpl_)
  {
    if (aset)
    {
      fgmod = Teuchos::rcp(new Epetra_Vector(*gactiven_));
      Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*wa,false,*invda,true,false,false,true);
      temp1->Multiply(false,*fa,*fgmod);
    }
  }
  //--------------------------------------------------------- FIFTH LINE
  Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);  // get global stick dofs

  // split the lagrange multiplier vector in stick and slip part
  Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
  Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
  Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

  LINALG::SplitVector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
  LINALG::SplitVector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
  Teuchos::RCP<Epetra_Vector> tempvec1;

  // fst: mutliply with linstickLM
  Teuchos::RCP<Epetra_Vector> fstmod;
  if (stickset)
  {
    fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    temp1->Multiply(false,*fa,*fstmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

    linstickLM_->Multiply(false, *zst, *tempvec1);
    fstmod->Update(-1.0,*tempvec1,1.0);
  }

  //--------------------------------------------------------- SIXTH LINE
  // fsl: mutliply with linslipLM
  Teuchos::RCP<Epetra_Vector> fslmod;
  Teuchos::RCP<Epetra_Vector> fslwmod;

  if (slipset)
  {
    fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
    Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    temp->Multiply(false,*fa,*fslmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

    linslipLM_->Multiply(false, *zsl, *tempvec1);

    fslmod->Update(-1.0,*tempvec1,1.0);

    if(wearimpl_)
    {
      fslwmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
      Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*wsla,false,*invda,true,false,false,true);
      temp2->Multiply(false,*fa,*fslwmod);
    }
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
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
    if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
      kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
      kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
      if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset)
    {
      smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);

      if(wearimpl_)
      {
        kgnmod = MORTAR::MatrixRowTransform(kgnmod,pgsdofrowmap_);
        kgmmod = MORTAR::MatrixRowTransform(kgmmod,pgsdofrowmap_);
        if (iset) kgimod = MORTAR::MatrixRowTransform(kgimod,pgsdofrowmap_);
        kgaamod = MORTAR::MatrixRowTransform(kgaamod,pgsdofrowmap_);
      }
    }

    //----------------------------------------------------------- FIFTH LINE
    if (stickset)
    {
      kstnmod = MORTAR::MatrixRowTransform(kstnmod,pgsdofrowmap_);
      kstmmod = MORTAR::MatrixRowTransform(kstmmod,pgsdofrowmap_);
      if (iset) kstimod = MORTAR::MatrixRowTransform(kstimod,pgsdofrowmap_);
      if (slipset) kstslmod = MORTAR::MatrixRowTransform(kstslmod,pgsdofrowmap_);
      kststmod = MORTAR::MatrixRowTransform(kststmod,pgsdofrowmap_);
      linstickDIS_ = MORTAR::MatrixRowTransform(linstickDIS_,pgsdofrowmap_);
    }

    //----------------------------------------------------------- SIXTH LINE
    if (slipset)
    {
      kslnmod = MORTAR::MatrixRowTransform(kslnmod,pgsdofrowmap_);
      kslmmod = MORTAR::MatrixRowTransform(kslmmod,pgsdofrowmap_);
      if (iset) kslimod = MORTAR::MatrixRowTransform(kslimod,pgsdofrowmap_);
      if (stickset) kslstmod = MORTAR::MatrixRowTransform(kslstmod,pgsdofrowmap_);
      kslslmod = MORTAR::MatrixRowTransform(kslslmod,pgsdofrowmap_);
      linslipDIS_ = MORTAR::MatrixRowTransform(linslipDIS_,pgsdofrowmap_);

      if(wearimpl_)
      {
        kslwnmod = MORTAR::MatrixRowTransform(kslwnmod,pgsdofrowmap_);
        kslwmmod = MORTAR::MatrixRowTransform(kslwmmod,pgsdofrowmap_);
        if (iset) kslwimod = MORTAR::MatrixRowTransform(kslwimod,pgsdofrowmap_);
        kslwaamod = MORTAR::MatrixRowTransform(kslwaamod,pgsdofrowmap_);
      }
    }
  }

  /********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                */
  /********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> kteffnew = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(),81,true,false,kteffmatrix->GetMatrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*ProblemDofs());

  //--------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  if (sset) kteffnew->Add(*kns,false,1.0,1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kmimod,false,1.0,1.0);
  if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (iset && aset) kteffnew->Add(*kiamod,false,1.0,1.0);

  //-------------------------------------------------------- FOURTH LINE

  // add a submatrices to kteffnew
  if (aset) kteffnew->Add(*smatrix_,false,1.0,1.0);

  if(wearimpl_)
  {
    //if (aset) kteffnew->Add(*wlinmatrix_,false,1.0,1.0);
    if (aset) kteffnew->Add(*kgnmod,false,-1.0,1.0);
    if (aset) kteffnew->Add(*kgmmod,false,-1.0,1.0);
    if (aset and iset) kteffnew->Add(*kgimod,false,-1.0,1.0);
    if (aset) kteffnew->Add(*kgaamod,false,-1.0,1.0);
    //if (aset and stickset) kteffnew->Add(*kgstmod,false,1.0,1.0);
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st submatrices to kteffnew
  if (stickset) kteffnew->Add(*kstnmod,false,1.0,1.0);
  if (stickset) kteffnew->Add(*kstmmod,false,1.0,1.0);
  if (stickset && iset) kteffnew->Add(*kstimod,false,1.0,1.0);
  if (stickset && slipset) kteffnew->Add(*kstslmod,false,1.0,1.0);
  if (stickset) kteffnew->Add(*kststmod,false,1.0,1.0);

  // add terms of linearization of sick condition to kteffnew
  if (stickset) kteffnew->Add(*linstickDIS_,false,-1.0,1.0);

  //--------------------------------------------------------- SIXTH LINE
  // add sl submatrices to kteffnew
  if (slipset) kteffnew->Add(*kslnmod,false,1.0,1.0);
  if (slipset) kteffnew->Add(*kslmmod,false,1.0,1.0);
  if (slipset && iset) kteffnew->Add(*kslimod,false,1.0,1.0);
  if (slipset) kteffnew->Add(*kslslmod,false,1.0,1.0);
  if (slipset && stickset) kteffnew->Add(*kslstmod,false,1.0,1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (slipset) kteffnew->Add(*linslipDIS_,false,-1.0,+1.0);

  if(wearimpl_)
  {
    //if (slipset) kteffnew->Add(*wlinmatrixsl_,false,1.0,1.0);
    if (slipset) kteffnew->Add(*kslwnmod,false,1.0,1.0);
    if (slipset) kteffnew->Add(*kslwmmod,false,1.0,1.0);
    if (slipset and iset) kteffnew->Add(*kslwimod,false,1.0,1.0);
    if (slipset) kteffnew->Add(*kslwaamod,false,1.0,1.0);
    //if (slipset and stickset) kteffnew->Add(*kslwstmod,false,1.0,1.0);
  }

  // FillComplete kteffnew (square)
  kteffnew->Complete();

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fimod,*fimodexp);
    feffnew->Update(1.0,*fimodexp,1.0);
  }

  //-------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  Teuchos::RCP<Epetra_Vector> fwexp;
  Teuchos::RCP<Epetra_Vector> fgmodexp;

  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*gact,*gexp);
    feffnew->Update(-1.0,*gexp,1.0);

    if(wearimpl_)
    {
      //commented due to incremental solution algorithm.
      fwexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fw,*fwexp);
      feffnew->Update(+1.0,*fwexp,1.0);

      fgmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fgmod,*fgmodexp);
      feffnew->Update(-1.0,*fgmodexp,+1.0);
    }
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fstmodexp;
  if (stickset)
  {
    fstmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fstmod,*fstmodexp);
    feffnew->Update(1.0,*fstmodexp,+1.0);
  }

  // add terms of linearization feffnew
   if (stickset)
   {
      Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*linstickRHS_,*linstickRHSexp);
      feffnew->Update(-1.0,*linstickRHSexp,1.0);
   }

  //--------------------------------------------------------- SIXTH LINE

  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fslmodexp;
  Teuchos::RCP<Epetra_Vector> fwslexp;
  Teuchos::RCP<Epetra_Vector> fslwmodexp;


  if (slipset)
  {
    fslmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fslmod,*fslmodexp);
    feffnew->Update(1.0,*fslmodexp,1.0);
  }

  if(wearimpl_)
  {
    if (slipset)
    {
      fslwmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fslwmod,*fslwmodexp);
      feffnew->Update(+1.0,*fslwmodexp,1.0);
    }
    //commented due to incremental solution algorithm
    if (slipset)
    {
      fwslexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fwsl,*fwslexp);
      feffnew->Update(-1.0,*fwslexp,1.0);
    }
  }

  if (slipset)
  {
    Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*linslipRHS_,*linslipRHSexp);
    feffnew->Update(-1.0,*linslipRHSexp,1.0);
  }

  // finally do the replacement
  kteff = kteffnew;
  feff = feffnew;

  return;
}

/*----------------------------------------------------------------------*
 | condense wear and lm. for discr. wear algorithm           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::CondenseWearDiscr(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                   Teuchos::RCP<Epetra_Vector>& feff,
                                                   Teuchos::RCP<Epetra_Vector>& gact)
{
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  INPAR::CONTACT::WearShape wearshapefcn = DRT::INPUT::IntegralValue<INPAR::CONTACT::WearShape>(Params(),"WEAR_SHAPEFCN");

  // double-check if this is a dual LM system
  if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
    dserror("Condensation only for dual LM");
  if (wearshapefcn != INPAR::CONTACT::wear_shape_dual)
    dserror("Condensation only for dual wear");

  // get stick map
  Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);

  /********************************************************************/
  /* (1a) Multiply Mortar matrices: m^ = inv(d) * m                    */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> invd = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix_));
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gsdofrowmap_,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);
  // we cannot use this check, as we deliberately replaced zero entries
  //if (err>0) dserror("ERROR: ReplaceDiagonalValues: Missing diagonal entry!");

  // do the multiplication mhat = inv(D) * M
  mhatmatrix_ = LINALG::MLMultiply(*invd,false,*mmatrix_,false,false,false,true);

  /********************************************************************/
  /* (2) Add contact stiffness terms to kteff                         */
  /********************************************************************/
  // transform if necessary
  if (ParRedist())
  {
    lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
    linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
  }

  kteff->UnComplete();
  kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
  kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
  kteff->Complete();

  /********************************************************************/
  /* (3a) Split kteff into 3x3 matrix blocks                          */
  /********************************************************************/

  // we want to split k into 3 groups s,m,n = 9 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> tempmap;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
  Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;

  // split into slave/master part + structure part
  Teuchos::RCP<LINALG::SparseMatrix> kteffmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(kteff);
  if (ParRedist())
  {
    // split and transform to redistributed maps
    LINALG::SplitMatrix2x2(kteffmatrix,pgsmdofrowmap_,gndofrowmap_,pgsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
    ksmsm = MORTAR::MatrixRowColTransform(ksmsm,gsmdofrowmap_,gsmdofrowmap_);
    ksmn  = MORTAR::MatrixRowTransform(ksmn,gsmdofrowmap_);
    knsm  = MORTAR::MatrixColTransform(knsm,gsmdofrowmap_);
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitMatrix2x2(kteffmatrix,gsmdofrowmap_,gndofrowmap_,gsmdofrowmap_,gndofrowmap_,ksmsm,ksmn,knsm,knn);
  }

  // further splits into slave part + master part
  LINALG::SplitMatrix2x2(ksmsm,gsdofrowmap_,gmdofrowmap_,gsdofrowmap_,gmdofrowmap_,kss,ksm,kms,kmm);
  LINALG::SplitMatrix2x2(ksmn,gsdofrowmap_,gmdofrowmap_,gndofrowmap_,tempmap,ksn,tempmtx1,kmn,tempmtx2);
  LINALG::SplitMatrix2x2(knsm,gndofrowmap_,tempmap,gsdofrowmap_,gmdofrowmap_,kns,knm,tempmtx1,tempmtx2);

  /********************************************************************/
  /* (4) Split feff into 3 subvectors                                 */
  /********************************************************************/

  // we want to split f into 3 groups s.m,n
  Teuchos::RCP<Epetra_Vector> fs, fm, fn;

  // temporarily we need the group sm
  Teuchos::RCP<Epetra_Vector> fsm;

  // do the vector splitting smn -> sm+n
  if (ParRedist())
  {
    // split and transform to redistributed maps
    LINALG::SplitVector(*ProblemDofs(),*feff,pgsmdofrowmap_,fsm,gndofrowmap_,fn);
    Teuchos::RCP<Epetra_Vector> fsmtemp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    LINALG::Export(*fsm,*fsmtemp);
    fsm = fsmtemp;
  }
  else
  {
    // only split, no need to transform
    LINALG::SplitVector(*ProblemDofs(),*feff,gsmdofrowmap_,fsm,gndofrowmap_,fn);
  }

  // abbreviations for slave and master set
  int sset = gsdofrowmap_->NumGlobalElements();
  int mset = gmdofrowmap_->NumGlobalElements();

  // we want to split fsm into 2 groups s,m
  fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

  // do the vector splitting sm -> s+m
  LINALG::SplitVector(*gsmdofrowmap_,*fsm,gsdofrowmap_,fs,gmdofrowmap_,fm);

  // store some stuff for static condensation of LM
  fs_   = fs;
  invd_ = invd;
  ksn_  = ksn;
  ksm_  = ksm;
  kss_  = kss;

  //--------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //--------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  // D^(-1)    ---->   T * D^(-1)
  // \hat{M}   ---->   T * \hat{M}
  //--------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_, invd_ and mhatmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> temp3 = LINALG::MLMultiply(*trafo_,false,*invd_,false,false,false,true);
    Teuchos::RCP<LINALG::SparseMatrix> temp4 = LINALG::MLMultiply(*trafo_,false,*mhatmatrix_,false,false,false,true);
    dmatrix_    = temp2;
    invd_       = temp3;
    mhatmatrix_ = temp4;
  }

  /********************************************************************/
  /* (5a) Split slave quantities into active / inactive, stick / slip  */
  /********************************************************************/
  // we want to split kssmod into 2 groups a,i = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;

  // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;

  // we will get the i rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gidofs;

  // do the splitting
  LINALG::SplitMatrix2x2(kss,gactivedofs_,gidofs,gactivedofs_,gidofs,kaa,kai,kia,kii);
  LINALG::SplitMatrix2x2(ksn,gactivedofs_,gidofs,gndofrowmap_,tempmap,kan,tempmtx1,kin,tempmtx2);
  LINALG::SplitMatrix2x2(ksm,gactivedofs_,gidofs,gmdofrowmap_,tempmap,kam,tempmtx1,kim,tempmtx2);
  LINALG::SplitMatrix2x2(kms,gmdofrowmap_,tempmap,gactivedofs_,gidofs,kma,kmi,tempmtx1,tempmtx2);

  // we want to split kaa into 2 groups sl,st = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> kslsl, kslst, kstsl, kstst, kast, kasl;

  // we want to split kan / kam / kai into 2 groups sl,st = 2 blocks
  Teuchos::RCP<LINALG::SparseMatrix> ksln, kstn, kslm, kstm, ksli, ksti;

  // some temporary Teuchos::RCPs
  Teuchos::RCP<Epetra_Map> temp1map;
  Teuchos::RCP<LINALG::SparseMatrix> temp1mtx4,temp1mtx5;

  // we will get the stick rowmap as a by-product
  Teuchos::RCP<Epetra_Map> gstdofs;

  LINALG::SplitMatrix2x2(kaa,gactivedofs_,gidofs,gstdofs,gslipdofs_,kast,kasl,temp1mtx4,temp1mtx5);

  // abbreviations for active and inactive set, stick and slip set
  int aset = gactivedofs_->NumGlobalElements();
  int iset = gidofs->NumGlobalElements();
  int stickset = gstdofs->NumGlobalElements();
  int slipset = gslipdofs_->NumGlobalElements();
  gidofs_ = gidofs;

  // we want to split fs into 2 groups a,i
  Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*gidofs));

  // do the vector splitting s -> a+i
  LINALG::SplitVector(*gsdofrowmap_,*fs,gactivedofs_,fa,gidofs,fi);

  // we want to split fa into 2 groups sl,st
  Teuchos::RCP<Epetra_Vector> fsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));
  Teuchos::RCP<Epetra_Vector> fst = Teuchos::rcp(new Epetra_Vector(*gstdofs));

  // do the vector splitting a -> sl+st
  if(aset)
    LINALG::SplitVector(*gactivedofs_,*fa,gslipdofs_,fsl,gstdofs,fst);

  /********************************************************************/
  /* (6) Isolate necessary parts from invd and mhatmatrix             */
  /********************************************************************/

  // active, stick and slip part of invd
  Teuchos::RCP<LINALG::SparseMatrix> invda, invdsl, invdst;
  LINALG::SplitMatrix2x2(invd_,gactivedofs_,gidofs,gactivedofs_,gidofs,invda,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gslipdofs_,gstdofs,invdsl,tempmtx1,tempmtx2,tempmtx3);
  LINALG::SplitMatrix2x2(invda,gactivedofs_,gidofs,gstdofs,gslipdofs_,invdst,tempmtx1,tempmtx2,tempmtx3);

  // for Implicit Wear
  Teuchos::RCP<LINALG::SparseMatrix> wa, wi, wsl, wst;          // split lin.matrix for Cgap
  Teuchos::RCP<LINALG::SparseMatrix> wsla, wsli, wslsl, wslst;  // split lin.matrix for Cslip

  // coupling part of dmatrix (only nonzero for 3D quadratic case!)
  Teuchos::RCP<LINALG::SparseMatrix> dai;
  LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,tempmtx1,dai,tempmtx2,tempmtx3);

   // do the multiplication dhat = invda * dai
  Teuchos::RCP<LINALG::SparseMatrix> dhat = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
  if (aset && iset) dhat = LINALG::MLMultiply(*invda,false,*dai,false,false,false,true);
  dhat->Complete(*gidofs,*gactivedofs_);

  // active part of mmatrix
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixa;
  LINALG::SplitMatrix2x2(mmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mmatrixa,tempmtx1,tempmtx2,tempmtx3);

  // do the multiplication mhataam = invda * mmatrixa
  // (this is only different from mhata for 3D quadratic case!)
  Teuchos::RCP<LINALG::SparseMatrix> mhataam = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,10));
  if (aset) mhataam = LINALG::MLMultiply(*invda,false,*mmatrixa,false,false,false,true);
  mhataam->Complete(*gmdofrowmap_,*gactivedofs_);

  // for the case without full linearization, we still need the
  // "classical" active part of mhat, which is isolated here
  Teuchos::RCP<LINALG::SparseMatrix> mhata;
  LINALG::SplitMatrix2x2(mhatmatrix_,gactivedofs_,gidofs,gmdofrowmap_,tempmap,mhata,tempmtx1,tempmtx2,tempmtx3);

  // scaling of invd and dai
  invda->Scale(1/(1-alphaf_));
  invdsl->Scale(1/(1-alphaf_));
  invdst->Scale(1/(1-alphaf_));
  dai->Scale(1-alphaf_);

  /********************************************************************/
  //                       WEARBLOCKS - BASIC
  /********************************************************************/
  // wcoeff: all blocks based on T-matrix have to be scaled with this
  // coefficent...
  double wcoeff = Params().get<double>("WEARCOEFF");

  /********************************************************************/
  /* (a) create inv E                                                */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> inve = Teuchos::rcp(new LINALG::SparseMatrix(*ematrix_));
  Teuchos::RCP<Epetra_Vector> diage = LINALG::CreateVector(*gslipn_,true);
  int erre = 0;

  // extract diagonal of inve into diage
  inve->ExtractDiagonalCopy(*diage);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diage->MyLength();++i)
    if ((*diage)[i]==0.0) (*diage)[i]=1.0;

  // scalar inversion of diagonal values
  erre = diage->Reciprocal(*diage);
  if (erre>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  erre = inve->ReplaceDiagonalValues(*diage);

  /********************************************************************/
  /* (b) build linedis + lintdis                                     */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> linetdis = Teuchos::rcp(new LINALG::SparseMatrix(*lintdis_));

  // scale T-part with wear coeff.
  linetdis->Scale(-wcoeff);
  linetdis->UnComplete();

  // add E-part
  linetdis->Add(*linedis_,false,1.0,1.0);

  //complete
  linetdis->Complete(*gsmdofrowmap_,*gslipn_);

  /********************************************************************/
  /* (c) Split linetdis into 4 matrix blocks                         */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> linte_m, linte_s, dummy1, dummy2;
  Teuchos::RCP<LINALG::SparseMatrix> linte_st, linte_sl, linte_i, linte_a;

  //split slave master...
  LINALG::SplitMatrix2x2(linetdis,gslipn_,gslipn_,gsdofrowmap_,gmdofrowmap_,linte_s,linte_m,dummy1,dummy2);

  //split active inactive
  LINALG::SplitMatrix2x2(linte_s,gslipn_,gslipn_,gactivedofs_,gidofs,linte_a,linte_i,dummy1,dummy2);

  //split stick slip
  LINALG::SplitMatrix2x2(linte_a,gslipn_,gslipn_,gslipdofs_,gstdofs,linte_sl,linte_st,dummy1,dummy2);

  /********************************************************************/
  /* (d) Split lintlm_ into 4 matrix blocks                           */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> lintlm_st, lintlm_sl, lintlm_i, lintlm_a;

  //split active inactive
  lintlm_->Scale(wcoeff);
  LINALG::SplitMatrix2x2(lintlm_,gslipn_,gslipn_,gactivedofs_,gidofs,lintlm_a,lintlm_i,dummy1,dummy2);

  //split stick slip
  LINALG::SplitMatrix2x2(lintlm_a,gslipn_,gslipn_,gslipdofs_,gstdofs,lintlm_sl,lintlm_st,dummy1,dummy2);

  /********************************************************************/
  /* (e) create Tn*D^-T                                               */
  /********************************************************************/
  Teuchos::RCP<LINALG::SparseMatrix> tndt_a;
  if (aset) tndt_a = LINALG::MLMultiply(*lintlm_a,false,*invda,true,false,false,true);

  Teuchos::RCP<LINALG::SparseMatrix> wear_dn, wear_dm, wear_di, wear_da;

  /********************************************************************/
  /* (f) delta_w = wear_dn*delta_dn + wear_dm*delta_dm + wear_da*delta_da + wear_di*delta_di + fw_z                                                  */
  /********************************************************************/
  // these are the matrices for the wear --> wear at slipnodes
  if(slipset)
  {
    // dn
    if (aset)
    {
      wear_dn = LINALG::MLMultiply(*inve,false,*tndt_a,false,false,false,true);
      wear_dn = LINALG::MLMultiply(*wear_dn,false,*kan,false,false,false,true);

      wear_dn->Scale(-1.0);
    }
    // dm
    if (aset)
    {
      wear_dm = LINALG::MLMultiply(*inve,false,*tndt_a,false,false,false,true);
      wear_dm = LINALG::MLMultiply(*wear_dm,false,*kam,false,false,false,true);

      wear_dm->Scale(-1.0);
    }
    // da
    if (aset)
    {
      wear_da = LINALG::MLMultiply(*inve,false,*tndt_a,false,false,false,true);
      wear_da = LINALG::MLMultiply(*wear_da,false,*kaa,false,false,false,false);

      Teuchos::RCP<LINALG::SparseMatrix> inve_te;
      inve_te = LINALG::MLMultiply(*inve,false,*linte_a,false,false,false,true);

      wear_da->Add(*inve_te,false,1.0,1.0);
      wear_da->Complete(*gactivedofs_,*gslipn_);
      wear_da->Scale(-1.0);
    }
    // di
    if (aset && iset)
    {
      wear_di = LINALG::MLMultiply(*inve,false,*tndt_a,false,false,false,true);
      wear_di = LINALG::MLMultiply(*wear_di,false,*kai,false,false,false,false);

      Teuchos::RCP<LINALG::SparseMatrix> inve_te;
      inve_te = LINALG::MLMultiply(*inve,false,*linte_i,false,false,false,true);

      wear_di->Add(*inve_te,false,1.0,1.0);
      wear_di->Complete(*gidofs,*gslipn_);
      wear_di->Scale(-1.0);
    }
  }

  // store to header_variables --> reused for recovering
  dnblock_ = wear_dn;
  dmblock_ = wear_dm;
  diblock_ = wear_di;
  dablock_ = wear_da;

  //*************************************************************
  // (g) prepare f_w
  //*************************************************************
  Teuchos::RCP<Epetra_Vector> dvec = Teuchos::rcp(new Epetra_Vector(*gslipn_));
  Teuchos::RCP<Epetra_Vector> fw_z = Teuchos::rcp(new Epetra_Vector(*gslipn_));
  Teuchos::RCP<Epetra_Vector> fw_wrhs = Teuchos::rcp(new Epetra_Vector(*gslipn_));

  if (slipset)
  {
    tndt_a->Multiply(false,*fa,*dvec);
    inve->Multiply(false,*dvec,*fw_z);

    // fa = -ra
    fw_z->Scale(-1.0);

    // this wear condrhs excludes the lm-part of the rhs
    WearCondRhs_->Scale(-1.0);
    inve->Multiply(false,*WearCondRhs_,*fw_wrhs);
    fw_z->Update(1.0,*fw_wrhs,1.0);

    fw_z->Scale(-1.0);

    fw_ = Teuchos::rcp(new Epetra_Vector(*gslipn_));
    fw_ = fw_z;
  }

  // ************************************************************
  // (h) prepare smatrixW and linslipW_
  // ************************************************************
  Teuchos::RCP<LINALG::SparseMatrix> smatrixW_sl,smatrixW_i,d2,d3;
  LINALG::SplitMatrix2x2(smatrixW_,gactiven_,gactiven_,gslipn_,gwinact_,smatrixW_sl,smatrixW_i,d2,d3);

  Teuchos::RCP<LINALG::SparseMatrix> linslipW_sl,linslipW_i,d5,d6;
  LINALG::SplitMatrix2x2(linslipW_,gslipt_,gslipt_,gslipn_,gwinact_,linslipW_sl,linslipW_i,d5,d6);

  /********************************************************************/
  /* (7) Build the final K blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // knn: nothing to do

  // knm: nothing to do

  // kns: nothing to do

  //-------------------------------------------------------- SECOND LINE
  // kmn: add T(mhataam)*kan
  Teuchos::RCP<LINALG::SparseMatrix> kmnmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  kmnmod->Add(*kmn,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhataam,true,*kan,false,false,false,true);
  kmnmod->Add(*kmnadd,false,1.0,1.0);
  kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());

  // kmm: add T(mhataam)*kam
  Teuchos::RCP<LINALG::SparseMatrix> kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
  kmmmod->Add(*kmm,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhataam,true,*kam,false,false,false,true);
  kmmmod->Add(*kmmadd,false,1.0,1.0);
  kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());

  // kmi: add T(mhataam)*kai
  Teuchos::RCP<LINALG::SparseMatrix> kmimod;
  if (iset)
  {
    kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmimod->Add(*kmi,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhataam,true,*kai,false,false,false,true);
    kmimod->Add(*kmiadd,false,1.0,1.0);
    kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
  }

  // kma: add T(mhataam)*kaa
  Teuchos::RCP<LINALG::SparseMatrix> kmamod;
  if (aset)
  {
    kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofrowmap_,100));
    kmamod->Add(*kma,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhataam,true,*kaa,false,false,false,true);
    kmamod->Add(*kmaadd,false,1.0,1.0);
    kmamod->Complete(kma->DomainMap(),kma->RowMap());
  }

  //--------------------------------------------------------- THIRD LINE
  //------------------- FOR 3D QUADRATIC CASE ----------------------------
  // kin: subtract T(dhat)*kan
  Teuchos::RCP<LINALG::SparseMatrix> kinmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
  kinmod->Add(*kin,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kinadd = LINALG::MLMultiply(*dhat,true,*kan,false,false,false,true);
  kinmod->Add(*kinadd,false,-1.0,1.0);
  kinmod->Complete(kin->DomainMap(),kin->RowMap());

  // kim: subtract T(dhat)*kam
  Teuchos::RCP<LINALG::SparseMatrix> kimmod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
  kimmod->Add(*kim,false,1.0,1.0);
  Teuchos::RCP<LINALG::SparseMatrix> kimadd = LINALG::MLMultiply(*dhat,true,*kam,false,false,false,true);
  kimmod->Add(*kimadd,false,-1.0,1.0);
  kimmod->Complete(kim->DomainMap(),kim->RowMap());

  // kii: subtract T(dhat)*kai
  Teuchos::RCP<LINALG::SparseMatrix> kiimod;
  if (iset)
  {
    kiimod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kiimod->Add(*kii,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kiiadd = LINALG::MLMultiply(*dhat,true,*kai,false,false,false,true);
    kiimod->Add(*kiiadd,false,-1.0,1.0);
    kiimod->Complete(kii->DomainMap(),kii->RowMap());
  }

  // kia: subtract T(dhat)*kaa
  Teuchos::RCP<LINALG::SparseMatrix> kiamod;
  if (iset && aset)
  {
    kiamod = Teuchos::rcp(new LINALG::SparseMatrix(*gidofs,100));
    kiamod->Add(*kia,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kiaadd = LINALG::MLMultiply(*dhat,true,*kaa,false,false,false,true);
    kiamod->Add(*kiaadd,false,-1.0,1.0);
    kiamod->Complete(kia->DomainMap(),kia->RowMap());
  }

  //--------------------------------------------------------- FOURTH LINE
  // WEAR STUFF
  Teuchos::RCP<LINALG::SparseMatrix> kgnw;
  Teuchos::RCP<LINALG::SparseMatrix> kgmw;
  Teuchos::RCP<LINALG::SparseMatrix> kgiw;
  Teuchos::RCP<LINALG::SparseMatrix> kgaw;

  if (aset && slipset) kgnw = LINALG::MLMultiply(*smatrixW_sl,false,*wear_dn,false,false,false,true);
  if (aset && slipset) kgmw = LINALG::MLMultiply(*smatrixW_sl,false,*wear_dm,false,false,false,true);
  if (aset && slipset && iset) kgiw = LINALG::MLMultiply(*smatrixW_sl,false,*wear_di,false,false,false,true);
  if (aset && slipset) kgaw = LINALG::MLMultiply(*smatrixW_sl,false,*wear_da,false,false,false,true);

  //--------------------------------------------------------- FIFTH  LINE
  // blocks for complementary conditions (stick nodes)

  // kstn: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstnmod;
  if (stickset)
  {
    kstnmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstnmod = LINALG::MLMultiply(*kstnmod,false,*kan,false,false,false,true);
  }

  // kstm: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstmmod;
  if(stickset)
  {
    kstmmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstmmod = LINALG::MLMultiply(*kstmmod,false,*kam,false,false,false,true);
  }

  // ksti: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstimod;
  if(stickset && iset)
  {
    kstimod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstimod = LINALG::MLMultiply(*kstimod,false,*kai,false,false,false,true);
  }

  // kstsl: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kstslmod;
  if(stickset && slipset)
  {
    kstslmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kstslmod = LINALG::MLMultiply(*kstslmod,false,*kasl,false,false,false,true);
  }

  // kststmod: multiply with linstickLM
  Teuchos::RCP<LINALG::SparseMatrix> kststmod;
  if (stickset)
  {
    kststmod = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    kststmod = LINALG::MLMultiply(*kststmod,false,*kast,false,false,false,true);
  }

  //--------------------------------------------------------- SIXTH LINE
  // blocks for complementary conditions (slip nodes)

  // ksln: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslnmod;
  if(slipset)
  {
    kslnmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslnmod = LINALG::MLMultiply(*kslnmod,false,*kan,false,false,false,true);
  }

  // kslm: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslmmod;
  if(slipset)
  {
    kslmmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslmmod = LINALG::MLMultiply(*kslmmod,false,*kam,false,false,false,true);
  }

  // ksli: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslimod;
  if (slipset && iset)
  {
    kslimod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslimod = LINALG::MLMultiply(*kslimod,false,*kai,false,false,false,true);
  }

  // kslsl: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslslmod;
  if(slipset)
  {
    kslslmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslslmod = LINALG::MLMultiply(*kslslmod,false,*kasl,false,false,false,true);
  }

  // slstmod: multiply with linslipLM
  Teuchos::RCP<LINALG::SparseMatrix> kslstmod;
  if (slipset && stickset)
  {
    kslstmod = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    kslstmod = LINALG::MLMultiply(*kslstmod,false,*kast,false,false,false,true);
  }

  // WEARSTUFF FOR THIS LINE
  Teuchos::RCP<LINALG::SparseMatrix> kslnw;
  Teuchos::RCP<LINALG::SparseMatrix> kslmw;
  Teuchos::RCP<LINALG::SparseMatrix> ksliw;
  Teuchos::RCP<LINALG::SparseMatrix> kslaw;

  if (slipset) kslnw = LINALG::MLMultiply(*linslipW_sl,false,*wear_dn,false,false,false,true);
  if (slipset) kslmw = LINALG::MLMultiply(*linslipW_sl,false,*wear_dm,false,false,false,true);
  if (slipset && iset) ksliw = LINALG::MLMultiply(*linslipW_sl,false,*wear_di,false,false,false,true);
  if (slipset) kslaw = LINALG::MLMultiply(*linslipW_sl,false,*wear_da,false,false,false,true);

  /********************************************************************/
  /* (8) Build the final f blocks                                     */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // fn: nothing to do

  //---------------------------------------------------------- SECOND LINE
  // fm: add alphaf * old contact forces (t_n)
  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Mold^T * zold to fit
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    Teuchos::RCP<Epetra_Vector> tempvecm2  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(mold_->RowMap()));
    if (mold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    mold_->Multiply(true,*zoldexp,*tempvecm2);
    if (mset) LINALG::Export(*tempvecm2,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    Teuchos::RCP<Epetra_Vector> tempvecm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->Multiply(true,*zold_,*tempvecm);
    fm->Update(alphaf_,*tempvecm,1.0);
  }

  // fs: prepare alphaf * old contact forces (t_n)
  Teuchos::RCP<Epetra_Vector> fsadd = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // for self contact, slave and master sets may have changed,
  // thus we have to export the product Dold^T * zold to fit
  if (IsSelfContact())
  {
    Teuchos::RCP<Epetra_Vector> tempvec  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    dold_->Multiply(true,*zoldexp,*tempvec);
    if (sset) LINALG::Export(*tempvec,*fsadd);
  }
  // if there is no self contact everything is ok
  else
  {
    dold_->Multiply(true,*zold_,*fsadd);
  }

  // fa: subtract alphaf * old contact forces (t_n)
  if (aset)
  {
    Teuchos::RCP<Epetra_Vector> faadd = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    LINALG::Export(*fsadd,*faadd);
    fa->Update(-alphaf_,*faadd,1.0);
  }

  // fm: add T(mhat)*fa
  Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (aset) mhataam->Multiply(true,*fa,*fmmod);
  fmmod->Update(1.0,*fm,1.0);

  //--------------------------------------------------------- THIRD LINE
  // fi: subtract alphaf * old contact forces (t_n)
  if (iset)
  {
    Teuchos::RCP<Epetra_Vector> fiadd = Teuchos::rcp(new Epetra_Vector(*gidofs));
    LINALG::Export(*fsadd,*fiadd);
    fi->Update(-alphaf_,*fiadd,1.0);
  }

  // fi: add T(dhat)*fa
  Teuchos::RCP<Epetra_Vector> fimod = Teuchos::rcp(new Epetra_Vector(*gidofs));
  if (aset) dhat->Multiply(true,*fa,*fimod);
  fimod->Update(1.0,*fi,-1.0);

  //--------------------------------------------------------- FOURTH LINE
  //
  Teuchos::RCP<Epetra_Vector> fw_g = Teuchos::rcp(new Epetra_Vector(*gactiven_));
  Teuchos::RCP<Epetra_Vector> fwi_g = Teuchos::rcp(new Epetra_Vector(*gactiven_));

  if(aset && slipset) smatrixW_sl->Multiply(false,*fw_z,*fw_g);

  if(aset && (iset || stickset)) smatrixW_i->Multiply(false,*inactiveWearRhs_,*fwi_g);

  //--------------------------------------------------------- FIFTH LINE
  Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);  // get global stick dofs

  // split the lagrange multiplier vector in stick and slip part
  Teuchos::RCP<Epetra_Vector> za = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> zi = Teuchos::rcp(new Epetra_Vector(*gidofs));
  Teuchos::RCP<Epetra_Vector> zst = Teuchos::rcp(new Epetra_Vector(*gstickdofs));
  Teuchos::RCP<Epetra_Vector> zsl = Teuchos::rcp(new Epetra_Vector(*gslipdofs_));

  LINALG::SplitVector(*gsdofrowmap_, *z_, gactivedofs_, za, gidofs, zi);
  LINALG::SplitVector(*gactivedofs_, *za, gstickdofs, zst, gslipdofs_, zsl);
  Teuchos::RCP<Epetra_Vector> tempvec1;

  // fst: mutliply with linstickLM
  Teuchos::RCP<Epetra_Vector> fstmod;
  if (stickset)
  {
    fstmod = Teuchos::rcp(new Epetra_Vector(*gstickt));
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*linstickLM_,false,*invdst,true,false,false,true);
    temp1->Multiply(false,*fa,*fstmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gstickt));

    linstickLM_->Multiply(false, *zst, *tempvec1);
    fstmod->Update(-1.0,*tempvec1,1.0);
  }

  //--------------------------------------------------------- SIXTH LINE
  // fsl: mutliply with linslipLM
  Teuchos::RCP<Epetra_Vector> fslmod;

  if (slipset)
  {
    fslmod = Teuchos::rcp(new Epetra_Vector(*gslipt_));
    Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*linslipLM_,false,*invdsl,true,false,false,true);
    temp->Multiply(false,*fa,*fslmod);

    tempvec1 = Teuchos::rcp(new Epetra_Vector(*gslipt_));

    linslipLM_->Multiply(false, *zsl, *tempvec1);

    fslmod->Update(-1.0,*tempvec1,1.0);
  }

  //WEAR
  Teuchos::RCP<Epetra_Vector> fw_sl = Teuchos::rcp(new Epetra_Vector(*gslipt_));
  if (slipset) linslipW_sl->Multiply(false,*fw_z,*fw_sl);

  Teuchos::RCP<Epetra_Vector> fwi_sl = Teuchos::rcp(new Epetra_Vector(*gslipt_));
  if (slipset && (iset || stickset)) linslipW_i->Multiply(false,*inactiveWearRhs_,*fwi_sl);

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
    //----------------------------------------------------------- FIRST LINE
    // nothing to do (ndof-map independent of redistribution)

    //---------------------------------------------------------- SECOND LINE
    kmnmod = MORTAR::MatrixRowTransform(kmnmod,pgmdofrowmap_);
    kmmmod = MORTAR::MatrixRowTransform(kmmmod,pgmdofrowmap_);
    if (iset) kmimod = MORTAR::MatrixRowTransform(kmimod,pgmdofrowmap_);
    if (aset) kmamod = MORTAR::MatrixRowTransform(kmamod,pgmdofrowmap_);

    //----------------------------------------------------------- THIRD LINE
    if (iset)
    {
      kinmod = MORTAR::MatrixRowTransform(kinmod,pgsdofrowmap_);
      kimmod = MORTAR::MatrixRowTransform(kimmod,pgsdofrowmap_);
      kiimod = MORTAR::MatrixRowTransform(kiimod,pgsdofrowmap_);
      if (aset) kiamod = MORTAR::MatrixRowTransform(kiamod,pgsdofrowmap_);
    }

    //---------------------------------------------------------- FOURTH LINE
    if (aset)
    {
      smatrix_ = MORTAR::MatrixRowTransform(smatrix_,pgsdofrowmap_);

      //WEAR
      if (slipset) kgnw = MORTAR::MatrixRowTransform(kgnw,pgsdofrowmap_);
      if (slipset) kgmw = MORTAR::MatrixRowTransform(kgmw,pgsdofrowmap_);
      if (slipset && iset) kgiw = MORTAR::MatrixRowTransform(kgiw,pgsdofrowmap_);
      if (slipset) kgaw = MORTAR::MatrixRowTransform(kgaw,pgsdofrowmap_);
    }

    //----------------------------------------------------------- FIFTH LINE
    if (stickset)
    {
      kstnmod = MORTAR::MatrixRowTransform(kstnmod,pgsdofrowmap_);
      kstmmod = MORTAR::MatrixRowTransform(kstmmod,pgsdofrowmap_);
      if (iset) kstimod = MORTAR::MatrixRowTransform(kstimod,pgsdofrowmap_);
      if (slipset) kstslmod = MORTAR::MatrixRowTransform(kstslmod,pgsdofrowmap_);
      kststmod = MORTAR::MatrixRowTransform(kststmod,pgsdofrowmap_);
      linstickDIS_ = MORTAR::MatrixRowTransform(linstickDIS_,pgsdofrowmap_);
    }

    //----------------------------------------------------------- SIXTH LINE
    if (slipset)
    {
      kslnmod = MORTAR::MatrixRowTransform(kslnmod,pgsdofrowmap_);
      kslmmod = MORTAR::MatrixRowTransform(kslmmod,pgsdofrowmap_);
      if (iset) kslimod = MORTAR::MatrixRowTransform(kslimod,pgsdofrowmap_);
      if (stickset) kslstmod = MORTAR::MatrixRowTransform(kslstmod,pgsdofrowmap_);
      kslslmod = MORTAR::MatrixRowTransform(kslslmod,pgsdofrowmap_);
      linslipDIS_ = MORTAR::MatrixRowTransform(linslipDIS_,pgsdofrowmap_);

      //WEAR
      kslnw = MORTAR::MatrixRowTransform(kslnw,pgsdofrowmap_);
      kslmw = MORTAR::MatrixRowTransform(kslmw,pgsdofrowmap_);
      if (iset) ksliw = MORTAR::MatrixRowTransform(ksliw,pgsdofrowmap_);
      kslaw = MORTAR::MatrixRowTransform(kslaw,pgsdofrowmap_);
    }
  }

  /********************************************************************/
  /* (10) Global setup of kteffnew (including contact)                */
  /********************************************************************/

  Teuchos::RCP<LINALG::SparseMatrix> kteffnew = Teuchos::rcp(new LINALG::SparseMatrix(*ProblemDofs(),81,true,false,kteffmatrix->GetMatrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*ProblemDofs());

  //--------------------------------------------------------- FIRST LINE
  // add n submatrices to kteffnew
  kteffnew->Add(*knn,false,1.0,1.0);
  kteffnew->Add(*knm,false,1.0,1.0);
  if (sset) kteffnew->Add(*kns,false,1.0,1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m submatrices to kteffnew
  kteffnew->Add(*kmnmod,false,1.0,1.0);
  kteffnew->Add(*kmmmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kmimod,false,1.0,1.0);
  if (aset) kteffnew->Add(*kmamod,false,1.0,1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i submatrices to kteffnew
  if (iset) kteffnew->Add(*kinmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kimmod,false,1.0,1.0);
  if (iset) kteffnew->Add(*kiimod,false,1.0,1.0);
  if (iset && aset) kteffnew->Add(*kiamod,false,1.0,1.0);

  //-------------------------------------------------------- FOURTH LINE
  // WEAR - STUFF
  // add a submatrices to kteffnew
  if (aset) kteffnew->Add(*smatrix_,false,1.0,1.0);
  if (aset && slipset) kteffnew->Add(*kgnw,false,1.0,1.0);
  if (aset && slipset) kteffnew->Add(*kgmw,false,1.0,1.0);
  if (aset && slipset && iset) kteffnew->Add(*kgiw,false,1.0,1.0);
  if (aset && slipset) kteffnew->Add(*kgaw,false,1.0,1.0);

  //--------------------------------------------------------- FIFTH LINE
  // add st submatrices to kteffnew
  if (stickset) kteffnew->Add(*kstnmod,false,1.0,1.0);
  if (stickset) kteffnew->Add(*kstmmod,false,1.0,1.0);
  if (stickset && iset) kteffnew->Add(*kstimod,false,1.0,1.0);
  if (stickset && slipset) kteffnew->Add(*kstslmod,false,1.0,1.0);
  if (stickset) kteffnew->Add(*kststmod,false,1.0,1.0);

  // add terms of linearization of sick condition to kteffnew
  if (stickset) kteffnew->Add(*linstickDIS_,false,-1.0,1.0);

  //--------------------------------------------------------- SIXTH LINE
  // add sl submatrices to kteffnew
  if (slipset) kteffnew->Add(*kslnmod,false,1.0,1.0);
  if (slipset) kteffnew->Add(*kslmmod,false,1.0,1.0);
  if (slipset && iset) kteffnew->Add(*kslimod,false,1.0,1.0);
  if (slipset) kteffnew->Add(*kslslmod,false,1.0,1.0);
  if (slipset && stickset) kteffnew->Add(*kslstmod,false,1.0,1.0);

  // add terms of linearization of slip condition to kteffnew and feffnew
  if (slipset) kteffnew->Add(*linslipDIS_,false,-1.0,+1.0);

  // WEAR - STUFF
  if (slipset) kteffnew->Add(*kslnw,false,-1.0,1.0);
  if (slipset) kteffnew->Add(*kslmw,false,-1.0,1.0);
  if (slipset && (iset)) kteffnew->Add(*ksliw,false,-1.0,1.0);
  if (slipset) kteffnew->Add(*kslaw,false,-1.0,1.0);

  // FillComplete kteffnew (square)
  kteffnew->Complete();

  /********************************************************************/
  /* (11) Global setup of feffnew (including contact)                 */
  /********************************************************************/

  //--------------------------------------------------------- FIRST LINE
  // add n subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fn,*fnexp);
  feffnew->Update(1.0,*fnexp,1.0);

  //-------------------------------------------------------- SECOND LINE
  // add m subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  LINALG::Export(*fmmod,*fmmodexp);
  feffnew->Update(1.0,*fmmodexp,1.0);

  //--------------------------------------------------------- THIRD LINE
  // add i subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fimodexp;
  if (iset)
  {
    fimodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fimod,*fimodexp);
    feffnew->Update(1.0,*fimodexp,1.0);
  }

  //-------------------------------------------------------- FOURTH LINE
  // add weighted gap vector to feffnew, if existing
  Teuchos::RCP<Epetra_Vector> gexp;
  Teuchos::RCP<Epetra_Vector> fwexp;
  Teuchos::RCP<Epetra_Vector> fwiexp;
  Teuchos::RCP<Epetra_Vector> fgmodexp;

  if (aset)
  {
    gexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*gact,*gexp);
    feffnew->Update(-1.0,*gexp,1.0);

    //DUE TO WEAR
    fwexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fw_g,*fwexp);
    feffnew->Update(-1.0,*fwexp,1.0);

    // DUE TO WEAR
    if(iset || stickset)
    {
      fwiexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*fwi_g,*fwiexp);
      feffnew->Update(+1.0,*fwiexp,1.0);
    }
  }

  //--------------------------------------------------------- FIFTH LINE
  // add st subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fstmodexp;
  if (stickset)
  {
    fstmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fstmod,*fstmodexp);
    feffnew->Update(1.0,*fstmodexp,+1.0);
  }

  // add terms of linearization feffnew
   if (stickset)
   {
      Teuchos::RCP<Epetra_Vector> linstickRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
      LINALG::Export(*linstickRHS_,*linstickRHSexp);
      feffnew->Update(-1.0,*linstickRHSexp,1.0);
   }

  //--------------------------------------------------------- SIXTH LINE

  // add a subvector to feffnew
  Teuchos::RCP<Epetra_Vector> fslmodexp;
  Teuchos::RCP<Epetra_Vector> fwslexp;
  Teuchos::RCP<Epetra_Vector> fwsliexp;
  Teuchos::RCP<Epetra_Vector> fslwmodexp;

  if (slipset)
  {
    fslmodexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fslmod,*fslmodexp);
    feffnew->Update(1.0,*fslmodexp,1.0);
  }

  if (slipset)
  {
    Teuchos::RCP<Epetra_Vector> linslipRHSexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*linslipRHS_,*linslipRHSexp);
    feffnew->Update(-1.0,*linslipRHSexp,1.0);
  }

  // DUE TO WEAR
  if (slipset)
  {
    fwslexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fw_sl,*fwslexp);
    feffnew->Update(+1.0,*fwslexp,1.0);
  }
  if (slipset && (iset || stickset))
  {
    fwsliexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fwi_sl,*fwsliexp);
    feffnew->Update(-1.0,*fwsliexp,1.0);
  }

  // finally do the replacement
  kteff = kteffnew;
  feff  = feffnew;

  return;
}


/*----------------------------------------------------------------------*
 | evaluate frictional wear contact (public)                 farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                     Teuchos::RCP<Epetra_Vector>& feff)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  kteff->Complete();

  // systemtype
  INPAR::CONTACT::SystemType systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

  // get wear shapefunction type
  INPAR::CONTACT::WearShape wearshapefcn =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearShape>(Params(),"WEAR_SHAPEFCN");

  /**********************************************************************/
  /* export weighted gap vector to gactiveN-map                         */
  /* Here, the additional wear is already included !!!                  */
  /**********************************************************************/
  Teuchos::RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
  if (gact->GlobalLength())
  {
    LINALG::Export(*g_,*gact);
    gact->ReplaceMap(*gactiven_);
  }

  /**********************************************************************/
  /* build global matrix t with tangent vectors of active nodes         */
  /* and global matrix s with normal derivatives of active nodes        */
  /* and global matrix linstick with derivatives of stick nodes         */
  /* and global matrix linslip with derivatives of slip nodes           */
  /* and inactive right-hand side with old lagrange multipliers (incr)  */
  /**********************************************************************/
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->AssembleT(*tmatrix_);
    interface_[i]->AssembleS(*smatrix_);
    interface_[i]->AssembleLinDM(*lindmatrix_,*linmmatrix_);
    interface_[i]->AssembleLinStick(*linstickLM_,*linstickDIS_,*linstickRHS_);
    interface_[i]->AssembleLinSlip(*linslipLM_,*linslipDIS_,*linslipRHS_);
    if (systype != INPAR::CONTACT::system_condensed)
      interface_[i]->AssembleInactiverhs(*inactiverhs_);

    //***************************************************
    // Assemble lin. for implicit wear algorithm
    if (wearimpl_)
    {    // assemble wear-specific matrices
      interface_[i]->AssembleLinWLm(*wlinmatrix_);
      interface_[i]->AssembleLinWLmSl(*wlinmatrixsl_);

#ifdef CONSISTENTSTICK
    interface_[i]->AssembleLinWLmSt(*wlinmatrixst_);
#endif
    }

    //***************************************************
    // Assemble all Entries for discrete wear algorithm
    if (weardiscr_)
    {
      // blocks for w-lines
      interface_[i]->AssembleTE(*twmatrix_,*ematrix_);
      interface_[i]->AssembleLinT_D(*lintdis_);
      interface_[i]->AssembleLinT_LM(*lintlm_);
      interface_[i]->AssembleLinE_D(*linedis_);

      // blocks for z-lines
      interface_[i]->AssembleLinG_W(*smatrixW_);
      interface_[i]->AssembleLinSlip_W(*linslipW_);

      // w-line rhs
      interface_[i]->AssembleInactiveWearRhs(*inactiveWearRhs_);
      interface_[i]->AssembleWearCondRhs(*WearCondRhs_);

      // for both-sided discrete wear
      if(wearbothdiscr_)
      {
        // blocks for w-lines
        interface_[i]->AssembleTE_Master(*twmatrixM_,*ematrixM_);
        interface_[i]->AssembleLinT_D_Master(*lintdisM_);
        interface_[i]->AssembleLinT_LM_Master(*lintlmM_);
        interface_[i]->AssembleLinE_D_Master(*linedisM_);

        // w-line rhs
        interface_[i]->AssembleInactiveWearRhs_Master(*inactiveWearRhsM_);
        interface_[i]->AssembleWearCondRhs_Master(*WearCondRhsM_);
      }
    }
  }// end interface loop

  /**********************************************************************/
  /* Complete matrices                                                  */
  /**********************************************************************/
  // FillComplete() global matrix T
  tmatrix_->Complete(*gactivedofs_,*gactivet_);

  // FillComplete() global matrix S
  smatrix_->Complete(*gsmdofrowmap_,*gactiven_);

  // FillComplete() global matrices LinD, LinM
  // (again for linD gsdofrowmap_ is sufficient as domain map,
  // but in the edge node modification case, master entries occur!)
  lindmatrix_->Complete(*gsmdofrowmap_,*gsdofrowmap_);
  linmmatrix_->Complete(*gsmdofrowmap_,*gmdofrowmap_);

  // FillComplete global Matrix linstickLM_, linstickDIS_
  Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);
  Teuchos::RCP<Epetra_Map> gstickdofs = LINALG::SplitMap(*gactivedofs_,*gslipdofs_);
  linstickLM_->Complete(*gstickdofs,*gstickt);
  linstickDIS_->Complete(*gsmdofrowmap_,*gstickt);

  // FillComplete global Matrix linslipLM_ and linslipDIS_
  linslipLM_->Complete(*gslipdofs_,*gslipt_);
  linslipDIS_->Complete(*gsmdofrowmap_,*gslipt_);

  /**********************************************************************/
  /* Complete impl wear-specific matrices                               */
  /**********************************************************************/
  if (wearimpl_)
  {
    // complete wear-specific matrices
    wlinmatrix_->Complete(*gsdofrowmap_,*gactiven_);
    wlinmatrixsl_->Complete(*gsdofrowmap_,*gslipt_);

#ifdef CONSISTENTSTICK
  wlinmatrixst_->Complete(*gsdofrowmap_,*gstickt);
#endif
  }

  /**********************************************************************/
  /* Complete discr wear-specific matrices                              */
  /**********************************************************************/
  if(weardiscr_)
  {
    if(sswear_)
    {
      if(wearshapefcn == INPAR::CONTACT::wear_shape_dual)
        ematrix_->Complete(*gactiven_,*gactiven_); // quadr. matrix --> for dual shapes --> diag.
      else if (wearshapefcn == INPAR::CONTACT::wear_shape_standard)
        ematrix_->Complete(*gsdofnrowmap_,*gactiven_); // quadr. matrix --> for dual shapes --> diag.
      else
        dserror("chosen shape fnc for wear not supported!");

      twmatrix_->Complete(*gsdofnrowmap_,*gactiven_);

      lintdis_->Complete(*gsmdofrowmap_,*gactiven_);
      lintlm_->Complete(*gsdofrowmap_,*gactiven_);
      linedis_->Complete(*gsmdofrowmap_,*gactiven_);
    }
    else
    {
      if(wearshapefcn == INPAR::CONTACT::wear_shape_dual)
        ematrix_->Complete(*gslipn_,*gslipn_); // quadr. matrix --> for dual shapes --> diag.
      else if (wearshapefcn == INPAR::CONTACT::wear_shape_standard)
        ematrix_->Complete(*gsdofnrowmap_,*gslipn_); // quadr. matrix --> for dual shapes --> diag.
      else
        dserror("ERROR: chosen shape fnc for wear not supported!");

      twmatrix_->Complete(*gsdofnrowmap_,*gslipn_);

      lintdis_->Complete(*gsmdofrowmap_,*gslipn_);
      lintlm_->Complete(*gsdofrowmap_,*gslipn_);
      linedis_->Complete(*gsmdofrowmap_,*gslipn_);
    }

    // if both-sided complete with all wear dofs!
    if(wearbothdiscr_)
    {
      smatrixW_->Complete(*galldofnrowmap_,*gactiven_);
      linslipW_->Complete(*galldofnrowmap_,*gslipt_);
    }
    else
    {
      smatrixW_->Complete(*gsdofnrowmap_,*gactiven_);
      linslipW_->Complete(*gsdofnrowmap_,*gslipt_);
    }

    if(wearbothdiscr_)
    {
      ematrixM_->Complete(*gmdofnrowmap_,*gmslipn_);
      twmatrixM_->Complete(*gsdofnrowmap_,*gmslipn_);

      lintdisM_->Complete(*gsmdofrowmap_,*gmslipn_);
      lintlmM_->Complete(*gsdofrowmap_,*gmslipn_);
      linedisM_->Complete(*gsmdofrowmap_,*gmslipn_);
    }
  }

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // LinD      ---->   T^(-T) * LinD
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify lindmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp1 = LINALG::MLMultiply(*invtrafo_,true,*lindmatrix_,false,false,false,true);
    lindmatrix_   = temp1;
  }

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL) ------ WEARDISCR
  // HERE THE LM ARE SOLVED ABSOLUTELY !!!
  //**********************************************************************
  //**********************************************************************
  if ((systype == INPAR::CONTACT::system_condensed) && weardiscr_==true )
  {
    CondenseWearDiscr(kteff,feff,gact);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE B: CONDENSED SYSTEM (DUAL) ------ WEARIMPL
  // HERE THE LM ARE SOLVED ABSOLUTELY !!!
  //**********************************************************************
  //**********************************************************************
  else if ((systype == INPAR::CONTACT::system_condensed) && weardiscr_==false )
  {
    CondenseWearImplExpl(kteff,feff,gact);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE C: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    PrepareSaddlePointSystem(kteff,feff);
  }
  // FD checks...
#ifdef WEARIMPLICITFDLM
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (!wearimpl_)
      dserror("Explicit wear algorithm: no FD check necessary!");

    interface_[i]->FDCheckWearDerivLm();
  }
#endif

#ifdef WEARIMPLICITFD
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    if (!wearimpl_)
      dserror("Explicit wear algorithm: no FD check necessary!");

    interface_[i]->FDCheckWearDeriv();
  }
#endif

#ifdef CONTACTFDGAP
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckGapDeriv();
    interface_[i]->FDCheckGapDeriv_W();
  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSLIPINCR
  // FD check of weighted gap g derivatives (non-penetr. condition)
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->FDCheckSlipIncrDerivTXI();
    if (Dim()==3)
      interface_[i]->FDCheckSlipIncrDerivTETA();

  }
#endif // #ifdef CONTACTFDGAP

#ifdef CONTACTFDSTICK

  if (gstickt->NumGlobalElements())
  {
    // FD check of stick condition
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckStickDeriv(*linstickLM_,*linstickDIS_);
    }
  }
#endif // #ifdef CONTACTFDSTICK

#ifdef CONTACTFDT_D

  if (!weardiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->FDCheckDerivT_D(*lintdis_);

#endif

#ifdef CONTACTFDT_D_MASTER

  if (!weardiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothdiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->FDCheckDerivT_D_Master(*lintdisM_);

#endif

#ifdef CONTACTFDE_D

  if (!weardiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->FDCheckDerivE_D(*linedis_);

#endif

#ifdef CONTACTFDE_D_MASTER

  if (!weardiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  if (!wearbothdiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

  // FD check of stick condition
  for (int i=0; i<(int)interface_.size(); ++i)
    interface_[i]->FDCheckDerivE_D_Master(*linedisM_);

#endif

#ifdef CONTACTFDSLIP

  if (gslipnodes_->NumGlobalElements())
  {
    // FD check of slip condition
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      interface_[i]->FDCheckSlipDeriv(*linslipLM_,*linslipDIS_,*linslipW_);
    }
  }
#endif // #ifdef CONTACTFDSLIP

#ifdef CONTACTFDMORTART

  if (!weardiscr_)
    dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTART -----------------------------------" << std::endl;
    twmatrix_->Complete();
    if( twmatrix_->NormOne() )
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->FDCheckMortarTDeriv();
    //twmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTART -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTARE

    if (!weardiscr_)
      dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARE -----------------------------------" << std::endl;
    ematrix_->Complete();
    if( ematrix_->NormOne() )
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->FDCheckMortarEDeriv();
    //ematrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARE -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTARE_MASTER

    if (!weardiscr_)
      dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    if (!wearbothdiscr_)
      dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARE_MASTER -----------------------------------" << std::endl;
    ematrixM_->Complete();
    if( ematrixM_->NormOne() )
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->FDCheckMortarE_Master_Deriv();
    //ematrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARE_MASTER -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD

#ifdef CONTACTFDMORTART_MASTER

    if (!weardiscr_)
      dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    if (!wearbothdiscr_)
      dserror("FD CHECK ONLY FOR DISCRETE WEAR!");

    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTART_MASTER -----------------------------------" << std::endl;
    twmatrixM_->Complete();
    if( twmatrixM_->NormOne() )
      for (int i=0; i<(int)interface_.size(); ++i)
        interface_[i]->FDCheckMortarT_Master_Deriv();
    //ematrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTART_MASTER -----------------------------------" << std::endl;
#endif // #ifdef CONTACTFDMORTARD

  return;
}

/*----------------------------------------------------------------------*
 | preparation for self-contact and assemble lind/m          farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::PrepareSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                                                   Teuchos::RCP<Epetra_Vector>& feff)
{
  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp2 = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrix_    = temp2;
  }

  // transform if necessary
  if (ParRedist())
  {
    lindmatrix_ = MORTAR::MatrixRowTransform(lindmatrix_,pgsdofrowmap_);
    linmmatrix_ = MORTAR::MatrixRowTransform(linmmatrix_,pgmdofrowmap_);
  }

  // add contact stiffness
  kteff->UnComplete();
  kteff->Add(*lindmatrix_,false,1.0-alphaf_,1.0);
  kteff->Add(*linmmatrix_,false,1.0-alphaf_,1.0);
  kteff->Complete();

  // for self contact, slave and master sets may have changed,
  // thus we have to export the products Dold^T * zold / D^T * z to fit
  // thus we have to export the products Mold^T * zold / M^T * z to fit
  if (IsSelfContact())
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecd  = Teuchos::rcp(new Epetra_Vector(dmatrix_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zexp  = Teuchos::rcp(new Epetra_Vector(dmatrix_->RowMap()));
    if (dmatrix_->RowMap().NumGlobalElements()) LINALG::Export(*z_,*zexp);
    dmatrix_->Multiply(true,*zexp,*tempvecd);
    LINALG::Export(*tempvecd,*fsexp);
    feff->Update(-(1.0-alphaf_),*fsexp,1.0);

    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecm  = Teuchos::rcp(new Epetra_Vector(mmatrix_->DomainMap()));
    mmatrix_->Multiply(true,*zexp,*tempvecm);
    LINALG::Export(*tempvecm,*fmexp);
    feff->Update(1.0-alphaf_,*fmexp,1.0);

    // add old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecdold  = Teuchos::rcp(new Epetra_Vector(dold_->DomainMap()));
    Teuchos::RCP<Epetra_Vector> zoldexp  = Teuchos::rcp(new Epetra_Vector(dold_->RowMap()));
    if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*zoldexp);
    dold_->Multiply(true,*zoldexp,*tempvecdold);
    LINALG::Export(*tempvecdold,*fsoldexp);
    feff->Update(-alphaf_,*fsoldexp,1.0);

    Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    Teuchos::RCP<Epetra_Vector> tempvecmold  = Teuchos::rcp(new Epetra_Vector(mold_->DomainMap()));
    mold_->Multiply(true,*zoldexp,*tempvecmold);
    LINALG::Export(*tempvecmold,*fmoldexp);
    feff->Update(alphaf_,*fmoldexp,1.0);
  }
  // if there is no self contact everything is ok
  else
  {
    // add contact force terms
    Teuchos::RCP<Epetra_Vector> fs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dmatrix_->Multiply(true,*z_,*fs);
    Teuchos::RCP<Epetra_Vector> fsexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fs,*fsexp);
    feff->Update(-(1.0-alphaf_),*fsexp,1.0);

    Teuchos::RCP<Epetra_Vector> fm = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mmatrix_->Multiply(true,*z_,*fm);
    Teuchos::RCP<Epetra_Vector> fmexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fm,*fmexp);
    feff->Update(1.0-alphaf_,*fmexp,1.0);

    ///////////////////////////////////////////////////////////////////
    //// FOR STATIC PROBLEMS --> alphaf_=0 !!! --> this is not needed!
    // add old contact forces (t_n)
    Teuchos::RCP<Epetra_Vector> fsold = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    dold_->Multiply(true,*zold_,*fsold);
    Teuchos::RCP<Epetra_Vector> fsoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fsold,*fsoldexp);
    feff->Update(-alphaf_,*fsoldexp,1.0);

    Teuchos::RCP<Epetra_Vector> fmold = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    mold_->Multiply(true,*zold_,*fmold);
    Teuchos::RCP<Epetra_Vector> fmoldexp = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
    LINALG::Export(*fmold,*fmoldexp);
    feff->Update(alphaf_,*fmoldexp,1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Solve linear system of saddle point type                  farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::SaddlePointSolve(LINALG::Solver& solver,
                  LINALG::Solver& fallbacksolver,
                  Teuchos::RCP<LINALG::SparseOperator> kdd,  Teuchos::RCP<Epetra_Vector> fd,
                  Teuchos::RCP<Epetra_Vector>  sold, Teuchos::RCP<Epetra_Vector> dirichtoggle,
                  int numiter)
{
  // get system type
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

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

  // initialize merged system (matrix, rhs, sol);
  Teuchos::RCP<Epetra_Map>           mergedmap   = Teuchos::null;
  if (!weardiscr_)
    mergedmap   = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
  else if (weardiscr_ and !wearbothdiscr_)
  {
    Teuchos::RCP<Epetra_Map> map_dummy = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
    mergedmap   = LINALG::MergeMap(map_dummy,gwdofrowmap_,false);
  }
  else
  {
    Teuchos::RCP<Epetra_Map> map_dummy = LINALG::MergeMap(ProblemDofs(),glmdofrowmap_,false);
    mergedmap   = LINALG::MergeMap(map_dummy,gwdofrowmap_,false);
    mergedmap   = LINALG::MergeMap(mergedmap,gwmdofrowmap_,false); // slave + master wear
  }

  Teuchos::RCP<LINALG::SparseMatrix> mergedmt    = Teuchos::null;
  Teuchos::RCP<Epetra_Vector>        mergedrhs   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedsol   = LINALG::CreateVector(*mergedmap);
  Teuchos::RCP<Epetra_Vector>        mergedzeros = LINALG::CreateVector(*mergedmap);

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> constrrhs = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // initialize transformed constraint matrices
  Teuchos::RCP<LINALG::SparseMatrix> trkdz, trkzd, trkzz;

  // Wear stuff for own discretization:
  // initialize wear r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> wearrhs  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> wearrhsM = Teuchos::null; // for master

  if (weardiscr_) wearrhs = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
  if (weardiscr_ and wearbothdiscr_) wearrhsM = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
  // initialize transformed constraint matrices
  Teuchos::RCP<LINALG::SparseMatrix> trkwd, trkwz, trkww, trkzw, trkdw; // slave
  Teuchos::RCP<LINALG::SparseMatrix> trkwmd, trkwmz, trkwmwm, trkzwm, trkdwm; // master
  // currently slave + master coeff
  double wcoeff  = Params().get<double>("WEARCOEFF");
  double wcoeffM = Params().get<double>("WEARCOEFF_MASTER");

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 1: FRICTIONLESS CONTACT ************************************
  if (!friction_)
  {
    dserror("WearLagrangeStrategy::SaddlePointSolve: Wear called without friction!");
  }

  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 2: Wear as pp quantity *************************************
  else if (!weardiscr_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
    kdz->Add(*dmatrix_,true,1.0-alphaf_,1.0);
    kdz->Add(*mmatrix_,true,-(1.0-alphaf_),1.0);
    kdz->Complete(*gsdofrowmap_,*gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = MORTAR::MatrixColTransformGIDs(kdz,glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkdz = MORTAR::MatrixRowTransform(trkdz,ProblemDofs());

    // build constraint matrix kzd
    Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gactiven_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzd->Add(*linstickDIS_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzd->Add(*linslipDIS_,false,1.0,1.0);
    kzd->Complete(*gdisprowmap_,*gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = MORTAR::MatrixRowTransformGIDs(kzd,glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkzd = MORTAR::MatrixColTransform(trkzd,ProblemDofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzz->Add(*linstickLM_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzz->Add(*linslipLM_,false,1.0,1.0);

    if(wearimpl_)
    {
      // add C-fnc. derivatives w.r.t. lm-values to kzz
      if (gactiven_->NumGlobalElements()) kzz->Add(*wlinmatrix_,false,1.0,1.0);
      if (gslipt_->NumGlobalElements()) kzz->Add(*wlinmatrixsl_,false,1.0,1.0);

#ifdef CONSISTENTSTICK
    if (gstickt->NumGlobalElements()) kzz->Add(*wlinmatrixst_,false,1.0,1.0);
#endif
    }
    kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (MatrixRowColTransform)
    trkzz = MORTAR::MatrixRowColTransformGIDs(kzz,glmdofrowmap_,glmdofrowmap_);

    /****************************************************************************************
   ***                RIGHT-HAND SIDE                   ***
   ****************************************************************************************/


  // We solve for the incremental Langrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) * z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gactiven_->NumGlobalElements())
    {
      LINALG::Export(*g_,*gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*gact,*gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linstickRHS_,*stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linslipRHS_,*slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0,*gactexp,1.0);
    constrrhs->Update(1.0,*stickexp,1.0);
    constrrhs->Update(1.0,*slipexp,1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs; // set constraint rhs vector
  }
  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 3: Wear as discrete variable *******************************
  else if (weardiscr_ and !wearbothdiscr_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
    kdz->Add(*dmatrix_,true,1.0-alphaf_,1.0);
    kdz->Add(*mmatrix_,true,-(1.0-alphaf_),1.0);
    kdz->Complete(*gsdofrowmap_,*gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = MORTAR::MatrixColTransformGIDs(kdz,glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkdz = MORTAR::MatrixRowTransform(trkdz,ProblemDofs());

    // build constraint matrix kzd
    Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gactiven_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzd->Add(*linstickDIS_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzd->Add(*linslipDIS_,false,1.0,1.0);
    kzd->Complete(*gdisprowmap_,*gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = MORTAR::MatrixRowTransformGIDs(kzd,glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkzd = MORTAR::MatrixColTransform(trkzd,ProblemDofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzz->Add(*linstickLM_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzz->Add(*linslipLM_,false,1.0,1.0);

    kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (MatrixRowColTransform)
    trkzz = MORTAR::MatrixRowColTransformGIDs(kzz,glmdofrowmap_,glmdofrowmap_);


    // ***************************************************************************************************
    // additional wear
    // ***************************************************************************************************
    // build wear matrix kwd
    Teuchos::RCP<LINALG::SparseMatrix> kwd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if(sswear_)
    {
      if (gactiven_->NumGlobalElements())  kwd->Add(*lintdis_,false,-wcoeff,1.0);
      if (gactiven_->NumGlobalElements())  kwd->Add(*linedis_,false,1.0,1.0);
      if (gactiven_->NumGlobalElements())  kwd->Complete(*gdisprowmap_,*gsdofnrowmap_);
    }
    else
    {
      if (gslipn_->NumGlobalElements())  kwd->Add(*lintdis_,false,-wcoeff,1.0);
      if (gslipn_->NumGlobalElements())  kwd->Add(*linedis_,false,1.0,1.0);
      if (gslipn_->NumGlobalElements())  kwd->Complete(*gdisprowmap_,*gsdofnrowmap_);
    }

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwd = MORTAR::MatrixRowTransformGIDs(kwd,gwdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkwd = MORTAR::MatrixColTransform(trkwd,ProblemDofs());

    // *********************************
    // build wear matrix kwz
    Teuchos::RCP<LINALG::SparseMatrix> kwz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if(sswear_)
    {
      if (gactiven_->NumGlobalElements()) kwz->Add(*lintlm_,false,-wcoeff,1.0);
    }
    else
    {
      if (gslipn_->NumGlobalElements()) kwz->Add(*lintlm_,false,-wcoeff,1.0);
    }

    kwz->Complete(*gsdofrowmap_,*gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwz = MORTAR::MatrixRowColTransformGIDs(kwz,gwdofrowmap_,glmdofrowmap_);

    // *********************************
    // build wear matrix kww
    Teuchos::RCP<LINALG::SparseMatrix> kww = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if(sswear_)
    {
      if (gactiven_->NumGlobalElements()) kww->Add(*ematrix_,false,1.0,1.0);
    }
    else
    {
      if (gslipn_->NumGlobalElements()) kww->Add(*ematrix_,false,1.0,1.0);
    }

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> onesw = Teuchos::rcp(new Epetra_Vector(*gwinact_));
    onesw->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiagw = Teuchos::rcp(new LINALG::SparseMatrix(*onesw));
    onesdiagw->Complete();
    // build constraint matrix kzz
    if (gwinact_->NumGlobalElements()) kww->Add(*onesdiagw,false,1.0,1.0);

    kww->Complete(*gsdofnrowmap_,*gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkww = MORTAR::MatrixRowColTransformGIDs(kww,gwdofrowmap_,gwdofrowmap_);

    // *********************************
    // build wear matrix kzw
    Teuchos::RCP<LINALG::SparseMatrix>  kzw = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gactiven_->NumGlobalElements()) kzw->Add(*smatrixW_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements())   kzw->Add(*linslipW_,false,1.0,1.0);
    kzw->Complete(*gsdofnrowmap_,*gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzw = MORTAR::MatrixRowColTransformGIDs(kzw,glmdofrowmap_,gwdofrowmap_);

    /****************************************************************************************
   ***                RIGHT-HAND SIDE                   ***
   ****************************************************************************************/


  // We solve for the incremental Langrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) * z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gactiven_->NumGlobalElements())
    {
      LINALG::Export(*g_,*gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*gact,*gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linstickRHS_,*stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linslipRHS_,*slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0,*gactexp,1.0);
    constrrhs->Update(1.0,*stickexp,1.0);
    constrrhs->Update(1.0,*slipexp,1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs; // set constraint rhs vector

    // ***************************************************************************************************
    // additional wear-rhs
    // ***************************************************************************************************

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    LINALG::Export(*WearCondRhs_, *WearCondRhsexp);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    LINALG::Export(*inactiveWearRhs_, *inactiveWearRhsexp);

    wearrhs->Update(1.0,*WearCondRhsexp,1.0);
    wearrhs->Update(1.0,*inactiveWearRhsexp,1.0);
    wearrhs->ReplaceMap(*gwdofrowmap_);

    wearrhs_=wearrhs;
  }
  //**********************************************************************
  // build matrix and vector blocks
  //**********************************************************************
  // *** CASE 4: Wear as discrete variable on both sides******************
  else if (weardiscr_ and wearbothdiscr_)
  {
    // global stick dof map
    Teuchos::RCP<Epetra_Map> gstickt = LINALG::SplitMap(*gactivet_,*gslipt_);

    // build constraint matrix kdz
    Teuchos::RCP<LINALG::SparseMatrix> kdz = Teuchos::rcp(new LINALG::SparseMatrix(*gdisprowmap_,100,false,true));
    kdz->Add(*dmatrix_,true,1.0-alphaf_,1.0);
    kdz->Add(*mmatrix_,true,-(1.0-alphaf_),1.0);
    kdz->Complete(*gsdofrowmap_,*gdisprowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
    trkdz = MORTAR::MatrixColTransformGIDs(kdz,glmdofrowmap_);

    // transform parallel row distribution of constraint matrix kdz
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkdz = MORTAR::MatrixRowTransform(trkdz,ProblemDofs());

    // build constraint matrix kzd
    Teuchos::RCP<LINALG::SparseMatrix> kzd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gactiven_->NumGlobalElements()) kzd->Add(*smatrix_,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzd->Add(*linstickDIS_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzd->Add(*linslipDIS_,false,1.0,1.0);
    kzd->Complete(*gdisprowmap_,*gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzd = MORTAR::MatrixRowTransformGIDs(kzd,glmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkzd = MORTAR::MatrixColTransform(trkzd,ProblemDofs());

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Map> gidofs = LINALG::SplitMap(*gsdofrowmap_,*gactivedofs_);
    Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(*gidofs));
    ones->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiag = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
    onesdiag->Complete();

    // build constraint matrix kzz
    Teuchos::RCP<LINALG::SparseMatrix> kzz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gidofs->NumGlobalElements())    kzz->Add(*onesdiag,false,1.0,1.0);
    if (gstickt->NumGlobalElements()) kzz->Add(*linstickLM_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements()) kzz->Add(*linslipLM_,false,1.0,1.0);

    kzz->Complete(*gsdofrowmap_,*gsdofrowmap_);

    // transform constraint matrix kzz to lmdofmap (MatrixRowColTransform)
    trkzz = MORTAR::MatrixRowColTransformGIDs(kzz,glmdofrowmap_,glmdofrowmap_);


    // ***************************************************************************************************
    // additional wear SLAVE
    // ***************************************************************************************************
    // build wear matrix kwd
    Teuchos::RCP<LINALG::SparseMatrix> kwd = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if (gslipn_->NumGlobalElements())  kwd->Add(*lintdis_,false,-wcoeff,1.0);
    if (gslipn_->NumGlobalElements())  kwd->Add(*linedis_,false,1.0,1.0);
    if (gslipn_->NumGlobalElements())  kwd->Complete(*gdisprowmap_,*gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwd = MORTAR::MatrixRowTransformGIDs(kwd,gwdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkwd = MORTAR::MatrixColTransform(trkwd,ProblemDofs());

    // *********************************
    // build wear matrix kwz
    Teuchos::RCP<LINALG::SparseMatrix> kwz = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if (gslipn_->NumGlobalElements()) kwz->Add(*lintlm_,false,-wcoeff,1.0);
    kwz->Complete(*gsdofrowmap_,*gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwz = MORTAR::MatrixRowColTransformGIDs(kwz,gwdofrowmap_,glmdofrowmap_);

    // *********************************
    // build wear matrix kww
    Teuchos::RCP<LINALG::SparseMatrix> kww = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofnrowmap_,100,false,true));
    if (gslipn_->NumGlobalElements()) kww->Add(*ematrix_,false,1.0,1.0);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> onesw = Teuchos::rcp(new Epetra_Vector(*gwinact_));
    onesw->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiagw = Teuchos::rcp(new LINALG::SparseMatrix(*onesw));
    onesdiagw->Complete();
    // build constraint matrix kzz
    if (gwinact_->NumGlobalElements()) kww->Add(*onesdiagw,false,1.0,1.0);

    kww->Complete(*gsdofnrowmap_,*gsdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkww = MORTAR::MatrixRowColTransformGIDs(kww,gwdofrowmap_,gwdofrowmap_);

    // FOR SLAVE AND MASTER
    // ********************************* S+M
    // build wear matrix kzw
    Teuchos::RCP<LINALG::SparseMatrix>  kzw = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,false,true));
    if (gactiven_->NumGlobalElements()) kzw->Add(*smatrixW_,false,1.0,1.0);
    if (gslipt_->NumGlobalElements())   kzw->Add(*linslipW_,false,1.0,1.0);
    kzw->Complete(*galldofnrowmap_,*gsdofrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkzw = MORTAR::MatrixRowColTransformGIDs(kzw,glmdofrowmap_,gwalldofrowmap_);

    // ***************************************************************************************************
    // additional wear MASTER
    // ***************************************************************************************************
    // build wear matrix kwmd
    Teuchos::RCP<LINALG::SparseMatrix> kwmd = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofnrowmap_,100,false,true));
    if (gmslipn_->NumGlobalElements())  kwmd->Add(*lintdisM_,false,-wcoeffM,1.0);
    if (gmslipn_->NumGlobalElements())  kwmd->Add(*linedisM_,false,1.0,1.0);
    if (gmslipn_->NumGlobalElements())  kwmd->Complete(*gdisprowmap_,*gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmd = MORTAR::MatrixRowTransformGIDs(kwmd,gwmdofrowmap_);

    // transform parallel column distribution of constraint matrix kzd
    // (only necessary in the parallel redistribution case)
    if (ParRedist()) trkwmd = MORTAR::MatrixColTransform(trkwmd,ProblemDofs());

    // *********************************
    // build wear matrix kwmz
    Teuchos::RCP<LINALG::SparseMatrix> kwmz = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofnrowmap_,100,false,true));
    if (gmslipn_->NumGlobalElements()) kwmz->Add(*lintlmM_,false,-wcoeffM,1.0);
    kwmz->Complete(*gsdofrowmap_,*gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmz = MORTAR::MatrixRowColTransformGIDs(kwmz,gwmdofrowmap_,glmdofrowmap_);

    // *********************************
    // build wear matrix kwmwm
    Teuchos::RCP<LINALG::SparseMatrix> kwmwm = Teuchos::rcp(new LINALG::SparseMatrix(*gmdofnrowmap_,100,false,true));
    if (gmslipn_->NumGlobalElements()) kwmwm->Add(*ematrixM_,false,1.0,1.0);

    // build unity matrix for inactive dofs
    Teuchos::RCP<Epetra_Vector> oneswm = Teuchos::rcp(new Epetra_Vector(*gwminact_));
    oneswm->PutScalar(1.0);
    Teuchos::RCP<LINALG::SparseMatrix> onesdiagwm = Teuchos::rcp(new LINALG::SparseMatrix(*oneswm));
    onesdiagwm->Complete();
    // build constraint matrix kzz
    if (gwminact_->NumGlobalElements()) kwmwm->Add(*onesdiagwm,false,1.0,1.0);

    kwmwm->Complete(*gmdofnrowmap_,*gmdofnrowmap_);

    // transform constraint matrix kzd to lmdofmap (MatrixRowTransform)
    trkwmwm = MORTAR::MatrixRowColTransformGIDs(kwmwm,gwmdofrowmap_,gwmdofrowmap_);

    /****************************************************************************************
    ***                                   RIGHT-HAND SIDE                                 ***
    ****************************************************************************************/


  // We solve for the incremental Langrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  //
  // r = r_effdyn,co = r_effdyn + a_f * B_co(d_(n)) * z_(n) + (1-a_f) * B_co(d^(i)_(n+1)) * z^(i)_(n+1)

    // export weighted gap vector
    Teuchos::RCP<Epetra_Vector> gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gactiven_->NumGlobalElements())
    {
      LINALG::Export(*g_,*gact);
      gact->ReplaceMap(*gactiven_);
    }
    Teuchos::RCP<Epetra_Vector> gactexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*gact,*gactexp);

    // export stick and slip r.h.s.
    Teuchos::RCP<Epetra_Vector> stickexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linstickRHS_,*stickexp);
    Teuchos::RCP<Epetra_Vector> slipexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*linslipRHS_,*slipexp);

    // export inactive rhs
    Teuchos::RCP<Epetra_Vector> inactiverhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*inactiverhs_, *inactiverhsexp);

    // build constraint rhs (1)
    constrrhs->Update(1.0, *inactiverhsexp, 1.0);

    // build constraint rhs
    constrrhs->Update(-1.0,*gactexp,1.0);
    constrrhs->Update(1.0,*stickexp,1.0);
    constrrhs->Update(1.0,*slipexp,1.0);
    constrrhs->ReplaceMap(*glmdofrowmap_);

    constrrhs_ = constrrhs; // set constraint rhs vector

    // ***************************************************************************************************
    // additional wear-rhs SLAVE
    // ***************************************************************************************************

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    LINALG::Export(*WearCondRhs_, *WearCondRhsexp);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    LINALG::Export(*inactiveWearRhs_, *inactiveWearRhsexp);

    wearrhs->Update(1.0,*WearCondRhsexp,1.0);
    wearrhs->Update(1.0,*inactiveWearRhsexp,1.0);
    wearrhs->ReplaceMap(*gwdofrowmap_);

    wearrhs_=wearrhs;
    // ***************************************************************************************************
    // additional wear-rhs Master
    // ***************************************************************************************************
    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> WearCondRhsexpM = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
    LINALG::Export(*WearCondRhsM_, *WearCondRhsexpM);

    // export inactive wear rhs
    Teuchos::RCP<Epetra_Vector> inactiveWearRhsexpM = Teuchos::rcp(new Epetra_Vector(*gmdofnrowmap_));
    LINALG::Export(*inactiveWearRhsM_, *inactiveWearRhsexpM);

    wearrhsM->Update(1.0,*WearCondRhsexpM,1.0);
    wearrhsM->Update(1.0,*inactiveWearRhsexpM,1.0);
    wearrhsM->ReplaceMap(*gwmdofrowmap_);

    wearmrhs_=wearrhsM;
  }
  else
    dserror("unknown wear algorithm!");

  //**********************************************************************
  // build and solve saddle point system
  //**********************************************************************
  if (systype==INPAR::CONTACT::system_saddlepoint)
  {
    // apply Dirichlet conditions to (0,0) and (0,1) blocks
    Teuchos::RCP<Epetra_Vector> zeros   = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(),true));
    Teuchos::RCP<Epetra_Vector> rhscopy = Teuchos::rcp(new Epetra_Vector(*fd));
    LINALG::ApplyDirichlettoSystem(stiffmt,sold,rhscopy,zeros,dirichtoggle);
    trkdz->ApplyDirichlet(dirichtoggle,false);

    // row map (equals domain map) extractor
    std::vector<Teuchos::RCP<const Epetra_Map> > mapvec;
    mapvec.push_back(ProblemDofs());
    mapvec.push_back(glmdofrowmap_);
    if(weardiscr_)
    {
      mapvec.push_back(gwdofrowmap_);
      if(wearbothdiscr_)
        mapvec.push_back(gwmdofrowmap_);
    }

    // set a helper flag for the CheapSIMPLE preconditioner (used to detect, if Teuchos::nullspace has to be set explicitely)
    // do we need this? if we set the Teuchos::nullspace when the solver is constructed?
    solver.Params().set<bool>("CONTACT",true); // for simpler precond

    // MODIFICATION OF SYSTEM:
    // =======================
    // build 2x2 soe by inserting wear blocks into lm blocks
    // this results to a well-suited block system for the iterative solvers
    Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > mat = Teuchos::null;

    if(weardiscr_)
    {
      Teuchos::RCP<LINALG::SparseMatrix> trkgd = Teuchos::null;
      Teuchos::RCP<LINALG::SparseMatrix> trkgg = Teuchos::null;

      if(!wearbothdiscr_)
      {
        // merged map ws + wm + z
        Teuchos::RCP<Epetra_Map>     gmap   = Teuchos::null;
        gmap   = LINALG::MergeMap(gwdofrowmap_,glmdofrowmap_,false);

        // row map (equals domain map) extractor
        LINALG::MapExtractor rowmapext(*mergedmap,gmap,ProblemDofs());
        LINALG::MapExtractor dommapext(*mergedmap,gmap,ProblemDofs());

        mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));

        trkgd = Teuchos::rcp(new LINALG::SparseMatrix(*gmap,100,false,true));
        trkgg = Teuchos::rcp(new LINALG::SparseMatrix(*gmap,100,false,true));

        trkgd->Add(*trkzd,false,1.0,1.0);
        trkgd->Add(*trkwd,false,1.0,1.0);

        trkgg->Add(*trkzz,false,1.0,1.0);
        trkgg->Add(*trkwz,false,1.0,1.0);
        trkgg->Add(*trkww,false,1.0,1.0);
        trkgg->Add(*trkzw,false,1.0,1.0);

        trkgd->Complete(*ProblemDofs(),*gmap);
        trkgg->Complete(*gmap,*gmap);

        trkdz->UnComplete();
        trkdz->Complete(*gmap,*ProblemDofs());

      }
      // BOTH_SIDED DISCRETE WEAR
      else
      {
        // merged map ws + wm + z
        Teuchos::RCP<Epetra_Map>     gmap   = Teuchos::null;
        Teuchos::RCP<Epetra_Map> map_dummyg = LINALG::MergeMap(gwdofrowmap_,gwmdofrowmap_,false);
        gmap   = LINALG::MergeMap(map_dummyg,glmdofrowmap_,false);

        // row map (equals domain map) extractor
        LINALG::MapExtractor rowmapext(*mergedmap,gmap,ProblemDofs());
        LINALG::MapExtractor dommapext(*mergedmap,gmap,ProblemDofs());

        mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));


        trkgd = Teuchos::rcp(new LINALG::SparseMatrix(*gmap,100,false,true));
        trkgg = Teuchos::rcp(new LINALG::SparseMatrix(*gmap,100,false,true));

        trkgd->Add(*trkzd,false,1.0,1.0);
        trkgd->Add(*trkwd,false,1.0,1.0);
        trkgd->Add(*trkwmd,false,1.0,1.0);

        trkgg->Add(*trkzz,false,1.0,1.0);
        trkgg->Add(*trkwz,false,1.0,1.0);
        trkgg->Add(*trkww,false,1.0,1.0);
        trkgg->Add(*trkzw,false,1.0,1.0);
        trkgg->Add(*trkwmz,false,1.0,1.0);
        trkgg->Add(*trkwmwm,false,1.0,1.0);

        trkgd->Complete(*ProblemDofs(),*gmap);
        trkgg->Complete(*gmap,*gmap);

        trkdz->UnComplete();
        trkdz->Complete(*gmap,*ProblemDofs());
      }


      mat->Assign(0,0,View,*stiffmt);
      mat->Assign(0,1,View,*trkdz);
      mat->Assign(1,0,View,*trkgd);
      mat->Assign(1,1,View,*trkgg);

    }
    // without wear unknowns...
    else
    {
      LINALG::MultiMapExtractor rowmapext(*mergedmap,mapvec);
      LINALG::MultiMapExtractor dommapext(*mergedmap,mapvec);

      // build block matrix for SIMPLER
      mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(dommapext,rowmapext,81,false,false));

      mat->Assign(0,0,View,*stiffmt);
      mat->Assign(0,1,View,*trkdz);
      mat->Assign(1,0,View,*trkzd);
      mat->Assign(1,1,View,*trkzz);
    }


    // matrix filling done
    mat->Complete();

    // we also need merged rhs here
    Teuchos::RCP<Epetra_Vector> fresmexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*fd,*fresmexp);
    mergedrhs->Update(1.0,*fresmexp,1.0);
    Teuchos::RCP<Epetra_Vector> constrexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*constrrhs,*constrexp);
    mergedrhs->Update(1.0,*constrexp,1.0);

    // add wear rhs
    if (weardiscr_)
    {
      Teuchos::RCP<Epetra_Vector> wearexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
      LINALG::Export(*wearrhs,*wearexp);
      mergedrhs->Update(1.0,*wearexp,1.0);

      if(wearbothdiscr_)
      {
        Teuchos::RCP<Epetra_Vector> wearexpM = Teuchos::rcp(new Epetra_Vector(*mergedmap));
        LINALG::Export(*wearrhsM,*wearexpM);
        mergedrhs->Update(1.0,*wearexpM,1.0);
      }
    }

    // apply Dirichlet B.C. to mergedrhs and mergedsol
    Teuchos::RCP<Epetra_Vector> dirichtoggleexp = Teuchos::rcp(new Epetra_Vector(*mergedmap));
    LINALG::Export(*dirichtoggle,*dirichtoggleexp);
    LINALG::ApplyDirichlettoSystem(mergedsol,mergedrhs,mergedzeros,dirichtoggleexp);

    // SIMPLER preconditioning solver call
    solver.Solve(mat->EpetraOperator(),mergedsol,mergedrhs,true,numiter==0);
  }

  //**********************************************************************
  // invalid system types
  //**********************************************************************
  else dserror("ERROR: Invalid system type in SaddlePointSolve");

  //**********************************************************************
  // extract results for displacement and LM increments
  //**********************************************************************
  Teuchos::RCP<Epetra_Vector> sollm = Teuchos::rcp(new Epetra_Vector(*glmdofrowmap_));
  Teuchos::RCP<Epetra_Vector> solw  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> solwm = Teuchos::null;

  if(weardiscr_) solw      = Teuchos::rcp(new Epetra_Vector(*gwdofrowmap_));
  if(wearbothdiscr_) solwm = Teuchos::rcp(new Epetra_Vector(*gwmdofrowmap_));

  LINALG::MapExtractor mapextd(*mergedmap,ProblemDofs(),glmdofrowmap_);
  LINALG::MapExtractor mapextlm(*mergedmap,glmdofrowmap_,glmdofrowmap_);
  mapextd.ExtractCondVector(mergedsol,sold);
  mapextlm.ExtractCondVector(mergedsol,sollm);
  sollm->ReplaceMap(*gsdofrowmap_);

  if(weardiscr_)
  {
    LINALG::MapExtractor mapextw(*mergedmap,gwdofrowmap_,glmdofrowmap_);
    mapextw.ExtractCondVector(mergedsol,solw);
    solw->ReplaceMap(*gsdofnrowmap_);
  }
  if(wearbothdiscr_)
  {
    LINALG::MapExtractor mapextwm(*mergedmap,gwmdofrowmap_,glmdofrowmap_);
    mapextwm.ExtractCondVector(mergedsol,solwm);
    solwm->ReplaceMap(*gmdofnrowmap_);
  }

  if (IsSelfContact())
  // for self contact, slave and master sets may have changed,
  // thus we have to reinitialize the LM vector map
  {
    zincr_ = Teuchos::rcp(new Epetra_Vector(*sollm));
    LINALG::Export(*z_, *zincr_);             // change the map of z_
    z_ = Teuchos::rcp(new Epetra_Vector(*zincr_));
    zincr_->Update(1.0, *sollm, 0.0);           // save sollm in zincr_
    z_->Update(1.0, *zincr_, 1.0);              // update z_
  }
  else
  {
    zincr_->Update(1.0, *sollm, 0.0);
    z_->Update(1.0, *zincr_, 1.0);

    if (weardiscr_)
    {
      wincr_->Update(1.0, *solw, 0.0);
      w_->Update(1.0, *wincr_, 1.0);

      if(wearbothdiscr_)
      {
        wmincr_->Update(1.0, *solwm, 0.0);
        wm_->Update(1.0, *wmincr_, 1.0);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------*
|  Output de-weighted wear vector                             farah 09/14|
*-----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::OutputWear()
{
  if(weardiscr_)
  {
    if(weartimescales_)
    {
      // multiply the wear with its normal direction and store in wear_vector
      // loop over all interfaces
      for (int i=0; i<(int)interface_.size(); ++i)
      {
        // FIRST: get the wear values and the normal directions for the interface
        // loop over all slave row nodes on the current interface
        for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->SlaveRowNodes()->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* frinode = dynamic_cast<FriNode*>(node);

          // be aware of problem dimension
          int dim = Dim();
          int numdof = frinode->NumDof();
          if (dim!=numdof) dserror("ERROR: Inconsistency Dim <-> NumDof");

          // nodal normal vector and wear
          double nn[3];
          double wear = 0.0;

          for (int j=0;j<3;++j)
            nn[j]=frinode->MoData().n()[j];

          if (abs(frinode->FriDataPlus().wcurr()[0] + frinode->FriDataPlus().waccu()[0])>1e-12)
            wear = frinode->FriDataPlus().wcurr()[0] + frinode->FriDataPlus().waccu()[0];
          else
            wear=0.0;

          // find indices for DOFs of current node in Epetra_Vector
          // and put node values (normal and tangential stress components) at these DOFs
          std::vector<int> locindex(dim);

          for (int dof=0;dof<dim;++dof)
          {
             locindex[dof] = (wearoutput_->Map()).LID(frinode->Dofs()[dof]);
            (*wearoutput_)[locindex[dof]] = wear * nn[dof];
          }
        }
      }

      // MASTER!!!
      // multiply the wear with its normal direction and store in wear_vector
      // loop over all interfaces
      if(wearbothdiscr_)
      {
        for (int i=0; i<(int)interface_.size(); ++i)
        {
          // FIRST: get the wear values and the normal directions for the interface
          // loop over all slave row nodes on the current interface
          for (int j=0; j<interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
          {
            int gid = interface_[i]->MasterRowNodes()->GID(j);
            DRT::Node* node = interface_[i]->Discret().gNode(gid);
            if (!node) dserror("ERROR: Cannot find node with gid %",gid);
            FriNode* frinode = dynamic_cast<FriNode*>(node);

            // be aware of problem dimension
            int dim = Dim();
            int numdof = frinode->NumDof();
            if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

            // nodal normal vector and wear
            double nn[3];
            double wear = 0.0;

            for (int j=0;j<3;++j)
              nn[j]=frinode->MoData().n()[j];

            if (abs(frinode->FriDataPlus().wcurr()[0] + frinode->FriDataPlus().waccu()[0])>1e-12)
              wear = frinode->FriDataPlus().wcurr()[0] + frinode->FriDataPlus().waccu()[0];
            else
              wear=0.0;

            // find indices for DOFs of current node in Epetra_Vector
            // and put node values (normal and tangential stress components) at these DOFs
            std::vector<int> locindex(dim);

            for (int dof=0;dof<dim;++dof)
            {
               locindex[dof] = (wearoutput2_->Map()).LID(frinode->Dofs()[dof]);
              (*wearoutput2_)[locindex[dof]] = wear * nn[dof];
            }
          }
        }
      }
    }
    else
    {
      // multiply the wear with its normal direction and store in wear_vector
      // loop over all interfaces
      for (int i=0; i<(int)interface_.size(); ++i)
      {
        // FIRST: get the wear values and the normal directions for the interface
        // loop over all slave row nodes on the current interface
        for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
        {
          int gid = interface_[i]->SlaveRowNodes()->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* frinode = dynamic_cast<FriNode*>(node);

          // be aware of problem dimension
          int dim = Dim();
          int numdof = frinode->NumDof();
          if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

          // nodal normal vector and wear
          double nn[3];
          double wear = 0.0;

          for (int j=0;j<3;++j)
            nn[j]=frinode->MoData().n()[j];

          if (abs(frinode->FriDataPlus().wcurr()[0])>1e-15)
            wear = frinode->FriDataPlus().wcurr()[0];
          else
            wear=0.0;

          // find indices for DOFs of current node in Epetra_Vector
          // and put node values (normal and tangential stress components) at these DOFs
          std::vector<int> locindex(dim);

          for (int dof=0;dof<dim;++dof)
          {
             locindex[dof] = (wearoutput_->Map()).LID(frinode->Dofs()[dof]);
            (*wearoutput_)[locindex[dof]] = wear * nn[dof];
          }
        }
      }

      // MASTER!!!
      // multiply the wear with its normal direction and store in wear_vector
      // loop over all interfaces
      if(wearbothdiscr_)
      {
        for (int i=0; i<(int)interface_.size(); ++i)
        {
          // FIRST: get the wear values and the normal directions for the interface
          // loop over all slave row nodes on the current interface
          for (int j=0; j<interface_[i]->MasterRowNodes()->NumMyElements(); ++j)
          {
            int gid = interface_[i]->MasterRowNodes()->GID(j);
            DRT::Node* node = interface_[i]->Discret().gNode(gid);
            if (!node) dserror("ERROR: Cannot find node with gid %",gid);
            FriNode* frinode = dynamic_cast<FriNode*>(node);

            // be aware of problem dimension
            int dim = Dim();
            int numdof = frinode->NumDof();
            if (dim!=numdof) dserror("ERROR: Inconsistency Dim <-> NumDof");

            // nodal normal vector and wear
            double nn[3];
            double wear = 0.0;

            for (int j=0;j<3;++j)
              nn[j]=frinode->MoData().n()[j];

            if (abs(frinode->FriDataPlus().wcurr()[0])>1e-12)
              wear = frinode->FriDataPlus().wcurr()[0];
            else
              wear=0.0;

            // find indices for DOFs of current node in Epetra_Vector
            // and put node values (normal and tangential stress components) at these DOFs
            std::vector<int> locindex(dim);

            for (int dof=0;dof<dim;++dof)
            {
               locindex[dof] = (wearoutput2_->Map()).LID(frinode->Dofs()[dof]);
              (*wearoutput2_)[locindex[dof]] = wear * nn[dof];
            }
          }
        }
      }
    }
  }
  //**********************************************************************************************
  //**********************************************************************************************
  //                                        weighted wear
  //**********************************************************************************************
  //**********************************************************************************************
  else
  {
    // only for dual Lagrange multiplier so far
    // diagonality of mortar matrix D is assumed
    INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
    if (shapefcn == INPAR::MORTAR::shape_standard)
      dserror("ERROR: Evaluation of wear only for dual shape functions so far.");

    // vectors
    Teuchos::RCP<Epetra_Vector> wear_vector = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    Teuchos::RCP<Epetra_Vector> real_weara = Teuchos::rcp(new Epetra_Vector(*gactivedofs_,true));
    Teuchos::RCP<Epetra_Vector> real_wear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_,true));

    // solver
    Teuchos::RCP<Teuchos::ParameterList> solver_params = Teuchos::rcp(new Teuchos::ParameterList());
    solver_params->set("solver","umfpack");
    solver_params->set("symmetric",false);
    LINALG::Solver solver(solver_params, Comm(), NULL);

    // multiply the wear with its normal direction and store in wear_vector
    // loop over all interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      // FIRST: get the wear values and the normal directions for the interface
      // loop over all slave row nodes on the current interface
      for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* frinode = dynamic_cast<FriNode*>(node);

        // be aware of problem dimension
        int dim = Dim();
        int numdof = frinode->NumDof();
        if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

        // nodal normal vector and wear
        double nn[3];
        double wear = 0.0;

        for (int j=0;j<3;++j)
          nn[j]=frinode->MoData().n()[j];
        wear = frinode->FriDataPlus().Wear();

        // find indices for DOFs of current node in Epetra_Vector
        // and put node values (normal and tangential stress components) at these DOFs
        std::vector<int> locindex(dim);

        for (int dof=0;dof<dim;++dof)
        {
           locindex[dof] = (wear_vector->Map()).LID(frinode->Dofs()[dof]);
          (*wear_vector)[locindex[dof]] = wear * nn[dof];
        }
      }
    }

    // extract active parts of D matrix
    // matrices, maps
    Teuchos::RCP<LINALG::SparseMatrix> daa,dai,dia,dii;
    Teuchos::RCP<Epetra_Map> gidofs;

    // ****************************************************************
    // split the matrix
    // why is here an empty gidofs map instead of a full map used???
    // ****************************************************************
    LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,gidofs,gactivedofs_,gidofs,daa,dai,dia,dii);

    // extract active parts of wear vector
    Teuchos::RCP<Epetra_Vector> wear_vectora = Teuchos::rcp(new Epetra_Vector(*gactivedofs_,true));
    Teuchos::RCP<Epetra_Vector> wear_vectori = Teuchos::rcp(new Epetra_Vector(*gidofs));

    // split the vector
    LINALG::SplitVector(*gsdofrowmap_,*wear_vector,gactivedofs_,wear_vectora,gidofs,wear_vectori);

    // approx. undo the weighting of the wear by solving
    // D * w = w~
    // dmatrix_ * real_wear = wear_
    if(gactivedofs_->NumGlobalElements())
      solver.Solve(daa->EpetraMatrix(), real_weara, wear_vectora, true);

    Teuchos::RCP<Epetra_Vector> real_wearexp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(*real_weara,*real_wearexp);
    real_wear->Update(1.0,*real_wearexp,0.0);

    // different wear coefficients on both sides...
    double wearcoeff_s = Params().get<double>("WEARCOEFF", 0.0);
    double wearcoeff_m = Params().get<double>("WEARCOEFF_MASTER", 0.0);
    if (wearcoeff_s<1e-12)
      dserror("wcoeff negative!!!");

    double fac         = wearcoeff_s / (wearcoeff_s+wearcoeff_m);

    // copy the local part of real_wear into wearoutput_
    for (int i=0; i<(int)gsdofrowmap_->NumMyElements(); ++i)
    {
     int gid = gsdofrowmap_->MyGlobalElements()[i];
     double tmp = (*real_wear)[real_wear->Map().LID(gid)];
     (*wearoutput_)[wearoutput_->Map().LID(gid)] = tmp*fac;
    }

    /**********************************************************************
     * Here the wearoutput_ - vector is the unweighted ("real") wearvector.
     * To calculate the wearvector for the master surface we transform
     * the slavewear vector via w_2~ = M^T * D^-1 * w~. In addition, we
     * unweight the resulting vector by D_2^-1*w_2~ and get the final
     * unweighted wear vector.
     **********************************************************************/
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(Params(),"BOTH_SIDED_WEAR")
        == INPAR::CONTACT::wear_both_map)
    {
      // different wear coefficients on both sides...
      double wearcoeff_s = Params().get<double>("WEARCOEFF", 0.0);
      double wearcoeff_m = Params().get<double>("WEARCOEFF_MASTER", 0.0);
      if (wearcoeff_s<1e-12)
        dserror("wcoeff negative!!!");

      double fac         = wearcoeff_m / (wearcoeff_s+wearcoeff_m);

      // extract involved parts of d2 matrix
      // matrices, maps - i: involved ; n: non-involved
      Teuchos::RCP<LINALG::SparseMatrix> d2ii,d2in,d2ni,d2nn;
      Teuchos::RCP<Epetra_Map> gndofs; // non-involved dofs

      LINALG::SplitMatrix2x2(d2matrix_,gminvolveddofs_,gndofs,gminvolveddofs_,gndofs,d2ii,d2in,d2ni,d2nn);

      Teuchos::RCP<Epetra_Vector> wear_master = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      Teuchos::RCP<Epetra_Vector> real_wear2 = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));

      // now we calc the weighted wear on the master surface
      mmatrix_->Multiply(true,*real_wear,*wear_master);

      // extract involved parts of wear vector
      Teuchos::RCP<Epetra_Vector> wear2_real    = Teuchos::rcp(new Epetra_Vector(*gminvolveddofs_,true));
      Teuchos::RCP<Epetra_Vector> wear2_vectori = Teuchos::rcp(new Epetra_Vector(*gminvolveddofs_,true));
      Teuchos::RCP<Epetra_Vector> wear2_vectorn = Teuchos::rcp(new Epetra_Vector(*gndofs));

      // split the vector
      LINALG::SplitVector(*gmdofrowmap_,*wear_master,gminvolveddofs_,wear2_vectori,gndofs,wear2_vectorn);

      if(gminvolveddofs_->NumGlobalElements())
      {
        solver.Solve(d2ii->EpetraMatrix(), wear2_real, wear2_vectori, true);
      }

      Teuchos::RCP<Epetra_Vector> real_wear2exp = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
      LINALG::Export(*wear2_real,*real_wear2exp);
      real_wear2->Update(1.0,*real_wear2exp,0.0);

      // copy the local part of real_wear into wearoutput_
      for (int i=0; i<(int)gmdofrowmap_->NumMyElements(); ++i)
      {
       int gid = gmdofrowmap_->MyGlobalElements()[i];
       double tmp = (*real_wear2)[real_wear2->Map().LID(gid)];
       (*wearoutput2_)[wearoutput2_->Map().LID(gid)] = -(tmp*fac); // negative sign because on other interface side
       //--> this Wear-vector (defined on master side) is along slave-side normal field!
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::DoWriteRestart(Teuchos::RCP<Epetra_Vector>& activetoggle,
                                                 Teuchos::RCP<Epetra_Vector>& sliptoggle,
                                                 Teuchos::RCP<Epetra_Vector>& weightedwear,
                                                 Teuchos::RCP<Epetra_Vector>& realwear,
                                                 bool forcedrestart)
{
  //TODO: extend this function to forcedrestart -- write output for
  // last converged wear... see contact_lagrange_strategy.cpp!

  // initalize
  activetoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  if (friction_)
    sliptoggle = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  if (wear_)
  {
    weightedwear = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));   // weighted
    realwear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));       // unweighted
  }

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoWriteRestart: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0; j<interface_[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      // set value active / inactive in toggle vector
      if (cnode->Active()) (*activetoggle)[dof]=1;

      // set value slip / stick in the toggle vector
      if (friction_)
      {
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (frinode->FriData().Slip()) (*sliptoggle)[dof]=1;
        if (wear_)
        {
          (*weightedwear)[dof] = frinode->FriDataPlus().Wear();
        }
      }
    }
  }

  if (friction_ and wear_)
    realwear=  wearoutput_;

  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::Recover(Teuchos::RCP<Epetra_Vector> disi)
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!IsInContact() && !WasInContact() && !WasInContactLastTimeStep()) return;

  // shape function and system types
  INPAR::MORTAR::ShapeFcn shapefcn = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(Params(),"SHAPEFCN");
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM");

  //**********************************************************************
  //**********************************************************************
  // CASE A: CONDENSED SYSTEM (DUAL) + WEAR DISCR (DUAL)
  //**********************************************************************
  //**********************************************************************
  if ((systype == INPAR::CONTACT::system_condensed) && weardiscr_)
  {
    // double-check if this is a dual LM system
    if (shapefcn != INPAR::MORTAR::shape_dual && shapefcn != INPAR::MORTAR::shape_petrovgalerkin)
      dserror("Condensation only for dual LM");

    // **********************************************
    // LAGR MULT RECOVERING
    // extract slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disis = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disis);

    // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> disim = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
    if (gmdofrowmap_->NumGlobalElements()) LINALG::Export(*disi, *disim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> disin = Teuchos::rcp(new Epetra_Vector(*gndofrowmap_));
    if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(invd_,gactivedofs_,tempmap,gactivedofs_,tempmap,invda,tempmtx1,tempmtx2,tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> invdmod = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    invdmod->Add(*invda,false,1.0,1.0);
    invdmod->Complete();

    // approximate update
    //invdmod->Multiply(false,*fs_,*z_);

    // full update
    //zincr_->Update(1.0,*z_,0.0); // z_i

    z_->Update(1.0,*fs_,0.0);
    Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    kss_->Multiply(false,*disis,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksm_->Multiply(false,*disim,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksn_->Multiply(false,*disin,*mod);
    z_->Update(-1.0,*mod,1.0);
    dold_->Multiply(true,*zold_,*mod);
    z_->Update(-alphaf_,*mod,1.0);
    Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
    invdmod->Multiply(true,*zcopy,*z_);
    z_->Scale(1/(1-alphaf_));

    //zincr_->Update(-1.0,*z_,1.0); // zi-zi+1

    // **********************************************
    // WEAR RECOVERING
    //wincr_ up to w_
    // extract active slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disia = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
    if (gactivedofs_->NumGlobalElements()) LINALG::Export(*disi, *disia);

    // extract inactive slave displacements from disi
    Teuchos::RCP<Epetra_Vector> disii = Teuchos::rcp(new Epetra_Vector(*gidofs_));
    if (gidofs_->NumGlobalElements()) LINALG::Export(*disi,*disii);

    // recovering just for slipnodes!
    if(gslipnodes_->NumGlobalElements()>0)
    {
      // rhs
      Teuchos::RCP<Epetra_Vector> fwexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
      LINALG::Export(*fw_,*fwexp);
      wincr_->Update(1.0, *fwexp, 0.0);

      Teuchos::RCP<Epetra_Vector> modw = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));

      // neutral part
      dnblock_ = MORTAR::MatrixRowTransformGIDs(dnblock_,gsdofnrowmap_);
      dnblock_->Multiply(false,*disin,*modw);
      wincr_->Update(1.0, *modw, 1.0);

      // master part
      dmblock_ = MORTAR::MatrixRowTransformGIDs(dmblock_,gsdofnrowmap_);
      dmblock_->Multiply(false,*disim,*modw);
      wincr_->Update(1.0, *modw, 1.0);

      //active part (stick and slip)
      dablock_ = MORTAR::MatrixRowTransformGIDs(dablock_,gsdofnrowmap_);
      dablock_->Multiply(false,*disia,*modw);
      wincr_->Update(1.0, *modw, 1.0);

      // inactive part
      if (gidofs_->NumGlobalElements()>0)
      {
        diblock_ = MORTAR::MatrixRowTransformGIDs(diblock_,gsdofnrowmap_);
        diblock_->Multiply(false,*disii,*modw);
        wincr_->Update(1.0, *modw, 1.0);
      }
    }
    else
    {
      wincr_->Scale(0.0);
    }
    //wear rhs  for inactive/stick nodes
    Teuchos::RCP<Epetra_Vector> wrhsexp = Teuchos::rcp(new Epetra_Vector(*gsdofnrowmap_));
    LINALG::Export(*inactiveWearRhs_,*wrhsexp);
    wincr_->Update(1.0, *wrhsexp, 1.0);
    w_->Update(1.0, *wincr_, 1.0);
  }
  //**********************************************************************
  //**********************************************************************
  // CASE B: CONDENSED SYSTEM (DUAL) + WEAR IMPLICIT/EXPLICIT
  //**********************************************************************
  //**********************************************************************
  else if ((systype == INPAR::CONTACT::system_condensed) && !weardiscr_)
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
    if (gndofrowmap_->NumGlobalElements()) LINALG::Export(*disi,*disin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(invd_,gactivedofs_,tempmap,gactivedofs_,tempmap,invda,tempmtx1,tempmtx2,tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> invdmod = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,10));
    invdmod->Add(*invda,false,1.0,1.0);
    invdmod->Complete();

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/

    // for self contact, slave and master sets may have changed,
    // thus we have to export the products Dold * zold and Mold^T * zold to fit
    if (IsSelfContact())
    {
      // approximate update
      //z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      //invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      z_->Update(1.0,*fs_,0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->Multiply(false,*disis,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksm_->Multiply(false,*disim,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksn_->Multiply(false,*disin,*mod);
      z_->Update(-1.0,*mod,1.0);
      Teuchos::RCP<Epetra_Vector> mod2 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
      if (dold_->RowMap().NumGlobalElements()) LINALG::Export(*zold_,*mod2);
      Teuchos::RCP<Epetra_Vector> mod3 = Teuchos::rcp(new Epetra_Vector((dold_->RowMap())));
      dold_->Multiply(true,*mod2,*mod3);
      Teuchos::RCP<Epetra_Vector> mod4 = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      if (gsdofrowmap_->NumGlobalElements()) LINALG::Export(*mod3,*mod4);
      z_->Update(-alphaf_,*mod4,1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->Multiply(true,*zcopy,*z_);
      z_->Scale(1/(1-alphaf_));
    }
    else
    {
      // approximate update
      //invdmod->Multiply(false,*fs_,*z_);

      // full update
      z_->Update(1.0,*fs_,0.0);
      Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
      kss_->Multiply(false,*disis,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksm_->Multiply(false,*disim,*mod);
      z_->Update(-1.0,*mod,1.0);
      ksn_->Multiply(false,*disin,*mod);
      z_->Update(-1.0,*mod,1.0);
      dold_->Multiply(true,*zold_,*mod);
      z_->Update(-alphaf_,*mod,1.0);
      Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
      invdmod->Multiply(true,*zcopy,*z_);
      z_->Scale(1/(1-alphaf_));
    }
  }
  //**********************************************************************
  //**********************************************************************
  // CASE C: SADDLE POINT SYSTEM
  //**********************************************************************
  //**********************************************************************
  else
  {
    // do nothing (z_ was part of solution already)
  }

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);
  if (weardiscr_)
  {
      StoreNodalQuantities(MORTAR::StrategyBase::wupdate);
  }
  if (wearbothdiscr_)
    StoreNodalQuantities(MORTAR::StrategyBase::wmupdate);

  return;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::RedistributeContact(Teuchos::RCP<Epetra_Vector> dis)
{
  // get out of here if parallel redistribution is switched off
  // or if this is a single processor (serial) job
  if (!ParRedist() || Comm().NumProc()==1) return;

  for (int i=0; i<(int)interface_.size();++i)
  {
    interface_[i]->IsRedistributed() = false;
  }

  // decide whether redistribution should be applied or not
  double taverage = 0.0;
  double eaverage = 0;
  bool doredist = false;
  double max_balance = Params().get<double>("MAX_BALANCE");

  //**********************************************************************
  // (1) static redistribution: ONLY at time t=0 or after restart
  // (both cases can be identified via empty unbalance vectors)
  //**********************************************************************
  if (WhichParRedist()==INPAR::MORTAR::parredist_static)
  {
    // this is the first time step (t=0) or restart
    if ((int)tunbalance_.size()==0 && (int)eunbalance_.size()==0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k=0;k<(int)tunbalance_.size();++k) taverage += tunbalance_[k];
      taverage/=(int)tunbalance_.size();
      for (int k=0;k<(int)eunbalance_.size();++k) eaverage += eunbalance_[k];
      eaverage/=(int)eunbalance_.size();

      // delete balance factors of last time step
      tunbalance_.resize(0);
      eunbalance_.resize(0);

      // no redistribution
      doredist = false;
    }
  }

  //**********************************************************************
  // (2) dynamic redistribution: whenever system is out of balance
  //**********************************************************************
  else if (WhichParRedist()==INPAR::MORTAR::parredist_dynamic)
  {
    // this is the first time step (t=0) or restart
    if ((int)tunbalance_.size()==0 && (int)eunbalance_.size()==0)
    {
      // do redistribution
      doredist = true;
    }

    // this is a regular time step (neither t=0 nor restart)
    else
    {
      // compute average balance factors of last time step
      for (int k=0;k<(int)tunbalance_.size();++k) taverage += tunbalance_[k];
      taverage/=(int)tunbalance_.size();
      for (int k=0;k<(int)eunbalance_.size();++k) eaverage += eunbalance_[k];
      eaverage/=(int)eunbalance_.size();

      // delete balance factors of last time step
      tunbalance_.resize(0);
      eunbalance_.resize(0);

      // decide on redistribution
      // -> (we allow a maximum value of the balance measure in the
      // system as defined in the input parameter MAX_BALANCE, i.e.
      // the maximum local processor workload and the minimum local
      // processor workload for mortar evaluation of all interfaces
      // may not differ by more than (MAX_BALANCE - 1.0)*100%)
      // -> (moreover, we redistribute if in the majority of iteration
      // steps of the last time step there has been an unbalance in
      // element distribution, i.e. if eaverage >= 0.5)
      if (taverage >= max_balance || eaverage >= 0.5)
        doredist = true;
    }
  }

  // print balance information to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "**********************************************************" << std::endl;
    if (taverage>0)
    {
      printf("Parallel balance (time): %e (limit %e) \n",taverage,max_balance);
      printf("Parallel balance (eles): %e (limit %e) \n",eaverage,0.5);
    }
    else
      printf("Parallel balance: t=0/restart \n");
    std::cout << "**********************************************************" << std::endl;
  }

  // get out of here if simulation is still in balance
  if (!doredist) return;

  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // set old and current displacement state
  // (needed for search within redistribution)
  SetState("displacement",dis);
  SetState("olddisplacement",dis);

  // global flag for redistribution
  bool anyinterfacedone = false;

  // parallel redistribution of all interfaces
  for (int i=0; i<(int)interface_.size();++i)
  {
    // redistribute optimally among procs
    bool done = interface_[i]->Redistribute(i+1);

    // if redistribution has really been performed
    // (the previous method might have found that there
    // are no "close" slave elements and thud redistribution
    // might not be necessary ->indicated by boolean)
    if (done)
    {
      // call fill complete again
      interface_[i]->FillComplete(maxdof_);

      // print new parallel distribution
      interface_[i]->PrintParallelDistribution(i+1);

      // re-create binary search tree
      interface_[i]->CreateSearchTree();

      // set bool for redistribution
      interface_[i]->IsRedistributed() = true;

      // set global flag to TRUE
      anyinterfacedone = true;
    }
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  if (anyinterfacedone)
  {
    Setup(true,false);
    SetupWear(true,false);
  }

  // time measurement
  Comm().Barrier();
  double t_end = Teuchos::Time::wallTime()-t_start;
  if (Comm().MyPID()==0)
    std::cout << "\nTime for parallel redistribution.........." << t_end << " secs\n" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::DoReadRestart(IO::DiscretizationReader& reader,
                                                Teuchos::RCP<Epetra_Vector> dis)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = DRT::INPUT::IntegralValue<int>(Params(),"RESTART_WITH_CONTACT");

  // set restart displacement state
  SetState("displacement", dis);
  SetState("olddisplacement", dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  InitMortar();
  InitEvalInterface();
  AssembleMortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (Dualquadslave3d())
  {
    // modify dmatrix_
    Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*dmatrix_,false,*invtrafo_,false,false,false,true);
    dmatrix_ = temp;
  }

  // read restart information on actice set and slip set (leave sets empty
  // if this is a restart with contact of a non-contact simulation run)
  Teuchos::RCP<Epetra_Vector> activetoggle =Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
  if (!restartwithcontact) reader.ReadVector(activetoggle,"activetoggle");

  // friction
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;
  Teuchos::RCP<Epetra_Vector> realwear;

  if (friction_)
  {
    sliptoggle =Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    if (!restartwithcontact) reader.ReadVector(sliptoggle,"sliptoggle");
  }

  // wear
  if (wear_)
  {
    weightedwear = Teuchos::rcp(new Epetra_Vector(*gsnoderowmap_));
    reader.ReadVector(weightedwear, "weightedwear");

    realwear = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    reader.ReadVector(realwear, "realwear");
  }

  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    // currently this only works safely for 1 interface
    //if (i>0) dserror("ERROR: DoReadRestart: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0; j<(interface_[i]->SlaveRowNodes())->NumMyElements(); ++j)
    {
      int gid = (interface_[i]->SlaveRowNodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof]==1)
      {
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // set value active / inactive in cnode
        cnode->Active()=true;

        if (friction_)
        {
          // set value stick / slip in cnode
          // set wear value
          if ((*sliptoggle)[dof]==1) dynamic_cast<CONTACT::FriNode*>(cnode)->FriData().Slip()=true;
          if (wear_)
            dynamic_cast<CONTACT::FriNode*>(cnode)->FriDataPlus().Wear() = (*weightedwear)[dof];
        }
      }
    }
  }

  if (friction_ and wear_)
    wearoutput_=realwear;


  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  zold_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  if (!restartwithcontact) reader.ReadVector(LagrMult(),"lagrmultold");
  if (!restartwithcontact) reader.ReadVector(LagrMultOld(),"lagrmultold");

  // Lagrange multiplier increment is always zero (no restart value to be read)
  zincr_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));

  // store restart information on Lagrange multipliers into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmcurrent);
  StoreNodalQuantities(MORTAR::StrategyBase::lmold);

  // TODO: same procedure for discrete wear...

  // only for Uzawa Augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  INPAR::CONTACT::SolvingStrategy st = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
  if (st == INPAR::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    if (!restartwithcontact) reader.ReadVector(LagrMultUzawa(),"lagrmultold");
    StoreNodalQuantities(MORTAR::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  StoreDM("old");

  if (friction_)
  {
    StoreNodalQuantities(MORTAR::StrategyBase::activeold);
    StoreToOld(MORTAR::StrategyBase::dm);
  }

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
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_, interface_[i]->ActiveNodes(), false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_, interface_[i]->ActiveDofs(), false);
    gactiven_ = LINALG::MergeMap(gactiven_, interface_[i]->ActiveNDofs(), false);
    gactivet_ = LINALG::MergeMap(gactivet_, interface_[i]->ActiveTDofs(), false);
    if (friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_, interface_[i]->SlipNodes(), false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_, interface_[i]->SlipDofs(), false);
      gslipt_ = LINALG::MergeMap(gslipt_, interface_[i]->SlipTDofs(), false);
    }
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_=true;
    wasincontact_=true;
    wasincontactlts_=true;
  }

  // evaluate relative movement (jump)
  // needed because it is not called in the predictor of the
  // lagrange multiplier strategy
  EvaluateRelMov();

  // reset unbalance factors for redistribution
  // (during restart the interface has been evaluated once)
  tunbalance_.resize(0);
  eunbalance_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::UpdateActiveSetSemiSmooth()
{
  // FIXME: Here we do not consider zig-zagging yet!

  // get out gof here if not in the semi-smooth Newton case
  // (but before doing this, check if there are invalid active nodes)
  bool semismooth = DRT::INPUT::IntegralValue<int>(Params(),"SEMI_SMOOTH_NEWTON");
  if (!semismooth)
  {
    // loop over all interfaces
    for (int i=0; i<(int)interface_.size(); ++i)
    {
      // loop over all slave nodes on the current interface
      for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
      {
        int gid = interface_[i]->SlaveRowNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
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
          dserror("ERROR: Active node %i without any segment/cell attached",cnode->Id());
      }
    }
    return;
  }

  // get input parameter ftype
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(Params(),"FRICTION");

  // read weighting factor cn
  // (this is necessary in semi-smooth Newton case, as the search for the
  // active set is now part of the Newton iteration. Thus, we do not know
  // the active / inactive status in advance and we can have a state in
  // which both the condition znormal = 0 and wgap = 0 are violated. Here
  // we have to weigh the two violations via cn!
  double cn = Params().get<double>("SEMI_SMOOTH_CN");

  // assume that active set has converged and check for opposite
  activesetconv_=true;

  // loop over all interfaces
  for (int i=0; i<(int)interface_.size(); ++i)
  {
    //if (i>0) dserror("ERROR: UpdateActiveSet: Double active node check needed for n interfaces!");

    // loop over all slave nodes on the current interface
    for (int j=0;j<interface_[i]->SlaveRowNodes()->NumMyElements();++j)
    {
      int gid = interface_[i]->SlaveRowNodes()->GID(j);
      DRT::Node* node = interface_[i]->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);

      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // get scaling factor
      double scalefac=1.;
      if (DRT::INPUT::IntegralValue<int>(scontact_,"LM_NODAL_SCALE")==true &&
                cnode->MoData().GetScale() != 0.)
        scalefac = cnode->MoData().GetScale();

      // compute weighted gap
      double wgap = (*g_)[g_->Map().LID(gid)]/scalefac;

      // compute normal part of Lagrange multiplier
      double nz = 0.0;
      double nzold = 0.0;
      for (int k=0;k<3;++k)
      {
        nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];
        nzold += cnode->MoData().n()[k] * cnode->MoData().lmold()[k];
      }

      // friction
      double ct = Params().get<double>("SEMI_SMOOTH_CT");
      std::vector<double> tz (Dim()-1,0);
      std::vector<double> tjump (Dim()-1,0);
      double euclidean = 0.0;

      if (friction_)
      {
        // static cast
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute tangential parts and of Lagrange multiplier and incremental jumps
        for (int i=0;i<Dim();++i)
        {
          tz[0] += frinode->CoData().txi()[i]*frinode->MoData().lm()[i];
          if(Dim()==3) tz[1] += frinode->CoData().teta()[i]*frinode->MoData().lm()[i];

          if (DRT::INPUT::IntegralValue<int>(Params(),"GP_SLIP_INCR")==false)
          {
            tjump[0] += frinode->CoData().txi()[i]*frinode->FriData().jump()[i];
            if(Dim()==3) tjump[1] += frinode->CoData().teta()[i]*frinode->FriData().jump()[i];
          }
        }

        if (DRT::INPUT::IntegralValue<int>(Params(),"GP_SLIP_INCR")==true)
        {
          tjump[0] = frinode->FriData().jump_var()[0];
          if(Dim()==3) tjump[1] = frinode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm |tz+ct.tjump|
        std::vector<double> sum (Dim()-1,0);
        sum[0] = tz[0]+ct*tjump[0];
        if (Dim()==3) sum[1] = tz[1]+ct*tjump[1];
        if (Dim()==2) euclidean = abs(sum[0]);
        if (Dim()==3) euclidean = sqrt(sum[0]*sum[0]+sum[1]*sum[1]);
       }


      // check nodes of inactive set *************************************
      if (cnode->Active()==false)
      {
        // check for fulfilment of contact condition
        //if (abs(nz) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact inactive node condition violated "
        //       <<  "for node ID: " << cnode->Id() << std::endl;

        // check for penetration and/or tensile contact forces
        if (nz - cn*wgap > 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap > 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = true;
          activesetconv_ = false;

          // friction
          if (friction_)
          {
            // nodes coming into contact
            dynamic_cast<FriNode*>(cnode)->FriData().Slip() = true;
          }
        }
      }

      // check nodes of active set ***************************************
      else
      {
        // check for fulfilment of contact condition
        //if (abs(wgap) > 1e-8)
        //  std::cout << "ERROR: UpdateActiveSet: Exact active node condition violated "
        //       << "for node ID: " << cnode->Id() << std::endl;

        // check for tensile contact forces and/or penetration
        if (nz - cn*wgap <= 0) // no averaging of Lagrange multipliers
        //if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
        {
          cnode->Active() = false;
          activesetconv_ = false;

          // friction
          if (friction_) dynamic_cast<FriNode*>(cnode)->FriData().Slip() = false;
        }

        // only do something for friction
        else
        {
          // friction tresca
          if (ftype == INPAR::CONTACT::friction_tresca)
          {
            FriNode* frinode = dynamic_cast<FriNode*>(cnode);

            // CAREFUL: friction bound is now interface-local (popp 08/2012)
            double frbound = interface_[i]->IParams().get<double>("FRBOUND");

            if(frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if(euclidean-frbound <= 0) {}
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
              if(euclidean-frbound > 0) {}
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

            if(frinode->FriData().Slip() == false)
            {
              // check (euclidean)-frbound <= 0
              if(euclidean-frcoeff*(nz-cn*wgap) <= 1e-10) {}
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
              if(euclidean-frcoeff*(nz-cn*wgap) > -1e-10) {}
              // do nothing (slip was correct)
              else
              {
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
  Comm().SumAll(&localcheck,&convcheck,1);

  // active set is only converged, if converged on all procs
  // if not, increase no. of active set steps too
  if (convcheck!=Comm().NumProc())
  {
    activesetconv_ = false;
    activesetsteps_ += 1;
  }

  // also update special flag for semi-smooth Newton convergence
  activesetssconv_ = activesetconv_;

  // update zig-zagging history (shift by one)
  if (zigzagtwo_!=Teuchos::null) zigzagthree_  = Teuchos::rcp(new Epetra_Map(*zigzagtwo_));
  if (zigzagone_!=Teuchos::null) zigzagtwo_    = Teuchos::rcp(new Epetra_Map(*zigzagone_));
  if (gactivenodes_!=Teuchos::null) zigzagone_ = Teuchos::rcp(new Epetra_Map(*gactivenodes_));

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_  = Teuchos::null;
  gactiven_     = Teuchos::null;
  gactivet_     = Teuchos::null;
  gslipnodes_   = Teuchos::null;
  gslipdofs_    = Teuchos::null;
  gslipt_       = Teuchos::null;

  // for both-sided wear
  gminvolvednodes_  = Teuchos::null;
  gminvolveddofs_   = Teuchos::null;

  // for wear with own discretization
  gslipn_  = Teuchos::null;
  gwinact_ = Teuchos::null;
  gmslipn_  = Teuchos::null;
  gwminact_ = Teuchos::null;
  gmslipnodes_ = Teuchos::null;
  gmactivenodes_=Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i=0;i<(int)interface_.size();++i)
  {
    interface_[i]->BuildActiveSet();
    gactivenodes_ = LINALG::MergeMap(gactivenodes_,interface_[i]->ActiveNodes(),false);
    gactivedofs_ = LINALG::MergeMap(gactivedofs_,interface_[i]->ActiveDofs(),false);
    gactiven_ = LINALG::MergeMap(gactiven_,interface_[i]->ActiveNDofs(),false);
    gactivet_ = LINALG::MergeMap(gactivet_,interface_[i]->ActiveTDofs(),false);

    // for both-sided wear
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(scontact_,"BOTH_SIDED_WEAR")
        == INPAR::CONTACT::wear_both_map)
    {
      gminvolvednodes_ = LINALG::MergeMap(gminvolvednodes_, interface_[i]->InvolvedNodes(), false);
      gminvolveddofs_  = LINALG::MergeMap(gminvolveddofs_, interface_[i]->InvolvedDofs(), false);
    }

    if(friction_)
    {
      gslipnodes_ = LINALG::MergeMap(gslipnodes_,interface_[i]->SlipNodes(),false);
      gslipdofs_ = LINALG::MergeMap(gslipdofs_,interface_[i]->SlipDofs(),false);
      gslipt_ = LINALG::MergeMap(gslipt_,interface_[i]->SlipTDofs(),false);
    }

    if(weardiscr_ and wearbothdiscr_)
    {
      interface_[i]->BuildActiveSetMaster();
      gmslipn_= LINALG::MergeMap(gmslipn_,interface_[i]->SlipMasterNDofs(),false);
      gmslipnodes_=LINALG::MergeMap(gmslipnodes_,interface_[i]->SlipMasterNodes(),false);
      gmactivenodes_=LINALG::MergeMap(gmactivenodes_,interface_[i]->ActiveMasterNodes(),false);
    }
  }// end interface loop

  if(weardiscr_)
  {
    gslipn_  =  LINALG::SplitMap(*gslipdofs_,*gslipt_);
    if(sswear_)
      gwinact_ =  LINALG::SplitMap(*gsdofnrowmap_,*gactiven_);
    else
      gwinact_ =  LINALG::SplitMap(*gsdofnrowmap_,*gslipn_);

    if(wearbothdiscr_)
      gwminact_ = LINALG::SplitMap(*gmdofnrowmap_,*gmslipn_);

  }

  //****************************************

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
  if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)
  {
    // frictionless contact
    if (ActiveSetSteps()>2)
    {
      if (zigzagtwo_!=Teuchos::null)
      {
        if (zigzagtwo_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 1;
        }
      }

      if (zigzagthree_!=Teuchos::null)
      {
        if (zigzagthree_->SameAs(*gactivenodes_))
        {
          // detect zig-zagging
          zigzagging = 2;
        }
      }
    }
  } // if (ftype != INPAR::CONTACT::friction_tresca && ftype != INPAR::CONTACT::friction_coulomb)

  // output to screen
  if (Comm().MyPID()==0)
  {
    if (zigzagging==1)
    {
      std::cout << "DETECTED 1-2 ZIG-ZAGGING OF ACTIVE SET................." << std::endl;
    }
    else if (zigzagging==2)
    {
      std::cout << "DETECTED 1-2-3 ZIG-ZAGGING OF ACTIVE SET................" << std::endl;
    }
    else
    {
      // do nothing, no zig-zagging
    }
  }

  // reset zig-zagging history
  if (activesetconv_==true)
  {
    zigzagone_  = Teuchos::null;
    zigzagtwo_  = Teuchos::null;
    zigzagthree_= Teuchos::null;
  }

  // output of active set status to screen
  if (Comm().MyPID()==0 && activesetconv_==false)
    std::cout << "ACTIVE CONTACT SET HAS CHANGED... CHANGE No. " << ActiveSetSteps()-1 << std::endl;

  // update flag for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_=true;
    wasincontact_=true;
  }
  else
    isincontact_=false;

  return;
}

/*----------------------------------------------------------------------*
 |  Update Wear rhs for seq. staggered partitioned sol.      farah 11/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::UpdateWearDiscretIterate(bool store)
{
  if (store)
  {
    StoreNodalQuantities(MORTAR::StrategyBase::wold);
    if (wearbothdiscr_)
      StoreNodalQuantities(MORTAR::StrategyBase::wmold);
  }
  else
  {
    // loop over all interfaces
    for(int i=0;i<(int)interface_.size();++i)
    {
      for (int j=0; j<(int)interface_[i]->SlaveColNodes()->NumMyElements(); ++j)
      {
        int gid = interface_[i]->SlaveColNodes()->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cnode = dynamic_cast<FriNode*>(node);

        cnode->FriDataPlus().wcurr()[0] = 0.0;
        cnode->FriDataPlus().wold()[0]  = 0.0;
        cnode->FriDataPlus().waccu()[0] = 0.0;
      }
      if(wearbothdiscr_)
      {
        const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(interface_[i]->MasterRowNodes()));

        for (int j=0; j<(int)masternodes->NumMyElements(); ++j)
        {
          int gid = masternodes->GID(j);
          DRT::Node* node = interface_[i]->Discret().gNode(gid);
          if (!node) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cnode = dynamic_cast<FriNode*>(node);

          cnode->FriDataPlus().wcurr()[0] = 0.0;
          cnode->FriDataPlus().wold()[0]  = 0.0;
          cnode->FriDataPlus().waccu()[0] = 0.0;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Update Wear for different time scales                    farah 12/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearLagrangeStrategy::UpdateWearDiscretAccumulation(bool wearaccumulation)
{
  if(wearaccumulation)
  {
    if(weartimescales_)
      StoreNodalQuantities(MORTAR::StrategyBase::wupdateT);
  }

  return;
}
