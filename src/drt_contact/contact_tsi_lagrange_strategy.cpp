/*---------------------------------------------------------------------*/
/*!
\file contact_tsi_lagrange_strategy.cpp

\brief a derived strategy handling the Lagrange multiplier based TSI contact

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/

#include "Epetra_SerialComm.h"
#include "contact_tsi_lagrange_strategy.H"
#include "contact_interface.H"
#include "contact_tsi_interface.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_multiply.H"
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "contact_tsi_interface.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                             seitz 08/15|
 *----------------------------------------------------------------------*/
CONTACT::CoTSILagrangeStrategy::CoTSILagrangeStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap,
    const Epetra_Map* NodeRowMap,
    Teuchos::ParameterList params,
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
    int dim,
    Teuchos::RCP<const Epetra_Comm> comm,
    double alphaf,
    int maxdof)
    : CoLagrangeStrategy(data_ptr,DofRowMap,NodeRowMap,params,interface,dim,comm,alphaf,maxdof),
      tsi_alpha_(1.)
{

  return;
}

/*------------------------------------------------------------------------*
 | Assign general thermo contact state                         seitz 08/15|
 *------------------------------------------------------------------------*/
void CONTACT::CoTSILagrangeStrategy::SetState(
    const enum MORTAR::StateType& statetype,
    const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case MORTAR::state_temperature:
    {
      for (int j=0;j<(int)interface_.size(); ++j)
      {
        DRT::Discretization& idiscr = interface_[j]->Discret();
        Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*idiscr.DofColMap(),false));
        LINALG::Export(vec,*global);

        for (int i=0;i<idiscr.NumMyColNodes();++i)
        {
          CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscr.lColNode(i));
          if (node==NULL) dserror("cast failed");
          std::vector<double> mytemp(1);
          std::vector<int> lm(1,node->Dofs()[0]);

          DRT::UTILS::ExtractMyValues(*global,mytemp,lm);
          if (node->HasCoTSIData()) // in case the interface has not been initialized yet
            node->CoTSIData().Temp() = mytemp[0];
        }
      }
      break;
    }
    case MORTAR::state_thermo_lagrange_multiplier:
    {
      for (int j=0;j<(int)interface_.size(); ++j)
      {
        DRT::Discretization& idiscr = interface_[j]->Discret();

        Teuchos::RCP<Epetra_Vector> global = Teuchos::rcp(new Epetra_Vector(*idiscr.DofColMap(),false));
        LINALG::Export(vec,*global);

        for (int i=0;i<idiscr.NumMyColNodes();++i)
        {
          CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscr.lColNode(i));
          std::vector<int> lm(1,node->Dofs()[0]);
          std::vector<double> myThermoLM(1,0.);
          DRT::UTILS::ExtractMyValues(*global,myThermoLM,lm);
          node->CoTSIData().ThermoLM() = myThermoLM[0];
        }
      }
      break;
    }
    default:
    {
      CONTACT::CoAbstractStrategy::SetState(statetype,vec);
      break;
    }
  }

  return;
}


void CONTACT::CoTSILagrangeStrategy::Evaluate(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>& combined_RHS,
    Teuchos::RCP<ADAPTER::Coupling> coupST,
    Teuchos::RCP<const Epetra_Vector> dis,
    Teuchos::RCP<const Epetra_Vector> temp,
    Teuchos::RCP<const LINALG::MapExtractor> str_dbc,
    Teuchos::RCP<const LINALG::MapExtractor> thr_dbc
    )
{
  if (thr_s_dofs_==Teuchos::null)
    thr_s_dofs_  = coupST->MasterToSlaveMap(gsdofrowmap_);

  // set the new displacements
  SetState(MORTAR::state_new_displacement, *dis);

  for (unsigned i=0;i<interface_.size();++i)
    interface_[i]->Initialize();

  // set new temperatures
  Teuchos::RCP<Epetra_Vector> temp2 = coupST()->SlaveToMaster(temp);
  SetState(MORTAR::state_temperature,*temp2);

  // error checks
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(Params(),"SYSTEM")
      !=INPAR::CONTACT::system_condensed)
    dserror("only condensed system implemented");

  // First, we need to evaluate all the interfaces
  InitMortar(); // initialize mortar matrices and vectors
  InitEvalInterface(); // evaluate mortar terms (integrate...)
  AssembleMortar();

  // get the relative movement for frictional contact
  EvaluateRelMov();

  // update active set
  UpdateActiveSetSemiSmooth();

  // init lin-matrices
  Initialize();

  // get the necessary maps on the thermo dofs
  Teuchos::RCP<Epetra_Map> gactive_themo_dofs = coupST->MasterToSlaveMap(gactivedofs_);
  Teuchos::RCP<Epetra_Map> master_thermo_dofs = coupST->MasterToSlaveMap(gmdofrowmap_);
  Teuchos::RCP<Epetra_Map> thr_act_dofs = coupST->MasterToSlaveMap(gactivedofs_);
  Teuchos::RCP<Epetra_Map> thr_m_dofs   = coupST->MasterToSlaveMap(gmdofrowmap_);
  Teuchos::RCP<Epetra_Map> thr_sm_dofs  = coupST->MasterToSlaveMap(gsmdofrowmap_);
  Teuchos::RCP<Epetra_Map> thr_all_dofs = Teuchos::rcp(new Epetra_Map(*coupST->SlaveDofMap()));

  // assemble the constraint lines for the active contact nodes
  Teuchos::RCP<LINALG::SparseMatrix> dcsdd =
      Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  LINALG::SparseMatrix dcsdT  (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  Teuchos::RCP<LINALG::SparseMatrix> dcsdLMc =
      Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  Teuchos::RCP<Epetra_Vector> rcsa = LINALG::CreateVector(*gactivedofs_,true);
  Teuchos::RCP<Epetra_Vector> g_all;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
    g_all = LINALG::CreateVector(*gsdofrowmap_,true);
  else
    g_all = LINALG::CreateVector(*gsnoderowmap_,true);

  // assemble linearization of heat conduction (thermal contact)
  LINALG::SparseMatrix dcTdd  (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix dcTdT  (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix dcTdLMc (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix dcTdLMt (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  Teuchos::RCP<Epetra_Vector> rcTa = LINALG::CreateVector(*gactivedofs_,true);

  // D and M matrix for the active nodes
  Teuchos::RCP<LINALG::SparseMatrix> dInv = Teuchos::rcp(new LINALG::SparseMatrix(*gsdofrowmap_,100,true,false));
  LINALG::SparseMatrix s  (*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);

  LINALG::SparseMatrix m_LinDissDISP            (*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix m_LinDissContactLM       (*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);

  // setup some linearizations
  LINALG::SparseMatrix linDcontactLM            (*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix linMcontactLM            (*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix linMdiss                 (*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix linMThermoLM             (*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix linDThermoLM             (*gsdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix m_LinDissContactLM_thrRow(*thr_m_dofs  ,100,true,false,LINALG::SparseMatrix::FE_MATRIX);

  // stick / slip linearization
  for (unsigned i=0;i<interface_.size();++i)
  {
    CONTACT::CoTSIInterface* tsi_interface = dynamic_cast<CONTACT::CoTSIInterface*>(&(*interface_[i]));
    if (!tsi_interface)
      dserror("in TSI contact, this should be a CoTSIInterface!");

    // linearized normal contact
    interface_[i]->AssembleS(s);
    interface_[i]->AssembleG(*g_all);

    // linearized tangential contact (friction)
    if (friction_)
    {
      tsi_interface->AssembleLinSlip (*dcsdLMc,*dcsdd,dcsdT,*rcsa);
      tsi_interface->AssembleLinStick(*dcsdLMc,*dcsdd,dcsdT,*rcsa);
    }
    else
    {
      tsi_interface->AssembleTN(dcsdLMc,Teuchos::null);
      tsi_interface->AssembleTNderiv(dcsdd,Teuchos::null);
      tsi_interface->AssembleTangrhs(*rcsa);
    }

    // linearized thermal contact (heat conduction)
    tsi_interface->AssembleLinConduct(dcTdd,dcTdT,dcTdLMt,dcTdLMc);

    tsi_interface->AssembleDM_linDiss(NULL,&m_LinDissDISP,NULL,&m_LinDissContactLM,1.);

    tsi_interface->AssembleLinDM(linDcontactLM,linMcontactLM);
    tsi_interface->AssembleLinDM_X(NULL,&linMdiss,1.,CONTACT::CoTSIInterface::LinDM_Diss,gsnoderowmap_);
    tsi_interface->AssembleLinDM_X(&linDThermoLM,&linMThermoLM,1.,CONTACT::CoTSIInterface::LinDM_ThermoLM,gsnoderowmap_);

  }

  // complete all those linearizations
  //                             colmap        rowmap
  linDcontactLM     .Complete(*gsmdofrowmap_,*gsdofrowmap_);
  linMcontactLM     .Complete(*gsmdofrowmap_,*gmdofrowmap_);
  m_LinDissDISP     .Complete(*gsmdofrowmap_,*gmdofrowmap_);
  linMdiss          .Complete(*gsmdofrowmap_,*gmdofrowmap_);
  linMThermoLM      .Complete(*gsmdofrowmap_,*gmdofrowmap_);
  linDThermoLM      .Complete(*gsmdofrowmap_,*gsdofrowmap_);
  m_LinDissContactLM.Complete(*gactivedofs_,*gmdofrowmap_);
  s                 .Complete(*gsmdofrowmap_,*gactivedofs_);

  dcsdd->Add(s,false,1.,-1.);
  dcsdLMc->Scale(-1.);
  dcsdT.Scale(-1.);
  rcsa->Scale(1.);

  // normal contact
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*gactivedofs_,true);
    if (gact->GlobalLength())
      LINALG::Export(*g_all,*gact);
  }
  else
  {
    gact = LINALG::CreateVector(*gactivenodes_,true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_all,*gact);
      if(gact->ReplaceMap(*gactiven_)) dserror("replaceMap went wrong");
    }
  }
  AddVector(*gact,*rcsa);
  rcsa->Norm2(&mech_contact_res_);

  // complete all the new matrix blocks
  // Note: since the contact interace assemled them, they are all based
  //       on displacement row and col maps. Hence, some still need to be transformed
  dcsdd->Complete(*gsmdofrowmap_,*gactivedofs_);
  dcsdT.Complete(*gsmdofrowmap_,*gactivedofs_);
  dcsdLMc->Complete(*gactivedofs_,*gactivedofs_);
  dcTdd.Complete(*gsmdofrowmap_,*gactivedofs_);
  dcTdT.Complete(*gsmdofrowmap_,*gactivedofs_);
  dcTdLMc.Complete(*gactivedofs_,*gactivedofs_);

  // get the seperate blocks of the 2x2 TSI block system
  // View mode!!! Since we actually want to add things there
  Teuchos::RCP<LINALG::SparseMatrix> kss = Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(0,0),LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> kst = Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(0,1),LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> kts = Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(1,0),LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> ktt = Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(1,1),LINALG::Copy));
  kss->UnComplete();
  kts->UnComplete();

  // split rhs
  Teuchos::RCP<Epetra_Vector> rs = Teuchos::rcp(new Epetra_Vector(*gdisprowmap_,true));
  Teuchos::RCP<Epetra_Vector> rt = Teuchos::rcp(new Epetra_Vector(*coupST->SlaveDofMap(),true));
  LINALG::Export(*combined_RHS,*rs);
  LINALG::Export(*combined_RHS,*rt);

  // we don't want the rhs but the residual
  rs->Scale(-1.);
  rt->Scale(-1.);

  // add last time step contact forces to rhs
  if (fscn_!=Teuchos::null) // in the first time step, we don't have any history of the
                            // contact force, after that, fscn_ should be initialized propperly
  {
    Epetra_Vector tmp(*gdisprowmap_);
    LINALG::Export(*fscn_,tmp);
    if (rs->Update(alphaf_,tmp,1.)!=0) // fscn already scaled with alphaf_ in update
      dserror("update went wrong");
  }

  if (ftcn_!=Teuchos::null)
  {
    Epetra_Vector tmp(*coupST->SlaveDofMap());
    LINALG::Export(*ftcn_,tmp);
    if (rt->Update((1.-tsi_alpha_),tmp,1.)!=0)
      dserror("update went wrong");
  }

  // map containing the inactive and non-contact structural dofs
  Teuchos::RCP<Epetra_Map> str_gni_dofs = LINALG::SplitMap(
      *LINALG::SplitMap(*gdisprowmap_,*gmdofrowmap_),*gactivedofs_);
  // map containing the inactive and non-contact thermal dofs
  Teuchos::RCP<Epetra_Map> thr_gni_dofs = coupST->MasterToSlaveMap(str_gni_dofs);

  // add to kss
  kss->Add(linDcontactLM,false,1.-alphaf_,1.);
  kss->Add(linMcontactLM,false,1.-alphaf_,1.);

  // transform and add to kts
  FSI::UTILS::MatrixRowTransform()(m_LinDissDISP,+tsi_alpha_,
    ADAPTER::CouplingMasterConverter(*coupST),*kts,true);
  FSI::UTILS::MatrixRowTransform()(linMdiss,-tsi_alpha_,     // this minus sign is there, since assemble linM does not actually
    ADAPTER::CouplingMasterConverter(*coupST),*kts,true);    // assemble the linearization of M but the negative linearization of M
  FSI::UTILS::MatrixRowTransform()(linMThermoLM,tsi_alpha_,
    ADAPTER::CouplingMasterConverter(*coupST),*kts,true);
  FSI::UTILS::MatrixRowTransform()(linDThermoLM,tsi_alpha_,
    ADAPTER::CouplingMasterConverter(*coupST),*kts,true);

  FSI::UTILS::MatrixRowTransform().operator ()(m_LinDissContactLM,1.,
    ADAPTER::CouplingMasterConverter(*coupST),m_LinDissContactLM_thrRow,false);
  m_LinDissContactLM_thrRow.Complete(*gactivedofs_,*thr_m_dofs);

  // complete the matrix blocks again, now that we have added
  // the additional displacement linearizations
  kss->Complete();
  kts->Complete(*gdisprowmap_,*coupST->SlaveDofMap());

  // split matrix blocks in 3 rows: Active, Master and (Inactive+others)
  Teuchos::RCP<LINALG::SparseMatrix> kss_ni, kss_m, kss_a,
                                     kst_ni, kst_m, kst_a,
                                     kts_ni, kts_m, kts_a,
                                     ktt_ni, ktt_m, ktt_a,
                                     dummy1,dummy2,dummy3;

  // temporary matrix
  Teuchos::RCP<LINALG::SparseMatrix> tmp;
  Teuchos::RCP<Epetra_Vector> tmpv;

  // an empty dummy map
  Teuchos::RCP<Epetra_Map> dummy_map1,dummy_map2;

  // ****************************************************
  // split kss block*************************************
  // ****************************************************
  // split first row
  LINALG::SplitMatrix2x2(kss,str_gni_dofs,dummy_map1,gdisprowmap_,dummy_map2,
      kss_ni,dummy1,tmp,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(tmp,gmdofrowmap_,dummy_map1,gdisprowmap_,dummy_map2,
      kss_m,dummy1,kss_a,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;
  tmp=Teuchos::null;
  // ****************************************************
  // split kss block*************************************
  // ****************************************************

  // ****************************************************
  // split kst block*************************************
  // ****************************************************
  // split first row
  LINALG::SplitMatrix2x2(kst,str_gni_dofs,dummy_map1,thr_all_dofs,dummy_map2,
      kst_ni,dummy1,tmp,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(tmp,gmdofrowmap_,dummy_map1,thr_all_dofs,dummy_map2,
      kst_m,dummy1,kst_a,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;
  tmp=Teuchos::null;
  // ****************************************************
  // split kst block*************************************
  // ****************************************************

  // ****************************************************
  // split kts block*************************************
  // ****************************************************
  // split first row
  LINALG::SplitMatrix2x2(kts,thr_gni_dofs,dummy_map1,gdisprowmap_,dummy_map2,
      kts_ni,dummy1,tmp,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(tmp,thr_m_dofs,dummy_map1,gdisprowmap_,dummy_map2,
      kts_m,dummy1,kts_a,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;
  tmp=Teuchos::null;
  // ****************************************************
  // split kts block*************************************
  // ****************************************************

  // ****************************************************
  // split ktt block*************************************
  // ****************************************************
  // split first row
  LINALG::SplitMatrix2x2(ktt,thr_gni_dofs,dummy_map1,thr_all_dofs,dummy_map2,
      ktt_ni,dummy1,tmp,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(tmp,thr_m_dofs,dummy_map1,thr_all_dofs,dummy_map2,
      ktt_m,dummy1,ktt_a,dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements()!=0
      ||
      dummy2->DomainMap().NumGlobalElements()!=0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1=Teuchos::null; dummy2=Teuchos::null; dummy_map1=Teuchos::null; dummy_map2=Teuchos::null;
  tmp=Teuchos::null;
  // ****************************************************
  // split ktt block*************************************
  // ****************************************************

  // ****************************************************
  // split rhs vectors***********************************
  // ****************************************************
  // split structural rhs
  Epetra_Vector rsni(*str_gni_dofs); LINALG::Export(*rs,rsni);
  Epetra_Vector rsm (*gmdofrowmap_); LINALG::Export(*rs,rsm);
  Teuchos::RCP<Epetra_Vector> rsa = Teuchos::rcp(new Epetra_Vector(*gactivedofs_)); LINALG::Export(*rs,*rsa);

  // split thermal rhs
  Epetra_Vector rtni(*thr_gni_dofs); LINALG::Export(*rt,rtni);
  Epetra_Vector rtm (*thr_m_dofs  ); LINALG::Export(*rt,rtm);
  Teuchos::RCP<Epetra_Vector> rta =Teuchos::rcp(new Epetra_Vector(*thr_act_dofs)); LINALG::Export(*rt,*rta);
  // ****************************************************
  // split rhs vectors***********************************
  // ****************************************************

  // D and M matrix for the active nodes
  Teuchos::RCP<LINALG::SparseMatrix> dInvA = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,100,true,false));
  Teuchos::RCP<LINALG::SparseMatrix> mA    = Teuchos::rcp(new LINALG::SparseMatrix(*gactivedofs_,100,true,false));

  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  LINALG::SplitMatrix2x2(dmatrix_,gactivedofs_,dummy_map1,gactivedofs_,dummy_map2,dInvA,dummy1,dummy2,dummy3);
  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  LINALG::SplitMatrix2x2(mmatrix_,gactivedofs_,dummy_map1,gmdofrowmap_,dummy_map2,mA,dummy1,dummy2,dummy3);

  // now we have added the additional linearizations.
  // if there are no active nodes, we can leave now

  if (gactivenodes_->NumGlobalElements()==0)
  {
    sysmat->Reset();
    sysmat->Assign(0,0,LINALG::Copy,*kss);
    sysmat->Assign(0,1,LINALG::Copy,*kst);
    sysmat->Assign(1,0,LINALG::Copy,*kts);
    sysmat->Assign(1,1,LINALG::Copy,*ktt);
    return;
  }


  // we need to add another term, since AssembleLinStick/Slip assumes that we solve
  // for the Lagrange multiplier increments. However, we solve for the LM directly.
  // We can do that, since the system is linear in the LMs.
  tmpv=Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  Teuchos::RCP<Epetra_Vector> tmpv2 = Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  LINALG::Export(*z_,*tmpv2);
  dcsdLMc->Multiply(false,*tmpv2,*tmpv);
  tmpv->Scale(-1.);
  AddVector(*tmpv,*rcsa);
  tmpv = Teuchos::null;
  tmpv2= Teuchos::null;

  dcTdLMt.Complete(*gsdofrowmap_,*gactivedofs_);
  LINALG::SparseMatrix test(*gactivedofs_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  Teuchos::RCP<LINALG::SparseMatrix> a1(&dcTdLMt,false);
  Teuchos::RCP<LINALG::SparseMatrix> a2;
  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  LINALG::SplitMatrix2x2(a1,gactivedofs_,dummy_map1,gactivedofs_,dummy_map2,a2,dummy1,dummy2,dummy3);
  dcTdLMt=*a2;

  dcTdLMt.Complete(*gactivedofs_,*gactivedofs_);
  dInvA->Complete(*gactivedofs_,*gactivedofs_);
  mA->Complete(*gmdofrowmap_,*gactivedofs_);

  LINALG::SparseMatrix dcTdLMc_thr(*thr_act_dofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix dcTdLMt_thr(*thr_act_dofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  FSI::UTILS::MatrixRowTransform()(dcTdLMc,1.,ADAPTER::CouplingMasterConverter(*coupST),dcTdLMc_thr,true);
  FSI::UTILS::MatrixRowColTransform()(dcTdLMt,1.,ADAPTER::CouplingMasterConverter(*coupST),
      ADAPTER::CouplingMasterConverter(*coupST),dcTdLMt_thr,true,false);
  dcTdLMc_thr.Complete(*gactivedofs_,*thr_act_dofs);
  dcTdLMt_thr.Complete(*thr_act_dofs,*thr_act_dofs);

  // invert D-matrix
  Epetra_Vector dDiag(*gactivedofs_);
  dInvA->ExtractDiagonalCopy(dDiag);
  if (dDiag.Reciprocal(dDiag)) dserror("inversion of diagonal D matrix failed");
  dInvA->ReplaceDiagonalValues(dDiag);

  // get dinv on thermal dofs
  Teuchos::RCP<LINALG::SparseMatrix> dInvaThr = Teuchos::rcp(new LINALG::SparseMatrix(*thr_act_dofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  FSI::UTILS::MatrixRowColTransform()(*dInvA,1.,ADAPTER::CouplingMasterConverter(*coupST),
      ADAPTER::CouplingMasterConverter(*coupST),*dInvaThr,false,false);
  dInvaThr->Complete(*thr_act_dofs,*thr_act_dofs);

  // save some matrix blocks for recovery
  dinvA_=dInvA;
  dinvAthr_=dInvaThr;
  kss_a_=kss_a;
  kst_a_=kst_a;
  kts_a_=kts_a;
  ktt_a_=ktt_a;
  rs_a_=rsa;
  rt_a_=rta;
  thr_act_dofs_=thr_act_dofs;

  // get dinv * M
  Teuchos::RCP<LINALG::SparseMatrix> dInvMa = LINALG::MLMultiply(*dInvA,false,*mA,false,false,false,true);

  // get dinv * M on the thermal dofs
  LINALG::SparseMatrix dInvMaThr(*thr_act_dofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  FSI::UTILS::MatrixRowColTransform()(*dInvMa,1.,ADAPTER::CouplingMasterConverter(*coupST),
      ADAPTER::CouplingMasterConverter(*coupST),dInvMaThr,false,false);
  dInvMaThr.Complete(*thr_m_dofs,*thr_act_dofs);

  // apply contact symmetry conditions
  if (constr_direction_==INPAR::CONTACT::constr_xyz)
  {
    double haveDBC=0;
    pgsdirichtoggle_->Norm1(&haveDBC);
    if (haveDBC>0.)
    {
      Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*gactivedofs_,true);
      dInvA->ExtractDiagonalCopy(*diag);
      Teuchos::RCP<Epetra_Vector> lmDBC=LINALG::CreateVector(*gactivedofs_,true);
      LINALG::Export(*pgsdirichtoggle_,*lmDBC);
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*gactivedofs_,true);
      tmp->Multiply(1.,*diag,*lmDBC,0.);
      diag->Update(-1.,*tmp,1.);
      dInvA->ReplaceDiagonalValues(*diag);
      dInvMa = LINALG::MLMultiply(*dInvA,false,*mA,false,false,false,true);
    }
  }

  // reset the tangent stiffness
  // (for the condensation we have constructed copies above)
  sysmat->Reset();
  sysmat->UnComplete();

  // need diagonal block kss with explicitdirichtlet_=true
  // to be able to apply dirichlet values for contact symmetry condition
  LINALG::SparseMatrix tmpkss(*gdisprowmap_,100,false,false,LINALG::SparseMatrix::FE_MATRIX);
  sysmat->Assign(0,0,LINALG::Copy,tmpkss);

  // get references to the blocks (just for convenience)
  LINALG::SparseMatrix& kss_new = sysmat->Matrix(0,0);
  LINALG::SparseMatrix& kst_new = sysmat->Matrix(0,1);
  LINALG::SparseMatrix& kts_new = sysmat->Matrix(1,0);
  LINALG::SparseMatrix& ktt_new = sysmat->Matrix(1,1);

  // reset rhs
  combined_RHS->Scale(0.);

  // **********************************************************************
  // **********************************************************************
  // BUILD CONDENSED SYSTEM
  // **********************************************************************
  // **********************************************************************

  // (1) add the blocks, we do nothing with (i.e. (Inactive+others))
  kss_new.Add(*kss_ni,false,1.,1.);
  kst_new.Add(*kst_ni,false,1.,1.);
  kts_new.Add(*kts_ni,false,1.,1.);
  ktt_new.Add(*ktt_ni,false,1.,1.);
  AddVector(rsni,*combined_RHS);
  AddVector(rtni,*combined_RHS);

  // (2) add the 'uncondensed' blocks (i.e. everything w/o a D^-1
  // (2)a actual stiffness blocks of the master-rows
  kss_new.Add(*kss_m,false,1.,1.);
  kst_new.Add(*kst_m,false,1.,1.);
  kts_new.Add(*kts_m,false,1.,1.);
  ktt_new.Add(*ktt_m,false,1.,1.);
  AddVector(rsm,*combined_RHS);
  AddVector(rtm,*combined_RHS);

  // (2)b active constraints in the active slave rows
  kss_new.Add(*dcsdd,false,1.,1.);

  FSI::UTILS::MatrixColTransform()(*gactivedofs_,*gsmdofrowmap_,dcsdT,1.,
      ADAPTER::CouplingMasterConverter(*coupST),kst_new,false,true);
  FSI::UTILS::MatrixRowTransform()(dcTdd,1.,ADAPTER::CouplingMasterConverter(*coupST),kts_new,true);
  FSI::UTILS::MatrixRowColTransform()(dcTdT,1.,ADAPTER::CouplingMasterConverter(*coupST),
      ADAPTER::CouplingMasterConverter(*coupST),ktt_new,true,true);
  AddVector(*rcsa,*combined_RHS);

  // (3) condensed parts
  // second row
  kss_new.Add(*LINALG::MLMultiply(*dInvMa,true,*kss_a,false,false,false,true),false,1.,1.);
  kst_new.Add(*LINALG::MLMultiply(*dInvMa,true,*kst_a,false,false,false,true),false,1.,1.);
  tmpv = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  dInvMa->Multiply(true,*rsa,*tmpv);
  AddVector(*tmpv,*combined_RHS);
  tmpv = Teuchos::null;

  // third row
  Teuchos::RCP<LINALG::SparseMatrix> wDinv = LINALG::MLMultiply(*dcsdLMc,false,*dInvA,true,false,false,true);
  kss_new.Add(*LINALG::MLMultiply(*wDinv,false,*kss_a,false,false,false,true),false,-1./(1.-alphaf_),1.);
  kst_new.Add(*LINALG::MLMultiply(*wDinv,false,*kst_a,false,false,false,true),false,-1./(1.-alphaf_),1.);
  tmpv=Teuchos::rcp(new Epetra_Vector(*gactivedofs_));
  wDinv->Multiply(false,*rsa,*tmpv);
  tmpv->Scale(-1./(1.-alphaf_));
  AddVector(*tmpv,*combined_RHS);
  tmpv = Teuchos::null;
  wDinv = Teuchos::null;

  // fourth row: no condensation. Terms already added in (1)

  // fifth row
  tmp = Teuchos::null;
  tmp = LINALG::MLMultiply(m_LinDissContactLM_thrRow,false,*dInvA,false,false,false,true);
  kts_new.Add(*LINALG::MLMultiply(*tmp,false,*kss_a,false,false,false,true),false,-tsi_alpha_/(1.-alphaf_),1.);
  ktt_new.Add(*LINALG::MLMultiply(*tmp,false,*kst_a,false,false,false,true),false,-tsi_alpha_/(1.-alphaf_),1.);
  tmpv=Teuchos::rcp(new Epetra_Vector(*thr_m_dofs));
  tmp->Multiply(false,*rsa,*tmpv);
  tmpv->Scale(-tsi_alpha_/(1.-alphaf_));
  AddVector(*tmpv,*combined_RHS);
  tmpv=Teuchos::null;

  kts_new.Add(*LINALG::MLMultiply(dInvMaThr,true,*kts_a,false,false,false,true),false,1.,1.);
  ktt_new.Add(*LINALG::MLMultiply(dInvMaThr,true,*ktt_a,false,false,false,true),false,1.,1.);
  tmpv=Teuchos::rcp(new Epetra_Vector(*thr_m_dofs));
  dInvMaThr.Multiply(true,*rta,*tmpv);
  AddVector(*tmpv,*combined_RHS);
  tmp=Teuchos::null;

  // sixth row
  Teuchos::RCP<LINALG::SparseMatrix> yDinv = LINALG::MLMultiply(dcTdLMc_thr,false,*dInvA,false,false,false,true);
  kts_new.Add(*LINALG::MLMultiply(*yDinv,false,*kss_a,false,false,false,true),false,-1./(1.-alphaf_),1.);
  ktt_new.Add(*LINALG::MLMultiply(*yDinv,false,*kst_a,false,false,false,true),false,-1./(1.-alphaf_),1.);
  tmpv=Teuchos::rcp(new Epetra_Vector(*thr_act_dofs));
  yDinv->Multiply(false,*rsa,*tmpv);
  tmpv->Scale(-1./(1.-alphaf_));
  AddVector(*tmpv,*combined_RHS);
  tmpv=Teuchos::null;

  Teuchos::RCP<LINALG::SparseMatrix> gDinv = LINALG::MLMultiply(dcTdLMt_thr,false,*dInvaThr,false,false,false,true);
  kts_new.Add(*LINALG::MLMultiply(*gDinv,false,*kts_a,false,false,false,true),false,-1./(tsi_alpha_),1.);
  ktt_new.Add(*LINALG::MLMultiply(*gDinv,false,*ktt_a,false,false,false,true),false,-1./(tsi_alpha_),1.);
  tmpv=Teuchos::rcp(new Epetra_Vector(*thr_act_dofs));
  gDinv->Multiply(false,*rta,*tmpv);
  tmpv->Scale(-1./tsi_alpha_);
  AddVector(*tmpv,*combined_RHS);

  // and were done with the system matrix
  sysmat->Complete();

  // we need to return the rhs, not the residual
  combined_RHS->Scale(-1.);

  return;

}


void CONTACT::CoTSILagrangeStrategy::AddVector(Epetra_Vector& src,Epetra_Vector& dst)
{

  // return if src has no elements
  if (src.GlobalLength()==0)
    return;

#ifdef DEBUG
  for (int i=0;i<src.Map().NumMyElements();++i)
    if ((dst.Map().LID(src.Map().GID(i)))<0)
      dserror("src is not a vector on a sub-map of dst");
#endif

  Epetra_Vector tmp = Epetra_Vector(dst.Map(),true);
  LINALG::Export(src,tmp);
  if(dst.Update(1.,tmp,1.)) dserror("vector update went wrong");
  return;
}

void CONTACT::CoTSILagrangeStrategy::RecoverCoupled(
    Teuchos::RCP<Epetra_Vector> sinc,
    Teuchos::RCP<Epetra_Vector> tinc,
    Teuchos::RCP<ADAPTER::Coupling> coupST)
{
  Teuchos::RCP<Epetra_Vector> z_old=Teuchos::null;
  if (z_!=Teuchos::null)
    z_old = Teuchos::rcp(new Epetra_Vector(*z_));
  Teuchos::RCP<Epetra_Vector> z_thr_old=Teuchos::null;
  if (z_thr_!=Teuchos::null)
    z_thr_old = Teuchos::rcp(new Epetra_Vector(*z_thr_));

  // recover contact LM
  if (gactivedofs_->NumGlobalElements()>0)
  {
    // do we have everything we need?
    if (   rs_a_ ==Teuchos::null
        || kss_a_==Teuchos::null
        || kst_a_==Teuchos::null
        || dinvA_==Teuchos::null
    )
      dserror("some data for LM recovery is missing");

    Epetra_Vector lmc_a_new(*gactivedofs_,false);
    Epetra_Vector tmp(*gactivedofs_,false);
    lmc_a_new.Update(1.,*rs_a_,0.);
    kss_a_->Multiply(false,*sinc,tmp);
    lmc_a_new.Update(1.,tmp,1.);
    kst_a_->Multiply(false,*tinc,tmp);
    lmc_a_new.Update(1.,tmp,1.);
    dinvA_->Multiply(false,lmc_a_new,tmp);
    tmp.Scale(-1./(1.-alphaf_));
    z_=Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    LINALG::Export(tmp,*z_);

    // recover thermo LM
    // do we have everything we need?
    if (   rt_a_    ==Teuchos::null
        || kts_a_   ==Teuchos::null
        || ktt_a_   ==Teuchos::null
        || dinvAthr_==Teuchos::null
    )
      dserror("some data for LM recovery is missing");

    Epetra_Vector lmt_a_new(*thr_act_dofs_,false);
    Epetra_Vector tmp2(*thr_act_dofs_,false);
    lmt_a_new.Update(1.,*rt_a_,0.);
    kts_a_->Multiply(false,*sinc,tmp2);
    lmt_a_new.Update(1.,tmp2,1.);
    ktt_a_->Multiply(false,*tinc,tmp2);
    lmt_a_new.Update(1.,tmp2,1.);
    dinvAthr_->Multiply(false,lmt_a_new,tmp2);
    tmp2.Scale(-1./(tsi_alpha_));
    z_thr_=Teuchos::rcp(new Epetra_Vector(*thr_s_dofs_));
    LINALG::Export(tmp2,*z_thr_);
  }

  else
  {
    z_=Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
    z_thr_=Teuchos::rcp(new Epetra_Vector(*thr_s_dofs_));
  }

  if (z_old!=Teuchos::null)
    {
      z_old->Update(-1.,*z_,1.);
      z_old->Norm2(&mech_contact_incr_);
    }
  if (z_thr_old!=Teuchos::null)
    {
      z_thr_old->Update(-1.,*z_thr_,1.);
      z_thr_old->Norm2(&thr_contact_incr_);
    }

  // store updated LM into nodes
  StoreNodalQuantities(MORTAR::StrategyBase::lmupdate,Teuchos::null);
  StoreNodalQuantities(MORTAR::StrategyBase::lmThermo,coupST);

  return;
};

void CONTACT::CoTSILagrangeStrategy::StoreNodalQuantities(
    MORTAR::StrategyBase::QuantityType type,Teuchos::RCP<ADAPTER::Coupling> coupST)
{
  Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;
  // start type switch
  switch (type)
  {
  case MORTAR::StrategyBase::lmThermo:
  {
    Teuchos::RCP<Epetra_Vector> tmp=Teuchos::rcp(new Epetra_Vector(*coupST->SlaveDofMap()));

    LINALG::Export(*z_thr_,*tmp);
    vectorglobal = z_thr_;
    vectorglobal = coupST->SlaveToMaster(tmp);
    Teuchos::RCP<Epetra_Map> sdofmap, snodemap;
    // loop over all interfaces
    for (int i = 0; i < (int) interface_.size(); ++i)
    {
      sdofmap  = interface_[i]->SlaveColDofs();
      snodemap = interface_[i]->SlaveColNodes();
      Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::null;
      vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofmap));
      if (vectorglobal != Teuchos::null)
        LINALG::Export(*vectorglobal, *vectorinterface);

      // loop over all slave nodes (column or row) on the current interface
      for (int j = 0; j < snodemap->NumMyElements(); ++j)
      {
        int gid = snodemap->GID(j);
        DRT::Node* node = interface_[i]->Discret().gNode(gid);
        if (!node)
          dserror("ERROR: Cannot find node with gid %", gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        cnode->CoTSIData().ThermoLM() = (*vectorinterface)[(vectorinterface->Map()).LID(cnode->Dofs()[0])];
      }
    }
    break;
  }
  default:
    CONTACT::CoAbstractStrategy::StoreNodalQuantities(type);
    break;
  }
}

void CONTACT::CoTSILagrangeStrategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  if (fscn_==Teuchos::null)
    fscn_ = Teuchos::rcp (new Epetra_Vector(*gsmdofrowmap_));
  fscn_->Scale(0.);

  if (ftcnp_==Teuchos::null)
    ftcnp_ = Teuchos::rcp (new Epetra_Vector(*coupST_->MasterToSlaveMap(gsmdofrowmap_)));
  ftcnp_->Scale(0.);

  Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_));
  dmatrix_->Multiply(false,*z_,*tmp);
  AddVector(*tmp,*fscn_);

  tmp = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  mmatrix_->Multiply(true,*z_,*tmp);
  tmp->Scale(-1.);
  AddVector(*tmp,*fscn_);

  CONTACT::CoAbstractStrategy::Update(dis);

  LINALG::SparseMatrix dThr(*coupST_->MasterToSlaveMap(gsdofrowmap_),100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  FSI::UTILS::MatrixRowColTransform()(*dmatrix_,1.,ADAPTER::CouplingMasterConverter(*coupST_),
      ADAPTER::CouplingMasterConverter(*coupST_),dThr,false,false);
  dThr.Complete();
  tmp = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterToSlaveMap(gsdofrowmap_)));
  if(dThr.Apply(*z_thr_,*tmp)!=0) dserror("apply went wrong");
  AddVector(*tmp,*ftcnp_);

  LINALG::SparseMatrix mThr(*coupST_->MasterToSlaveMap(gsdofrowmap_),100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  FSI::UTILS::MatrixRowColTransform()(*mmatrix_,1.,ADAPTER::CouplingMasterConverter(*coupST_),
      ADAPTER::CouplingMasterConverter(*coupST_),mThr,false,false);
  mThr.Complete(*coupST_->MasterToSlaveMap(gmdofrowmap_),*coupST_->MasterToSlaveMap(gsdofrowmap_));
  mThr.UseTranspose();
  tmp = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterToSlaveMap(gmdofrowmap_)));
  if (mThr.Multiply(true,*z_thr_,*tmp)!=0) dserror("multiply went wrong");
  tmp->Scale(-1.);
  AddVector(*tmp,*ftcnp_);

  LINALG::SparseMatrix m_LinDissContactLM(*gmdofrowmap_,100,true,false,LINALG::SparseMatrix::FE_MATRIX);
  for (unsigned i=0;i<interface_.size();++i)
    dynamic_cast<CONTACT::CoTSIInterface*>(&(*interface_[i]))->AssembleDM_linDiss(
        NULL,NULL,NULL,&m_LinDissContactLM,1.);
  m_LinDissContactLM.Complete(*gactivedofs_,*gmdofrowmap_);
  Teuchos::RCP<Epetra_Vector> z_act = Teuchos::rcp (new Epetra_Vector(*gactivedofs_));;
  LINALG::Export(*z_,*z_act);
  tmp = Teuchos::rcp(new Epetra_Vector(*gmdofrowmap_));
  if (m_LinDissContactLM.Multiply(false,*z_act,*tmp)!=0)
    dserror("multiply went wrong");
  Teuchos::RCP<Epetra_Vector>tmp2 = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterDofMap()));
  LINALG::Export(*tmp,*tmp2);
  Teuchos::RCP<Epetra_Vector>tmp3 = coupST_->MasterToSlave(tmp2);
  Teuchos::RCP<Epetra_Vector>tmp4 = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterToSlaveMap(gmdofrowmap_)));
  LINALG::Export(*tmp3,*tmp4);
  AddVector(*tmp4,*ftcnp_);

  ftcn_=ftcnp_;
}

void CONTACT::CoTSILagrangeStrategy::SetAlphafThermo(const Teuchos::ParameterList& tdyn)
{
  INPAR::THR::DynamicType dyn_type = DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn,"DYNAMICTYP");
  switch (dyn_type)
  {
  case INPAR::THR::dyna_genalpha:
    tsi_alpha_ = tdyn.sublist("GENALPHA").get<double>("ALPHA_F");
    break;
  case INPAR::THR::dyna_onesteptheta:
    tsi_alpha_ = tdyn.sublist("ONESTEPTHETA").get<double>("THETA");
    break;
  case INPAR::THR::dyna_statics:
    tsi_alpha_ = 1.;
    break;
  default:
    dserror("unknown thermal time integration type");
  }
  return;
}


/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoTSILagrangeStrategy::DoWriteRestart(
    std::map<std::string,Teuchos::RCP<Epetra_Vector> >& restart_vectors,
    bool forcedrestart) const
{
  CONTACT::CoAbstractStrategy::DoWriteRestart(restart_vectors,forcedrestart);

  if (fscn_!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
    LINALG::Export(*fscn_,*tmp);
    restart_vectors["last_contact_force"]=tmp;
  }
  if (ftcn_!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*coupST_->SlaveDofMap()));
    LINALG::Export(*ftcn_,*tmp);
  restart_vectors["last_thermo_force"]=coupST_->SlaveToMaster(tmp);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoTSILagrangeStrategy::DoReadRestart(
    IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis,
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  bool restartwithcontact = DRT::INPUT::IntegralValue<int>(Params(),
    "RESTART_WITH_CONTACT");

  CONTACT::CoAbstractStrategy::DoReadRestart(reader,dis);
  fscn_ = Teuchos::rcp(new Epetra_Vector(*gsmdofrowmap_));
  if (!restartwithcontact)
    reader.ReadVector(fscn_, "last_contact_force");

  Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterDofMap()));
  if (!restartwithcontact)
    reader.ReadVector(tmp,"last_thermo_force");
  ftcn_=coupST_->MasterToSlave(tmp);
  tmp = Teuchos::rcp(new Epetra_Vector(*coupST_->MasterToSlaveMap(gsmdofrowmap_)));
  LINALG::Export(*ftcn_,*tmp);
  ftcn_=tmp;
}
