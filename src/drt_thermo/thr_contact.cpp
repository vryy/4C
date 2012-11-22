/*----------------------------------------------------------------------*/
/*!
\file thr_contact.cpp
\brief Thermal contact routines for

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                                     06/11 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                    mgit 10/10 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_impl.H"
#include "thr_aux.H"
#include "thr_contact.H"

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"


/*----------------------------------------------------------------------*
 | constructor                                               mgit 06/11 |
 *----------------------------------------------------------------------*/
THR::ThermoContactMan::ThermoContactMan(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                                        Teuchos::RCP<DRT::Discretization> discretstruct,
                                        Teuchos::RCP<DRT::Discretization> discretthermo):
cmtman_(cmtman),
discretstruct_(discretstruct),
discretthermo_(discretthermo)
{
  // communicator
  comm_ = Teuchos::rcp(discretthermo_->Comm().Clone());

  // done so far
  return;
}

/*----------------------------------------------------------------------*
 | modify thermal system of equation towards thermal contact mgit 09/10 |
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::ApplyThermoContact(Teuchos::RCP<LINALG::SparseMatrix>& tang,
                                               Teuchos::RCP<Epetra_Vector>& feff,
                                               Teuchos::RCP<Epetra_Vector>& temp,
                                               double dt)
{
  
  //**********************************************************************
  // prepare / convert some maps 
  //**********************************************************************
 
  // complete stiffness matrix
  // (this is a prerequisite for the Split2x2 methods to be called later)
  tang->Complete();

  // convert maps (from structure discretization to thermo discretization)
  // slave-, active-, inactive-, master-, activemaster-, n- smdofs
  Teuchos::RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  ConvertMaps (sdofs,adofs,mdofs);

  // map of active and master dofs
  amdofs = LINALG::MergeMap(adofs,mdofs,false);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);

  // row map of thermal problem
  Teuchos::RCP<Epetra_Map> problemrowmap = Teuchos::rcp(new Epetra_Map(*(discretthermo_->DofRowMap(0))));
  
  // split problemrowmap in n+am
  ndofs = LINALG::SplitMap(*problemrowmap,*smdofs);

  //**********************************************************************
  // Modification of thermal system of equations towards contact
  // 1. Without Lagrange Multipliers (Laursen, Wriggers)
  // 2. With Lagrange Multipliers (Hüeber)
  //**********************************************************************

  // modifications only if there are active nodes
  if (adofs->NumGlobalElements()==0)
    return;

  // new matrices to build up 
  Teuchos::RCP<LINALG::SparseMatrix> tangmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(tang);
  Teuchos::RCP<LINALG::SparseMatrix> tangnew = Teuchos::rcp(new LINALG::SparseMatrix(*problemrowmap,81,true,false,tangmatrix->GetMatrixtype()));
  Teuchos::RCP<Epetra_Vector> feffnew = LINALG::CreateVector(*problemrowmap);

  // assemble Mortar Matrices D and M in thermo dofs for active nodes
  Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*sdofs,10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*sdofs,100));

  TransformDM(*dmatrix,*mmatrix,sdofs,mdofs);
  
  // FillComplete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*mdofs,*sdofs);

  // now the choice, with or without lagrange multipliers
  bool thermolagmult = DRT::INPUT::IntegralValue<int>(cmtman_->GetStrategy().Params(),"THERMOLAGMULT");

  if(thermolagmult == false)
  {
    //********************************************************************
    // 1. Modification without Lagrange Multipliers
    //********************************************************************

    // static cast of mortar strategy to contact strategy
    MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
    CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

    // get vector of contact interfaces
    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

    // this currently works only for one interface yet and for one heat
    // transfer coefficient
    // FIXGIT: The heat transfer coefficient should be a condition on
    // the single interfaces!!
    if (interface.size()>1)
      dserror("Error in TSI::Algorithm::AssembleThermContCondition: Only for one interface yet.");
   
    // heat transfer coefficient for slave and master surface
    double heattranss = interface[0]->IParams().get<double>("HEATTRANSSLAVE");
    double heattransm = interface[0]->IParams().get<double>("HEATTRANSMASTER");

    if (heattranss <= 0 or heattransm <= 0)
     dserror("Error: Choose realistic heat transfer parameter");

    double beta = heattranss*heattransm/(heattranss+heattransm);
    double delta = heattranss/(heattranss+heattransm);

    // assemble Matrix B in addition to D and M
    Teuchos::RCP<LINALG::SparseMatrix> bmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*mdofs,10));
    AssembleB(*bmatrix);
  
    // Fill Complete
    bmatrix->Complete();
    
    // assemble mechanical dissipation for master side
    Teuchos::RCP<Epetra_Vector> mechdissratemaster = LINALG::CreateVector(*mdofs,true);
    AssembleMechDissMaster(*mechdissratemaster,dt);
    
    // assemble mechanical dissipation for slave side
    Teuchos::RCP<Epetra_Vector> mechdissrateslave = LINALG::CreateVector(*sdofs,true);
    AssembleMechDissSlave(*mechdissrateslave,dt);
   
    // split matrices
    Teuchos::RCP<Epetra_Map> tmp;
    Teuchos::RCP<LINALG::SparseMatrix> dmatrixa,mmatrixa,mmatrixat,bmatrixa,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    LINALG::SplitMatrix2x2(dmatrix,adofs,idofs,sdofs,tmp,dmatrixa,tmp1,tmp2,tmp3);
    LINALG::SplitMatrix2x2(mmatrix,adofs,idofs,mdofs,tmp,mmatrixa,tmp4,tmp5,tmp6);

    // master node participating the active interface
    Teuchos::RCP<Epetra_Map> masterpar = Teuchos::rcp(new Epetra_Map(mmatrixa->ColMap()));
    Teuchos::RCP<Epetra_Map> masternotpar =  LINALG::SplitMap(*mdofs,*masterpar);
    
    // split matrices
    LINALG::SplitMatrix2x2(bmatrix,masterpar,masternotpar,mdofs,tmp,bmatrixa,tmp4,tmp5,tmp6);
    LINALG::SplitMatrix2x2(mmatrix,sdofs,tmp,masterpar,masternotpar,mmatrixat,tmp4,tmp5,tmp6);
    
    // FIXGIT, dserror
    // at the moment, it is not possible to extract from the matrices B 
    // and M the "active master nodes" part. Only segments with adjacent 
    // active slave nodes should contribute to B and M. This is not the case 
    // at the moment.
    if(masternotpar->NumGlobalElements()!= 0.0)
      dserror ("ERROR: Thermal contact without LM only possible for "
                "all masternodes belonging to the contact interface.");
    
    // modify left hand side --------------------------------------------- 
    tangnew->UnComplete();
    
    // add existing tang matrix
    tangnew->Add(*tang,false,1.0,1.0);
    
    // add matrices to master row nodes
    tangnew->Add(*bmatrixa,false,+beta,1.0);
    tangnew->Add(*mmatrixat,true,-beta,1.0);
 
    // add matrices to slave row nodes
    tangnew->Add(*mmatrixa,false,-beta,1.0);
    tangnew->Add(*dmatrixa,false,+beta,1.0);
    
    // already finished
    tangnew->Complete();
    
    // modify right hand side -------------------------------------------- 
    // slave and master temperature vectors  
    Teuchos::RCP<Epetra_Vector> fs, fm;
    LINALG::SplitVector(*problemrowmap,*temp,sdofs,fs,mdofs,fm);
    
    // add existing feff 
    feffnew->Update(1.0,*feff,1.0);
    
    // add B.T (master row)
    Teuchos::RCP <Epetra_Vector> BdotTemp = Teuchos::rcp(new Epetra_Vector(*masterpar));
    bmatrixa->Multiply(false,*fm,*BdotTemp);
    Teuchos::RCP<Epetra_Vector> BdotTempexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*BdotTemp,*BdotTempexp);
    feffnew->Update(+beta,*BdotTempexp,1.0);
    
    // add M(T).T (master row)
    Teuchos::RCP <Epetra_Vector> MTdotTemp = Teuchos::rcp(new Epetra_Vector(*masterpar));
    mmatrixat->Multiply(true,*fs,*MTdotTemp);
    Teuchos::RCP<Epetra_Vector> MTdotTempexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*MTdotTemp,*MTdotTempexp);
    feffnew->Update(-beta,*MTdotTempexp,1.0);

    // add M.T (slave row)
    Teuchos::RCP <Epetra_Vector> MdotTemp = Teuchos::rcp(new Epetra_Vector(*adofs));
    mmatrixa->Multiply(false,*fm,*MdotTemp);
    Teuchos::RCP<Epetra_Vector> MdotTempexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*MdotTemp,*MdotTempexp);
    feffnew->Update(-beta,*MdotTempexp,1.0);

    // add D.T (slave row)
    Teuchos::RCP <Epetra_Vector> DdotTemp = Teuchos::rcp(new Epetra_Vector(*adofs));
    dmatrixa->Multiply(false,*fs,*DdotTemp);
    Teuchos::RCP<Epetra_Vector> DdotTempexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*DdotTemp,*DdotTempexp);
    feffnew->Update(+beta,*DdotTempexp,1.0);
    
    // mechanical dissipation master side
    Teuchos::RCP<Epetra_Vector> mechdissratemasterexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*mechdissratemaster,*mechdissratemasterexp);
    feffnew->Update(-(1-delta),*mechdissratemasterexp,1.0);

    // mechanical dissipation slave side
    Teuchos::RCP<Epetra_Vector> mechdissrateslaveexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*mechdissrateslave,*mechdissrateslaveexp);
    feffnew->Update(-delta,*mechdissrateslaveexp,1.0);
    
    // already finished
  }
  else
  {  
    //********************************************************************
    // 2. Modification with the use of Lagrange Multipliers
    //********************************************************************
   
    // assemble Matrix A (in addition to D and M)
    Teuchos::RCP<LINALG::SparseMatrix> amatrix = Teuchos::rcp(new LINALG::SparseMatrix(*sdofs,10));
    TransformA(*amatrix,sdofs);
    
    // Fill Complete
    amatrix->Complete();
  
    // active part of dmatrix and mmatrix
    Teuchos::RCP<Epetra_Map> tmp;
    Teuchos::RCP<LINALG::SparseMatrix> dmatrixa,mmatrixa,amatrixa,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    LINALG::SplitMatrix2x2(dmatrix,adofs,idofs,adofs,idofs,dmatrixa,tmp1,tmp2,tmp3);
    LINALG::SplitMatrix2x2(mmatrix,adofs,idofs,mdofs,tmp,mmatrixa,tmp4,tmp5,tmp6);
    LINALG::SplitMatrix2x2(amatrix,adofs,idofs,sdofs,tmp,amatrixa,tmp4,tmp5,tmp6);
  
    // matrices from linearized thermal contact condition
    Teuchos::RCP<LINALG::SparseMatrix>  thermcontLM = Teuchos::rcp(new LINALG::SparseMatrix(*adofs,3));
    Teuchos::RCP<LINALG::SparseMatrix>  thermcontTEMP = Teuchos::rcp(new LINALG::SparseMatrix(*adofs,3));
    Teuchos::RCP<Epetra_Vector>         thermcontRHS = LINALG::CreateVector(*adofs,true);
  
    // assemble thermal contact condition
    AssembleThermContCondition(*thermcontLM,*thermcontTEMP,*thermcontRHS,*dmatrixa,*mmatrixa,*amatrixa,adofs,mdofs,temp,dt);
  
    // complete the matrices
    thermcontLM->Complete(*sdofs,*adofs);
    thermcontTEMP->Complete(*smdofs,*adofs);
    
    // store thermcondLM (needed for monolithic tsi with contact)
    thermcondLM_ = thermcontLM;
    
    // assemble mechanical dissipation for master side
    Teuchos::RCP<Epetra_Vector> mechdissrate;
    if(cmtman_->GetStrategy().Friction())
    {  
      mechdissrate = LINALG::CreateVector(*mdofs,true);
      AssembleMechDissMaster(*mechdissrate,dt);
    }  
  
    /**********************************************************************/
    /* Modification of the stiff matrix and rhs towards thermo contact    */
    /**********************************************************************/
  
    /**********************************************************************/
    /* Create inv(D)                                                      */
    /**********************************************************************/
    Teuchos::RCP<LINALG::SparseMatrix> invd = Teuchos::rcp(new LINALG::SparseMatrix(*dmatrix));
    Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*sdofs,true);
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
  
    // do the multiplication M^ = inv(D) * M
    Teuchos::RCP<LINALG::SparseMatrix> mhatmatrix;
    mhatmatrix = LINALG::MLMultiply(*invd,false,*mmatrix,false,false,false,true);
  
    /**********************************************************************/
    /* Split tang into 3x3 block matrix                                  */
    /**********************************************************************/
    // we want to split k into 3 groups s,m,n = 9 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kss, ksm, ksn, kms, kmm, kmn, kns, knm, knn;
  
    // temporarily we need the blocks ksmsm, ksmn, knsm
    // (FIXME: because a direct SplitMatrix3x3 is still missing!)
    Teuchos::RCP<LINALG::SparseMatrix> ksmsm, ksmn, knsm;
  
    // some temporary RCPs
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx2;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx3;
  
    // split into slave/master part + structure part
    LINALG::SplitMatrix2x2(tangmatrix,smdofs,ndofs,smdofs,ndofs,ksmsm,ksmn,knsm,knn);
  
    // further splits into slave part + master part
    LINALG::SplitMatrix2x2(ksmsm,sdofs,mdofs,sdofs,mdofs,kss,ksm,kms,kmm);
    LINALG::SplitMatrix2x2(ksmn,sdofs,mdofs,ndofs,tempmap,ksn,tempmtx1,kmn,tempmtx2);
    LINALG::SplitMatrix2x2(knsm,ndofs,tempmap,sdofs,mdofs,kns,knm,tempmtx1,tempmtx2);
  
    /**********************************************************************/
    /* Split feff into 3 subvectors                                       */
    /**********************************************************************/
    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fs, fm, fn;
  
    // temporarily we need the group sm
    Teuchos::RCP<Epetra_Vector> fsm;
  
    // do the vector splitting smn -> sm+n -> s+m+n
    LINALG::SplitVector(*problemrowmap,*feff,smdofs,fsm,ndofs,fn);
    LINALG::SplitVector(*smdofs,*fsm,sdofs,fs,mdofs,fm);
  
    /**********************************************************************/
    /* Split slave quantities into active / inactive                      */
    /**********************************************************************/
    // we want to split kssmod into 2 groups a,i = 4 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kaa, kai, kia, kii;
  
    // we want to split ksn / ksm / kms into 2 groups a,i = 2 blocks
    Teuchos::RCP<LINALG::SparseMatrix> kan, kin, kam, kim, kma, kmi;
  
    // do the splitting
    LINALG::SplitMatrix2x2(kss,adofs,idofs,adofs,idofs,kaa,kai,kia,kii);
    LINALG::SplitMatrix2x2(ksn,adofs,idofs,ndofs,tempmap,kan,tempmtx1,kin,tempmtx2);
    LINALG::SplitMatrix2x2(ksm,adofs,idofs,mdofs,tempmap,kam,tempmtx1,kim,tempmtx2);
    LINALG::SplitMatrix2x2(kms,mdofs,tempmap,adofs,idofs,kma,kmi,tempmtx1,tempmtx2);
  
    // we want to split fsmod into 2 groups a,i
    Teuchos::RCP<Epetra_Vector> fa = Teuchos::rcp(new Epetra_Vector(*adofs));
    Teuchos::RCP<Epetra_Vector> fi = Teuchos::rcp(new Epetra_Vector(*idofs));
  
    // do the vector splitting s -> a+i
    LINALG::SplitVector(*sdofs,*fs,adofs,fa,idofs,fi);
  
    // abbreviations for active and inactive set
    int aset = adofs->NumGlobalElements();
    int iset = idofs->NumGlobalElements();
  
    // active part of invd and mhatmatrix
    Teuchos::RCP<Epetra_Map> tmpmap;
    Teuchos::RCP<LINALG::SparseMatrix> invda,mhata;
    LINALG::SplitMatrix2x2(invd,sdofs,tmpmap,adofs,idofs,invda,tmp1,tmp2,tmp3);
    LINALG::SplitMatrix2x2(mhatmatrix,adofs,idofs,mdofs,tmpmap,mhata,tmp1,tmp2,tmp3);
    
    // store some stuff for static condensation of LM
    fs_   = fs;
    invd_ = invd;
    ksn_  = ksn;
    ksm_  = ksm;
    kss_  = kss;
  
    /**********************************************************************/
    /* Build the final K and f blocks                                     */
    /**********************************************************************/
    // knn: nothing to do
  
    // knm: nothing to do
  
    // kns: nothing to do
  
    // kmn: add T(mbaractive)*kan
    Teuchos::RCP<LINALG::SparseMatrix> kmnmod = Teuchos::rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmnmod->Add(*kmn,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmnadd = LINALG::MLMultiply(*mhata,true,*kan,false,false,false,true);
    kmnmod->Add(*kmnadd,false,1.0,1.0);
    kmnmod->Complete(kmn->DomainMap(),kmn->RowMap());
  
    // kmm: add T(mbaractive)*kam
    Teuchos::RCP<LINALG::SparseMatrix> kmmmod = Teuchos::rcp(new LINALG::SparseMatrix(*mdofs,100));
    kmmmod->Add(*kmm,false,1.0,1.0);
    Teuchos::RCP<LINALG::SparseMatrix> kmmadd = LINALG::MLMultiply(*mhata,true,*kam,false,false,false,true);
    kmmmod->Add(*kmmadd,false,1.0,1.0);
    kmmmod->Complete(kmm->DomainMap(),kmm->RowMap());
  
    // kmi: add T(mbaractive)*kai
    Teuchos::RCP<LINALG::SparseMatrix> kmimod;
    if (iset)
    {
      kmimod = Teuchos::rcp(new LINALG::SparseMatrix(*mdofs,100));
      kmimod->Add(*kmi,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmiadd = LINALG::MLMultiply(*mhata,true,*kai,false,false,false,true);
      kmimod->Add(*kmiadd,false,1.0,1.0);
      kmimod->Complete(kmi->DomainMap(),kmi->RowMap());
    }
  
    // kmi: add T(mbaractive)*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kmamod;
    if (aset)
    {
      kmamod = Teuchos::rcp(new LINALG::SparseMatrix(*mdofs,100));
      kmamod->Add(*kma,false,1.0,1.0);
      Teuchos::RCP<LINALG::SparseMatrix> kmaadd = LINALG::MLMultiply(*mhata,true,*kaa,false,false,false,true);
      kmamod->Add(*kmaadd,false,1.0,1.0);
      kmamod->Complete(kma->DomainMap(),kma->RowMap());
    }
  
    // kan: thermcontlm*invd*kan
    Teuchos::RCP<LINALG::SparseMatrix> kanmod;
    if (aset)
    {
      kanmod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
      kanmod = LINALG::MLMultiply(*kanmod,false,*kan,false,false,false,true);
      kanmod->Complete(kan->DomainMap(),kan->RowMap());
    }
  
    // kam: thermcontlm*invd*kam
    Teuchos::RCP<LINALG::SparseMatrix> kammod;
    if (aset)
    {
      kammod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
      kammod = LINALG::MLMultiply(*kammod,false,*kam,false,false,false,true);
      kammod->Complete(kam->DomainMap(),kam->RowMap());
    }
  
    // kai: thermcontlm*invd*kai
    Teuchos::RCP<LINALG::SparseMatrix> kaimod;
    if (aset && iset)
    {
      kaimod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
      kaimod = LINALG::MLMultiply(*kaimod,false,*kai,false,false,false,true);
      kaimod->Complete(kai->DomainMap(),kai->RowMap());
    }
  
    // kaa: thermcontlm*invd*kaa
    Teuchos::RCP<LINALG::SparseMatrix> kaamod;
    if (aset)
    {
      kaamod = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
      kaamod = LINALG::MLMultiply(*kaamod,false,*kaa,false,false,false,true);
      kaamod->Complete(kaa->DomainMap(),kaa->RowMap());
    }
  
    // Modifications towards rhs
    // FIXGIT: pay attention to genalpha
    // fm: add T(mbaractive)*fa
    Teuchos::RCP<Epetra_Vector> fmmod = Teuchos::rcp(new Epetra_Vector(*mdofs));
    mhata->Multiply(true,*fa,*fmmod);
    fmmod->Update(1.0,*fm,1.0);
  
    // fa: mutliply with thermcontlm
    Teuchos::RCP<Epetra_Vector> famod;
    {
      famod = Teuchos::rcp(new Epetra_Vector(*adofs));
      Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::MLMultiply(*thermcontLM,false,*invda,false,false,false,true);
      temp->Multiply(false,*fa,*famod);
    }
  
    /**********************************************************************/
    /* Global setup of tangnew, feffnew (including contact)              */
    /**********************************************************************/
  
    // add n submatrices to tangnew
    tangnew->Add(*knn,false,1.0,1.0);
    tangnew->Add(*knm,false,1.0,1.0);
    tangnew->Add(*kns,false,1.0,1.0);
  
    // add m submatrices to tangnew
    tangnew->Add(*kmnmod,false,1.0,1.0);
    tangnew->Add(*kmmmod,false,1.0,1.0);
    if (iset) tangnew->Add(*kmimod,false,1.0,1.0);
    if (aset) tangnew->Add(*kmamod,false,1.0,1.0);
  
    // add i submatrices to tangnew
    if (iset) tangnew->Add(*kin,false,1.0,1.0);
    if (iset) tangnew->Add(*kim,false,1.0,1.0);
    if (iset) tangnew->Add(*kii,false,1.0,1.0);
    if (iset) tangnew->Add(*kia,false,1.0,1.0);
  
    // add a submatrices to tangnew
    if (aset) tangnew->Add(*kanmod,false,-1.0,1.0);
    if (aset) tangnew->Add(*kammod,false,-1.0,1.0);
    if (aset && iset) tangnew->Add(*kaimod,false,-1.0,1.0);
    if (aset) tangnew->Add(*kaamod,false,-1.0,1.0);
  
    // add n subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fnexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*fn,*fnexp);
    feffnew->Update(1.0,*fnexp,1.0);
  
    // add m subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fmmodexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*fmmod,*fmmodexp);
    feffnew->Update(1.0,*fmmodexp,1.0);

    // add mechanical dissipation to feffnew, only in frictional case
    if (cmtman_->GetStrategy().Friction())
    {  
      Teuchos::RCP<Epetra_Vector> mechdissrateexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
      LINALG::Export(*mechdissrate,*mechdissrateexp);
      feffnew->Update(+1.0,*mechdissrateexp,1.0);
    }
      
    // add i subvector to feffnew
    Teuchos::RCP<Epetra_Vector> fiexp;
    if (iset)
    {
      fiexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
      LINALG::Export(*fi,*fiexp);
      feffnew->Update(1.0,*fiexp,1.0);
    }
  
    // add a subvector to feffnew
    Teuchos::RCP<Epetra_Vector> famodexp;
    if (aset)
    {
      famodexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
      LINALG::Export(*famod,*famodexp);
      feffnew->Update(-1.0,*famodexp,+1.0);
    }
  
    // add linearized thermo contact condition
    tangnew->Add(*thermcontTEMP,false,+1.0,+1.0);
  
    // add rhs of thermal contact condition to feffnew
    Teuchos::RCP<Epetra_Vector> thermcontRHSexp = Teuchos::rcp(new Epetra_Vector(*problemrowmap));
    LINALG::Export(*thermcontRHS,*thermcontRHSexp);
    feffnew->Update(+1.0,*thermcontRHSexp,1.0);
  
    // FillComplete tangnew (square)
    tangnew->Complete();
  }
  
  /**********************************************************************/
  /* Replace tang and feff by tangnew and feffnew                     */
  /**********************************************************************/
  tang = tangnew;
  feff = feffnew;

  // leave this place
  return;
}

/*----------------------------------------------------------------------*
 | convert maps form structure dofs to thermo dofs            mgit 04/10 |
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::ConvertMaps(Teuchos::RCP<Epetra_Map>& slavedofs,
                                        Teuchos::RCP<Epetra_Map>& activedofs,
                                        Teuchos::RCP<Epetra_Map>& masterdofs)
{

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes/dofs
    const Teuchos::RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // define local variables
    int slavecountnodes = 0;
    std::vector<int> myslavegids(slavenodes->NumMyElements());

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != comm_->MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      myslavegids[slavecountnodes] = (discretstruct_->Dof(1,node))[0];
      ++slavecountnodes;
    }

    // resize the temporary vectors
    myslavegids.resize(slavecountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gslavecountnodes;
    comm_->SumAll(&slavecountnodes,&gslavecountnodes,1);

    // create active node map and active dof map
    slavedofs = Teuchos::rcp(new Epetra_Map(gslavecountnodes,slavecountnodes,&myslavegids[0],0,*comm_));

    // active nodes/dofs
    const Teuchos::RCP<Epetra_Map> activenodes = interface[m]->ActiveNodes();

    // define local variables
    int countnodes = 0;
    std::vector<int> mynodegids(activenodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<activenodes->NumMyElements();++i)
    {
      int gid = activenodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != comm_->MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mynodegids[countnodes] = (discretstruct_->Dof(1,node))[0];
      ++countnodes;
    }

    // resize the temporary vectors
    mynodegids.resize(countnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gcountnodes;
    comm_->SumAll(&countnodes,&gcountnodes,1);

    // create active node map and active dof map
    activedofs = Teuchos::rcp(new Epetra_Map(gcountnodes,countnodes,&mynodegids[0],0,*comm_));

    // master nodes/dofs
    const Teuchos::RCP<Epetra_Map> masternodes = interface[m]->MasterRowNodes();

    // define local variables
    int mastercountnodes = 0;
    std::vector<int> mymastergids(masternodes->NumMyElements());

    // loop over all active nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = discretstruct_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != comm_->MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      mymastergids[mastercountnodes] = (discretstruct_->Dof(1,node))[0];
      ++mastercountnodes;
    }

    // resize the temporary vectors
    mymastergids.resize(mastercountnodes);

    // communicate countnodes, countdofs, countslipnodes and countslipdofs among procs
    int gmastercountnodes;
    comm_->SumAll(&mastercountnodes,&gmastercountnodes,1);

    // create active node map and active dof map
    masterdofs = Teuchos::rcp(new Epetra_Map(gmastercountnodes,mastercountnodes,&mymastergids[0],0,*comm_));
  }
  return;
}

/*----------------------------------------------------------------------*
 | transform mortar matrices in thermo dofs                   mgit 04/10 |
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::TransformDM(LINALG::SparseMatrix& dmatrix,
                                       LINALG::SparseMatrix& mmatrix,
                                       Teuchos::RCP<Epetra_Map>& slavedofs,
                                       Teuchos::RCP<Epetra_Map>& masterdofs)
{
  // static cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // dimension of the problem
  int dim = cstrategy.Dim();
  if (dim==2)
    dserror("In THR::TimIntImpl::TransformDM: Thermal problems only in 3D so far");
  
  int myrow;
  int numentries;
  
  /****************************************************** D-matrix ******/
  
  // mortar matrix D in structural dofs 
  Teuchos::RCP<Epetra_CrsMatrix> dstruct = (cstrategy.DMatrix())->EpetraMatrix();
  
  // row and column map of mortar matrices
  const Epetra_Map& rowmap = dstruct->RowMap();
  const Epetra_Map& colmap = dstruct->ColMap();
  
  // set of every dim-th dof (row map)
  std::set<int> rowsetred;
  for (myrow=0; myrow<dstruct->NumMyRows(); myrow=myrow+dim)
    rowsetred.insert(rowmap.GID(myrow));
  
  // map of every dim-th dof (row map)
  Teuchos::RCP<Epetra_Map> rowmapred = LINALG::CreateMap(rowsetred,*comm_);
  
  // mortar matrices in reduced structural dofs
  // this means no loss of information
  Teuchos::RCP<LINALG::SparseMatrix> dstructred = Teuchos::rcp(new LINALG::SparseMatrix(*rowmapred,10)) ;

  // loop over all rows of mortar matrix D 
  for (myrow=0; myrow<dstruct->NumMyRows(); myrow=myrow+dim)
  {
    double *Values;
    int *Indices;

    int err = dstruct->ExtractMyRowView(myrow, numentries, Values, Indices);
    if (err)
      dserror("ExtractMyRowView failed: err=%d", err);

     // row
     int row = rowmap.GID(myrow);
    
    // loop over entries of the row 
    for (int i=0;i<numentries; ++i)
    {
      // col and val
      int col  = colmap.GID(Indices[i]);
      double  val = Values[i];
      
      // assembly
      dstructred->Assemble(val,row,col); 
    }
  }
  
  // complete the matrix
  dstructred->Complete();
  
  // transform reduced structural D matrix in thermo dofs
  dmatrix=*(MORTAR::MatrixRowColTransformGIDs(dstructred,slavedofs,slavedofs));
  
//  dmatrix.Assemble(2.866034e-02,792,792); 
//  dmatrix.Assemble(2.866034e-02,810,810); 
//  dmatrix.Assemble(2.866034e-02,834,834); 
//  dmatrix.Assemble(2.866034e-02,843,843); 

  /****************************************************** M-matrix ******/
  
  // mortar matrix M in structural dofs 
  Teuchos::RCP<Epetra_CrsMatrix> mstruct = (cstrategy.MMatrix())->EpetraMatrix();
  
  // row and column map of mortar matrix M
  const Epetra_Map& rowmapm = mstruct->RowMap();
  const Epetra_Map& colmapm = mstruct->ColMap();
  const Epetra_Map& domainmapm = mstruct->DomainMap();
  
  // set of every dim-th dof (row map)
  rowsetred.clear();
  for (myrow=0; myrow<mstruct->NumMyRows(); myrow=myrow+dim)
    rowsetred.insert(rowmapm.GID(myrow));

  // set of every dim-th dof (domain map)
  std::set<int> domainsetred;
    for (int mycol=0; mycol<domainmapm.NumMyElements(); mycol=mycol+dim)
      domainsetred.insert(domainmapm.GID(mycol));
    
  // map of every dim-th dof (row map)
  rowmapred = LINALG::CreateMap(rowsetred,*comm_);
  
  // map of every dim-th dof (domain map)
  Teuchos::RCP<Epetra_Map> domainmapred = LINALG::CreateMap(domainsetred,*comm_);
  
  // mortar matrix M in reduced structural dofs
  // this means no loss of information
  Teuchos::RCP<LINALG::SparseMatrix> mstructred = Teuchos::rcp(new LINALG::SparseMatrix(*rowmapred,100)) ;

  // loop over all rows of mortar matrix M 
  for (myrow=0; myrow<mstruct->NumMyRows(); myrow=myrow+dim)
  {
    double *Values;
    int *Indices;

    int err = mstruct->ExtractMyRowView(myrow, numentries, Values, Indices);
    if (err)
      dserror("ExtractMyRowView failed: err=%d", err);

     // row
     int row = rowmapm.GID(myrow);
     
    // loop over entries of the row 
    for (int i=0;i<numentries; ++i)
    {
      // col and val
      int col  = colmapm.GID(Indices[i]);
      double  val = Values[i];
      
      // assembly
      mstructred->Assemble(val,row,col); 
    }
  }
  
  // complete the matrix
  mstructred->Complete(*domainmapred,*rowmapred);
  
  // transform reduced structural M matrix in thermo dofs
  mmatrix=*(MORTAR::MatrixRowColTransformGIDs(mstructred,slavedofs,masterdofs));
  
//  mmatrix.Assemble(4.597135e-03,792,           851);            
//  mmatrix.Assemble(8.293508e-03,792,           863);            
//  mmatrix.Assemble(-1.410774e-03,792,           871);          
//  mmatrix.Assemble(8.293508e-03, 792,           903);           
//  mmatrix.Assemble( 1.495921e-02,792,           909 );           
//  mmatrix.Assemble(-2.547689e-03,792,           913  );        
//  mmatrix.Assemble(-1.410774e-03,792,           931   );        
//  mmatrix.Assemble(-2.547689e-03,792,           937    );       
//  mmatrix.Assemble(4.339013e-04, 792,           941     );       
//  mmatrix.Assemble(-1.410774e-03,810,           871    );       
//  mmatrix.Assemble(8.293508e-03, 810,           879   );        
//  mmatrix.Assemble(4.597135e-03,  810,           887   );         
//  mmatrix.Assemble(-2.547689e-03, 810,           913   );        
//  mmatrix.Assemble( 1.495921e-02, 810,           917   );        
//  mmatrix.Assemble( 8.293508e-03, 810,           921   );         
//  mmatrix.Assemble(4.339013e-04,  810,           941   );         
//  mmatrix.Assemble(-2.547689e-03, 810,           945   );        
//  mmatrix.Assemble(-1.410774e-03, 810,           949   );        
//  mmatrix.Assemble(-1.410774e-03, 834,           931   );        
//  mmatrix.Assemble(-2.547689e-03, 834,           937   );        
//  mmatrix.Assemble( 4.339013e-04, 834,           941   );        
//  mmatrix.Assemble( 8.293508e-03, 834,           959   );        
//  mmatrix.Assemble( 1.495921e-02, 834,           965   );         
//  mmatrix.Assemble(-2.547689e-03, 834,           969   );        
//  mmatrix.Assemble( 4.597135e-03, 834,           987   );         
//  mmatrix.Assemble( 8.293508e-03, 834,           993   );         
//  mmatrix.Assemble(-1.410774e-03,  834,           997   );        
//  mmatrix.Assemble( 4.339013e-04, 843,           941   );        
//  mmatrix.Assemble( -2.547689e-03, 843,           945   );        
//  mmatrix.Assemble( -1.410774e-03, 843,           949  );         
//  mmatrix.Assemble( -2.547689e-03, 843,           969  );         
//  mmatrix.Assemble( 1.495921e-02,  843,           973  );          
//  mmatrix.Assemble( 8.293508e-03,  843,           977  );         
//  mmatrix.Assemble( -1.410774e-03, 843,           997  );         
//  mmatrix.Assemble( 8.293508e-03, 843,          1001  );         
//  mmatrix.Assemble( 4.597135e-03,  843,          1005  );          

  return;
}

/*----------------------------------------------------------------------*
 | transform matrix A in thermo dofs                          mgit 10/10 |
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::TransformA(LINALG::SparseMatrix& amatrix,
                                       Teuchos::RCP<Epetra_Map>& slavedofs)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // dimension of the problem
  int dim = cstrategy.Dim();
  if (dim==2)
    dserror("ERROR: THR::TimIntImpl::TransformA: Thermal problems only in 3D so far");
  
  int myrow;
  int numentries;
  
  /****************************************************** A-matrix ******/
  
  // matrix A in structural dofs 
  Teuchos::RCP<Epetra_CrsMatrix> astruct = (cstrategy.AMatrix())->EpetraMatrix();
  
  // row and column map of mortar matrices
  const Epetra_Map& rowmap = astruct->RowMap();
  const Epetra_Map& colmap = astruct->ColMap();
  
  // set of every dim-th dof (row map)
  std::set<int> rowsetred;
  for (myrow=0; myrow<astruct->NumMyRows(); myrow=myrow+dim)
    rowsetred.insert(rowmap.GID(myrow));
  
  // map of every dim-th dof (row map)
  Teuchos::RCP<Epetra_Map> rowmapred = LINALG::CreateMap(rowsetred,*comm_);
  
  // mortar matrices in reduced structural dofs
  // this means no loss of information
  Teuchos::RCP<LINALG::SparseMatrix> astructred = Teuchos::rcp(new LINALG::SparseMatrix(*rowmapred,10)) ;

  // loop over all rows of matrix A
  for (myrow=0; myrow<astruct->NumMyRows(); myrow=myrow+dim)
  {
    double *Values;
    int *Indices;

    int err = astruct->ExtractMyRowView(myrow, numentries, Values, Indices);
    if (err)
      dserror("ExtractMyRowView failed: err=%d", err);

     // row
     int row = rowmap.GID(myrow);
    
    // loop over entries of the row 
    for (int i=0;i<numentries; ++i)
    {
      // col and val
      int col  = colmap.GID(Indices[i]);
      double  val = Values[i];
      
      // assembly
      astructred->Assemble(val,row,col); 
    }
  }
  
  // complete the matrix
  astructred->Complete();
  
  // transform reduced structural A matrix in thermo dofs
  amatrix=*(MORTAR::MatrixRowColTransformGIDs(astructred,slavedofs,slavedofs));

  return;
}

/*----------------------------------------------------------------------*
 | assemble matrix B in thermo dofs                           mgit 10/10 |
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::AssembleB(LINALG::SparseMatrix& bmatrix)

{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // master nodes (full map)
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(interface[m]->MasterRowNodes()));

    // loop over all slave nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = 0;
      if(comm_->MyPID()==cnode->Owner())
        rowtemp = discretstruct_->Dof(1,nodeges)[0];

      /************************************************** B-matrix ******/
      std::set<int> bnodes;
      int mastergid=0;
      std::set<int>::iterator mcurr;
      int size = 0;
      std::vector<map<int,double> > bmap;

      bmap = cnode->GetB();
      bnodes = cnode->GetBNodes();
      size = bnodes.size();
      
      // flag if proc has B-matrix entries for this master node
      bool hasentries = false;
      if(size!=0)
        hasentries=true;
 
      if(hasentries)
        mcurr = bnodes.begin();
  
      // sum all entries of size
      int sizeglobal = 0;
      comm_->SumAll(&size,&sizeglobal,1);
      
      // the node should have entries only from one proc
      if(size!=sizeglobal)
      {  
        if (abs(size)>1e-12)
          dserror ("Error in AssembleB: Entries from more than one proc");
      }
      
      // loop over all according master nodes
      for (int l=0;l<sizeglobal;++l)
      {
        if (hasentries)
          mastergid=*mcurr;
        
        // global id of master node
        int masterglobal = 0;

        // sum all entries of mastergid
        comm_->SumAll(&mastergid,&masterglobal,1);
        
        DRT::Node* mnode = (interface[m]->Discret()).gNode(masterglobal);
        DRT::Node* mnodeges = discretstruct_->gNode(masterglobal);
        
        // temperature and displacement dofs
        int coltemp = 0;
        int coldis = 0;
        if(comm_->MyPID()==mnode->Owner())
        {
          CONTACT::CoNode* cmnode = static_cast<CONTACT::CoNode*>(mnode);
          coltemp = discretstruct_->Dof(1,mnodeges)[0];
          coldis = (cmnode->Dofs())[0];
        }

        // communicate temperature and displacement dof
        comm_->Broadcast(&coltemp,1,mnode->Owner());
        comm_->Broadcast(&coldis,1,mnode->Owner());

        // value from procs
        double val = 0;
        
        if(hasentries)
          val = bmap[0][coldis];
        
        // communicate the value and sum in valglobal
        double valglobal = 0;
        comm_->SumAll(&val,&valglobal,1);
         
        // do the assembly
        if (comm_->MyPID()==cnode->Owner() and (abs(valglobal)>1e-12) )
          bmatrix.Assemble(valglobal, rowtemp, coltemp);
        
        if(hasentries)
          ++mcurr;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mechanical dissipation for master nodes            mgit 08/10|
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::AssembleMechDissMaster(Epetra_Vector& mechdissrate,double dt)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");
  
  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // loop over master full nodes
    // master nodes are redundant on all procs and the entry of the
    // mechanical dissipation lies on proc which did the evaluation of 
    // mortar integrals and mechanical dissipation 

    // master nodes
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(interface[m]->MasterRowNodes()));

    // loop over all masternodes nodes of the interface
    for (int i=0;i<masternodes->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // mechanical dissipation to be assembled     
      double mechdissglobal = 0;
      
      // mechanical dissipation on proc
      double mechdissproc = 1/dt*cnode->MechDiss();
            
      // sum all entries to mechdissglobal
      comm_->SumAll(&mechdissproc,&mechdissglobal,1);
      
      // owner of master node does the assembly
      if(comm_->MyPID()==cnode->Owner())
      {
        // row dof of temperature
        DRT::Node* nodeges = discretstruct_->gNode(gid);
        int rowtemp = discretstruct_->Dof(1,nodeges)[0];

        Epetra_SerialDenseVector mechdissiprate(1);
        std::vector<int> dof(1);
        std::vector<int> owner(1);

        mechdissiprate(0) = mechdissglobal;
        dof[0] = rowtemp;
        owner[0] = cnode->Owner();

        // do assembly
        if(abs(mechdissiprate(0))>1e-12)
          LINALG::Assemble(mechdissrate, mechdissiprate, dof, owner);
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble mechanical dissipation for slave nodes            mgit 11/10|
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::AssembleMechDissSlave(Epetra_Vector& mechdissrate,double dt)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");
  
  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
     // slave nodes
    const Teuchos::RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
      DRT::Node* nodeges = discretstruct_->gNode(gid);

      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);

      // row dof of temperature
      int rowtemp = discretstruct_->Dof(1,nodeges)[0];

      Epetra_SerialDenseVector mechdissiprate(1);
      std::vector<int> dof(1);
      std::vector<int> owner(1);

      mechdissiprate(0) = 1/dt*cnode->MechDiss();
      dof[0] = rowtemp;
      owner[0] = cnode->Owner();

      // do assembly
      if(abs(mechdissiprate(0))>1e-12)
        LINALG::Assemble(mechdissrate, mechdissiprate, dof, owner);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assemble the thermal contact conditions for slave nodes    mgit 04/10 |
 *----------------------------------------------------------------------*/

void THR::ThermoContactMan::AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                                       LINALG::SparseMatrix& thermcontTEMP,
                                                       Epetra_Vector& thermcontRHS,
                                                       LINALG::SparseMatrix& dmatrix,
                                                       LINALG::SparseMatrix& mmatrix,
                                                       LINALG::SparseMatrix& amatrix,
                                                       Teuchos::RCP<Epetra_Map> activedofs,
                                                       Teuchos::RCP<Epetra_Map> masterdofs,
                                                       Teuchos::RCP<Epetra_Vector>& temp,
                                                       double dt)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet and for one heat
  // transfer coefficient
  // FIXGIT: The heat transfer coefficient should be a condition on
  // the single interfaces!!
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::AssembleThermContCondition: Only for one interface yet.");

  // heat transfer coefficient for slave and master surface
  double heattranss = interface[0]->IParams().get<double>("HEATTRANSSLAVE");
  double heattransm = interface[0]->IParams().get<double>("HEATTRANSMASTER");

  if (heattranss <= 0 or heattransm <= 0)
   dserror("Error: Choose realistic heat transfer parameter");

  double beta = heattranss*heattransm/(heattranss+heattransm);
  double delta = heattranss/(heattranss+heattransm);

  // with respect to Lagrange multipliers
  if (cmtman_->GetStrategy().Friction())
    thermcontLM.Add(amatrix,false,1.0,1.0);
  else
    thermcontLM.Add(dmatrix,false,1.0,1.0);
    
  // with respect to temperature
  thermcontTEMP.Add(dmatrix,false,-beta,1.0);
  thermcontTEMP.Add(mmatrix,false,+beta,1.0);

  Teuchos::RCP<Epetra_Vector> fa, fm, rest1, rest2;

  // row map of thermal problem
  Teuchos::RCP<Epetra_Map> problemrowmap = Teuchos::rcp(new Epetra_Map(*(discretthermo_->DofRowMap(0))));

  LINALG::SplitVector(*problemrowmap,*temp,activedofs,fa,masterdofs,fm);

  Teuchos::RCP <Epetra_Vector> DdotTemp = Teuchos::rcp(new Epetra_Vector(*activedofs));
  dmatrix.Multiply(false,*fa,*DdotTemp);
  thermcontRHS.Update(+beta,*DdotTemp,1.0);

  Teuchos::RCP <Epetra_Vector> MdotTemp = Teuchos::rcp(new Epetra_Vector(*activedofs));
  mmatrix.Multiply(false,*fm,*MdotTemp);
  thermcontRHS.Update(-beta,*MdotTemp,1.0);
  
//  // loop over active dofs
//  for (int i=0; i<cstrategy.ActiveRowNodes()->NumMyElements(); ++i)
//  {
//    int gid = cstrategy.ActiveRowNodes()->GID(i);
//    //cout << "GID " << gid << endl;
//    
//     DRT::Node* node = discretstruct_->gNode(gid);
//     
//     int row = (discretstruct_->Dof(1,node))[0];
//     
//     //cout << "ROW " << row << endl;
//     
//     int locid = (fa->Map()).LID(row);
//   
//     thermcontRHS[locid] = (-1)*(*fa)[locid]+5;
//     thermcontTEMP.Assemble(-1,row,row);
//  }

  // add mechanical dissipation only in the frictional case
  if (cmtman_->GetStrategy().Friction())
  {  
    // loop over all interfaces
    for (int m=0; m<(int)interface.size(); ++m)
    {
      // slave nodes (full map)
      const Teuchos::RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();
  
      // loop over all slave nodes of the interface
      for (int i=0;i<slavenodes->NumMyElements();++i)
      {
        int gid = slavenodes->GID(i);
        DRT::Node* node    = (interface[m]->Discret()).gNode(gid);
        DRT::Node* nodeges = discretstruct_->gNode(gid);
  
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CONTACT::FriNode* cnode = static_cast<CONTACT::FriNode*>(node);
  
        // row dof of temperature
        int rowtemp = 0;
        if(comm_->MyPID()==cnode->Owner())
          rowtemp = discretstruct_->Dof(1,nodeges)[0];
  
        Epetra_SerialDenseVector mechdissiprate(1);
        std::vector<int> dof(1);
        std::vector<int> owner(1);
  
        mechdissiprate(0) = (-1)*delta/dt*cnode->MechDiss();
        dof[0] = rowtemp;
        owner[0] = cnode->Owner();
  
        // do assembly
        if(abs(mechdissiprate(0))>1e-12 and cnode->Active())
          LINALG::Assemble(thermcontRHS, mechdissiprate, dof, owner);
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Initialize thermal LM                                       mgit 05/11|
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::InitializeThermLM(Teuchos::RCP<Epetra_Map> sthermdofs)
{
 
  // initialize thermal LM
  z_ = Teuchos::rcp(new Epetra_Vector(*sthermdofs));
  
  return;
}

/*----------------------------------------------------------------------*
 | Recovery method                                            mgit 05/11|
 *----------------------------------------------------------------------*/
void THR::ThermoContactMan::RecoverThermLM(Teuchos::RCP<Epetra_Vector> tempi)
{
 
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  // static cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);
 
  if (!cstrategy.IsInContact() && !cstrategy.WasInContact() && !cstrategy.WasInContactLastTimeStep())
    return;

  // necessary maps 
  // convert maps (from structure discretization to thermo discretization)
  // slave-, active-, inactive-, master-, activemaster-, n- smdofs
  Teuchos::RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  ConvertMaps (sdofs,adofs,mdofs);
 
  // check if contact contributions are present,
  // speed up in the case of no contact
  if(!(adofs->NumGlobalElements()))
    return;

  smdofs = LINALG::MergeMap(sdofs,mdofs,false);
 
  // row map of thermal problem
  Teuchos::RCP<Epetra_Map> problemrowmap = Teuchos::rcp(new Epetra_Map(*(discretthermo_->DofRowMap(0))));
  
  // split problemrowmap in n+am
  ndofs = LINALG::SplitMap(*problemrowmap,*smdofs);

    // double-check if this is a dual LM system
    //if (shapefcn!=INPAR::MORTAR::shape_dual) dserror("Condensation only for dual LM");
        
    // extract slave temperatures from tempi
    Teuchos::RCP<Epetra_Vector> tempis = Teuchos::rcp(new Epetra_Vector(*sdofs));
    if (sdofs->NumGlobalElements()) LINALG::Export(*tempi, *tempis);

   // extract master displacements from disi
    Teuchos::RCP<Epetra_Vector> tempim = Teuchos::rcp(new Epetra_Vector(*mdofs));
    if (mdofs->NumGlobalElements()) LINALG::Export(*tempi, *tempim);

    // extract other displacements from disi
    Teuchos::RCP<Epetra_Vector> tempin = Teuchos::rcp(new Epetra_Vector(*ndofs));
    if (ndofs->NumGlobalElements()) LINALG::Export(*tempi,*tempin);

    // condensation has been performed for active LM only,
    // thus we construct a modified invd matrix here which
    // only contains the active diagonal block
    // (this automatically renders the incative LM to be zero)
    Teuchos::RCP<LINALG::SparseMatrix> invda;
    Teuchos::RCP<Epetra_Map> tempmap;
    Teuchos::RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;
    LINALG::SplitMatrix2x2(invd_,adofs,tempmap,adofs,tempmap,invda,tempmtx1,tempmtx2,tempmtx3);
    Teuchos::RCP<LINALG::SparseMatrix> invdmod = Teuchos::rcp(new LINALG::SparseMatrix(*sdofs,10));
    invdmod->Add(*invda,false,1.0,1.0);
    invdmod->Complete();

    /**********************************************************************/
    /* Update Lagrange multipliers z_n+1                                  */
    /**********************************************************************/

    // full update
    z_->Update(1.0,*fs_,0.0);
    Teuchos::RCP<Epetra_Vector> mod = Teuchos::rcp(new Epetra_Vector(*sdofs));
    kss_->Multiply(false,*tempis,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksm_->Multiply(false,*tempim,*mod);
    z_->Update(-1.0,*mod,1.0);
    ksn_->Multiply(false,*tempin,*mod);
    z_->Update(-1.0,*mod,1.0);
    Teuchos::RCP<Epetra_Vector> zcopy = Teuchos::rcp(new Epetra_Vector(*z_));
    invdmod->Multiply(true,*zcopy,*z_);
    
  return;
}


/*----------------------------------------------------------------------*/
