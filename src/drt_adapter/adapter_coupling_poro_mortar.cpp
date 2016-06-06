/*---------------------------------------------------------------------*/
/*!
\file adapter_coupling_poro_mortar.cpp

\brief coupling adapter for porous meshtying

Masterthesis of h.Willmann under supervision of Anh-Tu Vuong
and Christoph Ager, Originates from ADAPTER::CouplingNonLinMortar

\level 3

\maintainer Anh-Tu Vuong, Christoph Ager

*/
/*---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  includes                                                  ager 10/15|
 *----------------------------------------------------------------------*/
//lib
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

//contact
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_interface.H"

//nurbs
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_knotvector.H"

//linalg
#include "../linalg/linalg_sparsematrix.H"

//inpar
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_structure.H"

//header
#include "adapter_coupling_poro_mortar.H"

/*----------------------------------------------------------------------*
 |  ctor                                                      ager 10/15|
 *----------------------------------------------------------------------*/
ADAPTER::CouplingPoroMortar::CouplingPoroMortar():
  CouplingNonLinMortar(),
  firstinit_(false),
  slavetype_(-1),
  mastertype_(-1)
{
  //empty...
}

/*----------------------------------------------------------------------*
 |  Read Mortar Condition                                     ager 10/15|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::ReadMortarCondition(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<DRT::Discretization>   slavedis,
    std::vector<int>                    coupleddof,
    const std::string&                  couplingcond,
    Teuchos::ParameterList&             input,
    std::map<int, DRT::Node*>& mastergnodes,
    std::map<int, DRT::Node*>& slavegnodes,
    std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
    std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements
    )
{
  //Call Base Class
  CouplingNonLinMortar::ReadMortarCondition(masterdis, slavedis, coupleddof, couplingcond, input, mastergnodes, slavegnodes, masterelements, slaveelements);

  //Set Problem Type to Poro
  input.set<int>("PROBTYPE", INPAR::CONTACT::poro);

  //porotimefac = 1/(theta*dt) --- required for derivation of structural displacements!
  const Teuchos::ParameterList& stru     = DRT::Problem::Instance()->StructuralDynamicParams();
  double porotimefac = 1/(stru.sublist("ONESTEPTHETA").get<double>("THETA") * stru.get<double>("TIMESTEP"));
  input.set<double> ("porotimefac", porotimefac);
  const Teuchos::ParameterList& porodyn     = DRT::Problem::Instance()->PoroelastDynamicParams();
  input.set<bool>("CONTACTNOPEN", DRT::INPUT::IntegralValue<int>(porodyn,"CONTACTNOPEN")); //used in the integrator
  if (!DRT::INPUT::IntegralValue<int>(porodyn,"CONTACTNOPEN"))
    dserror("Set CONTACTNOPEN for Poroelastic meshtying!");
}

/*----------------------------------------------------------------------*
 |  Add Mortar Elements                                        ager 10/15|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::AddMortarElements(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<DRT::Discretization>   slavedis,
    Teuchos::ParameterList&             input,
    std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
    std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements,
    Teuchos::RCP<CONTACT::CoInterface>& interface,
    int numcoupleddof
    )
{
  bool isnurbs = input.get<bool>("NURBS");

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  // We need to determine an element offset to start the numbering of the slave
  // mortar elements AFTER the master mortar elements in order to ensure unique
  // eleIDs in the interface discretization. The element offset equals the
  // overall number of master mortar elements (which is not equal to the number
  // of elements in the field that is chosen as master side).
  //
  // If masterdis==slavedis, the element numbering is right without offset
  int eleoffset = 0;
  if(masterdis.get()!=slavedis.get())
  {
    int nummastermtreles = masterelements.size();
    comm_->SumAll(&nummastermtreles,&eleoffset,1);
  }

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
                new CONTACT::CoElement(ele->Id(), ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), false,isnurbs));

    Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
    if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
    cele->PhysType() = MORTAR::MortarElement::PhysicalType::other;

    std::vector<Teuchos::RCP<DRT::Condition> > porocondvec;
    masterdis->GetCondition("PoroCoupling",porocondvec);
    for(unsigned int i=0;i<porocondvec.size();++i)
    {
      std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleitergeometry;
      for (eleitergeometry = porocondvec[i]->Geometry().begin();
          eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
      {
        if(faceele->ParentElement()->Id() == eleitergeometry->second->Id())
        {
          if (mastertype_==0)
            dserror("struct and poro master elements on the same processor - no mixed interface supported");
          cele->PhysType() = MORTAR::MortarElement::PhysicalType::poro;
          mastertype_=1;
          break;
        }
      }
    }
    if(cele->PhysType()==MORTAR::MortarElement::PhysicalType::other)
    {
      if (mastertype_==1)
        dserror("struct and poro master elements on the same processor - no mixed interface supported");
      cele->PhysType() = MORTAR::MortarElement::PhysicalType::structure;
      mastertype_=0;
    }

    cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());

    if(isnurbs)
    {
      Teuchos::RCP<DRT::NURBS::NurbsDiscretization> nurbsdis =
          Teuchos::rcp_dynamic_cast<DRT::NURBS::NurbsDiscretization>(masterdis);

      Teuchos::RCP<DRT::NURBS::Knotvector> knots =
          (*nurbsdis).GetKnotVector();
      std::vector<Epetra_SerialDenseVector> parentknots(dim);
      std::vector<Epetra_SerialDenseVector> mortarknots(dim - 1);

      Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
      double normalfac = 0.0;
      bool zero_size = knots->GetBoundaryEleAndParentKnots(
          parentknots,
          mortarknots,
          normalfac,
          faceele->ParentMasterElement()->Id(),
          faceele->FaceMasterNumber());

      // store nurbs specific data to node
      cele->ZeroSized() = zero_size;
      cele->Knots()     = mortarknots;
      cele->NormalFac() = normalfac;
    }

    interface->AddCoElement(cele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;

      Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
                  new CONTACT::CoElement(ele->Id(), ele->Owner(), ele->Shape(),
                      ele->NumNode(), ele->NodeIds(), true,isnurbs));

        Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
        if (faceele == Teuchos::null) dserror("Cast to FaceElement failed!");
        cele->PhysType() = MORTAR::MortarElement::PhysicalType::other;

        std::vector<Teuchos::RCP<DRT::Condition> > porocondvec;
        masterdis->GetCondition("PoroCoupling",porocondvec);

        for(unsigned int i=0;i<porocondvec.size();++i)
        {
          std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleitergeometry;
          for (eleitergeometry = porocondvec[i]->Geometry().begin();
              eleitergeometry != porocondvec[i]->Geometry().end(); ++eleitergeometry)
          {
            if(faceele->ParentElement()->Id() == eleitergeometry->second->Id())
            {
              if (slavetype_==0)
                dserror("struct and poro slave elements on the same processor - no mixed interface supported");
              cele->PhysType() = MORTAR::MortarElement::PhysicalType::poro;
              slavetype_=1;
              break;
            }
          }
        }
        if(cele->PhysType()==MORTAR::MortarElement::PhysicalType::other)
        {
          if (slavetype_==1)
            dserror("struct and poro slave elements on the same processor - no mixed interface supported");
          cele->PhysType() = MORTAR::MortarElement::PhysicalType::structure;
          slavetype_=0;
        }
        cele->SetParentMasterElement(faceele->ParentElement(), faceele->FaceParentNumber());

      if(isnurbs)
      {
        Teuchos::RCP<DRT::NURBS::NurbsDiscretization> nurbsdis =
            Teuchos::rcp_dynamic_cast<DRT::NURBS::NurbsDiscretization>(slavedis);

        Teuchos::RCP<DRT::NURBS::Knotvector> knots =
            (*nurbsdis).GetKnotVector();
        std::vector<Epetra_SerialDenseVector> parentknots(dim);
        std::vector<Epetra_SerialDenseVector> mortarknots(dim - 1);

        Teuchos::RCP<DRT::FaceElement> faceele = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(ele,true);
        double normalfac = 0.0;
        bool zero_size = knots->GetBoundaryEleAndParentKnots(
            parentknots,
            mortarknots,
            normalfac,
            faceele->ParentMasterElement()->Id(),
            faceele->FaceMasterNumber());

        // store nurbs specific data to node
        cele->ZeroSized() = zero_size;
        cele->Knots()     = mortarknots;
        cele->NormalFac() = normalfac;
      }

      interface->AddCoElement(cele);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  read mortar condition                                    Ager 02/16 |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::CreateStrategy(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<DRT::Discretization>   slavedis,
    Teuchos::ParameterList&             input,
    int numcoupleddof)
{
  // poro lagrange strategy:

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  //bools to decide which side is structural and which side is poroelastic to manage all 4 constellations
  // s-s, p-s, s-p, p-p
  bool poromaster = false;
  bool poroslave = false;
  bool structmaster = false;
  bool structslave = false;

  //wait for all processors to determine if they have poro or structural master or slave elements
  comm_->Barrier();
  int slaveTypeList[comm_->NumProc()];
  int masterTypeList[comm_->NumProc()];
  comm_->GatherAll(&slavetype_,&slaveTypeList[0],1);
  comm_->GatherAll(&mastertype_,&masterTypeList[0],1);
  comm_->Barrier();

  for(int i=0; i<comm_->NumProc();++i)
  {
    switch (slaveTypeList[i])
    {
      case -1:
        break;
      case 1:
        if(structslave) dserror("struct and poro slave elements on the same adapter - no mixed interface supported");
        //adjust dserror text, when more than one interface is supported
        poroslave=true;
        break;
      case 0:
        if(poroslave) dserror("struct and poro slave elements on the same adapter - no mixed interface supported");
        structslave=true;
        break;
      default:
        dserror("this cannot happen");
        break;
    }
  }

  for(int i=0; i<comm_->NumProc();++i)
  {
    switch (masterTypeList[i])
    {
      case -1:
        break;
      case 1:
        if(structmaster) dserror("struct and poro master elements on the same adapter - no mixed interface supported");
        //adjust dserror text, when more than one interface is supported
        poromaster=true;
        break;
      case 0:
        if(poromaster) dserror("struct and poro master elements on the same adapter - no mixed interface supported");
        structmaster=true;
        break;
      default:
        dserror("this cannot happen");
        break;
    }
  }

  const Teuchos::ParameterList& stru = DRT::Problem::Instance()->StructuralDynamicParams();
  double theta = stru.sublist("ONESTEPTHETA").get<double>("THETA");
  //what if problem is static ? there should be an error for previous line called in a dyna_statics problem
  //and not a value of 0.5
  //a proper disctinction is necessary if poro meshtying is expanded to other time integration strategies

  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(stru,"DYNAMICTYP")
      == INPAR::STR::dyna_statics)
  {
    theta = 1.0;
  }
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces;
  interfaces.push_back(interface_);
  double alphaf = 1.0 - theta;

  // build the correct data container
  Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr =
      Teuchos::rcp(new CONTACT::AbstractStratDataContainer());
  //create contact poro lagrange strategy for mesh tying
  porolagstrategy_ = Teuchos::rcp(new CONTACT::PoroLagrangeStrategy(
      data_ptr,
      masterdis->DofRowMap(),
      masterdis->NodeRowMap(),
      input,
      interfaces,
      dim,
      comm_,
      alphaf,
      numcoupleddof,
      poroslave,
      poromaster));

  porolagstrategy_->PoroMtInitialize();

  firstinit_ = true;

  return;
}


/*----------------------------------------------------------------------*
 |  complete interface (also print and parallel redist.)      Ager 02/16|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::CompleteInterface(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<CONTACT::CoInterface>& interface)
{
  // finalize the contact interface construction
  int maxdof = masterdis->DofRowMap()->MaxAllGID();
  interface->FillComplete(maxdof);

  //interface->CreateVolumeGhosting(*masterdis);

  // store old row maps (before parallel redistribution)
  slavedofrowmap_  = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowDofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->MasterRowDofs()));
  slavenoderowmap_ = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowNodes()));

  // print parallel distribution
  interface->PrintParallelDistribution(1);

  // TODO: is this possible?
  //**********************************************************************
  // PARALLEL REDISTRIBUTION OF INTERFACE
  //**********************************************************************
//  if (parredist && comm_->NumProc()>1)
//  {
//    // redistribute optimally among all procs
//    interface->Redistribute(1);
//
//    // call fill complete again
//    interface->FillComplete();
//
//    // print parallel distribution again
//    interface->PrintParallelDistribution(1);
//  }

  // store interface
  interface_ = interface;

  // create binary search tree
  interface_->CreateSearchTree();

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate blockmatrices for poro meshtying             ager 10/15    |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::EvaluatePoroMt(
    Teuchos::RCP<Epetra_Vector> fvel,
    Teuchos::RCP<Epetra_Vector> svel,
    Teuchos::RCP<Epetra_Vector> fpres,
    Teuchos::RCP<Epetra_Vector> sdisp,
    const Teuchos::RCP<DRT::Discretization> sdis,
    Teuchos::RCP<LINALG::SparseMatrix>& f,
    Teuchos::RCP<LINALG::SparseMatrix>& k_fs,
    Teuchos::RCP<Epetra_Vector>& frhs,
    ADAPTER::Coupling& coupfs,
    Teuchos::RCP<const Epetra_Map> fdofrowmap)
{
  // safety check
  CheckSetup();

  //write interface values into poro contact interface data containers for integration in contact integrator
  porolagstrategy_->SetState(MORTAR::state_fvelocity,*fvel);
  porolagstrategy_->SetState(MORTAR::state_svelocity,*svel);
  porolagstrategy_->SetState(MORTAR::state_fpressure,*fpres);
  porolagstrategy_->SetState(MORTAR::state_new_displacement,*sdisp);

  //store displacements of parent elements for deformation gradient determinant and its linearization
  porolagstrategy_->SetParentState("displacement",sdisp,sdis);

  interface_->Initialize();
  //in the end of Evaluate coupling condition residuals and linearizations are computed in contact integrator
  interface_->Evaluate();

  porolagstrategy_->PoroMtPrepareFluidCoupling();
  porolagstrategy_->PoroInitialize(coupfs,fdofrowmap,firstinit_);
  if(firstinit_) firstinit_=false;

  //do system matrix manipulations
  porolagstrategy_->EvaluatePoroNoPenContact(k_fs,f,frhs);
  return;
}// ADAPTER::CouplingNonLinMortar::EvaluatePoroMt()


/*----------------------------------------------------------------------*
 |  update poro meshtying quantities                      ager 10/15    |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::UpdatePoroMt()
{
  // safety check
  CheckSetup();

  porolagstrategy_->PoroMtUpdate();
  return;
}//ADAPTER::CouplingNonLinMortar::UpdatePoroMt()

/*----------------------------------------------------------------------*
 |  recover fluid coupling lagrange multiplier            ager 10/15    |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingPoroMortar::RecoverFluidLMPoroMt(Teuchos::RCP<Epetra_Vector> disi,Teuchos::RCP<Epetra_Vector> veli)
{
  // safety check
  CheckSetup();

  porolagstrategy_->RecoverPoroNoPen(disi,veli);
  return;
}//ADAPTER::CouplingNonLinMortar::RecoverFluidLMPoroMt
