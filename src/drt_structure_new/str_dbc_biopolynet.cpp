/*-----------------------------------------------------------*/
/*!
\file str_dbc_biopolynet.cpp

\brief dbc management for biopolymer network simulations

\maintainer Jonas Eichinger

\date Jun 14, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_dbc_biopolynet.H"

#include "str_timint_base.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beam3/beam3r.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_locsys.H"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>
#include <Teuchos_ParameterList.hpp>
#include "str_model_evaluator_browniandyn.H"


/*----------------------------------------------------------------------*
 | constructor                                  (public) eichinger 06/16|
 *----------------------------------------------------------------------*/
STR::DbcBioPolyNet::DbcBioPolyNet()
 :smdyn_ptr_(Teuchos::null),
  basisnodes_(0),
  basiselements_(0),
  useinitdbcset_(false)
{
  // empty constructor
} // DbcBioPolyNet()

/*----------------------------------------------------------------------*
 | initialize dirichlet statmech object         (public) eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::Init(const Teuchos::RCP<DRT::Discretization>&     discret_ptr_,
                            const Teuchos::RCP<Epetra_Vector>&           freact_ptr,
                            const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  // call Init from base class
  STR::Dbc::Init(discret_ptr_,freact_ptr,timint_ptr);

  smdyn_ptr_ = timint_ptr_->GetDataSDyn().GetMeshFreeDataPtr();

  // initial clearing of dbc management vectors
  dbcnodesets_.clear();
  // number of nodes which are part of the discretization already before starting adding crosslinkers
  basisnodes_    = discret_ptr_->NumGlobalNodes();
  // number of nodes which are part of the discretization already before starting adding crosslinkers
  basiselements_ = discret_ptr_->NumGlobalElements();

  // BC sanity checks
  DbcSanityCheck();

  // print dbc type
  BioPolyNetPrintBCType();

  isinit_ = true;
} // Init()

/*----------------------------------------------------------------------*
 |  print Boundary condition type to screen      (public) mueller 03/12 |
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::BioPolyNetPrintBCType()
{
  if(!DiscretPtr()->Comm().MyPID())
  {
    INPAR::STATMECH::DBCType dbctype = smdyn_ptr_->DbcType();
    switch(dbctype)
    {
      // standard DBC application
      case INPAR::STATMECH::dbctype_std:
        std::cout<<"- Conventional Input file based application of DBCs"<<std::endl;
      break;
      // shear with a fixed Dirichlet node set
      case INPAR::STATMECH::dbctype_shearfixed:
        std::cout<<"- DBCs for rheological measurements applied: fixed node set"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_shearfixeddel:
        std::cout<<"- DBCs for rheological measurements applied: fixed node set"<<std::endl;
      break;
      // shear with an updated Dirichlet node set (only DOF in direction of oscillation is subject to BC, others free)
      case INPAR::STATMECH::dbctype_sheartrans:
        std::cout<<"- DBCs for rheological measurements applied: transient node set"<<std::endl;
      break;
      // pin down and release individual nodes
      case INPAR::STATMECH::dbctype_pinnodes:
        std::cout<<"- Special DBCs pinning down selected nodes"<<std::endl;
      break;
      // apply affine shear deformation
      case INPAR::STATMECH::dbctype_affineshear:
        std::cout<<"- DBCs for affine shear deformation"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_affinesheardel:
        std::cout<<"- DBCs for affine shear deformation"<<std::endl;
      break;
      case INPAR::STATMECH::dbctype_movablesupport1d:
        std::cout<<"- DBCs 1D movable support"<<std::endl;
      break;
      // no DBCs at all
      case INPAR::STATMECH::dbctype_none:
        std::cout<<"- No application of DBCs (i.e. also no DBCs by Input file)"<<std::endl;
      break;
      // default: everything involving periodic boundary conditions
      default:
        dserror("Check your DBC type! %i", dbctype);
      break;
    }

  }
  return;
} //BioPolyNetPrintBCType()

/*----------------------------------------------------------------------*
 | apply Dirichlet Boundary Conditions            (public) mueller 03/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::ApplyDirichletBC(const double&               time,
                                        Teuchos::RCP<Epetra_Vector> dis,
                                        Teuchos::RCP<Epetra_Vector> vel,
                                        Teuchos::RCP<Epetra_Vector> acc,
                                        bool                        recreatemap
  )
{
  CheckInitSetup();
  // We have to rotate forward ...
  // ---------------------------------------------------------------------------
  if (!dis.is_null())
    RotateGlobalToLocal(dis,true);
  if (!vel.is_null())
    RotateGlobalToLocal(vel);
  if (!acc.is_null())
    RotateGlobalToLocal(acc);

  // Apply DBCs
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList p;
  p.set("total time",time);

  // Only new DBC application method implemented using the map extractor! Old version using toggle vectors discontinued
  if(dbcmap_ptr_ == Teuchos::null)
    dserror("Only new DBC application method implemented using the map extractor! Old version using toggle vectors discontinued!");
  recreatemap = true;

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  DiscretPtr()->ClearState();
  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"DBCTYPE");

  // standard DBC application
  if(dbctype == INPAR::STATMECH::dbctype_std)
    DiscretPtr()->EvaluateDirichlet(p, dis, vel, acc, Teuchos::null, dbcmap_ptr_);
  // default: everything involving periodic boundary conditions
  else
  EvaluateDirichletPeriodic(p, dis, vel, dbcmap_ptr_);

  DiscretPtr()->ClearState();

  // We have to rotate back into global Cartesian frame
  // ---------------------------------------------------------------------------
  if (dis != Teuchos::null)
    RotateLocalToGlobal(dis, true);
  if (vel != Teuchos::null)
    RotateLocalToGlobal(vel);
  if (acc != Teuchos::null)
    RotateLocalToGlobal(acc);

  return;
} // ApplyDirichletBC()

/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (private)        mueller 02/10|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::EvaluateDirichletPeriodic(Teuchos::ParameterList&     params,
                                                 Teuchos::RCP<Epetra_Vector> dis,
                                                 Teuchos::RCP<Epetra_Vector> vel,
                                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
#ifdef MEASURETIME
  const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

  // some preliminary checks
  if (!(DiscretPtr()->Filled())) dserror("FillComplete() was not called");
  if (!(DiscretPtr()->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");

  // retrieve DBC type
  INPAR::STATMECH::DBCType dbctype = smdyn_ptr_->DbcType();
  // vector holding Dirichlet increments
  Teuchos::RCP<Epetra_Vector> deltadbc = Teuchos::rcp(new Epetra_Vector(*(DiscretPtr()->DofRowMap()),true));
  // vector of DOF-IDs which are Dirichlet BCs
  Teuchos::RCP<std::set<int> > dbcgids = Teuchos::rcp(new std::set<int>());

  // take action according to DBC type in the input file
  switch(dbctype)
  {
    // shear with a fixed Dirichlet node set
    case INPAR::STATMECH::dbctype_shearfixed:
    // shear with a fixed Dirichlet node set and setting axial stiffness of disrupted elements close to zero
    case INPAR::STATMECH::dbctype_shearfixeddel:
    {
      DBCOscillatoryMotion(params,dis,vel,deltadbc);
      useinitdbcset_ = true;
    }
    break;
    // shear with an updated Dirichlet node set (only DOF in direction of oscillation is subject to BC, others free)
    case INPAR::STATMECH::dbctype_sheartrans:
    {
      DBCOscillatoryMotion(params,dis,vel,deltadbc);
    }
    break;
    // pin down and release individual nodes
    case INPAR::STATMECH::dbctype_pinnodes:
    {
      // apply DBCs from input file first
      DBCGetPredefinedConditions(params,dis,vel,dbcgids);
      // then, determine nodes that get inhibited by Crosslinkers
      DBCPinNodes();
    }
    break;
    // apply affine shear displacement to all nodes
    case INPAR::STATMECH::dbctype_affineshear:
    // apply affine shear displacement to all nodes
    case INPAR::STATMECH::dbctype_affinesheardel:
    {
      DBCAffineShear(params,dis,vel,deltadbc);
      // used here to store the node set that remains under DBCs after the initial affine displacement
      useinitdbcset_ = true;
    }
    break;
    // apply movable support to all upper face nodes of cut elements
    case INPAR::STATMECH::dbctype_movablesupport1d:
    {
      DBCMovableSupport1D(params,dis);
    }
    break;
    default:
      return;
  }
//------------------------------------set Dirichlet values
  DBCSetValues(dis, deltadbc, dbcgids);

//----create DBC and free map and build their common extractor
  DBCCreateMap(dbcgids, dbcmapextractor);

#ifdef MEASURETIME
  const double t_end = Teuchos::Time::wallTime();
  if(!DiscretPtr()->Comm().MyPID())
    std::cout<<"DBC Evaluation time: "<<t_end-t_start<<std::endl;
#endif // #ifdef MEASURETIME

  return;
} //EvaluateDirichletPeriodic()


/*----------------------------------------------------------------------*
 | Gather information on where DBCs are to be applied in the case of    |
 | viscoelastic measurements                   (private)  mueller 05/12 |
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCOscillatoryMotion(Teuchos::ParameterList&     params,
                                            Teuchos::RCP<Epetra_Vector> dis,
                                            Teuchos::RCP<Epetra_Vector> vel,
                                            Teuchos::RCP<Epetra_Vector> deltadbc)
{
  INPAR::STATMECH::DBCType dbctype = smdyn_ptr_->DbcType();
  // get the absolute time, time step size, time curve number and oscillation direction
  const double time = params.get<double>("total time", 0.0); // target time (i.e. timen_)
  double dt = (*timint_ptr_->GetDataGlobalState().GetDeltaTime())[0];
  int curvenumber = smdyn_ptr_->CurveNumber();
  int oscdir = smdyn_ptr_->DbcDispDir();
  double amp = smdyn_ptr_->ShearAmplitude();

  // time curve increment
  double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                       DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);

  std::vector<double> curvefac(2,0.0);
  if(vel!=Teuchos::null)
    curvefac = DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1);

//------------------------------------------------------gather dbc node set(s)
  // after the the very first entry into the following part, reenter only if the Dirichlet Node Set is dynamically updated
  if(!useinitdbcset_)
  {
    // note: we need column map displacements as we might need access to ghost nodes
    Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*DiscretPtr()->DofColMap(), true));
    LINALG::Export(*dis, *discol);

    // clear node sets, as they will be filled for the first time or updated with new node sets
    dbcnodesets_.clear();
    // vectors to manipulate DBC
    std::vector<int> oscillnodes;
    std::vector<int> fixednodes;
    oscillnodes.clear();
    fixednodes.clear();
    // toggle vector
    std::vector<bool> nodedbcstatus(DiscretPtr()->NumMyColNodes(), false);

    for(int lid=0; lid<DiscretPtr()->NumMyRowElements(); lid++)
    {
      // An element used to browse through local Row Elements
      DRT::Element* element = DiscretPtr()->lRowElement(lid);
      // skip element if it is a crosslinker element or in addition, in case of the Bead Spring model, Torsion3 elements
      if(element->Id() <= basisnodes_)
      {
        // number of translational DOFs
        int numdof = 3;
        // positions of nodes of an element with n nodes
        LINALG::SerialDenseMatrix coord(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode(), true);
        // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
        LINALG::SerialDenseMatrix cut(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode()-1,true);

//-------------------------------- obtain nodal coordinates of the current element
        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);

//-----------------------detect broken/fixed/free elements and fill position vector
        // determine existence and location of broken element
        if(CheckForBrokenElement(coord,cut))
        {
          // reduce the axial stiffness of the element drastically, close to zero in order to take this element out
          // only for Beam3r case
          if(element->ElementType()==DRT::ELEMENTS::Beam3rType::Instance() && dbctype == INPAR::STATMECH::dbctype_shearfixeddel)
          {
            for(int n=0; n<cut.N(); n++)
            {
              if(cut(2,n)>0.0)
              {
                //std::cout<<"Element "<<element->Id()<<" now has a reduced cross section"<<std::endl;
                dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSec(1.0e-9);
                dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSecShear(1.1e-9);
                break;
              }
            }
          }
          // loop over number of cuts (columns)
          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn     = DiscretPtr()->NodeColMap()->LID(element->Nodes()[n]->Id());
            int nodelidnp    = DiscretPtr()->NodeColMap()->LID(element->Nodes()[n+1]->Id());
            int noderowlidn  = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n]->Id());
            int noderowlidnp = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n+1]->Id());

            // only consider this node pair if it is not already subject to Dirichlet BCs
            if(nodedbcstatus.at(nodelidn)==false && nodedbcstatus.at(nodelidnp)==false)
            {
              // case 1: broken element (in z-dir); node_n+1 oscillates, node_n is fixed in dir. of oscillation
              if(cut(2,n)==1.0)
              {
                // update dbc status for these nodes (dofs)
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                // add GID of fixed node to fixed-nodes-vector (to be added to condition later)
                if(noderowlidn>-1)
                  fixednodes.push_back(element->Nodes()[n]->Id());
                // add GID of oscillating node to osc.-nodes-vector
                if(noderowlidnp>-1)
                {
                  oscillnodes.push_back(element->Nodes()[n+1]->Id());

                  // incremental Dirichlet displacement for an oscillating node (all DOFs except oscdir = 0.0)
                  (*deltadbc)[doflids.at(numdof*(n+1)+oscdir)] = (coord(2,n+1)/(*smdyn_ptr_->PeriodLength())[2])*amp*tcincrement;
                  if(vel!=Teuchos::null)
                    (*vel)[doflids.at(numdof*(n+1)+oscdir)] = (coord(2,n+1)/(*smdyn_ptr_->PeriodLength())[2])*amp*curvefac[1];
                }
              }// end of case 1
              if(cut(2,n)==2.0)// case 2: broken element (in z-dir); node_n oscillates, node_n+1 is fixed in dir. of oscillation
              {
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                if(noderowlidnp>-1)
                  fixednodes.push_back(element->Nodes()[n+1]->Id());
                // oscillating node
                if(noderowlidn>-1)
                {
                  oscillnodes.push_back(element->Nodes()[n]->Id());
                  (*deltadbc)[doflids.at(numdof*n+oscdir)] = (coord(2,n)/(*smdyn_ptr_->PeriodLength())[2])*amp*tcincrement;
                  if(vel!=Teuchos::null)
                    (*vel)[doflids.at(numdof*n+oscdir)] = (coord(2,n)/(*smdyn_ptr_->PeriodLength())[2])*amp*curvefac[1];
                }
              } // end of case 2
            }
          }
        }
      }
    }
//-----------------------------------------update node sets
    dbcnodesets_.push_back(oscillnodes);
    dbcnodesets_.push_back(fixednodes);
//---------check/set force sensors anew for each time step
//    UpdateForceSensors(dbcnodesets_[0], oscdir);
  }
  else // only ever use the fixed set of nodes after the first time dirichlet values were imposed
  {
    // dofs of fixednodes_ remain untouched since fixednodes_ dofs are 0.0
    for(int i=0; i<(int)dbcnodesets_[0].size(); i++)
    {
      int nodelid = DiscretPtr()->NodeRowMap()->LID(dbcnodesets_[0][i]);
      DRT::Node* oscnode = DiscretPtr()->lRowNode(nodelid);
      std::vector<int> dofnode = DiscretPtr()->Dof(oscnode);
      // oscillating node
      double tcincrement = 0.0;
      if(curvenumber>-1)
        tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                      DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
      (*deltadbc)[DiscretPtr()->DofRowMap()->LID(dofnode[oscdir])] = amp*tcincrement;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | pins down a particular given node by Dirichlet BCs                   |
 |                                           (private)     mueller 05/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCPinNodes()
{
  dserror("PinNodes() not yet implemented in new structure time integration.");
//  // gather node sets
//  // end nodes of horizontal filaments
//  std::vector<int> hfilfixednodes;
//  // end nodes of vertical filaments
//  std::vector<int> vfilfixednodes;
//  // regular nodes
//  std::vector<int> fixednodes;
//
//  // clear sets initially
//  dbcnodesets_.clear();
//  hfilfixednodes.clear();
//  vfilfixednodes.clear();
//  fixednodes.clear();
//
//  int maxnodegid = DiscretPtr()->NumGlobalNodes()-1;
//
//  for(int i=0; i<bspotstatus_->MyLength(); i++)
//  {
//    if((*bspotstatus_)[i]>-0.1)
//    {
//      int nodegid = bspotcolmap_->GID(i);
//      if(nodegid==0)
//      {
//        hfilfixednodes.push_back(nodegid);
//        continue;
//      }
//      // the very last node
//      if(nodegid==maxnodegid)
//      {
//        if((*filamentnumber_)[bspotcolmap_->LID(maxnodegid)]==0)
//          hfilfixednodes.push_back(nodegid);
//        else
//          vfilfixednodes.push_back(nodegid);
//        continue;
//      }
//      // in between end nodes and regular nodes
//      if(nodegid<maxnodegid && nodegid>0)
//      {
//        // look for end nodes
//        if((*filamentnumber_)[i-1]==(*filamentnumber_)[i] && (*filamentnumber_)[i+1]!=(*filamentnumber_)[i])
//        {
//          if((*filamentnumber_)[i]==0)
//            hfilfixednodes.push_back(bspotcolmap_->GID(i));
//          else
//            vfilfixednodes.push_back(bspotcolmap_->GID(i));
//        }
//        else
//          fixednodes.push_back(bspotcolmap_->GID(i));
//      }
//    }
//  }
//
//  // add to node sets
//  dbcnodesets_.push_back(hfilfixednodes);
//  dbcnodesets_.push_back(vfilfixednodes);
//  dbcnodesets_.push_back(fixednodes);

  return;
}

/*----------------------------------------------------------------------*
 | Apply affine shear displacement              (private)  mueller 05/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCAffineShear(Teuchos::ParameterList&     params,
                                      Teuchos::RCP<Epetra_Vector> dis,
                                      Teuchos::RCP<Epetra_Vector> vel,
                                      Teuchos::RCP<Epetra_Vector> deltadbc)
{
  /*====the following is mostly taken from DBCOscillatoryMotion()====
   * we apply an affine deformation to the whole network, i.e. all nodes
   * are displaced according to their vertical position in order to
   * shear the network by an angle gamma*/

  const double time = params.get<double>("total time", 0.0);
  double dt = params.get<double>("delta time", 0.01);
  int curvenumber = smdyn_ptr_->CurveNumber();
  int displacementdir = smdyn_ptr_->DbcDispDir();
  double amp = smdyn_ptr_->ShearAmplitude();
  double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                       DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);

  std::vector<double> curvefac(2,0.0);
  if(vel!=Teuchos::null)
    curvefac = DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1);

  // We need column map displacements as we might need access to ghost nodes
  Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*DiscretPtr()->DofColMap(), true));
  LINALG::Export(*dis, *discol);

  // reenter once the third time threshold is hit in order to release the
  if(!useinitdbcset_)
  {
    dbcnodesets_.clear();
    std::vector<int> sensornodes;
    std::vector<int> fixednodes;
    std::vector<int> affineshearnodes;
    sensornodes.clear();
    fixednodes.clear();
    affineshearnodes.clear();
    std::vector<bool> nodedbcstatus(DiscretPtr()->NumMyRowNodes(), false);

    // displaced and fixed nodes
    for(int lid=0; lid<DiscretPtr()->NumMyRowElements(); lid++)
    {
      DRT::Element* element = DiscretPtr()->lRowElement(lid);
      // skip element if it is a crosslinker element or in addition, in case of the Bead Spring model, Torsion3 elements
      if(element->Id() <= basisnodes_)
      {
        // number of translational DOFs
        int numdof = 3;
        LINALG::SerialDenseMatrix coord(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode(), true);
        // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
        LINALG::SerialDenseMatrix cut(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode()-1,true);

        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);

        // shifted elements
        if(CheckForBrokenElement(coord,cut))
        {
          // reduce the axial stiffness of the element drastically, close to zero in order to take this element out
          // only for Beam3r case
          for(int n=0; n<cut.N(); n++)
            if(element->ElementType()==DRT::ELEMENTS::Beam3rType::Instance() && cut(2,n)>0.0)
            {
              dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSec(1.0e-12);
              dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSecShear(1.1e-12);
              break;
            }

          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n]->Id());
            //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": nodelidn = "<<nodelidn<<"/"<<DiscretPtr()->NodeRowMap()->GID(nodelidn)<<std::endl;
            if(nodelidn>-1)
            {
              //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": nodelidn cut = "<<(int)cut(2,n)<<std::endl;
              if(nodedbcstatus.at(nodelidn)==false)
              {
                switch((int)(round(cut(2,n))))
                {
                  // case 1: broken element (in z-dir); node_n+1 displaced, node_n is fixed in dir. of displacement
                  case 1:
                  {
                    //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": case 1 n"<<std::endl;
                    nodedbcstatus.at(nodelidn) = true;
                    fixednodes.push_back(element->Nodes()[n]->Id());
                  }
                  break;
                  // case 2: broken element (in z-dir); node_n displaced, node_n+1 is fixed in dir. of displacement
                  case 2:
                  {
                    //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": case 2 n"<<std::endl;
                    nodedbcstatus.at(nodelidn) = true;
                    sensornodes.push_back(element->Nodes()[n]->Id());
                    (*deltadbc)[doflids.at(numdof*n+displacementdir)] =  (coord(2,n)/(*smdyn_ptr_->PeriodLength())[2])*amp*tcincrement;
                    if(vel!=Teuchos::null)
                      (*vel)[doflids.at(numdof*n+displacementdir)] = (coord(2,n)/(*smdyn_ptr_->PeriodLength())[2])*amp*curvefac[1];
                  }
                  break;
                }
              }
            }
            int nodelidnp = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n+1]->Id());
            //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": nodelidnp = "<<nodelidnp<<"/"<<DiscretPtr()->NodeRowMap()->GID(nodelidnp)<<std::endl;
            if(nodelidnp>-1)
            {
              //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": nodelidnp cut = "<<(int)cut(2,n)<<std::endl;
              if(nodedbcstatus.at(nodelidnp)==false)
              {
                switch((int)(round(cut(2,n))))
                {
                  // case 1: broken element (in z-dir); node_n+1 displaced, node_n is fixed in dir. of displacement
                  case 1:
                  {
                    //std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": case 1 np"<<std::endl;
                    nodedbcstatus.at(nodelidnp) = true;
                    sensornodes.push_back(element->Nodes()[n+1]->Id());
                    (*deltadbc)[doflids.at(numdof*(n+1)+displacementdir)] = (coord(2,n+1)/(*smdyn_ptr_->PeriodLength())[2])*amp*tcincrement;
                    if(vel!=Teuchos::null)
                      (*vel)[doflids.at(numdof*(n+1)+displacementdir)] = (coord(2,n+1)/(*smdyn_ptr_->PeriodLength())[2])*amp*curvefac[1];
                  }
                  break;
                  // case 2: broken element (in z-dir); node_n displaced, node_n+1 is fixed in dir. of displacement
                  case 2:
                  {
                   // std::cout<<"Proc "<<DiscretPtr()->Comm().MyPID()<<": case 2 np"<<std::endl;
                    nodedbcstatus.at(nodelidnp) = true;
                    fixednodes.push_back(element->Nodes()[n+1]->Id());
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
    // affine shear nodes (in analogy to prior DBC nodes
    for(int lid=0; lid<DiscretPtr()->NumMyRowElements(); lid++)
    {
      DRT::Element* element = DiscretPtr()->lRowElement(lid);
      if(element->Id() <= basisnodes_)
      {
        int numdof = 3;
        LINALG::SerialDenseMatrix coord(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode(), true);
        LINALG::SerialDenseMatrix cut(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode()-1,true);
        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);
        for(int n=0; n<cut.N(); n++)
        {
          for(int m=n; m<n+2; m++)
          {
            int nodelidm = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[m]->Id());
            if(nodelidm>-1)
            {
              if(nodedbcstatus.at(nodelidm)==false)
              {
                nodedbcstatus.at(nodelidm) = true;
                affineshearnodes.push_back(DiscretPtr()->NodeRowMap()->GID(nodelidm));
                (*deltadbc)[doflids.at(numdof*m+displacementdir)]     =  (coord(2,m)/(*smdyn_ptr_->PeriodLength())[2])*amp*tcincrement;
              }
            }
          }
        }
      }
    }
    dbcnodesets_.push_back(sensornodes);
    dbcnodesets_.push_back(fixednodes);

//    UpdateForceSensors(dbcnodesets_[0], displacementdir);
  }
  else // only ever use the fixed set of nodes after the first time dirichlet values were imposed
  {
    // top plate nodes
    for(int i=0; i<(int)dbcnodesets_[0].size(); i++)
    {
      int nodelid = DiscretPtr()->NodeRowMap()->LID(dbcnodesets_[0][i]);
      DRT::Node* displacednode = DiscretPtr()->lRowNode(nodelid);
      std::vector<int> dofnode = DiscretPtr()->Dof(displacednode);

      double znode = displacednode->X()[2] + (*discol)[DiscretPtr()->DofColMap()->LID(dofnode[2])];
      double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                           DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
      (*deltadbc)[DiscretPtr()->DofRowMap()->LID(dofnode[displacementdir])] = znode/(*smdyn_ptr_->PeriodLength())[2]*amp*tcincrement;
    }
    // all other network nodes that undergo affine deformation (hard coded for now)
    if(dbcnodesets_.size()==3)
    {
      for(int i=0; i<(int)dbcnodesets_[2].size(); i++)
      {
        int nodelid = DiscretPtr()->NodeRowMap()->LID(dbcnodesets_[2][i]);
        DRT::Node* affinenode = DiscretPtr()->lRowNode(nodelid);
        std::vector<int> dofnode = DiscretPtr()->Dof(affinenode);

        double znode = affinenode->X()[2] + (*dis)[DiscretPtr()->DofRowMap()->LID(dofnode[2])];
        double tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
                             DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
        (*deltadbc)[DiscretPtr()->DofRowMap()->LID(dofnode[displacementdir])] = znode/(*smdyn_ptr_->PeriodLength())[2]*amp*tcincrement;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Dirichlet condition constraining all but one translational           |
 | direction                                    (private)  mueller 09/13|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCMovableSupport1D(Teuchos::ParameterList&     params,
                                           Teuchos::RCP<Epetra_Vector> dis)
{
  // note in analogy to DBCOscillatoryMotion()
//  int freedir = smdyn_->get<int>("DBCDISPDIR",-1)-1;
  if(!useinitdbcset_)
  {
    Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*DiscretPtr()->DofColMap(), true));
    LINALG::Export(*dis, *discol);

    dbcnodesets_.clear();
    std::vector<int> movablenodes;
    std::vector<int> fixednodes;
    movablenodes.clear();
    fixednodes.clear();
    std::vector<bool> nodedbcstatus(DiscretPtr()->NumMyColNodes(), false);

    for(int lid=0; lid<DiscretPtr()->NumMyRowElements(); lid++)
    {
      DRT::Element* element = DiscretPtr()->lRowElement(lid);
      if(element->Id() <= basisnodes_)
      {
        int numdof = 3;
        LINALG::SerialDenseMatrix coord(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode(), true);
        LINALG::SerialDenseMatrix cut(numdof,(int)DiscretPtr()->lRowElement(lid)->NumNode()-1,true);
        std::vector<int> doflids;
        GetElementNodeCoords(element, discol, coord, &doflids);
        if(CheckForBrokenElement(coord,cut))
        {
//          if(element->ElementType()==DRT::ELEMENTS::Beam3rType::Instance())
//          {
//            for(int n=0; n<cut.N(); n++)
//            {
//              if(cut(2,n)>0.0)
//              {
//                dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSec(1.0e-12);
//                dynamic_cast<DRT::ELEMENTS::Beam3r*>(element)->SetCrossSecShear(1.1e-12);
//                break;
//              }
//            }
//          }
          for(int n=0; n<cut.N(); n++)
          {
            int nodelidn = DiscretPtr()->NodeColMap()->LID(element->Nodes()[n]->Id());
            int nodelidnp = DiscretPtr()->NodeColMap()->LID(element->Nodes()[n+1]->Id());
            int noderowlidn = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n]->Id());
            int noderowlidnp = DiscretPtr()->NodeRowMap()->LID(element->Nodes()[n+1]->Id());

            if(nodedbcstatus.at(nodelidn)==false && nodedbcstatus.at(nodelidnp)==false)
            {
              if(cut(2,n)==1.0)
              {
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                if(noderowlidn>-1)
                  fixednodes.push_back(element->Nodes()[n]->Id());
                if(noderowlidnp>-1)
                  movablenodes.push_back(element->Nodes()[n+1]->Id());
              }
              if(cut(2,n)==2.0)
              {
                nodedbcstatus.at(nodelidn) = true;
                nodedbcstatus.at(nodelidnp) = true;
                if(noderowlidnp>-1)
                  fixednodes.push_back(element->Nodes()[n+1]->Id());
                if(noderowlidn>-1)
                  movablenodes.push_back(element->Nodes()[n]->Id());
              }
            }
          }
        }
      }
    }
    dbcnodesets_.push_back(movablenodes);
    dbcnodesets_.push_back(fixednodes);
//    UpdateForceSensors(dbcnodesets_[0], freedir);
    useinitdbcset_ = true;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (private), mwgee 01/07, mueller 05/12 |
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DoDirichletConditionPredefined(DRT::Condition&              cond,
                                                      DRT::Discretization&         dis,
                                                      const bool                   usetime,
                                                      const double                 time,
                                                      Teuchos::RCP<Epetra_Vector>  systemvector,
                                                      Teuchos::RCP<Epetra_Vector>  systemvectord,
                                                      Teuchos::RCP<Epetra_Vector>  systemvectordd,
                                                      Teuchos::RCP<Epetra_Vector>  toggle,
                                                      Teuchos::RCP<std::set<int> > dbcgids)
{
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
  const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
  const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");
  const std::vector<double>* val    = cond.Get<std::vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // loop nodes to identify and evaluate load curves and spatial distributions
  // of Dirichlet boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    int nlid = dis.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = dis.lRowNode( nlid );

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = dis.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    for (unsigned j=0; j<total_numdf; ++j)
    {
      int onesetj = j % numdf;
      if ((*onoff)[onesetj]==0)
      {
        const int lid = (*systemvectoraux).Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (toggle!=Teuchos::null)
          (*toggle)[lid] = 0.0;
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
        continue;
      }
      const int gid = dofs[j];
      std::vector<double> value(deg+1,(*val)[onesetj]);

      // factor given by time curve
      std::vector<double> curvefac(deg+1, 1.0);
      int curvenum = -1;
      if (curve) curvenum = (*curve)[onesetj];
      if (curvenum>=0 && usetime)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
      else
        for (unsigned i=1; i<(deg+1); ++i) curvefac[i] = 0.0;

      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct)
      {
        funct_num = (*funct)[onesetj];
        if (funct_num>0)
          functfac =
            DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(onesetj,
                                                                  actnode->X(),
                                                                  time,
                                                                  &dis);
      }

      // apply factors to Dirichlet value
      for (unsigned i=0; i<deg+1; ++i)
      {
        value[i] *= functfac * curvefac[i];
      }

      // assign value
      const int lid = (*systemvectoraux).Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      if (systemvector != Teuchos::null)
        (*systemvector)[lid] = value[0];
      if (systemvectord != Teuchos::null)
        (*systemvectord)[lid] = value[1];
      if (systemvectordd != Teuchos::null)
        (*systemvectordd)[lid] = value[2];
      // set toggle vector
      if (toggle != Teuchos::null)
        (*toggle)[lid] = 1.0;
      // amend vector of DOF-IDs which are Dirichlet BCs
      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}

/*----------------------------------------------------------------------*
 |  fill system vector and toggle vector (private)         mueller  3/10|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DoDirichletConditionPeriodic(std::vector<int>*            nodeids,
                                                    std::vector<int>*            onoff,
                                                    Teuchos::RCP<Epetra_Vector>  dis,
                                                    Teuchos::RCP<Epetra_Vector>  deltadbc,
                                                    Teuchos::RCP<std::set<int> > dbcgids)
/*
 * This basically does the same thing as DoDirichletCondition() (to be found in drt_discret_evaluate.cpp),
 * but with the slight difference of taking current displacements into account.
 * Time curve values aren't added to the reference position(s) of the discretization as usual,
 * but to the latest known 0-position(s). These positions are calculated using the deltadbc
 * vector holding the latest incremental Dirichlet displacement. It is added to the displacement
 * at the end of the preceding time step.
 */
{
  /*/ "condition output"
  std::cout<<"Node Ids: ";
  for(int i=0; i<(int)nodeids->size(); i++)
    std::cout<<nodeids->at(i)<<" ";
  std::cout<<"onoff: ";
  for(int i=0; i<(int)DiscretPtr()->Dof(0,DiscretPtr()->gNode(nodeids->at(0))).size(); i++)
    std::cout<<onoff->at(i)<<" ";
  std::cout<<std::endl;*/

  // some checks for errors
  if (!nodeids) dserror("No Node IDs were handed over!");

  // get the condition properties
  const int nnode = nodeids->size();

  // loop over all nodes in condition
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!DiscretPtr()->NodeRowMap()->MyGID(nodeids->at(i))) continue;
    DRT::Node* actnode = DiscretPtr()->gNode(nodeids->at(i));
    if (!actnode) dserror("Cannot find global node %d",nodeids->at(i));

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = DiscretPtr()->Dof(0,actnode);
    const unsigned numdf = dofs.size();
    // loop over DOFs
    for (unsigned j=0; j<numdf; ++j)
    {
      // get the GID and the LID of the currently handled DOF
      const int gid = dofs[j];
      const int lid = (*dis).Map().LID(gid);

//---------------------------------------------Dirichlet Value Assignment
      if (onoff->at(j)==1)
      {
        // assign value
        if (lid<0) dserror("Global id %d not on this proc in system vector", dofs[j]);
        if (dis != Teuchos::null)
          (*dis)[lid] += (*deltadbc)[lid];
        // amend vector of DOF-IDs which are Dirichlet BCs
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }
      else // if DOF in question is not subject to DBCs (anymore)
      {
        if (lid<0)
          dserror("Global id %d not on this proc in system vector",dofs[j]);
        // get rid of entry in DBC map - if it exists
        if (dbcgids != Teuchos::null)
          (*dbcgids).erase(dofs[j]);
      }
    }  // loop over nodal DOFs
  }  // loop over nodes
  return;
}

/*----------------------------------------------------------------------*
 | Get DBCs defined in Input file and ad the DBCs DOFs to dbcgids       |            |
 |                                           (private)     mueller 05/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCGetPredefinedConditions(Teuchos::ParameterList&      params,
                                                  Teuchos::RCP<Epetra_Vector>  dis,
                                                  Teuchos::RCP<Epetra_Vector>  vel,
                                                  Teuchos::RCP<std::set<int> > dbcgids)
{
  //get absolute time
  double time = params.get<double>("total time", 0.0);
  bool usetime = true;
  if (time<0.0) usetime = false;

  // get the Dirichlet conditions given in the input file
  std::vector<DRT::Condition*> dirichlet;
  DiscretPtr()->GetCondition("Dirichlet",dirichlet);

  for (int i=0; i<(int)dirichlet.size(); i++)
  {
    if (dirichlet[i]->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletConditionPredefined(*dirichlet[i],*DiscretPtr(),usetime,time,dis,vel,Teuchos::null,Teuchos::null,dbcgids);
  }
  for (int i=0; i<(int)dirichlet.size(); i++)
  {
    if (dirichlet[i]->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletConditionPredefined(*dirichlet[i],*DiscretPtr(),usetime,time,dis,vel,Teuchos::null,Teuchos::null,dbcgids);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Set Dirichlet Values                 (private)         mueller  5/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCSetValues(Teuchos::RCP<Epetra_Vector>  dis,
                                    Teuchos::RCP<Epetra_Vector>  deltadbc,
                                    Teuchos::RCP<std::set<int> > dbcgids)
{
  // number of DOFs per node and DOF toggle vectors
  int numdof = (int)DiscretPtr()->Dof(DiscretPtr()->gNode(DiscretPtr()->NodeRowMap()->GID(0))).size();
  // DOF toggle vector
  std::vector<std::vector<int> > onoff((int)dbcnodesets_.size(),std::vector<int>(numdof,0));

  // manipulation of onoff vectors according to the different DBC types
  INPAR::STATMECH::DBCType dbctype = smdyn_ptr_->DbcType();
  switch(dbctype)
  {
    case INPAR::STATMECH::dbctype_shearfixed:
      // inhibit translational degrees of freedom
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
          onoff[i].at(j) = 1;
    break;
    case INPAR::STATMECH::dbctype_shearfixeddel:
      // inhibit translational degrees of freedom
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
          onoff[i].at(j) = 1;
    break;
    case INPAR::STATMECH::dbctype_sheartrans:
    {
      int dbcdispdir = smdyn_ptr_->DbcDispDir();
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        onoff[i].at(dbcdispdir) = 1;
    }
    break;
    case INPAR::STATMECH::dbctype_affineshear:
      // inhibit translational degrees of freedom
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
          onoff[i].at(j) = 1;
    break;
    case INPAR::STATMECH::dbctype_affinesheardel:
      // inhibit translational degrees of freedom
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
          onoff[i].at(j) = 1;
    break;
    case INPAR::STATMECH::dbctype_pinnodes:
      // DOF toggle vector for end nodes of the horizontal filament
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
          onoff[i].at(j) = 1;
    break;
    case INPAR::STATMECH::dbctype_movablesupport1d:
    {
      int freedir = smdyn_ptr_->DbcDispDir();
      // inhibit translational degrees of freedom
      for(int i=0; i<(int)dbcnodesets_.size(); i++)
        for(int j=0; j<3; j++)
        {
          if(i==0) //movable set
          {
            if(j!=freedir)
              onoff[i].at(j) = 1;
          }
          else // fixed set
            onoff[i].at(j) = 1;
        }
    }
    break;
    default:
      dserror("DBC values cannot be imposed for this DBC type: %d", dbctype);
    break;
  }

  // Apply DBC values
  for(int i=0; i<(int)dbcnodesets_.size(); i++)
    if(!dbcnodesets_[i].empty())
      DoDirichletConditionPeriodic(&dbcnodesets_[i], &onoff[i], dis, deltadbc, dbcgids);

  return;
}

/*----------------------------------------------------------------------*
 | Create Dirichlet DOF maps                    (private)  mueller 05/12|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DBCCreateMap(Teuchos::RCP<std::set<int> >       dbcgids,
                                    Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor)
{
  if (dbcmapextractor != Teuchos::null)
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids->size() > 0)
    {
      dbcgidsv.reserve(dbcgids->size());
      dbcgidsv.assign(dbcgids->begin(),dbcgids->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements, DiscretPtr()->DofRowMap()->IndexBase(), DiscretPtr()->DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *dbcmapextractor = LINALG::MapExtractor(*(DiscretPtr()->DofRowMap()), dbcmap);
  }
  return;
}

/*----------------------------------------------------------------------*
 | get a matrix with node coordinates and their DOF LIDs                |
 |                                                 (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::GetElementNodeCoords(DRT::Element* element,
                                            Teuchos::RCP<Epetra_Vector> discol,
                                            LINALG::SerialDenseMatrix& coord,
                                            std::vector<int>* lids)
{
  // clear LID vector just in case it was handed over non-empty
  lids->clear();
  for (int j=0; j<element->NumNode(); j++)
  {
    int nodegid = element->NodeIds()[j];
    DRT::Node* node = DiscretPtr()->lColNode(DiscretPtr()->NodeColMap()->LID(nodegid));
    for (int k=0; k<3; k++)
    {
      // obtain k-th spatial component of the reference position of the j-th node
      double referenceposition = node->X()[k];
      // get the GIDs of the node's DOFs
      std::vector<int> dofnode = DiscretPtr()->Dof(node);
      // store the displacement of the k-th spatial component
      double displacement = (*discol)[DiscretPtr()->DofColMap()->LID(dofnode[k])];
      // write updated components into coord (only translational
      coord(k, j) = referenceposition + displacement;
      // store current lid(s) (3 translational DOFs per node)
      if (lids != NULL)
        lids->push_back(DiscretPtr()->DofRowMap()->LID(dofnode[k]));
    }
  }
  return;
} // BioPolyNetManager::GetElementNodeCoords

/*----------------------------------------------------------------------*
 | Some simple sanity checks                   (private)  mueller 01/13 |
 *----------------------------------------------------------------------*/
void STR::DbcBioPolyNet::DbcSanityCheck()
{
  INPAR::STATMECH::DBCType dbctype = smdyn_ptr_->DbcType();
  switch(dbctype)
  {
    case INPAR::STATMECH::dbctype_std:
      break;
    case INPAR::STATMECH::dbctype_pinnodes:
      break;
    case INPAR::STATMECH::dbctype_none:
      break;
    // default: everything involving periodic boundary conditions and rheology
    default:
    {
      // "-1" because of counting convention
      int curvenumber = smdyn_ptr_->CurveNumber();
      int displacementdir = smdyn_ptr_->DbcDispDir();

      if(smdyn_ptr_->PeriodLength()->at(0)<=0.0)
        dserror("PERIODLENGTH = %f! This method works only in conjunction with periodic boundary conditions!", smdyn_ptr_->PeriodLength()->at(0));
      if(smdyn_ptr_->PeriodLength()->at(0)!=smdyn_ptr_->PeriodLength()->at(1) || smdyn_ptr_->PeriodLength()->at(0)!=smdyn_ptr_->PeriodLength()->at(2) || smdyn_ptr_->PeriodLength()->at(1)!=smdyn_ptr_->PeriodLength()->at(2))
        dserror("Check PERIODLENGTH in your input file! This method only works for a cubic boundary volume");
      if(displacementdir>2 || displacementdir<0 || (curvenumber<0 && dbctype != INPAR::STATMECH::dbctype_movablesupport1d))
        dserror("In case of imposed DBC values, please define the BioPolyNet parameters DBCDISPDIR={1,2,3} and/or CURVENUMBER correctly");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | check for broken element                        (public) mueller 3/10|
 *----------------------------------------------------------------------*/
bool STR::DbcBioPolyNet::CheckForBrokenElement(LINALG::SerialDenseMatrix& coord,
                                             LINALG::SerialDenseMatrix& cut)
{
  // empty cut just in case it was handed over non-empty
  cut.Zero();
  bool broken = false;
  int ndim = coord.M();
  // flag "broken" signals a broken element, flag "cut" hints at location of nodes 0 and 1 (of the two-noded beam)
  for (int dof = 0; dof < ndim; dof++)
    //loop through columns of "cut"
    for (int n = 0; n < cut.N(); n++)
    {
      // broken element with node_n close to "0.0"-boundary
      if (fabs(coord(dof, n + 1) - smdyn_ptr_->PeriodLength()->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        // set value for the spatial component in question at n-th cut
        cut(dof, n) = 1.0;
      }
      else if (fabs(coord(dof, n + 1) + smdyn_ptr_->PeriodLength()->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        cut(dof, n) = 2.0;
      }
    }
  return broken;
}// CheckForBrokenElement()
