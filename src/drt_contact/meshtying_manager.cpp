/*!----------------------------------------------------------------------
\file meshtying_manager.cpp
\brief BACI implementation of main class to control all meshtying

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

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "meshtying_manager.H"
#include "meshtying_lagrange_strategy.H"
#include "meshtying_penalty_strategy.H"
#include "meshtying_defines.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"

#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_knotvector.H"
/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/08|
 *----------------------------------------------------------------------*/
CONTACT::MtManager::MtManager(DRT::Discretization& discret, double alphaf) :
MORTAR::ManagerBase(),
discret_(discret)
{
  // overwrite base class communicator
  comm_ = Teuchos::rcp(Discret().Comm().Clone());

  // welcome message

  // create some local variables (later to be stored in strategy)
  int dim = DRT::Problem::Instance()->NDim();
  if (dim!= 2 && dim!=3) dserror("ERROR: Meshtying problem must be 2D or 3D");
  std::vector<Teuchos::RCP<MORTAR::MortarInterface> > interfaces;
  Teuchos::ParameterList mtparams;

  // read and check meshtying input parameters
  if(Comm().MyPID()==0)
  {
    std::cout << "Checking meshtying input parameters...........";
    fflush(stdout);
  }
  ReadAndCheckInput(mtparams);
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;

  // check for FillComplete of discretization
  if (!Discret().Filled()) dserror("Discretization is not fillcomplete");

  // let's check for meshtying boundary conditions in discret
  // and detect groups of matching conditions
  // for each group, create a contact interface and store it
  if(Comm().MyPID()==0)
  {
    std::cout << "Building meshtying interface(s)...............";
    fflush(stdout);
  }

  std::vector<DRT::Condition*> contactconditions(0);
  Discret().GetCondition("Mortar",contactconditions);

  // there must be more than one meshtying condition
  if ((int)contactconditions.size()<2)
    dserror("Not enough contact conditions in discretization");

  // find all pairs of matching meshtying conditions
  // there is a maximum of (conditions / 2) groups
  std::vector<int> foundgroups(0);
  int numgroupsfound = 0;

  // maximum dof number in discretization
  // later we want to create NEW Lagrange multiplier degrees of
  // freedom, which of course must not overlap with displacement dofs
  int maxdof = Discret().DofRowMap()->MaxAllGID();

  for (int i=0; i<(int)contactconditions.size(); ++i)
  {
    // initialize vector for current group of conditions and temp condition
    std::vector<DRT::Condition*> currentgroup(0);
    DRT::Condition* tempcond = NULL;

    // try to build meshtying group around this condition
    currentgroup.push_back(contactconditions[i]);
    const std::vector<int>* group1v = currentgroup[0]->Get<std::vector<int> >("Interface ID");
    if (!group1v) dserror("Contact Conditions does not have value 'Interface ID'");
    int groupid1 = (*group1v)[0];
    bool foundit = false;

    for (int j=0; j<(int)contactconditions.size(); ++j)
    {
      if (j==i) continue; // do not detect contactconditions[i] again
      tempcond = contactconditions[j];
      const std::vector<int>* group2v = tempcond->Get<std::vector<int> >("Interface ID");
      if (!group2v) dserror("Contact Conditions does not have value 'Interface ID'");
      int groupid2 = (*group2v)[0];
      if (groupid1 != groupid2) continue; // not in the group
      foundit = true; // found a group entry
      currentgroup.push_back(tempcond); // store it in currentgroup
    }

    // now we should have found a group of conds
    if (!foundit) dserror("Cannot find matching contact condition for id %d",groupid1);

    // see whether we found this group before
    bool foundbefore = false;
    for (int j=0; j<numgroupsfound; ++j)
      if (groupid1 == foundgroups[j])
      {
        foundbefore = true;
        break;
      }

    // if we have processed this group before, do nothing
    if (foundbefore) continue;

    // we have not found this group before, process it
    foundgroups.push_back(groupid1);
    ++numgroupsfound;

    // find out which sides are Master and Slave
    bool hasslave = false;
    bool hasmaster = false;
    std::vector<const std::string*> sides((int)currentgroup.size());
    std::vector<bool> isslave((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      sides[j] = currentgroup[j]->Get<std::string>("Side");
      if (*sides[j] == "Slave")
      {
        hasslave = true;
        isslave[j] = true;
      }
      else if (*sides[j] == "Master")
      {
        hasmaster = true;
        isslave[j] = false;
      }
      else
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
    }

    if (!hasslave) dserror("Slave side missing in contact condition group!");
    if (!hasmaster) dserror("Master side missing in contact condition group!");

    // find out which sides are initialized as Active
    std::vector<const std::string*> active((int)currentgroup.size());
    std::vector<bool> isactive((int)currentgroup.size());

    for (int j=0;j<(int)sides.size();++j)
    {
      active[j] = currentgroup[j]->Get<std::string>("Initialization");
      if (*sides[j] == "Slave")
      {
        // slave sides must be initialized as "Active"
        if (*active[j] == "Active")        isactive[j] = true;
        else if (*active[j] == "Inactive") dserror("ERROR: Slave side must be active for meshtying!");
        else                               dserror("ERROR: Unknown contact init qualifier!");
      }
      else if (*sides[j] == "Master")
      {
        // master sides must NOT be initialized as "Active" as this makes no sense
        if (*active[j] == "Active")        dserror("ERROR: Master side cannot be active!");
        else if (*active[j] == "Inactive") isactive[j] = false;
        else                               dserror("ERROR: Unknown contact init qualifier!");
      }
      else
        dserror("ERROR: MtManager: Unknown contact side qualifier!");
    }

    // create an empty meshtying interface and store it in this Manager
    // (for structural meshtying we currently choose redundant master storage)
    INPAR::MORTAR::RedundantStorage redundant = DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(mtparams,"REDUNDANT_STORAGE");
//    if (redundant != INPAR::MORTAR::redundant_master)
//      dserror("ERROR: MtManager: Meshtying requires redundant master storage");
    interfaces.push_back(Teuchos::rcp(new MORTAR::MortarInterface(groupid1,Comm(),dim,mtparams,redundant)));

    // get it again
    Teuchos::RCP<MORTAR::MortarInterface> interface = interfaces[(int)interfaces.size()-1];

    // note that the nodal ids are unique because they come from
    // one global problem discretization conatining all nodes of the
    // contact interface
    // We rely on this fact, therefore it is not possible to
    // do meshtying between two distinct discretizations here

    //-------------------------------------------------- process nodes
    for (int j=0;j<(int)currentgroup.size();++j)
    {
      // get all nodes and add them
      const std::vector<int>* nodeids = currentgroup[j]->Nodes();
      if (!nodeids) dserror("Condition does not have Node Ids");
      for (int k=0; k<(int)(*nodeids).size(); ++k)
      {
        int gid = (*nodeids)[k];
        // do only nodes that I have in my discretization
        if (!Discret().NodeColMap()->MyGID(gid)) continue;
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("Cannot find node with gid %",gid);

        // create MortarNode object
        Teuchos::RCP<MORTAR::MortarNode> mtnode = Teuchos::rcp(new MORTAR::MortarNode(node->Id(),node->X(),
                                                              node->Owner(),
                                                              Discret().NumDof(node),
                                                              Discret().Dof(node),
                                                              isslave[j]));
        //-------------------
        // get nurbs weight!
        if(mtparams.get<bool>("NURBS")==true)
        {
          DRT::NURBS::ControlPoint* cp =
            dynamic_cast<DRT::NURBS::ControlPoint* > (node);

          mtnode->NurbsW() = cp->W();
        }

//#ifdef CONTACTCONSTRAINTXYZ
        // Check, if this node (and, in case, which dofs) are in the contact symmetry condition
        std::vector<DRT::Condition*> contactSymconditions(0);
        Discret().GetCondition("mrtrsym",contactSymconditions);

        for (unsigned j=0; j<contactSymconditions.size(); j++)
          if (contactSymconditions.at(j)->ContainsNode(node->Id()))
          {
            const std::vector<int>*    onoff  = contactSymconditions.at(j)->Get<std::vector<int> >("onoff");
            for (unsigned k=0; k<onoff->size(); k++)
              if (onoff->at(k)==1)
                mtnode->DbcDofs()[k]=true;
          }
//#endif

        // note that we do not have to worry about double entries
        // as the AddNode function can deal with this case!
        interface->AddMortarNode(mtnode);
      }
    }

    //----------------------------------------------- process elements
    int ggsize = 0;
    for (int j=0;j<(int)currentgroup.size();++j)
    {
      // get elements from condition j of current group
      std::map<int,Teuchos::RCP<DRT::Element> >& currele = currentgroup[j]->Geometry();

      // elements in a boundary condition have a unique id
      // but ids are not unique among 2 distinct conditions
      // due to the way elements in conditions are build.
      // We therefore have to give the second, third,... set of elements
      // different ids. ids do not have to be continous, we just add a large
      // enough number ggsize to all elements of cond2, cond3,... so they are
      // different from those in cond1!!!
      // note that elements in ele1/ele2 already are in column (overlapping) map
      int lsize = (int)currele.size();
      int gsize = 0;
      Comm().SumAll(&lsize,&gsize,1);


      std::map<int,Teuchos::RCP<DRT::Element> >::iterator fool;
      for (fool=currele.begin(); fool != currele.end(); ++fool)
      {
        Teuchos::RCP<DRT::Element> ele = fool->second;
        Teuchos::RCP<MORTAR::MortarElement> mtele = Teuchos::rcp(new MORTAR::MortarElement(ele->Id()+ggsize,
                                                                   ele->Owner(),
                                                                   ele->Shape(),
                                                                   ele->NumNode(),
                                                                   ele->NodeIds(),
                                                                   isslave[j],
                                                                   mtparams.get<bool>("NURBS")));
        //------------------------------------------------------------------
        // get knotvector, normal factor and zero-size information for nurbs
        if(mtparams.get<bool>("NURBS")==true)
        {
          DRT::NURBS::NurbsDiscretization* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discret_));

          Teuchos::RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();
          std::vector<Epetra_SerialDenseVector> parentknots(dim);
          std::vector<Epetra_SerialDenseVector> mortarknots(dim-1);

          double normalfac = 0.0;
          bool zero_size=knots->GetBoundaryEleAndParentKnots(parentknots,
                                                             mortarknots,
                                                             normalfac,
                                                             ele->ParentMasterElement()->Id(),
                                                             ele->FaceMasterNumber());

          // store nurbs specific data to node
          mtele->ZeroSized() = zero_size;
          mtele->Knots()     = mortarknots;
          mtele->NormalFac() = normalfac;
        }

        mtele->IsHermite() = DRT::INPUT::IntegralValue<int>(mtparams,"HERMITE_SMOOTHING") ;

        interface->AddMortarElement(mtele);
      } // for (fool=ele1.start(); fool != ele1.end(); ++fool)

      ggsize += gsize; // update global element counter
    }

    //-------------------- finalize the meshtying interface construction
    interface->FillComplete(maxdof);

  } // for (int i=0; i<(int)contactconditions.size(); ++i)
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;

  //**********************************************************************
  // create the solver strategy object
  // and pass all necessary data to it
  if(Comm().MyPID()==0)
  {
    std::cout << "Building meshtying strategy object............";
    fflush(stdout);
  }

  INPAR::CONTACT::SolvingStrategy stype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(mtparams,"STRATEGY");
  if (stype == INPAR::CONTACT::solution_lagmult)
    strategy_ = Teuchos::rcp(new MtLagrangeStrategy(Discret(),mtparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_penalty)
    strategy_ = Teuchos::rcp(new MtPenaltyStrategy(Discret(),mtparams,interfaces,dim,comm_,alphaf,maxdof));
  else if (stype == INPAR::CONTACT::solution_uzawa)
    strategy_ = Teuchos::rcp(new MtPenaltyStrategy(Discret(),mtparams,interfaces,dim,comm_,alphaf,maxdof));
  else
    dserror("Unrecognized strategy");
  if(Comm().MyPID()==0) std::cout << "done!" << std::endl;
  //**********************************************************************

  //**********************************************************************
  // parallel redistribution of all interfaces
  GetStrategy().RedistributeMeshtying();
  //**********************************************************************

  // create binary search tree
  for (int i=0; i<(int)interfaces.size();++i)
    interfaces[i]->CreateSearchTree();

  // print parameter list to screen
  if (Comm().MyPID()==0)
  {
    std::cout << "\ngiven parameters in list '" << GetStrategy().Params().name() << "':\n";
    std::cout << GetStrategy().Params() << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  read and check input parameters (public)                  popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::MtManager::ReadAndCheckInput(Teuchos::ParameterList& mtparams)
{
  // read parameter lists from DRT::Problem
  const Teuchos::ParameterList& mortar    = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList& meshtying = DRT::Problem::Instance()->ContactDynamicParams();
  const Teuchos::ParameterList& wearlist  = DRT::Problem::Instance()->WearParams();

  // read Problem Type and Problem Dimension from DRT::Problem
  const PROBLEM_TYP problemtype = DRT::Problem::Instance()->ProblemType();
  int dim = DRT::Problem::Instance()->NDim();
  std::string distype = DRT::Problem::Instance()->SpatialApproximation();

  // get mortar information
  std::vector<DRT::Condition*> mtcond(0);
  std::vector<DRT::Condition*> ccond(0);

  Discret().GetCondition("Mortar", mtcond);
  Discret().GetCondition("Contact",ccond);

  bool onlymeshtying       = false;
  bool meshtyingandcontact = false;

  // check for case
  if(mtcond.size()!=0 and ccond.size()!=0)
    meshtyingandcontact = true;

  if(mtcond.size()!=0 and ccond.size()==0)
    onlymeshtying = true;


  // *********************************************************************
  // invalid parameter combinations
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") == INPAR::CONTACT::solution_penalty &&
                                                 meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") == INPAR::CONTACT::solution_uzawa &&
                                                 meshtying.get<double>("PENALTYPARAM") <= 0.0)
    dserror("Penalty parameter eps = 0, must be greater than 0");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") == INPAR::CONTACT::solution_uzawa &&
                                                   meshtying.get<int>("UZAWAMAXSTEPS") <  2)
    dserror("Maximum number of Uzawa / Augmentation steps must be at least 2");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") == INPAR::CONTACT::solution_uzawa &&
                                               meshtying.get<double>("UZAWACONSTRTOL") <= 0.0)
    dserror("Constraint tolerance for Uzawa / Augmentation scheme must be greater than 0");

  if (onlymeshtying &&
            DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(meshtying,"FRICTION") != INPAR::CONTACT::friction_none)
    dserror("Friction law supplied for mortar meshtying");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") == INPAR::CONTACT::solution_lagmult &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_standard &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(meshtying,"SYSTEM") == INPAR::CONTACT::system_condensed)
    dserror("Condensation of linear system only possible for dual Lagrange multipliers");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") == INPAR::MORTAR::parredist_dynamic and
      onlymeshtying)
    dserror("ERROR: Dynamic parallel redistribution not possible for meshtying");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none &&
                                                   mortar.get<int>("MIN_ELEPROC") <  0)
    dserror("Minimum number of elements per processor for parallel redistribution must be >= 0");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_DUAL_CONSISTENT")==true &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") != INPAR::CONTACT::solution_lagmult&&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") != INPAR::MORTAR::shape_standard)
    dserror("ERROR: Consistent dual shape functions in boundary elements only for Lagrange multiplier strategy.");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_DUAL_CONSISTENT")==true &&
       DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements &&
       (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_dual ||
        DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_petrovgalerkin   ))

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true &&
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(meshtying,"STRATEGY") != INPAR::CONTACT::solution_lagmult)
    dserror("ERROR: Nodal scaling of Lagrange multipliers only for Lagrange multiplier strategy.");

  if(DRT::INPUT::IntegralValue<int>(mortar,"LM_NODAL_SCALE")==true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar,"INTTYPE") == INPAR::MORTAR::inttype_elements)
   dserror("ERROR: Nodal scaling of Lagrange multipliers not for purely element-based integration.");


  // *********************************************************************
  // not (yet) implemented combinations
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true && dim == 3)
    dserror("ERROR: Crosspoints / edge node modification not yet implemented for 3D");

  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_lin)
    dserror("ERROR: Crosspoints and linear LM interpolation for quadratic FE not yet compatible");

  if (DRT::INPUT::IntegralValue<int>(mortar,"CROSSPOINTS") == true &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(mortar,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    dserror("ERROR: Crosspoints and parallel redistribution not yet compatible");

  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_petrovgalerkin and
      onlymeshtying)
    dserror("Petrov-Galerkin approach makes no sense for meshtying");

  // *********************************************************************
  // 3D quadratic mortar (choice of interpolation and testing fcts.)
  // *********************************************************************
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_pwlin &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_dual)
    dserror("No pwlin approach (for LM) implemented for quadratic meshtying with DUAL shape fct.");

#ifndef MORTARTRAFO
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(mortar,"LM_QUAD") == INPAR::MORTAR::lagmult_lin &&
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(mortar,"LM_SHAPEFCN") == INPAR::MORTAR::shape_dual)
    dserror("Linear approach (for LM) for quadratic meshtying with DUAL shape fct. requires MORTARTRAFO");
#endif // #ifndef MORTARTRAFO

  // *********************************************************************
  // element-based vs. segment-based mortar integration
  // *********************************************************************
  INPAR::MORTAR::IntType inttype = DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(mortar, "INTTYPE");

  if (inttype == INPAR::MORTAR::inttype_segments && mortar.get<int>("NUMGP_PER_DIM") != 0)
    dserror("It is not possible to choose a Gauss rule with NUMGP_PER_DIM for segment-based integration."
            "\nSegment-based integration always uses pre-defined default values (5 GP per segment in 2D "
            "\nand 7 GP per triangular cell in 3D). Ask a 'contact person' if you want to change this.");

  if (inttype == INPAR::MORTAR::inttype_elements && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  if (inttype == INPAR::MORTAR::inttype_elements_BS && mortar.get<int>("NUMGP_PER_DIM") <= 0)
    dserror("Invalid Gauss point number NUMGP_PER_DIM for element-based integration with boundary segmentation."
            "\nPlease note that the value you have to provide only applies to the element-based integration"
            "\ndomain, while pre-defined default values will be used in the segment-based boundary domain.");

  if ((problemtype!=prb_tfsi_aero &&
      (inttype == INPAR::MORTAR::inttype_elements || inttype == INPAR::MORTAR::inttype_elements_BS) &&
      mortar.get<int>("NUMGP_PER_DIM") <= 1))
    dserror(
        "Invalid Gauss point number NUMGP_PER_DIM for element-based integration.");

  // *********************************************************************
  // warnings
  // *********************************************************************
  if (mortar.get<double>("SEARCH_PARAM") == 0.0 && Comm().MyPID()==0)
    std::cout << ("Warning: Meshtying search called without inflation of bounding volumes\n") << std::endl;

  // get parameter lists
  mtparams.setParameters(mortar);
  mtparams.setParameters(meshtying);
  mtparams.setParameters(wearlist);

  // geometrically decoupled elements cannot be given via input file
  mtparams.set<bool>("GEO_DECOUPLED", false);

  // *********************************************************************
  // predefined params for meshtying and contact
  // *********************************************************************
  if(meshtyingandcontact)
  {
    // set options for mortar coupling
    mtparams.set<std::string>("SEARCH_ALGORITHM","Binarytree");
    mtparams.set<double>("SEARCH_PARAM", 0.1);
    mtparams.set<std::string>("SEARCH_USE_AUX_POS", "no");
    mtparams.set<std::string>("PARALLEL_REDIST","static");
    mtparams.set<std::string>("LM_SHAPEFCN","dual");
    mtparams.set<std::string>("REDUNDANT_STORAGE","Master");
    mtparams.set<std::string>("SYSTEM","condensed");
    mtparams.set<bool>("NURBS",false);
  }

  // *********************************************************************
  // smooth interfaces
  // *********************************************************************
  // NURBS PROBLEM?
  if(distype=="Nurbs")
    mtparams.set<bool>("NURBS",true);
  else
    mtparams.set<bool>("NURBS",false);

  if (DRT::INPUT::IntegralValue<int>(mtparams,"HERMITE_SMOOTHING") == true)
    dserror("Hermite smoothing only for mortar contact!");

  mtparams.setName("CONTACT DYNAMIC / MORTAR COUPLING");

  // no parallel redistribution in the serial case
  if (Comm().NumProc()==1) mtparams.set<std::string>("PARALLEL_REDIST","None");

  return true;
}

/*----------------------------------------------------------------------*
 |  write restart information for meshtying (public)          popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart)
{
  // write restart information for meshtying
  Teuchos::RCP<Epetra_Map> problemdofs = GetStrategy().ProblemDofs();
  Teuchos::RCP<Epetra_Vector> lagrmultoldexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*(GetStrategy().LagrMultOld()),*lagrmultoldexp);
  output.WriteVector("mt_lagrmultold",lagrmultoldexp);

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for meshtying (public)           popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::ReadRestart(IO::DiscretizationReader& reader,
                                     Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero)
{
  // this is meshtying, thus we need zeros for restart
  // let strategy object do all the work
  GetStrategy().DoReadRestart(reader, zero);

  return;
}

/*----------------------------------------------------------------------*
 |  write interface tractions for postprocessing (public)     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::MtManager::PostprocessTractions(IO::DiscretizationWriter& output)
{
  // evaluate interface tractions
  Teuchos::RCP<Epetra_Map> problem = GetStrategy().ProblemDofs();
  Teuchos::RCP<Epetra_Vector> traction = Teuchos::rcp(new Epetra_Vector(*(GetStrategy().LagrMultOld())));
  Teuchos::RCP<Epetra_Vector> tractionRescaled = Teuchos::rcp(new Epetra_Vector(*(GetStrategy().LagrMultOldRescaled())));
  Teuchos::RCP<Epetra_Vector> tractionexp = Teuchos::rcp(new Epetra_Vector(*problem));
  LINALG::Export(*tractionRescaled, *tractionexp);

  // evaluate slave and master forces
  Teuchos::RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(GetStrategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(GetStrategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcslaveexp = Teuchos::rcp(new Epetra_Vector(*problem));
  Teuchos::RCP<Epetra_Vector> fcmasterexp = Teuchos::rcp(new Epetra_Vector(*problem));
  GetStrategy().DMatrix()->Multiply(true, *traction, *fcslave);
  GetStrategy().MMatrix()->Multiply(true, *traction, *fcmaster);
  LINALG::Export(*fcslave, *fcslaveexp);
  LINALG::Export(*fcmaster, *fcmasterexp);

  // write to output
  output.WriteVector("interfacetraction",tractionexp);
  output.WriteVector("slaveforces",fcslaveexp);
  output.WriteVector("masterforces",fcmasterexp);

  return;
}

