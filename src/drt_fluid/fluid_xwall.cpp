/*----------------------------------------------------------------------*/
/*!
\file fluid_xwall.cpp
\brief XWall

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_xwall.H"
#include "fluidimplicitintegration.H"
#include "../drt_fluid_ele/fluid_ele_xwall.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/drt_periodicbc.H"

#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 04/14 |
 *----------------------------------------------------------------------*/
FLD::XWall::XWall(
    Teuchos::RCP<DRT::Discretization>      dis,
    int                           nsd,
    Teuchos::RCP<Teuchos::ParameterList>& params):
    discret_(dis),
    params_(params)
{

  // get the processor ID from the communicator
  myrank_  = discret_->Comm().MyPID();

  if(myrank_==0)
  {
    std::cout << "\nWall modeling with a Spalding's law enrichment" << std::endl;
  }

  //some exclusions and safety checks:
  if (nsd!=3)
    dserror("Only 3D problems considered in xwall modelling!");
  if(DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") != INPAR::FLUID::timeint_afgenalpha&&
      DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") != INPAR::FLUID::timeint_npgenalpha)
    dserror("Use Af-Genalpha for time integration in combination with xwall wall modeling. There would be additional updates necessary otherwise");

  if(params_->get<std::string>("predictor","steady_state_predictor")!= "steady_state")
    dserror("The meshtying framework does only support a steady-state predictor");

  std::string tauwtype = params_->sublist("WALL MODEL").get<std::string>("Tauw_Type","constant");

  if(tauwtype=="constant")
    tauwtype_=INPAR::FLUID::constant;
  else if(tauwtype=="mean_between_steps")
    tauwtype_=INPAR::FLUID::mean_between_steps;
  else if(tauwtype=="mean_iter")
    tauwtype_=INPAR::FLUID::mean_iter;
  else if(tauwtype=="between_steps")
    tauwtype_=INPAR::FLUID::between_steps;
  else if(tauwtype=="fix_point_iter_with_step_control")
    tauwtype_=INPAR::FLUID::fix_point_iter_with_step_control;
  else if(tauwtype=="fully_linearized")
    tauwtype_=INPAR::FLUID::fully_linearized;
  else
    dserror("unknown Tauw_Type");

  std::string tauwcalctype = params_->sublist("WALL MODEL").get<std::string>("Tauw_Calc_Type","residual");

  if(tauwcalctype=="residual")
    tauwcalctype_=INPAR::FLUID::residual;
  else if(tauwcalctype=="spalding")
    tauwcalctype_=INPAR::FLUID::spalding;
  else if(tauwcalctype=="gradient")
    tauwcalctype_=INPAR::FLUID::gradient;
  else if(tauwcalctype=="gradient_to_residual")
    tauwcalctype_=INPAR::FLUID::gradient_to_residual;
  else
    dserror("unknown Tauw_Calc_Type");

  constant_tauw_=params_->sublist("WALL MODEL").get<double>("C_Tauw") ;

  min_tauw_=params_->sublist("WALL MODEL").get<double>("Min_Tauw") ;

  fac_=params_->sublist("WALL MODEL").get<double>("Inc_Tauw") ;

  //get gauss points
  gp_norm_=params_->sublist("WALL MODEL").get<int>("GP_Wall_Normal") ;
  gp_par_=params_->sublist("WALL MODEL").get<int>("GP_Wall_Parallel") ;

  // compute initial pressure
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id==-1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  dens_ = actmat->density_;
  visc_ = actmat->viscosity_/dens_;//here I want to have the kinematic viscosity


  std::string projectiontype = params_->sublist("WALL MODEL").get<std::string>("Projection","No");

  if(projectiontype=="onlyl2projection")
  {
    projconstr_=false;
    proj_=true;
  }
  else if(projectiontype=="l2projectionwithcontinuityconstraint")
  {
    projconstr_=true;
    proj_=false;
  }
  else if(projectiontype=="No")
  {
    projconstr_=false;
    proj_=false;
  }
  else
    dserror("unknown projection type");

  penalty_param_ =params_->sublist("WALL MODEL").get<double>("Penalty_Param") ;

  std::string blendingtype = params_->sublist("WALL MODEL").get<std::string>("Blending_Type","none");

  if(blendingtype=="none")
    blendingtype_=INPAR::FLUID::none;
  else if(blendingtype=="ramp_function")
    blendingtype_=INPAR::FLUID::ramp_function;
  else if(blendingtype=="tauw_transformation")
    blendingtype_=INPAR::FLUID::tauw_transformation;
  else
    dserror("unknown Blending_Type");
  if(myrank_==0)
    std::cout << "\nConsider changing the element formulation such that y+ lives on the standard finite element space! this would improve tauw_transformation blending" << std::endl;

  inctauwnorm_=0.0;

  smooth_res_aggregation_ = DRT::INPUT::IntegralValue<int>(params_->sublist("WALL MODEL"),"SMOOTH_TAUW");

  switch_step_ = params_->sublist("WALL MODEL").get<int>("Switch_Step");
  if(tauwcalctype_ == INPAR::FLUID::gradient_to_residual && switch_step_ < 2)
    dserror("provide reasonable Switch_Step if you want to use gradient_to_residual");

  if(smooth_res_aggregation_ && tauwcalctype_ == INPAR::FLUID::gradient)
    dserror("smoothing of tauw works only for residual-based tauw, as the residual is smoothed");

  SepEnr_=Teuchos::null;


  //output:
  if(myrank_==0)
  {
    std::cout << "\nXWall settings: " << std::endl;
    std::cout << "Tau_w is updated with:        " << tauwtype << std::endl;
    std::cout << "Tau_w is calculated with:     " << tauwcalctype << std::endl;
    std::cout << "Switching from grad to res:   " << switch_step_ << std::endl;
    std::cout << "Constant tau_w:               " << constant_tauw_ << std::endl;
    std::cout << "Minimum tau_w (clipping):     " << min_tauw_ << std::endl;
    std::cout << "Increment of tau_w:           " << fac_ << std::endl;
    std::cout << "Gauss rule:                   normal:  " << gp_norm_ << "  parallel:  " << gp_par_ << "  overall:  " << gp_norm_*gp_par_*gp_par_ << std::endl;
    std::cout << "Enriched DOFs l2-projected:   "<< projectiontype << std::endl;
    std::cout << "Penalty parameter:            " << penalty_param_ << std::endl;
    std::cout << "Blending method:              " << blendingtype << std::endl;
    std::cout << "Smooth tau_w:                 " << smooth_res_aggregation_ << std::endl;
    std::cout << "Solver for tau_w smoothing:   " << params_->sublist("WALL MODEL").get<int>("ML_SOLVER") << std::endl;
    std::cout << "Solver for projection:        " << params_->sublist("WALL MODEL").get<int>("PROJECTION_SOLVER") << std::endl;
    std::cout << std::endl;
    std::cout << "WARNING: ramp functions are used to treat fluid Neumann inflow conditions" << std::endl;
    std::cout << "WARNING: ramp functions are used to treat fluid Mortar coupling conditions" << std::endl;
    std::cout << "WARNING: face element with enrichment not implemented" << std::endl;
  }

  Setup();
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetXWallParams(Teuchos::ParameterList& eleparams)
{

  //makes sure to use the xwall template
//  eleparams.set<std::string>("impltype","xwall");

  //params required for the shape functions
  eleparams.set("walldist",wdist_);
  eleparams.set("tauw",tauw_);
  eleparams.set("inctauw",inctauw_);
  eleparams.set("xwalltoggle",xwalltoggle_);

  eleparams.set("mk",mkstate_);

  eleparams.set<int>("gpnorm",gp_norm_);
  eleparams.set<int>("gppar",gp_par_);

  return;
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 08/14 |
 |  Used for the xwdiscret_, which is redistributed                     |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetXWallParamsXWDis(Teuchos::ParameterList& eleparams)
{

  //params required for the shape functions
  eleparams.set("walldist",wdistxwdis_);
  eleparams.set("tauw",tauwxwdis_);
  eleparams.set("inctauw",inctauwxwdis_);
  eleparams.set("xwalltoggle",xwalltogglexwdis_);

  eleparams.set("mk",mkxwstate_);

  eleparams.set<int>("gpnorm",gp_norm_);
  eleparams.set<int>("gppar",gp_par_);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup XWall                                                bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::Setup()
{
  if(myrank_==0)
    std::cout << "Setup: "<< std::endl;

  InitXWallMaps();

  InitWallDist();

  SetupXWallDis();

  tauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()),true));
  inctauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()),true));

  //initialize for first call or for constant tauw setting
  tauw_->PutScalar(constant_tauw_);

  wdistxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()),true));
  LINALG::Export(*walldist_,*wdistxwdis_);

  tauwxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()),true));
  inctauwxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()),true));

  //initialize for first call or for constant tauw setting
  tauwxwdis_->PutScalar(constant_tauw_);

  InitToggleVector();

  SetupL2Projection();

  {
    //the value for linear elements is 1/3
    //initialize just in case
    mkxwstate_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->ElementColMap()),true));
    mkxwstate_->PutScalar(0.33333333333);
    mkstate_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementColMap()),true));
    mkstate_->PutScalar(0.33333333333);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Setup basic xwall map and dirichlet map                    bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitXWallMaps()
{
  if(myrank_==0)
    std::cout << "- build xwall maps...                                   " ;

  //build row vector of all xwall nodes
  {
    std::vector<int> rowvec;          // node row map
    for (int i=0;i<discret_->NodeRowMap()->NumMyElements();++i)
    {
      int xwallgid = discret_->NodeRowMap()->GID(i);
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if(!xwallnode) dserror("Cannot find node");

      bool enriched=false;

      //check if one of the surrounding elements is xwall element
      DRT::Element** surrele = xwallnode->Elements();
      for (int k=0;k<xwallnode->NumElement();++k)
      {
        DRT::ELEMENTS::FluidXWall* xwallele=dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

        if(xwallele)
          enriched=true;
      }

      if(enriched)
        rowvec.push_back(xwallgid);
    }

    xwallrownodemap_ = Teuchos::rcp(new Epetra_Map(-1,(int)rowvec.size(),&rowvec[0],0,discret_->Comm()));
  }

  //get Dirichlet conditions
  std::vector<DRT::Condition*> dircond;
  discret_->GetCondition("FluidStressCalc",dircond);

  if (not dircond.empty())
  {
    std::vector<int> testcollect;
    int count = 0;
    for (unsigned numcond=0;numcond<dircond.size();++numcond)
    {
      const std::vector<int>* test= dircond[numcond]->Nodes();
      int j=0;
      for( std::vector<int>::const_iterator i = (*test).begin(); i != (*test).end(); ++i)
      {
        ++count;
        testcollect.push_back((*test)[j]);
        ++j;
      }
    }

    int gcount;
    (discret_->Comm()).SumAll(&count,&gcount,1);
    dircolnodemap_ = Teuchos::rcp(new Epetra_Map(gcount,count,&testcollect[0],0,discret_->Comm()));
  } // end loop this conditions
  else
    dserror("You need DESIGN FLUID STRESS CALC SURF CONDITIONS for xwall");


  //map is of course not unique as it is a column map
//  if(dircolnodemap_->UniqueGIDs())
//    dserror("Map resulting from DESIGN FLUID STRESS CALC SURF CONDITIONS not unique, probably node specified on two conditions?");

  if(myrank_==0)
    std::cout << xwallrownodemap_->NumGlobalElements()<<  " XWall nodes initialized!" << std::endl;
  if(xwallrownodemap_->NumGlobalElements()==0)
    dserror("No XWall elements found");
  return;
}


/*----------------------------------------------------------------------*
 |  Calculate wall distance and coupling matrix of for tauw    bk 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitWallDist()
{
  if(myrank_==0)
    std::cout << "- calculate wall distance...                            " ;

  tauwcouplingmattrans_ = Teuchos::rcp(new LINALG::SparseMatrix(*xwallrownodemap_,2,false,false));
  wdist_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()),true));
  walldist_ = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));

  //build a new discretization which lies on all procs
  Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(discret_->Comm().Clone());

  //this is very expensive in terms of memory
  //we will delete it as soon as we are ready here
  Teuchos::RCP<DRT::Discretization> commondis=Teuchos::rcp(new DRT::Discretization((std::string)"Commondis",newcomm));

  // loop over all column nodes of underlying problem discret and add
  for (int i=0;i<(discret_->NodeColMap())->NumMyElements();++i)
  {
    DRT::Node* node = discret_->lColNode(i);
    if (!node) dserror("Cannot find node with lid %",i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    commondis->AddNode(newnode);
  }
  // loop over all column elements of underlying problem discret and add
  for (int i=0;i<(discret_->ElementColMap())->NumMyElements();++i)
  {
    DRT::Element* node = discret_->lColElement(i);
    if (!node) dserror("Cannot find ele with lid %",i);
    Teuchos::RCP<DRT::Element> newnode = Teuchos::rcp(node->Clone());
    commondis->AddElement(newnode);
  }

  Teuchos::RCP<Epetra_Map> testrednodecolmap = LINALG::AllreduceEMap(*(discret_->NodeRowMap()));
  commondis->ExportColumnNodes(*testrednodecolmap);
  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (commondis->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new discretization
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentDofSet(discret_,parallel));

  commondis->ReplaceDofSet(newdofset);
  commondis->FillComplete(true,false,false);

  newdofset=Teuchos::null;

  //also build a fully overlapping map of enriched nodes here:
  std::vector<int> colvec;          // node col map
  for (int i=0;i<(commondis->NodeColMap())->NumMyElements();++i)
  {
    int gid =commondis->NodeColMap()->GID(i);
    DRT::Node* xwallnode = commondis->lColNode(i);
    if (!xwallnode) dserror("Cannot find node with lid %",i);
    int enriched=0;

    DRT::Element** surrele = xwallnode->Elements();
    for (int k=0;k<xwallnode->NumElement();++k)
    {
      DRT::ELEMENTS::FluidXWall* xwallele=dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

      if(xwallele)
        enriched=1;
    }
    int genriched=0;
    (commondis->Comm()).SumAll(&enriched,&genriched,1);
    if(genriched>0)
      colvec.push_back(gid);
  }
  int count=(int)colvec.size();

  xwallcolnodemap_ = Teuchos::rcp(new Epetra_Map(count,count,&colvec[0],0,discret_->Comm()));

  for (int j=0; j<xwallcolnodemap_->NumMyElements();++j)
  {
    int xwallgid = xwallcolnodemap_->GID(j);

    DRT::Node* xwallnode = commondis->gNode(xwallgid);
    if(!xwallnode) dserror("Cannot find node");

    double mydist=1.0E10;
    double gdist=1.0E10;
    int mygid=0;


    for (int i=0; i<dircolnodemap_->NumMyElements(); ++i)
    {
      int gid = dircolnodemap_->GID(i);

      if (discret_->NodeRowMap()->MyGID(gid))
      {
        DRT::Node* node = discret_->gNode(gid);

        if (!node) dserror("ERROR: Cannot find wall node with gid %",gid);

        double newdist=sqrt(((xwallnode->X())[0]-(node->X())[0])*((xwallnode->X())[0]-(node->X())[0])+((xwallnode->X())[1]-(node->X())[1])*((xwallnode->X())[1]-(node->X())[1])+((xwallnode->X())[2]-(node->X())[2])*((xwallnode->X())[2]-(node->X())[2]));
        if(newdist < mydist)
        {
          mydist=newdist;
          mygid=gid;
        }
      }
    }

    discret_->Comm().MinAll(&mydist,&gdist,1);

    //now write this value in the node based vector
    if (xwallrownodemap_->MyGID(xwallgid))
    {
      int err =     walldist_->ReplaceGlobalValues(1,&gdist,&xwallgid);
      if(err>0)
        dserror("global row not on proc");
      else if(err<0)
        dserror("wrong vector index");
    }

    //this is the processor that knows the respective node
    if(mydist==gdist)
    {
      tauwcouplingmattrans_->Assemble(1.0,mygid,xwallgid);
    }
  }

  LINALG::Export(*walldist_,*wdist_);
  tauwcouplingmattrans_->Complete();

  double mean=0.0;
  walldist_->MeanValue(&mean);

  if(myrank_==0)
    std::cout << "the mean distance from the wall of all XWall nodes is: "<< mean << "... ";

  if(myrank_==0)
    std::cout << "done!  " << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Build a node-based toggle vector (on/off=1.0/0.0, 0.7)     bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitToggleVector()
{
  if(myrank_==0)
    std::cout << "- build enriched/blending toggle vector for elements... ";

  xwalltoggle_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()),true));
  xwalltogglexwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()),true));
  xtoggleloc_ = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));
  int count =0;
  for (int j=0; j<xwallrownodemap_->NumMyElements();++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (discret_->NodeRowMap()->MyGID(xwallgid)) //just in case
    {
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if(!xwallnode) dserror("Cannot find node");

      bool fullyenriched=true;

      // Neumann inflow
      std::vector<DRT::Condition*> inflcond;
      xwallnode->GetCondition("FluidNeumannInflow",inflcond);
      if(not inflcond.empty())
        fullyenriched=false;

      // Mortar interface
      std::vector<DRT::Condition*> mortarcond;
      xwallnode->GetCondition("Mortar",mortarcond);
      if(not mortarcond.empty())
        fullyenriched=false;

      //get all surrounding elements
      DRT::Element** surrele = xwallnode->Elements();
      for (int k=0;k<xwallnode->NumElement();++k)
      {
        DRT::ELEMENTS::FluidXWall* xwallele=dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

        if(!xwallele)
          fullyenriched=false;
      }

      if(fullyenriched==true)
      {
        int err =     xtoggleloc_->ReplaceMyValue(j,0,1.0);
        if (err!=0)
          dserror("something went wrong");
      }
      else
      {
        if(blendingtype_!=INPAR::FLUID::ramp_function)
        {
          int err =     xtoggleloc_->ReplaceMyValue(j,0,0.7);
          if (err!=0)
            dserror("something went wrong");
        }
        count++;
      }
    }
  }

  LINALG::Export(*xtoggleloc_,*xwalltoggle_);
  LINALG::Export(*xtoggleloc_,*xwalltogglexwdis_);

  int gcount;
  (discret_->Comm()).SumAll(&count,&gcount,1);
  if(myrank_==0)
    std::cout << gcount << " blending nodes identified... ";

  if(myrank_==0)
    std::cout << "done!  " << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Build a discretization only including xwall nodes/elements bk 07/14 |
 |  and redistribute                                                    |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetupXWallDis()
{
  //build a new discretization
  Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(discret_->Comm().Clone());

  xwdiscret_=Teuchos::rcp(new DRT::Discretization((std::string)"xwalldis",newcomm));

  // loop over all xwall row nodes and add
  for (int i=0;i<(discret_->NodeColMap())->NumMyElements();++i)
  {
    int gid = (discret_->NodeColMap())->GID(i);

    if (xwallcolnodemap_->MyGID(gid))
    {
      DRT::Node* node = discret_->lColNode(i);
      if (!node) dserror("Cannot find node with lid %",i);
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      xwdiscret_->AddNode(newnode);
    }
  }
  // loop over all column elements of underlying problem discret and add
  for (int i=0;i<discret_->NumMyColElements();++i)
  {
    DRT::Element* ele = discret_->lColElement(i);
    if (!ele) dserror("Cannot find ele with lid %",i);

    DRT::ELEMENTS::FluidXWall* xwallele=dynamic_cast<DRT::ELEMENTS::FluidXWall*>(ele);

    if(!xwallele)
      continue;

    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    xwdiscret_->AddElement(newele);
  }

  // make all conditions known to the child discretization
  // i.e. periodic boundary conditions, dirichlet conditions, ...
  {
    // get all conditions types prescribed in the input file
    std::vector<std::string> allcond;
    discret_->GetConditionNames(allcond);
    //loop all conditions types
    for (unsigned numcond=0;numcond<allcond.size();++numcond)
    {
      // get condition
      std::vector<DRT::Condition*> actcond;
      discret_->GetCondition(allcond[numcond],actcond);
      // loop all condition of the current type
      for (unsigned numactcond=0;numactcond<actcond.size();++numactcond)
      {
        // finally set condition
        xwdiscret_->SetCondition(allcond[numcond],Teuchos::rcp(new DRT::Condition(*actcond[numactcond])));
      }
    }
  }


  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (xwdiscret_->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new discretization
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentDofSet(discret_,parallel));

  xwdiscret_->ReplaceDofSet(newdofset);
  xwdiscret_->FillComplete(true,true,true);

  //redistribute and treat periodic bc if parallel
  if(parallel)
  {
  //redistribute
#if defined(PARALLEL) && defined(PARMETIS)

    Teuchos::RCP<Epetra_Map> elemap=Teuchos::rcp( new Epetra_Map(*xwdiscret_->ElementRowMap()));
    Teuchos::RCP<Epetra_Map> rownodes;
    Teuchos::RCP<Epetra_Map> colnodes;
    Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(discret_->Comm().Clone());
    DRT::UTILS::PartUsingParMetis(xwdiscret_, elemap, rownodes, colnodes, comm,true);
    // rebuild of the system with new maps
    xwdiscret_->Redistribute(*rownodes,*colnodes,false,false);

#else
#if defined(PARALLEL)
    dserror("require PARMETIS not METIS");
#endif
#endif

    PeriodicBoundaryConditions pbc(xwdiscret_,false);
    pbc.UpdateDofsForPeriodicBoundaryConditions();
    xwdiscret_->ReplaceDofSet(newdofset);
    xwdiscret_->FillComplete(true,true,true);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Setup matrix, vectors and solver for projection            bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetupL2Projection()
{

  //create matrix for projection
  if(proj_)
  {
    //build the dof maps of p and u_y

    std::vector<int> enrdf;          // enriched dofs

    for (int i=0; i<xwdiscret_->DofRowMap()->NumMyElements(); ++i)
    {
      int gdfid=xwdiscret_->DofRowMap()->GID(i);
      if(gdfid%8>3&&gdfid%8<7)
        enrdf.push_back(gdfid);
    }

    enrdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)enrdf.size(),&enrdf[0],0,xwdiscret_->Comm()));

    massmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*enrdofrowmap_,108,false,true));

    incveln_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
    incvelnp_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
    incaccn_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));

    stateveln_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
    statevelnp_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
    stateaccn_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));

    mergedmap_=Teuchos::null;
    lagrdofrowmap_=Teuchos::null;
  }
  else if (projconstr_)
  {
    //build the dof maps of p and u_y

    std::vector<int> enrdf;          // enriched dofs
    std::vector<int> lagrdf;          // lagrange dofs

    for (int i=0; i<xwdiscret_->DofRowMap()->NumMyElements(); ++i)
    {
      int gdfid=xwdiscret_->DofRowMap()->GID(i);
      if(gdfid%8>3&&gdfid%8<7)
        enrdf.push_back(gdfid);
      if(gdfid%8>6)
        lagrdf.push_back(gdfid);
    }

    enrdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)enrdf.size(),&enrdf[0],0,xwdiscret_->Comm()));
    lagrdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)lagrdf.size(),&lagrdf[0],0,xwdiscret_->Comm()));
    mergedmap_=LINALG::MergeMap(*enrdofrowmap_,*lagrdofrowmap_,false);

    massmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*mergedmap_,108,false,true));

    incveln_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
    incvelnp_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
    incaccn_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));

    stateveln_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
    statevelnp_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
    stateaccn_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
  }
  else
  {
    massmatrix_=Teuchos::null;
    incveln_=Teuchos::null;
    incvelnp_=Teuchos::null;
    incaccn_=Teuchos::null;
    stateveln_=Teuchos::null;
    statevelnp_=Teuchos::null;
    stateaccn_=Teuchos::null;
    enrdofrowmap_ =Teuchos::null;
    lagrdofrowmap_ =Teuchos::null;
  }

    //setup solver
  if(proj_ || projconstr_)
  {
    const int solvernumber = params_->sublist("WALL MODEL").get<int>("PROJECTION_SOLVER");
    if(solvernumber<0)
      dserror("provide a solver number for l2-projection");
    // get solver parameter list of linear solver
    const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(solvernumber);
    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

    solver_ =
              Teuchos::rcp(new LINALG::Solver(solverparams,
              xwdiscret_->Comm(),
              DRT::Problem::Instance()->ErrorFile()->Handle()));

    if(solvertype != INPAR::SOLVER::umfpack)
    {
      if(solvertype != INPAR::SOLVER::belos&&myrank_==0)
          std::cout << "\nUse Belos as solver because it can handle several right hand sides at once!\n"<< std::endl;
      const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
      //watch out: only ILU might work right now because of compute nullspace might not work...? ... test!
      switch (prectyp)
      {
      case INPAR::SOLVER::azprec_ML:
      case INPAR::SOLVER::azprec_MLfluid:
      case INPAR::SOLVER::azprec_MLAPI:
      case INPAR::SOLVER::azprec_MLfluid2:
      case INPAR::SOLVER::azprec_MueLuAMG_sym:
      case INPAR::SOLVER::azprec_MueLuAMG_nonsym:
      {
        if(projconstr_||proj_)
        { //has 3 dofs, velocity dofs
          //BUT: enriched nodes have 8 dofs, so we have to calculate our own nullspace for 3 dofs
          // store nv and np at unique location in solver parameter list
          solver_->Params().sublist("NodalBlockInformation").set("nv",3);
          solver_->Params().sublist("NodalBlockInformation").set("np",0);
          solver_->Params().sublist("NodalBlockInformation").set("numdf",3);
          solver_->Params().sublist("NodalBlockInformation").set("dimns",3);
          // get the belos list and see whether we use downwinding
          Teuchos::ParameterList& beloslist = solver_->Params().sublist("Belos Parameters");

          beloslist.set<int>("downwinding nv",3);
          beloslist.set<int>("downwinding np",0);

          Teuchos::ParameterList* mllist_ptr = NULL;
          mllist_ptr = &((solver_->Params()).sublist("ML Parameters"));
          Teuchos::ParameterList& mllist = *mllist_ptr; //solveparams.sublist("ML Parameters");

          Teuchos::RCP<std::vector<double> > ns = mllist.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);
          if (ns != Teuchos::null) break;

          // no, we have not previously computed the nullspace
          // -> compute nullspace
          ns = Teuchos::null;

          mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);
          // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
          mllist.set<bool>("ML validate parameter list",false);

          mllist.set("PDE equations",3);
          mllist.set("null space: dimension",3);
          mllist.set("null space: type","pre-computed");
          mllist.set("null space: add default vectors",false);

          // allocate dimns times the local length of the rowmap
          const int lrows = enrdofrowmap_->NumMyElements();
          ns = Teuchos::rcp(new std::vector<double>(3*lrows));
          double* nullsp = &((*ns)[0]);
          mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
          mllist.set("null space: vectors",nullsp);

          //now compute the nullspace
          double* mode[6];
          for (int i=0; i<3; ++i)
            mode[i] = &((*ns)[i*lrows]);

          for (int i=0; i<xwdiscret_->NumMyRowNodes(); ++i)
          {
            const unsigned int ndof = 3;
            DRT::Node* actnode = xwdiscret_->lRowNode(i);
            if(!actnode)
              dserror("cannot find node");
            std::vector<int> dofs = xwdiscret_->Dof(0,actnode);
            std::vector<int> actdofs;
            //only dof 4...6 (enriched dofs)
            actdofs.push_back(dofs.at(4));
            actdofs.push_back(dofs.at(5));
            actdofs.push_back(dofs.at(6));

            for (unsigned j=0; j<ndof; ++j)
            {
              const int dof = actdofs.at(j);

              const int lid = enrdofrowmap_->LID(dof);
              if (lid<0) dserror("Cannot find dof %i",dof);

              for (unsigned k=0; k<ndof; ++k)
              {
                if (k%ndof == j%ndof)
                  mode[k%ndof][lid] = 1.0;
                else
                  mode[k%ndof][lid] = 0.0;
              }
            } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
          } // for (int i=0; i<NumMyRowNodes(); ++i)
        }
      }
      break;
      case INPAR::SOLVER::azprec_ILU:
      case INPAR::SOLVER::azprec_ILUT:
      case INPAR::SOLVER::azprec_SymmGaussSeidel:
      case INPAR::SOLVER::azprec_Jacobi:
        // do nothing
      break;
      default:
        dserror("You have to choose ML, MueLu or ILU preconditioning");
      break;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Routine to update Tauw                                     bk 07/14 |
 *----------------------------------------------------------------------*/
  void FLD::XWall::UpdateTauW(int step,Teuchos::RCP<Epetra_Vector>   trueresidual, int itnum,Teuchos::RCP<Epetra_Vector>   accn,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   veln, FluidImplicitTimeInt& fluid)
{
  Teuchos::RCP<Epetra_Vector> newtauw = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));

  if(((tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step >= switch_step_) || (tauwcalctype_ != INPAR::FLUID::gradient_to_residual && step > 1)) && smooth_res_aggregation_)
  {
    //filtering matrix for wall shear stress
    if(SepEnr_==Teuchos::null)
    {
      if(myrank_==0)
        std::cout << "build filtering matrix for residual via aggregation:" << std::endl;

      MLAPI::Init();
      const int scale_sep_solvernumber = params_->sublist("WALL MODEL").get<int>("ML_SOLVER");
      if(scale_sep_solvernumber==-1)
        dserror("Set a solver number in ML_SOLVER for smooting of tauw via aggregation");
        Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(scale_sep_solvernumber),
                                              discret_->Comm(),
                                              DRT::Problem::Instance()->ErrorFile()->Handle()));

      Teuchos::ParameterList& mlparams = solver->Params().sublist("ML Parameters");
      // compute the null space,
      discret_->ComputeNullSpaceIfNecessary(solver->Params(),true);

      // get nullspace parameters
      double* nullspace = mlparams.get("null space: vectors",(double*)NULL);
      if (!nullspace) dserror("No nullspace supplied in parameter list");
      int nsdim = mlparams.get("null space: dimension",1);
      if(nsdim!=4)
        dserror("Wrong Nullspace dimension for XWall");
      int lrowdofs = discret_->DofRowMap()->NumMyElements();
    //  std::cout << "lrowdofs  " << lrowdofs << std::endl;
    //std::cout << "check the nullspace for mfs" << std::endl;
      for (int j=0; j<discret_->NodeRowMap()->NumMyElements();++j)
      {
        int xwallgid = discret_->NodeRowMap()->GID(j);

        if (not discret_->NodeRowMap()->MyGID(xwallgid)) //just in case
          dserror("not on proc");
        {
          DRT::Node* xwallnode = discret_->gNode(xwallgid);
          if(!xwallnode) dserror("Cannot find node");

          int firstglobaldofid=discret_->Dof(xwallnode,0);
          int firstlocaldofid=discret_->DofRowMap()->LID(firstglobaldofid);

          std::vector<DRT::Condition*> nodedircond;
          xwallnode->GetCondition("FluidStressCalc",nodedircond);

          if(not nodedircond.empty())
          {
            //these nodes are wall nodes, so aggregate them
            nullspace[firstlocaldofid]=1.0;
            nullspace[lrowdofs+firstlocaldofid+1]=1.0;
            nullspace[lrowdofs*2+firstlocaldofid+2]=1.0;
            nullspace[lrowdofs*3+firstlocaldofid+3]=1.0;
          }
          else
          {
            //set everything to zero
            nullspace[firstlocaldofid]=0.0;
            nullspace[lrowdofs+firstlocaldofid+1]=0.0;
            nullspace[lrowdofs*2+firstlocaldofid+2]=0.0;
            nullspace[lrowdofs*3+firstlocaldofid+3]=0.0;
          }

          if(discret_->NumDof(xwallnode)>5)
          {
            nullspace[firstlocaldofid+4]=0.0;
            nullspace[lrowdofs+firstlocaldofid+5]=0.0;
            nullspace[lrowdofs*2+firstlocaldofid+6]=0.0;
            nullspace[lrowdofs*3+firstlocaldofid+7]=0.0;
          }
        }
      }

      // get plain aggregation Ptent
      Teuchos::RCP<Epetra_CrsMatrix> crsPtent;
      MLAPI::GetPtent(*fluid.SystemMatrix()->EpetraMatrix(),mlparams,nullspace,crsPtent);
      LINALG::SparseMatrix Ptent(crsPtent);

      // compute scale-separation matrix: S = I - Ptent*Ptent^T
      SepEnr_ = LINALG::Multiply(Ptent,false,Ptent,true);

      SepEnr_->Complete();
    }

    Teuchos::RCP<Epetra_Vector> fsresidual = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
    SepEnr_->Multiply(false,*trueresidual,*fsresidual);
    trueresidual->Update(1.0,*fsresidual,0.0);
  }
  switch (tauwtype_)
  {
    case INPAR::FLUID::constant: //works
    {
      if(step==1)
        CalcMK();
      //tauw_ is constant and inctauw_ is zero
      return;
    }
    break;
    case INPAR::FLUID::mean_between_steps:
    {
      dserror("calculation of tauw not supported right now");
    }
    break;
    case INPAR::FLUID::mean_iter: //works
    {
      dserror("calculation of tauw not supported right now");
    }
    break;
    case INPAR::FLUID::between_steps:
    {
      inctauw_->PutScalar(0.0);

      if(itnum==0)//in between steps
        CalcTauW(step,trueresidual, velnp,fluid.CalcWallShearStresses());
      else
        return;
    }
    break;
    case INPAR::FLUID::fix_point_iter_with_step_control:
    {
      dserror("calculation of tauw not supported right now");
    }
    break;
    case INPAR::FLUID::fully_linearized:
    {
      dserror("not yet implemented");
    }
    break;
    default: dserror("unknown tauwtype_");
    break;
  }

  LINALG::Export(*tauw_,*newtauw);

  double actmean=-1.0;
  newtauw->MeanValue(&actmean);

  double min=-1.0;
  newtauw->MinValue(&min);
  if(min<1e-10)
    dserror("tauw is zero");
  double max=-1.0;
  newtauw->MaxValue(&max);

  //convergence check, works only if we don't take the mean
  newtauw->PutScalar(0.0);
  LINALG::Export(*inctauw_,*newtauw);
  newtauw->Norm2(&inctauwnorm_);
  //rescale inctauw to full increment
  if(fac_>1.0e-8)
    inctauwnorm_/=fac_;

  if(myrank_==0)
    std::cout << "  min:  "<<min << "  max:  "<< max << "  mean-applied:  "<< actmean<< "  inc norm2:  "<< inctauwnorm_ ;

  //also project the other vectors
  //later I could also run all three at the same time (with the same matrix)

  if(proj_)
  {
    if(myrank_==0)
      std::cout << "  L2-project... ";
    if(tauwtype_==INPAR::FLUID::between_steps)
    {
      L2ProjectVector(veln,Teuchos::null,accn);

      //at the beginning of this time step they are equal -> calculate only one of them
      velnp->Update(1.0,*veln,0.0);
    }
    else
      L2ProjectVector(veln,velnp,accn);

    if(myrank_==0)
      std::cout << "done!" << std::endl;
  }
  else if(projconstr_)
  {
    if(myrank_==0)
      std::cout << "  L2-project wcc... ";
    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(9), "SOLVER");

    if(solvertype == INPAR::SOLVER::umfpack||solvertype == INPAR::SOLVER::aztec_msr)
    {
      L2ProjectVectorWithContinuityConstraint(veln,Teuchos::null,Teuchos::null);
      L2ProjectVectorWithContinuityConstraint(accn,Teuchos::null,Teuchos::null);
      //at the beginning of this time step they are equal -> calculate only one of them
      velnp->Update(1.0,*veln,0.0);
    }
    else if(tauwtype_==INPAR::FLUID::between_steps)
    {
      L2ProjectVectorWithContinuityConstraint(veln,Teuchos::null,accn);
      //at the beginning of this time step they are equal -> calculate only one of them
      velnp->Update(1.0,*veln,0.0);
    }
    else
      L2ProjectVectorWithContinuityConstraint(veln,velnp,accn);

    if(myrank_==0)
      std::cout << "done!" << std::endl;
  }
  else
    std::cout << std::endl;

  CalcMK();

  return;
}

/*----------------------------------------------------------------------*
 |  Routines to calculate Tauw                                 bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::CalcTauW(int step, Teuchos::RCP<Epetra_Vector>   trueresidual,Teuchos::RCP<Epetra_Vector>   velnp, Teuchos::RCP<Epetra_Vector>   wss)
{

  Teuchos::RCP<Epetra_Vector> newtauw = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));
  Teuchos::RCP<Epetra_Vector> newtauw2 = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));
  Teuchos::RCP<Epetra_Vector> tauw=Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()),true));

  if(tauwcalctype_ == INPAR::FLUID::gradient_to_residual && switch_step_ == step && myrank_ == 0)
    std::cout << "\n switching from gradient to residual \n" << std::endl;

  if(tauwcalctype_==INPAR::FLUID::spalding)
  {
    Teuchos::RCP<Epetra_Vector> velnpxw = Teuchos::rcp(new Epetra_Vector(*(discret_->DofColMap()),true));

    LINALG::Export(*velnp,*velnpxw);

    for (int j=0; j<dircolnodemap_->NumMyElements();++j)
    {
      int xwallgid = dircolnodemap_->GID(j);

      if (discret_->NodeRowMap()->MyGID(xwallgid)) //just in case
      {
        DRT::Node* xwallnode = discret_->gNode(xwallgid);
        if(!xwallnode) dserror("Cannot find node");

        DRT::Element** surrele = xwallnode->Elements();

        for (int k=0;k<xwallnode->NumElement();++k)
        {
          DRT::Node** surrnodes=surrele[k]->Nodes();

          for(int i=0; i<surrele[k]->NumNode();i++)
          {
            //now I need the velocity at this node,
            //the distance from the wall
            //and the spalding's law

            if(!surrnodes[i])
              dserror("can't find surrounding node");

            std::vector<DRT::Condition*> nodedircond;
            surrnodes[i]->GetCondition("FluidStressCalc",nodedircond);

            if(nodedircond.empty())
            {
              int firstgdofid=discret_->Dof(0,surrnodes[i],0);

              int firstldofid=velnpxw->Map().LID(firstgdofid);
        //      std::cout << der2psiold/der2psinew << std::endl;
              double velx =(*velnpxw)[firstldofid];
              double vely =(*velnpxw)[firstldofid+1];
              double velz =(*velnpxw)[firstldofid+2];
              double u=sqrt(velx*velx+vely*vely+velz*velz); //we approximate the wall-parallel velocity
                                                            //by its norm

              int ylid=wdist_->Map().LID(surrnodes[i]->Id());
              double y=(*wdist_)[ylid];

              double k=0.409836066;
              double B=5.17;

              int tauwlid=tauw_->Map().LID(surrnodes[i]->Id());
              double tauw=(*tauw_)[tauwlid];//use old value for a start...
              if(tauw<1e-6)
                tauw=1.0;
              double yp=0.0;
              double up=0.0;
              double derup=0.0;
              double utau=0.0;
              double derutau=0.0;
              double inc=100.0;
              double fn=0.0;
              double dfn=0.0;

              int count=0;
              while (abs(inc)>1e-8)
              {

                utau=sqrt(tauw/dens_);
                derutau=1.0/(2.0*dens_*tauw);

                yp=utau*y/visc_;
                up=u/utau;

                derup=-u/(utau*utau)*derutau;
                fn=-yp+up+exp(-k*B)*(exp(k*up)-1.0-k*up-(k*up)*(k*up)*0.5-(k*up)*(k*up)*(k*up)/6.0);
                dfn=-y/visc_*derutau+derup+exp(-k*B)*(k*exp(k*up)*derup-k*derup-k*k*up*derup-k*k*k*up*up*derup*0.5);

                if(abs(dfn)<1.0e-10||abs(dfn)>1e100||std::isinf(dfn)!=0)
                {
                  inc=-2.0e-8;//dummy value
                }
                else
                  inc=fn/dfn;

                if(tauw-inc<0.0|| count>50)
                  inc/=2.0; //damping for robustness

                tauw-=inc;
                if(tauw<1.0e-5)
                  tauw+=0.001*abs(inc); //it is problematic if tauw approaches zero, thus go a step away of zero

                tauw=abs(tauw);//just in case, because it can happen, that it becomes negative during newton iteration

                count++;
                //this converges very slowly... might there be a mistake somewhere?
                if(count>100000)
                {
                  std::cout << "Tauw does not converge... u= " << u << "  tauw= " << tauw << "  inc= " << inc<< "  dfn= " << dfn << std::endl;
                  inc=0.0;
                }
              }

              if(tauw < min_tauw_)
                tauw = min_tauw_;

              //store in vector
              int err =     newtauw->ReplaceGlobalValue(xwallgid,0,tauw);
              if(err==1)
                dserror("global row not on proc");
              else if(err==-1)
                dserror("wrong vector index");
            }
          }
        }
      }
    }
  }
  else if(tauwcalctype_==INPAR::FLUID::residual || (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step >= switch_step_))
  {
    for(int lnodeid=0;lnodeid<dircolnodemap_->NumMyElements();lnodeid++)
    {
      int gid = dircolnodemap_->GID(lnodeid);
      //continue only if on this proc
      if (discret_->NodeRowMap()->MyGID(gid))
      {
        DRT::Node* node = discret_->gNode(gid);
         if (!node) dserror("ERROR: Cannot find off wall node with gid %",gid);

         int firstglobaldofid=discret_->Dof(0,node,0);
         int firstlocaldofid=wss->Map().LID(firstglobaldofid);

         if (firstlocaldofid < 0)
           dserror("localdofid not found in map for given globaldofid");
         double forcex=(*wss)[firstlocaldofid];
         double forcey=(*wss)[firstlocaldofid+1];
         double forcez=(*wss)[firstlocaldofid+2];

         double tauw=sqrt(forcex*forcex+forcey*forcey+forcez*forcez);

         //this is necessary since we are dividing by tauw on element level
         //also, the shape functions become singular, if tauw==0
         if(tauw < min_tauw_)
           tauw = min_tauw_;
         //store in vector
         int err =     newtauw->ReplaceGlobalValue(gid,0,tauw);
         if(err!=0)
           dserror("something went wrong during replacemyvalue");
      }
    }
  }
  else if(tauwcalctype_==INPAR::FLUID::gradient || (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step < switch_step_))
  {
    //necessary to set right state (the maps of the state vector and discretization have to be equal)
    Teuchos::RCP<Epetra_Vector> statevel = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()),true));
    Teuchos::RCP<Epetra_Vector> newtauwxwdis = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeRowMap()),true));
    LINALG::Export(*velnp,*statevel);

    xwdiscret_->SetState("vel",statevel);

    Teuchos::RCP<Epetra_Vector> timesvec = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeRowMap()),true));

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    // get number of elements
    const int numele = xwdiscret_->NumMyColElements();

    // loop column elements: vector
    for (int i=0; i<numele; ++i)
    {
      DRT::Element* actele = xwdiscret_->lColElement(i);

      const int numnode = actele->NumNode();

      // get element location vector and ownerships
      actele->LocationVector(*xwdiscret_,lm,lmowner,lmstride);

      elevector1.Size(numnode);
      elevector2.Size(numnode);
      // set action in order to project element void fraction to nodal void fraction
      Teuchos::ParameterList params;
      SetXWallParamsXWDis(params);
      params.set<int>("action",FLD::tauw_via_gradient);
      params.set("newtauw",newtauwxwdis);
      // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
      //elevector1 has to be NULL here, because I am assuming a dof-based vector otherwise
      actele->Evaluate(params,*xwdiscret_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      // get element location vector for nodes
      lm.resize(numnode);
      lmowner.resize(numnode);

      DRT::Node** nodes = actele->Nodes();
      for(int n=0; n<numnode; ++n)
      {
        lm[n] = nodes[n]->Id();
        lmowner[n] = nodes[n]->Owner();
      }

      // assembling into node maps
      LINALG::Assemble(*newtauwxwdis,elevector1,lm,lmowner);
      LINALG::Assemble(*timesvec,elevector2,lm,lmowner);
    } //end element loop

    xwdiscret_->ClearState();

    //scale with times:
    for (int l=0;l<(xwdiscret_->NodeRowMap())->NumMyElements();l++)
    {
      double sumnewtauw=(*newtauwxwdis)[l];
      double timesfac=(*timesvec)[l];
      double newtauwsc=1.0;
      //we only have to do something, if we assembled at least once the same value
      if(timesfac>0.5)
      {
        newtauwsc=sumnewtauw/timesfac;

        if(newtauwsc < min_tauw_)
          newtauwsc = min_tauw_;
        int err  =     newtauwxwdis->ReplaceMyValue(l,0,newtauwsc);
        if(err!=0)
          dserror("something went wrong during replacemyvalue");
      }
    }
    LINALG::Export(*newtauwxwdis,*newtauw);
  }
  else
    dserror("unknown tauwcalctype_");

  tauw->Update(1.0,*tauw_,0.0);


  inctauw_->Update(1.0,*tauw,0.0);
  tauw->PutScalar(0.0);

  tauwcouplingmattrans_->Multiply(true,*newtauw,*newtauw2);


  double meansp=0.0;
  newtauw2->MeanValue(&meansp);

  //here i can modify newtauw2 to get rid of ramp functions!
  BlendingViaModificationOfTauw(newtauw2);

  LINALG::Export(*newtauw2,*tauw);
  inctauw_->Update(fac_,*tauw,-fac_); //now this is the increment (new-old)

  tauw_->Update(1.0,*inctauw_,1.0);

  LINALG::Export(*inctauw_,*newtauw2);
  LINALG::Export(*newtauw2,*inctauwxwdis_);
  LINALG::Export(*tauw_,*newtauw2);
  LINALG::Export(*newtauw2,*tauwxwdis_);

  if(meansp<2.0e-9)
    dserror("Average wall shear stress is zero. You probably forgot to specify approprite DESIGN FLUID STRESS CALC SURF CONDITIONS where the stress should be calculated.");

  if(myrank_==0)
  std::cout << "tauw mean:  "  << meansp ;

  return;
}


/*----------------------------------------------------------------------*
 |  Modify Tauw so that ramp functions are not necessary       bk 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::BlendingViaModificationOfTauw(Teuchos::RCP<Epetra_Vector>   newtauw2)
{
  if(blendingtype_==INPAR::FLUID::tauw_transformation)
  {
    Teuchos::RCP<Epetra_Vector> yp = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_,true));
    double sumyp=0.0;
    int count=0;
    for (int j=0; j<xwallrownodemap_->NumMyElements();++j)
    {
      double toggle=(*xtoggleloc_)[j];
      if(toggle>0.5&&toggle<0.9)
      {
        double wdist=(*walldist_)[j];
        double tauw=(*newtauw2)[j];

        double utau=sqrt(tauw/dens_);
        double yplus=wdist*utau/visc_;
        sumyp+=yplus;
        count++;
        int err  =     yp->ReplaceMyValue(j,0,yplus);
        if(err!=0)
          dserror("cannot replace my value");
      }
    }
    double max=0.0;
    yp->MaxValue(&max);
    double sumypall;
    int countall;
    (discret_->Comm()).SumAll(&count,&countall,1);
    (discret_->Comm()).SumAll(&sumyp,&sumypall,1);
    double ypint=0.0*max+1.0*sumypall/(double)countall;

    //now calculate respective tauw and replace in newtauw2
    for (int j=0; j<xwallrownodemap_->NumMyElements();++j)
    {
      double toggle=(*xtoggleloc_)[j];
      if(toggle>0.5&&toggle<0.9)
      {
        double wdist=(*walldist_)[j];
        if(wdist<1e-9)
          dserror("blending nodes may not be wall nodes if transformation is used");
        double tauw=dens_*ypint*ypint*visc_*visc_/(wdist*wdist);
        int err  =     newtauw2->ReplaceMyValue(j,0,tauw);
        if(err!=0)
          dserror("cannot replace my value");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  L2-project enriched dofs of vector                         bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::L2ProjectVector(Teuchos::RCP<Epetra_Vector>   veln,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   accn)
{

  if(not veln->Map().SameAs(*discret_->DofRowMap()))
    dserror("input map is not the dof row map of the fluid discretization");

  massmatrix_->Zero();

  incveln_->PutScalar(0.0);
  if(accn!=Teuchos::null)
    incaccn_->PutScalar(0.0);
  if(velnp!=Teuchos::null)
    incvelnp_->PutScalar(0.0);

  LINALG::Export(*veln,*stateveln_);
  if(accn!=Teuchos::null)
    LINALG::Export(*accn,*stateaccn_);
  if(velnp!=Teuchos::null)
    LINALG::Export(*velnp,*statevelnp_);

  //number of right hand sides during solving
  //is the number of velocity components that is solved for
  //3 since we are in 3D
  int numberofrhs=0;
  if(velnp==Teuchos::null&&accn==Teuchos::null)
    numberofrhs=1;
  else if(velnp==Teuchos::null||accn==Teuchos::null)
    numberofrhs=2;
  else
    numberofrhs=3;

  xwdiscret_->SetState("veln",stateveln_);
  if(accn!=Teuchos::null)
    xwdiscret_->SetState("accn",stateaccn_);
  if(velnp!=Teuchos::null)
    xwdiscret_->SetState("velnp",statevelnp_);

  // set action in order to project nodal enriched values to new shape functions
  Teuchos::ParameterList params;
  SetXWallParamsXWDis(params);
  params.set<int>("action",FLD::xwall_l2_projection);

  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhsassemble = Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_, numberofrhs));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevectordummy;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get number of elements
  const int numele = xwdiscret_->NumMyColElements();

  // loop column elements
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = xwdiscret_->lColElement(i);

    const int numnode = actele->NumNode();
    const int numdf=3;

    // get element location vector and ownerships
    actele->LocationVector(*xwdiscret_,lm,lmowner,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elematrix1.Shape(numnode*numdf,numnode*numdf);
    // Reshape element matrices and vectors and initialize to zero
    elematrix2.Shape(numnode*numdf,numberofrhs);//we have 3 right hand sides for now: 3 velocity components

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    actele->Evaluate(params,*xwdiscret_,lm,elematrix1,elematrix2,elevectordummy,elevector2,elevector3);

    // get element location vector for enriched dofs
    std::vector<int> lmassemble;
    std::vector<int> lmownerassemble;
    lmassemble.resize(numnode*numdf);
    lmownerassemble.resize(numnode*numdf);

    for(int n=0; n<numnode; ++n)
    {
      for(int df=4;df<7;++df)
      {
        lmassemble[n*numdf+df-4] = lm[n*8+df];
        lmownerassemble[n*numdf+df-4] = lmowner[n*8+df];
      }
    }

    // assembling into node maps
    massmatrix_->Assemble(actele->Id(),elematrix1,lmassemble,lmownerassemble);

    // assembling into node maps
    // assemble numberofrhs entries in rhs vector sequentially
    elevector1.Size(numnode*numdf);
    for(int n=0; n<numberofrhs; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for(int idf=0; idf<numnode*numdf; ++idf)
        elevector1(idf) = elematrix2(idf,n);
      // assemble into nth vector of MultiVector
      LINALG::Assemble(*rhsassemble,n,elevector1,lmassemble,lmownerassemble);
    }
  } //end element loop

  xwdiscret_->ClearState();
  // finalize the matrix
  massmatrix_->Complete();

  // solution vector
  Teuchos::RCP<Epetra_MultiVector> resultvec = Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_,numberofrhs));

  // solve for 1, 2 or 3 right hand sides at the same time --> thanks to Belos
  solver_->Solve(massmatrix_->EpetraOperator(), resultvec, rhsassemble, true, true, Teuchos::null);

  //now copy result in original vector: the result is an increment of the velocity/ acceleration
  LINALG::Export(*((*resultvec)(0)),*incveln_);
  if(numberofrhs>1)
  LINALG::Export(*((*resultvec)(1)),*incaccn_);
  if(numberofrhs>2)
    LINALG::Export(*((*resultvec)(2)),*incvelnp_);

  veln->Update(1.0,*incveln_,1.0);
  if(accn!=Teuchos::null)
    accn->Update(1.0,*incaccn_,1.0);
  if(velnp!=Teuchos::null)
    velnp->Update(1.0,*incvelnp_,1.0);


  return;
}

/*----------------------------------------------------------------------*
 |  L2-project enriched dofs of vector                         bk 07/14 |
 |  with a penalty continuity constraint                                |
 *----------------------------------------------------------------------*/
void FLD::XWall::L2ProjectVectorWithContinuityConstraint(Teuchos::RCP<Epetra_Vector>   veln,Teuchos::RCP<Epetra_Vector>   velnp,Teuchos::RCP<Epetra_Vector>   accn)
{

  if(not veln->Map().SameAs(*discret_->DofRowMap()))
    dserror("input map is not the dof row map of the fluid discretization");

  massmatrix_->Zero();

  incveln_->PutScalar(0.0);
  if(accn!=Teuchos::null)
    incaccn_->PutScalar(0.0);
  if(velnp!=Teuchos::null)
    incvelnp_->PutScalar(0.0);

  LINALG::Export(*veln,*stateveln_);
  if(accn!=Teuchos::null)
    LINALG::Export(*accn,*stateaccn_);
  if(velnp!=Teuchos::null)
    LINALG::Export(*velnp,*statevelnp_);

  //number of right hand sides during solving
  //is the number of velocity components that is solved for
  //3 since we are in 3D
  int numberofrhs=0;
  if(velnp==Teuchos::null&&accn==Teuchos::null)
    numberofrhs=1;
  else if(velnp==Teuchos::null||accn==Teuchos::null)
    numberofrhs=2;
  else
    numberofrhs=3;

  xwdiscret_->SetState("veln",stateveln_);
  if(accn!=Teuchos::null)
    xwdiscret_->SetState("accn",stateaccn_);
  if(velnp!=Teuchos::null)
    xwdiscret_->SetState("velnp",statevelnp_);

  // set action in order to project nodal enriched values to new shape functions
  Teuchos::ParameterList params;
  SetXWallParamsXWDis(params);
  params.set<int>("action",FLD::xwall_l2_projection_with_continuity_constraint);

  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhsassemble = Teuchos::rcp(new Epetra_MultiVector(*mergedmap_, numberofrhs));

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevectordummy;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get number of elements
  const int numele = xwdiscret_->NumMyColElements();

  // loop column elements
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = xwdiscret_->lColElement(i);
    DRT::ELEMENTS::FluidXWall* xwallele=dynamic_cast<DRT::ELEMENTS::FluidXWall*>(actele);
    if(!xwallele)//not a xwall element and the node row maps don't know it's nodes
      dserror("must be xwallele");
    const int numnode = actele->NumNode();
    const int numdf=4;

    // get element location vector and ownerships
    actele->LocationVector(*xwdiscret_,lm,lmowner,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elematrix1.Shape(numnode*numdf,numnode*numdf);
    // Reshape element matrices and vectors and initialize to zero
    elematrix2.Shape(numnode*numdf,numberofrhs);//we have 3 right hand sides for now: 3 velocity components

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    actele->Evaluate(params,*xwdiscret_,lm,elematrix1,elematrix2,elevectordummy,elevector2,elevector3);

    // get element location vector for enriched dofs
    std::vector<int> lmassemble;
    std::vector<int> lmownerassemble;
    lmassemble.resize(numnode*numdf);
    lmownerassemble.resize(numnode*numdf);

    for(int n=0; n<numnode; ++n)
    {
      for(int df=4;df<8;++df)
      {
        lmassemble[n*numdf+df-4] = lm[n*8+df];
        lmownerassemble[n*numdf+df-4] = lmowner[n*8+df];
      }
    }

    // assembling into node maps
    massmatrix_->Assemble(actele->Id(),elematrix1,lmassemble,lmownerassemble);

    // assembling into node maps
    // assemble numberofrhs entries in rhs vector sequentially
    elevector1.Size(numnode*numdf);
    for(int n=0; n<numberofrhs; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for(int idf=0; idf<numnode*numdf; ++idf)
        elevector1(idf) = elematrix2(idf,n);
      // assemble into nth vector of MultiVector
//      std::cout << elevector1 << std::endl;
      LINALG::Assemble(*rhsassemble,n,elevector1,lmassemble,lmownerassemble);
    }
  } //end element loop

  xwdiscret_->ClearState();
  // finalize the matrix
  massmatrix_->Complete();

  // we want to split M into 2 groups = 4 blocks
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix_solve, B, ktemp1, ktemp2;
  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhssolve = Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_, numberofrhs));

  // split into enriched and lagrange multiplyer dof's
  LINALG::SplitMatrix2x2(massmatrix_,enrdofrowmap_,lagrdofrowmap_,enrdofrowmap_,lagrdofrowmap_,massmatrix_solve,ktemp1,B,ktemp2);

  Teuchos::RCP<LINALG::SparseMatrix> kbtb_add = MLMultiply(*B,true,*B,false,false,false,true);

  // we want to split f into 2 groups
  Teuchos::RCP<Epetra_Vector> fm1, fb1,fb1_add;
  // do the vector splitting
  LINALG::SplitVector(*mergedmap_,*((*rhsassemble)(0)),enrdofrowmap_,fm1,lagrdofrowmap_,fb1);

  fb1_add=Teuchos::rcp(new Epetra_Vector(*enrdofrowmap_,true));
  B->Multiply(true,*fb1,*fb1_add);

  //determine penalty parameter
  double normconstr=0.0;
  fb1_add->Norm2(&normconstr);
  double normeq=0.0;
  fm1->Norm2(&normeq);

  double k=0.0;
  if(normconstr>1e-8)
    k= penalty_param_*normeq/normconstr;
  //else use 0.0, no continuity constraint
  if(myrank_==0)
    std::cout << "k = " << k << " ... " ;

//  std::cout << *B << std::endl;
  massmatrix_solve->UnComplete();
  massmatrix_solve->Add(*kbtb_add,false,k,1.0);

  massmatrix_solve->Complete();

  fm1->Update(k,*fb1_add,1.0);

  ((*rhssolve)(0))->Update(1.0,*fm1,0.0);

  // do the vector splitting smn -> sm+n
  if(numberofrhs>1)
  {
    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fm2, fb2,fb2_add;
    fb2_add=Teuchos::rcp(new Epetra_Vector(*enrdofrowmap_,true));
    LINALG::SplitVector(*mergedmap_,*(*rhsassemble)(1),enrdofrowmap_,fm2,lagrdofrowmap_,fb2);
    B->Multiply(true,*fb2,*fb2_add);
    fm2->Update(k,*fb2_add,1.0);
    ((*rhssolve)(1))->Update(1.0,*fm2,0.0);
  }

  // do the vector splitting smn -> sm+n
  if(numberofrhs>2)
  {
    // we want to split f into 3 groups s.m,n
    Teuchos::RCP<Epetra_Vector> fm2, fb2,fb2_add;
    fb2_add=Teuchos::rcp(new Epetra_Vector(*enrdofrowmap_,true));
    LINALG::SplitVector(*mergedmap_,*(*rhsassemble)(2),enrdofrowmap_,fm2,lagrdofrowmap_,fb2);
    B->Multiply(true,*fb2,*fb2_add);
    fm2->Update(k,*fb2_add,1.0);
    ((*rhssolve)(2))->Update(1.0,*fm2,0.0);
  }

  // solution vector
  Teuchos::RCP<Epetra_MultiVector> resultvec = Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_,numberofrhs));

  // solve for 1, 2 or 3 right hand sides at the same time --> thanks to Belos
  solver_->Solve(massmatrix_solve->EpetraOperator(), resultvec, rhssolve, true, true, Teuchos::null);

  //now copy result in original vector: the result is an increment of the velocity/ acceleration
  LINALG::Export(*((*resultvec)(0)),*incveln_);
  if(numberofrhs>1)
  LINALG::Export(*((*resultvec)(1)),*incaccn_);
  if(numberofrhs>2)
    LINALG::Export(*((*resultvec)(2)),*incvelnp_);

  veln->Update(1.0,*incveln_,1.0);
  if(accn!=Teuchos::null)
    accn->Update(1.0,*incaccn_,1.0);
  if(velnp!=Teuchos::null)
    velnp->Update(1.0,*incvelnp_,1.0);

  return;
}

/*----------------------------------------------------------------------*
 |  Adapt ML Nullspace for MFS aggregation                     bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::AdaptMLNullspace(const Teuchos::RCP<LINALG::Solver>&   solver)
{
  // extract the ML parameters:
  Teuchos::ParameterList&  mlparams = solver->Params().sublist("ML Parameters");

  // get nullspace parameters
  double* nullspace = mlparams.get("null space: vectors",(double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int nsdim = mlparams.get("null space: dimension",1);
  if(nsdim!=4)
    dserror("Wrong Nullspace dimension for XWall");
  int lrowdofs = discret_->DofRowMap()->NumMyElements();
//  std::cout << "lrowdofs  " << lrowdofs << std::endl;
//std::cout << "check the nullspace for mfs" << std::endl;
  for (int j=0; j<xwallrownodemap_->NumMyElements();++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (not discret_->NodeRowMap()->MyGID(xwallgid)) //just in case
      dserror("not on proc");
    {
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if(!xwallnode) dserror("Cannot find node");

      int firstglobaldofid=discret_->Dof(xwallnode,0);
      int firstlocaldofid=discret_->DofRowMap()->LID(firstglobaldofid);

      nullspace[firstlocaldofid+4]=0.0;
      nullspace[lrowdofs+firstlocaldofid+5]=0.0;
      nullspace[lrowdofs*2+firstlocaldofid+6]=0.0;
      nullspace[lrowdofs*3+firstlocaldofid+7]=0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Zero Enriched dofs of row vector                           bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::ZeroEnrDofs(Teuchos::RCP<Epetra_Vector>   vec)
{
  for (int j=0; j<xwallrownodemap_->NumMyElements();++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (not discret_->NodeRowMap()->MyGID(xwallgid)) //just in case
      dserror("not on proc");
    {
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if(!xwallnode) dserror("Cannot find node");

      int firstgdofid=discret_->Dof(xwallnode,0);

      int err =     vec->ReplaceGlobalValue(firstgdofid+4,0,0.0);
          err +=     vec->ReplaceGlobalValue(firstgdofid+5,0,0.0);
          err +=     vec->ReplaceGlobalValue(firstgdofid+6,0,0.0);
          err +=     vec->ReplaceGlobalValue(firstgdofid+7,0,0.0);

       if (err!=0)
         dserror("something went wrong");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Calculate MK for residual-based stabilization parameter    bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::CalcMK()
{
  Teuchos::RCP<Epetra_MultiVector> mkxw=Teuchos::rcp(new Epetra_MultiVector(*(xwdiscret_->ElementRowMap()),1,true));

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",FLD::xwall_calc_mk);

  SetXWallParamsXWDis(eleparams);

  xwdiscret_->EvaluateScalars(eleparams, mkxw);


  Teuchos::RCP<Epetra_Vector> mkxwv=Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->ElementRowMap()),true));
  Teuchos::RCP<Epetra_Vector> mkv=Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()),true));

  //export
  LINALG::Export(*((*mkxw)(0)),*mkxwv);
  LINALG::Export(*mkxwv,*mkxwstate_);
  LINALG::Export(*mkxwv,*mkv);
  LINALG::Export(*mkv,*mkstate_);


  return;
} // end CalcMK

/*----------------------------------------------------------------------*
 |  Filter residual for output of wall shear stress            bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::FilterResVectorForOutput(Teuchos::RCP<Epetra_Vector>   residual)
{
  if(SepEnr_!=Teuchos::null)
  {
  Teuchos::RCP<Epetra_Vector> fsresidual = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
  SepEnr_->Multiply(false,*residual,*fsresidual);
  residual->Update(1.0,*fsresidual,0.0);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Write enriched dofs in standard dofs for output            bk 09/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::XWall::GetOutputVector(Teuchos::RCP<Epetra_Vector>   vel)
{
  Teuchos::RCP<Epetra_Vector> velenr = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
  for (int i=0;i<xwallrownodemap_->NumMyElements();++i)
  {
    int xwallgid = xwallrownodemap_->GID(i);
    DRT::Node* xwallnode = discret_->gNode(xwallgid);
    if(!xwallnode) dserror("Cannot find node");

    int firstglobaldofid=discret_->Dof(xwallnode,0);
    int firstlocaldofid=discret_->DofRowMap()->LID(firstglobaldofid);

    int err  = velenr->ReplaceMyValue(firstlocaldofid,0,(*vel)[firstlocaldofid+4]);
        err += velenr->ReplaceMyValue(firstlocaldofid+1,0,(*vel)[firstlocaldofid+5]);
        err += velenr->ReplaceMyValue(firstlocaldofid+2,0,(*vel)[firstlocaldofid+6]);
        err += velenr->ReplaceMyValue(firstlocaldofid+3,0,(*vel)[firstlocaldofid+7]);
    if(err!=0)
      dserror("error during replacemyvalue");
  }
  return velenr;
}
