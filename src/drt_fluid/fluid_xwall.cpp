/*-----------------------------------------------------------*/
/*!

\brief implementation of enrichment-based wall model for turbulence

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_xwall.H"
#include "fluidimplicitintegration.H"
#include "fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_xwall.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_utils_rebalancing.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_fluid_turbulence/drt_transfer_turb_inflow.H"
#include "../drt_io/io.H"

#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 04/14 |
 *----------------------------------------------------------------------*/
FLD::XWall::XWall(Teuchos::RCP<DRT::Discretization> dis, int nsd,
    Teuchos::RCP<Teuchos::ParameterList>& params, Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
    Teuchos::RCP<FLD::UTILS::StressManager> wssmanager)
    : discret_(dis), params_(params), mystressmanager_(wssmanager), iter_(0)
{
  // get the processor ID from the communicator
  myrank_ = discret_->Comm().MyPID();

  if (myrank_ == 0)
  {
    std::cout << "\nWall modeling with a Spalding's law enrichment" << std::endl;
  }

  // some exclusions and safety checks:
  if (nsd != 3) dserror("Only 3D problems considered in xwall modelling!");
  if (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") !=
          INPAR::FLUID::timeint_afgenalpha &&
      DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo") !=
          INPAR::FLUID::timeint_npgenalpha)
    dserror(
        "Use Af-Genalpha for time integration in combination with xwall wall modeling. There would "
        "be additional updates necessary otherwise");

  if (params_->get<std::string>("predictor", "steady_state_predictor") != "steady_state")
    dserror("The meshtying framework does only support a steady-state predictor");

  std::string tauwtype = params_->sublist("WALL MODEL").get<std::string>("Tauw_Type", "constant");

  if (tauwtype == "constant")
    tauwtype_ = INPAR::FLUID::constant;
  else if (tauwtype == "between_steps")
    tauwtype_ = INPAR::FLUID::between_steps;
  else
    dserror("unknown Tauw_Type");

  std::string tauwcalctype =
      params_->sublist("WALL MODEL").get<std::string>("Tauw_Calc_Type", "residual");

  if (tauwcalctype == "residual")
    tauwcalctype_ = INPAR::FLUID::residual;
  else if (tauwcalctype == "gradient")
    tauwcalctype_ = INPAR::FLUID::gradient;
  else if (tauwcalctype == "gradient_to_residual")
    tauwcalctype_ = INPAR::FLUID::gradient_to_residual;
  else
    dserror("unknown Tauw_Calc_Type");

  constant_tauw_ = params_->sublist("WALL MODEL").get<double>("C_Tauw");

  min_tauw_ = params_->sublist("WALL MODEL").get<double>("Min_Tauw");

  fac_ = params_->sublist("WALL MODEL").get<double>("Inc_Tauw");

  // get gauss points
  gp_norm_ = params_->sublist("WALL MODEL").get<int>("GP_Wall_Normal");
  gp_norm_ow_ = params_->sublist("WALL MODEL").get<int>("GP_Wall_Normal_Off_Wall");
  gp_par_ = params_->sublist("WALL MODEL").get<int>("GP_Wall_Parallel");

  // compute initial pressure
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id == -1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
  dens_ = actmat->density_;
  visc_ = actmat->viscosity_ / dens_;  // here I want to have the kinematic viscosity


  std::string projectiontype = params_->sublist("WALL MODEL").get<std::string>("Projection", "No");

  if (projectiontype == "onlyl2projection")
    proj_ = true;
  else if (projectiontype == "No")
    proj_ = false;
  else
    dserror("unknown projection type");

  std::string blendingtype =
      params_->sublist("WALL MODEL").get<std::string>("Blending_Type", "none");

  if (blendingtype == "none")
    blendingtype_ = INPAR::FLUID::none;
  else if (blendingtype == "ramp_function")
    blendingtype_ = INPAR::FLUID::ramp_function;
  else
    dserror("unknown Blending_Type");

  inctauwnorm_ = 0.0;

  const int mlsmooth =
      (DRT::Problem::Instance()->FluidDynamicParams()).get<int>("WSS_ML_AGR_SOLVER");
  if (mlsmooth != -1)
    smooth_res_aggregation_ = true;
  else
    smooth_res_aggregation_ = false;

  switch_step_ = params_->sublist("WALL MODEL").get<int>("Switch_Step");
  if (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && switch_step_ < 2)
    dserror("provide reasonable Switch_Step if you want to use gradient_to_residual");

  if (smooth_res_aggregation_ && tauwcalctype_ == INPAR::FLUID::gradient)
    dserror("smoothing of tauw works only for residual-based tauw, as the residual is smoothed");

  fix_residual_on_inflow_ = DRT::INPUT::IntegralValue<int>(
      params_->sublist("WALL MODEL"), "Treat_Tauw_on_Dirichlet_Inflow");

  // output:
  if (myrank_ == 0)
  {
    std::cout << "\nXWall settings: " << std::endl;
    std::cout << "Tau_w is updated with:        " << tauwtype << std::endl;
    std::cout << "Tau_w is calculated with:     " << tauwcalctype << std::endl;
    std::cout << "Switching from grad to res:   " << switch_step_ << std::endl;
    std::cout << "Constant tau_w:               " << constant_tauw_ << std::endl;
    std::cout << "Minimum tau_w (clipping):     " << min_tauw_ << std::endl;
    std::cout << "Increment of tau_w:           " << fac_ << std::endl;
    std::cout << "Gauss rule:                   normal:  " << gp_norm_ << "  parallel:  " << gp_par_
              << "  overall:  " << gp_norm_ * gp_par_ * gp_par_ << std::endl;
    std::cout << "Gauss rule (off-wall):        normal:  " << gp_norm_ow_
              << "  parallel:  " << gp_par_ << "  overall:  " << gp_norm_ow_ * gp_par_ * gp_par_
              << std::endl;
    std::cout << "Enriched DOFs l2-projected:   " << projectiontype << std::endl;
    std::cout << "Blending method:              " << blendingtype << std::endl;
    std::cout << "Smooth tau_w:                 " << smooth_res_aggregation_ << std::endl;
    std::cout << "Solver for tau_w smoothing:   "
              << (DRT::Problem::Instance()->FluidDynamicParams()).get<int>("WSS_ML_AGR_SOLVER")
              << std::endl;
    std::cout << "Solver for projection:        "
              << params_->sublist("WALL MODEL").get<int>("PROJECTION_SOLVER") << std::endl;
    std::cout << "Fix Dirichlet inflow:         " << fix_residual_on_inflow_ << std::endl;
    std::cout << std::endl;
    std::cout << "WARNING: ramp functions are used to treat fluid Mortar coupling conditions"
              << std::endl;
    std::cout << "WARNING: face element with enrichment not implemented" << std::endl;
    std::cout << "WARNING: stabilization terms for Neumann inflow only on standard space"
              << std::endl;
  }

  Setup();

  turbulent_inflow_condition_ =
      Teuchos::rcp(new TransferTurbulentInflowConditionNodal(discret_, dbcmaps));

  if (turbulent_inflow_condition_->IsActive())
  {
    oldtauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    oldinctauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
  }
  else
  {
    oldtauw_ = Teuchos::null;
    oldinctauw_ = Teuchos::null;
  }
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetXWallParams(Teuchos::ParameterList& eleparams)
{
  // params required for the shape functions
  eleparams.set("walldist", wdist_);
  eleparams.set("tauw", tauw_);
  eleparams.set("inctauw", inctauw_);
  eleparams.set("xwalltoggle", xwalltoggle_);

  eleparams.set("mk", mkstate_);

  eleparams.set("gpnorm", gp_norm_);
  eleparams.set("gpnormow", gp_norm_ow_);
  eleparams.set("gppar", gp_par_);

  return;
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 08/14 |
 |  Used for the xwdiscret_, which is redistributed                     |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetXWallParamsXWDis(Teuchos::ParameterList& eleparams)
{
  // params required for the shape functions
  eleparams.set("walldist", wdistxwdis_);
  eleparams.set("tauw", tauwxwdis_);
  eleparams.set("inctauw", inctauwxwdis_);
  eleparams.set("xwalltoggle", xwalltogglexwdis_);

  eleparams.set("mk", mkxwstate_);

  eleparams.set("gpnorm", gp_norm_);
  eleparams.set("gpnormow", gp_norm_ow_);
  eleparams.set("gppar", gp_par_);
  return;
}

/*----------------------------------------------------------------------*
 |  Setup XWall                                                bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::Setup()
{
  if (myrank_ == 0) std::cout << "Setup: " << std::endl;

  InitXWallMaps();

  InitWallDist();

  SetupXWallDis();

  tauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()), true));
  inctauw_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()), true));

  // initialize for first call or for constant tauw setting
  tauw_->PutScalar(constant_tauw_);

  wdistxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()), true));
  LINALG::Export(*walldist_, *wdistxwdis_);

  tauwxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()), true));
  inctauwxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()), true));

  // initialize for first call or for constant tauw setting
  tauwxwdis_->PutScalar(constant_tauw_);

  InitToggleVector();

  SetupL2Projection();

  {
    // the value for linear elements is 1/3
    // initialize just in case
    mkxwstate_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->ElementColMap()), true));
    mkxwstate_->PutScalar(0.33333333333);
    mkstate_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementColMap()), true));
    mkstate_->PutScalar(0.33333333333);

    restart_wss_ = Teuchos::null;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Setup basic xwall map and dirichlet map                    bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitXWallMaps()
{
  if (myrank_ == 0) std::cout << "- build xwall maps...                                   ";

  // build row vector of all xwall nodes
  {
    std::vector<int> rowvec;  // node row map
    for (int i = 0; i < discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      int xwallgid = discret_->NodeRowMap()->GID(i);
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if (!xwallnode) dserror("Cannot find node");

      bool enriched = false;

      // check if one of the surrounding elements is xwall element
      DRT::Element** surrele = xwallnode->Elements();
      for (int k = 0; k < xwallnode->NumElement(); ++k)
      {
        DRT::ELEMENTS::FluidXWall* xwallele = dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

        if (xwallele) enriched = true;
      }

      if (enriched) rowvec.push_back(xwallgid);
    }

    xwallrownodemap_ =
        Teuchos::rcp(new Epetra_Map(-1, (int)rowvec.size(), &rowvec[0], 0, discret_->Comm()));
  }

  // get Dirichlet conditions
  std::vector<DRT::Condition*> dircond;
  discret_->GetCondition("FluidStressCalc", dircond);

  if (not dircond.empty())
  {
    std::vector<int> testcollect;
    int count = 0;
    for (unsigned numcond = 0; numcond < dircond.size(); ++numcond)
    {
      const std::vector<int>* test = dircond[numcond]->Nodes();
      int j = 0;
      for (std::vector<int>::const_iterator i = (*test).begin(); i != (*test).end(); ++i)
      {
        ++count;
        testcollect.push_back((*test)[j]);
        ++j;
      }
    }

    int gcount;
    (discret_->Comm()).SumAll(&count, &gcount, 1);
    dircolnodemap_ =
        Teuchos::rcp(new Epetra_Map(gcount, count, &testcollect[0], 0, discret_->Comm()));
  }  // end loop this conditions
  else
    dserror("You need DESIGN FLUID STRESS CALC SURF CONDITIONS for xwall");


  // map is of course not unique as it is a column map
  //  if(dircolnodemap_->UniqueGIDs())
  //    dserror("Map resulting from DESIGN FLUID STRESS CALC SURF CONDITIONS not unique, probably
  //    node specified on two conditions?");

  if (myrank_ == 0)
    std::cout << xwallrownodemap_->NumGlobalElements() << " XWall nodes initialized!" << std::endl;
  if (xwallrownodemap_->NumGlobalElements() == 0) dserror("No XWall elements found");
  return;
}


/*----------------------------------------------------------------------*
 |  Calculate wall distance and coupling matrix of for tauw    bk 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitWallDist()
{
  if (myrank_ == 0) std::cout << "- calculate wall distance...                            ";

  tauwcouplingmattrans_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*xwallrownodemap_, 2, false, false));
  wdist_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()), true));
  walldist_ = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));

  // build a new discretization which lies on all procs
  Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(discret_->Comm().Clone());

  // this is very expensive in terms of memory
  // we will delete it as soon as we are ready here
  Teuchos::RCP<DRT::Discretization> commondis =
      Teuchos::rcp(new DRT::Discretization((std::string) "Commondis", newcomm));

  // loop over all column nodes of underlying problem discret and add
  for (int i = 0; i < (discret_->NodeColMap())->NumMyElements(); ++i)
  {
    DRT::Node* node = discret_->lColNode(i);
    if (!node) dserror("Cannot find node with lid %", i);
    Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
    commondis->AddNode(newnode);
  }
  // loop over all column elements of underlying problem discret and add
  for (int i = 0; i < (discret_->ElementColMap())->NumMyElements(); ++i)
  {
    DRT::Element* node = discret_->lColElement(i);
    if (!node) dserror("Cannot find ele with lid %", i);
    Teuchos::RCP<DRT::Element> newnode = Teuchos::rcp(node->Clone());
    commondis->AddElement(newnode);
  }

  Teuchos::RCP<Epetra_Map> testrednodecolmap = LINALG::AllreduceEMap(*(discret_->NodeRowMap()));
  commondis->ExportColumnNodes(*testrednodecolmap);

  // do not assign any dofs to save memory
  // only the nodes are needed in this discretization
  commondis->FillComplete(false, false, false);

  // also build a fully overlapping map of enriched nodes here:
  std::vector<int> colvec;  // node col map
  for (int i = 0; i < (commondis->NodeColMap())->NumMyElements(); ++i)
  {
    int gid = commondis->NodeColMap()->GID(i);
    DRT::Node* xwallnode = commondis->lColNode(i);
    if (!xwallnode) dserror("Cannot find node with lid %", i);
    int enriched = 0;

    DRT::Element** surrele = xwallnode->Elements();
    for (int k = 0; k < xwallnode->NumElement(); ++k)
    {
      DRT::ELEMENTS::FluidXWall* xwallele = dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

      if (xwallele) enriched = 1;
    }
    int genriched = 0;
    (commondis->Comm()).SumAll(&enriched, &genriched, 1);
    if (genriched > 0) colvec.push_back(gid);
  }
  int count = (int)colvec.size();

  xwallcolnodemap_ = Teuchos::rcp(new Epetra_Map(count, count, &colvec[0], 0, discret_->Comm()));

  for (int j = 0; j < xwallcolnodemap_->NumMyElements(); ++j)
  {
    int xwallgid = xwallcolnodemap_->GID(j);

    DRT::Node* xwallnode = commondis->gNode(xwallgid);
    if (!xwallnode) dserror("Cannot find node");

    double mydist = 1.0E10;
    double gdist = 1.0E10;
    int mygid = 0;


    for (int i = 0; i < dircolnodemap_->NumMyElements(); ++i)
    {
      int gid = dircolnodemap_->GID(i);

      if (discret_->NodeRowMap()->MyGID(gid))
      {
        DRT::Node* node = discret_->gNode(gid);

        if (!node) dserror("ERROR: Cannot find wall node with gid %", gid);

        double newdist =
            sqrt(((xwallnode->X())[0] - (node->X())[0]) * ((xwallnode->X())[0] - (node->X())[0]) +
                 ((xwallnode->X())[1] - (node->X())[1]) * ((xwallnode->X())[1] - (node->X())[1]) +
                 ((xwallnode->X())[2] - (node->X())[2]) * ((xwallnode->X())[2] - (node->X())[2]));
        if (newdist < mydist)
        {
          mydist = newdist;
          mygid = gid;
        }
      }
    }

    discret_->Comm().MinAll(&mydist, &gdist, 1);

    // now write this value in the node based vector
    if (xwallrownodemap_->MyGID(xwallgid))
    {
      int err = walldist_->ReplaceGlobalValues(1, &gdist, &xwallgid);
      if (err > 0)
        dserror("global row not on proc");
      else if (err < 0)
        dserror("wrong vector index");
    }

    // this is the processor that knows the respective node
    if (mydist == gdist)
    {
      tauwcouplingmattrans_->Assemble(1.0, mygid, xwallgid);
    }
  }

  LINALG::Export(*walldist_, *wdist_);
  tauwcouplingmattrans_->Complete();

  double mean = 0.0;
  walldist_->MeanValue(&mean);

  if (myrank_ == 0)
    std::cout << "the mean distance from the wall of all XWall nodes is: " << mean << "... ";

  if (myrank_ == 0) std::cout << "done!  " << std::endl;

  commondis = Teuchos::null;
  return;
}

/*----------------------------------------------------------------------*
 |  Build a node-based toggle vector (on/off=1.0/0.0, 0.7)     bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::InitToggleVector()
{
  if (myrank_ == 0) std::cout << "- build enriched/blending toggle vector for elements... ";

  xwalltoggle_ = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()), true));
  xwalltogglexwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()), true));
  xtoggleloc_ = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  int count = 0;
  for (int j = 0; j < xwallrownodemap_->NumMyElements(); ++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (discret_->NodeRowMap()->MyGID(xwallgid))  // just in case
    {
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if (!xwallnode) dserror("Cannot find node");

      bool fullyenriched = true;

      // Mortar interface
      std::vector<DRT::Condition*> mortarcond;
      xwallnode->GetCondition("Mortar", mortarcond);
      if (not mortarcond.empty()) fullyenriched = false;

      // get all surrounding elements
      DRT::Element** surrele = xwallnode->Elements();
      for (int k = 0; k < xwallnode->NumElement(); ++k)
      {
        DRT::ELEMENTS::FluidXWall* xwallele = dynamic_cast<DRT::ELEMENTS::FluidXWall*>(surrele[k]);

        if (!xwallele) fullyenriched = false;
      }

      if (fullyenriched == true)
      {
        int err = xtoggleloc_->ReplaceMyValue(j, 0, 1.0);
        if (err != 0) dserror("something went wrong");
      }
      else
      {
        if (blendingtype_ != INPAR::FLUID::ramp_function)
        {
          int err = xtoggleloc_->ReplaceMyValue(j, 0, 0.7);
          if (err != 0) dserror("something went wrong");
        }
        count++;
      }
    }
  }

  LINALG::Export(*xtoggleloc_, *xwalltoggle_);
  LINALG::Export(*xtoggleloc_, *xwalltogglexwdis_);

  int gcount;
  (discret_->Comm()).SumAll(&count, &gcount, 1);
  if (myrank_ == 0) std::cout << gcount << " blending nodes identified... ";

  if (myrank_ == 0) std::cout << "done!  " << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Build a discretization only including xwall nodes/elements bk 07/14 |
 |  and redistribute                                                    |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetupXWallDis()
{
  // build a new discretization
  Teuchos::RCP<Epetra_Comm> newcomm = Teuchos::rcp(discret_->Comm().Clone());

  xwdiscret_ = Teuchos::rcp(new DRT::Discretization((std::string) "xwalldis", newcomm));

  // loop over all xwall row nodes and add
  for (int i = 0; i < (discret_->NodeColMap())->NumMyElements(); ++i)
  {
    int gid = (discret_->NodeColMap())->GID(i);

    if (xwallcolnodemap_->MyGID(gid))
    {
      DRT::Node* node = discret_->lColNode(i);
      if (!node) dserror("Cannot find node with lid %", i);
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(node->Clone());
      xwdiscret_->AddNode(newnode);
    }
  }
  // loop over all column elements of underlying problem discret and add
  for (int i = 0; i < discret_->NumMyColElements(); ++i)
  {
    DRT::Element* ele = discret_->lColElement(i);
    if (!ele) dserror("Cannot find ele with lid %", i);

    DRT::ELEMENTS::FluidXWall* xwallele = dynamic_cast<DRT::ELEMENTS::FluidXWall*>(ele);

    if (!xwallele) continue;

    Teuchos::RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    xwdiscret_->AddElement(newele);
  }

  // make all conditions known to the child discretization
  // i.e. periodic boundary conditions, dirichlet conditions, ...
  {
    // get all conditions types prescribed in the input file
    std::vector<std::string> allcond;
    discret_->GetConditionNames(allcond);
    // loop all conditions types
    for (unsigned numcond = 0; numcond < allcond.size(); ++numcond)
    {
      // get condition
      std::vector<DRT::Condition*> actcond;
      discret_->GetCondition(allcond[numcond], actcond);
      // loop all condition of the current type
      for (unsigned numactcond = 0; numactcond < actcond.size(); ++numactcond)
      {
        // finally set condition
        xwdiscret_->SetCondition(
            allcond[numcond], Teuchos::rcp(new DRT::Condition(*actcond[numactcond])));
      }
    }
  }


  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (xwdiscret_->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new discretization
  Teuchos::RCP<DRT::DofSet> newdofset =
      Teuchos::rcp(new DRT::TransparentDofSet(discret_, parallel));

  xwdiscret_->ReplaceDofSet(newdofset);
  xwdiscret_->FillComplete(true, true, true);

  // redistribute and treat periodic bc if parallel
  if (parallel)
  {
    // redistribute
#if defined(PARALLEL) && defined(HAVE_PARMETIS)

    Teuchos::RCP<Epetra_Map> elemap = Teuchos::rcp(new Epetra_Map(*xwdiscret_->ElementRowMap()));
    Teuchos::RCP<Epetra_Map> rownodes;
    Teuchos::RCP<Epetra_Map> colnodes;
    Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(discret_->Comm().Clone());
    DRT::UTILS::ComputeRebalancedNodeMaps(
        xwdiscret_, elemap, rownodes, colnodes, comm, true, comm->NumProc());
    // rebuild of the system with new maps
    xwdiscret_->Redistribute(*rownodes, *colnodes, false, false);

#else
#if defined(PARALLEL)
    dserror("require PARMETIS not METIS");
#endif
#endif

    PeriodicBoundaryConditions pbc(xwdiscret_, false);
    pbc.UpdateDofsForPeriodicBoundaryConditions();
    xwdiscret_->ReplaceDofSet(newdofset);
    xwdiscret_->FillComplete(true, true, true);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Setup matrix, vectors and solver for projection            bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::SetupL2Projection()
{
  // create matrix for projection
  if (proj_)
  {
    // build the dof maps of p and u_y

    std::vector<int> enrdf;  // enriched dofs

    for (int i = 0; i < xwdiscret_->NodeRowMap()->NumMyElements(); ++i)
    {
      int gid = xwdiscret_->NodeRowMap()->GID(i);
      // continue only if on this proc
      if (xwdiscret_->NodeRowMap()->MyGID(gid))
      {
        DRT::Node* node = xwdiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find off wall node with gid %", gid);
        // make sure that periodic nodes are not assembled twice
        std::vector<DRT::Condition*> periodiccond;
        node->GetCondition("SurfacePeriodic", periodiccond);
        // make sure that slave periodic bc are not included
        bool includedofs = true;
        if (not periodiccond.empty())
        {
          for (unsigned numcondper = 0; numcondper < periodiccond.size(); ++numcondper)
          {
            const std::string* mymasterslavetoggle =
                periodiccond[numcondper]->Get<std::string>("Is slave periodic boundary condition");
            if (*mymasterslavetoggle == "Slave")
            {
              includedofs = false;
            }
          }
        }
        if (includedofs)
        {
          int firstglobaldofid = xwdiscret_->Dof(0, node, 0);
          enrdf.push_back(firstglobaldofid + 4);
          enrdf.push_back(firstglobaldofid + 5);
          enrdf.push_back(firstglobaldofid + 6);
        }
      }
    }

    enrdofrowmap_ =
        Teuchos::rcp(new Epetra_Map(-1, (int)enrdf.size(), &enrdf[0], 0, xwdiscret_->Comm()));

    massmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*enrdofrowmap_, 108, false, true));

    incveln_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
    incvelnp_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
    incaccn_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));

    stateveln_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()), true));
    statevelnp_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()), true));
    stateaccn_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()), true));

    mergedmap_ = Teuchos::null;
    lagrdofrowmap_ = Teuchos::null;
  }
  else
  {
    massmatrix_ = Teuchos::null;
    incveln_ = Teuchos::null;
    incvelnp_ = Teuchos::null;
    incaccn_ = Teuchos::null;
    stateveln_ = Teuchos::null;
    statevelnp_ = Teuchos::null;
    stateaccn_ = Teuchos::null;
    enrdofrowmap_ = Teuchos::null;
    lagrdofrowmap_ = Teuchos::null;
  }

  // setup solver
  if (proj_)
  {
    const int solvernumber = params_->sublist("WALL MODEL").get<int>("PROJECTION_SOLVER");
    if (solvernumber < 0) dserror("provide a solver number for l2-projection");
    // get solver parameter list of linear solver
    const Teuchos::ParameterList& solverparams =
        DRT::Problem::Instance()->SolverParams(solvernumber);
    const int solvertype =
        DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

    solver_ = Teuchos::rcp(new LINALG::Solver(
        solverparams, xwdiscret_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

    if (solvertype != INPAR::SOLVER::umfpack)
    {
      if (solvertype != INPAR::SOLVER::belos && myrank_ == 0)
        std::cout
            << "\nUse Belos as solver because it can handle several right hand sides at once!\n"
            << std::endl;
      const int prectyp =
          DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams, "AZPREC");
      // watch out: only ILU might work right now because of compute nullspace might not work...?
      // ... test!
      switch (prectyp)
      {
        case INPAR::SOLVER::azprec_ML:
        case INPAR::SOLVER::azprec_MLfluid:
        case INPAR::SOLVER::azprec_MLAPI:
        case INPAR::SOLVER::azprec_MLfluid2:
        case INPAR::SOLVER::azprec_MueLuAMG_sym:
        case INPAR::SOLVER::azprec_MueLuAMG_nonsym:
        {
          if (proj_)
          {  // has 3 dofs, velocity dofs
            // BUT: enriched nodes have 8 dofs, so we have to calculate our own nullspace for 3 dofs
            // store nv and np at unique location in solver parameter list
            solver_->Params().sublist("NodalBlockInformation").set("nv", 3);
            solver_->Params().sublist("NodalBlockInformation").set("np", 0);
            solver_->Params().sublist("NodalBlockInformation").set("numdf", 3);
            solver_->Params().sublist("NodalBlockInformation").set("dimns", 3);
            // get the belos list and see whether we use downwinding
            Teuchos::ParameterList& beloslist = solver_->Params().sublist("Belos Parameters");

            beloslist.set<int>("downwinding nv", 3);
            beloslist.set<int>("downwinding np", 0);

            Teuchos::ParameterList* mllist_ptr = NULL;
            mllist_ptr = &((solver_->Params()).sublist("ML Parameters"));
            Teuchos::ParameterList& mllist = *mllist_ptr;  // solveparams.sublist("ML Parameters");

            Teuchos::RCP<std::vector<double>> ns =
                mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
            if (ns != Teuchos::null) break;

            // no, we have not previously computed the nullspace
            // -> compute nullspace
            ns = Teuchos::null;

            mllist.set<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
            // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
            mllist.set<bool>("ML validate parameter list", false);

            mllist.set("PDE equations", 3);
            mllist.set("null space: dimension", 3);
            mllist.set("null space: type", "pre-computed");
            mllist.set("null space: add default vectors", false);

            // allocate dimns times the local length of the rowmap
            const int lrows = enrdofrowmap_->NumMyElements();
            ns = Teuchos::rcp(new std::vector<double>(3 * lrows));
            double* nullsp = &((*ns)[0]);
            mllist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
            mllist.set("null space: vectors", nullsp);

            // now compute the nullspace
            double* mode[6];
            for (int i = 0; i < 3; ++i) mode[i] = &((*ns)[i * lrows]);

            for (int i = 0; i < xwdiscret_->NumMyRowNodes(); ++i)
            {
              const unsigned int ndof = 3;
              DRT::Node* actnode = xwdiscret_->lRowNode(i);
              if (!actnode) dserror("cannot find node");
              std::vector<int> dofs = xwdiscret_->Dof(0, actnode);
              std::vector<int> actdofs;
              // only dof 4...6 (enriched dofs)
              actdofs.push_back(dofs.at(4));
              actdofs.push_back(dofs.at(5));
              actdofs.push_back(dofs.at(6));

              for (unsigned j = 0; j < ndof; ++j)
              {
                const int dof = actdofs.at(j);

                const int lid = enrdofrowmap_->LID(dof);
                if (lid < 0) dserror("Cannot find dof %i", dof);

                for (unsigned k = 0; k < ndof; ++k)
                {
                  if (k % ndof == j % ndof)
                    mode[k % ndof][lid] = 1.0;
                  else
                    mode[k % ndof][lid] = 0.0;
                }
              }  // for (int j=0; j<actnode->Dof().NumDof(); ++j)
            }    // for (int i=0; i<NumMyRowNodes(); ++i)
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
void FLD::XWall::UpdateTauW(int step, Teuchos::RCP<Epetra_Vector> trueresidual, int itnum,
    Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> veln)
{
  iter_ = 0;

  TransferAndSaveTauw();

  Teuchos::RCP<Epetra_Vector> newtauw = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  // calculate wall stress
  Teuchos::RCP<Epetra_Vector> wss;
  if (restart_wss_ != Teuchos::null)
  {
    wss = restart_wss_;
  }
  else if (((tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step >= switch_step_) ||
               (tauwcalctype_ != INPAR::FLUID::gradient_to_residual && step > 1)))
  {
    if (mystressmanager_ == Teuchos::null) dserror("wssmanager not available in xwall");
    // fix nodal forces on dirichlet inflow surfaces if desired
    wss = mystressmanager_->GetPreCalcWallShearStresses(FixDirichletInflow(trueresidual));
  }
  switch (tauwtype_)
  {
    case INPAR::FLUID::constant:  // works
    {
      if (step == 1) CalcMK();
      // tauw_ is constant and inctauw_ is zero
      return;
    }
    break;
    case INPAR::FLUID::between_steps:
    {
      inctauw_->PutScalar(0.0);

      if (itnum == 0)  // in between steps
        CalcTauW(step, velnp, wss);
      else
        return;
    }
    break;
    default:
      dserror("unknown tauwtype_");
      break;
  }

  LINALG::Export(*tauw_, *newtauw);

  double actmean = -1.0;
  newtauw->MeanValue(&actmean);

  double min = -1.0;
  newtauw->MinValue(&min);
  if (min < 1e-10) dserror("tauw is zero");
  double max = -1.0;
  newtauw->MaxValue(&max);

  // convergence check, works only if we don't take the mean
  newtauw->PutScalar(0.0);
  LINALG::Export(*inctauw_, *newtauw);
  newtauw->Norm2(&inctauwnorm_);
  // rescale inctauw to full increment
  if (fac_ > 1.0e-8) inctauwnorm_ /= fac_;

  if (myrank_ == 0)
    std::cout << "  min:  " << min << "  max:  " << max << "  mean-applied:  " << actmean
              << "  inc norm2:  " << inctauwnorm_;

  // also project the other vectors
  // later I could also run all three at the same time (with the same matrix)

  if (proj_)
  {
    if (myrank_ == 0) std::cout << "  L2-project... ";
    if (tauwtype_ == INPAR::FLUID::between_steps)
    {
      L2ProjectVector(veln, Teuchos::null, accn);

      // at the beginning of this time step they are equal -> calculate only one of them
      velnp->Update(1.0, *veln, 0.0);
    }
    else
      L2ProjectVector(veln, velnp, accn);

    if (myrank_ == 0) std::cout << "done!" << std::endl;
  }
  else
    std::cout << std::endl;

  CalcMK();

  // destruct vector so that we don't use it next time
  if (restart_wss_ != Teuchos::null) restart_wss_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  Routines to calculate Tauw                                 bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::CalcTauW(
    int step, Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> wss)
{
  Teuchos::RCP<Epetra_Vector> newtauw = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> newtauw2 = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> tauw =
      Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap()), true));

  if (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && switch_step_ == step && myrank_ == 0)
    std::cout << "\n switching from gradient to residual \n" << std::endl;


  if (tauwcalctype_ == INPAR::FLUID::residual ||
      (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step >= switch_step_))
  {
    for (int lnodeid = 0; lnodeid < dircolnodemap_->NumMyElements(); lnodeid++)
    {
      int gid = dircolnodemap_->GID(lnodeid);
      // continue only if on this proc
      if (discret_->NodeRowMap()->MyGID(gid))
      {
        DRT::Node* node = discret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find off wall node with gid %", gid);

        int firstglobaldofid = discret_->Dof(0, node, 0);
        int firstlocaldofid = wss->Map().LID(firstglobaldofid);

        if (firstlocaldofid < 0) dserror("localdofid not found in map for given globaldofid");
        double forcex = (*wss)[firstlocaldofid];
        double forcey = (*wss)[firstlocaldofid + 1];
        double forcez = (*wss)[firstlocaldofid + 2];

        double tauw = sqrt(forcex * forcex + forcey * forcey + forcez * forcez);

        // this is necessary since we are dividing by tauw on element level
        // also, the shape functions become singular, if tauw==0
        if (tauw < min_tauw_) tauw = min_tauw_;
        // store in vector
        int err = newtauw->ReplaceGlobalValue(gid, 0, tauw);
        if (err != 0) dserror("something went wrong during replacemyvalue");
      }
    }
  }
  else if (tauwcalctype_ == INPAR::FLUID::gradient ||
           (tauwcalctype_ == INPAR::FLUID::gradient_to_residual && step < switch_step_))
  {
    // necessary to set right state (the maps of the state vector and discretization have to be
    // equal)
    Teuchos::RCP<Epetra_Vector> statevel =
        Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->DofRowMap()), true));
    Teuchos::RCP<Epetra_Vector> newtauwxwdis =
        Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeRowMap()), true));
    LINALG::Export(*velnp, *statevel);

    xwdiscret_->SetState("vel", statevel);

    Teuchos::RCP<Epetra_Vector> timesvec =
        Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeRowMap()), true));

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

    // set action in order to project element void fraction to nodal void fraction
    Teuchos::ParameterList params;

    SetXWallParamsXWDis(params);

    params.set<int>("action", FLD::tauw_via_gradient);

    // loop column elements: vector
    for (int i = 0; i < numele; ++i)
    {
      DRT::Element* actele = xwdiscret_->lColElement(i);

      const int numnode = actele->NumNode();

      // get element location vector and ownerships
      actele->LocationVector(*xwdiscret_, lm, lmowner, lmstride);

      elevector1.Size(numnode);
      elevector2.Size(numnode);

      // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
      // elevector1 has to be NULL here, because I am assuming a dof-based vector otherwise
      actele->Evaluate(
          params, *xwdiscret_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);

      // get element location vector for nodes
      lm.resize(numnode);
      lmowner.resize(numnode);

      DRT::Node** nodes = actele->Nodes();
      for (int n = 0; n < numnode; ++n)
      {
        lm[n] = nodes[n]->Id();
        lmowner[n] = nodes[n]->Owner();
      }

      // assembling into node maps
      LINALG::Assemble(*newtauwxwdis, elevector1, lm, lmowner);
      LINALG::Assemble(*timesvec, elevector2, lm, lmowner);
    }  // end element loop

    xwdiscret_->ClearState();

    // scale with times:
    for (int l = 0; l < (xwdiscret_->NodeRowMap())->NumMyElements(); l++)
    {
      double sumnewtauw = (*newtauwxwdis)[l];
      double timesfac = (*timesvec)[l];
      double newtauwsc = 1.0;
      // we only have to do something, if we assembled at least once the same value
      if (timesfac > 0.5)
      {
        newtauwsc = sumnewtauw / timesfac;

        if (newtauwsc < min_tauw_) newtauwsc = min_tauw_;
        int err = newtauwxwdis->ReplaceMyValue(l, 0, newtauwsc);
        if (err != 0) dserror("something went wrong during replacemyvalue");
      }
    }
    LINALG::Export(*newtauwxwdis, *newtauw);
  }
  else
    dserror("unknown tauwcalctype_");

  tauw->Update(1.0, *tauw_, 0.0);


  inctauw_->Update(1.0, *tauw, 0.0);
  tauw->PutScalar(0.0);

  tauwcouplingmattrans_->Multiply(true, *newtauw, *newtauw2);
  double meansp = 0.0;
  newtauw2->MeanValue(&meansp);

  LINALG::Export(*newtauw2, *tauw);
  inctauw_->Update(fac_, *tauw, -fac_);  // now this is the increment (new-old)

  tauw_->Update(1.0, *inctauw_, 1.0);

  OverwriteTransferredValues();

  LINALG::Export(*inctauw_, *newtauw2);
  LINALG::Export(*newtauw2, *inctauwxwdis_);
  LINALG::Export(*tauw_, *newtauw2);
  LINALG::Export(*newtauw2, *tauwxwdis_);

  if (meansp < 2.0e-9)
    dserror(
        "Average wall shear stress is zero. You probably forgot to specify approprite DESIGN FLUID "
        "STRESS CALC SURF CONDITIONS where the stress should be calculated.");

  if (myrank_ == 0) std::cout << "tauw mean:  " << meansp;

  return;
}

/*----------------------------------------------------------------------*
 |  L2-project enriched dofs of vector                         bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::L2ProjectVector(Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> accn)
{
  if (not veln->Map().SameAs(*discret_->DofRowMap()))
    dserror("input map is not the dof row map of the fluid discretization");

  massmatrix_->Zero();

  incveln_->PutScalar(0.0);
  if (accn != Teuchos::null) incaccn_->PutScalar(0.0);
  if (velnp != Teuchos::null) incvelnp_->PutScalar(0.0);

  LINALG::Export(*veln, *stateveln_);
  if (accn != Teuchos::null) LINALG::Export(*accn, *stateaccn_);
  if (velnp != Teuchos::null) LINALG::Export(*velnp, *statevelnp_);

  // number of right hand sides during solving
  // is the number of velocity components that is solved for
  // 3 since we are in 3D
  int numberofrhs = 0;
  if (velnp == Teuchos::null && accn == Teuchos::null)
    numberofrhs = 1;
  else if (velnp == Teuchos::null || accn == Teuchos::null)
    numberofrhs = 2;
  else
    numberofrhs = 3;

  xwdiscret_->SetState("veln", stateveln_);
  if (accn != Teuchos::null) xwdiscret_->SetState("accn", stateaccn_);
  if (velnp != Teuchos::null) xwdiscret_->SetState("velnp", statevelnp_);

  // set action in order to project nodal enriched values to new shape functions
  Teuchos::ParameterList params;
  SetXWallParamsXWDis(params);
  params.set<int>("action", FLD::xwall_l2_projection);

  // create empty right hand side
  Teuchos::RCP<Epetra_MultiVector> rhsassemble =
      Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_, numberofrhs));

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
  for (int i = 0; i < numele; ++i)
  {
    DRT::Element* actele = xwdiscret_->lColElement(i);

    const int numnode = actele->NumNode();
    const int numdf = 3;

    // get element location vector and ownerships
    actele->LocationVector(*xwdiscret_, lm, lmowner, lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elematrix1.Shape(numnode * numdf, numnode * numdf);
    // Reshape element matrices and vectors and initialize to zero
    elematrix2.Shape(
        numnode * numdf, numberofrhs);  // we have 3 right hand sides for now: 3 velocity components

    // call the element specific evaluate method (elemat1 = mass matrix, elemat2 = rhs)
    actele->Evaluate(
        params, *xwdiscret_, lm, elematrix1, elematrix2, elevectordummy, elevector2, elevector3);

    // get element location vector for enriched dofs
    std::vector<int> lmassemble;
    std::vector<int> lmownerassemble;
    lmassemble.resize(numnode * numdf);
    lmownerassemble.resize(numnode * numdf);

    for (int n = 0; n < numnode; ++n)
    {
      for (int df = 4; df < 7; ++df)
      {
        lmassemble[n * numdf + df - 4] = lm[n * 8 + df];
        lmownerassemble[n * numdf + df - 4] = lmowner[n * 8 + df];
      }
    }

    // assembling into node maps
    massmatrix_->Assemble(actele->Id(), elematrix1, lmassemble, lmownerassemble);

    // assembling into node maps
    // assemble numberofrhs entries in rhs vector sequentially
    elevector1.Size(numnode * numdf);
    for (int n = 0; n < numberofrhs; ++n)
    {
      // copy results into Serial_DenseVector for assembling
      for (int idf = 0; idf < numnode * numdf; ++idf) elevector1(idf) = elematrix2(idf, n);
      // assemble into nth vector of MultiVector
      LINALG::Assemble(*rhsassemble, n, elevector1, lmassemble, lmownerassemble);
    }
  }  // end element loop

  xwdiscret_->ClearState();
  // finalize the matrix
  massmatrix_->Complete();

  // solution vector
  Teuchos::RCP<Epetra_MultiVector> resultvec =
      Teuchos::rcp(new Epetra_MultiVector(*enrdofrowmap_, numberofrhs));

  // solve for 1, 2 or 3 right hand sides at the same time --> thanks to Belos
  solver_->Solve(massmatrix_->EpetraOperator(), resultvec, rhsassemble, true, true, Teuchos::null);

  // now copy result in original vector: the result is an increment of the velocity/ acceleration
  LINALG::Export(*((*resultvec)(0)), *incveln_);
  if (numberofrhs > 1) LINALG::Export(*((*resultvec)(1)), *incaccn_);
  if (numberofrhs > 2) LINALG::Export(*((*resultvec)(2)), *incvelnp_);

  veln->Update(1.0, *incveln_, 1.0);
  if (accn != Teuchos::null) accn->Update(1.0, *incaccn_, 1.0);
  if (velnp != Teuchos::null) velnp->Update(1.0, *incvelnp_, 1.0);


  return;
}

/*----------------------------------------------------------------------*
 |  Adapt ML Nullspace for MFS aggregation                     bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::AdaptMLNullspace(const Teuchos::RCP<LINALG::Solver>& solver)
{
  // extract the ML parameters:
  Teuchos::ParameterList& mlparams = solver->Params().sublist("ML Parameters");

  // get nullspace parameters
  double* nullspace = mlparams.get("null space: vectors", (double*)NULL);
  if (!nullspace) dserror("No nullspace supplied in parameter list");
  int nsdim = mlparams.get("null space: dimension", 1);
  if (nsdim != 4) dserror("Wrong Nullspace dimension for XWall");
  int lrowdofs = discret_->DofRowMap()->NumMyElements();
  //  std::cout << "lrowdofs  " << lrowdofs << std::endl;
  // std::cout << "check the nullspace for mfs" << std::endl;
  for (int j = 0; j < xwallrownodemap_->NumMyElements(); ++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (not discret_->NodeRowMap()->MyGID(xwallgid))  // just in case
      dserror("not on proc");
    {
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if (!xwallnode) dserror("Cannot find node");

      int firstglobaldofid = discret_->Dof(xwallnode, 0);
      int firstlocaldofid = discret_->DofRowMap()->LID(firstglobaldofid);

      nullspace[firstlocaldofid + 4] = 0.0;
      nullspace[lrowdofs + firstlocaldofid + 5] = 0.0;
      nullspace[lrowdofs * 2 + firstlocaldofid + 6] = 0.0;
      nullspace[lrowdofs * 3 + firstlocaldofid + 7] = 0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Calculate MK for residual-based stabilization parameter    bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::CalcMK()
{
  Teuchos::RCP<Epetra_MultiVector> mkxw =
      Teuchos::rcp(new Epetra_MultiVector(*(xwdiscret_->ElementRowMap()), 1, true));

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", FLD::xwall_calc_mk);

  SetXWallParamsXWDis(eleparams);

  xwdiscret_->EvaluateScalars(eleparams, mkxw);


  Teuchos::RCP<Epetra_Vector> mkxwv =
      Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->ElementRowMap()), true));
  Teuchos::RCP<Epetra_Vector> mkv =
      Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()), true));

  // export
  LINALG::Export(*((*mkxw)(0)), *mkxwv);
  LINALG::Export(*mkxwv, *mkxwstate_);
  LINALG::Export(*mkxwv, *mkv);
  LINALG::Export(*mkv, *mkstate_);


  return;
}  // end CalcMK

/*----------------------------------------------------------------------*
 |  Write enriched dofs in standard dofs for output            bk 09/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::XWall::GetOutputVector(Teuchos::RCP<Epetra_Vector> vel)
{
  Teuchos::RCP<Epetra_Vector> velenr =
      Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
  for (int i = 0; i < xwallrownodemap_->NumMyElements(); ++i)
  {
    int xwallgid = xwallrownodemap_->GID(i);
    DRT::Node* xwallnode = discret_->gNode(xwallgid);
    if (!xwallnode) dserror("Cannot find node");

    int firstglobaldofid = discret_->Dof(xwallnode, 0);
    int firstlocaldofid = discret_->DofRowMap()->LID(firstglobaldofid);

    int err = velenr->ReplaceMyValue(firstlocaldofid, 0, (*vel)[firstlocaldofid + 4]);
    err += velenr->ReplaceMyValue(firstlocaldofid + 1, 0, (*vel)[firstlocaldofid + 5]);
    err += velenr->ReplaceMyValue(firstlocaldofid + 2, 0, (*vel)[firstlocaldofid + 6]);
    err += velenr->ReplaceMyValue(firstlocaldofid + 3, 0, (*vel)[firstlocaldofid + 7]);
    if (err != 0) dserror("error during replacemyvalue");
  }
  return velenr;
}

/*----------------------------------------------------------------------*
 |  transfer tauw for turbulent inflow channel                 bk 09/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::XWall::GetTauw() { return tauw_; }

/*----------------------------------------------------------------------*
 |  transfer tauw for turbulent inflow channel                 bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::TransferAndSaveTauw()
{
  if (turbulent_inflow_condition_->IsActive())
  {
    LINALG::Export(*tauw_, *oldtauw_);
    LINALG::Export(*inctauw_, *oldinctauw_);
    turbulent_inflow_condition_->Transfer(oldtauw_, oldtauw_, 0.0);
    turbulent_inflow_condition_->Transfer(oldinctauw_, oldinctauw_, 0.0);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  transfer tauw for turbulent inflow channel                 bk 09/14 |
 *----------------------------------------------------------------------*/
void FLD::XWall::OverwriteTransferredValues()
{
  if (turbulent_inflow_condition_->IsActive())
  {
    Teuchos::RCP<Epetra_Vector> inctauwtmp =
        Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    LINALG::Export(*inctauw_, *inctauwtmp);
    Teuchos::RCP<Epetra_Vector> tauwtmp =
        Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    LINALG::Export(*tauw_, *tauwtmp);

    for (int i = 0; i < discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      int xwallgid = discret_->NodeRowMap()->GID(i);
      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if (!xwallnode) dserror("Cannot find node");
      std::vector<DRT::Condition*> nodecloudstocouple;
      xwallnode->GetCondition("TransferTurbulentInflow", nodecloudstocouple);
      if (not nodecloudstocouple.empty())
      {
        // usually we will only have one condition in nodecloudstocouple
        // but it doesn't hurt if there are several ones
        for (std::vector<DRT::Condition*>::iterator cond = nodecloudstocouple.begin();
             cond != nodecloudstocouple.end(); ++cond)
        {
          const std::string* mytoggle = (*cond)->Get<std::string>("toggle");
          if (*mytoggle == "slave")
          {
            inctauwtmp->ReplaceMyValue(i, 0, (*oldinctauw_)[i]);
            tauwtmp->ReplaceMyValue(i, 0, (*oldtauw_)[i]);
          }
        }
      }
    }

    LINALG::Export(*inctauwtmp, *inctauw_);
    LINALG::Export(*tauwtmp, *tauw_);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Read Restart                                               bk 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XWall::ReadRestart(IO::DiscretizationReader& reader)
{
  Teuchos::RCP<Epetra_Vector> tauw =
      Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
  reader.ReadVector(tauw, "xwall_tauw");
  LINALG::Export(*tauw, *tauw_);
  LINALG::Export(*tauw, *tauwxwdis_);

  restart_wss_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
  reader.ReadVector(restart_wss_, "wss");
  return;
}



/*----------------------------------------------------------------------*
 |  treat Dirichlet inflow                                     bk 04/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::XWall::FixDirichletInflow(Teuchos::RCP<Epetra_Vector> trueresidual)
{
  // copy for safety reasons
  Teuchos::RCP<Epetra_Vector> fixedtrueresidual =
      Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
  fixedtrueresidual->Update(1.0, *trueresidual, 0.0);

  // fix nodal forces on dirichlet inflow surfaces
  if (fix_residual_on_inflow_)
  {
    Teuchos::RCP<Epetra_Vector> res =
        Teuchos::rcp(new Epetra_Vector(*(discret_->DofColMap()), true));
    LINALG::Export(*trueresidual, *res);
    for (int j = 0; j < xwallrownodemap_->NumMyElements(); ++j)
    {
      int xwallgid = xwallrownodemap_->GID(j);

      DRT::Node* xwallnode = discret_->gNode(xwallgid);
      if (!xwallnode) dserror("Cannot find node");
      std::vector<DRT::Condition*> periodiccond;
      xwallnode->GetCondition("SurfacePeriodic", periodiccond);

      bool includedofs = true;
      if (not periodiccond.empty())
      {
        for (unsigned numcondper = 0; numcondper < periodiccond.size(); ++numcondper)
        {
          const std::string* mymasterslavetoggle =
              periodiccond[numcondper]->Get<std::string>("Is slave periodic boundary condition");
          if (*mymasterslavetoggle == "Slave")
          {
            includedofs = false;
          }
        }
      }
      if (includedofs)
      {
        if (discret_->NodeRowMap()->MyGID(xwallgid))
        {
          std::vector<DRT::Condition*> dircond;
          xwallnode->GetCondition("Dirichlet", dircond);

          std::vector<DRT::Condition*> stresscond;
          xwallnode->GetCondition("FluidStressCalc", stresscond);

          int numdf = discret_->NumDof(xwallnode);

          if ((not dircond.empty()) && (not stresscond.empty()) && numdf > 5)
          {
            bool isuglydirnode = false;
            for (unsigned numcond = 0; numcond < dircond.size(); ++numcond)
            {
              const std::vector<int>* flag = dircond[numcond]->Get<std::vector<int>>("onoff");

              if ((*flag)[4] or (*flag)[5] or (*flag)[6]) isuglydirnode = true;
            }

            if (isuglydirnode)
            {
              //
              // the new node has to be on these as well
              //  std::vector<DRT::Condition*> dircond;
              //    discret_->GetCondition("FluidStressCalc",dircond);
              DRT::Element** surrele = xwallnode->Elements();

              // loop over all surrounding elements and find indices of node k, l which is closes
              // while fulfilling all criteria
              double founddist = 1e9;
              int foundk = -1;
              int foundl = -1;
              for (int k = 0; k < (xwallnode->NumElement()); ++k)  // loop over elements
              {
                DRT::Node** test = surrele[k]->Nodes();
                for (int l = 0; l < surrele[k]->NumNode(); ++l)  // loop over nodes of element
                {
                  // it has to be on fluidstresscalc
                  // it may not be a dirichlet inflow node
                  // get Dirichlet conditions
                  std::vector<DRT::Condition*> stresscond;
                  test[l]->GetCondition("FluidStressCalc", stresscond);
                  int numdf = discret_->NumDof(test[l]);
                  if (not stresscond.empty() and numdf > 5)
                  {
                    std::vector<DRT::Condition*> dircond;
                    test[l]->GetCondition("Dirichlet", dircond);
                    bool isuglydirnode = false;
                    if (dircond.empty())
                    {
                      test[l]->GetCondition("FSICoupling", dircond);
                      if (dircond.empty())
                        dserror("this should be a Dirichlet or fsi coupling node");
                    }
                    else
                    {
                      for (unsigned numcond = 0; numcond < dircond.size(); ++numcond)
                      {
                        const std::vector<int>* flag =
                            dircond[numcond]->Get<std::vector<int>>("onoff");
                        if ((*flag)[4] or (*flag)[5] or (*flag)[6]) isuglydirnode = true;
                      }
                    }

                    if (not isuglydirnode)
                    {
                      const double* x = test[l]->X();
                      double dist = abs(x[0] - xwallnode->X()[0]) + abs(x[1] - xwallnode->X()[1]) +
                                    abs(x[2] - xwallnode->X()[2]);
                      if (founddist > dist)
                      {
                        founddist = dist;
                        foundk = k;
                        foundl = l;
                      }
                    }
                  }
                }
              }

              if (foundk < 0 or foundl < 0) dserror("haven't found required node");

              DRT::Node** test = surrele[foundk]->Nodes();

              int firstglobaldofidtoreplace = discret_->Dof(xwallnode, 0);
              int secondglobaldofidtoreplace = discret_->Dof(xwallnode, 0) + 1;
              int thirdglobaldofidtoreplace = discret_->Dof(xwallnode, 0) + 2;
              int firstglobaldofidnewvalue = discret_->Dof(test[foundl], 0);

              int firstlocaldofidnewvalue = discret_->DofColMap()->LID(firstglobaldofidnewvalue);
              // half because the area is half on a boundary node compared to an inner node
              double newvalue1 = 0.5 * (*res)[firstlocaldofidnewvalue];
              double newvalue2 = 0.5 * (*res)[firstlocaldofidnewvalue + 1];
              double newvalue3 = 0.5 * (*res)[firstlocaldofidnewvalue + 2];

              int err =
                  fixedtrueresidual->ReplaceGlobalValues(1, &newvalue1, &firstglobaldofidtoreplace);
              err = fixedtrueresidual->ReplaceGlobalValues(
                  1, &newvalue2, &secondglobaldofidtoreplace);
              err =
                  fixedtrueresidual->ReplaceGlobalValues(1, &newvalue3, &thirdglobaldofidtoreplace);
              if (err != 0) dserror("something wrong");
            }
          }
        }
      }
    }
  }
  return fixedtrueresidual;
}



/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 01/15 |
 *----------------------------------------------------------------------*/
FLD::XWallAleFSI::XWallAleFSI(Teuchos::RCP<DRT::Discretization> dis, int nsd,
    Teuchos::RCP<Teuchos::ParameterList>& params, Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
    Teuchos::RCP<FLD::UTILS::StressManager> wssmanager, Teuchos::RCP<Epetra_Vector> dispnp,
    Teuchos::RCP<Epetra_Vector> gridv)
    : XWall(dis, nsd, params, dbcmaps, wssmanager), mydispnp_(dispnp), mygridv_(gridv)
{
  incwdistxwdis_ = Teuchos::rcp(new Epetra_Vector(*(xwdiscret_->NodeColMap()), true));
  return;
}

void FLD::XWallAleFSI::UpdateWDistWALE()
{
  // save old one for projection
  incwdistxwdis_->Update(1.0, *wdistxwdis_, 0.0);

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> z = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));

  // fill vectors with coords
  for (int j = 0; j < xwallrownodemap_->NumMyElements(); ++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (not discret_->NodeRowMap()->MyGID(xwallgid))  // just in case
      dserror("not on proc");
    DRT::Node* xwallnode = discret_->gNode(xwallgid);
    if (!xwallnode) dserror("Cannot find node");

    int firstglobaldofid = discret_->Dof(xwallnode, 0);
    int firstlocaldofid = discret_->DofRowMap()->LID(firstglobaldofid);

    int err = x->ReplaceMyValue(j, 0, (xwallnode->X())[0] + (*mydispnp_)[firstlocaldofid]);
    err += y->ReplaceMyValue(j, 0, (xwallnode->X())[1] + (*mydispnp_)[firstlocaldofid + 1]);
    err += z->ReplaceMyValue(j, 0, (xwallnode->X())[2] + (*mydispnp_)[firstlocaldofid + 2]);
    if (err > 0) dserror("something wrong");
  }

  Teuchos::RCP<Epetra_Vector> wdistx = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> wdisty = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));
  Teuchos::RCP<Epetra_Vector> wdistz = Teuchos::rcp(new Epetra_Vector(*xwallrownodemap_, true));

  // project coordinates of the closest wall node to the node itself
  tauwcouplingmattrans_->Multiply(true, *x, *wdistx);
  tauwcouplingmattrans_->Multiply(true, *y, *wdisty);
  tauwcouplingmattrans_->Multiply(true, *z, *wdistz);

  // get delta
  wdistx->Update(-1.0, *x, 1.0);
  wdisty->Update(-1.0, *y, 1.0);
  wdistz->Update(-1.0, *z, 1.0);

  // fill vectors with coords
  for (int j = 0; j < xwallrownodemap_->NumMyElements(); ++j)
  {
    int xwallgid = xwallrownodemap_->GID(j);

    if (not discret_->NodeRowMap()->MyGID(xwallgid))  // just in case
      dserror("not on proc");
    DRT::Node* xwallnode = discret_->gNode(xwallgid);
    if (!xwallnode) dserror("Cannot find node");
    double x = (*wdistx)[j];
    double y = (*wdisty)[j];
    double z = (*wdistz)[j];
    double newwdist = sqrt(x * x + y * y + z * z);
    int err = walldist_->ReplaceMyValue(j, 0, newwdist);
    if (err > 0) dserror("something wrong");
  }

  LINALG::Export(*walldist_, *wdist_);
  LINALG::Export(*walldist_, *wdistxwdis_);
  // save old one for projection
  incwdistxwdis_->Update(1.0, *wdistxwdis_, -1.0);

  double mean = 0.0;
  walldist_->MeanValue(&mean);

  if (myrank_ == 0)
    std::cout << "the new mean distance from the wall of all XWall nodes is: " << mean << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWallAleFSI::SetXWallParams(Teuchos::ParameterList& eleparams)
{
  XWall::SetXWallParams(eleparams);
  return;
}

/*----------------------------------------------------------------------*
 |  Set params required to build the shape functions           bk 08/14 |
 |  Used for the xwdiscret_, which is redistributed                     |
 *----------------------------------------------------------------------*/
void FLD::XWallAleFSI::SetXWallParamsXWDis(Teuchos::ParameterList& eleparams)
{
  XWall::SetXWallParamsXWDis(eleparams);
  // params required for the shape functions
  eleparams.set("incwalldist", incwdistxwdis_);
  Teuchos::RCP<Epetra_Vector> xwdisdispnp = LINALG::CreateVector(*(xwdiscret_->DofRowMap()), true);
  LINALG::Export(*mydispnp_, *xwdisdispnp);
  Teuchos::RCP<Epetra_Vector> xwdisgridv = LINALG::CreateVector(*(xwdiscret_->DofRowMap()), true);
  LINALG::Export(*mygridv_, *xwdisgridv);

  xwdiscret_->SetState("dispnp", xwdisdispnp);
  xwdiscret_->SetState("gridv", xwdisgridv);

  return;
}

/*----------------------------------------------------------------------*
 |  Routine to update Tauw                                     bk 07/14 |
 *----------------------------------------------------------------------*/
void FLD::XWallAleFSI::UpdateTauW(int step, Teuchos::RCP<Epetra_Vector> trueresidual, int itnum,
    Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> veln)
{
  UpdateWDistWALE();
  FLD::XWall::UpdateTauW(step, trueresidual, itnum, accn, velnp, veln);
  if (tauwtype_ == INPAR::FLUID::constant)
  {
    if (proj_)
    {
      if (myrank_ == 0) std::cout << "  L2-project... ";

      L2ProjectVector(veln, Teuchos::null, accn);

      // at the beginning of this time step they are equal -> calculate only one of them
      velnp->Update(1.0, *veln, 0.0);

      if (myrank_ == 0) std::cout << "done!" << std::endl;
    }
    else
      dserror("projection required for ale case even with constant tauw, since wdist is updating");
  }

  return;
}
