/*!----------------------------------------------------------------------
\file fluidwomersleycondition.cpp
\brief evaluation of womersley flow bc

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/


#include <stdio.h>

#include "fluid_volumetric_surfaceFlow_condition.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_ana.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 09/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::FluidVolumetricSurfaceFlowWrapper(RCP<DRT::Discretization> actdis,
                                                                                 IO::DiscretizationWriter& output,
                                                                                 double dta) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{


  //--------------------------------------------------------------------
  // extract the womersley boundary dof
  //--------------------------------------------------------------------

  // get the fluid map extractor
  womersley_mp_extractor_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  womersley_mp_extractor_->Setup(*actdis);

  // creat the boundary map vector
  //  womersleybc_ = LINALG::CreateVector(*(womersley_mp_extractor_->WomersleyCondMap()),true);
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  womersleyvbc_ = LINALG::CreateVector(*dofrowmap,true);


  // Get the surfaces to whome the Womersley flow profile must be applied
  std::vector<DRT::Condition*> womersleycond;
  discret_->GetCondition("VolumetricSurfaceFlowCond",womersleycond);
  int num_of_wom_conds = womersleycond.size();

  // Get the lines which define the surrounding nodes of the womersley surface
  std::vector<DRT::Condition*> womersley_border_nodes_cond;
  discret_->GetCondition("VolumetricFlowBorderNodesCond",womersley_border_nodes_cond);
  int num_of_borders   = womersley_border_nodes_cond.size();

  //--------------------------------------------------------------------
  // Make sure that both each surface has one and only one border
  //--------------------------------------------------------------------
  if (num_of_wom_conds != num_of_borders)
  {
    dserror("Each Womersley surface condition must have one and only one border condition");
    exit(0);
  }
  // Check if each surface has it's corresponding border
  for (unsigned int i = 0; i < womersleycond.size(); i++)
  {
    bool ConditionIsWrong = true;
    // get the Womersley surface ID
    int surfID = womersleycond[i]->GetInt("ConditionID");

    // loop over all of the border conditions
    for (unsigned int j = 0; j< womersley_border_nodes_cond.size(); j++)
    {
      // get the border ID
      int lineID = womersley_border_nodes_cond[j]->GetInt("ConditionID");
      if(lineID == surfID)
      {
        // Since the condition is ok then create the corresponding the condition
        RCP<FluidVolumetricSurfaceFlowBc> fvsf_bc = Teuchos::rcp(new FluidVolumetricSurfaceFlowBc(discret_, output_, dta,"VolumetricSurfaceFlowCond","VolumetricFlowBorderNodesCond", surfID, i, j));
        bool inserted = fvsf_map_.insert( std::make_pair( surfID, fvsf_bc ) ).second;
        if ( !inserted )
        {
          dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
          exit(0);
        }
        ConditionIsWrong = false;
        break;
      }
    }

    // if a surface womersley doesn't have a correspondiong border defined!
    if (ConditionIsWrong)
    {
      dserror("Each Womersley surface condition must have one and only one border condition");
      exit(1);
    }
  }

  return;
} // end FluidVolumetricSurfaceFlowWrapper


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Output (public)                                         ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::Output( IO::DiscretizationWriter&  output)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    mapiter->second->FluidVolumetricSurfaceFlowBc::Output(output,"VolumetricSurfaceFlowCond",mapiter->first);
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                    ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::ReadRestart( IO::DiscretizationReader& reader)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    mapiter->second->FluidVolumetricSurfaceFlowBc::ReadRestart(reader,"VolumetricSurfaceFlowCond",mapiter->first);
  }

  return;
}//FluidVolumetricSurfaceFlowWrapper::ReadRestart


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the velocities of the dof and the map          ismail 09/10 |
 | extractor of boundary condition                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::EvaluateVelocities(RCP<Epetra_Vector> velocities,
                                                                       double             time)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {

    double flowrate =  mapiter->second->FluidVolumetricSurfaceFlowBc::EvaluateFlowrate("VolumetricSurfaceFlowCond",time);

    mapiter->second->FluidVolumetricSurfaceFlowBc::EvaluateVelocities(flowrate,"VolumetricSurfaceFlowCond",time);

    discret_->SetState("velaf",velocities);
    mapiter->second->FluidVolumetricSurfaceFlowBc::CorrectFlowRate("VolumetricSurfaceFlowCond",FLD::calc_flowrate,time);

    mapiter->second->FluidVolumetricSurfaceFlowBc::SetVelocities(velocities);
  }

  return;

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the dof map of the womersley conditions        ismail 10/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::EvaluateCondMap(RCP<Epetra_Map> &  bcmap )
{

  //std::cout<<*(womersley_mp_extractor_->WomersleyCondMap())<<std::endl;
  bcmap =   Teuchos::rcp(new Epetra_Map(*(womersley_mp_extractor_->VolumetricSurfaceFlowCondMap())));
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the map extractor of the womersley conditions  ismail 10/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::EvaluateMapExtractor(RCP<FLD::UTILS::MapExtractor>&  mapextractor )
{

  //std::cout<<*(womersley_mp_extractor_->WomersleyCondMap())<<std::endl;
  mapextractor =   Teuchos::rcp(new FLD::UTILS::MapExtractor(*(womersley_mp_extractor_)));
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Destructor dtor (public)                               ismail 09/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidVolumetricSurfaceFlowWrapper::~FluidVolumetricSurfaceFlowWrapper()
{
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 09/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidVolumetricSurfaceFlowBc::FluidVolumetricSurfaceFlowBc(RCP<DRT::Discretization>  actdis,
                                                                       IO::DiscretizationWriter& output,
                                                                       double                    dta,
                                                                       std::string                    ds_condname,
                                                                       std::string                    dl_condname,
                                                                       int                       condid,
                                                                       int                       surf_numcond,
                                                                       int                       line_numcond ) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // get dof row map
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = actdis->DofRowMap();

  // -------------------------------------------------------------------
  // get condition
  // -------------------------------------------------------------------
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(ds_condname,conditions);

  // -------------------------------------------------------------------
  // some standard initialized variables
  // -------------------------------------------------------------------
  condid_    = condid;
  condnum_s_ = surf_numcond;
  condnum_l_ = line_numcond;
  correction_factor_ = 0.0;
  dta_ = dta;

  // get the cycle period size
  period_      =   (conditions[surf_numcond])->GetDouble("Period");

  // get the polynomial order of the profile
  order_       =   (conditions[surf_numcond])->GetInt("Order");

  // get the number of harmonics
  n_harmonics_ =   (conditions[surf_numcond])->GetInt("Harmonics");

  // get the profile type
  flowprofile_type_  = (*(conditions[surf_numcond])->Get<std::string>("ConditionType"));

  // get the prebiasing flag
  prebiasing_flag_   = (*(conditions[surf_numcond])->Get<std::string>("prebiased"));

  // -------------------------------------------------------------------
  // calculate the center of mass and varage normal of the surface
  // condition
  // -------------------------------------------------------------------
  RCP<std::vector<double> > cmass  = Teuchos::rcp(new std::vector<double>);
  RCP<std::vector<double> > normal = Teuchos::rcp(new std::vector<double>);
  this->CenterOfMassCalculation(cmass,normal,ds_condname);

  // get the normal
  normal_ =  Teuchos::rcp(new std::vector<double>(*normal));
  std::string normal_info = *(conditions[surf_numcond])->Get<std::string>("NORMAL");
  if(normal_info == "SelfEvaluateNormal")
  {
    if(!myrank_)
    {
      std::cout<<"Normal is automatically evaluated"<<std::endl;
    }
    vnormal_ =  Teuchos::rcp(new std::vector<double>(*normal));
  }
  else if (normal_info ==  "UsePrescribedNormal")
  {
    if(!myrank_)
    {
      std::cout<<"Normal is manually setup"<<std::endl;
    }
    vnormal_ =  Teuchos::rcp(new std::vector<double>);
    (*vnormal_)[0] = (conditions[surf_numcond])->GetDouble("n1");
    (*vnormal_)[1] = (conditions[surf_numcond])->GetDouble("n2");
    (*vnormal_)[2] = (conditions[surf_numcond])->GetDouble("n3");
  }
  else
  {
    dserror("[%s]: is not a defined normal evaluation type",normal_info.c_str());
    exit(1);
  }

  // get the center of mass
  std::string c_mass_info = *(conditions[surf_numcond])->Get<std::string>("CenterOfMass");
  if(c_mass_info == "SelfEvaluateCenterOfMass")
  {
    if(!myrank_)
    {
      std::cout<<"Center of mass is automatically evaluated"<<std::endl;
    }
    cmass_ =  Teuchos::rcp(new std::vector<double>(*cmass));
  }
  else if (c_mass_info ==  "UsePrescribedCenterOfMass")
  {
    if(!myrank_)
    {
      std::cout<<"Center of mass is manually setup"<<std::endl;
    }
    normal_ =  Teuchos::rcp(new std::vector<double>);
    (*cmass_)[0] = (conditions[surf_numcond])->GetDouble("c1");
    (*cmass_)[1] = (conditions[surf_numcond])->GetDouble("c2");
    (*cmass_)[2] = (conditions[surf_numcond])->GetDouble("c3");
  }
  else
  {
    dserror("[%s]: is not a defined center-of-mass evaluation type",normal_info.c_str());
    exit(1);
  }


  // check if the condition surface is a inlet or outlet
  std::string flow_dir = *(conditions[surf_numcond])->Get<std::string>("FlowType");
  if (flow_dir == "InFlow")
  {
    flow_dir_ = -1.0;
  }
  else if (flow_dir == "OutFlow")
  {
    flow_dir_ =  1.0;
  }
  else
  {
    dserror("[%s]: is not a defined flow-direction-type",normal_info.c_str());
    exit(1);
  }

  // check if the flow is with correction
  std::string corr_flag = *(conditions[surf_numcond])->Get<std::string>("CorrectionFlag");
  correct_flow_ = (corr_flag == "WithCorrection");

  // -------------------------------------------------------------------
  // create the flow rates vector
  // -------------------------------------------------------------------
  int num_steps = int(period_/dta) + 1;

  flowrates_ = Teuchos::rcp(new std::vector<double>(num_steps,0.0));

  if (prebiasing_flag_=="PREBIASED"||prebiasing_flag_=="FORCED")
  {
    for(unsigned int i=0;i<flowrates_->size();i++)
    {
      (*flowrates_)[i] = this->EvaluateFlowrate(ds_condname,dta*double(i));
    }
  }

  // -------------------------------------------------------------------
  // get the node row maps of the condition node
  // -------------------------------------------------------------------
  // evaluate the surface node row map
  //  cond_surfnoderowmap_ = DRT::UTILS::ConditionNodeRowMap(*discret_,"WomersleySurfaceFlowCond");
  const std::string womersley_cond_name = ds_condname;
  this->BuildConditionNodeRowMap(discret_, womersley_cond_name, condid_, condnum_s_,cond_surfnoderowmap_);

  // evaluate the line node row map
  //  cond_linenoderowmap_ = DRT::UTILS::ConditionNodeRowMap(*discret_,"WomersleyBorderNodesCond");
  const std::string border_cond_name = dl_condname;
  this->BuildConditionNodeRowMap(discret_, womersley_cond_name, condid_, condnum_l_,cond_linenoderowmap_);

  // -------------------------------------------------------------------
  // get the dof row maps of the condition node
  // -------------------------------------------------------------------
  // evaluate the surface dof row map
  this->BuildConditionDofRowMap(discret_, womersley_cond_name, condid_, condnum_s_,cond_dofrowmap_);

  // -------------------------------------------------------------------
  // calculate the normalized  of mass of the surface condition
  // -------------------------------------------------------------------
  this->EvalLocalNormalizedRadii(ds_condname,dl_condname);

  // -------------------------------------------------------------------
  // create cond_velocities and codition traction velocity terms
  // -------------------------------------------------------------------
  cond_velocities_  = LINALG::CreateVector(*cond_dofrowmap_,true);
  cond_traction_vel_= LINALG::CreateVector(*dofrowmap,true);

  // -------------------------------------------------------------------
  // Evaluate the area of the design surface.
  // This will also return the viscosity and density of the fluid
  // -------------------------------------------------------------------
  area_ = this->Area(density_, viscosity_,ds_condname,condid_);

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Calculate center of mass (public)                       ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::CenterOfMassCalculation(RCP<std::vector<double> > coords,
                                                                       RCP<std::vector<double> > normal,
                                                                       std::string                    ds_condname)
{
  // Evaluate center of mass
  *coords = std::vector<double>(3,0.0);
  *normal = std::vector<double>(3,0.0);

  // fill the list od element evaluation
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",FLD::center_of_mass_calc);
  eleparams.set<double>("total area", 0.0);
  eleparams.set<RCP<std::vector<double> > >("center of mass", coords);
  eleparams.set<RCP<std::vector<double> > >("normal", normal);

  const std::string condstring(ds_condname);

  discret_->EvaluateCondition(eleparams,cond_velocities_,condstring,condid_);

  // get center of mass in parallel case
  std::vector<double> par_coord (3,0.0);
  for (unsigned int i = 0; i< par_coord.size();i++)
  {
    // get the actaul area on the local processor
    double act_area = eleparams.get<double>("total area");

    // get the actaul coordinate values on the local processor
    double act_val = (*coords)[i] * act_area;
    double act_n_val = (*normal)[i] * act_area;

    // define the parallel values that will be summed ove all of the processors
    double par_area  = 0.0;
    double par_val   = 0.0;
    double par_n_val = 0.0;

    // Summ all of the local processor values in one values
    discret_->Comm().SumAll(&act_val   ,&par_val   ,1);
    discret_->Comm().SumAll(&act_n_val ,&par_n_val ,1);
    discret_->Comm().SumAll(&act_area  ,&par_area  ,1);

    // finally evaluate the actual center of mass and avarage normal
    (*coords)[i] = par_val/par_area;
    (*normal)[i] = par_n_val/par_area;
  }


#if 1
  // Print out the results
  if (myrank_ == 0)
  {
    // print the center of mass
    printf("center of mass of cond(%d) is evaluated:\n",condid_);
    for (unsigned int i = 0; i< coords->size();i++)
    {
      printf ("| %f |",(*coords)[i]);
    }
    printf("\n");

    // print the avarage normal
    printf("avarage normal of cond(%d) surface is:\n",condid_);
    for (unsigned int i = 0; i< coords->size();i++)
    {
      printf ("| %f |",(*normal)[i]);
    }
    printf("\n");
  }
#endif

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Calculate local normalized radii (public)               ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::EvalLocalNormalizedRadii(
  std::string ds_condname,
  std::string dl_condname)
//RCP<std::vector<double> > center_of_mass,
//RCP<std::vector<double> > avg_normal,
//int condid_)
{
  //--------------------------------------------------------------------
  // get the processor rank
  //--------------------------------------------------------------------
  int myrank = discret_->Comm().MyPID();

  local_radii_   = LINALG::CreateVector(*cond_surfnoderowmap_,true);

  border_radii_  = LINALG::CreateVector(*cond_surfnoderowmap_,true);
  //--------------------------------------------------------------------
  // get all of the border nodes
  //--------------------------------------------------------------------
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(dl_condname,conditions);

  DRT::Condition* cond = conditions[condnum_l_];

  //--------------------------------------------------------------------
  // get the nodes of the condition
  //--------------------------------------------------------------------
  const std::vector<int>* border_nodes = cond->Nodes();

  //--------------------------------------------------------------------
  // Create a map of the border nodes
  //--------------------------------------------------------------------
  std::map <int, std::vector<double> > border_nodes_coords;
  for (unsigned int i =0 ; i<border_nodes->size();i++)
  {
    std::vector<double> coords(3,0.0);
    int node_num = (*border_nodes)[i];

    bool inserted = border_nodes_coords.insert( std::make_pair( node_num, coords ) ).second;
    if ( !inserted )
    {
      dserror("There are more than one node of the same number. something went wrong");
      exit(0);
    }
  }

  //--------------------------------------------------------------------
  // find the coordinates of each border node and get it parallel
  //--------------------------------------------------------------------
  for (std::map<int, std::vector<double> >::iterator it = border_nodes_coords.begin(); it!=border_nodes_coords.end(); it++)
  {
    // define the coordinate of a border node
    std::vector<double>  xyze(3,0.0);



    // if a node doesn't exist on the proc, then use the (0,0,0) origin
    for (unsigned int i=0; i<xyze.size(); i++)
    {
      xyze[i] = 0.0;
    }

    // get the coordinate of a node, if it exists on the current proc
    if ( discret_->HaveGlobalNode(it->first))
    {
      // check if the node is not a gohst node
      if (discret_->gNode(it->first)->Owner() == myrank)
      {
        // get the coordinate of a node
        const double* x = discret_->gNode(it->first)->X();
        for (unsigned int i=0; i<xyze.size(); i++)
        {
          xyze[i] = x[i];
        }
      }
    }

    //------------------------------------------------------------------
    // distribute the coordinates to all of the processors
    //------------------------------------------------------------------
    for (unsigned int i=0; i<xyze.size(); i++)
    {
      // define the actual value
      double act_val = xyze[i];

      // define the parallel value
      double par_val   = 0.0;

      // Summ all of the local processor values in one values
      discret_->Comm().SumAll(&act_val, &par_val, 1);

      // set the coordinate of the nodes on the local processor
      (it->second)[i] = par_val;
    }
  }

  // -------------------------------------------------------------------
  // loop over each node and compare its distance to the
  // [CenterOfMass BorderNodes)
  // -------------------------------------------------------------------

  // get the dimension of the node
  const int dim = 3;

  // define the vector between center-of-mass and current-node
  LINALG::Matrix< (dim) ,1>  c_cnd(true);

  // define the vector between center-of-mass and border-node
  LINALG::Matrix< (dim) ,1>  c_bnd(true);

  // define the vector that is the nearest from right
  LINALG::Matrix< (dim) ,1>  v_right(true);

  // define the vector that is the nearest from left
  LINALG::Matrix< (dim) ,1>  v_left(true);


  // define a direction vector perpendicular to the vector
  // of center-of-mass and current-node and to the surface normal.
  // This vector is also used to define whether a certain vector
  // is in the [0,pi] or [pi,2pi] awy from the
  // "center-of-mass and current-node" vector.
  LINALG::Matrix< (dim) ,1>  dir_vec(true);

  for (int lid = 0; lid<cond_surfnoderowmap_->NumMyElements(); lid++)
  {
    // get the global id of the current node
    int gid = cond_surfnoderowmap_->GID(lid);

    double R_r = 0.0;
    double R_l = 0.0;
    double R   = 0.0;
    // check if the node exists on the current processor
    if (discret_->HaveGlobalNode(gid))
    {
      // check if the node is not a gohst node
      if (discret_->gNode(gid)->Owner() == myrank)
      {

        double border_raduis = 0.0;
        const double * curr_xyze = discret_->gNode(gid)->X();

        //----------------------------------------------------------------
        // loop over all of the border nodes
        //----------------------------------------------------------------
        double diff_error_l = 2.0;
        double diff_error_r = 2.0;

        bool isBorderNode = false;
        for (std::map<int, std::vector<double> >::iterator it = border_nodes_coords.begin(); it!=border_nodes_coords.end(); it++)
        {
          isBorderNode = false;

          const std::vector <double> bord_xyze = it->second;

          //--------------------------------------------------------------
          // build the cener-of-mass to current-node vector
          // build the cener-of-mass to border-node vector
          //--------------------------------------------------------------
          for (int index = 0; index< dim; index++)
          {
            c_cnd(index) = curr_xyze[index] - (*cmass_)[index];
            c_bnd(index) = bord_xyze[index] - (*cmass_)[index];
          }

          // calculate the raduis of the border node
          R             = c_cnd.Norm2();
          border_raduis = c_bnd.Norm2();

          if(it->first == gid)
          {
            isBorderNode = true;
            break;
          }

          // normalize the two vectors
          c_cnd.Scale(1.0/c_cnd.Norm2());
          c_bnd.Scale(1.0/c_bnd.Norm2());

          //--------------------------------------------------------------
          // find the closest two vectors by calculating the norm of the
          // difference between the two vectors
          //--------------------------------------------------------------
          LINALG::Matrix< (dim) ,1> diff = c_cnd;
          diff -= c_bnd;

          //--------------------------------------------------------------
          // evaluate the direction vector = normal_vec X ref_vec
          //--------------------------------------------------------------
          dir_vec(0) = c_cnd(1)*(*normal_)[2] - c_cnd(2)*(*normal_)[1];
          dir_vec(1) = c_cnd(2)*(*normal_)[0] - c_cnd(0)*(*normal_)[2];
          dir_vec(2) = c_cnd(0)*(*normal_)[1] - c_cnd(1)*(*normal_)[0];

          // if the boundary is from the left
          if(dir_vec.Dot(c_bnd)> 0)
          {
            if (diff.Norm2() <= diff_error_l)
            {
              diff_error_l = diff.Norm2();
              R_l = border_raduis;
              v_left = c_bnd;
            }
          }
          else
          {
            if (diff.Norm2() <= diff_error_r)
            {
              diff_error_r = diff.Norm2();
              R_r = border_raduis;
              v_right= c_bnd;
            }
          }
        }
        // ---------------------------------------------------------------
        // calculate the local raduis of the node
        //  1- Calculate the angle between R_right and [center  curr-node]
        //  2- Calculate the angle between R_left and R_right
        //  3- Calculate the local Raduis by interpolation
        //  P.S: intersection method might be more accurate. But since
        //       some nodes might be slightly out of the plane, such a
        //       method might fail.
        // ---------------------------------------------------------------

        if (!isBorderNode)
        {

	  double cos_rl = v_right.Dot(v_left);
	  double cos_r  = v_right.Dot(c_cnd );

	  double angle_rl = 0.0;
          double angle_r  = 0.0;
	  double border_raduis = 0.0;

	  if (cos_r>=1.0 || cos_r<=-1.0)
          {
            border_raduis = R_r;
          }
	  else if (cos_rl>=1.0 || cos_rl<=-1.0)
          {
            border_raduis = R_l;
          }
	  else
          {
            angle_rl = acos(v_right.Dot(v_left));
            angle_r  = acos(v_right.Dot(c_cnd ));
            border_raduis = R_r + (R_l - R_r)*angle_r/angle_rl;
          }

          // update local raduis
          R /= border_raduis;
        }
        else
        {
          R = 1.0;
        }
        // update local raduis
        local_radii_->ReplaceGlobalValues(1,&R,&gid);

        // update border raduis
        border_radii_->ReplaceGlobalValues(1,&border_raduis,&gid);
      }
    }
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Calculate local normalized radii (public)               ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::BuildConditionNodeRowMap(
  RCP<DRT::Discretization>      dis,
  const std::string                  condname,
  int                           condid,
  int                           condnum  ,
  RCP<Epetra_Map> &             cond_noderowmap)
{

  //--------------------------------------------------------------------
  // get the processor rank
  //--------------------------------------------------------------------
  int myrank = dis->Comm().MyPID();

  //--------------------------------------------------------------------
  // read in the condition
  //--------------------------------------------------------------------
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(condname,conditions);

  DRT::Condition* cond = conditions[condnum];

  //--------------------------------------------------------------------
  // get the nodes of the condition
  //--------------------------------------------------------------------
  const std::vector<int>* nodes = cond->Nodes();
  std::vector<int>  nodeids;

  //--------------------------------------------------------------------
  // check which nodes belong to the current proc
  //--------------------------------------------------------------------
  for (unsigned int i = 0; i< nodes->size(); i++)
  {
    int Id = (*nodes)[i];
    if (dis->HaveGlobalNode(Id))
    {
      // check if the node is not a gohst node
      if (dis->gNode(Id)->Owner() == myrank)
      {
        nodeids.push_back(Id);
      }
    }
  }

  //--------------------------------------------------------------------
  // create the node row map of the nodes on the current proc
  //--------------------------------------------------------------------
  cond_noderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,dis->Comm()));

}//BuildConditionNodeRowMap

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Calculate local normalized radii (public)               ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::BuildConditionDofRowMap(
  RCP<DRT::Discretization>      dis,
  const std::string                  condname,
  int                           condid,
  int                           condnum  ,
  RCP<Epetra_Map> &             cond_dofrowmap)
{

  //--------------------------------------------------------------------
  // get the processor rank
  //--------------------------------------------------------------------
  int myrank = dis->Comm().MyPID();

  //--------------------------------------------------------------------
  // read in the condition
  //--------------------------------------------------------------------
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(condname,conditions);

  DRT::Condition* cond = conditions[condnum];

  //--------------------------------------------------------------------
  // get the nodes of the condition
  //--------------------------------------------------------------------
  const std::vector<int>* nodes = cond->Nodes();
  std::vector<int>  dofids;

  //--------------------------------------------------------------------
  // check which nodes belong to the current proc
  //--------------------------------------------------------------------
  for (unsigned int i = 0; i< nodes->size(); i++)
  {
    int Id = (*nodes)[i];
    if (dis->HaveGlobalNode(Id))
    {
      // check if the node is not a gohst node
      DRT::Node * node = dis->gNode(Id);
      if (node->Owner() == myrank)
      {
        for (int ldof = 0; ldof<3; ldof++)
        {
          dofids.push_back(dis->Dof(node, ldof));
        }
      }
    }
  }

  //--------------------------------------------------------------------
  // create the node row map of the nodes on the current proc
  //--------------------------------------------------------------------
  cond_dofrowmap = Teuchos::rcp(new Epetra_Map(-1,dofids.size(),&dofids[0],0,dis->Comm()));

}//FluidVolumetricSurfaceFlowBc::BuildConditionDofRowMap


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Output (public)                                         ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::Output(
  IO::DiscretizationWriter&  output,
  std::string                     ds_condname,
  int                        condnum )
{

  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!

  std::stringstream stream_rad, stream_brad, stream1, stream2, stream3, stream4;


  // write the flowrates of the previous period
  stream_rad <<ds_condname<<"_radii"<<condnum;
  output.WriteVector(stream_rad.str(),local_radii_);
  // write the flowrates of the previous period
  stream_brad <<ds_condname<<"_border_radii"<<condnum;
  output.WriteVector(stream_brad.str(),border_radii_);

  // write the flowrates of the previous period
  stream1 << ds_condname<<"_flowrates"<<condnum;
  output.WriteRedundantDoubleVector(stream1.str(),flowrates_);

  // write the time step
  stream2 << ds_condname<<"_dt"<<condnum;
  output.WriteDouble(stream2.str(), dta_);

  // write the condition velocities
  stream3 << ds_condname<<"_velocities"<<condnum;
  output.WriteVector(stream3.str(),cond_velocities_);

  stream4 << ds_condname<<"_traction_vel_component"<<condnum;
  output.WriteVector(stream4.str(),cond_traction_vel_);
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                    ismail 11/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::ReadRestart(
  IO::DiscretizationReader&  reader,
  std::string                     ds_condname,
  int                        condnum )
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!
  std::stringstream stream_rad, stream_brad, stream1, stream2, stream3, stream4;


  // write the flowrates of the previous period
  stream_rad <<ds_condname<<"_radii"<<condnum;
  reader.ReadVector(local_radii_,stream_rad.str());

  // write the flowrates of the previous period
  stream_brad <<ds_condname<<"_border_radii"<<condnum;
  reader.ReadVector(border_radii_,stream_brad.str());

  // old time step size
  stream2<< ds_condname<<"_dt"<<condnum;
  double odta = reader.ReadDouble(stream2.str());

  // write the condition velocities
  stream3 << ds_condname<<"_velocities"<<condnum;
  reader.ReadVector(cond_velocities_,stream3.str());

  stream4 << ds_condname<<"_traction_vel_component"<<condnum;
  reader.ReadVector(cond_traction_vel_,stream4.str());




  // get time step of the current problems
  double ndta = dta_;

  // -------------------------------------------------------------------
  // Read in the flowrates values and the flowrates position
  // -------------------------------------------------------------------
  stream1 << ds_condname<<"_flowrates"<<condnum;

  // read in flowrates
  reader.ReadRedundantDoubleVector(flowrates_ ,stream1.str());

  // read in the flowrates' position
  flowratespos_ = 0;

  // evaluate the new flowrate vector and position if time stepping
  // is changed
  if (odta != ndta)
  {
    // Get old flowrates Vector size
    int oQSize = (int)flowrates_->size();

    // Calculate new flowrates Vector size
    int nQSize = (int)(double(oQSize)*odta/ndta);

    // evaluate the new flowrates vector
    int nq_pos = 0;
    RCP<std::vector<double> > nq = Teuchos::rcp(new std::vector<double>(nQSize,0.0));
    this->Interpolate(flowrates_,nq,flowratespos_,nq_pos,period_);

    // store new values in class
    flowratespos_ = nq_pos;
    flowrates_    = nq;
  }

#if 0
  int myrank = discret_->Comm().MyPID();
  if (!myrank)
  {
    for(int i=0;i<flowrates_->size();i++)
    {
      std::cout<<"Flowrate"<<condnum<<": "<<(*flowrates_)[i]<<" time:"<<double(i)*odta<<std::endl;
    }
  }
#endif

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluates the Velocities (public)                       ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::EvaluateVelocities(
  double             flowrate,
  std::string             ds_condname,
  double             time)
{
  //--------------------------------------------------------------------
  // get the processor rank
  //--------------------------------------------------------------------
  //  int myrank = discret_->Comm().MyPID();

  double time_in_a_period = fmod(time,period_);//time - period_*floor(time/period_);
  // get the flowrate position
  int position = int(time_in_a_period/dta_+0.5);

  // insert flowrate into the flowrates vector
  (*flowrates_)[position] = flowrate;

  if(time_in_a_period < dta_)
  {
    (*flowrates_)[0] = flowrate;
    (*flowrates_)[flowrates_->size()-1] = flowrate;
  }

#if 0
  std::cout<<"condition("<<condid_<<"): has position: "<<position<<std::endl;
  std::cout<<"condition("<<condid_<<"): has area: "<<area_<<std::endl;
  if (!myrank)
  {
    for(int i=0;i<flowrates_->size();i++)
    {
      std::cout<<"Flowrate"<<condid_<<": "<<(*flowrates_)[i]<<" time: "<<double(i)*dta_<<std::endl;
    }
  }
#endif

  RCP<Teuchos::ParameterList>  params = Teuchos::rcp(new Teuchos::ParameterList);

  params->set<int>("Number of Harmonics",n_harmonics_);
  // condition id
  params->set<int>("Condition ID", condid_);
  // history of flowrates at the outlet
  params->set<RCP<std::vector<double> > >("Flowrates", flowrates_);
  // the velocity position
  params->set<int>("Velocity Position",position);
  // the flow type
  params->set<std::string>("flowrate type",flowprofile_type_);
  // time
  params->set<double>("time",time_in_a_period);
  // period of a cycle
  params->set<double>("period",period_);
  // polynomial order
  params->set<int>("polynomial order",order_);
  // flow direction
  params->set<double>("flow direction",flow_dir_);
  // surface area
  params->set<double>("area",area_);


  this->Velocities( discret_,
                    cond_velocities_,
                    cond_surfnoderowmap_,
                    local_radii_,
                    border_radii_,
                    vnormal_,
                    params);


}//FLD::UTILS::FluidWomersleyBc::EvaluateVelocities

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluates the Velocity componets of the traction        ismail 05/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::ResetTractionVelocityComp()
{
  cond_traction_vel_->Scale(0.0);
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluates the Velocity componets of the traction        ismail 05/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::EvaluateTractionVelocityComp(
  std::string             condname,
  double             flowrate,
  int                condid_,
  double             time,
  double             theta,
  double             dta)
{
  Teuchos::ParameterList eleparams;
  //  double norm2= 0.0;

  eleparams.set("thsl",theta*dta);
  eleparams.set("condition velocities",cond_velocities_);
  eleparams.set("condition dofrowmap",cond_dofrowmap_);

  discret_->SetState("velaf",cond_velocities_);

  eleparams.set("flowrate",flowrate);
  eleparams.set("area",area_);

#if 1
    eleparams.set<int>("action",FLD::traction_velocity_component);
    eleparams.set("velocities",cond_velocities_);
    discret_->EvaluateCondition(eleparams,cond_traction_vel_,condname,condid_);
#else
    eleparams.set<int>("action",FLD::traction_Uv_integral_component);
    discret_->EvaluateCondition(eleparams,cond_traction_vel_,condname,condid_);
#endif
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluates the Flowrate  (public)                        ismail 04/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
double FLD::UTILS::FluidVolumetricSurfaceFlowBc::EvaluateFlowrate (
  //  RCP<Epetra_Vector> bcdof,
  std::string             ds_condname,
  double             time)
{
  // -------------------------------------------------------------------
  // get curve information
  // -------------------------------------------------------------------

  // get condition
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(ds_condname,conditions);
  DRT::Condition* condition = conditions[condnum_s_];

  // get curve and curve_factor
  const std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
  double curvefac = 1.0;
  const std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

  // evaluate the current flowrate value
  double flowrate = 0.0;
  if((*curve)[0]>=0)
  {
    curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
    flowrate = (*vals)[0]*curvefac;
  }

  return flowrate;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluates the Velocities (public)                       ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::Velocities(
  RCP<DRT::Discretization>  disc,
  RCP<Epetra_Vector>        bcdof,
  RCP<Epetra_Map>           cond_noderowmap,
  RCP<Epetra_Vector>        local_radii,
  RCP<Epetra_Vector>        border_radii,
  RCP<std::vector<double> > normal,
  RCP<Teuchos::ParameterList>        params)


{

  //--------------------------------------------------------------------
  // get the processor rank
  //--------------------------------------------------------------------
  int myrank = disc->Comm().MyPID();


  // -------------------------------------------------------------------
  // Read in all of the parameters
  // -------------------------------------------------------------------
  // number of harmonics
  int n_harmonics = params->get<int>("Number of Harmonics");

  // condition id
  int condid      = params->get<int>("Condition ID");
  // history of avarage velocities at the outlet
  RCP<std::vector<double> > flowrates;
  flowrates            = params->get<RCP<std::vector<double> > >("Flowrates");
  RCP<std::vector<double> > velocities = Teuchos::rcp(new std::vector<double>(*flowrates));

  // the velocity position
  int velocityposition = params->get<int>("Velocity Position");
  // the flow type
  std::string  flowType     = params->get<std::string>("flowrate type");
  // time
  double time          = params->get<double>("time");
  // period of a cycle
  double period        = params->get<double>("period");
  // polynomial order
  int order            = params->get<int>("polynomial order");
  // flow direction
  //  double flow_dir      = params->get<double>("flow direction");

  // surface area
  double area          = params->get<double>("area");

  // -------------------------------------------------------------------
  // Get the volumetric flowrates Fourier coefficients
  // -------------------------------------------------------------------
  RCP<std::vector<std::complex<double> > > Vn;
  std::vector<double> Bn;

  //
  int cycle_num = 0;
  double in_cycle_time = 0.0;
  double time2= time- dta_;
  cycle_num = int(floor(time2/period));
  in_cycle_time = time2 - period*double(cycle_num);
  flowratespos_ = int((in_cycle_time/dta_));

  for(unsigned int i = 0; i< velocities->size(); i++)
  {
    (*velocities)[i] /= area*flow_dir_;
  }

  if (flowType == "WOMERSLEY")
  {
    if (n_harmonics< 1 )
    {
      dserror("The number of Womersley harmonics is %d (less than 1)",n_harmonics);
    }

    this->DFT(velocities,Vn,flowratespos_);

    Bn  = std::vector<double>(n_harmonics,0.0);

    double rl = real((*Vn)[0]);
    double im = imag((*Vn)[0]);

    Bn[0] = 0.5*rl;

    for(unsigned int k = 1; k<Bn.size(); k++)
    {
      //    double Mk   = sqrt(norm((*Qn)[k]));
      rl = real((*Vn)[k]);
      im = imag((*Vn)[k]);
      const double Mk = sqrt(rl*rl + im*im);
      const double Phik = atan2(-imag((*Vn)[k]),real((*Vn)[k]));


      // Bn[k] = Mk* cos(2.0*M_PI*double(k)*time/period - Phik);
      Bn[k] = Mk* cos(2.0*M_PI*double(k) - Phik);
    }
  }



  // -------------------------------------------------------------------
  // evaluate the avarage velocity and apply it to the design surface
  // -------------------------------------------------------------------
  // loop over all of the nodes
  for (int lid = 0; lid<cond_noderowmap->NumMyElements(); lid++)
  {
    int gid = cond_noderowmap->GID(lid);

    // check if the node exists on the current processor
    if (disc->HaveGlobalNode(gid))
    {
      const DRT::Node* node =  disc->gNode(gid);

      // check if the node is not a gohst node
      if (node->Owner() == myrank)
      {
        // loop over the dof of a map
        // eval the velocity of a dof
        double velocity = 0.0;
        double r = (*local_radii)[cond_noderowmap->LID(gid)];

        //------------------------------------------------------------
        // Check for the velocity profile type
        //------------------------------------------------------------

        // check for the polynomial type
        if (flowType == "POLYNOMIAL")
        {
          if (order_!=0)
          {
            velocity  = (1.0+2.0/double(order_))*(*velocities)[velocityposition];
            velocity *= this->PolynomailVelocity(r, order);
          }
          else
          {
            velocity  = (*velocities)[velocityposition];
          }
        }
        // else check for Womersley type
        else if (flowType == "WOMERSLEY")
        {
          double R = (*border_radii)[cond_noderowmap->LID(gid)];

          // first calculate the parabolic profile of the 0th
          // harmonics
          if (order_!=0)
          {
            velocity = (1.0+2.0/double(order))*Bn[0]*this->PolynomailVelocity(r, order);
          }
          else
          {
            velocity = Bn[0];
          }
          //velocity = 0.0;
          double velocity_wom = 0.0;
          double rl = real((*Vn)[0]);
          double im = imag((*Vn)[0]);
          for (unsigned int k = 1; k < Bn.size(); k++)
          {
            //    double Mk   = sqrt(norm((*Qn)[k]));
            rl = real((*Vn)[k]);
            im = imag((*Vn)[k]);
            const double Mk = sqrt(rl*rl + im*im);
            const double Phik = atan2(-imag((*Vn)[k]),real((*Vn)[k]));
            Bn[k] = Mk* cos(2.0*M_PI*double(k) - Phik);
            // Bn[k] = Mk* cos(2.0*M_PI*double(k)*time/period - Phik);

            velocity_wom += this->WomersleyVelocity(r,R,Mk,Phik,k,period);
          }
          velocity += velocity_wom;
        }
        else
        {
          dserror("[%s] in cond (%d): No such profile is defined. Please correct the input file ",flowType.c_str(),condid);
        }
        velocity *= flow_dir_;
        for (unsigned int ldof = 0; ldof< normal->size(); ldof++)
        {
          // get the global dof from using the local one
          int gdof = disc->Dof(node, ldof);

          //------------------------------------------------------------
          // Apply the velocity in the normal direction
          //------------------------------------------------------------
          double Vdof =  velocity * (*normal)[ldof];
          //          std::cout<<"Velocity["<<gdof<<"]: "<<Vdof<<std::endl;
          bcdof->ReplaceGlobalValues(1,&Vdof,&gdof);
        }
      }
    }
  }

}//FLD::UTILS::FluidVolumetricSurfaceFlowBc::Velocities



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Corrects the Flow Rate   (public)                       ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::CorrectFlowRate
(std::string             ds_condname,
 FLD::BoundaryAction action,
 double             time,
 bool               force_correction)
{

  if(!force_correction)
  {
    if(!correct_flow_)
    {
      return;
    }
  }

  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  int myrank  = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // calculate flow rate
  // -------------------------------------------------------------------
  discret_->SetState("velaf",cond_velocities_);

  double actflowrate = this->FlowRateCalculation(time,ds_condname,action,condid_);

  //  std::cout<<*cond_velocities_<<std::endl;

  // -------------------------------------------------------------------
  // evaluate the wanted flow rate
  // -------------------------------------------------------------------

  // get curve information
  double time_in_a_period = fmod(time,period_);//time - period_*floor(time/period_);

  // get the flowrate position
  int position = int(time_in_a_period/dta_+0.5);

  // insert flowrate into the flowrates vector
  double flowrate = (*flowrates_)[position];
  //  double flowrate = this->FlowRateCalculation(time,ds_condname,action,condid_);

  if(myrank == 0 )
  {
    double flow_error = 0.0;
    flow_error = actflowrate - flowrate;
    printf("DIR(%f): Flow_estimated = %f : Flow_wanted = %f : Flow_correction = %f\n",flow_dir_,actflowrate,flowrate,flow_error);
  }

  // loop over all of the nodes
  RCP<Epetra_Vector> correction_velnp = LINALG::CreateVector(*cond_dofrowmap_,true);

  RCP<Teuchos::ParameterList>  params = Teuchos::rcp(new Teuchos::ParameterList);

  params->set<int>("Number of Harmonics",0);
  // condition id
  params->set<int>("Condition ID", condid_);
  // history of avarage velocities at the outlet
  RCP<std::vector<double> > flowrates = Teuchos::rcp(new std::vector<double> );
  flowrates->push_back(1.0*area_);
  params->set<RCP<std::vector<double> > >("Flowrates", flowrates);
  // the velocity position
  params->set<int>("Velocity Position",0);
  // the flow type
  params->set<std::string>("flowrate type","POLYNOMIAL");
  // time
  params->set<double>("time",time_in_a_period);
  // period of a cycle
  params->set<double>("period",period_);
  // polynomial order
  params->set<int>("polynomial order",order_);
  // flow direction
  params->set<double>("flow direction",flow_dir_);
  // surface area
  params->set<double>("area",area_);

  this->Velocities( discret_,
                    correction_velnp,
                    cond_surfnoderowmap_,
                    local_radii_,
                    border_radii_,
                    vnormal_,
                    params);




  discret_->SetState("velaf",correction_velnp);
  double corrective_flowrate = this->FlowRateCalculation(time,ds_condname,action,condid_);
  //  discret_->ClearState();

  if (action == FLD::calc_flowrate)
  {
    correction_factor_ = (flowrate  - actflowrate)/(corrective_flowrate);

    double correction = 0.0;
    for (int lid = 0; lid<correction_velnp->MyLength(); lid++)
    {
      int gid = correction_velnp->Map().GID(lid);
      correction = correction_factor_*(*correction_velnp)[lid];

      int bc_lid = cond_velocities_->Map().LID(gid);
      (*cond_velocities_)[bc_lid] += correction;
    }
  }
  else
  {
    correction_factor_ = flowrate/(corrective_flowrate);
    correction_factor_ = sqrt(correction_factor_);
    double correction = 0.0;
    for (int lid = 0; lid<correction_velnp->MyLength(); lid++)
    {
      int gid = correction_velnp->Map().GID(lid);
      correction = correction_factor_*(*correction_velnp)[lid];

      int bc_lid = cond_velocities_->Map().LID(gid);
      (*cond_velocities_)[bc_lid] = correction;
    }
  }
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply velocities         (public)                       ismail 10/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::SetVelocities(RCP<Epetra_Vector> velocities)
{
  for (int lid = 0; lid <cond_velocities_->MyLength();lid++)
  {
    int gid    =  cond_velocities_->Map().GID(lid);
    double val = (*cond_velocities_)[lid];

    velocities->ReplaceGlobalValues(1,&val,&gid);
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Flow rate calculation                                      ac 03/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
  modified by chfoe 04/08

  Calculate the flow rate accros an impedance boundary surface

  Flow rates are
  (1) calculated for single element surfaces
  (2) added up over the elements of the single procs
  (3) communicated and added over the procs
  (4) and finally stored within the vector 'flowrates_'

  The vector of the flowrates holds the flow rate history of the
  very last cycle!

*/
double FLD::UTILS::FluidVolumetricSurfaceFlowBc::FlowRateCalculation(
  double time,
  std::string ds_condname,
  FLD::BoundaryAction action,
  int    condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",action);
  eleparams.set("total time",time);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> flowrates = LINALG::CreateVector(*dofrowmap,true);

  const std::string condstring(ds_condname);
  discret_->EvaluateCondition(eleparams,flowrates,condstring,condid);

  double local_flowrate = 0.0;
  for (int i=0; i < dofrowmap->NumMyElements(); i++)
  {
    local_flowrate +=((*flowrates)[i]);
  }

  double flowrate = 0.0;
  dofrowmap->Comm().SumAll(&local_flowrate,&flowrate,1);

  return flowrate;

}//FluidImplicitTimeInt::FlowRateCalculation


double FLD::UTILS::FluidVolumetricSurfaceFlowBc::PressureCalculation(
  double time,
  std::string ds_condname,
  std::string action,
  int    condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  double pressure = 0.0;
  eleparams.set("action",action);
  eleparams.set("Inlet integrated pressure",pressure);
  eleparams.set("total time",time);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> flowrates = LINALG::CreateVector(*dofrowmap,true);

  const std::string condstring(ds_condname);
  discret_->EvaluateCondition(eleparams,flowrates,condstring,condid);

  pressure = eleparams.get<double>("Inlet integrated pressure");
  std::cout<<"avg pressure: "<<pressure/area_<<std::endl;

  return pressure/area_;

}//FluidImplicitTimeInt::PressureCalculation

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Parabolic velocity at certain raduis and time         mueller 04/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
double FLD::UTILS::FluidVolumetricSurfaceFlowBc::PolynomailVelocity(double r,
                                                                    int order)
{
  return (1.0-pow(r,double(order)));

}// FLD::UTILS::FluidVolumetricSurfaceFlowBc::PolynomailVelocity


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Womersley velocity at certain raduis and time         mueller 04/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
double FLD::UTILS::FluidVolumetricSurfaceFlowBc::WomersleyVelocity(double r,
                                                                   double R,
                                                                   double Bn,
                                                                   double Phi,
                                                                   //complex<double> Bn,
                                                                   int    n,
                                                                   double time)
{

  //--------------------------------------------------------------------
  // define some variables
  //--------------------------------------------------------------------

  // velocity
  std::complex<double> velocity;

  // imaginary unit
  std::complex<double> i =  std::complex<double>(0.0,1.0);

  // exp^( i*w*t)
  double constexp = 2.0*M_PI*double(n)*time/period_-Phi;
  double realpart = cos(constexp);
  double imagpart = sin(constexp);
  std::complex<double> eiwt_phi (realpart,imagpart);

  // Jo_z
  std::complex<double> Jo_z;

  // J1_z
  std::complex<double> J1_z;

  // Jo_rz
  std::complex<double> Jo_rz;

  //--------------------------------------------------------------------
  // evaluate the nth harmonic velocity
  //--------------------------------------------------------------------

  // Womersley number
  //  double          alpha = R*sqrt(2.0*M_PI*double(n)/period_/viscosity_);
  double          alpha = R*sqrt(2.0*M_PI*double(n)/period_/(viscosity_/density_));



  // Bessel variable
  std::complex<double>     z = alpha*pow(i,1.5);

  // bessel functions
  Jo_z  = this->BesselJ01(z  ,false);
  J1_z  = this->BesselJ01(z  ,true );
  Jo_rz = this->BesselJ01(z*r,false);

  // velocity
  velocity =  (Bn)*(z*(Jo_z - Jo_rz)/(z*Jo_z - 2.0*J1_z))*eiwt_phi;

  // return the real part of the Womersley velocity
  return real(velocity);

}// FLD::UTILS::FluidVolumetricSurfaceFlowBc::PolynomailVelocity


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Womersley: Bessel functions of order 0 and 1          mueller 04/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
std::complex<double> FLD::UTILS::FluidVolumetricSurfaceFlowBc::BesselJ01(std::complex<double> z, bool order)
{
  // DESCRIPTION:
  // Bessel functions of order 0 (order==false) or 1 (order==true) are calculated for
  // a given argument z

  const int end = 70;
  std::complex<double> J(0.0,0.0);
  std::complex<double> Jmine(0.0,0.0);
  double fac = 1.0;

  int alpha = 1;
  if(order==false)
  {
    alpha = 0;
  }
  else
  {
    alpha = 1;
  }

  for(int m=0;m<end;m++)
  {
    double fac   = 1.0;
    double gamma = 1.0;
    for (int k=1;k<=m;k++)
    {
      fac *= (double)(k);
    }
    for (int k=1;k<=m+alpha;k++)
    {
      gamma *= (double)(k);
    }

    Jmine += std::pow(z*std::complex<double>(0.5),double(alpha))
      *  pow(-std::complex<double>(0.25)*z*z,double(m))
      /(std::complex<double> (fac) * std::complex<double> (gamma));

  }
#if 1
  // Bessel function of the first kind and order 0
  if(order==false)
  {
    for(int m=0;m<end;m++)
    {
      for(int k=2;k<=m;k++)
	fac *= (double)(k);
      J += (std::complex<double>)((double)(std::pow(-1.0,(double)(m)))/std::pow(fac,2.0))*
      pow((z/(std::complex<double>)(2.0)),(std::complex<double>)(2*m));
      fac = 1.0;
    }
    if(z == std::complex<double>(0.0,0.0))
      J= std::complex<double>(1.0,0.0);
  }
  // Bessel function of the first kind and order 1
  else
  {
    for(int m=0;m<end;m++)
    {
      for(int k=2;k<=m;k++)
	fac *= (double)(k);
      J += (std::complex<double>)((std::pow(-1.0,(double)(m)))/((double)(m+1)*std::pow(fac,2.0)))*
      std::pow((z/std::complex<double>(2.0)),(std::complex<double>)(2*m+1));
      fac = 1.0;
    }
  }
#endif

  if (fabs(real(Jmine)-  real(J))>0.001*fabs(real(J)) || fabs(imag(Jmine)- imag(J))>0.001*fabs(imag(J)))
  {
    dserror("J(%d) problems... [%f,%f]\t[%f,%f]",alpha,real(Jmine),imag(Jmine),real(J),imag(J));
  }

  return J;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Area calculation                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!

*/
double FLD::UTILS::FluidVolumetricSurfaceFlowBc::Area(
  double& density,
  double& viscosity,
  std::string  ds_condname,
  int     condid )
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",FLD::calc_area);
  eleparams.set<double>("area",0.0);
  eleparams.set<double>("viscosity",0.0);
  eleparams.set<double>("density",0.0);

  const std::string condstring(ds_condname);

  discret_->EvaluateCondition(eleparams,condstring,condid);

  double actarea = eleparams.get<double>("area");
  density = eleparams.get<double>("density");
  viscosity = eleparams.get<double>("viscosity");

  // find the lowest proc number that knows the material data
  int numproc = discret_->Comm().NumProc();
  int theproc = -1;   // the lowest proc that has the desired information
  std::vector<double> alldens(numproc);

  discret_->Comm().GatherAll( &density,&(alldens[0]),1 );
  for(int i=0; i<numproc; i++)
    if( alldens[i] > 0.0 )
    {
      theproc = i;
      break;
    }
  if(theproc < 0)
    dserror("Something parallel went terribly wrong!");

  // do the actual communication of density ...
  discret_->Comm().Broadcast(&density,1,theproc);
  // ... and viscosity
  discret_->Comm().Broadcast(&viscosity,1,theproc);

  // get total area in parallel case
  double pararea = 0.0;
  discret_->Comm().SumAll(&actarea,&pararea,1);

  if (myrank_ == 0)
  {
    std::cout << "Volumetric surface flow rate condition Id: " << condid << " area = " << pararea << std::endl;
  }
  return pararea;
}//FLD::UTILS::FluidVolumetricSurfaceFlowBc::Area

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Womersley: Discrete Fourier Transfomation              ismail 10/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::DFT(RCP<std::vector<double> >   f,
                                                   RCP<std::vector<std::complex<double> > > & F,
                                                   int starting_pos)
{
  //--------------------------------------------------------------------
  // Initialise the Fourier values
  //--------------------------------------------------------------------
  F =Teuchos::rcp(new std::vector<std::complex<double> >(f->size(),0.0));

  const double N  = double(f->size());
  const int fsize = f->size();

  //--------------------------------------------------------------------
  // Compute the Fourier values
  //--------------------------------------------------------------------
  for (int k = 0; k < fsize; k++)
  {
    (*F) [k] = std::complex<double>(0.0,0.0);
    for (int n = 0; n< fsize; n++)
    {
      int pos = 0;
      if (starting_pos - n >= 0)
      {
        pos = (starting_pos - n);
      }
      else
      {
        pos = fsize - (n - starting_pos);
      }

      double rl = (*f)[pos]*2.0/N*( cos(2.0*M_PI*double(k)*double(fsize-1 - n)/N));
      double im = (*f)[pos]*2.0/N*(-sin(2.0*M_PI*double(k)*double(fsize-1 - n)/N));

      (*F) [k] += std::complex<double> (rl,im);
    }

  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Interpolate values of Vector1 to fit in Vector2         ismail 04/10|
 |                                                                      |
 | V ^                                                                  |
 |   |                              +                                   |
 |   |                            , .                                   |
 |   |                          ,   .                                   |
 |   |                        ,     .                                   |
 |   |                      ,+      .                                   |
 |   |                    ,  .      .                                   |
 |   |                  //   .      .                                   |
 |   |                ,      .      .                                   |
 |   |              ,+       .      .                                   |
 |   |            ,  .       .      .                                   |
 |   |          ,    .       .      .                                   |
 |   |        +      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   +--------+------o---//--o------+----------->                       |
 |            T(i)   .              T(i+1)      t                       |
 |            .      .              .                                   |
 |            t(m)   t(m+j)         t(m+k)                              |
 |                                                                      |
 |                                                                      |
 |  (T) is the time step of the original vector                         |
 |  (t) is the time step of the new vector                              |
 |  1 - Loop over all intervals (T(i) and T(i+1))                       |
 |  2 - Check if V2 has any time steps between T(i) and T(i+1)          |
 |  3 - Using linear interpolation check get the value of V2 at t(m+j)  |
 |                                                                      |
 | *The advantage of this method is that is works for finer and coarser |
 |  interpolations.                                                     |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidVolumetricSurfaceFlowBc::Interpolate(RCP<std::vector<double> > V1,
                                                           RCP<std::vector<double> > V2,
                                                           int index1,
                                                           int & index2,
                                                           double period)
{
  // Get size of V1 and V2
  int n1 = V1->size();
  int n2 = V2->size();

  double TotalTime = period;

  // Get time step size of V1 and V2
  double dt1 = (TotalTime)/double(n1 - 1);
  double dt2 = (TotalTime)/double(n2 - 1);

  // defining some necessary variables
  double t1_1, t1_2;
  double v1_1, v1_2;

  // define t (time step of V2) and k (index of V2)
  double t = 0.0;
  int k    = 0;
  for (int i = 0; i< n1-1; i++)
  {
    // -----------------------------------------------------------------
    // Get V1 values at T(i) and T(i+1)
    // -----------------------------------------------------------------
    v1_1 = (*V1)[i];
    v1_2 = (*V1)[i+1];

    // -----------------------------------------------------------------
    // Calculate T(i) and T(i+1)
    // -----------------------------------------------------------------
    t1_1 = double(i)*dt1;
    t1_2 = double(i+1)*dt1;

    // -----------------------------------------------------------------
    // Evaluate V2 values between T(i) and  T(i+1)
    // -----------------------------------------------------------------
    while (t < t1_2)
    {
      // Evaluate value of V2 using Interpolation
      (*V2)[k] = (t1_2 - t)/dt1*v1_1 + (t - t1_1)/dt1*v1_2;
      // Increment k
      k++;
      // Increment t
      t += dt2;
    }
  }

  // -------------------------------------------------------------------
  // Finally resolve the last step where V2(n2) = V1(n1)
  // -------------------------------------------------------------------
  (*V2)[V2->size()-1] = (*V1)[V1->size()-1];


  // -------------------------------------------------------------------
  // Get the index of V2
  // where t = t     => dt1*index1 = dt2*index2
  //                              dt1
  //                 => index2 = ----- index1
  //                              dt2
  // -------------------------------------------------------------------
  index2 = int(double(index1)*(dt1/dt2));

}//FLD::UTILS::FluidVolumetricSurfaceFlowBc::interpolate

void FLD::UTILS::FluidVolumetricSurfaceFlowBc::UpdateResidual(RCP<Epetra_Vector>  residual )
{
  //  double norm2;
  //  cond_traction_vel_->Norm2(&norm2);
  //  std::cout<<"DIR("<<flow_dir_<<"): Residual Norm: "<<norm2<<std::endl;
  //  residual->Update(1.0,*cond_traction_vel_,1.0);
  residual->Update(1.0,*cond_traction_vel_,1.0);

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 04/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::TotalTractionCorrector::TotalTractionCorrector(RCP<DRT::Discretization> actdis,
                                                                   IO::DiscretizationWriter& output,
                                                                   double dta) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{


  //--------------------------------------------------------------------
  // extract the womersley boundary dof
  //--------------------------------------------------------------------

  // get the fluid map extractor
  traction_mp_extractor_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  traction_mp_extractor_->Setup(*actdis);

  // creat the boundary map vector
  // womersleybc_ = LINALG::CreateVector(*(womersley_mp_extractor_->WomersleyCondMap()),true);
  //  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  //  tractionvbc_ = LINALG::CreateVector(*dofrowmap,true);
  //  ExtractVolumetricSurfaceFlowCondVector(LINALG::CreateVector(*dofrowmap,true)


  // Get the surfaces to whome the traction flow profile must be applied
  std::vector<DRT::Condition*> tractioncond;
  discret_->GetCondition("TotalTractionCorrectionCond",tractioncond);
  int num_of_tr_conds = tractioncond.size();

  // Get the lines which define the surrounding nodes of the traction surface
  std::vector<DRT::Condition*> traction_border_nodes_cond;
  discret_->GetCondition("TotalTractionCorrectionBorderNodesCond",traction_border_nodes_cond);
  int num_of_borders   = traction_border_nodes_cond.size();

  //--------------------------------------------------------------------
  // Make sure that both each surface has one and only one border
  //--------------------------------------------------------------------
  if (num_of_tr_conds != num_of_borders)
  {
    dserror("Each Womersley surface condition must have one and only one border condition");
    exit(0);
  }
  // Check if each surface has it's corresponding border
  for (unsigned int i = 0; i < tractioncond.size(); i++)
  {
    bool ConditionIsWrong = true;
    // get the traction surface ID
    int surfID = tractioncond[i]->GetInt("ConditionID");

    // loop over all of the border conditions
    for (unsigned int j = 0; j< traction_border_nodes_cond.size(); j++)
    {
      // get the border ID
      int lineID = traction_border_nodes_cond[j]->GetInt("ConditionID");
      if(lineID == surfID)
      {
        // Since the condition is ok then create the corresponding the condition
        RCP<FluidVolumetricSurfaceFlowBc> fvsf_bc = Teuchos::rcp(new FluidVolumetricSurfaceFlowBc(discret_, output_, dta,"TotalTractionCorrectionCond","TotalTractionCorrectionBorderNodesCond", surfID, i, j));
        bool inserted = fvsf_map_.insert( std::make_pair( surfID, fvsf_bc ) ).second;
        if ( !inserted )
        {
          dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
          exit(0);
        }

        ConditionIsWrong = false;
        break;
      }
    }

    // if a surface traction doesn't have a correspondiong border defined!
    if (ConditionIsWrong)
    {
      dserror("Each Total traction correction surface condition must have one and only one border condition");
      exit(1);
    }
  }

  return;
} // end TotalTractionCorrector


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the velocities of the dof and the map          ismail 04/11 |
 | extractor of boundary condition                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::EvaluateVelocities(RCP<Epetra_Vector> velocities,
                                                            double             time,
                                                            double             theta,
                                                            double             dta)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    double flowrate = 0.0;

    if (mapiter->second->PrebiasingFlag()=="FORCED")
    {
      flowrate = mapiter->second->EvaluateFlowrate("TotalTractionCorrectionCond",time);
    }
    else
    {
      discret_->SetState("velaf",velocities);
      flowrate =    mapiter->second->FluidVolumetricSurfaceFlowBc::FlowRateCalculation(time,"TotalTractionCorrectionCond",FLD::calc_flowrate,mapiter->first);
    }

    mapiter->second->FluidVolumetricSurfaceFlowBc::EvaluateVelocities(flowrate,"TotalTractionCorrectionCond",time);

    discret_->SetState("velaf",velocities);

    //    mapiter->second->FluidVolumetricSurfaceFlowBc::CorrectFlowRate("TotalTractionCorrectionCond","calculate Uv integral component",time,true);
    mapiter->second->FluidVolumetricSurfaceFlowBc::CorrectFlowRate("TotalTractionCorrectionCond",FLD::calc_flowrate,time,true);

    mapiter->second->FluidVolumetricSurfaceFlowBc::ResetTractionVelocityComp();

    discret_->SetState("velaf",velocities);
    mapiter->second->FluidVolumetricSurfaceFlowBc::EvaluateTractionVelocityComp("TotalTractionCorrectionCond",flowrate,mapiter->first,time,theta,dta);
  }

  return;

}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the dof map of the womersley conditions        ismail 04/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::EvaluateCondMap(RCP<Epetra_Map> &  bcmap )
{

  //std::cout<<*(womersley_mp_extractor_->WomersleyCondMap())<<std::endl;
  bcmap =   Teuchos::rcp(new Epetra_Map(*(traction_mp_extractor_->TotalTractionCorrectionCondMap())));
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Evaluate the map extractor of the womersley conditions  ismail 04/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::EvaluateMapExtractor(RCP<FLD::UTILS::MapExtractor>&  mapextractor )
{

  //std::cout<<*(womersley_mp_extractor_->WomersleyCondMap())<<std::endl;
  mapextractor =   Teuchos::rcp(new FLD::UTILS::MapExtractor(*(traction_mp_extractor_)));
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Update residual                                         ismail 04/11 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::UpdateResidual(RCP<Epetra_Vector>  residual )
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    mapiter->second->FluidVolumetricSurfaceFlowBc::UpdateResidual(residual);
  }
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Output (public)                                         ismail 04/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::Output( IO::DiscretizationWriter&  output)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    mapiter->second->FluidVolumetricSurfaceFlowBc::Output(output,"TotalTractionCorrectionCond",mapiter->first);
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                    ismail 04/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::TotalTractionCorrector::ReadRestart( IO::DiscretizationReader& reader)
{
  std::map<const int, RCP<class FluidVolumetricSurfaceFlowBc> >::iterator mapiter;

  for (mapiter = fvsf_map_.begin(); mapiter != fvsf_map_.end(); mapiter++ )
  {
    mapiter->second->FluidVolumetricSurfaceFlowBc::ReadRestart(reader,"TotalTractionCorrectionCond",mapiter->first);
  }

  return;
}//FluidVolumetricSurfaceFlowWrapper::ReadRestart



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Distructor (public)                                    ismail 04/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

FLD::UTILS::TotalTractionCorrector::~TotalTractionCorrector()
{
  return;
}

