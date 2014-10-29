/*----------------------------------------------------------------------*/
/*!
\file mc_var_thickness_manager.H
\brief class to modify wall thickness of cardiovascular structures on the fly

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 28915276
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "mc_var_thickness_manager.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"
#include "../linalg/linalg_utils.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fs3i/biofilm_fsi_utils.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"


/*----------------------------------------------------------------------*/
/* constructor                                               jb 05/14   */
/*----------------------------------------------------------------------*/
UQ::MCVarThicknessManager::MCVarThicknessManager(Teuchos::RCP<DRT::Discretization> discret, const int my_thickness_field_id):
discret_(discret)
{
  if (not discret_->Filled() || not discret_->HaveDofs())
      dserror("Discretisation is not complete or has no dofs!");

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  //check if quasistatic analysis
  if(sdyn.get<std::string>("DYNAMICTYP")!= "Statics")
    dserror ("Structure with ale only for quasistatic analysis so in new sti so far.");

  // create empty disscretization and add it to the problem
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = Teuchos::rcp(new DRT::Discretization("ale",DRT::Problem::Instance()->GetNPGroup()->LocalComm()));

  aledis->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(aledis)));

  DRT::Problem::Instance()->AddDis("ale",aledis);

  if (!aledis->Filled()) aledis->FillComplete();

  // we use the structure discretization as layout for the ale discretization
  if (discret_->NumGlobalNodes()==0) dserror("Structure discretization is empty!");
  // clone ale mesh from structure discretization
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(discret_,aledis);
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else
    dserror("Reading an ALE mesh from the input file is not supported for this problem type.");

  // create ale
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale
    = Teuchos::rcp(new ALE::AleBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), DRT::Problem::Instance()->GetDis("ale")));
  ale_ = ale->AleFieldrcp();

  int my_dummy_seed =1 ;
  randomfield_=CreateRandomField(my_thickness_field_id, my_dummy_seed);

  // get the conditions for the current evaluation
  std::vector<DRT::Condition* > uncert_surface;
  discret_->GetCondition("UncertainSurface",uncert_surface);

  // check wether length of condition is one
    if (uncert_surface.size()!=1)
      dserror("Uncertain Surface currently only implemented for 1 condition only");

    my_uncert_nodeids_ = uncert_surface[0]->Nodes();
    my_uncert_nnodes_=(int)my_uncert_nodeids_->size();

    // center of element
    std::vector<double> nodepos(3,0.0);

    // loop nodes in condition
    for (int i=0; i<my_uncert_nnodes_; ++i)
    {
      // do only nodes in my row map
      int nlid = discret_->NodeRowMap()->LID((*my_uncert_nodeids_)[i]);
      if (nlid < 0) continue;
        DRT::Node* actnode = discret_->lRowNode( nlid );

        for(int dim=0; dim<3; ++dim)
        {
          nodepos[dim] = actnode->X()[dim];
        }

        my_uncertain_nodes_org_pos_.insert( std::pair<int,std::vector<double> >((*my_uncert_nodeids_)[i],nodepos) );

    }
    LINALG::GatherAll(my_uncertain_nodes_org_pos_,discret_->Comm());

    // check whether we have at least one ALE DBC specified in the input file
    // get the conditions for the current evaluation
     std::vector<DRT::Condition* > ale_dbc;
     discret_->GetCondition("ALEDirichlet",ale_dbc);
     if (ale_dbc.size()==0)
       dserror("No ALE dirichlet BCs specified! Specify at least one ALE DBC if you want to simulate uncertain surfaces");

     // store initial geometry in dofrowmap layout
     org_geom_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()),true));
     const int nummyrownode = (discret_->NodeRowMap())->NumMyElements();

     // loop over all nodes
     for(int i=0;i<nummyrownode;i++)
     {
       // get current node
       DRT::Node* mynode = discret_->lRowNode(i);

       std::vector<int> globaldofs = discret_->Dof(0,mynode);

       // determine number of space dimensions
       const int numdim = DRT::Problem::Instance()->NDim();

       for (int j=0; j<numdim; j++)
       {
         const int lid = org_geom_->Map().LID(globaldofs[j]);
         if (lid<0)
           dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",org_geom_->Comm().MyPID(),globaldofs[j]);
         (*org_geom_)[lid] = (mynode->X())[j];
       }

     }
     ComputeNormals();

     // input parameters
     const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

     initial_wall_thickness_= mlmcp.get<double>("INITIALTHICKNESS"); ;



}


/*----------------------------------------------------------------------*
 |  Compute new realizations of random field and set thickness jb 07/14  |
 *----------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::SetUpThickness(unsigned int myseed, double para_cont_parameter, bool reuse_rf)
{
  if(!reuse_rf)
    CreateNewRealizationOfRandomField(myseed);
  SetThickness(para_cont_parameter);
}


/*----------------------------------------------------------------------*/
/* Compute thickness and adjust geometry                      jb 05/14  */
/*----------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::SetThickness(double para_cont_parameter)
{
  double local_wall_thick;
  // clear map
  my_uncertain_nodes_pos_.clear();
  std::vector<double> old_node_pos(3,0.0);
  std::vector<double> new_node_pos(3,0.0);

  // loop nodes in condition
   for (int i=0; i<my_uncert_nnodes_; ++i)
   {
     // do only nodes in my row map
     int nlid = discret_->NodeRowMap()->LID((*my_uncert_nodeids_)[i]);
     if (nlid < 0) continue;
       std::vector<double> node_normal = my_uncertain_nodes_normals_.at((*my_uncert_nodeids_)[i]);
       old_node_pos=my_uncertain_nodes_org_pos_.at((*my_uncert_nodeids_)[i]);
       local_wall_thick=randomfield_->EvalFieldAtLocation(old_node_pos,para_cont_parameter,false,false)-initial_wall_thickness_;
       //IO::cout << "HACK " << IO::endl;
       //l/ocal_wall_thick=20;
       // deal with circular quasi 2D field
       if (randomfield_->Dimension()==2)
       {
         double radius =25.0;
         //special hack here assuming circular geometry with r=25 mm
            double phi= acos(old_node_pos[0]/radius);
            //compute x coord
            old_node_pos.at(0)=phi*radius;
            old_node_pos.at(1)=old_node_pos.at(2);
           // ele_center_temp.push_back(ele_center[2]);
            //ele_center_temp.push_back(ele_center[2]);
       }
       //if(old_node_pos.at(2)>-453)
       {
        if(local_wall_thick>-initial_wall_thickness_)
        {
          new_node_pos.at(0)=node_normal.at(0)*(local_wall_thick);
          new_node_pos.at(1)=node_normal.at(1)*(local_wall_thick);
          new_node_pos.at(2)=node_normal.at(2)*(local_wall_thick);
        }
        else
        {
          dserror("The wall thickness you're trying to obtain is negative. Please adjust your input parameters");
          new_node_pos.at(0)=0.0;
          new_node_pos.at(1)=0.0;
          new_node_pos.at(2)=0.0;
        }
        my_uncertain_nodes_pos_.insert(std::pair<int,std::vector<double> >((*my_uncert_nodeids_)[i],new_node_pos) );
       }
   }

   LINALG::GatherAll(my_uncertain_nodes_pos_,discret_->Comm());
   ALEStep(my_uncertain_nodes_pos_);

}

/*----------------------------------------------------------------------*
 |  Create Random field based on input data                   jb 07/14  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<UQ::RandomField> UQ::MCVarThicknessManager::CreateRandomField(int random_field_id, unsigned int myseed)
{
  const Teuchos::ParameterList& rfp = DRT::Problem::Instance()->RandomFieldParams(random_field_id);
  // before calling the constructor make a quick safety check whether this random field was activated in
  // the input file or if the section contains only the default parameters
  bool active = DRT::INPUT::IntegralValue<int>(rfp ,"ACTIVE");
  if (!active)
    dserror("Trying to setup random field that is not active");

  // call constructor based on type of random field
  INPAR::MLMC::CalcMethod calcm = DRT::INPUT::IntegralValue<INPAR::MLMC::CalcMethod>(rfp,"CALC_METHOD");
  Teuchos::RCP<RandomField> rf = Teuchos::null;

  switch(calcm)
  {
    case INPAR::MLMC::calc_m_fft:
      rf = Teuchos::rcp(new RandomFieldSpectral(myseed,discret_,rfp ));
      break;
    case INPAR::MLMC::calc_m_cos:
      rf = Teuchos::rcp(new RandomFieldSpectral(myseed,discret_,rfp ));
      break;
    case INPAR::MLMC::calc_m_fourier:
      rf = Teuchos::rcp(new RandomFieldFourier(myseed,discret_,rfp));
      break;
    default:
      dserror("Unknown simulation method for RF choose fft or cos or fourier");
      break;
  }
  return rf;
}

/*---------------------------------------------------------------------------------------------------*
 * Set ALE displacement conditions to move the uncertain surface
 *---------------------------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::ModifyConditions( const std::map<int, std::vector<double> >& ale_bc_nodes)
{
  AddConditions( ale_->Discretization(), ale_bc_nodes );
  ale_->Discretization()->FillComplete();

  // just adding the Conditions alone do not get the work done
  // we need to build the Dirichlet boundary condition maps again for these conditions to be implemented
  ale_->SetupDBCMapEx( false );
}

/*------------------------------------------------------------------------------------*
 * Add the given node ids and Displacement conditions to discretization    sudhakar 06/14
 * Type of the condition added is PointDirichlet
 * For each node, a separate condition is added. This is helpful to identify
 * these conditions later and delete them at appropriate time
 *------------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::AddConditions( Teuchos::RCP<DRT::Discretization> discret,
                                       const std::map<int, std::vector<double> >& ale_bc_nodes )
{
  std::map<int, std::vector<double> >::const_iterator iter;
  int id = 101;
  for( iter = ale_bc_nodes.begin(); iter != ale_bc_nodes.end(); iter++ )
  {
    Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp( new DRT::Condition( id++, DRT::Condition::PointDirichlet, false, DRT::Condition::Point ) );
    std::vector<int> onoff(3,1);

    //cond->SetConditionType(  DRT::Condition::PointDirichlet );
    cond->Add( "Node Ids", iter->first );
    cond->Add("onoff",onoff);
    cond->Add("val",iter->second);

    discret->SetCondition( "Dirichlet", cond );
  }
}


/*-----------------------------------------------------------------------------------------*
 * Delete all Dirichlet conditions that has only one associated node
 * These are the conditions that we added to move the uncertain surface around
 *-----------------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::DeleteConditions( Teuchos::RCP<DRT::Discretization> discret )
{
  std::vector<std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator> del;

  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator conit;
  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >& allcondn = discret->GetAllConditions();
  for( conit = allcondn.begin(); conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<DRT::Condition> cond = conit->second;
    if( cond->Nodes()->size() == 1 )
      del.push_back( conit );
  }

  for( unsigned i=0; i< del.size(); i++ )
  {
    conit = del[i];
    allcondn.erase( conit );
  }
}

/*----------------------------------------------------------------------------------*
 * Perform all operations of ALE step
 *----------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::ALEStep( const std::map<int, std::vector<double> >& ale_bc_nodes)
{
  ModifyConditions( ale_bc_nodes );

  ALESolve();

}

/*----------------------------------------------------------------------------------*
 * Build ALE system matrix and solve the system
 *----------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::ALESolve()
{
  ale_->BuildSystemMatrix();

  ale_->Solve();

  Teuchos::RCP<const Epetra_Vector> ale_final = ale_->Dispnp();

  // convert displacement in ALE field into structural field
  Teuchos::RCP<Epetra_Vector> ale_str = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap())));
  for( int i=0; i<ale_str->MyLength(); i++ )
    (*ale_str)[i] = (*ale_final)[i];

  FS3I::Biofilm::UTILS::updateMaterialConfigWithALE_Disp( discret_, ale_str );
  //ale_->Discretization()->ClearDiscret();
  ale_->Reset();

  DeleteConditions( ale_->Discretization() );
  ale_->Discretization()->FillComplete();

  discret_->FillComplete();

}

/*---------------------------------------------------------------------------------*
 |  Compute new realizations of random field                 jb 07/14              |
 *---------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::CreateNewRealizationOfRandomField(unsigned int myseed)
{
    randomfield_->CreateNewSample(myseed+174368);
}

/*---------------------------------------------------------------------------------*
 | Compute Surface normal for all node in uncertain surface condition   jb 07/14   |
 *---------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::ComputeNormals()
{
  ADAPTER::CouplingMortar mytemp;
  mytemp.SetupForUQAbuseNormalCalculation(discret_,discret_->Comm());
  mytemp.Interface()->EvaluateNodalNormals(my_uncertain_nodes_normals_);
  LINALG::GatherAll(my_uncertain_nodes_normals_, discret_->Comm());
}

/*-----------------------------------------------------------------------------------*
 * Reset geometry to its original state                                       jb 07/14
 *----------------------------------------------------------------------------------*/
void UQ::MCVarThicknessManager::ResetGeometry()
{
  const int numnode = (discret_->NodeColMap())->NumMyElements();

  //Create Vector which holds all col-displacments of processor
  Teuchos::RCP<Epetra_Vector> coldisp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofColMap())));

  //Export row-displacments to col-displacements
  LINALG::Export(*org_geom_, *coldisp);

  const Epetra_Vector& gvector =*coldisp;

  // loop over all nodes
  for (int index = 0; index < numnode; ++index)
  {
    // get current node
    DRT::Node* mynode = discret_->lColNode(index);

    std::vector<int> globaldofs = discret_->Dof(0,mynode);
    std::vector<double> nvector(globaldofs.size());

    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    for (int i=0; i<numdim; ++i)
    {
      const int lid = gvector.Map().LID(globaldofs[i]);

      if (lid<0)
        dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",gvector.Comm().MyPID(),globaldofs[i]);
      nvector[i] += gvector[lid];
    }

    mynode->SetPos(nvector);
  }

  return;
}

/*-----------------------------------------------------------------------------------*
 * Evaluate underlying random field at a specific location                       jb 07/14
 *----------------------------------------------------------------------------------*/
double UQ::MCVarThicknessManager::EvalThicknessAtLocation(std::vector<double> myloc,double para_cont_parameter)
{
  return randomfield_->EvalFieldAtLocation(myloc,para_cont_parameter,false,false)-initial_wall_thickness_;
}
