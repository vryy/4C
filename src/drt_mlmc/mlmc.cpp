/*!----------------------------------------------------------------------
\file  mlmc.cpp
\brief Class for performing Multi Level Monte Carlo (MLMC)analysis of structure


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef HAVE_FFTW

#include "../drt_adapter/ad_str_structure.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>
#include "mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../global_full/global_inp_control.H"
#include "../drt_io/io_hdf.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/matpar_bundle.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "mc_mat_par_manager.H"
#include "mc_var_thickness_manager.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"


#include "../drt_comm/comm_utils.H"
//for file output
#include <fstream>
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_element.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::UQ::MLMC::MLMC(Teuchos::RCP<DRT::Discretization> dis)
  : discret_(dis)
{

  // get coarse and fine discretizations
  actdis_coarse_ = DRT::Problem::Instance(0)->GetDis("structure");

  // set degrees of freedom in the discretization
  if (not actdis_coarse_->Filled() || not actdis_coarse_->HaveDofs()) actdis_coarse_->FillComplete();

  filename_ = DRT::Problem::Instance()->OutputControlFile()->FileName();

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // get number of timesteps
  tsteps_ = sdyn.get<int>("NUMSTEP");
  // input parameters multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // Get number of Newton iterations
  num_newton_it_ = mlmcp.get<int>("ITENODEINELE");
  // Get convergence tolerance
  convtol_    = mlmcp.get<double>("CONVTOL");

  // calculate difference to lower level yes/no
  calc_diff_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"DIFF_TO_LOWER_LEVEL");

  // prolongate results yes/no
  prolongate_res_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"PROLONGATERES");

  // value for stochmat blending in parameter continuation case
  cont_num_maxtrials_ = mlmcp.get<int>("CONTNUMMAXTRIALS");

  // get starting random seed
  start_random_seed_ = mlmcp.get<int>("INITRANDOMSEED");

  // get name of lower level outputfiles
  filename_lower_level_ =  mlmcp.get<std::string>("OUTPUT_FILE_OF_LOWER_LEVEL");

  // get numerb of current level
  num_level_ = mlmcp.get<int>("LEVELNUMBER");

  //write statistics every write_stat_ steps
  write_stats_ = mlmcp.get<int>("WRITESTATS");

  stoch_wall_thickness_=DRT::INPUT::IntegralValue<int>(mlmcp ,"RANDOMGEOMETRY");

  double word;
  std::istringstream bsdampingstream(Teuchos::getNumericStringParameter(mlmcp,"OUTPUT_ELEMENT_IDS"));
  while (bsdampingstream >> word)
	  AllMyOutputEleIds_.push_back(int (word));

  if(AllMyOutputEleIds_.front()== -1)
	  IO::cout << RED_LIGHT "No elements specified for output " END_COLOR << IO::endl;

  InEleRange_ = 1.0 + 10e-3;


  // controlling parameter
  start_run_ = mlmcp.get<int>("START_RUN");
  int numruns = mlmcp.get<int>("NUMRUNS");
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  int NNestedGroups = problem->GetNPGroup()->NumGroups();
  int i = problem->GetNPGroup()->GroupId();

  numruns_pergroup_= int(ceil(numruns/NNestedGroups));
  start_run_  += (i)*numruns_pergroup_;

  numb_run_ =  start_run_;    // counter of how many runs were made monte carlo

  // init number of wall elements
  tot_num_wall_elements_ = 0;
  SetupEvalDisAtEleCenters(AllMyOutputEleIds_);
  SetupEvalPeakWallStress();


  reduced_output_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"REDUCED_OUTPUT");

  write_rv_to_file_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"WRITE_RV_TO_FILE");

  // meshfile name to be written to controlfile in prolongated results
  std::stringstream meshfilename_helper1;
  std::string meshfilename_helper2;
  meshfilename_helper1 << filename_ << "_prolongated_run_" << start_run_ ;
  meshfilename_helper2 = meshfilename_helper1.str();
  // strip path from name
  std::string::size_type pos = meshfilename_helper2.find_last_of('/');
  if (pos==std::string::npos)
	  meshfilename_ = meshfilename_helper2;
  else
	  meshfilename_ = meshfilename_helper2.substr(pos+1);


  // parameter continuation parameter
  // init variables to store restart information for parameter continuation
  cont_step_ = Teuchos::rcp(new int);
  cont_time_  = Teuchos::rcp(new double);
  cont_disn_ = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  cont_disn_init_ = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  cont_veln_ = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  cont_accn_ = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  // this should also have row map layout
  // we need to of those because we need to store and keep eleementdata of initial run with constant beta
  cont_elementdata_init_ = Teuchos::rcp(new std::vector<char>);
  cont_elementdata_ = Teuchos::rcp(new std::vector<char>);
  cont_nodedata_init_ = Teuchos::rcp(new std::vector<char>);
  cont_nodedata_ = Teuchos::rcp(new std::vector<char>);


  //set up managers to create stochasticity
  my_matpar_manager_ = Teuchos::rcp(new STR::UQ::MCMatParManager(actdis_coarse_));
  // init field with some seed
  my_matpar_manager_->SetupRandomFields(2);

  my_matpar_manager_->NumPhysStochParams();

  if(stoch_wall_thickness_)
  {
    my_thickness_manager_ = Teuchos::rcp(new MCVarThicknessManager(actdis_coarse_,my_matpar_manager_->NumPhysStochParams()+1));
  }
  else
    my_thickness_manager_=Teuchos::null;

    //my_thickness_manager_->SetUpThickness(2,1.6,false);

  //init stuff that is only needed when we want to prolongate the results to a finer mesh,
  // and hence have a fine discretization
  if(prolongate_res_ )
  {
    // Get finest Grid problem instance

    actdis_fine_ = DRT::Problem::Instance(1)->GetDis("structure");
    // set degrees of freedom in the discretization
    if (not actdis_fine_->Filled()) actdis_fine_->FillComplete();
    // Get coarse Grid problem instance

    output_control_fine_ = Teuchos::rcp(new IO::OutputControl(actdis_fine_->Comm(), "structure", "Polynomial", filename_, filename_, 3, 0, 20, 1));
    output_fine_ = actdis_fine_->Writer();
    output_fine_->SetOutput(output_control_fine_);

    // init vectors to store mean stresses and displacements
    mean_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    mean_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    mean_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    // init vectors to store standard dev of  stresses and displacements
    var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    // init vectors to calc standard dev of  stresses and displacements
    delta_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    delta_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    delta_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    m2_var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    m2_var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    m2_var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    m2_helper_var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    m2_helper_var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    m2_helper_var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));


    // same vectors for difference between two levels
    // init vectors to store mean stresses and displacements
    diff_mean_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    diff_mean_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_mean_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    // init vectors to store standard dev of  stresses and displacements
    diff_var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    diff_var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    // init vectors to calc standard dev of  stresses and displacements
    diff_delta_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    diff_delta_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_delta_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_m2_var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    diff_m2_var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_m2_var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_m2_helper_var_disp_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    diff_m2_helper_var_stress_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    diff_m2_helper_var_strain_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));


    disp_lower_level_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));
    stress_lower_level_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
    strain_lower_level_ = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
  }
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::UQ::MLMC::Integrate()
{

  // init vector to store displacemnet

  //Teuchos::RCP<Epetra_Vector> dis_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  Teuchos::RCP<const Epetra_Vector> dis_coarse = Teuchos::null;

  const int myrank = discret_->Comm().MyPID();

  //measure time
  Epetra_Time timer(discret_->Comm());
  //double t1 = timer.ElapsedTime();
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  //int numruns = mlmcp.get<int>("NUMRUNS")+start_run_;
  // nested par hack
  int numruns =numruns_pergroup_+start_run_;


  // get initial random seed from inputfile
  unsigned int random_seed= mlmcp.get<int>("INITRANDOMSEED");

  // init variables to store restart information for parameter continuation
  Teuchos::RCP<int> cont_step = Teuchos::rcp(new int);
  Teuchos::RCP<double> cont_time  = Teuchos::rcp(new double);
  Teuchos::RCP<Epetra_Vector> cont_disn = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  Teuchos::RCP<Epetra_Vector> cont_veln = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  Teuchos::RCP<Epetra_Vector> cont_accn = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  // this should also have row map layout
  Teuchos::RCP<std::vector<char> > cont_elementdata = Teuchos::rcp(new std::vector<char>);

  const double t0 = Teuchos::Time::wallTime();
  do
  {
    if (myrank == 0)
    {
      IO::cout << GREEN_LIGHT "================================================================================" << IO::endl;
      IO::cout << "                            MULTILEVEL MONTE CARLO                              " << IO::endl;
      IO::cout << "                              RUN: " << numb_run_ << "  of  " << numruns  << IO::endl;
      IO::cout << "================================================================================" END_COLOR << IO::endl;
    }

    if (myrank == 0 )
    {
      IO::cout << GREEN_LIGHT " RESET PRESTRESS " END_COLOR << IO::endl;
    }
    ResetPrestress();
    const double t1 = Teuchos::Time::wallTime();
    my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),1.0,false);

    if(stoch_wall_thickness_)
    {
      my_thickness_manager_->SetUpThickness((random_seed+(unsigned int)numb_run_),1.0,false);
    }
    const double t2 = Teuchos::Time::wallTime();
    discret_->Comm().Barrier();

    discret_->Writer()->NewResultFile(filename_,(numb_run_));
    //discret_->Writer()->WriteMesh(0, 0.01);


    // get input lists
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    //store time
    double t3 =0;
    double t4 =0;
    // major switch to different time integrators
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
    {
      case INPAR::STR::dyna_statics:
      case INPAR::STR::dyna_genalpha:
      case INPAR::STR::dyna_onesteptheta:
      case INPAR::STR::dyna_gemm:
      case INPAR::STR::dyna_expleuler:
      case INPAR::STR::dyna_centrdiff:
      case INPAR::STR::dyna_ab2:
      case INPAR::STR::dyna_euma:
      case INPAR::STR::dyna_euimsto:
      {
        // instead of calling dyn_nlnstructural_drt(); here build the adapter here so that we have acces to the results
        // What follows is basicaly a copy of whats usually in dyn_nlnstructural_drt();
        // access the structural discretization
        Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
        // create an adapterbase and adapter
        ADAPTER::StructureBaseAlgorithm adapterbase(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
        ADAPTER::Structure& structadaptor = adapterbase.StructureField();

        // do restart
        const int restart = DRT::Problem::Instance()->Restart();
        if (restart)
        {
          structadaptor.ReadRestart(restart);
        }

        t3 = Teuchos::Time::wallTime();
        // do we want parameter continuation (for now we use the reset presstress flag)
        structadaptor.Integrate();

        t4 = Teuchos::Time::wallTime();

        dis_coarse= Teuchos::rcp(new const Epetra_Vector(*(structadaptor.Dispn())));
        // test results
        DRT::Problem::Instance()->AddFieldTest(structadaptor.CreateFieldTest());
        DRT::Problem::Instance()->TestAll(structadaptor.DofRowMap()->Comm());

        if(!reduced_output_)
        {
          Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(structadaptor.DofRowMap()->Comm());
          Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
        }

        // time to go home...
      }
        break;
      default:
        dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
        break;
    }

    Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_disp =Teuchos::rcp(new std::vector <std::vector<double > >);
    Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_stresses =Teuchos::rcp(new std::vector <std::vector<double > >);
    Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_strains =Teuchos::rcp(new std::vector <std::vector<double > >);
    Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_mat_params =Teuchos::rcp(new std::vector <std::vector<double > >);

    EvalDisAtEleCenters(dis_coarse,INPAR::STR::stress_2pk,INPAR::STR::strain_gl,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
    ExportEleDataAndWriteToFile(OutputMap_,INPAR::STR::stress_2pk,INPAR::STR::strain_gl,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
    // clear vectors as they are filled via pushback
    my_output_elements_c_disp->clear();
    my_output_elements_c_stresses->clear();
    my_output_elements_c_strains->clear();
    my_output_elements_mat_params->clear();
    EvalDisAtEleCenters(dis_coarse,INPAR::STR::stress_cauchy  ,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
    ExportEleDataAndWriteToFile(OutputMap_,INPAR::STR::stress_cauchy,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);

    // Evaluate Peak stress and strain vallues
    EvalPeakWallStress(dis_coarse, INPAR::STR::stress_cauchy,INPAR::STR::strain_ea );


    // write random variables to file for regression
    if(write_rv_to_file_)
    my_matpar_manager_->WriteRandomVariablesToFile(filename_,numb_run_);

    if (numb_run_-start_run_== 0 &&  prolongate_res_)
    {
     // not parallel
      SetupProlongatorParallel();
    }
    if (calc_diff_)
    {
      ReadResultsFromLowerLevel();
    }
    if ( prolongate_res_)
    {
      ProlongateResults();
      // write statoutput evey now and then
      if(numb_run_% write_stats_ == 0)
        WriteStatOutput();
    }

    // reset geometry to initial read in geometry
    if(stoch_wall_thickness_)
    {
      my_thickness_manager_->ResetGeometry();
    }

    numb_run_++;
    const double t5 = Teuchos::Time::wallTime();
    if(!reduced_output_)
    {
      // plot some time measurements
      if(!discret_->Comm().MyPID())
          {
            IO::cout<<"\n=================Time  Measurement================"<<IO::endl;
            IO::cout<<"Setup Stochastic Material:\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<IO::endl;
            IO::cout<<"Forward Solve  :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<IO::endl;
            IO::cout<<"Total Wall Time After " << numb_run_ << " runs:\t"<<std::setprecision(4)<<t5-t0<<"\ts"<<IO::endl;
            IO::cout<<"=================================================="<<IO::endl;
          }
    }


    } while (numb_run_< numruns);
  const double t6 = Teuchos::Time::wallTime();
  IO::cout<<"\n=================Time  Measurement================"<<IO::endl;
  IO::cout<<"Total runtime:\t"<<std::setprecision(4)<<t6-t0<<"\ts"<<IO::endl;
  IO::cout<<"=================================================="<<IO::endl;
  if( prolongate_res_ )
  {
    WriteStatOutput();
  }
  return;
}
/*----------------------------------------------------------------------*/
/* Nice and clean version with for not resetting prestress*/
void STR::UQ::MLMC::IntegrateNoReset()
{
  // init vector to store displacemnet
  Teuchos::RCP<const Epetra_Vector> dis_coarse = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> dis_coarse2 = Teuchos::null;

  const int myrank = discret_->Comm().MyPID();

  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // vector to store parameter continuation info in
  Teuchos::RCP<std::vector <int> > paramcont_info  = Teuchos::rcp(new std::vector <int>(3,0));

  // nested par hack
  int numruns =numruns_pergroup_+start_run_;

  // get initial random seed from inputfile
  unsigned int random_seed= mlmcp.get<int>("INITRANDOMSEED");

  // get initial number of continuation steps from input
  unsigned int num_cont_steps= mlmcp.get<int>("NUMCONTSTEPS");

  //double time_step_size = 0.0 ;


  //store time init with zero just to be sure
  t0_ = 0.0;
  t1_ = 0.0;
  t2_ = 0.0;
  t3_ = 0.0;
  t4_ = 0.0;
  t0_ = Teuchos::Time::wallTime();
  do
  {
    if (myrank == 0)
    {
      IO::cout << GREEN_LIGHT "================================================================================" << IO::endl;
      IO::cout << "                            MULTILEVEL MONTE CARLO                              " << IO::endl;
      IO::cout << "                          RUN: " << numb_run_ << "  of  " << numruns  << IO::endl;
      IO::cout << "================================================================================" END_COLOR << IO::endl;
    }
    ResetPrestress();
    discret_->Comm().Barrier();

    discret_->Writer()->NewResultFile(filename_,(numb_run_));
    //discret_->Writer()->WriteMesh(0, 0.0);


    // major switch to different time integrators
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
    {
      case INPAR::STR::dyna_statics:
      case INPAR::STR::dyna_genalpha:
      case INPAR::STR::dyna_onesteptheta:
      case INPAR::STR::dyna_gemm:
      case INPAR::STR::dyna_expleuler:
      case INPAR::STR::dyna_centrdiff:
      case INPAR::STR::dyna_ab2:
      case INPAR::STR::dyna_euma:
      case INPAR::STR::dyna_euimsto:
      {
        // instead of calling dyn_nlnstructural_drt(); here build the adapter here so that we have acces to the results
        // What follows is basicaly a copy of whats usually in dyn_nlnstructural_drt();
        // create an adapterbase and adapter
        // access the structural discretization
        Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
        ADAPTER::StructureBaseAlgorithm adapterbase(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
        ADAPTER::Structure& structadaptor = adapterbase.StructureField();


        // do restart
        const int restart = DRT::Problem::Instance()->Restart();
        if (restart)
        {
          structadaptor.ReadRestart(restart);
        }
        /// Try some nasty stuff
        if (numb_run_-start_run_== 0)
        {
        
          // for first run do normal integration
          my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),0.0,false);
          structadaptor.Integrate();
          // get all information
          structadaptor.GetRestartData(cont_step_,cont_time_,cont_disn_init_,cont_veln_,cont_accn_,cont_elementdata_init_,cont_nodedata_init_);
          // Get timestep size
          //time_step_size=structadaptor.Dt();
        }
        else
        {
          // compute solution with parameter continuation
          int error = ParameterContinuation( num_cont_steps, random_seed,false, structadaptor);
          // if not converged
          int num_trials = 1;
          int num_cont_steps_new = num_cont_steps;
          while( error && num_trials < cont_num_maxtrials_ )
          {
            // double number of cont steps
            num_cont_steps_new = num_cont_steps_new*2;
            error = ParameterContinuation( num_cont_steps_new, random_seed, true, structadaptor);
            num_trials++;
            // check for error
            if(num_trials==cont_num_maxtrials_  && error)
              dserror("your problem is stupid, go away");
          }
          // how many step did I plan to use
          paramcont_info->at(0)=num_cont_steps;
          paramcont_info->at(1)=num_cont_steps_new;
          paramcont_info->at(2)=num_trials;


        }


        // test results
        DRT::Problem::Instance()->AddFieldTest(structadaptor.CreateFieldTest());
        DRT::Problem::Instance()->TestAll(structadaptor.DofRowMap()->Comm());

        if(!reduced_output_)
        {
          Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(structadaptor.DofRowMap()->Comm());
          Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
        }

        // time to go home...
      }
      break;
      default:
        dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
        break;
      } // eof switch case

    // check wether we have elements to compute the stresses at
    if(AllMyOutputEleIds_.size()>0)
    {
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_disp =Teuchos::rcp(new std::vector <std::vector<double > >);
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_stresses =Teuchos::rcp(new std::vector <std::vector<double > >);
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_strains =Teuchos::rcp(new std::vector <std::vector<double > >);
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_mat_params =Teuchos::rcp(new std::vector <std::vector<double > >);

      EvalDisAtEleCenters(cont_disn_,INPAR::STR::stress_2pk,INPAR::STR::strain_gl,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
      ExportEleDataAndWriteToFile(OutputMap_,INPAR::STR::stress_2pk,INPAR::STR::strain_gl,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
      EvalDisAtEleCenters(cont_disn_,INPAR::STR::stress_cauchy  ,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
      ExportEleDataAndWriteToFile(OutputMap_,INPAR::STR::stress_cauchy,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);

      // Evaluate Peak stress and strain vallues
      EvalPeakWallStress(cont_disn_, INPAR::STR::stress_cauchy,INPAR::STR::strain_ea );
      // write parameter continuation data to file
      WriteParamContInfoToFile(paramcont_info);
    }

    // write random variables to file for regression
    if(write_rv_to_file_)
    my_matpar_manager_->WriteRandomVariablesToFile(filename_,numb_run_);

    t5_ = Teuchos::Time::wallTime();
    numb_run_++;
    if(!reduced_output_)
    {
      // plot some time measurements
      if(!discret_->Comm().MyPID())
      {
        IO::cout<<"\n=================Time  Measurement================"<<IO::endl;
        IO::cout<<"Setup Stochastic Material:\t"<<std::setprecision(4)<<t3_-t2_<<"\ts"<<IO::endl;
        IO::cout<<"Forward Solve  :\t"<<std::setprecision(4)<<t4_-t1_-(t3_-t2_)<<"\ts"<<IO::endl;
        IO::cout<<"Total Wall Time After " << numb_run_ << " runs:\t"<<std::setprecision(4)<<t5_-t0_<<"\ts"<<IO::endl;
        IO::cout<<"=================================================="<<IO::endl;
      }
    }

  }while (numb_run_< numruns); // eof outer MC loop
 t6_ = Teuchos::Time::wallTime();
  IO::cout<<"\n=================Time  Measurement================"<<IO::endl;
  IO::cout<<"Total runtime:\t"<<std::setprecision(4)<<t6_-t0_<<"\ts"<<IO::endl;
  IO::cout<<"=================================================="<<IO::endl;
  return;
}

//---------------------------------------------------------------------------------------------
int STR::UQ::MLMC::ParameterContinuation(unsigned int num_cont_steps, unsigned int random_seed, bool re_use_rf, ADAPTER::Structure&  structadaptor)
{
  int error = 0;
  t1_ = Teuchos::Time::wallTime();
  // init some variables
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  *cont_elementdata_=*cont_elementdata_init_;
  *cont_nodedata_=*cont_nodedata_init_;
  *cont_disn_=*cont_disn_init_;

  // Get timestep size
  const double time_step_size=structadaptor.Dt();

  // for loop to encapsulate parameter continuation scheme
  for(unsigned int k =0 ; k<num_cont_steps ; k++)
  {
    double gamma=1.0/num_cont_steps+(k*(1.0/num_cont_steps));
    // set new outputfile
    std::stringstream filename_helper1;
    std::string filename_helper2;
    std::string filename_paracont;
    filename_helper1 << filename_ << "_gamma_" << gamma ;
    filename_helper2 = filename_helper1.str();
    // strip path from name
    std::string::size_type pos = filename_helper2.find_last_of('/');
    if (pos==std::string::npos)
      filename_paracont = filename_helper2;
    else
      filename_paracont = filename_helper2.substr(pos+1);

    IO::cout << "filename_paracont " << filename_paracont << IO::endl;
    IO::cout << "filename_ " << filename_ << IO::endl;

    discret_->Writer()->NewResultFile(filename_paracont,(numb_run_));


    IO::cout << "Gamma is  " << gamma << IO::endl;
    // if we have prestress we need possibly a two step process
    // get prestress type
    INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");

     if(pstype==INPAR::STR::prestress_mulf)
     {
       // get prestress time
       double pstime = sdyn.get<double>("PRESTRESSTIME");
       double endtime;

       // compute the number of steps we did in prestressing mode
       int num_prestress_steps = (int)(pstime/time_step_size);
       // store maxtime
       endtime= structadaptor.GetTimeEnd();
       // alter timemax to compute only one timestep
       structadaptor.SetTimeEnd(pstime);
       // set restart without displacement (use acc instead which are zero in statics)
       structadaptor.SetRestart(num_prestress_steps-1,(num_prestress_steps-1)*time_step_size,cont_accn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
       // set new material parameters
       // Setup Material Parameters in each element based on deterministic value
       //measure time only if we really compute the random field (k=0)
       if(!re_use_rf)
       {
         t2_ = Teuchos::Time::wallTime();
         my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),gamma,re_use_rf);
         t3_ = Teuchos::Time::wallTime();
         // we only want to setupt once hence
         re_use_rf=true;
       }
       else if(re_use_rf)
       {
         my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),gamma,re_use_rf);
       }
       error=structadaptor.Integrate();
       if(error)
       {
         // reset endtime
         structadaptor.SetTimeEnd(endtime);
         return error;
       }

       // create dummy rcp to throw away
       Teuchos::RCP<int> garbage1 = Teuchos::rcp(new int);
       Teuchos::RCP<double> garbage2  = Teuchos::rcp(new double);
       // get the newly computed prestress deformation gradient
       // use acc for disp again which will be overritten by the accc directly after
       structadaptor.GetRestartData(garbage1,garbage2,cont_accn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
       // write new prestress data together with old displament data to discretization
       structadaptor.SetRestart((*(cont_step_)-1),(*cont_time_)-time_step_size,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
       // reset old maxtime
       structadaptor.SetTimeEnd(endtime);
       my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),gamma,true);
       error = structadaptor.Integrate();
       if(error)
         return error;
       structadaptor.GetRestartData(garbage1,garbage2,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
     }
     else // no prestress
     {
       // set prestress state and time, disp acc vel are al set to zero
       structadaptor.SetRestart((*(cont_step_)-1),(*cont_time_)-time_step_size,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
       //measure time only if we really compute the random field (k=0)
       if(!k)
         t2_ = Teuchos::Time::wallTime();
       my_matpar_manager_->SetUpStochMats((random_seed+(unsigned int)numb_run_),gamma,(bool)k);
       if(!k)
         t3_ = Teuchos::Time::wallTime();
       error = structadaptor.Integrate();
       if(error)
         return error;
       structadaptor.GetRestartData(cont_step_,cont_time_,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_,cont_nodedata_);
     }

  }// eof parameter continuation scheme
    t4_ = Teuchos::Time::wallTime();
    return 0;
}


// Setup Material Parameters in each element based on deterministic value
void STR::UQ::MLMC::SetupStochMatDet(double value)
{
  // flag have init stochmat??
  int stochmat_flag=0;
  // Get parameters from stochastic matlaw
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if (myrank == 0)
    IO::cout << "No. material laws considered: " << (int)mats.size() << IO::endl;
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); curr++)
  {
    const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
    case INPAR::MAT::m_aaaneohooke_stopro:
    {
      dserror("aaaneohook_stopro no longer supported");
    }
    break;
    case INPAR::MAT::m_aaaneohooke:
    {
      stochmat_flag=1;
      MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
      if (!params) dserror("Cannot cast material parameters");
      // get epetra vector in elecol layout
      Teuchos::RCP<Epetra_Vector> my_param = Teuchos::rcp(new Epetra_Vector(*discret_->ElementColMap(),true));
      my_param->PutScalar(value);
      // set values
      params->SetParameter(params->beta ,my_param);
    }
    break;
    default:
    {
      IO::cout << "MAT CURR " << actmat->Type() << "not stochastic" << IO::endl;
      break;
    }
    }
  }
  // EOF loop over mats
  if(!stochmat_flag)// ignore unknown materials ?
  {
    dserror("No stochastic material supplied");
  }

}

void STR::UQ::MLMC::ResetPrestress()
{
  // Reset Presstress possibly still present in Discretization
  // Get prestress parameter
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // get prestress type
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    switch(pstype)
    {
    case INPAR::STR::prestress_none:
    {
      IO::cout << "nothing to do no prestressing used " << IO::endl;
    }
    break;
    case INPAR::STR::prestress_mulf:
    {
      Teuchos::ParameterList p;
      // action for elements
      p.set("action","calc_struct_reset_all");
      discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
     }
    break;
    case INPAR::STR::prestress_id:
    {
      dserror("MLMC and pressstressing with ID do not go great together");
    }
    break;
    default:
      dserror("Unknown type of prestressing");
      break;
    }
}


void STR::UQ::MLMC::HelperFunctionOutput(Teuchos::RCP< Epetra_MultiVector> stress,Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp)
{
  // assamble name for outputfil
  std::stringstream outputfile;
  outputfile << filename_ << "_statistics_output_" << start_run_ << ".txt";
  std::string name = outputfile.str();;
  /// file to write output
  std::ofstream File;
  if (numb_run_ == 0 || numb_run_ == start_run_)
  {
    File.open(name.c_str(),std::ios::out);
    File << "run id   "<< "xdisp node 24 " << "S_xx node 436 "<< std::endl;
    File.close();
  }
  // reopen in append mode
  File.open(name.c_str(),std::ios::app);
  File << numb_run_ << "    "<< (*disp)[0][72]<< "    " << (*stress)[0][435] << std::endl;
  File.close();
}

void STR::UQ::MLMC::EvalDisAtNodes(Teuchos::RCP<const Epetra_Vector> disp )
{

  const int myrank = actdis_coarse_->Comm().MyPID();

  // build map that lives only on proc 0
  // nodes for fine discretization
   int node[5] = {1528, 3905, 7864, 10832, 13720};
   int dofs[15]={4584, 4585, 4586,11715, 11716, 11717,23592,23593, 23594,32496, 32497, 32498,41160, 41161, 41162};
   // nodes for coarse sicretization
   //int node[5] = {112, 257, 526, 728, 910};
   // set entries to zero for testing


  int numglobalelements =5;
  int numglobalelements_dof =15;
  int nummyelements;
  int nummyelements_dof;

  if (actdis_coarse_->Comm().MyPID()==0)
  {
   nummyelements = 5;
   nummyelements_dof = 15;
  }
  else
  {
   nummyelements = 0;
   nummyelements_dof = 0;
  }

  Epetra_Map output_node_map(numglobalelements,nummyelements,&node[0],0,actdis_coarse_->Comm());
  Epetra_Map output_dof_map(numglobalelements_dof,nummyelements_dof,&dofs[0],0,actdis_coarse_->Comm());

  INPAR::STR::StressType iostress =INPAR::STR::stress_2pk; //stress_none;
  INPAR::STR::StrainType iostrain= INPAR::STR::strain_gl; // strain_none;
  Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > plstrain = Teuchos::rcp(new std::vector<char>());

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  p.set("action","calc_struct_stress");
  p.set("stress", stress);
  p.set("plstrain",plstrain);
  p.set("strain", strain);

  p.set<int>("iostress", iostress);
  p.set<int>("iostrain", iostrain);

  Teuchos::RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  Teuchos::RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  // disp is passed to the function no need for reading in results anymore
  actdis_coarse_->SetState("displacement",disp);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alrigth lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");
  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
 p.set("gpstrainmap", gpstrainmap);


  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();


  // st action to calc poststresse
  p.set("action","postprocess_stress");
  // Multivector to store poststresses
  // WE need the collumn Map here because we need all the ghosted nodes to calculate stresses at the nodes
  Teuchos::RCP<Epetra_MultiVector> poststress =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeColMap()),6,true));

  p.set("gpstressmap", gpstressmap);
  p.set("poststress", poststress);
  p.set("stresstype","ndxyz");

  //actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();

  // again for strains
  p.set("action","postprocess_stress");
  // Multivector to store poststrains
  Teuchos::RCP<Epetra_MultiVector> poststrain =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeColMap()),6,true));
  p.set("poststress", poststrain);
  p.set("gpstressmap", gpstrainmap);
  p.set("stresstype","ndxyz");

  actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();

  // also get element stresses
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", gpstressmap);
  Teuchos::RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
  p.set("poststress",elestress);
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }
  IO::cout << "Result stress" << *elestress << IO::endl;
  dserror("stop here");
  // and strains
  p.set("gpstressmap", gpstrainmap);
  Teuchos::RCP<Epetra_MultiVector> elestrains = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
  p.set("poststress",elestrains);
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  // assamble name for outputfile
  std::stringstream outputfile2;
  outputfile2 << filename_ << "_statistics_output_" << start_run_ << ".txt";
  std::string name2 = outputfile2.str();;
  // file to write output

  Teuchos::RCP<Epetra_Vector> output_disp = LINALG::CreateVector(output_dof_map,false);
  Teuchos::RCP<Epetra_MultiVector> output_stress =  Teuchos::rcp(new Epetra_MultiVector(output_node_map,6,true));
  Teuchos::RCP<Epetra_MultiVector> output_strain =  Teuchos::rcp(new Epetra_MultiVector(output_node_map,6,true));
  LINALG::Export(*poststress,*output_stress);
  LINALG::Export(*poststrain,*output_strain);
  LINALG::Export(*disp,*output_disp);


  if(myrank==0)
  {
    std::ofstream File;
    if (numb_run_ == 0 || numb_run_ == start_run_)
    {
      File.open(name2.c_str(),std::ios::out);
      if (File.is_open())
      {
        //if(error == -1)
        //dserror("Unable to open Statistics output Filename");
        File << "run id   "<< "disp x  disp y disp z  node   " << node[0] << " "
               << "disp x  disp y disp z  node   " << node[1]<< " "
               << "disp x  disp y disp z  node   " << node[2]<< " "
               << "disp x  disp y disp z  node   " << node[3]<< " "
               << "disp x  disp y disp z  node   " << node[4]<< " "
               // stresses
               << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[0] << " "
               << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[1] << " "
               << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[2] << " "
               << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[3] << " "
               << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[4] << " "
              // << "stress xx  stress yy stress zz stress xy stress yz stress xz node   " << node[5] << " "
               // strains
               << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[0] << " "
               << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[1] << " "
               << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[2] << " "
               << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[3] << " "
               << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[4] << " "
              // << "strain xx  strain yy strain zz strain xy strain yz strain xz node   " << node[5] << " "
               << std::endl;
        File.close();
      }
      //
      else
      {
      dserror("Unable to open statistics output file");
      }
    }
    // reopen in append mode
    File.open(name2.c_str(),std::ios::app);
    File << numb_run_ ;
    for(int i=0;i<15;i++)
    {
      File <<  "  " << (*output_disp)[i];
    }
    for (int i=0;i<5;i++)
    {
      for(int j=0;j<6;j++)
      {
        File << " " << (*output_stress)[j][i];
      }
    }
    for (int i=0;i<5;i++)
    {
      for(int j=0;j<6;j++)
      {
        File << " " << (*output_strain)[j][i];
      }
    }
    File << std::endl;
    File.close();
  }
  actdis_coarse_->Comm().Barrier();
}


void STR::UQ::MLMC::SetupEvalDisAtEleCenters(std::vector <int> AllOutputEleIds)
{
  //const int myrank = actdis_coarse_->Comm().MyPID();
  for (unsigned int i=0; i<AllOutputEleIds.size(); i++)
  {
    if(actdis_coarse_->ElementRowMap()->LID(AllOutputEleIds[i]) > -1)
    {
      my_output_elements_.push_back(AllOutputEleIds[i]);
    }
  }

  int NumGlobalMapElements =AllOutputEleIds.size();
  int NumMyElements;

  if (actdis_coarse_->Comm().MyPID()==0)
  {
    NumMyElements = AllOutputEleIds.size();
  }
  else
  {
    NumMyElements= 0;
  }

  IO::cout << "Proc "<< actdis_coarse_->Comm().MyPID() << " NumOutputele " << my_output_elements_.size() << IO::endl;

  OutputMap_ = Teuchos::rcp(new Epetra_Map (NumGlobalMapElements,NumMyElements,&(AllOutputEleIds[0]),0,actdis_coarse_->Comm()));
}

void STR::UQ::MLMC::SetupEvalPeakWallStress()
{
  int  num_my_wall_elements =0;
 // Now loop over all elements in EleRowMap

  //TODO define better criterion for peak stress evaluation than hard coded material law
  for(int i=0;i<actdis_coarse_->ElementRowMap()->NumMyElements();i++)
  {
    // check wether element has AAANeohookeStopro Material as indicator for AAA wall
    if(actdis_coarse_->lRowElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      num_my_wall_elements++;
      my_wall_elements_.push_back(actdis_coarse_->lRowElement(i)->Id());
    }
    // check wether element has AAANeohookeStopro Material as indicator for AAA wall
    else if(actdis_coarse_->lRowElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke)
    {
      num_my_wall_elements++;
      my_wall_elements_.push_back(actdis_coarse_->lRowElement(i)->Id());
    }
  }//EOF Loop over Elements
  // Broadcast the information of how many wall elements we have to all other procs
  actdis_coarse_->Comm().SumAll(&num_my_wall_elements,&tot_num_wall_elements_,1);
}

void STR::UQ::MLMC::EvalPeakWallStress(Teuchos::RCP<const Epetra_Vector> disp, INPAR::STR::StressType iostress,INPAR::STR::StrainType iostrain )
{
  // Eval all wall elements on proc
  Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_disp =Teuchos::rcp(new std::vector <std::vector<double > >);
  Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_stresses =Teuchos::rcp(new std::vector <std::vector<double > >);
  Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_strains =Teuchos::rcp(new std::vector <std::vector<double > >);
  Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_mat_params =Teuchos::rcp(new std::vector <std::vector<double > >);
  EvalDisAtEleCenters(disp,
      iostress,
      iostrain,
      &my_wall_elements_,
      my_output_elements_c_disp,
      my_output_elements_c_stresses,
      my_output_elements_c_strains,
      my_output_elements_mat_params);

  // compute von Mises stress
  Teuchos::RCP< std::vector <double> > vonMisesStress = Teuchos::rcp(new std::vector<double > );
  Teuchos::RCP< std::vector <double> > vonMisesStrain = Teuchos::rcp(new std::vector<double > );
  Teuchos::RCP< std::vector <double> > DispMag = Teuchos::rcp(new std::vector<double > );
  Teuchos::RCP< std::vector <double> > Beta = Teuchos::rcp(new std::vector<double > );
  Teuchos::RCP< std::vector <double> > Youngs = Teuchos::rcp(new std::vector<double > );

  CalcVonMises(my_output_elements_c_stresses,vonMisesStress);
  CalcVonMises(my_output_elements_c_strains,vonMisesStrain);

  // compute disp magnitude
  for(unsigned int i=0;i<my_output_elements_c_disp->size();i++)
  {
    DispMag->push_back(sqrt( pow(my_output_elements_c_disp->at(i).at(0),2)+pow(my_output_elements_c_disp->at(i).at(1),2)+pow(my_output_elements_c_disp->at(i).at(2),2) ));
    // since my_output_elements_mat_params is the same size we use the same for loop here
    Beta->push_back(my_output_elements_mat_params->at(i).at(1));
    Youngs->push_back(my_output_elements_mat_params->at(i).at(0));
  }

  //calculate maximum and 99% quantile across procs
  // temp variables to store peak values in
  Teuchos::RCP<std::pair <int , double > > temp_peak = Teuchos::rcp(new std::pair <int , double  >);
  Teuchos::RCP<std::pair <int , double > > temp_quantile99 = Teuchos::rcp(new std::pair <int , double > );
  // list to store the peak values in
  Teuchos::RCP<std::list<std::pair <int , double > > > peak_vals = Teuchos::rcp(new std::list<std::pair <int , double  > >);
  Teuchos::RCP<std::list<std::pair <int , double > > > quantile99_vals = Teuchos::rcp(new std::list<std::pair <int , double  > >);


  ComputePeakAndQuantile(vonMisesStress,&my_wall_elements_,temp_peak,temp_quantile99);
  peak_vals->push_back(*temp_peak);
  quantile99_vals->push_back(*temp_quantile99);

  ComputePeakAndQuantile(vonMisesStrain,&my_wall_elements_,temp_peak,temp_quantile99);
  peak_vals->push_back(*temp_peak);
  quantile99_vals->push_back(*temp_quantile99);

  ComputePeakAndQuantile(DispMag,&my_wall_elements_,temp_peak,temp_quantile99);
  peak_vals->push_back(*temp_peak);
  quantile99_vals->push_back(*temp_quantile99);

  ComputePeakAndQuantile(Youngs,&my_wall_elements_,temp_peak,temp_quantile99);
   peak_vals->push_back(*temp_peak);
   quantile99_vals->push_back(*temp_quantile99);

  ComputePeakAndQuantile(Beta,&my_wall_elements_,temp_peak,temp_quantile99);
  peak_vals->push_back(*temp_peak);
  quantile99_vals->push_back(*temp_quantile99);

  ExportPeakStressDataAndWriteToFile(peak_vals,quantile99_vals);

}

void STR::UQ::MLMC::CalcVonMises(Teuchos::RCP<std::vector <std::vector<double > > > input_components , Teuchos::RCP< std::vector <double> > output_vM )
{
  for(unsigned int i=0;i<input_components->size();i++)
  {
    output_vM->push_back(sqrt( pow(input_components->at(i).at(0),2)+pow(input_components->at(i).at(1),2)+pow(input_components->at(i).at(2),2)
                                -input_components->at(i).at(0)*input_components->at(i).at(1)
                                -input_components->at(i).at(1)*input_components->at(i).at(2)
                                -input_components->at(i).at(2)*input_components->at(i).at(0)
                                +3*(
                                    pow(input_components->at(i).at(3),2)
                                   +pow(input_components->at(i).at(4),2)
                                  +pow(input_components->at(i).at(5),2))));
  }

}

void STR::UQ::MLMC::ComputePeakAndQuantile(Teuchos::RCP<std::vector <double > >  values , std::vector <int>* ele_ids, Teuchos::RCP<std::pair<int, double> > peak,Teuchos::RCP<std::pair<int, double> > quantile99)
{
  // Build a list< pair< int eleGID, double some quantity e.g. vM Stress > >
  // and sort it
  std::list< std::pair< int, double > > my_quantity;


  //loop over all element ids
  std::pair< int, double> thispair;
  for(unsigned int i=0;i<ele_ids->size();i++)
  {
    thispair.first  = ele_ids->at(i);
    thispair.second = values->at(i);
    my_quantity.push_back( thispair);
  }

  // sort the the listin ascending order by the estimated distance
  // this is the STL sorting, which is pretty fast
  my_quantity.sort(MyComparePairs);

  //now we need to communicate top ten percent
  int num_top_ten = (int) (tot_num_wall_elements_*0.1+0.5);
  thispair.first  = -1;
  thispair.second = -10000000;
  //resize to top ten fill up the remainder if needed with small number and negative ele id
  my_quantity.resize(num_top_ten,thispair);

  // gather does not like lists, thus
  std::vector<int> my_quantity_id_vec(0);
  std::vector<double> my_quantity_value_vec(0);

  // put data in separate vectors because Gather doesnt like vectors of pairs
  for (std::list< std::pair< int, double > >::iterator it = my_quantity.begin(); it != my_quantity.end(); it++)
  {
    my_quantity_id_vec.push_back(it->first);
    my_quantity_value_vec.push_back(it->second);
  }

  // vector to gather all the data in
  std::vector<double>  my_quantity_values_gathered_vec(0);
  std::vector<int>  my_quantity_id_gathered_vec(0);

  // gather informatio across all procs
  const int my_numprocs = actdis_coarse_->Comm().NumProc();

  // information how many processors participate in total
  std::vector<int> allproc(actdis_coarse_->Comm().NumProc());
  for (int i=0;i<actdis_coarse_->Comm().NumProc();++i) allproc[i] = i;

  LINALG::Gather<int>(my_quantity_id_vec,my_quantity_id_gathered_vec,my_numprocs,&allproc[0], actdis_coarse_->Comm());
  LINALG::Gather<double>(my_quantity_value_vec,my_quantity_values_gathered_vec,my_numprocs,&allproc[0], actdis_coarse_->Comm());

  //put data back into list for sorting
  std::list< std::pair< int, double > > my_quantity_gathered;
  for (unsigned int i=0; i<my_quantity_id_gathered_vec.size();i++)
  {
    std::pair<int,double> mytemp;
    mytemp.first=my_quantity_id_gathered_vec[i];
    mytemp.second=my_quantity_values_gathered_vec[i];
    my_quantity_gathered.push_back(mytemp);
  }
  // sort the list
  my_quantity_gathered.sort(MyComparePairs);

  peak->first=my_quantity_gathered.begin()->first;
  peak->second=my_quantity_gathered.begin()->second;

   std::list< std::pair< int ,double> >::iterator it = my_quantity_gathered.begin();
   std::advance( it,(int)(tot_num_wall_elements_*0.01));
   quantile99->first=it->first;
   quantile99->second=it->second;
}


void STR::UQ::MLMC::EvalDisAtEleCenters(Teuchos::RCP<const Epetra_Vector> disp,
      INPAR::STR::StressType iostress,
      INPAR::STR::StrainType iostrain,
      std::vector <int> * output_elements,
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_disp,
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_stresses,
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_c_strains,
      Teuchos::RCP<std::vector <std::vector<double > > > my_output_elements_mat_params)
{

  // we need the displacements in vector with colmap layout
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("displacement",disp);
  // get the displacements in col map layout
  Teuchos::RCP<const Epetra_Vector> disp_colmap = actdis_coarse_->GetState("displacement");
  for(unsigned int i = 0; i<output_elements->size(); i++)
  {
    std::vector<double>  my_c_disp (3, 0.0);
    int myNumNodes= actdis_coarse_->gElement(output_elements->at(i))->NumNode();
    const int* myNodeIds = actdis_coarse_->gElement(output_elements->at(i))->NodeIds();
    for (int k=0; k< myNumNodes ; k++)
    {
      const DRT::Node* node = actdis_coarse_->gNode(myNodeIds[k]);
      std::vector <int> myDofsPerNode = actdis_coarse_->Dof(node);
      for (unsigned int l=0; l<myDofsPerNode.size(); l++)
      {
        my_c_disp[l]+=1./myNumNodes*(*disp_colmap)[disp_colmap->Map().LID(myDofsPerNode[l])];
      }
    }
    my_output_elements_c_disp->push_back(my_c_disp);

    std::vector<double> mat_params (2, 0.0);
    // get the mat parameters
    if(actdis_coarse_->gElement(output_elements->at(i))->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->gElement(output_elements->at(i))->Material().get());

      mat_params[0]=(aaa_stopro->Youngs());
      mat_params[1]=(aaa_stopro->Beta());
    }
    if(actdis_coarse_->gElement(output_elements->at(i))->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke)
    {
      // Get LID
      int mylid=actdis_coarse_->gElement(output_elements->at(i))->LID();
      mat_params[0]=actdis_coarse_->gElement(output_elements->at(i))->Material()->Parameter()->GetParameter(0,mylid);
      mat_params[1]=actdis_coarse_->gElement(output_elements->at(i))->Material()->Parameter()->GetParameter(2,mylid);
    }
    // else ,no need initialized
    my_output_elements_mat_params->push_back(mat_params);
  }
  // Now we need to get ele stresses and strains
  Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > plstrain = Teuchos::rcp(new std::vector<char>());

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  p.set("action","calc_struct_stress");
  p.set("stress", stress);
  p.set("plstrain",plstrain);
  p.set("strain", strain);
  p.set<int>("iostress", iostress);
  p.set<int>("iostrain", iostrain);

  Teuchos::RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  Teuchos::RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  // disp is passed to the function no need for reading in results anymore
  actdis_coarse_->SetState("displacement",disp);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alright lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");
  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstrainmap", gpstrainmap);
  //actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Evaluate2(p,output_elements,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();
  // also get element stresses
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", gpstressmap);
  Teuchos::RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
  p.set("poststress",elestress);
  Evaluate2(p,output_elements,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  //actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }
  // and strains
  p.set("gpstressmap", gpstrainmap);
  Teuchos::RCP<Epetra_MultiVector> elestrains = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
  p.set("poststress",elestrains);
  Evaluate2(p,output_elements,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  //actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  // now we have stresses and strains in elestress and elestrain in  Element rowmap layout
  // hence loop over the elements
  for(unsigned int i = 0; i<output_elements->size(); i++)
  {
    std::vector<double>  element_c_stresses;
    std::vector<double>  element_c_strains;
    for(int k = 0;k<6 ; k++)
    {
      element_c_stresses.push_back((*elestress)[k][(elestress->Map().LID(output_elements->at(i)))]);
      element_c_strains.push_back((*elestrains)[k][(elestrains->Map().LID(output_elements->at(i)))]);
    }
    my_output_elements_c_stresses->push_back(element_c_stresses);
    my_output_elements_c_strains->push_back(element_c_strains);
  }
}

void STR::UQ::MLMC::ExportPeakStressDataAndWriteToFile(
      Teuchos::RCP<std::list <std::pair <int , double > > >  peak_values,
      Teuchos::RCP<std::list <std::pair <int , double > > >  quantile99_values)
{
  if (actdis_coarse_->Comm().MyPID()==0)
  {
    // assemble name for outputfile
    std::stringstream outputfile2;
    outputfile2 << filename_ << "_statistics_output_" << start_run_ << "_max_values" << ".txt";
    std::string name = outputfile2.str();;
    // file to write output
    std::ofstream File;
    if (numb_run_ == 0 || numb_run_ == start_run_)
    {
      File.open(name.c_str(),std::ios::out);
      if (File.is_open())
      {
        File << "runid "<< "vMstress EID vMstrain EID Phi EID disp EID youngs EID beta EID vMstress99 EID vMstrain99 EID Phi99 EID disp99 EID youngs99 EID beta99 EID "
            << std::endl;
        File.close();
      }
      else
      {
        dserror("Unable to open statistics output file");
      }
     }

    // reopen in append mode
    File.open(name.c_str(),std::ios::app);
    File << numb_run_ ;
    for (std::list< std::pair< int, double > >::iterator it = peak_values->begin(); it != peak_values->end(); it++)
    {
      File << " " << it->second << " "<< it->first;
    }
    for (std::list< std::pair< int, double > >::iterator it = quantile99_values->begin(); it != quantile99_values->end(); it++)
    {
      File << " " << it->second << " "<< it->first;
    }
    File << std::endl;
    File.close();
    }

  actdis_coarse_->Comm().Barrier();
}

void STR::UQ::MLMC::ExportEleDataAndWriteToFile(Teuchos::RCP<const Epetra_Map> OutputMap,
      INPAR::STR::StressType iostress,
      INPAR::STR::StrainType iostrain,
      std::vector <int> * output_elements,
      Teuchos::RCP<std::vector <std::vector<double > > >  my_output_elements_c_disp,
      Teuchos::RCP<std::vector <std::vector<double > > >  my_output_elements_c_stresses,
      Teuchos::RCP<std::vector <std::vector<double > > >  my_output_elements_c_strains,
      Teuchos::RCP<std::vector <std::vector<double > > >  my_output_elements_mat_params)
{
  // now we need to setup a map
  // set up map <int,DRT::CONTAINER>
  std::map<int,Teuchos::RCP <DRT::Container> >my_output_element_map;
  for(unsigned int i = 0; i<output_elements->size(); i++)
  {
    // put all the stuff into container
    Teuchos::RCP<DRT::Container> mycontainer = Teuchos::rcp(new DRT::Container);
    mycontainer->Add("stresses",(my_output_elements_c_stresses->at(i)));
    mycontainer->Add("strains",(my_output_elements_c_strains->at(i)));
    mycontainer->Add("disp",(my_output_elements_c_disp->at(i)));
    mycontainer->Add("mat_params",(my_output_elements_mat_params->at(i)));
    my_output_element_map.insert( std::pair <int,Teuchos::RCP< DRT::Container> >(output_elements->at(i),mycontainer) );
  }
  // build exporter
  DRT::Exporter myexporter(*(actdis_coarse_->ElementRowMap()),*OutputMap,actdis_coarse_->Comm());
  // this actually transfers everything to proc 0
  myexporter.Export(my_output_element_map);

  if (actdis_coarse_->Comm().MyPID()==0)
  {
    std::map<int,Teuchos::RCP <DRT::Container> >::iterator myit;
    for ( myit=my_output_element_map.begin() ; myit != my_output_element_map.end(); myit++ )
    {
      // get back all the data
     const std::vector<double> * stresses = myit->second()->Get< std::vector <double> >("stresses");
     const std::vector<double> * strains = myit->second()->Get< std::vector <double> >("strains");
     const std::vector<double> * mat_params = myit->second()->Get< std::vector <double> >("mat_params");
     const std::vector<double> * disp = myit->second()->Get< std::vector <double> >("disp");
     // IO::cout << "my_output_element_map first " << myit->first << *(myit->second)  << IO::endl;


    // assamble name for outputfile
    std::stringstream outputfile2;
    std::string stresstype;
    std::string straintype;
    if (iostress==INPAR::STR::stress_cauchy)
      stresstype="cauchy";
    else if (iostress==INPAR::STR::stress_2pk)
      stresstype="piola2";
    else
      dserror("unknown stresstype");
    if (iostrain==INPAR::STR::strain_ea)
      straintype="ea";
    else if (iostrain==INPAR::STR::strain_gl)
      straintype="gl";
    else
          dserror("unknown straintype");

    outputfile2 << filename_ << "_statistics_output_" << start_run_ << "_stress_" << stresstype << "_strain_" << straintype << "_EleId_"<< myit->first << ".txt";
    std::string name2 = outputfile2.str();;
    // file to write output
    std::ofstream File;
    if (numb_run_ == 0 || numb_run_ == start_run_)
    {
      File.open(name2.c_str(),std::ios::out);
      if (File.is_open())
      {
        File << "runid "<< "dispx dispy dispz stressxx stressyy stresszz stressxy stressyz stressxz strainxx  strainyy strainzz strainxy strainyz strainxz mat_param1 mat_param2"
            << std::endl;
        File.close();
      }
      else
      {
        dserror("Unable to open statistics output file");
      }
     }
    // reopen in append mode
    File.open(name2.c_str(),std::ios::app);
    File << numb_run_ ;
    for(unsigned int i=0;i<disp->size();i++)
    {
      File <<  " " << (*disp)[i];
    }
    for(unsigned int i=0;i<stresses->size();i++)
    {
      File <<  " " << (*stresses)[i];
    }
    for(unsigned int i=0;i<stresses->size();i++)
    {
      File <<  " " << (*strains)[i];
    }
    for(unsigned int i=0;i<mat_params->size();i++)
    {
       File << " " << (*mat_params)[i];
    }
    File << std::endl;
    File.close();
    }
  }
  actdis_coarse_->Comm().Barrier();
}
void STR::UQ::MLMC::WriteParamContInfoToFile(
      Teuchos::RCP<std::vector <int> > paramcont_info)
{
  if (actdis_coarse_->Comm().MyPID()==0)
  {
    // assemble name for outputfile
    std::stringstream outputfile2;
    outputfile2 << filename_ << "_parametercontinuation_output_" << start_run_ << ".txt";
    std::string name = outputfile2.str();;
    // file to write output
    std::ofstream File;
    if (numb_run_ == 0 || numb_run_ == start_run_)
    {
      File.open(name.c_str(),std::ios::out);
      if (File.is_open())
      {
        File << "runid "<< "num_cont_steps_planned num_cont_steps_needed  num_trials"
            << std::endl;
        File.close();
      }
      else
      {
        dserror("Unable to open statistics output file");
      }
     }

    // reopen in append mode
    File.open(name.c_str(),std::ios::app);
    File << numb_run_ ;
    for (std::vector <int>::iterator it = paramcont_info->begin(); it != paramcont_info->end(); it++)
    {
      File << " " << *it;
    }
    File << std::endl;
    File.close();
    }

  actdis_coarse_->Comm().Barrier();
}

/*----------------------------------------------------------------------*
 |  evaluate (public) / basically copy of the Evaluate function of
 |  the discretization                                                   |
 *----------------------------------------------------------------------*/
void STR::UQ::MLMC::Evaluate2(
                        Teuchos::ParameterList&              params,
                        std::vector <int>  * eval_elements,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                        Teuchos::RCP<Epetra_Vector>          systemvector1,
                        Teuchos::RCP<Epetra_Vector>          systemvector2,
                        Teuchos::RCP<Epetra_Vector>          systemvector3)
{
  DRT::AssembleStrategy strategy( 0, 0, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3 );

  TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate");

  if (!actdis_coarse_->Filled()) dserror("FillComplete() was not called");
  if (!actdis_coarse_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // call the element's register class preevaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate PreEvaluate");
    DRT::ParObjectFactory::Instance().PreEvaluate(*actdis_coarse_,params,
                                             strategy.Systemmatrix1(),
                                             strategy.Systemmatrix2(),
                                             strategy.Systemvector1(),
                                             strategy.Systemvector2(),
                                             strategy.Systemvector3());
  }
  // in the orginal function we have la(dofsets_.size()); here
  DRT::Element::LocationArray la(1);
  //my_output_elements_
  for (unsigned int i=0; i<eval_elements->size(); ++i)
  {
    // only evaluate the necessary elements for efficiency
    DRT::Element* actele = actdis_coarse_->gElement(eval_elements->at(i));

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate LocationVector");
    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(*actdis_coarse_,la,false);
    }

    {
    TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate Resize");

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    strategy.ClearElementStorage( la[row].Size(), la[col].Size() );
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate elements");
    // call the element evaluate method
    int err = actele->Evaluate(params,*actdis_coarse_,la,
                               strategy.Elematrix1(),
                               strategy.Elematrix2(),
                               strategy.Elevector1(),
                               strategy.Elevector2(),
                               strategy.Elevector3());
    if (err) dserror("Proc %d: Element %d returned err=%d",actdis_coarse_->Comm().MyPID(),actele->Id(),err);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR("DRT::Discretization::Evaluate assemble");
      int eid = actele->Id();
      strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
      strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
      strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
      strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
      strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );

    }
  } // for (int i=0; i<numcolele; ++i)

  return;
}



#endif // FFTW

