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
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "gen_randomfield.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_comm/comm_utils.H"
//for file output
#include <fstream>
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_element.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::MLMC::MLMC(Teuchos::RCP<DRT::Discretization> dis,
                Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    output_(output)
    //sti_(Teuchos::null)
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

  // use deterministic value yes/no
  use_det_value_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"USEDETVALUE");

   // get det value if needed
  det_value_ =mlmcp.get<double>("DETVALUE");

  // get starting random seed
  start_random_seed_ = mlmcp.get<int>("INITRANDOMSEED");

  // get name of lower level outputfiles
  filename_lower_level_ =  mlmcp.get<std::string>("OUTPUT_FILE_OF_LOWER_LEVEL");

  // get numerb of current level
  num_level_ = mlmcp.get<int>("LEVELNUMBER");

  //write statistics every write_stat_ steps
  write_stats_ = mlmcp.get<int>("WRITESTATS");

  // get OutputElements

  double word;
  std::istringstream bsdampingstream(Teuchos::getNumericStringParameter(mlmcp,"OUTPUT_ELEMENT_IDS"));
  while (bsdampingstream >> word)
	  AllMyOutputEleIds_.push_back(int (word));

  if(AllMyOutputEleIds_.front()== -1)
	  IO::cout << RED_LIGHT "No elements specified for output " END_COLOR << IO::endl;

  // In element critirion xsi_i < 1 + eps  eps = MLMCINELETOL
  InEleRange_ = 1.0 + 10e-3;
  //ReadInParameters();


  // controlling parameter
  start_run_ = mlmcp.get<int>("START_RUN");
  int numruns = mlmcp.get<int>("NUMRUNS");
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  int NNestedGroups = problem->GetNPGroup()->NumGroups();
  int i = problem->GetNPGroup()->GroupId();

  numruns_pergroup_= int(ceil(numruns/NNestedGroups));
  start_run_  += (i)*numruns_pergroup_;

  numb_run_ =  start_run_;//+numruns_pergroup_;     // counter of how many runs were made monte carlo

  // init number of wall elements
  tot_num_wall_elements_ = 0;
  SetupEvalDisAtEleCenters(AllMyOutputEleIds_);
  SetupEvalPeakWallStress();


  reduced_output_ = DRT::INPUT::IntegralValue<int>(mlmcp ,"REDUCED_OUTPUT");


  // store rf at all element centers
  rf_values_ = Teuchos::rcp(new std::vector<double>(discret_->NumMyColElements(),0.0));


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
    output_fine_ = Teuchos::rcp(new IO::DiscretizationWriter(actdis_fine_,output_control_fine_));

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
void STR::MLMC::Integrate()
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
    SetupStochMat((random_seed+(unsigned int)numb_run_));
    const double t2 = Teuchos::Time::wallTime();
    discret_->Comm().Barrier();


    if(!reduced_output_)
    {
      output_->NewResultFile(filename_,(numb_run_));
      output_->WriteMesh(1, 0.01);
    }

    else
    {
      IO::cout << "ATTENTION NO NEW RESULTFILE CREATED" << IO::endl;
    }
    // get input lists
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    //store time
    double t3 = 0;
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
    EvalDisAtEleCenters(dis_coarse,INPAR::STR::stress_cauchy  ,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);
    ExportEleDataAndWriteToFile(OutputMap_,INPAR::STR::stress_cauchy,INPAR::STR::strain_ea,&my_output_elements_,my_output_elements_c_disp,my_output_elements_c_stresses,my_output_elements_c_strains,my_output_elements_mat_params);

    // Evaluate Peak stress and strain vallues
    EvalPeakWallStress(dis_coarse, INPAR::STR::stress_cauchy,INPAR::STR::strain_ea );


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

  //ReadResultsFromLowerLevel(12);
  return;
}
/*----------------------------------------------------------------------*/
/* Nice and clean version with for not resetting prestress*/
void STR::MLMC::IntegrateNoReset()
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

    if(!reduced_output_)
    {
        output_->NewResultFile(filename_,(numb_run_));
        output_->WriteMesh(1, 0.01);
    }

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
          SetupStochMatDet(3.8);
          structadaptor.Integrate();
          // get all information
          structadaptor.GetRestartData(cont_step_,cont_time_,cont_disn_init_,cont_veln_,cont_accn_,cont_elementdata_init_);
          // Get timestep size
          //time_step_size=structadaptor.GetTimeStepSize();
        }
        else
        {
          // compute solution with parameter continuation
          int error = ParameterContinuation( num_cont_steps, random_seed,false, structadaptor);
          // if not converged
          int num_trials = 1;
          int max_trials = 4;
          int num_cont_steps_new = num_cont_steps;
          while( error && num_trials < max_trials)
          {
            // double number of cont steps
            num_cont_steps_new = num_cont_steps_new*2;
            error = ParameterContinuation( num_cont_steps_new, random_seed, true, structadaptor);
            num_trials++;
            // check for error
            if(num_trials==max_trials && error)
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
int STR::MLMC::ParameterContinuation(unsigned int num_cont_steps, unsigned int random_seed, bool re_use_rf, ADAPTER::Structure&  structadaptor)
{
  int error = 0;
  t1_ = Teuchos::Time::wallTime();
  // init some variables
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  *cont_elementdata_=*cont_elementdata_init_;
  *cont_disn_=*cont_disn_init_;

  // Get timestep size
  double time_step_size=structadaptor.GetTimeStepSize();

  // for loop to encapsulate parameter continuation scheme
  for(unsigned int k =0 ; k<num_cont_steps ; k++)
  {
     double gamma=1.0/num_cont_steps+(k*(1.0/num_cont_steps));
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
       structadaptor.SetRestart(num_prestress_steps-1,(num_prestress_steps-1)*time_step_size,cont_accn_,cont_veln_,cont_accn_,cont_elementdata_);
       // set new material parameters
       // Setup Material Parameters in each element based on deterministic value
       //measure time only if we really compute the random field (k=0)
       if(!re_use_rf)
       {
         t2_ = Teuchos::Time::wallTime();
         BlendStochMat((random_seed+(unsigned int)numb_run_), false , 1.5 ,gamma);
         t3_ = Teuchos::Time::wallTime();
         // we only want to setupt once hence
         re_use_rf=true;
       }
       else if(re_use_rf)
       {
         BlendStochMat((random_seed+(unsigned int)numb_run_), true , 1.5 ,gamma);
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
       structadaptor.GetRestartData(garbage1,garbage2,cont_accn_,cont_veln_,cont_accn_,cont_elementdata_);
       // write new prestress data together with old displament data to discretization
       structadaptor.SetRestart((*(cont_step_)-1),(*cont_time_)-time_step_size,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_);
       // reset old maxtime
       structadaptor.SetTimeEnd(endtime);
       //BlendStochMat((random_seed+(unsigned int)numb_run_),(bool)(k+1),4.61,gamma);
       BlendStochMat((random_seed+(unsigned int)numb_run_),true,1.5,gamma);
       error = structadaptor.Integrate();
       if(error)
         return error;
       structadaptor.GetRestartData(garbage1,garbage2,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_);
     }
     else // no prestress
     {
       // set prestress state and time, disp acc vel are al set to zero
       structadaptor.SetRestart((*(cont_step_)-1),(*cont_time_)-time_step_size,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_);
       //measure time only if we really compute the random field (k=0)
       if(!k)
         t2_ = Teuchos::Time::wallTime();
       BlendStochMat((random_seed+(unsigned int)numb_run_), (bool)k , 4.61,gamma);
       if(!k)
         t3_ = Teuchos::Time::wallTime();
       error = structadaptor.Integrate();
       if(error)
         return error;
       structadaptor.GetRestartData(cont_step_,cont_time_,cont_disn_,cont_veln_,cont_accn_,cont_elementdata_);
     }

  }// eof parameter continuation scheme
    t4_ = Teuchos::Time::wallTime();
    return 0;
}
//---------------------------------------------------------------------------------------------
void STR::MLMC::SetupProlongatorParallel()
{
  dserror("no longer maintianed");
/*  // number of outliers
  int num_outliers = 0;
  const Epetra_Map* rmap_disp = NULL;
  const Epetra_Map* dmap_disp = NULL;
  const Epetra_Map* rmap_stress = NULL;
  const Epetra_Map* dmap_stress = NULL;
  rmap_disp= (actdis_fine_->DofRowMap());
  dmap_disp=(actdis_coarse_->DofRowMap());
  //maps for stress prolongator
  rmap_stress = (actdis_fine_->NodeRowMap());
  dmap_stress = (actdis_coarse_->NodeRowMap());
  int bg_ele_id;
  // store location of node
  double xsi[3] = {0.0, 0.0,0.0};
  prolongator_disp_crs_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap_disp,*dmap_disp,8,false));
  prolongator_stress_crs_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap_stress,*dmap_stress,8,false));

  //loop over nodes of fine discretization on this proc
  Teuchos::RCP<Epetra_Vector> node_vector = Teuchos::rcp(new Epetra_Vector(*(actdis_fine_->NodeRowMap()),true));
  // loop over dofs
  for (int i=0; i< actdis_fine_->NumMyRowNodes(); i++  )
  //  for (int i=0; i<2; i++  )
  {
    DRT::Node* node = actdis_fine_->lRowNode(i);

    // Get background element and local coordinates
    num_outliers += FindBackgroundElement(*node, actdis_coarse_, &bg_ele_id, xsi);

    // Get element
    DRT::Element* bg_ele = actdis_coarse_->gElement(bg_ele_id);

    Epetra_SerialDenseVector shape_fcts(bg_ele->NumNode());
    DRT::UTILS::shape_function_3D(shape_fcts,xsi[0],xsi[1],xsi[2],bg_ele->Shape());


    int* rows = NULL;
    int* cols = NULL;
    double* values = NULL;
    int numColumns = bg_ele->NumNode();
    // insert rows dof wise
    int numRows   = 1;


    cols = new int[numColumns];
    rows = new int[numRows];
    values = new double[numColumns];

    // fill prolongators
    // DIM = 3
    for (int j = 0; j<3 ; j++)
    {
      int index = 0;
      for (int k=0; k< bg_ele->NumNode(); k ++)
      {
        // store global indices in cols
       // (*rcp_Array)[index]=dmap_disp->GID((bg_ele->Nodes()[k]->Id()*3)+j);
        cols[index]=dmap_disp->GID((bg_ele->Nodes()[k]->Id()*3)+j);
        rows[0]=i*3+j;
        values[index]=shape_fcts[k];
        index++;
      }
      int err=  prolongator_disp_crs_->InsertGlobalValues(1,rows,8,cols,values,Epetra_FECrsMatrix::COLUMN_MAJOR);
      if (err != 0)
      {
        dserror("Could not insert global values");
      }

    } // loop j
    // stress prolongator
    for (int k=0; k< bg_ele->NumNode(); k ++)
    {
      // store global indices in cols
      cols[k]=dmap_stress->GID(bg_ele->Nodes()[k]->Id());
      values[k]=shape_fcts[k];
    }
    rows[0]=i;
    int err=  prolongator_stress_crs_->InsertGlobalValues(1,rows,8,cols,values,Epetra_FECrsMatrix::COLUMN_MAJOR);
    if (err != 0)
    {
      dserror("Could not insert global values");
    }

    delete [] cols;
    delete [] rows;
    delete [] values;

    rows = NULL;
    cols = NULL;
    values = NULL;

  } // End of loop over nodes of fine discretzation

  // Assembly
  prolongator_disp_crs_->GlobalAssemble(*dmap_disp,*rmap_disp,true);
  prolongator_stress_crs_->GlobalAssemble(*dmap_stress,*rmap_stress,true);
  IO::cout << "################################################### " << IO::endl;
  IO::cout << "   SUCCESSFULLY INITIALIZED  PROLONGATOR" << IO::endl;
  IO::cout <<  num_outliers << " Nodes do not lie within a background element " << IO::endl;
  IO::cout << "################################################### " << IO::endl; */
}
//---------------------------------------------------------------------------------------------
void STR::MLMC::SetupProlongator()
{
  // This functions calculates the prolongtators for the displacement and the nodal stresses

  // 3D Problem
  int num_columns_prolongator_disp = actdis_fine_->NumGlobalNodes()*3;
  int num_columns_prolongator_stress = actdis_fine_->NumGlobalNodes();

  double xsi[3] = {0.0, 0.0,0.0};

  // loop over nodes of fine dis
  int num_nodes;
  int bg_ele_id;
  num_nodes = actdis_fine_->NumGlobalNodes();

  // init prolongators
  IO::cout << "num_columns proongator "  << num_columns_prolongator_disp << IO::endl;
  IO::cout << "num_columns proongator stress  "  << num_columns_prolongator_stress << IO::endl;
  prolongator_disp_ = Teuchos::rcp(new Epetra_MultiVector(*actdis_coarse_->DofRowMap(),num_columns_prolongator_disp,true));
  prolongator_stress_ = Teuchos::rcp(new Epetra_MultiVector(*actdis_coarse_->NodeRowMap(),num_columns_prolongator_stress,true));

  for (int i = 0; i < num_nodes ; i++)
  {

    // Get node
    DRT::Node* node = actdis_fine_->gNode(i);

    // Get background element and local coordinates
    FindBackgroundElement(*node, actdis_coarse_, &bg_ele_id, xsi);
    // Get element
    DRT::Element* bg_ele = actdis_coarse_->gElement(bg_ele_id);

    Epetra_SerialDenseVector shape_fcts(bg_ele->NumNode());
    DRT::UTILS::shape_function_3D(shape_fcts,xsi[0],xsi[1],xsi[2],bg_ele->Shape());

    // fill prolongators add values to prolongator
    for (int j = 0; j<3 ; j++)
    {

      for (int k=0; k< bg_ele->NumNode(); k ++)
      {
        //prolongator_disp_ is dof based
        (*prolongator_disp_)[(i*3)+j][(bg_ele->Nodes()[k]->Id()*3)+j]= shape_fcts[k];
        // prolongator stress_ is node based
        (*prolongator_stress_)[i][(bg_ele->Nodes()[k]->Id())]= shape_fcts[k];
      }
    }
  } // End of loop over nodes of fine discretzation

}
//---------------------------------------------------------------------------------------------
void STR::MLMC::ProlongateResults()
{
  IO::cout << "Prolongating Resuts " << IO::endl;
  // To avoid messing with the timeintegration we read in the results of the coarse discretization here
  // Get coarse Grid problem instance
  std::stringstream name;
  std::string filename_helper;
  name << filename_ << "_run_"<< numb_run_;
  filename_helper = name.str();

  RCP<IO::InputControl> inputcontrol = Teuchos::rcp(new IO::InputControl(filename_helper, false));

  IO::DiscretizationReader input_coarse(actdis_coarse_, inputcontrol,tsteps_);
  // Vector for displacements
  Teuchos::RCP<Epetra_Vector> dis_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));

  // read in displacements
  input_coarse.ReadVector(dis_coarse, "displacement");

  Teuchos::RCP<Epetra_MultiVector> dis_fine = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));

  // Try new and shiny prolongator based on crs Matrix
  int error = prolongator_disp_crs_->Multiply(false,*dis_coarse,*dis_fine);
  if(error!=0)
  {
    dserror("stuff went wrong");
  }


  // create new resultfile for prolongated results
  std::stringstream name_prolong;
  std::string filename_helper_prolong;
  name_prolong << filename_ << "_prolongated";
  filename_helper_prolong = name_prolong.str();

  output_fine_->NewResultFile(filename_helper_prolong ,(numb_run_));


  output_fine_->WriteMesh(1, 0.01, meshfilename_);
  // exception for first run
  if(numb_run_==start_run_)
    output_fine_->WriteMesh(1, 0.01);

  output_fine_->NewStep( 1, 0.01);


  // Write interpolated displacement to file
  output_fine_->WriteVector("displacement", dis_fine, output_fine_->dofvector);
  Teuchos::RCP<Epetra_Vector> dis_fine_single = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));

  // transfer to Epetra_Vector
  for( int i = 0;i< dis_fine->MyLength(); i++)
  {
   (*dis_fine_single)[i]= (*dis_fine)[0][i];
  }

  //#####################################################################################
  //
  //                  prolongate stresses based on interpolated displacement field
  //
  //#####################################################################################

  // set up parameters for Evaluation
  double timen         = 0.9;  // params_.get<double>("total time"             ,0.0);
  double dt            = 0.1; //params_.get<double>("delta time"             ,0.01);
  double alphaf        = 0.459; // params_.get<double>("alpha f"                ,0.459);
  INPAR::STR::StressType iostress =INPAR::STR::stress_2pk; //stress_none;
  INPAR::STR::StrainType iostrain= INPAR::STR::strain_gl; // strain_none;
  RCP<Epetra_Vector>    zeros_ = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));
  RCP<Epetra_Vector>    dis_ = dis_fine_single;
  RCP<Epetra_Vector>    vel_ = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements

  p.set("action","calc_struct_stress");
  // other parameters that might be needed by the elements
  p.set("total time",timen);
  p.set("delta time",dt);
  p.set("alpha f",alphaf);

  Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
  // plastic strains need to be init as well
  Teuchos::RCP<std::vector<char> > plstrain = Teuchos::rcp(new std::vector<char>());
  p.set("stress", stress);
  p.set("plstrain",plstrain);
  //

  p.set<int>("iostress", iostress);
  p.set("strain", strain);
  p.set<int>("iostrain", iostrain);
  // set vector values needed by elements
  p.set<double>("random test",5.0);
  actdis_fine_->ClearState();
  actdis_fine_->SetState("residual displacement",zeros_);
  actdis_fine_->SetState("displacement",dis_);
  actdis_fine_->SetState("velocity",vel_);

  // Evaluate Stresses based on interpolated displacements
  //actdis_fine_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  //actdis_fine_->ClearState();
  // Write to file
  // interpolated stresses from disp field
  //output_fine_->WriteVector("gauss_2PK_stresses_xyz",*stress,*(actdis_fine_->ElementRowMap()));




  //#####################################################################################
  //
  //                  prolongate stresses based on interpolated nodal stress field
  //
  //#####################################################################################

  // use same parameter list as above
  RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  //RCP<Epetra_Vector>    dis_ = dis_fine_single;
  RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  actdis_coarse_->SetState("displacement",dis_coarse);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alrigth lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");

  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstrainmap", gpstrainmap);

  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // st action to calc poststresse
  p.set("action","postprocess_stress");
  // Multivector to store poststresses
  RCP<Epetra_MultiVector> poststress =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeRowMap()),6,true));
  // for fine diskretization as well
  RCP<Epetra_MultiVector> poststress_fine =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));

  p.set("poststress", poststress);
  p.set("stresstype","ndxyz");

  actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();

  // again for strains
  p.set("action","postprocess_stress");
  // Multivector to store poststrains
  RCP<Epetra_MultiVector> poststrain =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeRowMap()),6,true));
  // for fine diskretization as well
  RCP<Epetra_MultiVector> poststrain_fine =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
  p.set("poststress", poststrain);
  p.set("gpstressmap", gpstrainmap);
  p.set("stresstype","ndxyz");

  actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();


  // try new and shiny crs prolongator
  int error2 = prolongator_stress_crs_->Multiply(false,*poststress,*poststress_fine);
  if(error2!=0)
  {
    dserror("stuff went wrong");
  }
  // strains as well
  int error3 = prolongator_stress_crs_->Multiply(false,*poststrain,*poststrain_fine);
  if(error3!=0)
  {
    dserror("stuff went wrong");
  }

  IO::cout << " before prolongated gauss 2pk stresses " << IO::endl;
  output_fine_->WriteVector("prolongated_gauss_2PK_stresses_xyz", poststress_fine, output_fine_->nodevector);
  output_fine_->WriteVector("prolongated_gauss_GL_strains_xyz", poststrain_fine, output_fine_->nodevector);
  // do some statistics
  if(calc_diff_)
  {
    CalcDifferenceToLowerLevel(poststress_fine,poststrain_fine, dis_fine);
    // Write Difference between two Discretizations to File
    output_fine_->WriteVector("diff_to_ll_displacement", disp_lower_level_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_to_ll_prolongated_gauss_2PK_stresses_xyz", stress_lower_level_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_to_ll_prolongated_gauss_GL_strains_xyz", strain_lower_level_, output_fine_->nodevector);
  }
  CalcStatStressDisp(poststress_fine,poststrain_fine, dis_fine);
  // write some output to another file
  HelperFunctionOutputTube(poststress_fine,poststrain_fine,dis_fine);
}
//-----------------------------------------------------------------------------------
int STR::MLMC::FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi)
{
  bool outlier = false;
  Teuchos::RCP<Epetra_Vector> element_vector = Teuchos::rcp(new Epetra_Vector(*(background_dis->ElementRowMap()),true));
  // loop over discretization

  double pseudo_distance = 0.0;
  double min_pseudo_distance = 2.0;
  int back_ground_ele_id = 0;
  double background_xsi[3] = {0,0,0};
  //for (int i=0; i< 2 && min_pseudo_distance > InEleRange_ ; i++  )
  for (int i=0; i< element_vector->MyLength() && min_pseudo_distance > InEleRange_ ; i++  )
  {
    int globel_ele_id= background_dis->ElementRowMap()->GID(i);
    DRT::Element* ele = background_dis->gElement(globel_ele_id);
    //inEle = CheckIfNodeInElement(node, *ele);
    pseudo_distance = CheckIfNodeInElement(node, *ele, xsi);
    //IO::cout << " pseudo_distance  " << pseudo_distance  << IO::endl;
    // check if node is in Element
    if (pseudo_distance < min_pseudo_distance)
    {
      min_pseudo_distance = pseudo_distance;
      back_ground_ele_id = globel_ele_id;
      background_xsi[0]=xsi[0];
      background_xsi[1]=xsi[1];
      background_xsi[2]=xsi[2];
    }
  } // end of loop over elements
  // Debug
  if(min_pseudo_distance < InEleRange_)
  {
    //IO::cout << "found background element Ele ID is " <<  back_ground_ele_id << IO::endl;
    //IO::cout << "Local Coordinates are "<< "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << IO::endl;

  }
  else
  {
   //IO::cout << "did not find background element, closest element is: Ele ID: " << back_ground_ele_id << IO::endl;
   //IO::cout << "Local Coordinates are "<< "xsi_0 " << background_xsi[0] << " xsi_1 " << background_xsi[1] << " xsi_2 " << background_xsi[2] << IO::endl;
   // dserror("stop right here");
   // write closest element xsi* into  pointer
    xsi[0]=background_xsi[0];
    xsi[1]=background_xsi[1];
    xsi[2]=background_xsi[2];
    outlier =true;
  }
  *bg_ele_id = back_ground_ele_id;
  return outlier;


}

//-----------------------------------------------------------------------------------
double STR::MLMC::CheckIfNodeInElement(DRT::Node& node, DRT::Element& ele, double* xsi)
{
  // init speudo distance, which is essentially largest values of xsi[i]
  double pseudo_distance = 0.0;
  //local/element coordinates
  xsi[0] = xsi[1] = xsi[2] = 0.0 ;
  // Replace later with problem dimension or element type
  //if (Dim()==3)

  if (ele.Shape()==DRT::Element::hex8)
    {
      // function f (vector-valued)
      double f[3] = {0.0, 0.0, 0.0};
      LINALG::Matrix<3,1> b;
      // gradient of f (df/dxsi[0], df/dxsi[1], df/dxsi[2])
      LINALG::Matrix<3,3> df;
      //convergeence check
      double conv = 0.0;

      for (int k=0;k<num_newton_it_;++k)
        {
          EvaluateF(f,node,ele,xsi);
          conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
          //IO::cout << "Iteration " << k << ": -> |f|=" << conv << IO::endl;
          if (conv <= convtol_) break;

          EvaluateGradF(df,node,ele,xsi);

          // solve dxsi = - inv(df) * f
          df.Invert();


          // update xsi
          xsi[0] += -df(0,0)*f[0] - df(1,0)*f[1] - df(2,0)*f[2];
          xsi[1] += -df(0,1)*f[0] - df(1,1)*f[1] - df(2,1)*f[2];
          xsi[2] += -df(0,2)*f[0] - df(1,2)*f[1] - df(2,2)*f[2];
        }

      // Newton iteration unconverged
      if (conv > convtol_)
      {
        dserror("ERROR: CheckIfNodeInElement: Newton unconverged for NodeID %i "
                   "and ElementID %i", node.Id(), ele.Id());
      }
      // Newton iteration converged
      // find largest value of xsi[i] and return as pseudo_distance
      for (int i = 0; i < 3; i++)
      {
        if (fabs(xsi[i]) > pseudo_distance)
        {
          pseudo_distance = fabs(xsi[i]);
        }
      }
      return pseudo_distance;

    }
    else
    {
      dserror("CheckIfNodeInElement only implememted for hex8 Elements");
      return pseudo_distance;
    }
}
//-----------------------------------------------------------------------------------
bool STR::MLMC::EvaluateF(double* f,DRT::Node& node, DRT::Element& ele,const double* xsi)
{

  LINALG::Matrix<8,1> funct ;
  //DRT::Element::DiscretizationType distype = ele.DiscretizationType;
  const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
  DRT::UTILS::shape_function_3D(funct, xsi[0], xsi[1], xsi[2],distype);

  //LINALG::Matrix<ele.NumNode(),ele->numdim> xrefe;  // material coord. of element
  LINALG::Matrix<8,3> xrefe;  // material coord. of element
    for (int i=0; i<ele.NumNode(); ++i){
      const double* x =ele.Nodes()[i]->X();
      xrefe(i,0) = x[0];
      xrefe(i,1) = x[1];
      xrefe(i,2) = x[2];
    }
    // Calc Difference in Location
    LINALG::Matrix<1,3> point;
    point.MultiplyTN(funct, xrefe);
    f[0]=point(0,0)-node.X()[0];
    f[1]=point(0,1)-node.X()[1];
    f[2]=point(0,2)-node.X()[2];
    return true;
}


//-----------------------------------------------------------------------------------
bool STR::MLMC::EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi)
{
  //static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  LINALG::Matrix<3,8> deriv1;
  const DRT::Element::DiscretizationType distype = DRT::Element::hex8;

  DRT::UTILS::shape_function_3D_deriv1(deriv1 ,xsi[0], xsi[1], xsi[2],distype);

  LINALG::Matrix<8,3> xrefe;  // material coord. of element
  int NUMNOD_SOH8 = 8;
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x =ele.Nodes()[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  fgrad.MultiplyNN(deriv1,xrefe);

  return true;
}
// Read in Stresses and Displacements from corresponding Run on lower Level
void STR::MLMC::ReadResultsFromLowerLevel()
{
  std::stringstream name_helper;
  // assamble filename and pathfor cluster jobs
  //name_helper << "../../level"<<num_level_-1<< "/START_RUN_"<< start_run_ <<"/"<<filename_lower_level_<<"_prolongated"<< "_run_" << numb_run_ ;

  name_helper << filename_lower_level_ <<"_prolongated"<< "_run_" << numb_run_ ;
  std::string name_ll = name_helper.str();
  // check if bool should be true or false
  RCP<IO::InputControl> inputcontrol = Teuchos::rcp(new IO::InputControl(name_ll, false));
  IO::DiscretizationReader input_coarse(actdis_fine_, inputcontrol,1);

  // read in displacements
  input_coarse.ReadMultiVector(disp_lower_level_, "displacement");
  // read in prolongated GP Stresses
  input_coarse.ReadMultiVector(stress_lower_level_, "prolongated_gauss_2PK_stresses_xyz");
  // read in prolongated GP Strain
   input_coarse.ReadMultiVector(strain_lower_level_, "prolongated_gauss_GL_strains_xyz");
}

void STR::MLMC::CalcDifferenceToLowerLevel(RCP< Epetra_MultiVector> stress, RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp)
{
  // store difference not in new vectors to save memory
  disp_lower_level_->Update(1.0,*disp,-1.0);
  stress_lower_level_->Update(1.0,*stress,-1.0);
  strain_lower_level_->Update(1.0,*strain,-1.0);

}

// Setup Material Parameters in each element based on deterministic value
void STR::MLMC::SetupStochMatDet(double value)
{
  // flag have init stochmat??
   int stochmat_flag=0;
   // Get parameters from stochastic matlaw
   const int myrank = discret_->Comm().MyPID();

   // loop all materials in problem
   const std::map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
   if (myrank == 0)
     IO::cout << "No. material laws considered: " << (int)mats.size() << IO::endl;
   std::map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
   for (curr=mats.begin(); curr != mats.end(); curr++)
   {
     const RCP<MAT::PAR::Material> actmat = curr->second;
     switch(actmat->Type())
     {
        case INPAR::MAT::m_aaaneohooke_stopro:
        {
          stochmat_flag=1;
          MAT::PAR::AAAneohooke_stopro* params = dynamic_cast<MAT::PAR::AAAneohooke_stopro*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
        }
        break;
       default:
       {
        IO::cout << "MAT CURR " << actmat->Type() << "not stochastic" << IO::endl;
        break;
       }

     }
   } // EOF loop over mats
   if (!stochmat_flag)// ignore unknown materials ?
   {
     dserror("No stochastic material supplied");
   }
   // loop over all elements
    for (int i=0; i< (discret_->NumMyColElements()); i++)
    {
      if(discret_->lColElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
      {
        MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->lColElement(i)->Material().get());
        aaa_stopro->Init(value,"BETA");
      }
    } // EOF loop elements

}
// Setup Material Parameters in each element based on Random Field
void STR::MLMC::SetupStochMat(unsigned int random_seed)
{
  // Variables for Random field
  double stoch_mat_par;
  // element center
  std::vector<double> ele_c_location;

  // flag have init stochmat??
  int stochmat_flag=0;
  // Get parameters from stochastic matlaw
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const std::map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if (myrank == 0)
    IO::cout << "No. material laws considered: " << (int)mats.size() << IO::endl;
  std::map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); curr++)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
       case INPAR::MAT::m_aaaneohooke_stopro:
       {
         stochmat_flag=1;
         MAT::PAR::AAAneohooke_stopro* params = dynamic_cast<MAT::PAR::AAAneohooke_stopro*>(actmat->Parameter());
         if (!params) dserror("Cannot cast material parameters");
       }
       break;
      default:
      {
       IO::cout << "MAT CURR " << actmat->Type() << "not stochastic" << IO::endl;
       break;
      }

    }
  } // EOF loop over mats
  if (!stochmat_flag)// ignore unknown materials ?
          {
            //IO::cout << "Mat Type   " << actmat->Type() << IO::endl;
            dserror("No stochastic material supplied");
          }
  // get elements on proc use col map to init ghost elements as well
  Teuchos::RCP<Epetra_Vector> my_ele = Teuchos::rcp(new Epetra_Vector(*discret_->ElementColMap(),true));
  // for doing quick deterministic computations instead of using a field use det_value instead
  // we do this so we can use our custom element output while doing deterministic simulation

  if (!use_det_value_)
  {
    if (numb_run_-start_run_== 0 )
    {
      random_field_ = Teuchos::rcp(new GenRandomField(random_seed,discret_));
      if (0)
      {
        // little testing loop
        // Get size
        //int size = random_field_->SizePerDim();
        //int  NumCosTerms = random_field_->NumberOfCosTerms();
        //int dim =random_field_->Dimension();
        //int mysize = pow(size,double(dim));
        //IO::cout << "stopped  calculation of sample psd " << IO::endl;
        //Teuchos::RCP<Teuchos::Array <double> > sample_psd= Teuchos::rcp( new Teuchos::Array<double>(mysize,0.0));
        //Teuchos::RCP<Teuchos::Array <double> > average_psd= Teuchos::rcp( new Teuchos::Array<double>(mysize,0.0));
        //random_field_->GetPSDFromSample3D(sample_psd);
        //random_field_->GetPSDFromSample(sample_psd);
        //random_field_->WriteSamplePSDToFile(sample_psd);
        // eof testing
        //dserror("stop here");
  //      Teuchos::RCP<Teuchos::Array <double> > average_psd=
  //          Teuchos::rcp( new Teuchos::Array<double>(size*size,0.0));
  //      Teuchos::RCP<Teuchos::Array <double> > sample_psd=
  //             Teuchos::rcp( new Teuchos::Array<double>(size*size,0.0));
  //      for (int i =0 ; i<500; i++)
  //      {
  //        random_field_->CreateNewSample(random_seed+i);
  //        std::vector<double> ele_center (3,1.0);
  //        // calc average psd
  //        //random_field_->GetPSDFromSample3D(sample_psd);
  //        random_field_->GetPSDFromSample(sample_psd);
  //        for (int j=0;j<size *size;j++)
  //        {
  //          (*average_psd)[j]+=1/500.*(*sample_psd)[j];
  //        }
  //        IO::cout << "Run " << i << " of 1000" << IO::endl;
  //      }
  //      ofstream File;
  //      File.open("Psd_average.txt",ios::app);
  //      for (int j=0;j<size*size;j++)
  //      {
  //          File << (*average_psd)[j] << IO::endl;
  //      }
  //      File.close();
  //      dserror("DONE stop now");
      }
    }
    else
    {
      random_field_->CreateNewSample(random_seed);
    }
  }
  // loop over all elements
  for (int i=0; i< (discret_->NumMyColElements()); i++)
  {
    if(discret_->lColElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->lColElement(i)->Material().get());
      std::vector<double> ele_center;
      ele_center = discret_->lColElement(i)->ElementCenterRefeCoords();
      if(use_det_value_)
      {
        stoch_mat_par=det_value_;
        IO::cout << "WARNING NOT USING RANDOM FIELD BUT A DET_VALUE OF " << det_value_<<"  INSTEAD" << IO::endl;
      }
      //if(numb_run_-start_run_== 0)
      //if(use_det_value_)
      //{
       // stoch_mat_par=4.61;
       // IO::cout << "WARNING NOT USING RANDOM FIELD BUT 4.61 INSTEAD" << IO::endl;
     // }
      else
      {
        // get dim of field
        if(random_field_->Dimension()==2)
        {
          //special hack here assuming circular geometry with r=25 mm
          double phi= acos(ele_center[0]/25);
          //compute x coord
          ele_center[0]=phi*25;
          ele_center[1]=ele_center[2];

          //IO::cout << "No Spherical Field" << IO::endl;
          //ele_center[1]=ele_center[2];
        }
        if (i==0)
          stoch_mat_par = random_field_->EvalFieldAtLocation(ele_center,false,true);
        else
          stoch_mat_par = random_field_->EvalFieldAtLocation(ele_center,false,false);
      }
        aaa_stopro->Init(stoch_mat_par,"BETA");

    }
  } // EOF loop elements

}
// Compute Stochmat parameter as linear combination
// beta= (1-gamma)beta_old+ gamma* beta_new, where beta_old is used for all elements
void STR::MLMC::BlendStochMat(unsigned int random_seed, bool reuse_rf_values, double beta_old, double gamma)
{
  // quick check wether we can reuse rf_values
  if(reuse_rf_values)
  {
    if(rf_values_->size() == 0)
      dserror("random field has not been initializes");
  }
  // element center
  std::vector<double> ele_c_location;

  // flag have init stochmat??
  int stochmat_flag=0;
  // Get parameters from stochastic matlaw
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const std::map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if (myrank == 0)
    IO::cout << "No. material laws considered: " << (int)mats.size() << IO::endl;
  std::map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); curr++)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
       case INPAR::MAT::m_aaaneohooke_stopro:
       {
         stochmat_flag=1;
         MAT::PAR::AAAneohooke_stopro* params = dynamic_cast<MAT::PAR::AAAneohooke_stopro*>(actmat->Parameter());
         if (!params) dserror("Cannot cast material parameters");
       }
       break;
      default:
      {
       IO::cout << "MAT CURR " << actmat->Type() << "not stochastic" << IO::endl;
       break;
      }

    }
  } // EOF loop over mats
  if (!stochmat_flag)// ignore unknown materials ?
          {
            //IO::cout << "Mat Type   " << actmat->Type() << IO::endl;
            dserror("No stochastic material supplied");
          }
  // get elements on proc use col map to init ghost elements as well

  Teuchos::RCP<Epetra_Vector> my_ele = Teuchos::rcp(new Epetra_Vector(*discret_->ElementColMap(),true));
  // for doing quick deterministic computations instead of using a field use det_value instead
  // we do this so we can use our custom element output while doing deterministic simulation
  if (!reuse_rf_values)
  {
    // first run uses deterministic value for all elements hence create rf here in first real run
    if (numb_run_-start_run_== 1 )
    {
      random_field_ = Teuchos::rcp(new GenRandomField(random_seed,discret_));
    }
    else
    {
      random_field_->CreateNewSample(random_seed);
    }
  }

  // loop over all elements
  for (int i=0; i< (discret_->NumMyColElements()); i++)
  {
    if(discret_->lColElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->lColElement(i)->Material().get());
      std::vector<double> ele_center;
      ele_center = discret_->lColElement(i)->ElementCenterRefeCoords();
      // get dim of field
      if(random_field_->Dimension()==2)
      {
        //special hack here assuming circular geometry with r=25 mm
        double phi= acos(ele_center[0]/25);
        //compute x coord
        ele_center[0]=phi*25;
        ele_center[1]=ele_center[2];
        //IO::cout << "No Spherical Field" << IO::endl;
        //ele_center[1]=ele_center[2];
      }
      // store all values in
      if(!reuse_rf_values)
      {
        rf_values_->at(i)=random_field_->EvalFieldAtLocation(ele_center,false,false);

      }
      // compute blended stopro parameter
      double blended = beta_old*(1-gamma) + (gamma)*rf_values_->at(i);
      aaa_stopro->Init(blended,"BETA");

    }
  } // EOF loop elements
}

// calculate some statistics
void STR::MLMC::CalcStatStressDisp(RCP< Epetra_MultiVector> curr_stress,RCP< Epetra_MultiVector> curr_strain,RCP<Epetra_MultiVector> curr_disp)
{
  // in order to avoid saving the stresses and displacements for each step an online
  // algorithm to compute the std deviation is needed .

  // Such an algorithm can be found here:

  // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

  // This it how it goes
  //def online_variance(data):
   //   n = 0
   //   mean = 0
   //   M2 = 0
   //
   //   for x in data:
    //      n = n + 1
    //      delta = x - mean
     //     mean = mean + delta/n
      //    M2 = M2 + delta*(x - mean)  # This expression uses the new value of mean
   //
     // variance_n = M2/n
     /// variance = M2/(n - 1)
    //  return variance*/
  // since numb_run_ does not start from zero anymore we need
  int n = numb_run_-start_run_+1;
  // calc mean and variance for displacement
  delta_disp_->Update(1.0,*curr_disp,-1.0,*mean_disp_,0.0);
  mean_disp_->Update(1.0/n,*delta_disp_,1.0);
  m2_helper_var_disp_->Update(1.0,*curr_disp,-1.0,*mean_disp_,0.0);
  m2_var_disp_->Multiply(1.0,*delta_disp_,*m2_helper_var_disp_ ,1.0);
  var_disp_->Update(1.0/(numb_run_),*m2_var_disp_,0.0);
  // do the same for stresses
  delta_stress_->Update(1.0,*curr_stress,-1.0,*mean_stress_,0.0);
  mean_stress_->Update(1.0/n,*delta_stress_,1.0);
  m2_helper_var_stress_->Update(1.0,*curr_stress,-1.0,*mean_stress_,0.0);
  m2_var_stress_->Multiply(1.0,*delta_stress_,*m2_helper_var_stress_ ,1.0);
  var_stress_->Update(1.0/(numb_run_),*m2_var_stress_,0.0);
  // and for strains
  delta_strain_->Update(1.0,*curr_strain,-1.0,*mean_strain_,0.0);
  mean_strain_->Update(1.0/n,*delta_strain_,1.0);
  m2_helper_var_strain_->Update(1.0,*curr_strain,-1.0,*mean_strain_,0.0);
  m2_var_strain_->Multiply(1.0,*delta_strain_,*m2_helper_var_strain_ ,1.0);
  var_strain_->Update(1.0/(numb_run_),*m2_var_strain_,0.0);

  // quick check if we need difference stats
  if (calc_diff_)
  {
  // calc mean and variance for difference between levels
  diff_delta_disp_->Update(1.0,*disp_lower_level_,-1.0,*diff_mean_disp_,0.0);
  diff_mean_disp_->Update(1.0/n,*diff_delta_disp_,1.0);
  diff_m2_helper_var_disp_->Update(1.0,*disp_lower_level_,-1.0,*diff_mean_disp_,0.0);
  diff_m2_var_disp_->Multiply(1.0,*diff_delta_disp_,*diff_m2_helper_var_disp_ ,1.0);
  diff_var_disp_->Update(1.0/(numb_run_),*diff_m2_var_disp_,0.0);
  // do the same for stresses
  diff_delta_stress_->Update(1.0,*stress_lower_level_,-1.0,*diff_mean_stress_,0.0);
  diff_mean_stress_->Update(1.0/n,*diff_delta_stress_,1.0);
  diff_m2_helper_var_stress_->Update(1.0,*stress_lower_level_,-1.0,*diff_mean_stress_,0.0);
  diff_m2_var_stress_->Multiply(1.0,*diff_delta_stress_,*diff_m2_helper_var_stress_ ,1.0);
  diff_var_stress_->Update(1.0/(numb_run_),*diff_m2_var_stress_,0.0);
  // do the same for strains
  diff_delta_strain_->Update(1.0,*strain_lower_level_,-1.0,*diff_mean_strain_,0.0);
  diff_mean_strain_->Update(1.0/n,*diff_delta_strain_,1.0);
  diff_m2_helper_var_strain_->Update(1.0,*strain_lower_level_,-1.0,*diff_mean_strain_,0.0);
  diff_m2_var_strain_->Multiply(1.0,*diff_delta_strain_,*diff_m2_helper_var_strain_ ,1.0);
  diff_var_strain_->Update(1.0/(numb_run_),*diff_m2_var_strain_,0.0);
  }



}
void STR::MLMC::WriteStatOutput()
{  //
  std::stringstream name_helper;
  name_helper << filename_ << "_statistics";
  output_fine_->NewResultFile(name_helper.str(),numb_run_);
  output_fine_->WriteMesh(1, 0.01);
  output_fine_->NewStep( 1, 0.01);
  output_fine_->WriteVector("mean_displacements", mean_disp_, output_fine_->dofvector);
  output_fine_->WriteVector("variance_displacements", var_disp_, output_fine_->dofvector);
  output_fine_->WriteVector("mean_gauss_2PK_stresses_xyz", mean_stress_, output_fine_->nodevector);
  output_fine_->WriteVector("variance_gauss_2PK_stresses_xyz", var_stress_, output_fine_->nodevector);
  output_fine_->WriteVector("mean_gauss_GL_strain_xyz", mean_strain_, output_fine_->nodevector);
  output_fine_->WriteVector("variance_gauss_GL_strain_xyz", var_strain_, output_fine_->nodevector);
  // write stats with respect to lower level
  if (calc_diff_)
  {
    output_fine_->WriteVector("diff_mean_displacements", diff_mean_disp_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_variance_displacements", diff_var_disp_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_mean_gauss_2PK_stresses_xyz", diff_mean_stress_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_variance_gauss_2PK_stresses_xyz", diff_var_stress_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_mean_gauss_GL_strain_xyz", diff_mean_strain_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_variance_gauss_GL_strain_xyz", diff_var_strain_, output_fine_->nodevector);
   }
}
void STR::MLMC::ResetPrestress()
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


void STR::MLMC::HelperFunctionOutput(RCP< Epetra_MultiVector> stress,RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp)
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
void STR::MLMC::HelperFunctionOutputTube(RCP< Epetra_MultiVector> stress,RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp)
{

  // assamble name for outputfil
  std::stringstream outputfile;
  outputfile << filename_ << "_statistics_output_" << start_run_ << ".txt";
  std::string name = outputfile.str();;
  /// file to write output
  // paraview ids of nodes
  int node[5] = {566, 1764, 3402,5194,6510};
  double disp_mag[5];
  double stress_mag[5];
  //double strain_mag[5];
  std::ofstream File;
  if (numb_run_ == 0 || numb_run_ == start_run_)
  {
    File.open(name.c_str(),std::ios::out);
    if (File.is_open())
    {
      //if(error == -1)
      //dserror("Unable to open Statistics output Filename");
      File << "run id   "<< "disp mag node   " << node[0] << "S_mag node   " << node[0] << "disp mag node   " << node[1] << "S_mag node   " << node[1]
       << "disp mag node   " << node[2] << "S_mag node   " << node[2] << "disp mag node   " << node[3] << "S_mag node   "
       << node[3] << "disp mag node   " << node[4] << "S_mag node   " << node[4]<< std::endl;
      File.close();
    }
    //
    else
    {
    dserror("Unable to open statistics output file");
    }
  }
  for(int i =0; i<5 ; i++)
  {
    //calc mag disp
    disp_mag[i]=sqrt(pow((*disp)[0][node[i]*3],2)+pow((*disp)[0][node[i]*3+1],2)+pow((*disp)[0][node[i]*3+2],2));
    stress_mag[i]=sqrt(pow((*stress)[0][node[i]],2)+pow((*stress)[1][node[i]],2)+ pow((*stress)[2][node[i]],2)+ pow((*stress)[3][node[i]],2)+
        pow((*stress)[4][node[i]],2)+pow((*stress)[5][node[i]],2));

  }
  // reopen in append mode
  File.open(name.c_str(),std::ios::app);
  File << numb_run_ << "    "<< disp_mag[0]<< "    " << stress_mag[0] << "    "<< disp_mag[1]<< "    " << stress_mag[1] << "    "<< disp_mag[2]<< "    " << stress_mag[2]
    << "    "<< disp_mag[3]<< "    " << stress_mag[3] << "    "<< disp_mag[4]<< "    " << stress_mag[4] <<  std::endl;
  File.close();
}

void STR::MLMC::EvalDisAtNodes(Teuchos::RCP<const Epetra_Vector> disp )
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

  RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  // disp is passed to the function no need for reading in results anymore
  actdis_coarse_->SetState("displacement",disp);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alrigth lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");
  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
 p.set("gpstrainmap", gpstrainmap);


  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();


  // st action to calc poststresse
  p.set("action","postprocess_stress");
  // Multivector to store poststresses
  // WE need the collumn Map here because we need all the ghosted nodes to calculate stresses at the nodes
  RCP<Epetra_MultiVector> poststress =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeColMap()),6,true));

  p.set("gpstressmap", gpstressmap);
  p.set("poststress", poststress);
  p.set("stresstype","ndxyz");

  //actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();

  // again for strains
  p.set("action","postprocess_stress");
  // Multivector to store poststrains
  RCP<Epetra_MultiVector> poststrain =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeColMap()),6,true));
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
  RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
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
  RCP<Epetra_MultiVector> elestrains = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
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

  RCP<Epetra_Vector> output_disp = LINALG::CreateVector(output_dof_map,false);
  RCP<Epetra_MultiVector> output_stress =  Teuchos::rcp(new Epetra_MultiVector(output_node_map,6,true));
  RCP<Epetra_MultiVector> output_strain =  Teuchos::rcp(new Epetra_MultiVector(output_node_map,6,true));
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


void STR::MLMC::SetupEvalDisAtEleCenters(std::vector <int> AllOutputEleIds)
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

void STR::MLMC::SetupEvalPeakWallStress()
{
  int  num_my_wall_elements =0;
 // Now loop over all elements in EleRowMap

  for(int i=0;i<actdis_coarse_->ElementRowMap()->NumMyElements();i++)
  {
    // check wether element has AAANeohookeStopro Material as indicator for AAA wall
    if(actdis_coarse_->lRowElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      num_my_wall_elements++;
      my_wall_elements_.push_back(actdis_coarse_->lRowElement(i)->Id());
    }
  }//EOF Loop over Elements
  // Broadcast the information of how many wall elements we have to all other procs
  actdis_coarse_->Comm().SumAll(&num_my_wall_elements,&tot_num_wall_elements_,1);
}

void STR::MLMC::EvalPeakWallStress(Teuchos::RCP<const Epetra_Vector> disp, INPAR::STR::StressType iostress,INPAR::STR::StrainType iostrain )
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

  CalcVonMises(my_output_elements_c_stresses,vonMisesStress);
  CalcVonMises(my_output_elements_c_strains,vonMisesStrain);

  // compute disp magnitude
  for(unsigned int i=0;i<my_output_elements_c_disp->size();i++)
  {
    DispMag->push_back(sqrt( pow(my_output_elements_c_disp->at(i).at(0),2)+pow(my_output_elements_c_disp->at(i).at(1),2)+pow(my_output_elements_c_disp->at(i).at(2),2) ));
    // since my_output_elements_mat_params is the same size we use the same for loop here
    Beta->push_back(my_output_elements_mat_params->at(i).at(1));
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

  ComputePeakAndQuantile(Beta,&my_wall_elements_,temp_peak,temp_quantile99);
  peak_vals->push_back(*temp_peak);
  quantile99_vals->push_back(*temp_quantile99);

  ExportPeakStressDataAndWriteToFile(peak_vals,quantile99_vals);

}

void STR::MLMC::CalcVonMises(Teuchos::RCP<std::vector <std::vector<double > > > input_components , Teuchos::RCP< std::vector <double> > output_vM )
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

void STR::MLMC::ComputePeakAndQuantile(Teuchos::RCP<std::vector <double > >  values , std::vector <int>* ele_ids, Teuchos::RCP<std::pair<int, double> > peak,Teuchos::RCP<std::pair<int, double> > quantile99)
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

   // write to screen
   //if(!actdis_coarse_->Comm().MyPID())
   //{
   //for (std::list< std::pair< int, double > >::iterator it = my_quantity_gathered.begin(); it != my_quantity_gathered.end(); it++)
     //{
        // IO::cout << " " << it->second << " "<< it->first << IO::endl;
     //}
   //}
   std::list< std::pair< int ,double> >::iterator it = my_quantity_gathered.begin();
   std::advance( it,(int)(tot_num_wall_elements_*0.01));
   quantile99->first=it->first;
   quantile99->second=it->second;
}


void STR::MLMC::EvalDisAtEleCenters(Teuchos::RCP<const Epetra_Vector> disp,
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

    std::vector<double> mat_params ;
    // get the mat parameters
    if(actdis_coarse_->gElement(output_elements->at(i))->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
    {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->gElement(output_elements->at(i))->Material().get());

      mat_params.push_back(aaa_stopro->Youngs());
      mat_params.push_back(aaa_stopro->Beta());
    }
    else
    {
      mat_params.push_back(0.0);
      mat_params.push_back(0.0);
    }
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

  RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  // disp is passed to the function no need for reading in results anymore
  actdis_coarse_->SetState("displacement",disp);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alright lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");
  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstrainmap", gpstrainmap);
  //actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Evaluate2(p,output_elements,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();
  // also get element stresses
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", gpstressmap);
  RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
  p.set("poststress",elestress);
  Evaluate2(p,output_elements,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  //actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }
  // and strains
  p.set("gpstressmap", gpstrainmap);
  RCP<Epetra_MultiVector> elestrains = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->ElementRowMap()),6));
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

void STR::MLMC::ExportPeakStressDataAndWriteToFile(
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
        File << "runid "<< "vMstress EID vMstrain EID Phi EID disp EID beta EID vMstress99 EID vMstrain99 EID Phi99 EID disp99 EID beta99 EID "
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

void STR::MLMC::ExportEleDataAndWriteToFile(Teuchos::RCP<const Epetra_Map> OutputMap,
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
    for(int i=0;i<2;i++)
    {
       File << " " << (*mat_params)[i];
    }
    File << std::endl;
    File.close();
    }
    //dsfdsf
  }
  actdis_coarse_->Comm().Barrier();
}
void STR::MLMC::WriteParamContInfoToFile(
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
void STR::MLMC::Evaluate2(
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

