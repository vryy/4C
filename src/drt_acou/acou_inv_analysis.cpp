/*!----------------------------------------------------------------------
\file acou_inv_analysis.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "acou_inv_analysis.H"
#include "acou_impl_euler.H"
#include "acou_impl_trap.H"
#include "acou_impl_dirk.H"
#include "acou_impl_bdf.H"
#include "acou_inv_resulttest.H"
#include "acou_ele_action.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/material.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "Epetra_SerialDenseSolver.h"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
ACOU::InvAnalysis::InvAnalysis(Teuchos::RCP<DRT::Discretization> scatradis,
                               Teuchos::RCP<DRT::DiscretizationHDG> acoudis,
                               const Teuchos::ParameterList& scapara,
                               const Teuchos::ParameterList& scasolvpara,
                               Teuchos::RCP<Teuchos::ParameterList> acoupara,
                               Teuchos::RCP<LINALG::Solver> acousolv,
                               Teuchos::RCP<IO::DiscretizationWriter> acouout
)
  : scatra_discret_(scatradis),
    acou_discret_(acoudis),
    scatraalgo_(Teuchos::null),
    acoualgo_(Teuchos::null),
    scatraparams_(scapara),
    scatrasolverparams_(scasolvpara),
    acouparams_(acoupara),
    acousolver_(acousolv),
    acououtput_(acouout),
    dyna_(DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*acouparams_,"TIMEINT")),
    phys_(DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*acouparams_,"PHYSICAL_TYPE")),
    myrank_(acoudis->Comm().MyPID()),
    error_(1.0e6),
    tol_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL")),
    iter_(0),
    max_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_MAX_RUN")),
    max_ls_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_LS_MAX_RUN")),
    output_count_(0),
    fdcheck_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("FDCHECK"))),
    alpha_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("ALPHA_MUA")),
    beta_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("BETA_MUA")),
    calcacougrad_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("INV_TOL_GRAD_YN"))),
    nm_(0),
    dtacou_(acouparams_->get<double>("TIMESTEP")),
    J_(0.0),
    normdiffp_(0.0),
    tstart_(Teuchos::Time::wallTime())
{
  if(phys_==INPAR::ACOU::acou_viscous) dserror("inverse analysis for now only implemented for lossless fluid");

  scatra_output_ = scatra_discret_->Writer();

  // get the output name
  name_ = DRT::Problem::Instance()->OutputControlFile()->FileName();

  // vectors we need
  adjoint_phi_0_ = LINALG::CreateVector(*(acou_discret_->NodeRowMap()),true);
  phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  adjoint_w_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  if(calcacougrad_)
    tol_grad_ = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL_GRAD");
  else
    tol_grad_ = 0.0;

  // initialize acou_rhs_: we need a vector with the nodes of the boundary where
  // the pressure is monitored-> Read the monitor file and create a vector with
  // corresponding nodes OR take the boundary where absorbing bcs are prescribed!
  // we need this map extractor thing!
  // we deal with NODES here, not with DOFS

  std::string condname = "Absorbing";
  std::vector<DRT::Condition*> absorbingBC;
  acou_discret_->GetCondition(condname,absorbingBC);
  if(absorbingBC.size()==0)
    dserror("you have to use absorbing line conditions for inverse analysis!");
  const std::vector<int> abcnodes = *(absorbingBC[0]->Nodes());

  abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, abcnodes.size(), &abcnodes[0], 0, acou_discret_->Comm()));
  abcnodes_mapex_ = Teuchos::rcp(new LINALG::MapExtractor(*(acou_discret_->NodeRowMap()),abcnodes_map_,true));

  // determine the number of vectors for monitoring
  // this is naive: later on, we won't be able to store everything at once and we have to implement
  // a smarter approach to reduce storage requirements
  int numvec = acouparams_->get<int>("NUMSTEP");
  int oderso = acouparams_->get<double>("MAXTIME")/dtacou_;
  if ( oderso < numvec)
    numvec = oderso+1;
  t_ = numvec * dtacou_;

  acou_rhs_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec));
  acou_rhsm_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec));

  // get measured values
  // open monitor file and read it
  nnodes_ = 0;
  {
    char* foundit = NULL;
    std::string monitorfilename = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE");
    if (monitorfilename=="none.monitor") dserror("No monitor file provided");
    // insert path to monitor file if necessary
    if (monitorfilename[0]!='/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos!=std::string::npos)
      {
        std::string path = filename.substr(0,pos+1);
        monitorfilename.insert(monitorfilename.begin(), path.begin(), path.end());
      }
    }

    FILE* file = fopen(monitorfilename.c_str(),"rb");
    if (file==NULL) dserror("Could not open monitor file %s",monitorfilename.c_str());

    char buffer[150000];
    fgets(buffer,150000,file);
    // read steps
    foundit = strstr(buffer,"steps"); foundit += strlen("steps");
    nsteps_ = strtol(foundit,&foundit,10);
    timesteps_.resize(nsteps_);
    // read nnodes
    foundit = strstr(buffer,"nnodes"); foundit += strlen("nnodes");
    nnodes_ = strtol(foundit,&foundit,10);
    // read nodes
    nodes_.resize(nnodes_);
    for (unsigned int i=0; i<nnodes_; ++i)
    {
      fgets(buffer,150000,file);
      foundit = buffer;
      nodes_[i] = strtol(foundit,&foundit,10) - 1;

      //if (!myrank_) printf("Monitored node %d ",nodes_[i]);
      //if (!myrank_) printf("\n");
    }
    // read in measured curve
    {
      mcurve_ = Epetra_SerialDenseVector(nnodes_*nsteps_);

      //if (!myrank_) printf("nsteps %d nnode %d\n",nsteps_,nnodes_);

      // read comment lines
      foundit = buffer;
      fgets(buffer,150000,file);
      while(strstr(buffer,"#"))
        fgets(buffer,150000,file);

      // read in the values for each node
      unsigned int count = 0;
      for (int i=0; i<nsteps_; ++i)
      {
        // read the time step
        timesteps_[i] = strtod(foundit,&foundit);
        for (unsigned int j=0; j<nnodes_; ++j)
          mcurve_[count++] = strtod(foundit,&foundit);
        fgets(buffer,150000,file);
        foundit = buffer;
      }
      if (count != nnodes_*nsteps_) dserror("Number of measured pressure values wrong on input");
    }
  }

  if (nnodes_ != abcnodes.size()) dserror("For now implemented only when all boundary nodes are monitored");
  if (nodes_ != abcnodes) dserror("And please provide the correct order (feel free to reimplement)");

  // every proc knows mcurve_, now, we want to write mcurve_ to a Epetra_MultiVector in the same form as acou_rhs_
  // with the same parallel distribution!
  // and we want to interpolate measured values in case the monitored time step size is not the same as the one for the simulation
  double eps = dtacou_/1000.0;


  if( timesteps_[0] != 0.0 )
    dserror("your measured values have to start at time 0.0");
  else if( timesteps_[0] == 0.0 && timesteps_[1] == dtacou_ ) // the standard case
    for(unsigned int i=0; i<nnodes_; ++i)
    {
      if( acou_discret_->HaveGlobalNode(nodes_[i]) )
        for(int j=0; j<nsteps_; ++j)
          acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,mcurve_(i+j*nnodes_)); // the proc who has this row, writes the value
    }
  else // we have to interpolate!
  {
    if( numvec < nsteps_ )
    {
      dserror("set your time step size smaller, at least to %14f or implement here",timesteps_[1]-timesteps_[0]);
    }
    else
    {

      for(unsigned int i=0; i<nnodes_; ++i)
      {
        if( acou_discret_->HaveGlobalNode(nodes_[i]) )
        {
          for(int j=0; j<numvec; ++j)
          {
            double actualt = j * dtacou_; // we need values for this time
            int timeval = 0;
            // find next higher and next lower value
            while(actualt>timesteps_[timeval]-eps)
            {
              timeval++;
            }
            // timesteps_[timeval] has the next higher point in time
            // now interpolate from this and the value before
            if(timeval == 0)
            {
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,0.0);
            }
            else if(actualt<timesteps_[timeval]+eps && actualt>timesteps_[timeval]-eps) // then this is more or less it
            {
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,mcurve_(i+timeval*nnodes_));
            }
            else
            {
              double value = mcurve_(i+(timeval-1)*nnodes_) + (mcurve_(i+(timeval)*nnodes_)-mcurve_(i+(timeval-1)*nnodes_)) * (actualt - timesteps_[timeval-1]) / (timesteps_[timeval]-timesteps_[timeval-1]);
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,value);
            }
          } // for(int j=0; j<numvec; ++j)
        } // if( acou_discret_->HaveGlobalNode(nodes_[i]) )
      } // for(unsigned int i=0; i<nnodes_; ++i)
    } // else ** if( numvec < nsteps_ )
  } // else ** if( timesteps_[0] == dtacou_ || (timesteps_[0]==0.0 && timesteps_[1] = dtacou_) )

  // we also need vectors containing node based material parameters
  // the build process of these vectors needs communication and we only want to do this once,
  // that's why it's done in the constructor!

  // allocate vectors
  node_c_    = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);
  node_rho_  = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);
  node_mu_a_ = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);

  // now: fill the vectors
  int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd) // cannot loop scatra nodes, because they don't necessarily start with gid 0
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou) )
    {
      opti_node = scatra_discret_->gNode(nd+maxnodeidacou);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;
    // here, every proc knows the owner and the gid of the optical node

    // get the corresponding acoustical node
    DRT::Node* acou_node = NULL;
    int myacounodeowner = -1;
    int acounodeowner = -1;

    if ( acou_discret_->HaveGlobalNode(nd) )
    {
      acou_node = acou_discret_->gNode(nd);
      myacounodeowner = acou_node->Owner();
      if ( myacounodeowner != myrank_ ) myacounodeowner = -1;
    }
    acou_discret_->Comm().MaxAll(&myacounodeowner,&acounodeowner,1);

    // the acoustical node shall calculate its speed of sound and density from neighboring elements
    int loc_numacouele = 0;
    double loc_c = 0.0;
    double loc_rho = 0.0;
    // every proc checks its elements, whether they are adjacent to this node, if so, add values to c, rho and numacouele
    for(int rael = 0; rael<acou_discret_->NumMyRowElements(); ++rael)
    {
      DRT::Element* racouele = acou_discret_->lRowElement(rael);
      const int* nodeids = racouele->NodeIds();
      int numnode = racouele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd )
        {
          const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(racouele->Material().get());
          loc_c += actmat->SpeedofSound();
          loc_rho += actmat->Density();
          loc_numacouele++;
        }
    }
    int glo_numacouele = 0;
    double glo_c = 0.0;
    double glo_rho = 0.0;
    // communicate to global values
    acou_discret_->Comm().SumAll(&loc_numacouele,&glo_numacouele,1);
    acou_discret_->Comm().SumAll(&loc_c,&glo_c,1);
    acou_discret_->Comm().SumAll(&loc_rho,&glo_rho,1);
    glo_c /= double(glo_numacouele);
    glo_rho /= double(glo_numacouele);

    // now, every proc knows the speed of sound and density of this node -> write them to vector
    int nodelid = scatra_discret_->NodeRowMap()->LID(nd+maxnodeidacou);
    if( nodelid >= 0 ) // only on owning proc
    {
      node_c_->ReplaceMyValue(nodelid,0,glo_c);
      node_rho_->ReplaceMyValue(nodelid,0,glo_rho);
    }

    // we have to do the same procedure for the absorption coefficient with the scatra discretization
    int loc_numoptiele = 0;
    double loc_mu_a = 0.0;
    for(int roel = 0; roel<scatra_discret_->NumMyRowElements(); ++roel)
    {
      DRT::Element* roptele = scatra_discret_->lRowElement(roel);
      const int* nodeids = roptele->NodeIds();
      int numnode = roptele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd+maxnodeidacou )
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(roptele->Material().get());
          loc_mu_a += actmat->ReaCoeff();
          loc_numoptiele++;
        }
    }
    int glo_numoptiele = 0;
    double glo_mu_a = 0.0;
    scatra_discret_->Comm().SumAll(&loc_numoptiele,&glo_numoptiele,1);
    scatra_discret_->Comm().SumAll(&loc_mu_a,&glo_mu_a,1);
    glo_mu_a /= double(glo_numoptiele);

    // write this value to vecot
    if( nodelid >= 0 ) // only on owning proc
      node_mu_a_->ReplaceMyValue(nodelid,0,glo_mu_a);
  } // for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)

  // cout<<"map"<<endl;
  // abcnodes_map_->Print(cout);
/*  cout<<"acou rhsm"<<endl;
  acou_rhsm_->Print(cout);
  for(int i=0; i<acou_rhsm_->NumVectors(); ++i)
  {
    cout<< double(i)*dtacou_ <<" "<<acou_rhsm_->operator ()(i)->operator [](0)<<" "<<acou_rhsm_->operator ()(i)->operator [](1)<<" "<<acou_rhsm_->operator ()(i)->operator [](2)<<" "<<acou_rhsm_->operator ()(i)->operator [](3)<<endl;
  }
*/

  // read material parameters from input file
  ReadInParameters();

  // number of material parameters
  np_ = p_.Length();
  G_.Resize(np_);
  G_.Scale(0.0);
  H_.Reshape(np_,np_);
  H_.Scale(0.0);
  // initialize the BFGS approximation of the Hessian to identity (there are smarter possibilities, feel free to implement)
  for(int i=0; i<np_; ++i)
    H_(i,i) = 1.0;

} // ACOU::InvAnalysis::InvAnalysis(...)

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::Integrate()
{
  // Solve the standard problem
  SolveStandardProblem();

  // fitting loop
  do
  {
    if(!myrank_)
    {
      std::cout<<std::endl;
      std::cout<<"*********************************************************************************"<<std::endl;
      std::cout<<"iteration "<<iter_+1<<" of maximal "<<max_iter_<<" iterations "<<std::endl;
      std::cout<<"*********************************************************************************"<<std::endl;
      std::cout<<std::endl;
    }

    // Calculate the error and the value of the objective function
    CalculateObjectiveFunctionValue();

    // Solve the adjoint problem
    SolveAdjointAcouProblem();
    SolveAdjointOptiProblem();

    // Build the gradient
    CalculateGradient();

    // Update the approximation to the inverse Hessian
    UpdateHessian();

    // Update the sought parameters
    UpdateParameters();

    // Tell all the elements their new material parameters
    UpdateMaterial();

    // Calculate some useful numbers
    CalculateStatsAndService();

    // Output some useful user information, like time consume, solution advance, ...
    OutputStats();

  } while ( ConvergenceCheck() );

  if( error_<=tol_ && !myrank_ )
  {
    std::cout<<std::endl;
    std::cout<<"*********************************************************************************"<<std::endl;
    std::cout<<"optimization terminated successfully, objective function fulfills given tolerance"<<std::endl;
    std::cout<<"*********************************************************************************"<<std::endl;
    std::cout<<std::endl;
  }
  else if(iter_ == max_iter_ && !myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"<<std::endl;
    std::cout<<"optimization terminated unsuccessfully, maximal number of iterations reached"<<std::endl;
    std::cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"<<std::endl;
    std::cout<<std::endl;
  }
  else if(G_.Norm2()<= tol_grad_ && !myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*"<<std::endl;
    std::cout<<"optimization terminated successfully, norm of gradient fulfills given tolerance"<<std::endl;
    std::cout<<"*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*"<<std::endl;
    std::cout<<std::endl;
  }
  else if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"<<std::endl;
    std::cout<<"optimization terminated, norm of differnce of parameters (almost) zero (<1e-14)"<<std::endl;
    std::cout<<"*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*"<<std::endl;
    std::cout<<std::endl;
  }

  return;
} // void ACOU::InvAnalysis::Integrate()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::SolveStandardProblem()
{
  // create scatra algorithm
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams_,"VELOCITYFIELD");

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(scatrasolverparams_,scatra_discret_->Comm(),DRT::Problem::Instance()->ErrorFile()->Handle()));
  Teuchos::RCP<Teuchos::ParameterList> scatratimeparams= Teuchos::rcp(new Teuchos::ParameterList(scatraparams_));
  scatratimeparams->set<double>   ("TIMESTEP"    ,1.0);
  scatratimeparams->set<double>   ("MAXTIME"     ,1.0);
  scatratimeparams->set<int>      ("NUMSTEP"     ,1);
  scatratimeparams->set           ("RESTARTEVRY" ,1);
  scatratimeparams->set           ("UPRES"       ,1);

  Teuchos::RCP<Teuchos::ParameterList> extraparams = Teuchos::rcp(new Teuchos::ParameterList());
  extraparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
  extraparams->set<bool>("isale",false);
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  extraparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  extraparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
  extraparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  extraparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

  scatra_discret_->Comm().Barrier();
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SCALAR TRANSPORT PROBLEM - OPTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatra_discret_->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      std::string outname = name_;
      outname.append("_invforward_opti");
      scatra_output_->NewResultFile(outname,output_count_);
      output_count_++;
      scatra_output_->WriteMesh(0,0.0);

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      scatraalgo_ = Teuchos::rcp(new SCATRA::TimIntStationary(scatra_discret_, solver, scatratimeparams, extraparams, scatra_output_));

      scatraalgo_->Init();
      scatraalgo_->SetVelocityField();

      scatraalgo_->TimeLoop();

      // output of elemental reaction coefficient
      scatraalgo_->DiscWriter().WriteVector("rea_coeff",ElementMatVec());

      phi_ = scatraalgo_->Phinp();

      break;
    }
    default:
      dserror("unknown velocity field type for transport of passive scalar in problem type Acoustics");
      break;
  }
  bool meshconform = DRT::INPUT::IntegralValue<bool>(*acouparams_,"MESHCONFORM");

  double pulse = acouparams_->get<double>("PULSEDURATION");
  acouparams_->set<bool>("adjoint",false);

  std::string outname = name_;
  outname.append("_invforward_acou");
  acououtput_->NewResultFile(outname,output_count_);
  output_count_++;

  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_trapezoidal:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplTrap(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_bdf2:
  case INPAR::ACOU::acou_bdf3:
  case INPAR::ACOU::acou_bdf4:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplBDF(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_dirk23:
  case INPAR::ACOU::acou_dirk33:
  case INPAR::ACOU::acou_dirk34:
  case INPAR::ACOU::acou_dirk54:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplDIRK(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }
  acoualgo_->SetInitialPhotoAcousticField(pulse,phi_,scatra_discret_, meshconform);

  // we have to call a slightly changed routine, which fills our history vector which we need for the adjoint problem
  // acoualgo_->Integrate();
  acoualgo_->Integrate(acou_rhs_,abcnodes_mapex_);

//  if(iter_==0 && output_count_<3) // do you want this output?
//  {
//    std::cout<<"measured values"<<std::endl;
//    std::cout.precision(15);
//    for(int j=0; j<acou_rhs_->NumVectors(); ++j)
//    {
//      std::cout<<dtacou_*double(j)<<" ";
//      for(int i=0; i<acou_rhs_->Map().NumGlobalElements(); ++i)
//        std::cout<<" "<<acou_rhs_->operator ()(j)->operator [](i);
//      std::cout<<std::endl;
//    }
//  }

  return;
} // void ACOU::InvAnalysis::SolveStandardProblem()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateObjectiveFunctionValue()
{
  //      alpha  ||      ||2      alpha  ||     ||2       1   ||         ||2
  // J = ------- || mu_a ||    + ------- ||  D  ||    + ----- ||  p - p  ||
  //        2    ||      ||L2O      2    ||     ||L2O     2   ||       m ||L2O

  // first two terms analog to the functions in the calculation of the gradient
  J_ = 0.0;
  INPAR::SCATRA::ScaTraType scatype = DRT::INPUT::IntegralValue<INPAR::SCATRA::ScaTraType>(scatraparams_,"SCATRATYPE");
  {
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("scatratype",scatype);
    eleparams.set<int>("action",SCATRA::calc_integr_objf);

    eleparams.set<double>("integr_objf",0.0);

    scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    double loc_integr_objf = eleparams.get<double>("integr_objf");
    double glo_integr_objf = 0.0;

    scatra_discret_->Comm().SumAll(&loc_integr_objf,&glo_integr_objf,1);

    J_ += glo_integr_objf;
    J_ *= alpha_; // regularization parameter (0.5 already in element routine)
  }

  Teuchos::RCP<Epetra_Vector> node_mu_a_col = LINALG::CreateVector(*(scatra_discret_->NodeColMap()),true);
  LINALG::Export(*node_mu_a_,*node_mu_a_col);

  Teuchos::RCP<Epetra_Vector> rea_grad_ele = LINALG::CreateVector(*(scatra_discret_->ElementRowMap()),true);
  // now, we do not have to do all communication manually!
  for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = scatra_discret_->lRowElement(i);
    double rea_ele = static_cast<const MAT::ScatraMat*>(ele->Material().get())->ReaCoeff();
    const int* nodeids = ele->NodeIds();
    int numnode = ele->NumNode();
    double grad = rea_ele;
    for(int j=0; j<numnode; ++j)
    {
      double loc_rea_node = node_mu_a_col->operator [](scatra_discret_->NodeColMap()->LID(nodeids[j]));
      grad -= loc_rea_node/double(numnode);
    }
    rea_grad_ele->operator [](i) = (grad);
  }
  double reg_grad = 0.0;
  rea_grad_ele->Dot(*rea_grad_ele,&reg_grad);
  J_ += reg_grad * beta_ / 2.0;


  // now the last term  0.5 || p - p_m ||^2_L2G
  {
    Teuchos::RCP<Epetra_MultiVector> temp = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,acou_rhs_->NumVectors()));;
    temp->Update(1.0,*acou_rhs_,0.0);
    temp->Update(-1.0,*acou_rhsm_,1.0);

    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",ACOU::calc_integr_objf);
    eleparams.set<double>("integr_objf",0.0);
    eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

    eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",temp);
    eleparams.set<double>("dt",dtacou_);

    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    acou_discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
      acou_discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,condname);

    double loc_integr_objf = eleparams.get<double>("integr_objf");
    double glo_integr_objf = 0.0;

    acou_discret_->Comm().SumAll(&loc_integr_objf,&glo_integr_objf,1);

    error_ = glo_integr_objf;

    error_ /= 2.0 ; //* t_;
    J_ += error_;
  }

  if(!myrank_)
  {
    std::cout.precision(15);
    std::cout<<"objective function value "<<J_<<std::endl;
  }
  return;
} // void ACOU::InvAnalysis::CalculateObjectiveFunctionValue()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::SolveAdjointAcouProblem()
{
  // again, we need the acoustical time integrator

  // the integration is exactly the same (due to some smart rearrangements of variables)
  // but we had to set the flag such that the integrator accounts for the right hand side values
  // from the difference of measured and simulated values, also these values have to be fed into
  // the parameter list

  acouparams_->set<bool>("adjoint",true); // but this time, we set the flag that we are doing the adjoint problem
  // set list of monitored nodes
  Teuchos::RCP<std::vector<int> > nodes_rcp= Teuchos::rcp(new std::vector<int> (nodes_.size()));
  for(unsigned int i=0; i<nodes_.size(); ++i)
    (*nodes_rcp)[i] = nodes_[i];
  acouparams_->set<Teuchos::RCP<std::vector<int> > >("monitorednodes",nodes_rcp);

  // set p-p_m
/*  cout<<"measured"<<endl;
  acou_rhsm_->Print(cout);
  cout<<"simulated"<<endl;
  acou_rhs_->Print(cout);
  cout<<"simulated"<<endl;
  for(int i=0; i<acou_rhs_->NumVectors(); ++i)
  {
    cout<< double(i)*dtacou_ <<" "<<acou_rhs_->operator ()(i)->operator [](0)<<" "<<acou_rhs_->operator ()(i)->operator [](1)<<" "<<acou_rhs_->operator ()(i)->operator [](2)<<" "<<acou_rhs_->operator ()(i)->operator [](3)<<endl;
  }
*/
  acou_rhs_->Update(-1.0,*acou_rhsm_,1.0);
  // acou_rhs_->Scale(1.0/t_);

  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",acou_rhs_);

  std::string outname = name_;
  outname.append("_invadjoint_acou");
  acououtput_->NewResultFile(outname,output_count_);
  output_count_++;

  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_trapezoidal:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplTrap(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_bdf2:
  case INPAR::ACOU::acou_bdf3:
  case INPAR::ACOU::acou_bdf4:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplBDF(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_dirk23:
  case INPAR::ACOU::acou_dirk33:
  case INPAR::ACOU::acou_dirk34:
  case INPAR::ACOU::acou_dirk54:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplDIRK(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }

  // here the initial field is zero everywhere
  acoualgo_->SetInitialZeroField();

  // integrate the adjoint problem
  acoualgo_->Integrate();

  // give me phi(0) which is needed for the source term of the adjoint optical problem
  adjoint_phi_0_->PutScalar(0.0);
  acoualgo_->NodalPressureField(adjoint_phi_0_);

  return;
} // void ACOU::InvAnalysis::SolveAdjointAcouProblem()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::SolveAdjointOptiProblem()
{
  // what we do here:
  // 1. create the source term or rhs vector for the scatra problem
  // 2. solve the scatra problem
  // 3. store the solution of the scatra problem

  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SCALAR TRANSPORT PROBLEM - ADJOINT OPTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  Teuchos::RCP<LINALG::SparseMatrix> sysmatscatra = scatraalgo_->SystemMatrix(); // this matrix and the algorithm should still exist

  // also, the system matrix has already the applied Dirichlets (1.0 on diagonals)
  // so, we can just use the systemmatrix and our own right hand side vector
  // which is built from the adjoint acoustical problem

  // calculate this rhs vector!!!
  // acou node based vector -> opti dof based vector
  Teuchos::RCP<Epetra_Vector> rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  std::string condname = "Dirichlet";
  std::vector<DRT::Condition*> dirichlets;
  scatra_discret_->GetCondition(condname,dirichlets);

  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd) )
    {
      opti_node = scatra_discret_->gNode(nd);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;

    double loc_value = 0.0;
    if(acou_discret_->NodeRowMap()->LID(nd)>-1)
    {
      loc_value = adjoint_phi_0_->operator [](acou_discret_->NodeRowMap()->LID(nd));
    }
    double glo_value = 0.0;
    acou_discret_->Comm().SumAll(&loc_value,&glo_value,1);

    // ok, we got the value, we still need c, rho and mu_a, but they are stored on the nodemap of the scatra dis, so this should not be a problem
    if(scatra_discret_->Comm().MyPID() == optnodeowner)
    {
      int lnodeid = node_mu_a_->Map().LID(nd);
      double mu_a = node_mu_a_->operator [](lnodeid);
      double c    = node_c_->operator [](lnodeid);
      int dofgid = scatra_discret_->Dof(opti_node,0);
      int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
      int err = rhsvec->ReplaceMyValue(doflid,0,-glo_value/c/c*mu_a);
      if (err) dserror("could not replace local vector entry");
    }
  }

  //this vector contains dof based values for the right hand side vector
  // to get a proper rhs vec we have to integrate it!!!!!!!
  Teuchos::ParameterList eleparams;

  INPAR::SCATRA::ScaTraType scatype = DRT::INPUT::IntegralValue<INPAR::SCATRA::ScaTraType>(scatraparams_,"SCATRATYPE");
  eleparams.set<int>("scatratype",scatype);
  scatra_discret_->SetState("rhsnodebasedvals",rhsvec);
  eleparams.set<int>("action",SCATRA::calc_integr_pat_rhsvec);
  Teuchos::RCP<Epetra_Vector> b = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,b,Teuchos::null,Teuchos::null);

  for(int nd=0; nd<scatra_discret_->NumMyRowNodes(); ++nd)
  {
    DRT::Node* opti_node = scatra_discret_->lRowNode(nd);
    int nodegid = opti_node->Id();
    for(unsigned int i=0; i<dirichlets.size(); ++i)
    {
      if (dirichlets[i]->ContainsNode(nodegid))
      {
        int dofgid = scatra_discret_->Dof(opti_node,0);
        int err = b->ReplaceGlobalValue(dofgid,0,0.0);
        if (err) dserror("could not replace global vector entry");
      }
    } // for(unsigned int i=0; i<dirichlets.size(); ++i)
  } // for(int nd=0; nd<scatra_discret_->NumMyRowNodes(); ++nd)


  // the rhsvec is ready for take off
  // solve the system now!
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(scatrasolverparams_,scatra_discret_->Comm(),DRT::Problem::Instance()->ErrorFile()->Handle()));
  solver->Solve(sysmatscatra->EpetraOperator(),adjoint_w_,b,true,true);

  std::string outname = name_;
  outname.append("_invadjoint_opti");
  scatra_output_->NewResultFile(outname,output_count_);
  output_count_++;
  scatraalgo_->DiscWriter().WriteMesh(0,0.0);
  scatraalgo_->DiscWriter().NewStep(1,1.0);
  scatraalgo_->DiscWriter().WriteElementData(true);
  scatraalgo_->DiscWriter().WriteVector("phinp",adjoint_w_);
  return;
} // void ACOU::InvAnalysis::SolveAdjointOptiProblem()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateGradient()
{
  // gradient consists of two parts: eintries with respect to the absorption coefficient
  // and entries with respect to the diffusion coefficient, both are treated in a similar way
  // and the formula to calculate the entries is given by

  /*
   * delta J = alpha mu_a delta mu_a + alpha D delta D + phi w delta mu_a - grad phi grad w delta D
   *
   * G_i_mu_a = alpha mu_a N_mu_a + phi w N_mu_a
   *
   * G_i_D = alpha D N_D - grad phi grad w N_D
   *
   * with N_mu_a and N_D the shape functions for the material parameters and hence 1 on every element
   * with this material and zero elsewhere (we need elements which belong to one material)
   *
   */

  G_.Scale(0.0);

  Teuchos::RCP<Epetra_Vector> node_mu_a_col = LINALG::CreateVector(*(scatra_discret_->NodeColMap()),true);
  LINALG::Export(*node_mu_a_,*node_mu_a_col);

  Teuchos::RCP<Epetra_Vector> rea_grad_ele = LINALG::CreateVector(*(scatra_discret_->ElementRowMap()),true);
  // now, we do not have to do all communication manually!
  for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = scatra_discret_->lRowElement(i);
    double rea_ele = static_cast<const MAT::ScatraMat*>(ele->Material().get())->ReaCoeff();
    const int* nodeids = ele->NodeIds();
    int numnode = ele->NumNode();
    double grad = rea_ele;
    for(int j=0; j<numnode; ++j)
    {
      double loc_rea_node = node_mu_a_col->operator [](scatra_discret_->NodeColMap()->LID(nodeids[j]));
      grad -= loc_rea_node/double(numnode);
    }
    rea_grad_ele->operator [](i) = (grad);
  }
  Teuchos::RCP<Epetra_Vector> rea_grad_ele_i = LINALG::CreateVector(*(scatra_discret_->ElementRowMap()),true);
  for(int mats=0; mats<nm_; mats++)
  {
    Epetra_SerialDenseVector tempparam(np_);
    tempparam.Scale(0.0);
    tempparam(mats) = 1.0;
    SetParameters(tempparam);
    Teuchos::RCP<Epetra_Vector> node_mu_a_col = LINALG::CreateVector(*(scatra_discret_->NodeColMap()),true);
    LINALG::Export(*node_mu_a_,*node_mu_a_col);
    for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
    {
      DRT::Element* ele = scatra_discret_->lRowElement(i);
      double rea_ele = static_cast<const MAT::ScatraMat*>(ele->Material().get())->ReaCoeff();
      const int* nodeids = ele->NodeIds();
      int numnode = ele->NumNode();
      double grad = rea_ele;
      for(int j=0; j<numnode; ++j)
      {
        double loc_rea_node = node_mu_a_col->operator [](scatra_discret_->NodeColMap()->LID(nodeids[j]));
        grad -= loc_rea_node/double(numnode);
      }
      rea_grad_ele_i->operator [](i) = (grad);
    }
    rea_grad_ele_i->Dot(*rea_grad_ele,&G_(mats));
  }
  G_.Scale(beta_);
  // this implementation depends heavily on the number of considered materials
  // if we got much less materials than elements we act differently than if we
  // got approximately as many materials as elements

  // in this case we loop all materials, set one material to one, and all others
  // to zero and evaluate the materials, or so to say, integrate the terms as described above
  Teuchos::ParameterList eleparams;

  INPAR::SCATRA::ScaTraType scatype = DRT::INPUT::IntegralValue<INPAR::SCATRA::ScaTraType>(scatraparams_,"SCATRATYPE");
  eleparams.set<int>("scatratype",scatype);
  scatra_discret_->SetState("phi",phi_);
  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou) )
    {
      opti_node = scatra_discret_->gNode(nd+maxnodeidacou);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;

    double loc_value = 0.0;
    if(acou_discret_->NodeRowMap()->LID(nd)>-1)
    {
      loc_value = adjoint_phi_0_->operator [](acou_discret_->NodeRowMap()->LID(nd));
    }
    double glo_value = 0.0;
    acou_discret_->Comm().SumAll(&loc_value,&glo_value,1);

    // ok, we got the value, we still need c, rho and mu_a, but they are stored on the nodemap of the scatra dis, so this should not be a problem
    if(scatra_discret_->Comm().MyPID() == optnodeowner)
    {
      int lnodeid = node_mu_a_->Map().LID(nd+maxnodeidacou);
      double c    = node_c_->operator [](lnodeid);
      int dofgid = scatra_discret_->Dof(opti_node,0);
      int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
      int err = psi->ReplaceMyValue(doflid,0,glo_value/c/c);
      if (err) dserror("could not replace local vector entry");
    }
  }
  scatra_discret_->SetState("psi",psi);

  Teuchos::RCP<Epetra_Vector> chi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  chi->Update(1.0,*adjoint_w_,0.0);
  scatra_discret_->SetState("w",chi);

  for(int mats=0; mats<nm_; mats++)
  {
    Epetra_SerialDenseVector tempparam(np_);
    tempparam.Scale(0.0);
    tempparam(mats) = 1.0;
    SetParameters(tempparam);

    eleparams.set<bool>("signum_mu",p_(mats)<0.0);
    eleparams.set<int>("action",SCATRA::calc_integr_grad_reac);
    eleparams.set<double>("regular_mu",alpha_*p_(mats));

    eleparams.set<double>("integr_grad_reac",0.0);

    scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    double loc_integr_grad_reac = eleparams.get<double>("integr_grad_reac");
    double glo_integr_grad_reac = 0.0;

    scatra_discret_->Comm().SumAll(&loc_integr_grad_reac,&glo_integr_grad_reac,1);

    G_(mats) += glo_integr_grad_reac;

  } // for(int mats=0; mats<nm_; mats++)
  scatra_discret_->ClearState();

  if(fdcheck_)
    FD_GradientCheck();

  return;
} // void ACOU::InvAnalysis::CalculateGradient()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::FD_GradientCheck()
{
  if(!myrank_)
  {
    std::cout<<"ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"<<std::endl;
    std::cout<<"ooooooooooooooooooooo FD Check oooooooooooooooooooooooooo"<<std::endl;
    std::cout<<"ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"<<std::endl;
  }

  Epetra_SerialDenseVector tempparams(np_);
  Epetra_SerialDenseVector FD_G(np_);

  double J_before = J_;
  Epetra_SerialDenseVector J_afterwards(np_);

  double delta = 1.0e-6;
  for(int i=0; i<np_; ++i) // loop parameters
  {
    tempparams = p_;
    // perturb ith parameter
    tempparams(i) += delta; // * p_(i) ;

    if(!myrank_)
    {
      std::cout<<std::endl;
      std::cout<<"ooooooooooooooooooooo FD Check iteration "<<i+1<<" of "<<np_<<std::endl;
      std::cout.precision(15);
      std::cout<<"ooooooooooooooooooooo parameter before ";
      for(int j=0; j<np_; ++j)
        std::cout<<p_(j)<<" ";
      std::cout<<" and disturbed ";
      for(int j=0; j<np_; ++j)
        std::cout<<tempparams(j)<<" ";
      std::cout<<std::endl;
    }

    SetParameters(tempparams);
    SolveStandardProblem();
    CalculateObjectiveFunctionValue();
    J_afterwards(i) = J_;
    FD_G(i) = (J_ - J_before) / delta;
  } // for(int i=0; i<np_; ++i) // loop parameters

  // output of the outcome of the check
  if(!myrank_)
  {
    std::cout<<"GRADIENT ACCORDING TO ADJOINT ANALYSIS"<<std::endl;
    for(int i=0; i<np_; ++i)
      std::cout<<G_(i)<<" ";
    std::cout<<endl;
    std::cout<<"GRADIENT ACCORDING TO FINITE DIFFERNECES"<<std::endl;
    for(int i=0; i<np_; ++i)
      std::cout<<FD_G(i)<<" ";
    std::cout<<endl;
    std::cout<<"VERHAELTNIS "<<FD_G(0)/G_(0)<<" DIFFERENZ "<<FD_G(0)-G_(0)<<" RELATIVER FEHLER "<<(FD_G(0)-G_(0))/FD_G(0)<<std::endl;
  }

  // we pretend this never happened
  SetParameters(p_);
  J_ = J_before;

  return;
} // void ACOU::InvAnalysis::FD_GradientCheck()


/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::UpdateHessian()
{
  if(iter_==0)
  {
//    if(p_.Norm2()<1.0)
//      H_.Scale(1.0);///G_.Norm2());
//    else
//      H_.Scale(p_.Norm2()/G_.Norm2()/10.0);
    return;
  }
  Epetra_SerialDenseVector s(np_);
  s = pm_;
  s.Scale(-1.0);
  s += p_;
  Epetra_SerialDenseVector y(np_);
  y = Gm_;
  y.Scale(-1.0);
  y += G_;

  Hm_ = H_;

  double rho = 1.0 / y.Dot(s);
  if (rho < 0.0 && !myrank_) // corresponds to y^T . s < 0
  {
    std::cout<<"WARNING: no guarantee for curvature condition to be satisfied, skipping Hessian Update"<<std::endl;
    return;
  }
  Epetra_SerialDenseMatrix temp(np_,np_);
  temp.Scale(0.0);
  for(int i=0; i<np_; ++i)
    temp(i,i) = 1.0;
  temp.Multiply('N','T',-rho,s,y,1.0);

  Epetra_SerialDenseMatrix temp2(np_,np_);
  temp2.Scale(0.0);
  for(int i=0; i<np_; ++i)
    temp2(i,i) = 1.0;
  temp2.Multiply('N','T',-rho,y,s,1.0);

  H_.Multiply('N','T',rho,s,s,0.0);

  Epetra_SerialDenseMatrix tempH(np_,np_);
  tempH.Multiply('N','N',1.0,Hm_,temp2,0.0);
  H_.Multiply('N','N',1.0,temp,tempH,1.0);

  return;
} // void ACOU::InvAnalysis::UpdateHessian()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::UpdateParameters()
{
  pm_ = p_;

//  if(iter_==0)
//    H_.Scale(1.0/G_.Norm2());//*p_.Norm2()/10.0);

  // BFGS
  Epetra_SerialDenseVector d(np_); // search direction
  d.Multiply('N','N',-1.0,H_,G_,0.0);

  // descent direction?
  if(d.Dot(G_)>=0)
  {
    H_.Print(std::cout);
    dserror("no descent direction");
  }
  if(!myrank_)
  {
    std::cout<<"search direction for line search"<<std::endl;
    d.Print(std::cout);
  }
  p_ = LineSearch(d);

  return;
} // void ACOU::InvAnalysis::UpdateParameters()


/*----------------------------------------------------------------------*/
Epetra_SerialDenseVector ACOU::InvAnalysis::LineSearch(Epetra_SerialDenseVector d)
{
  Epetra_SerialDenseVector pnew(d.M());

  // This is a backtracking line search, starting with initial alpha and reducing it,
  // as long as the sufficient decrease condition is not yet satisfied

  double alpha_0 = 1.0;// / d.Norm2();
  double c = 1.0e-4;
  double rho = 0.2;

  // does alpha_0 already satisfy the Armijo (sufficient decrease) condition?
  double alpha = alpha_0;
  double J_before = J_;
  double condition = J_before;// + c * alpha * d.Dot(G_);

  int count = 0;

  do{
    count++;

    if(!myrank_)
      std::cout<<"*************** line search iteration "<<count<<" of maximal "<<max_ls_iter_<<" line search iterations, condition "<<condition<<" step length "<<alpha<<std::endl;

    pnew = d;
    pnew.Scale(alpha);
    pnew += p_;

    SetParameters(pnew);
    SolveStandardProblem();
    CalculateObjectiveFunctionValue();
    alpha *= rho;

  } while(J_ >= condition && count<max_ls_iter_);

  if(count == max_ls_iter_ && J_ >= condition)
    dserror("line search did not succeed in %d iterations",max_ls_iter_);

  return pnew;
} // void ACOU::InvAnalysis::LineSearch()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::UpdateMaterial()
{
  SetParameters(p_);
  return;
} // void ACOU::InvAnalysis::UpdateMaterial()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateStatsAndService()
{
  // update iteration count
  iter_++;

  // update gradient
  Gm_ = G_;

  // calculate the norm of the difference of the given parameters, for the curious user
  Epetra_SerialDenseVector temp(np_);
  temp = pm_;
  temp.Scale(-1.0);
  temp += p_;
  normdiffp_ = temp.Norm2();

  return;
} // void ACOU::InvAnalysis::CalculateStatsAndService()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::OutputStats()
{
  double normB = 0.0;
  {
    Epetra_SerialDenseSolver inverseH;
    Epetra_SerialDenseMatrix B(np_,np_);
    for(int i=0; i<H_.M(); ++i)
      for(int j=0; j<H_.M(); ++j)
        B(i,j) = H_(i,j);
    inverseH.SetMatrix(B);
    inverseH.Invert();
    normB = B.NormInf();
  }
  double normH = H_.NormInf();

  if (!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"*** objective function value:             "<<J_<<std::endl;
    std::cout<<"*** error value:                          "<<error_<<std::endl;
    std::cout<<"*** norm of the gradient:                 "<<G_.Norm2()<<std::endl;
    std::cout<<"*** norm of the difference of parameters: "<<normdiffp_<<std::endl;
    std::cout<<"*** output count:                         "<<output_count_<<std::endl;
    std::cout<<"*** norm H "<<normH<<" norm B "<<normB<<" cond "<<normH*normB<<std::endl;
    std::cout<<"*** simulation time since start           "<<Teuchos::Time::wallTime()-tstart_<<std::endl;
    std::cout<<"*** parameters: ";
    for(int i=0; i<np_; ++i)
      std::cout<<p_(i)<<" ";
    std::cout<<std::endl;
  }
  return;
} // void ACOU::InvAnalysis::OutputStats()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::ConvergenceCheck()
{
  // check convergence criteria
  bool unconverged;

  if(calcacougrad_)
    unconverged = error_>tol_ && iter_<max_iter_ && G_.Norm2()>tol_grad_;
  else
    unconverged = error_>tol_ && iter_<max_iter_;

  if(std::abs(normdiffp_)<1.0e-14)
    unconverged = false;

  return unconverged;
} // bool ACOU::InvAnalysis::ConvergenceCheck()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ReadInParameters()
{
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance(0)->Materials()->Map();
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const Teuchos::RCP<MAT::PAR::Material> material = curr->second;
    switch(material->Type())
    {
    case INPAR::MAT::m_scatra:
    {
      MAT::PAR::ScatraMat* params = dynamic_cast<MAT::PAR::ScatraMat*>(material->Parameter());
      if (!params) dserror("Cannot cast material parameters");
      const int j = p_.Length();
      p_.Resize(j+1);
      p_[j]   = params->reacoeff_;
      nm_++;
      break;
    }
    case INPAR::MAT::m_acousticmat:
    {
      break;
    }
    default:
      dserror("photoacoustic inverse analysis only implemented for scatra mat and acou mat");
      break;
    }
  }
  return;
} // void ACOU::InvAnalysis::ReadInParameters()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::SetParameters(Epetra_SerialDenseVector p_cur)
{
  acou_discret_->Comm().Broadcast(&p_cur[0],np_,0);

  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance(0)->Materials()->Map();
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator curr;
  int j=0;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const Teuchos::RCP<MAT::PAR::Material> material = curr->second;
    switch(material->Type())
    {
    case INPAR::MAT::m_scatra:
    {
      MAT::PAR::ScatraMat* params = dynamic_cast<MAT::PAR::ScatraMat*>(material->Parameter());
      if (!params) dserror("Cannot cast material parameters");
      params->SetReaCoeff(p_cur(j));
      j += 1;
      break;
    }
    case INPAR::MAT::m_acousticmat:
    {
      break;
    }
    default:
      dserror("photoacoustic inverse analysis only implemented for scatra mat and acou mat");
      break;
    }
  }

  // we have to recalculate node_mu_a_ every time we update parameters
  // now: fill the vectors
  int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd) // cannot loop scatra nodes, because they don't necessarily start with gid 0
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou) )
    {
      opti_node = scatra_discret_->gNode(nd+maxnodeidacou);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;

    int nodelid = scatra_discret_->NodeRowMap()->LID(nd+maxnodeidacou);

    int loc_numoptiele = 0;
    double loc_mu_a = 0.0;
    for(int roel = 0; roel<scatra_discret_->NumMyRowElements(); ++roel)
    {
      DRT::Element* roptele = scatra_discret_->lRowElement(roel);
      const int* nodeids = roptele->NodeIds();
      int numnode = roptele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd+maxnodeidacou )
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(roptele->Material().get());
          loc_mu_a += actmat->ReaCoeff();
          loc_numoptiele++;
        }
    }
    int glo_numoptiele = 0;
    double glo_mu_a = 0.0;
    scatra_discret_->Comm().SumAll(&loc_numoptiele,&glo_numoptiele,1);
    scatra_discret_->Comm().SumAll(&loc_mu_a,&glo_mu_a,1);
    glo_mu_a /= double(glo_numoptiele);

    // write this value to vecot
    if( nodelid >= 0 ) // only on owning proc
      node_mu_a_->ReplaceMyValue(nodelid,0,glo_mu_a);
  }

  return;
} // void ACOU::InvAnalysis::SetParameters(Epetra_SerialDenseVector p_cur)

const Teuchos::RCP<Epetra_Vector> ACOU::InvAnalysis::ElementMatVec()
{
  Teuchos::RCP<Epetra_Vector> outvec = LINALG::CreateVector(*(scatra_discret_->ElementRowMap()),true);

  for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
  {
    DRT::Element* optele = scatra_discret_->lRowElement(i);
    const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(optele->Material().get());
    double mu_a = actmat->ReaCoeff();
    outvec->ReplaceMyValue(i,0,mu_a);
  }

  return outvec;
} // const Teuchos::RCP<Epetra_Vector> ACOU::InvAnalysis::ElementMatVec()

Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()
{
  return Teuchos::rcp(new AcouInvResultTest(*this));
} // Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()
