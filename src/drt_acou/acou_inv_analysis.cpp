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
#include "acou_ele.H"
#include "acou_impl_euler.H"
#include "acou_impl_trap.H"
#include "acou_impl_dirk.H"
#include "acou_impl_bdf.H"
#include "acou_inv_resulttest.H"
#include "acou_ele_action.H"
#include "pat_matpar_manager.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/material.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inv_analysis/invana_utils.H"
#include "../drt_inv_analysis/regularization_base.H"
#include "../drt_inv_analysis/regularization_tikhonov.H"
#include "../drt_inv_analysis/regularization_totalvariation.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
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
    opti_(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")),
    myrank_(acoudis->Comm().MyPID()),
    error_(1.0e6),
    tol_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL")),
    iter_(0),
    max_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_MAX_RUN")),
    max_ls_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_LS_MAX_RUN")),
    output_count_(0),
    fdcheck_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("FDCHECK"))),
    ls_c_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("LS_DECREASECOND")),
    ls_rho_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("LS_STEPLENGTHRED")),
    ls_gd_scal_(1.0),
    calcacougrad_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("INV_TOL_GRAD_YN"))),
    scalegradele_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("ELE_SCALING"))),
    normgrad_(1.0e6),
    normgrad0_(1.0e6),
    dtacou_(acouparams_->get<double>("TIMESTEP")),
    J_(0.0),
    normdiffp_(0.0),
    tstart_(Teuchos::Time::wallTime())
{
  // some checks, if everything is alright
  if(ls_rho_>=1.0) dserror("LS_STEPLENGTHRED has to be smaller than 1.0");
  if(ls_c_>1.0) dserror("LS_DECREASECOND is usually chosen in between 0.0 and 0.01, decrease it!");
  if(not acou_discret_->Filled() || not acou_discret_->HaveDofs()) dserror("acoustical discretization is not complete or has no dofs");
  if(not scatra_discret_->Filled() || not scatra_discret_->HaveDofs()) dserror("scatra discretization is not complete or has no dofs");

  // set up of the output
  scatra_output_ = scatra_discret_->Writer();
  name_ = DRT::Problem::Instance()->OutputControlFile()->FileName();

  // vectors we need
  adjoint_phi_0_ = LINALG::CreateVector(*(acou_discret_->NodeRowMap()),true);
  phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  adjoint_w_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  // tolerance for the gradient, if used
  if(calcacougrad_) tol_grad_ = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL_GRAD");
  else tol_grad_ = 0.0;

  // read monitor file, create multivector and map for measured values
  ReadMonitor(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE"));

  // allocate vectors
  node_mu_a_ = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);

  // and fill them
  ComputeNodeBasedVectors();

  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
    {
      //matman_ = Teuchos::rcp(new STR::INVANA::MatParManagerPerElement(scatra_discret_));
      matman_ = Teuchos::rcp(new ACOU::PatMatParManagerPerElement(scatra_discret_,scalegradele_));
    }
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
    {
      //matman_ = Teuchos::rcp(new STR::INVANA::MatParManagerUniform(scatra_discret_));
      matman_ = Teuchos::rcp(new ACOU::PatMatParManagerUniform(scatra_discret_));
    }
      break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
    break;
  }
  matman_->Setup();

  // create regularization manager
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
  case INPAR::INVANA::stat_inv_reg_none:
    break;
  case INPAR::INVANA::stat_inv_reg_tikhonov:
  {
    regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov(invp));
    break;
  }
  case INPAR::INVANA::stat_inv_reg_totalvariation:
  {
    regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation(invp));
    break;
  }
  default:
    dserror("no valid regularization type provided");
    break;
  }
  if (regman_!=Teuchos::null)
  {
    regman_->Init(scatra_discret_,matman_->GetConnectivityData());
    regman_->Setup();
  }

  ssize_ =  invp.get<int>("SIZESTORAGE");
  ssize_ *= matman_->NumVectors();
  actsize_ = 0;

  sstore_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, matman_->ParamLayoutMap().get(), true));
  ystore_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, matman_->ParamLayoutMap().get(), true));

  objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(),true));
  objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(),true));
  step_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(), true));
  d_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(), true));

} // ACOU::InvAnalysis::InvAnalysis(...)

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ReadMonitor(std::string monitorfilename)
{
  // initialize acou_rhs_: we need a vector with the nodes of the boundary where
  // the pressure is monitored-> Read the monitor file and create a vector with
  // corresponding nodes OR take the boundary where absorbing bcs are prescribed!
  // we need this map extractor thing!
  // we deal with NODES here, not with DOFS

  std::string condname = "PressureMonitor";
  std::vector<DRT::Condition*> pressuremon;
  acou_discret_->GetCondition(condname,pressuremon);
  if(pressuremon.size()==0)
    dserror("you have to use pressure monitor conditions for inverse analysis!");
  const std::vector<int> pressuremonnodes = *(pressuremon[0]->Nodes());
  std::vector<int> pressuremonnodesunique;

  // create unique map
  acou_discret_->Comm().Barrier();
  for(int i=0; i<acou_discret_->Comm().NumProc(); ++i)
  {
    if(acou_discret_->Comm().MyPID() == i)
    {
      for(unsigned int j=0; j<pressuremonnodes.size(); ++j)
      {
        if(acou_discret_->HaveGlobalNode(pressuremonnodes[j]))
        {
          if(acou_discret_->gNode(pressuremonnodes[j])->Owner()==int(i))
            pressuremonnodesunique.push_back(pressuremonnodes[j]);
        }
      }
    }
    acou_discret_->Comm().Barrier();
  }
  acou_discret_->Comm().Barrier();


//abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, abcnodes.size(), &abcnodes[0], 0, acou_discret_->Comm()));
  abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, pressuremonnodesunique.size(), &pressuremonnodesunique[0], 0, acou_discret_->Comm()));
  abcnodes_mapex_ = Teuchos::rcp(new LINALG::MapExtractor(*(acou_discret_->NodeRowMap()),abcnodes_map_,true));

  // determine the number of vectors for monitoring
  // this is naive: later on, we won't be able to store everything at once and we have to implement
  // a smarter approach to reduce storage requirements
  int numvec = acouparams_->get<int>("NUMSTEP");
  int oderso = acouparams_->get<double>("MAXTIME")/dtacou_;
  if ( oderso < numvec)
    numvec = oderso;

  acou_rhs_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec+1,true));
  acou_rhsm_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec+1,true));

  unsigned int nsteps = 0;
  Epetra_SerialDenseVector mcurve;
  std::vector<double>  timesteps;

  // get measured values
  // open monitor file and read it
  unsigned int nnodes = 0;
  {
    char* foundit = NULL;

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
    nsteps = strtol(foundit,&foundit,10);
    timesteps.resize(nsteps);
    // read nnodes
    foundit = strstr(buffer,"nnodes"); foundit += strlen("nnodes");
    nnodes = strtol(foundit,&foundit,10);
    // read nodes
    nodes_.resize(nnodes);
    for (unsigned int i=0; i<nnodes; ++i)
    {
      fgets(buffer,150000,file);
      foundit = buffer;
      nodes_[i] = strtol(foundit,&foundit,10) - 1;

      //if (!myrank_) printf("Monitored node %d ",nodes_[i]);
      //if (!myrank_) printf("\n");
    }
    // read in measured curve
    {
      mcurve = Epetra_SerialDenseVector(nnodes*nsteps);

      //if (!myrank_) printf("nsteps %d nnode %d\n",nsteps_,nnodes);

      // read comment lines
      foundit = buffer;
      fgets(buffer,150000,file);
      while(strstr(buffer,"#"))
        fgets(buffer,150000,file);

      // read in the values for each node
      unsigned int count = 0;
      for (unsigned int i=0; i<nsteps; ++i)
      {
        // read the time step
        timesteps[i] = strtod(foundit,&foundit);
        for (unsigned int j=0; j<nnodes; ++j)
          mcurve[count++] = strtod(foundit,&foundit);
        fgets(buffer,150000,file);
        foundit = buffer;
      }
      if (count != nnodes*nsteps) dserror("Number of measured pressure values wrong on input");
    }
  }

  double eps = dtacou_/1000.0;
  if((numvec-1)*dtacou_>timesteps[nsteps-1]+eps) dserror("You want to simulate till %.15f but your monitor file only provides values till %.15f! Fix it!",(numvec-1)*dtacou_,timesteps[nsteps-1]);
  if (nnodes != pressuremonnodes.size()) dserror("Given number of nodes in boundary condition and monitor file don't match!");
  if (nodes_ != pressuremonnodes) dserror("And please provide the correct order (feel free to reimplement)");

//  for(int i=0; i<acou_discret_->Comm().NumProc(); ++i)
//  {
//    acou_discret_->Comm().Barrier();
//    if(i == acou_discret_->Comm().MyPID())
//    {
//      std::cout<<"mypid "<<acou_discret_->Comm().MyPID()<<std::endl;
//      mcurve.Print(std::cout);std::cout<<std::endl;
//    }
//    acou_discret_->Comm().Barrier();
//  }

  // every proc knows mcurve, now, we want to write mcurve to a Epetra_MultiVector in the same form as acou_rhs_
  // with the same parallel distribution!
  // and we want to interpolate measured values in case the monitored time step size is not the same as the one for the simulation
  acou_rhsm_->PutScalar(0.0);
  if( timesteps[0] != 0.0 )
    dserror("your measured values have to start at time 0.0");
  else if( timesteps[0] == 0.0 && timesteps[1] == dtacou_ ) // the standard case
    for(unsigned int i=0; i<nnodes; ++i)
    {
      if( acou_discret_->HaveGlobalNode(nodes_[i]) )
        for(unsigned int j=0; j<nsteps; ++j)
          acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,mcurve(i+j*nnodes)); // the proc who has this row, writes the value
    }
  else // we have to interpolate!
  {
    if( numvec < int(nsteps) )
    {
      dserror("set your time step size smaller, at least to %14f or implement here",timesteps[1]-timesteps[0]);
    }
    else
    {

      for(unsigned int i=0; i<nnodes; ++i)
      {
        if( acou_discret_->HaveGlobalNode(nodes_[i]) )
        {
          for(int j=0; j<numvec; ++j)
          {
            double actualt = j * dtacou_; // we need values for this time
            int timeval = 0;

            // find next higher and next lower value
            while(actualt>timesteps[timeval]-eps)
            {
              timeval++;
            }

            // timesteps_[timeval] has the next higher point in time
            // now interpolate from this and the value before
            if(timeval == 0)
            {
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,0.0);
            }
            else if(actualt<timesteps[timeval]+eps && actualt>timesteps[timeval]-eps) // then this is more or less it
            {
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,mcurve(i+timeval*nnodes));
            }
            else
            {
              double value = mcurve(i+(timeval-1)*nnodes) + (mcurve(i+(timeval)*nnodes)-mcurve(i+(timeval-1)*nnodes)) * (actualt - timesteps[timeval-1]) / (timesteps[timeval]-timesteps[timeval-1]);
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,value);
            }
          } // for(int j=0; j<numvec; ++j)
        } // if( acou_discret_->HaveGlobalNode(nodes_[i]) )
        acou_discret_->Comm().Barrier();
      } // for(unsigned int i=0; i<nnodes; ++i)
    } // else ** if( numvec < nsteps_ )
  } // else ** if( timesteps_[0] == dtacou_ || (timesteps_[0]==0.0 && timesteps_[1] = dtacou_) )

  return;
} // void ACOU::InvAnalysis::ReadMonitor(...)

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ComputeNodeBasedVectors()
{
  // fill the vectors
  int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd) // cannot loop scatra nodes, because they don't necessarily start with gid 0
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou+1) )
    {
      opti_node = scatra_discret_->gNode(nd+maxnodeidacou+1);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;
    // here, every proc knows the owner and the gid of the optical node

    // now, every proc knows the speed of sound and density of this node -> write them to vector
    int nodelid = scatra_discret_->NodeRowMap()->LID(nd+maxnodeidacou+1);

    // we have to do the same procedure for the absorption coefficient with the scatra discretization
    int loc_numoptiele = 0;
    double loc_mu_a = 0.0;
    for(int roel = 0; roel<scatra_discret_->NumMyRowElements(); ++roel)
    {
      DRT::Element* roptele = scatra_discret_->lRowElement(roel);
      const int* nodeids = roptele->NodeIds();
      int numnode = roptele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd+maxnodeidacou+1 )
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(roptele->Material().get());
          loc_mu_a += actmat->ReaCoeff(scatra_discret_->ElementColMap()->LID(roptele->Id()));
          loc_numoptiele++;
        }
    }
    int glo_numoptiele = 0;
    double glo_mu_a = 0.0;
    scatra_discret_->Comm().SumAll(&loc_numoptiele,&glo_numoptiele,1);
    scatra_discret_->Comm().SumAll(&loc_mu_a,&glo_mu_a,1);
    glo_mu_a /= double(glo_numoptiele);

    // write this value to vect
    if( nodelid >= 0 ) // only on owning proc
      node_mu_a_->ReplaceMyValue(nodelid,0,glo_mu_a);
  } // for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::Integrate()
{
  // Solve the standard problem
  SolveStandardProblem();

  bool success = true;

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

    // Update the sought parameters
    success = UpdateParameters();

    // Calculate some useful numbers
    CalculateStatsAndService();

    // Output some useful user information, like time consume, solution advance, ...
    OutputStats();

  } while ( ConvergenceCheck() && success );

  if( success == false && !myrank_ )
  {
    std::cout<<std::endl;
    std::cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"<<std::endl;
    std::cout<<"optimization terminated unsuccessfully, line search failed"<<std::endl;
    std::cout<<"-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"<<std::endl;
    std::cout<<std::endl;
  }
  else if( error_<=tol_ && !myrank_ )
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
  else if(normgrad_<= tol_grad_ && !myrank_)
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
  acoualgo_->SetInitialPhotoAcousticField(phi_,scatra_discret_,meshconform);

  // we have to call a slightly changed routine, which fills our history vector which we need for the adjoint problem
  acou_rhs_->Scale(0.0);

  double norminf = 0.0;
  matman_->GetMatParams()->NormInf(&norminf);
  if(norminf>1.0e-5) // we do not have to solve when the maximum reaction coefficient is zero, just a waste of time!
    acoualgo_->Integrate(acou_rhs_,abcnodes_mapex_);

//    if(iter_==0 && output_count_<3) // do you want this output?
//    {
//      std::cout<<"measured values"<<std::endl;
//      std::cout.precision(15);
//      for(int j=0; j<acou_rhs_->NumVectors(); ++j)
//      {
//        std::cout<<dtacou_*double(j)<<" ";
//        for(int i=0; i<acou_rhs_->Map().NumGlobalElements(); ++i)
//          std::cout<<" "<<acou_rhs_->operator ()(j)->operator [](i);
//        std::cout<<std::endl;
//      }
//    }

  return;
} // void ACOU::InvAnalysis::SolveStandardProblem()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateObjectiveFunctionValue()
{
  J_ = 0.0;

  // contribution from difference between measured and simulated values
  Epetra_MultiVector tempvec = Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors());
  tempvec.Update(1.0,*acou_rhsm_,0.0);
  tempvec.Update(1.0,*acou_rhs_,-1.0);

  tempvec.Multiply(1.0,tempvec,tempvec,0.0);
  acou_discret_->Comm().Barrier();

  Epetra_SerialDenseVector normvec(acou_rhsm_->NumVectors());
  tempvec.Norm1(normvec.Values());

  J_ += error_ = 0.5 * normvec.Norm1();

  // contribution from regularization
  if(regman_ != Teuchos::null)
    regman_->Evaluate(*(matman_->GetParams()),&J_);

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

  // acou_rhs_ has to be scaled with weighting (adjoint of the weighting)
  Teuchos::RCP<Epetra_Vector> touchcountvec = LINALG::CreateVector(*abcnodes_map_);
  acoualgo_->FillTouchCountVec(touchcountvec);
  acou_rhs_->Multiply(1.0,*touchcountvec,*acou_rhs_,0.0);

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
  acoualgo_->NodalPsiField(adjoint_phi_0_);

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
  Teuchos::RCP<Epetra_Vector> rhsvec;
  bool meshconform = DRT::INPUT::IntegralValue<bool>(*acouparams_,"MESHCONFORM");
  if(meshconform)
  {
    rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

    int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
    for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
    {
      // get node and owner
      int myoptnodeowner = -1;
      int optnodeowner = -1;
      DRT::Node* opti_node = NULL;
      if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou+1) )
      {
        opti_node = scatra_discret_->gNode(nd+maxnodeidacou+1);
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
        int lnodeid = node_mu_a_->Map().LID(nd+maxnodeidacou+1);
        double mu_a = node_mu_a_->operator [](lnodeid);
        int dofgid = scatra_discret_->Dof(opti_node,0);
        int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
        int err = rhsvec->ReplaceMyValue(doflid,0,-glo_value*mu_a);
        if (err) dserror("could not replace local vector entry");
      }
    }
  } // if(meshconform)
  else
  {
    rhsvec = CalculateNonconformRhsvec(adjoint_phi_0_);
    for(int i=0; i<node_mu_a_->MyLength(); ++i)
    {
      double mu_a = node_mu_a_->operator [](i);
      int dofgid = scatra_discret_->Dof(scatra_discret_->lRowNode(i),0);
      int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
      rhsvec->operator [](doflid) *= -mu_a;
    }
  } // else ** if(meshconform)

  Teuchos::ParameterList eleparams;
  scatra_discret_->SetState("rhsnodebasedvals",rhsvec);
  eleparams.set<int>("action",SCATRA::calc_integr_pat_rhsvec);
  Teuchos::RCP<Epetra_Vector> b = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,b,Teuchos::null,Teuchos::null);


  std::string condname = "Dirichlet";
  std::vector<DRT::Condition*> dirichlets;
  scatra_discret_->GetCondition(condname,dirichlets);
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
Teuchos::RCP<Epetra_Vector> ACOU::InvAnalysis::CalculateNonconformRhsvec(Teuchos::RCP<Epetra_Vector> acounodevec)
{
  // this function is similar to the mapping in void ACOU::AcouImplicitTimeInt::SetInitialPhotoAcousticField
  // just the other way round

  Teuchos::RCP<Epetra_Vector> rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  // export input vector to column map
  Teuchos::RCP<Epetra_Vector> acounodeveccol = Teuchos::rcp(new Epetra_Vector(*(acou_discret_->NodeColMap())),true);
  LINALG::Export(*acounodevec,*acounodeveccol);

  int numdim = DRT::Problem::Instance()->NDim();
  int minoptnodegid = scatra_discret_->NodeRowMap()->MinAllGID();
  for(int optnd=0; optnd<scatra_discret_->NumGlobalNodes(); ++optnd)
  {
    DRT::Node* optnode = NULL;
    int myoptnodeowner = -1;
    if(scatra_discret_->HaveGlobalNode(optnd+minoptnodegid))
    {
      optnode = scatra_discret_->gNode(optnd+minoptnodegid);
      myoptnodeowner = optnode->Owner();
      if(myoptnodeowner != myrank_) myoptnodeowner = -1;
    }
    int optnodeowner = -1;
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);

    double optnodecoords[numdim];
    if(myrank_==optnodeowner)
    {
      for(int d=0; d<numdim; ++d)
        optnodecoords[d] = optnode->X()[d];
    }
    scatra_discret_->Comm().Broadcast(&optnodecoords[0],numdim,optnodeowner);

    double r = 0.0;
    for(int acouel=0; acouel<acou_discret_->NumMyRowElements(); ++acouel)
    {
      DRT::Element* ele = acou_discret_->lRowElement(acouel);
      // get the nodes of this element, and then check if acoustical node is inside
      if(ele->Shape()==DRT::Element::quad4)
      {
        double acounodecoords[4][numdim];
        double minmaxvals[2][numdim];
        for(int j=0; j<numdim; ++j)
        {
          minmaxvals[0][j] = 1.0e6; // minvals
          minmaxvals[1][j] = -1.0e6; // maxvals
        }
        for(int nd=0;nd<4;++nd) // quad4 has 4 nodes
          for(int d=0;d<numdim;++d)
          {
            acounodecoords[nd][d] = ele->Nodes()[nd]->X()[d];
            if(acounodecoords[nd][d] < minmaxvals[0][d]) minmaxvals[0][d]=acounodecoords[nd][d];
            if(acounodecoords[nd][d] > minmaxvals[1][d]) minmaxvals[1][d]=acounodecoords[nd][d];
          }
        // check, if acoustical node is in bounding box
        bool inside = true;
        for(int d=0;d<numdim;++d)
          if(optnodecoords[d]>minmaxvals[1][d]+5.0e-5 || optnodecoords[d]<minmaxvals[0][d]-5.0e-5)
            inside=false;
        if(inside)
        {
          // solve for xi by local Newton
          LINALG::Matrix<2,1> F(true);
          LINALG::Matrix<2,2> dFdxi(true);
          LINALG::Matrix<2,1> xi(true);
          LINALG::Matrix<2,1> deltaxi(true);
          double deltaxinorm = 0.0;
          int count = 0;
          do{
            count++;
            F(0) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * acounodecoords[0][0]
                 + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * acounodecoords[1][0]
                 + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * acounodecoords[2][0]
                 + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * acounodecoords[3][0]  - optnodecoords[0];
            F(1) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * acounodecoords[0][1]
                 + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * acounodecoords[1][1]
                 + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * acounodecoords[2][1]
                 + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * acounodecoords[3][1]  - optnodecoords[1] ;

            dFdxi(0,0) = - 0.25 * (1. - xi(1)) * acounodecoords[0][0]
                         + 0.25 * (1. - xi(1)) * acounodecoords[1][0]
                         + 0.25 * (1. + xi(1)) * acounodecoords[2][0]
                         - 0.25 * (1. + xi(1)) * acounodecoords[3][0] ;
            dFdxi(0,1) = - 0.25 * (1. - xi(0)) * acounodecoords[0][0]
                         - 0.25 * (1. + xi(0)) * acounodecoords[1][0]
                         + 0.25 * (1. + xi(0)) * acounodecoords[2][0]
                         + 0.25 * (1. - xi(0)) * acounodecoords[3][0] ;
            dFdxi(1,0) = - 0.25 * (1. - xi(1)) * acounodecoords[0][1]
                         + 0.25 * (1. - xi(1)) * acounodecoords[1][1]
                         + 0.25 * (1. + xi(1)) * acounodecoords[2][1]
                         - 0.25 * (1. + xi(1)) * acounodecoords[3][1] ;
            dFdxi(1,1) = - 0.25 * (1. - xi(1)) * acounodecoords[0][1]
                         - 0.25 * (1. + xi(1)) * acounodecoords[1][1]
                         + 0.25 * (1. + xi(1)) * acounodecoords[2][1]
                         + 0.25 * (1. - xi(1)) * acounodecoords[3][1] ;

            LINALG::FixedSizeSerialDenseSolver<2,2,1> inverser;
            inverser.SetMatrix(dFdxi);
            inverser.SetVectors(deltaxi,F);
            inverser.Solve();

            deltaxinorm = deltaxi.Norm2();
            xi.Update(-1.0,deltaxi,1.0);
          } while ( deltaxinorm > 1.0e-8 && count < 10 );
          if(!(count == 10 || xi.NormInf()>1.0+0.15))
          {
            // get the values!
            double values[4] = {0};
            for(int nd=0;nd<4;++nd)
            {
              int lid = acou_discret_->NodeColMap()->LID(ele->Nodes()[nd]->Id());
              if(lid<0)
                dserror("node of element not on this processor");
              else
                values[nd] = acounodeveccol->operator [](lid);
            }
            r = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * values[0]
              + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * values[1]
              + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * values[2]
              + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * values[3];
          }
        } // if(inside)
      }
      else dserror("up to now only implemented for quad4");
    } // for(int acouel=0; acouel<acou_discret_->NumMyRowElements(); ++acouel)
    // one processor might provide a value

    double glob_p_min = 0.0;
    scatra_discret_->Comm().MinAll(&r,&glob_p_min,1);
    double glob_p_max = 0.0;
    scatra_discret_->Comm().MaxAll(&r,&glob_p_max,1);
    // take higher absolut values
    double glob_p = 0.0;
    if(std::abs(glob_p_min)>std::abs(glob_p_max))
      glob_p = glob_p_min;
    else
      glob_p = glob_p_max;

    // set p value in node based vector
    if(myrank_==optnodeowner && glob_p != 0.0)
    {
      int dof = scatra_discret_->Dof(optnode,0);
      int lid = scatra_discret_->DofRowMap()->LID(dof);
      if(lid<0) dserror("cannot find dof for node %d ",optnd);

      int err = rhsvec->ReplaceMyValue(lid,0,glob_p);
      if (err) dserror("could not replace local vector entry");
    }
  } // for(int optnd=0; optnd<scatra_discret_->NumGlobalNodes(); ++optnd)


  return rhsvec;
} // void ACOU::InvAnalysis::CalculateNonconformRhsvec()


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

  //zero out gradient vector initially
  objgrad_->Scale(0.0);

  // contribution from regularization
  if(regman_ != Teuchos::null)
    regman_->EvaluateGradient(*(matman_->GetParams()),objgrad_);

  scatra_discret_->SetState("dual phi",adjoint_w_);
  scatra_discret_->SetState("phi",phi_);

  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  bool meshconform = DRT::INPUT::IntegralValue<bool>(*acouparams_,"MESHCONFORM");
  if(meshconform)
  {
    int maxnodeidacou = acou_discret_->NodeRowMap()->MaxAllGID();
    for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
    {
      // get node and owner
      int myoptnodeowner = -1;
      int optnodeowner = -1;
      DRT::Node* opti_node = NULL;
      if( scatra_discret_->HaveGlobalNode(nd+maxnodeidacou+1) )
      {
        opti_node = scatra_discret_->gNode(nd+maxnodeidacou+1);
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
        int dofgid = scatra_discret_->Dof(opti_node,0);
        int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
        int err = psi->ReplaceMyValue(doflid,0,glo_value);
        if (err) dserror("could not replace local vector entry");
      }
    }
  }
  else
    psi = CalculateNonconformRhsvec(adjoint_phi_0_);

  scatra_discret_->SetState("psi",psi);
  matman_->Evaluate(0.0,objgrad_,true);

  INVANA::MVNorm(*objgrad_,*matman_->ParamLayoutMapUnique(),2,&normgrad_);
  if( iter_ == 0)
    normgrad0_ = normgrad_;

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


  if(!myrank_)
    std::cout<<"Gradient according to adjoint analysis"<<std::endl;
  objgrad_->Print(std::cout);

  double J_before = J_;
  Epetra_MultiVector tempvec = Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors());
  tempvec.Update(1.0,*acou_rhs_,0.0);

  int numparams = matman_->NumVectors();

  Teuchos::RCP<Epetra_MultiVector> objgradFD = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(),true));

  double perturba = 1.0e-3;
  double perturbb = 1.0e-4;

  Teuchos::RCP<Epetra_MultiVector> perturb = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(),true));
  perturb->Update(1.0,*(matman_->GetParams()),0.0);

  //keep a copy of the current parameters to reset after perturbation:
  Teuchos::RCP<Epetra_MultiVector> pcurr = Teuchos::rcp(new Epetra_MultiVector(*(matman_->ParamLayoutMap()), matman_->NumVectors(),true));
  pcurr->Update(1.0,*(matman_->GetParams()),0.0);

  double pn=0.0;
  double p=0.0;
  double dp=0.0;
  for (int i=0; i<numparams; i++)
  {
    for (int j=0; j<matman_->ParamLayoutMap()->NumGlobalElements(); j++)
    {
      int lid = matman_->GetParams()->Map().LID(j+matman_->ParamLayoutMap()->MinAllGID());
      if(lid > -1)   p = (*(*(matman_->GetParams()))(i))[lid];
      else           p = 0.0;

      pn = p+p*perturba+perturbb;
      std::cout<<"p "<<p<<" disturbed "<<pn<<std::endl;
      perturb->ReplaceGlobalValue(j+matman_->ParamLayoutMap()->MinAllGID(),i,pn);
      perturb->Print(std::cout);
      matman_->ReplaceParams(*perturb);

      SolveStandardProblem();
      perturb->Update(1.0,*pcurr,0.0);

      CalculateObjectiveFunctionValue();
      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      objgradFD->ReplaceGlobalValue(j+matman_->ParamLayoutMap()->MinAllGID(),i,dp);
    }
  }

  // we pretend this never happened
  matman_->ReplaceParams(*pcurr);
  J_ = J_before;

  if(!myrank_)
    std::cout<<"Gradient according to finite difference check"<<std::endl;
  objgradFD->Print(std::cout);

  acou_rhs_->Update(1.0,tempvec,0.0);

  return;
} // void ACOU::InvAnalysis::FD_GradientCheck()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::UpdateParameters()
{
  double J_before = J_;

  if(iter_==0 || opti_ == INPAR::ACOU::inv_gd)
    d_->Update(-1.0,*objgrad_,0.0);
  else
  {
    StoreVectors();
    ComputeDirection();
  }

  if(iter_==0 || opti_ == INPAR::ACOU::inv_gd) // in the first iteration the lbfgs does not know how long the step should be, gradient descent should always do scaled step length
    d_->Scale(1.0/normgrad0_*double(objgrad_->Map().NumGlobalElements()));

  std::cout<<"search direction "<<std::endl;d_->Print(std::cout);
  bool success = LineSearch();

  if (success == false && opti_ == INPAR::ACOU::inv_lbfgs && iter_ != 0)
  {
    // in this case we try a steepest descent step!
    d_->Update(-1.0/normgrad0_,*objgrad_,0.0);

    if(!myrank_)
      std::cout<<"line search failed, we try steepest descent from now on!"<<std::endl;

    // reset objective function
    J_ = J_before;

    // and in every following step we do steepest descent
    opti_ = INPAR::ACOU::inv_gd;
    success = LineSearch();
  }

  return success;
} // void ACOU::InvAnalysis::UpdateParameters()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::StoreVectors()
{
  if (iter_*matman_->NumVectors()<=ssize_) // we have "<=" since we do not store the first run
    actsize_+=matman_->NumVectors();

  Epetra_MultiVector s(*(matman_->ParamLayoutMap()), (matman_->NumVectors()),true);

  //push back s
  s.Update(1.0,*(matman_->GetParams()),-1.0,*(matman_->GetParamsOld()),0.0);
  for (int i=0; i<s.NumVectors(); i++)
    sstore_->UpdateSteps(*s(i));

  // push back y
  s.Update(1.0,*objgrad_,-1.0,*objgrad_o_,0.0);
  for (int i=0; i<s.NumVectors(); i++)
    ystore_->UpdateSteps(*s(i));

  return;
} // void ACOU::InvAnalysis::StoreVectors()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ComputeDirection()
{

  d_->Update(1.0,*objgrad_,0.0);
  std::vector<double> alpha;

  // loop steps
  for (int i=0; i>-actsize_; i-=matman_->NumVectors())
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;

    int ind=0;
    for (int j=matman_->NumVectors(); j>0; j--)
    {
      INVANA::MVDotProduct(*(*ystore_)(i-j+1),*(*sstore_)(i-j+1),*matman_->ParamLayoutMapUnique(),&a);
      INVANA::MVDotProduct(*(*sstore_)(i-j+1),*(*d_)(ind),*matman_->ParamLayoutMapUnique(),&b);
      ind++;
      aa += a;
      bb += b;
    }
    alpha.push_back(1/aa*bb);

    ind=0;
    for (int j=matman_->NumVectors(); j>0; j--)
    {
      (*d_)(ind)->Update(-1.0*alpha.back(), *(*ystore_)(i-j+1),1.0 );
      ind++;
    }
  }

  // Some scaling of the initial hessian might come in here but has not been proven to be effective
  // altough they say so

  for (int i=-actsize_+1; i<=0; i+=matman_->NumVectors())
  {
    double a=0.0;
    double b=0.0;
    double aa=0.0;
    double bb=0.0;
    double beta=0.0;

    for (int j=0; j<matman_->NumVectors(); j++)
    {
      INVANA::MVDotProduct(*(*ystore_)(i+j),*(*sstore_)(i+j),*matman_->ParamLayoutMapUnique(),&a);
      INVANA::MVDotProduct(*(*ystore_)(i+j),*(*d_)(j),*matman_->ParamLayoutMapUnique(),&b);
      aa += a;
      bb += b;
    }

    beta=1/aa*bb;
    double alphac=alpha.back();
    alpha.pop_back();

    for (int j=0; j<matman_->NumVectors(); j++)
      (*d_)(j)->Update(alphac-beta, *(*sstore_)(i+j),1.0 );
  }

  // we do minimization not maximization
  d_->Scale(-1.0);

  return;
} // void ACOU::InvAnalysis::ComputeDirection()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::LineSearch()
{
  //Epetra_SerialDenseVector pnew(d.M());

  // This is a backtracking line search, starting with initial alpha and reducing it,
  // as long as the sufficient decrease condition is not yet satisfied

  double alpha_0 = 1.0;

  // does alpha_0 already satisfy the Armijo (sufficient decrease) condition?
  double alpha = alpha_0 * ls_gd_scal_;
  double J_before = J_;
  double dotproduct = 0.0;
  INVANA::MVDotProduct(*objgrad_,*d_,*matman_->ParamLayoutMapUnique(),&dotproduct);

  double condition = J_before + ls_c_ * alpha * dotproduct;

  double norm_d = 0.0;
  d_->Norm2(&norm_d);

  int count = 0;
  bool success = true;

  do{
    count++;

    if(!myrank_)
      std::cout<<"*************** line search iteration "<<count<<" of maximal "<<max_ls_iter_<<" line search iterations, condition "<<condition<<" step length "<<alpha<<" norm direction "<<norm_d<<std::endl;

    step_->Update(alpha,*d_,0.0);
    matman_->UpdateParams(step_);

    SolveStandardProblem();
    CalculateObjectiveFunctionValue();
    if( J_ < condition )
    {
      if( opti_ == INPAR::ACOU::inv_gd && count == 1)// we only needed one step within the line search
         ls_gd_scal_ /= ls_rho_;
      else if ( opti_ == INPAR::ACOU::inv_gd )
        ls_gd_scal_ *= ls_rho_;
      return success;
    }
    alpha *= ls_rho_;

    matman_->ResetParams();

  } while(J_ >= condition && count<max_ls_iter_);

  success = false;

  return success;
} // void ACOU::InvAnalysis::LineSearch()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateStatsAndService()
{
  // update iteration count
  iter_++;

  // update gradient
  objgrad_o_->Update(1.0,*objgrad_,0.0);

  // calculate the norm of the difference of the given parameters, for the curious user
  INVANA::MVNorm(*step_,*matman_->ParamLayoutMapUnique(),0,&normdiffp_);

  // update node based material parameter vector
  ComputeNodeBasedVectors();

  return;
} // void ACOU::InvAnalysis::CalculateStatsAndService()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::OutputStats()
{

  if (!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"*** objective function value:             "<<J_<<std::endl;
    std::cout<<"*** error value:                          "<<error_<<std::endl;
    std::cout<<"*** norm of the gradient:                 "<<normgrad_<<std::endl;
    std::cout<<"*** norm of the difference of parameters: "<<normdiffp_<<std::endl;
    std::cout<<"*** output count:                         "<<output_count_<<std::endl;
    std::cout<<"*** simulation time since start           "<<Teuchos::Time::wallTime()-tstart_<<std::endl;
    std::cout<<"*** parameters:                           "<<std::endl;
  }
  matman_->GetParams()->Print(std::cout);

  return;
} // void ACOU::InvAnalysis::OutputStats()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::ConvergenceCheck()
{
  // check convergence criteria
  bool unconverged;

  if(calcacougrad_)
    unconverged = error_>tol_ && iter_<max_iter_ && normgrad_>tol_grad_;
  else
    unconverged = error_>tol_ && iter_<max_iter_;

  if(std::abs(normdiffp_)<1.0e-14)
    unconverged = false;

  return unconverged;
} // bool ACOU::InvAnalysis::ConvergenceCheck()

/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_MultiVector> ACOU::InvAnalysis::ElementMatVec()
{
  return matman_->GetMatParams();
} // const Teuchos::RCP<Epetra_Vector> ACOU::InvAnalysis::ElementMatVec()

/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()
{
  return Teuchos::rcp(new AcouInvResultTest(*this));
} // Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()
