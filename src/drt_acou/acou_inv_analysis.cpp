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
#include "acou_expl.H"
#include "acou_impl_euler.H"
#include "acou_inv_resulttest.H"
#include "acou_ele_action.H"
#include "pat_matpar_manager.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_timestepping/timintmstep.H"
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
#include "../drt_inv_analysis/regularization_tvdtikh.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "Epetra_SerialDenseSolver.h"

#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*/
ACOU::InvAnalysis::InvAnalysis(Teuchos::RCP<DRT::Discretization> scatradis,
                               Teuchos::RCP<DRT::DiscretizationHDG> acoudis,
                               Teuchos::RCP<Teuchos::ParameterList> scatrapara,
                               Teuchos::RCP<LINALG::Solver> scatrasolv,
                               Teuchos::RCP<Teuchos::ParameterList> acoupara,
                               Teuchos::RCP<LINALG::Solver> acousolv,
                               Teuchos::RCP<IO::DiscretizationWriter> acouout
)
  : scatra_discret_(scatradis),
    acou_discret_(acoudis),
    scatraalgo_(Teuchos::null),
    acoualgo_(Teuchos::null),
    scatraparams_(scatrapara),
    scatrasolver_(scatrasolv),
    acouparams_(acoupara),
    acousolver_(acousolv),
    acououtput_(acouout),
    dyna_(DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*acouparams_,"TIMEINT")),
    phys_(DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*acouparams_,"PHYSICAL_TYPE")),
    opti_(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")),
    myrank_(acoudis->Comm().MyPID()),
    error_(1.0e6),
    error_start_(0.0),
    tol_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL")),
    iter_(0),
    max_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_MAX_RUN")),
    max_ls_iter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_LS_MAX_RUN")),
    output_count_(0),
    last_acou_out_(0),
    fdcheck_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("FDCHECK"))),
    ls_c_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("LS_DECREASECOND")),
    ls_rho_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("LS_STEPLENGTHRED")),
    ls_gd_scal_(1.0),
    backprojection_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("BACKPROJECTION"))),
    meshconform_(DRT::INPUT::IntegralValue<bool>(*acouparams_,"MESHCONFORM")),
    calcacougrad_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("INV_TOL_GRAD_YN"))),
    scalegradele_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("ELE_SCALING"))),
    segmentation_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("SEGMENTATION"))),
    penaltyparam_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("EQUALITYPENALTY")),
    adaptreguval_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("ADAPTREGUPARAMS")),
    reguweight_(1.0),
    scatra_normgrad_(1.0e6),
    acou_normgrad_(1.0e6),
    scatra_normgrad0_(1.0e6),
    acou_normgrad0_(1.0e6),
    dtacou_(acouparams_->get<double>("TIMESTEP")),
    J_(0.0),
    J_start_(0.0),
    normdiffp_(0.0),
    tstart_(Teuchos::Time::wallTime())
{
  // some checks, if everything is alright
  if(ls_rho_>=1.0) dserror("LS_STEPLENGTHRED has to be smaller than 1.0");
  if(ls_c_>1.0) dserror("LS_DECREASECOND is usually chosen in between 0.0 and 0.01, decrease it!");
  if(not acou_discret_->Filled() || not acou_discret_->HaveDofs()) dserror("acoustical discretization is not complete or has no dofs");
  if(not scatra_discret_->Filled() || not scatra_discret_->HaveDofs()) dserror("scatra discretization is not complete or has no dofs");
  if(segmentation_ && (meshconform_ == false)) dserror("segmentation requires conform mesh");
  if(segmentation_ && (opti_==INPAR::ACOU::inv_lbfgs)) dserror("segmentation only for gradientdescent not LBFGS");

  // create necessary extra parameter list for scatra
  scatraextraparams_ = Teuchos::rcp(new Teuchos::ParameterList());
  scatraextraparams_->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
  scatraextraparams_->set<bool>("isale",false);
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  scatraextraparams_->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  scatraextraparams_->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
  scatraextraparams_->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  scatraextraparams_->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

  // set up of the output
  scatraoutput_ = scatra_discret_->Writer();
  name_ = DRT::Problem::Instance()->OutputControlFile()->FileName();

  // vectors we need
  adjoint_psi_ = LINALG::CreateVector(*(acou_discret_->NodeRowMap()),true);
  phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  adjoint_phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  // tolerance for the gradient, if used
  if(calcacougrad_) tol_grad_ = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL_GRAD");
  else tol_grad_ = 0.0;

  // read monitor file, create multivector and map for measured values
  ReadMonitor(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE"));

  // allocate vectors
  node_mu_a_ = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);

  // and fill them
  ComputeNodeBasedVectors();

  // check for acoustic optimization
  {
    std::istringstream acoupstream(Teuchos::getNumericStringParameter(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"ACOUPARAMLIST"));
    std::string word; acoupstream>>word;

    if(word=="none")
      acouopt_ = false;
    else
      acouopt_ = true;
    acouparams_->set<bool>("acouopt",acouopt_);
  }
  if(segmentation_ && (acouopt_==false)) dserror("acoustic gradient has to be calculated for segmentation");

  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams(); // TODO override PARAMLIST for acoustic matman
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
    {
      scatra_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerPerElement(scatra_discret_,scalegradele_));
      if(acouopt_) acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerPerElement(acou_discret_));
    }
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
    {
      scatra_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerUniform(scatra_discret_));
      if(acouopt_) acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerUniform(acou_discret_));
    }
    break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
    break;
  }
  scatra_matman_->Init(invp);
  scatra_matman_->Setup();
  if(acouopt_)
  {
    Teuchos::ParameterList acouinvp(invp);
    acouinvp.set("PARAMLIST",Teuchos::getNumericStringParameter(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"ACOUPARAMLIST"));
    //acouinvp.set<INPAR::INVANA::StatInvMetaParamsType>("METAPARAMS",DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMetaParamsType>(invp,"METAPARAMS"));
    acou_matman_->Init(acouinvp);
    acou_matman_->Setup();
  }

  // create regularization manager
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
  case INPAR::INVANA::stat_inv_reg_none:
    break;
  case INPAR::INVANA::stat_inv_reg_tikhonov:
  {
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
    if(acouopt_) acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
    break;
  }
  case INPAR::INVANA::stat_inv_reg_totalvariation:
  {
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
    if(acouopt_) acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
    break;
  }
  case INPAR::INVANA::stat_inv_reg_tvdtikh:
  {
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
    if(acouopt_) acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
    break;
  }
  default:
    dserror("no valid regularization type provided");
    break;
  }
  if (scatra_regman_!=Teuchos::null)
  {
    scatra_regman_->Init(scatra_discret_,scatra_matman_->GetConnectivityData());
    scatra_regman_->Setup(invp);
    if(acouopt_)
    {
      acou_regman_->Init(acou_discret_,acou_matman_->GetConnectivityData());
      acou_regman_->Setup(invp);
    }
  }

  if(opti_==INPAR::ACOU::inv_lbfgs)
  {
    ssize_ =  invp.get<int>("SIZESTORAGE");
    ssize_ *= scatra_matman_->NumVectors();
    actsize_ = 0;

    scatra_sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, scatra_matman_->ParamLayoutMap().get(), true));
    scatra_ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, scatra_matman_->ParamLayoutMap().get(), true));
    if(acouopt_)
    {
      acou_sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, acou_matman_->ParamLayoutMap().get(), true));
      acou_ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(-ssize_+1, 0, acou_matman_->ParamLayoutMap().get(), true));
    }
  }

  // allocate the required quantities
  // scatra
  scatra_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
  scatra_objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
  scatra_step_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(), true));
  scatra_d_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(), true));
  // and acou
  if(acouopt_)
  {
    acou_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
    acou_objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
    acou_step_ = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(), true));
    acou_d_ = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(), true));
  }

  // check for segmentation
  if(segmentation_)
  {
    segmentationsolver_ = Teuchos::rcp(new Segmentation(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("SEGMENTATIONMATS"),penaltyparam_,scatra_discret_,acou_discret_));
    k_                  = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),segmentationsolver_->nummats,true));
    k_d_                = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),segmentationsolver_->nummats,true));
    k_objgrad_          = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),segmentationsolver_->nummats,true));
    segmentationsolver_->InitParameters(k_,scatra_matman_,acou_matman_,scatra_step_,acou_step_);
  }

  // create line search
  linesearch_ = Teuchos::rcp(new PATLineSearch(Teuchos::rcp(this,false),DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("SUCCESSIVEUPDATE"))));

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
  const std::vector<int> pressuremonmics = *(pressuremon[0]->Nodes());
  std::vector<int> pressuremonmicsunique;

  nodes_.resize(pressuremonmics.size());
  for(unsigned int i=0; i<pressuremonmics.size(); ++i)
    nodes_[i] = pressuremonmics[i];

  // create unique map
  acou_discret_->Comm().Barrier();
  for(int i=0; i<acou_discret_->Comm().NumProc(); ++i)
  {
    if(acou_discret_->Comm().MyPID() == i)
    {
      for(unsigned int j=0; j<pressuremonmics.size(); ++j)
      {
        if(acou_discret_->HaveGlobalNode(pressuremonmics[j]))
        {
          if(acou_discret_->gNode(pressuremonmics[j])->Owner()==int(i))
            pressuremonmicsunique.push_back(pressuremonmics[j]);
        }
      }
    }
    acou_discret_->Comm().Barrier();
  }
  acou_discret_->Comm().Barrier();

  //abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, abcnodes.size(), &abcnodes[0], 0, acou_discret_->Comm()));
  abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, pressuremonmicsunique.size(), &pressuremonmicsunique[0], 0, acou_discret_->Comm()));
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
  unsigned int nmics = 0;

  std::vector<std::vector<double> > meascoords;
  Epetra_SerialDenseVector mcurve;

  // get measured values
  // open monitor file and read it
  {
    // check for monitor file
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
    char* foundit = NULL;

    // read steps
    foundit = strstr(buffer,"steps"); foundit += strlen("steps");
    nsteps = strtol(foundit,&foundit,10);
    timesteps_.resize(nsteps);

    // read mics
    foundit = strstr(buffer,"mics"); foundit += strlen("mics");
    nmics = strtol(foundit,&foundit,10);

    // read measurement coordinates for every microphone
    meascoords.resize(nmics);
    for (unsigned int i=0; i<nmics; ++i)
    {
      meascoords[i].resize(3);
      fgets(buffer,150000,file);
      foundit = buffer;
      for(int j=0; j<3; ++j)
        meascoords[i][j] = strtod(foundit,&foundit);
    }

    // read in measured curve
    {
      mcurve = Epetra_SerialDenseVector(nmics*nsteps);

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
        timesteps_[i] = strtod(foundit,&foundit);
        for (unsigned int j=0; j<nmics; ++j)
          mcurve[count++] = strtod(foundit,&foundit);
        fgets(buffer,150000,file);
        foundit = buffer;
      }
      if (count != nmics*nsteps) dserror("Number of measured pressure values wrong on input");
    }
  }

  // interpolation
  // vector for the interpolated data
  nod_curve_interpol_=Epetra_SerialDenseVector(pressuremonmicsunique.size()*nsteps);

  if((pressuremonmicsunique.size())!=0)
  {
    unsigned int i=0, j=0, l=0;
    unsigned int must_set_curve=1;
    double help;
    double distance[nmics];
    double epsilon = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("EPSILON");

    //if the user doesn't want to give an espilon as input, we'll calculate it individually
    if(epsilon==-1.0)
      epsilon=GetEpsilon(pressuremonmicsunique.size());

    //interpolation-loop for every single nod
    for (i=0; i<pressuremonmicsunique.size(); ++i)
    {
      const double* nod_coords = acou_discret_->gNode(pressuremonmicsunique[i])->X();
      unsigned int M_1=0, M_2=0;

      for(j=0;j<nmics;j++)
      {
        distance[j] = Delta(meascoords[j][0],meascoords[j][1],meascoords[j][2],nod_coords[0],nod_coords[1],nod_coords[2]);
        // if the nod is in an epsilon bubble of any of the microphones, the measured curve of this microphone and the nod's curve should be equal
        if(distance[j]<=epsilon)
        {
          for(l=0;l<nsteps;l++)
            nod_curve_interpol_[i*nsteps+l]=mcurve[j+l*nmics];
          must_set_curve=0;
        }
      }
      // finds those two microphones, that are the nearest ones to the actual nod
      if(must_set_curve)
      {
        help=distance[0];
        for(j=0;j<nmics;j++)
        {
          if(distance[j]<help)
          {
            help=distance[j];
            M_1=j;
          }
        }
        if((M_1+1)==nmics)
        {
          help=distance[M_1-1];
          M_2=M_1-1;
        }
        else
        {
          help=distance[M_1+1];
          M_2=M_1+1;
        }
        for(j=0;j<nmics;j++)
        {
          if(j==M_1)
            ++j;
          if(distance[j]<help&&j<nmics)
          {
            help=distance[j];
            M_2=j;
          }
        }
        Interpol(nod_coords,meascoords, M_1,M_2, nmics,i, nsteps, mcurve,nod_curve_interpol_);
      }
    }
  }

  double eps = dtacou_/1000.0;
  if((numvec-1)*dtacou_>timesteps_[nsteps-1]+eps) dserror("You want to simulate till %.15f but your monitor file only provides values till %.15f! Fix it!",(numvec-1)*dtacou_,timesteps_[nsteps-1]);

  // every proc knows mcurve, now, we want to write mcurve to a Epetra_MultiVector in the same form as acou_rhs_
  // with the same parallel distribution!
  // and we want to interpolate measured values in case the monitored time step size is not the same as the one for the simulation
  acou_rhsm_->PutScalar(0.0);

  if( timesteps_[0] != 0.0 )
    dserror("your measured values have to start at time 0.0");
  else if( timesteps_[0] == 0.0 && timesteps_[1] == dtacou_ ) // the standard case
  {
    for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
      if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
        for(unsigned int j=0; j<nsteps; ++j)
          acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nod_curve_interpol_(i*nsteps+j)); // the proc who has this row, writes the value
  }
  else // we have to interpolate!
  {
    if( numvec < int(nsteps) )
    {
      dserror("set your time step size smaller, at least to %14f or implement here",timesteps_[1]-timesteps_[0]);
    }
    else
    {
      for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
      {
        if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
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
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,0.0);
            }
            else if(actualt<timesteps_[timeval]+eps && actualt>timesteps_[timeval]-eps) // then this is more or less it
            {
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nod_curve_interpol_(i*nsteps+timeval));
            }
            else
            {
              double value = nod_curve_interpol_(i*nsteps+(timeval-1)) + (nod_curve_interpol_(i*nsteps+(timeval))-nod_curve_interpol_(i*nsteps+(timeval-1))) * (actualt - timesteps_[timeval-1]) / (timesteps_[timeval]-timesteps_[timeval-1]);
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,value);
            }
          } // for(int j=0; j<numvec; ++j)
        } // if( acou_discret_->HaveGlobalNode(nodes_[i]) )
      } // for(unsigned int i=0; i<nnodes; ++i)
      acou_discret_->Comm().Barrier();
    } // else ** if( numvec < nsteps_ )
  } // else ** if( timesteps_[0] == dtacou || (timesteps_[0]==0.0 && timesteps_[1] = dtacou) )

  return;
} // void ACOU::InvAnalysis::ReadMonitor(...)

/*----------------------------------------------------------------------*/
double ACOU::InvAnalysis::Delta(double coord_M_x,double coord_M_y, double coord_M_z, double coord_N_x,double coord_N_y, double coord_N_z)
{
  double distance = sqrt((coord_M_x-coord_N_x)*(coord_M_x-coord_N_x)+(coord_M_y-coord_N_y)*(coord_M_y-coord_N_y)+(coord_M_z-coord_N_z)*(coord_M_z-coord_N_z));
  return distance;
}

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::Interpol(const double nod_coords[3],std::vector<std::vector<double> > mic_coords, unsigned int mic_1, unsigned int mic_2, int nmic, unsigned int nod, int timesteps_, Epetra_SerialDenseVector& curve, Epetra_SerialDenseVector& inter_curve)
{
  double d1;
  double D1;
  double d2;
  double D2;

  d1=Delta(mic_coords[mic_1][0],mic_coords[mic_1][1],mic_coords[mic_1][2],nod_coords[0],nod_coords[1],nod_coords[2]);
  d2=Delta(mic_coords[mic_2][0],mic_coords[mic_2][1],mic_coords[mic_2][2],nod_coords[0],nod_coords[1],nod_coords[2]);
  D2=d2/(d1+d2);
  D1=d1/(d1+d2);

  for(int i=0;i<timesteps_;i++)
  {
  inter_curve[nod*timesteps_+i]=D2*curve[mic_1+i*nmic]+D1*curve[mic_2+i*nmic];
  }
}

/*----------------------------------------------------------------------*/
double ACOU::InvAnalysis::GetEpsilon(int nnodes)
{
  double min_dis[nnodes];
  double dist;
  double min_abs;
  double eps;

  //creates a vector which contains the distance of every nod to its nearest neighbor
  for(int i=0;i<abcnodes_map_->NumMyElements();i++)
  {
    const double* nc= acou_discret_->gNode(abcnodes_map_->GID(i))->X();//acou_discret_->gNode(nodes_[i])->X();
    int iter=0;
    for(int j=0; j<nnodes; j++)
    {
      if(j==i)
        j++;
      if(j==nnodes)
        break;
      const double* ncc=acou_discret_->gNode(abcnodes_map_->GID(j))->X();
      dist = sqrt((nc[0]-ncc[0])*(nc[0]-ncc[0])+(nc[1]-ncc[1])*(nc[1]-ncc[1])+(nc[2]-ncc[2])*(nc[2]-ncc[2]));
      if(iter==0)
        min_dis[i]=dist;
      else if(dist<min_dis[i])
        min_dis[i]=dist;
      ++iter;
    }
  }

  //searches for the (absolute) smallest distance
  min_abs = min_dis[0];
  for(int i=0; i<nnodes; i++)
    if(min_abs>min_dis[i])
      min_abs=min_dis[i];

  eps = min_abs*10.0;

  return eps;
}

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ComputeNodeBasedVectors()
{
  // fill the vectors
  int minnodeidscatra = scatra_discret_->NodeRowMap()->MinAllGID();

  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd) // cannot loop scatra nodes, because they don't necessarily start with gid 0
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+minnodeidscatra) )
    {
      opti_node = scatra_discret_->gNode(nd+minnodeidscatra);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;
    // here, every proc knows the owner and the gid of the optical node

    // now, every proc knows the speed of sound and density of this node -> write them to vector
    int nodelid = scatra_discret_->NodeRowMap()->LID(nd+minnodeidscatra);

    // we have to do the same procedure for the absorption coefficient with the scatra discretization
    int loc_numoptiele = 0;
    double loc_mu_a = 0.0;
    for(int roel = 0; roel<scatra_discret_->NumMyRowElements(); ++roel)
    {
      DRT::Element* roptele = scatra_discret_->lRowElement(roel);
      const int* nodeids = roptele->NodeIds();
      int numnode = roptele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd+minnodeidscatra )
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
  // Estimate values for mu_a
  if(backprojection_)
    EstimateMua(0);

  // Solve the standard problem
  SolveStandardProblemScatra();

  if(backprojection_)
  {
    EstimateMua(1);
    SolveStandardProblemScatra(); //to output the actual distribution
  }
  SolveStandardProblemAcou();

  // Calculate the error and the value of the objective function
  CalculateObjectiveFunctionValue(true);

  // save the initial objective function value
  J_start_ = J_;
  error_start_ = error_;

  // Solve the adjoint problem
  SolveAdjointAcouProblem();
  SolveAdjointOptiProblem();

  // Build the gradient
  CalculateGradient();

  // init
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
  else if(scatra_normgrad_<= tol_grad_ && !myrank_)
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
void ACOU::InvAnalysis::SolveStandardProblemScatra()
{
  // output for user
  scatra_discret_->Comm().Barrier();
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SCALAR TRANSPORT PROBLEM - OPTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  // create and run scatra algorithm
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*scatraparams_,"VELOCITYFIELD");
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
      scatraoutput_->NewResultFile(outname,output_count_);
      output_count_++;
      scatraoutput_->WriteMesh(0,0.0);

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      scatraalgo_ = Teuchos::rcp(new SCATRA::TimIntStationary(scatra_discret_, scatrasolver_, scatraparams_, scatraextraparams_, scatraoutput_));

      scatraalgo_->Init();
      scatraalgo_->SetVelocityField();

      scatraalgo_->TimeLoop();

      // output of elemental reaction coefficient
      OutputReactionAndDiffusion();

      // store the solution vector
      phi_ = scatraalgo_->Phinp();

      break;
    }
    default:
      dserror("unknown velocity field type for transport of passive scalar in problem type Acoustics");
      break;
  }
  return;
}

void ACOU::InvAnalysis::OutputReactionAndDiffusion()
{
  // build the two vectors
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> diffvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    double diff = actele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
    diffvec->operator [](i) = diff;
  }
  scatraoutput_->WriteVector("rea_coeff",reacvec);
  scatraoutput_->WriteVector("diff_coeff",diffvec);
  if(segmentation_)
  {
    for(unsigned int m=0; m<segmentationsolver_->nummats; ++m)
    {
      std::ostringstream nameos;
      nameos<< "k_" << m+1;
      Teuchos::RCP<Epetra_Vector> kvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
      for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
        kvec->operator [](i) = 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;
      scatraoutput_->WriteVector(nameos.str(),kvec);
    }
  }
  return;
}

void ACOU::InvAnalysis::SolveStandardProblemAcou()
{
  acouparams_->set<bool>("adjoint",false);

  std::string outname = name_;
  outname.append("_invforward_acou");
  acououtput_->NewResultFile(outname,output_count_);
  last_acou_out_ = output_count_;
  output_count_++;

  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }
  acoualgo_->SetInitialPhotoAcousticField(phi_,scatra_discret_,meshconform_);

  // we have to call a slightly changed routine, which fills our history vector which we need for the adjoint problem
  acou_rhs_->Scale(0.0);

  // determine if we have to perform this evaluation
  double maxval = 0.0;
  scatra_matman_->GetParams()->MaxValue(&maxval);
  if(scatra_matman_->NumVectors()==0)
    maxval = 1.0; // in case only acoustic parameters are optimized!

  // only necessary for values which are nonzero
  if(maxval>1.0e-9) // we do not have to solve when the maximum reaction coefficient is zero, just a waste of time!
    acoualgo_->Integrate(acou_rhs_,abcnodes_mapex_);

  return;
} // void ACOU::InvAnalysis::SolveStandardProblem()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateObjectiveFunctionValue(bool first)
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

  error_ = 0.5 * normvec.Norm1();

  // contribution from regularization
  if(scatra_regman_ != Teuchos::null)
  {
    scatra_regman_->Evaluate(*(scatra_matman_->GetParams()),&J_);
    if(acouopt_) acou_regman_->Evaluate(*(acou_matman_->GetParams()),&J_);
  }

  // contribution from error
  if( adaptreguval_==0.0 || scatra_regman_==Teuchos::null)
    J_ += error_;
  else
  {
    if(first)
    {
      reguweight_ = adaptreguval_ * error_/J_;
      if(iter_==0)
        reguweight_ = adaptreguval_;
    }
    J_ *= reguweight_;
    J_ += error_;
    if(!myrank_ && first)
      std::cout<<"set reg weight to "<<reguweight_<<std::endl;
  }

  // contribution from penalty term
  if(segmentation_)
  {
    Teuchos::RCP<Epetra_Vector> kvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
    kvec->PutScalar(-1.0);
    for(unsigned int m=0; m<segmentationsolver_->nummats; ++m)
      for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
        kvec->operator [](i) += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;
    double val = 0.0;
    kvec->Norm2(&val);
    J_ += 0.5 * penaltyparam_ * val * val;
  }

  if(!myrank_)
  {
    std::cout.precision(15);
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<std::endl;
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
  if(acouopt_)
  {
    acouparams_->set<int>("outputcount",last_acou_out_);
    acouparams_->set<std::string>("name",name_);
  }

  // build difference
  Teuchos::RCP<Epetra_MultiVector> tempvec = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors(),true));
  tempvec->Update(1.0,*acou_rhs_,0.0);
  tempvec->Update(-1.0,*acou_rhsm_,1.0);
  //acou_rhs_->Update(-1.0,*acou_rhsm_,1.0);

  // acou_rhs_ has to be scaled with weighting (adjoint of the mappint (|R))
  Teuchos::RCP<Epetra_Vector> touchcountvec = LINALG::CreateVector(*abcnodes_map_);
  acoualgo_->FillTouchCountVec(touchcountvec);
  tempvec->Multiply(1.0,*touchcountvec,*tempvec,0.0);

  // set the difference between measured and simulated values
  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",tempvec);

  // prepare the output
  std::string outname = name_;
  outname.append("_invadjoint_acou");
  acououtput_->NewResultFile(outname,output_count_);
  output_count_++;

  // create the acoustic algorithm
  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
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

  // give me psi which is needed for the source term of the adjoint optical problem
  adjoint_psi_->PutScalar(0.0);
  acoualgo_->NodalPsiField(adjoint_psi_);

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

  // consider both cases
  if(meshconform_)
  {
    rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

    int minnodeidscatra = scatra_discret_->NodeRowMap()->MinAllGID();
    for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
    {
      // get node and owner
      int myoptnodeowner = -1;
      int optnodeowner = -1;
      DRT::Node* opti_node = NULL;
      if( scatra_discret_->HaveGlobalNode(nd+minnodeidscatra) )
      {
        opti_node = scatra_discret_->gNode(nd+minnodeidscatra);
        myoptnodeowner = opti_node->Owner();
        if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
      }
      scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
      if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
        continue;

      double loc_value = 0.0;
      if(acou_discret_->NodeRowMap()->LID(nd)>-1)
      {
        loc_value = adjoint_psi_->operator [](acou_discret_->NodeRowMap()->LID(nd));
      }
      double glo_value = 0.0;
      acou_discret_->Comm().SumAll(&loc_value,&glo_value,1);

      // ok, we got the value, we still need c, rho and mu_a, but they are stored on the nodemap of the scatra dis, so this should not be a problem
      if(scatra_discret_->Comm().MyPID() == optnodeowner)
      {
        int lnodeid = node_mu_a_->Map().LID(nd+minnodeidscatra);
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
    rhsvec = CalculateNonconformRhsvec(adjoint_psi_);
    for(int i=0; i<node_mu_a_->MyLength(); ++i)
    {
      double mu_a = node_mu_a_->operator [](i);
      int dofgid = scatra_discret_->Dof(scatra_discret_->lRowNode(i),0);
      int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
      rhsvec->operator [](doflid) *= -mu_a;
    }
  } // else ** if(meshconform)

  // perform the element integration
  Teuchos::ParameterList eleparams;
  scatra_discret_->SetState("rhsnodebasedvals",rhsvec);
  eleparams.set<int>("action",SCATRA::calc_integr_pat_rhsvec);
  Teuchos::RCP<Epetra_Vector> b = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,b,Teuchos::null,Teuchos::null);

  // consider Dirichlet boundaries in the right hand side vector
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
    }
  }

  // solve the system now!
  scatrasolver_->Solve(sysmatscatra->EpetraOperator(),adjoint_phi_,b,true,true);

  // output the solution
  std::string outname = name_;
  outname.append("_invadjoint_opti");
  scatraoutput_->NewResultFile(outname,output_count_);
  output_count_++;
  scatraoutput_->WriteMesh(0,0.0);
  scatraoutput_->NewStep(1,1.0);
  scatraoutput_->WriteElementData(true);
  scatraoutput_->WriteVector("phinp",adjoint_phi_);
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
  /* THE SCATRA GRADIENT PART */
  {
    // gradient consists of two parts: eintries with respect to the absorption coefficient
    // and entries with respect to the diffusion coefficient, both are treated in a similar way

    //zero out gradient vector initially
    scatra_objgrad_->Scale(0.0);

    // contribution from regularization
    if(scatra_regman_ != Teuchos::null)
      scatra_regman_->EvaluateGradient(*(scatra_matman_->GetParams()),scatra_objgrad_);
    scatra_objgrad_->Scale(reguweight_); // weight is 1.0 in case it is not adapted

    scatra_discret_->SetState("adjoint phi",adjoint_phi_);
    scatra_discret_->SetState("phi",phi_);

    Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

    if(meshconform_)
    {
      int minscatranodegid = scatra_discret_->NodeRowMap()->MinAllGID();
      for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
      {
        // get node and owner
        int myoptnodeowner = -1;
        int optnodeowner = -1;
        DRT::Node* opti_node = NULL;
        if( scatra_discret_->HaveGlobalNode(nd+minscatranodegid) )
        {
          opti_node = scatra_discret_->gNode(nd+minscatranodegid);
          myoptnodeowner = opti_node->Owner();
          if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
        }
        scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
        if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
          continue;

        double loc_value = 0.0;
        if(acou_discret_->NodeRowMap()->LID(nd)>-1)
        {
          loc_value = adjoint_psi_->operator [](acou_discret_->NodeRowMap()->LID(nd));
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
      psi = CalculateNonconformRhsvec(adjoint_psi_);
    // set the psi vector
    scatra_discret_->SetState("psi",psi);

    // do the actual evaluation
    scatra_matman_->AddEvaluate(0.0,scatra_objgrad_);
    scatra_matman_->Finalize(scatra_objgrad_);

    INVANA::MVNorm(*scatra_objgrad_,*scatra_matman_->ParamLayoutMapUnique(),2,&scatra_normgrad_);
    scatra_discret_->ClearState();
    //TODO
    //scatra_objgrad_->Scale(1.0/scatra_normgrad_);
    //INVANA::MVNorm(*scatra_objgrad_,*scatra_matman_->ParamLayoutMapUnique(),2,&scatra_normgrad_);
  }

  /* THE ACOUSTIC GRADIENT PART */
  if(acouopt_)
  {
    acou_objgrad_->Scale(0.0);

    // contribution from regularization
    if(acou_regman_ != Teuchos::null)
    {
      acou_regman_->EvaluateGradient(*(acou_matman_->GetParams()),acou_objgrad_);
      acou_objgrad_->Scale(reguweight_);
    }

    // normal contribution
    acou_matman_->AddEvaluate(0.0,acou_objgrad_);
    acou_matman_->Finalize(acou_objgrad_);
    INVANA::MVNorm(*acou_objgrad_,*acou_matman_->ParamLayoutMapUnique(),2,&acou_normgrad_);
    //TODO
    //acou_objgrad_->Scale(1.0/acou_normgrad_);
    //INVANA::MVNorm(*acou_objgrad_,*acou_matman_->ParamLayoutMapUnique(),2,&acou_normgrad_);
    //std::cout<<scatra_normgrad_<<" "<<acou_normgrad_<<std::endl;
  }

  // set initial norm of gradient
  if( iter_ == 0)
  {
    scatra_normgrad0_ = scatra_normgrad_;
    acou_normgrad0_ = acou_normgrad_;
    std::cout<<"Gradient Norm: "<<scatra_normgrad0_<<" "<<acou_normgrad0_<<std::endl;
  }

  if(segmentation_)
    segmentationsolver_->CalculateGradient(
                             scatra_matman_->GetParams(),
                             acou_matman_->GetParams(),
                             scatra_objgrad_,
                             acou_objgrad_,
                             k_,
                             k_objgrad_);

  //TODO
  /*if(iter_==0)
  for(int i=0; i<f_objgrad_->MyLength(); ++i)
    if(f_objgrad_->operator [](i)>0.0)
      f_objgrad_->operator [](i)= 0.0;
*/
  // check gradient if required
  if(fdcheck_)
    FD_GradientCheck();

  return;
} // void ACOU::InvAnalysis::CalculateGradient()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::FD_GradientCheck()
{
  if(myrank_)
    dserror("The Gradient Check works only for serial calculations");

  std::cout<<"ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"<<std::endl;
  std::cout<<"ooooooooooooooooooooo FD Check oooooooooooooooooooooooooo"<<std::endl;
  std::cout<<"ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"<<std::endl;

  double J_before = J_;
  Epetra_MultiVector tempvec = Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors());
  tempvec.Update(1.0,*acou_rhs_,0.0);

  /* THE SCATRA GRADIENT PART */
  if(0)
  {
    std::cout<<"Scatra gradient according to adjoint analysis"<<std::endl;
    scatra_objgrad_->Print(std::cout);
    int numparams = scatra_matman_->NumVectors();

    Teuchos::RCP<Epetra_MultiVector> objgradFD = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));

    double perturba = 1.0e-3;
    double perturbb = 0.0; //1.0e-4;

    Teuchos::RCP<Epetra_MultiVector> perturb = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
    perturb->Update(1.0,*(scatra_matman_->GetParams()),0.0);

    //keep a copy of the current parameters to reset after perturbation:
    Teuchos::RCP<Epetra_MultiVector> pcurr = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
    pcurr->Update(1.0,*(scatra_matman_->GetParams()),0.0);

    double pn=0.0;
    double p=0.0;
    double dp=0.0;
    for (int i=0; i<numparams; i++)
    {
      for (int j=0; j<scatra_matman_->ParamLayoutMap()->NumGlobalElements(); j++)
      {
        int lid = scatra_matman_->GetParams()->Map().LID(j+scatra_matman_->ParamLayoutMap()->MinAllGID());
        if(lid > -1)   p = (*(*(scatra_matman_->GetParams()))(i))[lid];
        else           p = 0.0;

        pn = p+p*perturba+perturbb;
        std::cout<<"i "<<i<<" j "<<j<<" p "<<p<<" disturbed "<<pn<<std::endl;
        perturb->ReplaceGlobalValue(j+scatra_matman_->ParamLayoutMap()->MinAllGID(),i,pn);
        //perturb->Print(std::cout);
        scatra_matman_->ReplaceParams(*perturb);

        SolveStandardProblemScatra();
        SolveStandardProblemAcou();

        perturb->Update(1.0,*pcurr,0.0);

        CalculateObjectiveFunctionValue(false);
        dp=(J_before-J_)/(p-pn);
        std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
        objgradFD->ReplaceGlobalValue(j+scatra_matman_->ParamLayoutMap()->MinAllGID(),i,dp);
      }
    }
    std::cout<<"Scatra gradient according to finite difference check"<<std::endl;
    objgradFD->Print(std::cout);

    // we pretend this never happened
    scatra_matman_->ReplaceParams(*pcurr);
  }

  /* THE ACOUSTIC GRADIENT PART */
  if(0)   //if(acouopt_)
  {
    std::cout<<"Acoustic gradient according to adjoint analysis"<<std::endl;
    acou_objgrad_->Print(std::cout);
    int numparams = acou_matman_->NumVectors();

    Teuchos::RCP<Epetra_MultiVector> objgradFD = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));

    double perturba = 1.0e-3;
    double perturbb = 0.0; //1.0e-4;

    Teuchos::RCP<Epetra_MultiVector> perturb = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
    perturb->Update(1.0,*(acou_matman_->GetParams()),0.0);

    //keep a copy of the current parameters to reset after perturbation:
    Teuchos::RCP<Epetra_MultiVector> pcurr = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
    pcurr->Update(1.0,*(acou_matman_->GetParams()),0.0);

    double pn=0.0;
    double p=0.0;
    double dp=0.0;
    for (int i=0; i<numparams; i++)
    {
      for (int j=0; j<acou_matman_->ParamLayoutMap()->NumGlobalElements(); j++)
      {
        int lid = acou_matman_->GetParams()->Map().LID(j+acou_matman_->ParamLayoutMap()->MinAllGID());
        if(lid > -1)   p = (*(*(acou_matman_->GetParams()))(i))[lid];
        else           p = 0.0;

        pn = p+p*perturba+perturbb;
        std::cout<<"i "<<i<<" j "<<j<<" p "<<p<<" disturbed "<<pn<<std::endl;
        perturb->ReplaceGlobalValue(j+acou_matman_->ParamLayoutMap()->MinAllGID(),i,pn);
        //perturb->Print(std::cout);
        acou_matman_->ReplaceParams(*perturb);

        SolveStandardProblemScatra();
        SolveStandardProblemAcou();

        perturb->Update(1.0,*pcurr,0.0);

        CalculateObjectiveFunctionValue(false);
        dp=(J_before-J_)/(p-pn);
        std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
        objgradFD->ReplaceGlobalValue(j+acou_matman_->ParamLayoutMap()->MinAllGID(),i,dp);
      }
    }
    std::cout<<"Acoustic gradient according to finite difference check"<<std::endl;
    objgradFD->Print(std::cout);

    // we pretend this never happened
    acou_matman_->ReplaceParams(*pcurr);
  }

  /* THE GRADIENT WITH RESPECT TO THE CONCENTRATIONS */
  if(segmentation_)
  {
    std::cout<<"Concentration gradient according to adjoint analysis"<<std::endl;
    k_objgrad_->Print(std::cout);

    Teuchos::RCP<Epetra_MultiVector> objgradFD = Teuchos::rcp(new Epetra_MultiVector(k_objgrad_->Map(),k_objgrad_->NumVectors(),true));

    double perturba = 1.0e-3;
    double perturbb = 0.0; //1.0e-4;

    Teuchos::RCP<Epetra_MultiVector> perturb = Teuchos::rcp(new Epetra_MultiVector(k_objgrad_->Map(),k_objgrad_->NumVectors(),true));
    perturb->Update(1.0,*k_,0.0);

    //keep a copy of the current parameters to reset after perturbation:
    Teuchos::RCP<Epetra_MultiVector> pcurr = Teuchos::rcp(new Epetra_MultiVector(k_objgrad_->Map(),k_objgrad_->NumVectors(),true));
    pcurr->Update(1.0,*k_,0.0);

    double pn=0.0;
    double p=0.0;
    double dp=0.0;
    for(unsigned int m=0; m<segmentationsolver_->nummats; ++m)
    {
      for(int j=0; j<k_objgrad_->MyLength(); ++j)
      {
        p = pcurr->operator ()(m)->operator [](j);
        pn = p+p*perturba+perturbb;

        std::cout<<"m "<<m<<" j "<<j<<" p "<<p<<" disturbed "<<pn<<std::endl;
        perturb->operator ()(m)->operator [](j) = pn;
        k_->Update(1.0,*perturb,0.0);
        segmentationsolver_->UpdateParameters(perturb,scatra_matman_,acou_matman_,scatra_step_,acou_step_);
        scatra_matman_->ReplaceParams(*scatra_step_);
        acou_matman_->ReplaceParams(*acou_step_);

        SolveStandardProblemScatra();
        SolveStandardProblemAcou();

        perturb->Update(1.0,*pcurr,0.0);

        CalculateObjectiveFunctionValue(false);
        dp=(J_before-J_)/(p-pn);
        std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn meta"<<p-pn<<" val "<<dp<<std::endl;
        objgradFD->operator ()(m)->operator [](j) = dp;
      }
    }
    std::cout<<"Concentration gradient according to finite difference check"<<std::endl;
    objgradFD->Print(std::cout);
  }

  // we pretend this never happened
  J_ = J_before;
  acou_rhs_->Update(1.0,tempvec,0.0);
  dserror("stop it");
  return;
} // void ACOU::InvAnalysis::FD_GradientCheck()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::UpdateParameters()
{
  double J_before = J_;
  Teuchos::RCP<Epetra_MultiVector> optparams_before = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
  Teuchos::RCP<Epetra_MultiVector> acouparams_before;
  if(acouopt_)
    acouparams_before = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
  optparams_before->Update(1.0,*scatra_matman_->GetParams(),0.0);
  if(acouopt_) acouparams_before->Update(1.0,*acou_matman_->GetParams(),0.0);

  // first iteration is steepest descent
  if(segmentation_)
  {
    k_d_->Update(-1.0,*k_objgrad_,0.0);
    double norms[k_d_->NumVectors()];
    k_d_->Norm1(norms);
    double norm = 0.0;
    for(unsigned int i=0;i<segmentationsolver_->nummats;++i)
      norm+=norms[i];
    norm/=double(segmentationsolver_->nummats);
    norm/=double(k_d_->GlobalLength());
    norm*=10.0;
    k_d_->Scale(1./norm);
  }
  else if(iter_==0 || opti_ == INPAR::ACOU::inv_gd)
  {
    scatra_d_->Update(-1.0,*scatra_objgrad_,0.0);
    if(acouopt_)
      acou_d_->Update(-1.0,*acou_objgrad_,0.0);
  }
  else
  {
    // else the inverse hessian is applied
    StoreVectors();
    ComputeDirection();
  }

  if( (segmentation_== false && (iter_==0 || opti_ == INPAR::ACOU::inv_gd)) ) // in the first iteration the lbfgs does not know how long the step should be, gradient descent should always do scaled step length
  {
    //scatra_d_->Scale(1.0/scatra_normgrad0_*double(scatra_objgrad_->Map().NumGlobalElements()));
    double max_val = 0.0;
    scatra_d_->MaxValue(&max_val);
    if(max_val != 0.0)
      scatra_d_->Scale(1.0/max_val/20.0);

    if(acouopt_)
    {
      //acou_d_->Scale(1.0/acou_normgrad0_*double(acou_objgrad_->Map().NumGlobalElements()));
      double min_val = 0.0;
      acou_d_->MinValue(&min_val);
      acou_d_->MaxValue(&max_val);
      if(min_val<0.0)
      {
        std::cout<<"MIN VAL IN ACOU D WARNING "<<min_val<<std::endl;
        acou_d_->Scale(1.0/std::abs(min_val)*0.2);
      }
      else if(max_val != 0.0)
        acou_d_->Scale(1.0/max_val/20.0);
    }
  }

  // output for small problems
  if(scatra_d_->GlobalLength()<30)
  {
    std::cout<<"search direction "<<std::endl;
    scatra_d_->Print(std::cout);
    if(acouopt_)
      acou_d_->Print(std::cout);
  }

  // update gradient
  scatra_objgrad_o_->Update(1.0,*scatra_objgrad_,0.0);
  if(acouopt_)
    acou_objgrad_o_->Update(1.0,*acou_objgrad_,0.0);

  linesearch_->Init();
  bool success = linesearch_->Run();

  if (success == false && opti_ == INPAR::ACOU::inv_lbfgs && iter_ != 0)
  {
    // in this case we try a steepest descent step!
    scatra_d_->Update(-1.0/scatra_normgrad0_,*scatra_objgrad_,0.0);
    if(acouopt_)
      acou_d_->Update(-1.0/acou_normgrad0_,*acou_objgrad_,0.0);

    if(!myrank_)
      std::cout<<"line search failed, we try steepest descent from now on!"<<std::endl;

    // reset objective function
    J_ = J_before;

    // and in every following step we do steepest descent
    opti_ = INPAR::ACOU::inv_gd;
    linesearch_->Init();
    success = linesearch_->Run();
  }
  return success;
} // void ACOU::InvAnalysis::UpdateParameters()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::StoreVectors()
{
  if (iter_*scatra_matman_->NumVectors()<=ssize_) // we have "<=" since we do not store the first run
    actsize_+=scatra_matman_->NumVectors();

  /* THE SCATRA PART */
  {
    Epetra_MultiVector s(*(scatra_matman_->ParamLayoutMap()), (scatra_matman_->NumVectors()),true);

    //push back s
    s.Update(1.0,*(scatra_matman_->GetParams()),-1.0,*(scatra_matman_->GetParamsOld()),0.0);
    for (int i=0; i<s.NumVectors(); i++)
      scatra_sstore_->UpdateSteps(*s(i));

    // push back y
    s.Update(1.0,*scatra_objgrad_,-1.0,*scatra_objgrad_o_,0.0);
    for (int i=0; i<s.NumVectors(); i++)
      scatra_ystore_->UpdateSteps(*s(i));
  }

  /* THE ACOUSTIC PART */
  if(acouopt_)
  {
    Epetra_MultiVector s(*(acou_matman_->ParamLayoutMap()), (acou_matman_->NumVectors()),true);

    //push back s
    s.Update(1.0,*(acou_matman_->GetParams()),-1.0,*(acou_matman_->GetParamsOld()),0.0);
    for (int i=0; i<s.NumVectors(); i++)
      acou_sstore_->UpdateSteps(*s(i));

    // push back y
    s.Update(1.0,*acou_objgrad_,-1.0,*acou_objgrad_o_,0.0);
    for (int i=0; i<s.NumVectors(); i++)
      acou_ystore_->UpdateSteps(*s(i));
  }

  return;
} // void ACOU::InvAnalysis::StoreVectors()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::ComputeDirection()
{
  /* THE SCATRA PART */
  {
    scatra_d_->Update(1.0,*scatra_objgrad_,0.0);
    std::vector<double> alpha;

    // loop steps
    for (int i=0; i>-actsize_; i-=scatra_matman_->NumVectors())
    {
      double a=0.0;
      double b=0.0;
      double aa=0.0;
      double bb=0.0;

      int ind=0;
      for (int j=scatra_matman_->NumVectors(); j>0; j--)
      {
        INVANA::MVDotProduct(*(*scatra_ystore_)(i-j+1),*(*scatra_sstore_)(i-j+1),*scatra_matman_->ParamLayoutMapUnique(),&a);
        INVANA::MVDotProduct(*(*scatra_sstore_)(i-j+1),*(*scatra_d_)(ind),*scatra_matman_->ParamLayoutMapUnique(),&b);
        ind++;
        aa += a;
        bb += b;
      }
      alpha.push_back(1/aa*bb);

      ind=0;
      for (int j=scatra_matman_->NumVectors(); j>0; j--)
      {
        (*scatra_d_)(ind)->Update(-1.0*alpha.back(), *(*scatra_ystore_)(i-j+1),1.0 );
        ind++;
      }
    }

    // Some scaling of the initial hessian might come in here but has not been proven to be effective
    // altough they say so

    for (int i=-actsize_+1; i<=0; i+=scatra_matman_->NumVectors())
    {
      double a=0.0;
      double b=0.0;
      double aa=0.0;
      double bb=0.0;
      double beta=0.0;

      for (int j=0; j<scatra_matman_->NumVectors(); j++)
      {
        INVANA::MVDotProduct(*(*scatra_ystore_)(i+j),*(*scatra_sstore_)(i+j),*scatra_matman_->ParamLayoutMapUnique(),&a);
        INVANA::MVDotProduct(*(*scatra_ystore_)(i+j),*(*scatra_d_)(j),*scatra_matman_->ParamLayoutMapUnique(),&b);
        aa += a;
        bb += b;
      }

      beta=1/aa*bb;
      double alphac=alpha.back();
      alpha.pop_back();

      for (int j=0; j<scatra_matman_->NumVectors(); j++)
        (*scatra_d_)(j)->Update(alphac-beta, *(*scatra_sstore_)(i+j),1.0 );
    }

    // we do minimization not maximization
    scatra_d_->Scale(-1.0);
  }
  /* THE ACOUSTIC PART */
  if(acouopt_)
  {
    acou_d_->Update(1.0,*acou_objgrad_,0.0);
    std::vector<double> alpha;

    // loop steps
    for (int i=0; i>-actsize_; i-=acou_matman_->NumVectors())
    {
      double a=0.0;
      double b=0.0;
      double aa=0.0;
      double bb=0.0;

      int ind=0;
      for (int j=acou_matman_->NumVectors(); j>0; j--)
      {
        INVANA::MVDotProduct(*(*acou_ystore_)(i-j+1),*(*acou_sstore_)(i-j+1),*acou_matman_->ParamLayoutMapUnique(),&a);
        INVANA::MVDotProduct(*(*acou_sstore_)(i-j+1),*(*acou_d_)(ind),*acou_matman_->ParamLayoutMapUnique(),&b);
        ind++;
        aa += a;
        bb += b;
      }
      alpha.push_back(1/aa*bb);

      ind=0;
      for (int j=acou_matman_->NumVectors(); j>0; j--)
      {
        (*acou_d_)(ind)->Update(-1.0*alpha.back(), *(*acou_ystore_)(i-j+1),1.0 );
        ind++;
      }
    }

    // Some scaling of the initial hessian might come in here but has not been proven to be effective
    // altough they say so

    for (int i=-actsize_+1; i<=0; i+=acou_matman_->NumVectors())
    {
      double a=0.0;
      double b=0.0;
      double aa=0.0;
      double bb=0.0;
      double beta=0.0;

      for (int j=0; j<acou_matman_->NumVectors(); j++)
      {
        INVANA::MVDotProduct(*(*acou_ystore_)(i+j),*(*acou_sstore_)(i+j),*acou_matman_->ParamLayoutMapUnique(),&a);
        INVANA::MVDotProduct(*(*acou_ystore_)(i+j),*(*acou_d_)(j),*acou_matman_->ParamLayoutMapUnique(),&b);
        aa += a;
        bb += b;
      }

      beta=1/aa*bb;
      double alphac=alpha.back();
      alpha.pop_back();

      for (int j=0; j<acou_matman_->NumVectors(); j++)
        (*acou_d_)(j)->Update(alphac-beta, *(*acou_sstore_)(i+j),1.0 );
    }

    // we do minimization not maximization
    acou_d_->Scale(-1.0);
  }

  return;
} // void ACOU::InvAnalysis::ComputeDirection()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::CalculateStatsAndService()
{
  // update iteration count
  iter_++;

  // calculate the norm of the difference of the given parameters, for the curious user
  INVANA::MVNorm(*scatra_step_,*scatra_matman_->ParamLayoutMapUnique(),0,&normdiffp_);
  if(acouopt_)
  {
    double acoucontrib = 0.0;
    INVANA::MVNorm(*acou_step_,*acou_matman_->ParamLayoutMapUnique(),0,&acoucontrib);
    normdiffp_ += acoucontrib;
  }

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
    std::cout<<"*** relative objective function value:    "<<J_/J_start_<<std::endl;
    std::cout<<"*** error value:                          "<<error_<<std::endl;
    std::cout<<"*** relative error value:                 "<<error_/error_start_<<std::endl;
    if(adaptreguval_!=0.0)
      std::cout<<"*** regularization weight                 "<<reguweight_<<std::endl;
    std::cout<<"*** norm of the scatra gradient:          "<<scatra_normgrad_<<std::endl;
    if(acouopt_)
      std::cout<<"*** norm of the acoustic gradient:        "<<acou_normgrad_<<std::endl;
    std::cout<<"*** norm of the difference of parameters: "<<normdiffp_<<std::endl;
    std::cout<<"*** output count:                         "<<output_count_<<std::endl;
    std::cout<<"*** simulation time since start [h]:      "<<(Teuchos::Time::wallTime()-tstart_)/(60.0*60.0)<<std::endl;
    std::cout<<"*** parameters:                           "<<std::endl;
  }
  if(scatra_matman_->GetParams()->GlobalLength()<30)
  {
    scatra_matman_->GetParams()->Print(std::cout);
    if(acouopt_)
      acou_matman_->GetParams()->Print(std::cout);
  }

  return;
} // void ACOU::InvAnalysis::OutputStats()

/*----------------------------------------------------------------------*/
bool ACOU::InvAnalysis::ConvergenceCheck()
{
  // check convergence criteria
  bool unconverged;

  if(calcacougrad_)
    unconverged = error_>tol_ && iter_<max_iter_ && scatra_normgrad_>tol_grad_;
  else
    unconverged = error_>tol_ && iter_<max_iter_;

  if(std::abs(normdiffp_)<1.0e-14)
    unconverged = false;

  return unconverged;
} // bool ACOU::InvAnalysis::ConvergenceCheck()

/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_MultiVector> ACOU::InvAnalysis::ElementMatVec()
{
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
  }
  return reacvec;
} // const Teuchos::RCP<Epetra_Vector> ACOU::InvAnalysis::ElementMatVec()

/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()
{
  return Teuchos::rcp(new AcouInvResultTest(*this));
} // Teuchos::RCP<DRT::ResultTest> ACOU::InvAnalysis::CreateFieldTest()

/*----------------------------------------------------------------------*/
void ACOU::InvAnalysis::EstimateMua(int count)
{
  // backprojection-algorithm

  // create unique map
  std::vector<int> pressuremonmicsunique;
  acou_discret_->Comm().Barrier();
  for(int i=0; i<acou_discret_->Comm().NumProc(); ++i)
  {
    if(acou_discret_->Comm().MyPID() == i)
    {
      for(unsigned int j=0; j<nodes_.size(); ++j)
      {
        if(acou_discret_->HaveGlobalNode(nodes_[j]))
        {
          if(acou_discret_->gNode(nodes_[j])->Owner()==int(i))
            pressuremonmicsunique.push_back(nodes_[j]);
        }
      }
    }
    acou_discret_->Comm().Barrier();
  }
  acou_discret_->Comm().Barrier();

  //*************************************************************************
  // Calculate middle point of monitor line

  // load number of nods in the circle around the measurement geometry
  int N=pressuremonmicsunique.size();
  double middlepoint[3];
  double nod_coords[N][3];
  for(int nd=0;nd<N;++nd)
  {
    const double* nod_coords_help = acou_discret_->gNode(pressuremonmicsunique[nd])->X();
    nod_coords[nd][0]=nod_coords_help[0];
    nod_coords[nd][1]=nod_coords_help[1];
    nod_coords[nd][2]=nod_coords_help[2];
  }


  LINALG::FixedSizeSerialDenseSolver<3,3,1> inverseCoeff;
  LINALG::Matrix<3,1> RHS(true);
  LINALG::Matrix<3,3> LHS(true);

  if(N>3)
  {
    for(unsigned int i=0;i<3;++i)
    {
      int m= int(N/3)*(i+1)-1;
      LHS(i,0)=1;
      LHS(i,1)= -nod_coords[m][0];
      LHS(i,2)= -nod_coords[m][1];
      RHS(i,0)=-(nod_coords[m][0]*nod_coords[m][0]+nod_coords[m][1]*nod_coords[m][1]);
    }
    inverseCoeff.SetMatrix(LHS);
    inverseCoeff.SetVectors(RHS,RHS);
    int err = inverseCoeff.Solve();
    if(err != 0) dserror("Inversion of matrix in light evaluation failed with error-code %d",err);
    //std::cout<<"Der Mittelpunkt liegt bei x= "<<RHS(1,0)/2<<" und bei y= "<<RHS(2,0)/2<<std::endl;

    middlepoint[0]=RHS(1,0)/2;
    middlepoint[1]=RHS(2,0)/2;
    middlepoint[2]=0.0;
  }

  double gmiddlepoint[3];
  int procnum =acou_discret_->Comm().NumProc();
  scatra_discret_->Comm().SumAll(&middlepoint[0],&gmiddlepoint[0],1);
  scatra_discret_->Comm().SumAll(&middlepoint[1],&gmiddlepoint[1],1);
  gmiddlepoint[0]=gmiddlepoint[0]/procnum;
  gmiddlepoint[1]=gmiddlepoint[1]/procnum;
  gmiddlepoint[2]=0;
  for(int check=0;check<1;check++)
  {
    if(fabs(gmiddlepoint[check]-middlepoint[check])>(1e-10))
      dserror("Error: Your monitor-line-geometry is not a circle, backprojection only works with circles in 2D");
  }

  //***********************************************************************

  // occupied area of a detector
  double deltaS=1.0;
  // speed of sound in water
  double c=1484.0;
  // Grueneisen-parameter
  double gamma=1.0;
  double phi[scatra_discret_->NumMyRowElements()];

  //***********************************************************************
  // phi-correction:

  if(count)
  {
    for(int a=0;a<scatra_discret_->NumMyRowElements();a++)
    {
      DRT::Element* myElement=scatra_discret_->lRowElement(a);
      std::vector<int> eledofs(myElement->NumNode());
      DRT::Node** nodes = myElement->Nodes();

      for(int nd=0; nd<myElement->NumNode();++nd)
        eledofs[nd]=scatra_discret_->DofRowMap()->LID(scatra_discret_->Dof(nodes[nd])[0]);

      double phiK[myElement->NumNode()];
      int c=0;
      for(unsigned int i=0;i<eledofs.size();i++)
      {
        if(eledofs[i]<0)
        {
          c++;
          phiK[i]=0;
        }
        if(eledofs[i]>=0)
          phiK[i]=phi_->operator[](eledofs[i]);
      }
      double help=0.0;
      for(int j=0;j<myElement->NumNode();j++)
        help=help+phiK[j];

      phi[a]=help/(myElement->NumNode()-c);
    }
  }
  else
  {
    for(int j=0;j<scatra_discret_->NumMyRowElements();j++)
      phi[j]=1.0;
  }
  //******************************************************************

  double r_e[3]={0.0,0.0,0.0};
  int mingid=scatra_discret_->ElementRowMap()->MinAllGID();
  for(int e=0; e<scatra_discret_->NumGlobalElements(); ++e)
  {
    DRT::Element* optele = NULL;
    int myopteleowner = -1;
    int opteleowner = -1;
    if(scatra_discret_->HaveGlobalElement(e+mingid))
    {
      double x_sum=0.0;
      double y_sum=0.0;
      double z_sum=0.0;
      myopteleowner = scatra_discret_->Comm().MyPID();
      optele = scatra_discret_->gElement(e+mingid);
      //load r_e (element middle-pointing vector)
      //load the nod IDs of the analyzed element
      const int* nodeids = optele->NodeIds();
      int numnode = optele->NumNode();

      //load the coordinates of every nod for the current element
      for(int loc_pos=0;loc_pos<numnode;loc_pos++)
      {
        const double* nodcoords = scatra_discret_->gNode(nodeids[loc_pos])->X();
        x_sum=x_sum+nodcoords[0];
        y_sum=y_sum+nodcoords[1];
        z_sum=z_sum+nodcoords[2];
      }

      //create the element's middle-pointing vector by building an average of coord's sums
      r_e[0]=x_sum/numnode;
      r_e[1]=y_sum/numnode;
      r_e[2]=z_sum/numnode;
    }
    scatra_discret_->Comm().MaxAll(&myopteleowner,&opteleowner,1);
    scatra_discret_->Comm().Broadcast(r_e,3,opteleowner);


    // every proc calculates its contribution to the sums
    //----------------------------------------------------------------------------------------------
    // build normal on measurement surface
    double normal[N][3];
    for(int i=0;i<N;++i)
    {
      double f=nod_coords[i][0]-gmiddlepoint[0];
      double s=nod_coords[i][1]-gmiddlepoint[1];
      double t=nod_coords[i][2]-gmiddlepoint[2];

      double a=1/sqrt(f*f+s*s+t*t);
      normal[i][0]=a*f;
      normal[i][1]=a*s;
      normal[i][2]=a*t;
    }

    //----------------------------------------------------------------------------------------------
    // build norm and quadratic norm of (d_i-r)
    double norm[N];
    double quadnorm[N];
    for(int i=0;i<N;i++)
      norm[i]=sqrt((nod_coords[i][0]-r_e[0])*(nod_coords[i][0]-r_e[0])+(nod_coords[i][1]-r_e[1])*(nod_coords[i][1]-r_e[1])+(nod_coords[i][2]-r_e[2])*(nod_coords[i][2]-r_e[2]));

    for(int i=0;i<N;i++)
      quadnorm[i]=(nod_coords[i][0]-r_e[0])*(nod_coords[i][0]-r_e[0])+(nod_coords[i][1]-r_e[1])*(nod_coords[i][1]-r_e[1])+(nod_coords[i][2]-r_e[2])*(nod_coords[i][2]-r_e[2]);

    // build scalar product of the normal and (r_e-d_i)
    double scalar[N];
    double difference[N][3];
    for(int i=0;i<N;i++)
    {
      difference[i][0]=(-nod_coords[i][0]+r_e[0]);
      difference[i][1]=(-nod_coords[i][1]+r_e[1]);
      difference[i][2]=(-nod_coords[i][2]+r_e[2]);
    }
    for(int i=0;i<N;i++)
      scalar[i]=difference[i][0]*normal[i][0]+difference[i][1]*normal[i][1]+difference[i][2]*normal[i][2];

    //build omega_i
    double omega[N];
    for(int i=0;i<N;i++)
      omega[i]=(deltaS/quadnorm[i])*(scalar[i]/norm[i]);

    //build b_i
    double b[N];
    double p[N];
    double tol;
    unsigned int ts[N];
    double calc_ts[N];
    double derivp[N];

    for(int i=0;i<N;i++)
      calc_ts[i]=norm[i]/c;
    int z=timesteps_.size();

    // search for the corresponding time-step and pressure in the measured curve
    for(int i=0;i<N;i++)
    {
      double lastone=1.0;
      tol=timesteps_[1]-timesteps_[0];
      for(unsigned int j=0;j<timesteps_.size();j++)
      {
        if((abs(timesteps_[j]-calc_ts[i])<tol)&&(abs(timesteps_[j]-calc_ts[i])<lastone))
        {
          p[i]=nod_curve_interpol_[z*i+j];
          ts[i]=j;
          lastone=abs(timesteps_[j]-calc_ts[i]);
        }
      }
    }

    // calculate the derivative of p via a finite difference approach
    for(int i=0;i<N;i++)
      derivp[i]=(-1/60*nod_curve_interpol_[(ts[i]-3)+z*i]+3/20*nod_curve_interpol_[(ts[i]-2)+z*i]-3/4*nod_curve_interpol_[(ts[i]-1)+z*i]+1/60*nod_curve_interpol_[(ts[i]+3)+z*i]-3/20*nod_curve_interpol_[(ts[i]+2)+z*i]+3/4*nod_curve_interpol_[(ts[i]+1)+z*i])/(c*(timesteps_[ts[i]]-timesteps_[ts[i]-1]));

    // calculate b_i as difference of 2*p-2*t_queer*derivative(p)
    for(int i=0;i<N;i++)
      b[i]=2*p[i]-2*norm[i]*derivp[i];

    //build sums
    double enumerator_sum=0.0;
    double denominator_sum=0.0;
    for(int i=0;i<N;i++)
    {
      enumerator_sum=enumerator_sum+(omega[i]*b[i]);
      denominator_sum=denominator_sum+omega[i];
    }

    double genumerator_sum=0.0;
    double gdenominator_sum=0.0;

    // sum over all processors
    scatra_discret_->Comm().SumAll(&enumerator_sum,&genumerator_sum,1);
    scatra_discret_->Comm().SumAll(&denominator_sum,&gdenominator_sum,1);

    double mu_estimated = 0.0;
    if(scatra_discret_->HaveGlobalElement(e+mingid))
    {
      int a=scatra_discret_->ElementRowMap()->LID(e+mingid);
      if(a>=0)
      {
        mu_estimated=(-1.0/(gamma*phi[a]))*(genumerator_sum/gdenominator_sum);
        scatra_step_->operator ()(0)->operator [](a) = mu_estimated;
      }
    }
  }
  scatra_matman_->ReplaceParams(*scatra_step_);

  return;
}

/*----------------------------------------------------------------------*/
ACOU::PATLineSearch::PATLineSearch(Teuchos::RCP<InvAnalysis> problem, bool successive)
{
  invproblem_ = problem;
  itermax_ = problem->max_ls_iter_;
  acouopt_ = problem->acouopt_;
  alpha_max_ = 15.0;
  normacou_ = 0.0;

  c1_ =  1.0e-12;
  c2_ = 0.9;

  // TODO: validate
  successive_ = successive;
  segmentation_ = problem->segmentation_;

  if(segmentation_)
    k_o_ = Teuchos::rcp(new Epetra_MultiVector(problem->k_->Map(),problem->k_->NumVectors()));
}

/*----------------------------------------------------------------------*/
void ACOU::PATLineSearch::Init()
{
  // initial objective function value
  J_0_ = invproblem_->J_;

  // calculation of the initial  norm for scatra
  INVANA::MVDotProduct(*invproblem_->scatra_objgrad_,*invproblem_->scatra_d_,*invproblem_->scatra_matman_->ParamLayoutMapUnique(),&normsca_0_);

  // calculation of the initial  norm for acou
  if(acouopt_) INVANA::MVDotProduct(*invproblem_->acou_objgrad_,*invproblem_->acou_d_,*invproblem_->acou_matman_->ParamLayoutMapUnique(),&normacou_0_);
  else         normacou_0_ = 0.0;

  // set initial trial step length
  alpha_i_ = 1.0;
  alpha_im1_ = 0.0;

  if(segmentation_)
  {
    k_o_->Update(1.0,*invproblem_->k_,0.0);
    INVANA::MVDotProduct(*invproblem_->k_objgrad_,*invproblem_->k_d_,*invproblem_->scatra_discret_->ElementColMap(),&normk_0_);
    normk_0_ *= 3.0;
  }

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::Run()
{
  // Do it!
  bool success = false;

  if(segmentation_)                                  success = RunSegmentation();
  else if(successive_)                               success = RunSuccessive();
  else                                               success = RunSimultaneous();

  return success;
}


/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::RunSuccessive()
{
  // Do it!

  // perform the optical optimization, break if succeeded
  for(int i=0; i<itermax_; ++i)
  {
    if(!(invproblem_->myrank_))
      std::cout<<"*************** optical line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_ "<<invproblem_->J_<<", normsca_0_ "<<normsca_0_<<", normacou_0_ "<<normacou_0_<<std::endl;

    // save previous objective function value
    J_im1_ = invproblem_->J_;

    // update the optical parameters
    if(i>0) invproblem_->scatra_matman_->ResetParams();
    invproblem_->scatra_step_->Update(alpha_i_,*invproblem_->scatra_d_,0.0);
    invproblem_->scatra_matman_->UpdateParams(invproblem_->scatra_step_);

    // do the forward problem to be able to evaluate the objective function
    invproblem_->SolveStandardProblemScatra();
    invproblem_->SolveStandardProblemAcou();

    // evaluate the objective function
    invproblem_->CalculateObjectiveFunctionValue(i==0);

    if( invproblem_->J_> J_0_ + c1_ * alpha_i_ * normsca_0_ || (invproblem_->J_ >= J_im1_ && i>0) )
    {
      alpha_x_ = Zoom(alpha_im1_, alpha_i_, J_0_, J_im1_, normsca_0_, false, true, false);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 1 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 1 not met"<<std::endl;

    // evaluate phi'(alpha_i)
    //{
    // Solve the adjoint problem
    invproblem_->SolveAdjointAcouProblem();
    invproblem_->SolveAdjointOptiProblem();

    // Build the gradient
    invproblem_->CalculateGradient();

    INVANA::MVDotProduct(*invproblem_->scatra_objgrad_,*invproblem_->scatra_d_,*invproblem_->scatra_matman_->ParamLayoutMapUnique(),&normsca_);
    //}

    if(std::abs(normsca_) <= -c2_*normsca_0_)
    {
      alpha_x_ = alpha_i_;
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 2 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 2 not met"<<std::endl;

    if(normsca_ >= 0)
    {
      alpha_x_ = Zoom(alpha_i_,alpha_im1_, J_0_, invproblem_->J_, normsca_, true, true, false);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 3 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 3 not met"<<std::endl;

    alpha_im1_ = alpha_i_;
    alpha_i_ *= 2.0;
    //alpha_i_ = (alpha_i_+alpha_max_)/2.0;
  }

  bool scasucc;
  if(alpha_x_ != 0.0)
  {
    scasucc = true;
    if(!invproblem_->myrank_)
      std::cout<<"*************** optical line search succeeded"<<std::endl;
  }
  else
  {
    scasucc = false;
    if(!invproblem_->myrank_)
      std::cout<<"*************** optical line search failed"<<std::endl;
    // return false; // not necessart, it is enough if one of them succeeds
  }

  alpha_im1_ = 0.0;
  alpha_i_ = 1.0;
  J_0_ = invproblem_->J_;

  if(acouopt_)
  {
    // use the new gradient for the acoustical update!
    //{
    if(invproblem_->opti_ != INPAR::ACOU::inv_gd)
      dserror("not yet implemented for BFGS");
    invproblem_->acou_d_->Update(-1.0,*invproblem_->acou_objgrad_,0.0);
    invproblem_->acou_d_->Scale(1.0/invproblem_->acou_normgrad0_*double(invproblem_->acou_objgrad_->Map().NumGlobalElements()));
    invproblem_->acou_objgrad_o_->Update(1.0,*invproblem_->acou_objgrad_,0.0);
    //}

    // perform the  acoustical optimization, break if succeeded
    for(int i=0; i<itermax_; ++i)
    {
      if(!(invproblem_->myrank_))
        std::cout<<"*************** acoustical line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_ "<<invproblem_->J_<<", normsca_0_ "<<normsca_0_<<", normacou_0_ "<<normacou_0_<<std::endl;

      // save previous objective function value
      J_im1_ = invproblem_->J_;

      // update the acoustical parameters
      if(i>0) invproblem_->acou_matman_->ResetParams();
      invproblem_->acou_step_->Update(alpha_i_,*invproblem_->acou_d_,0.0);
      invproblem_->acou_matman_->UpdateParams(invproblem_->acou_step_);

      // do the forward problem to be able to evaluate the objective function
      invproblem_->SolveStandardProblemScatra();
      invproblem_->SolveStandardProblemAcou();

      // evaluate the objective function
      invproblem_->CalculateObjectiveFunctionValue(i==0);

      if( invproblem_->J_> J_0_ + c1_ * alpha_i_ * normacou_0_ || (invproblem_->J_ >= J_im1_ && i>0) )
      {
        alpha_x_ = Zoom(alpha_im1_, alpha_i_, J_0_, J_im1_, normacou_0_, false, false, true);
        if(!invproblem_->myrank_)
          std::cout<<"*************** condition 1 met"<<std::endl;
        break;
      }
      else if(!invproblem_->myrank_)
        std::cout<<"*************** condition 1 not met"<<std::endl;

      // evaluate phi'(alpha_i)
      //{
      // Solve the adjoint problem
      invproblem_->SolveAdjointAcouProblem();
      invproblem_->SolveAdjointOptiProblem();

      // Build the gradient
      invproblem_->CalculateGradient();

      INVANA::MVDotProduct(*invproblem_->acou_objgrad_,*invproblem_->acou_d_,*invproblem_->acou_matman_->ParamLayoutMapUnique(),&normacou_);
      //}

      if(std::abs(normacou_) <= -c2_*normacou_0_)
      {
        alpha_x_ = alpha_i_;
        if(!invproblem_->myrank_)
          std::cout<<"*************** condition 2 met"<<std::endl;
        break;
      }
      else if(!invproblem_->myrank_)
        std::cout<<"*************** condition 2 not met: normacou_ "<<normacou_<<", normacou_0_ "<<normacou_0_<<std::endl;

      if(normacou_ >= 0)
      {
        alpha_x_ = Zoom(alpha_i_,alpha_im1_, J_0_, invproblem_->J_, normacou_0_, true, false, true);
        if(!invproblem_->myrank_)
          std::cout<<"*************** condition 3 met"<<std::endl;
        break;
      }
      else if(!invproblem_->myrank_)
        std::cout<<"*************** condition 3 not met"<<std::endl;

      alpha_im1_ = alpha_i_;
      alpha_i_ *= 2.0;
      //alpha_i_ = (alpha_i_+alpha_max_)/2.0;
    }

    if(alpha_x_ != 0.0)
    {
      if(!invproblem_->myrank_)
        std::cout<<"*************** acoustical line search succeeded"<<std::endl;
      return true;
    }
    else
    {
      if(!invproblem_->myrank_)
        std::cout<<"*************** acoustical line search failed"<<std::endl;
      if(scasucc == false)
        return false;
      else
        return true; // is is enough when one of them succeeds
    }
  }
  // in case you come here, you failed (otherwise return earlier) or you only had scatra parameter
  if(acouopt_==false)
    return scasucc;
  else
    return false;
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::RunSimultaneous()
{
  // Do it!

  // perform the optimization, break if succeeded
  for(int i=0; i<itermax_; ++i)
  {
    if(!(invproblem_->myrank_))
      std::cout<<"*************** simultaneous line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_ "<<invproblem_->J_<<", normsca_0_ "<<normsca_0_<<", normacou_0_ "<<normacou_0_<<std::endl;

    // save previous objective function value
    J_im1_ = invproblem_->J_;

    // update the optical parameters
    if(i>0) invproblem_->scatra_matman_->ResetParams();
    invproblem_->scatra_step_->Update(alpha_i_,*invproblem_->scatra_d_,0.0);
    invproblem_->scatra_matman_->UpdateParams(invproblem_->scatra_step_);

    if(acouopt_)
    {
      // update the acoustical parameters
      if(i>0) invproblem_->acou_matman_->ResetParams();
      invproblem_->acou_step_->Update(alpha_i_,*invproblem_->acou_d_,0.0);
      invproblem_->acou_matman_->UpdateParams(invproblem_->acou_step_);
    }
    // do the forward problem to be able to evaluate the objective function
    invproblem_->SolveStandardProblemScatra();
    invproblem_->SolveStandardProblemAcou();

    // evaluate the objective function
    invproblem_->CalculateObjectiveFunctionValue(i==0);

    if( invproblem_->J_> J_0_ + c1_ * alpha_i_ * (normsca_0_+normacou_0_) || (invproblem_->J_ >= J_im1_ && i>0) )
    {
      alpha_x_ = Zoom(alpha_im1_, alpha_i_, J_0_, J_im1_, (normsca_0_+normacou_0_), false, true, true);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 1 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 1 not met"<<std::endl;

    // evaluate phi'(alpha_i)
    //{
    // Solve the adjoint problem
    invproblem_->SolveAdjointAcouProblem();
    invproblem_->SolveAdjointOptiProblem();

    // Build the gradient
    invproblem_->CalculateGradient();

    INVANA::MVDotProduct(*invproblem_->scatra_objgrad_,*invproblem_->scatra_d_,*invproblem_->scatra_matman_->ParamLayoutMapUnique(),&normsca_);
    if(acouopt_) INVANA::MVDotProduct(*invproblem_->acou_objgrad_,*invproblem_->acou_d_,*invproblem_->acou_matman_->ParamLayoutMapUnique(),&normacou_);
    //}

    if(std::abs(normsca_+normacou_) <= -c2_*(normsca_0_+normacou_0_))
    {
      alpha_x_ = alpha_i_;
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 2 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 2 not met: norm_ "<<normacou_+normsca_<<", norm_0_ "<<normacou_0_+normsca_0_<<std::endl;

    if(normsca_+normacou_ >= 0)
    {
      alpha_x_ = Zoom(alpha_i_,alpha_im1_, J_0_, invproblem_->J_, normsca_0_+normacou_0_, true, true, true);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 3 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 3 not met: norm_ "<<normacou_+normsca_<<", norm_0_ "<<normacou_0_+normsca_0_<<std::endl;

    alpha_im1_ = alpha_i_;
    alpha_i_ *= 2.0;
    //alpha_i_ = (alpha_i_+alpha_max_)/2.0;
  }

  if(alpha_x_ != 0.0)
  {
    if(!invproblem_->myrank_)
      std::cout<<"*************** simultaneous line search succeeded"<<std::endl;
    return true;
  }
  else
  {
    if(!invproblem_->myrank_)
      std::cout<<"*************** simultaneous line search failed"<<std::endl;
    return false;
  }
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::RunSegmentation()
{
  // Do it!
  //invproblem_->f_d_->Print(std::cout);

  std::cout<<"RUNSEGMENTATION "<<J_0_<<" "<<normk_0_<<std::endl;

  // perform the optimization, break if succeeded
  for(int i=0; i<itermax_; ++i)
  {
    if(!(invproblem_->myrank_))
      std::cout<<"*************** segmentation line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_ "<<invproblem_->J_<<", normk_0_ "<<normk_0_<<std::endl;

    // save previous objective function value
    J_im1_ = invproblem_->J_;

    // update parameters
    invproblem_->k_->Update(1.0,*k_o_,0.0);
    invproblem_->k_->Update(alpha_i_,*invproblem_->k_d_,1.0);
    invproblem_->segmentationsolver_->UpdateParameters(invproblem_->k_,invproblem_->scatra_matman_,invproblem_->acou_matman_,invproblem_->scatra_step_,invproblem_->acou_step_);
    invproblem_->scatra_matman_->ReplaceParams(*invproblem_->scatra_step_);
    invproblem_->acou_matman_->ReplaceParams(*invproblem_->acou_step_);

    // do the forward problem to be able to evaluate the objective function
    invproblem_->SolveStandardProblemScatra();
    invproblem_->SolveStandardProblemAcou();

    // evaluate the objective function
    invproblem_->CalculateObjectiveFunctionValue(i==0);

    // check first condition (sufficient decrease)
    if( invproblem_->J_> J_0_ + c1_ * alpha_i_ * normk_0_ || (invproblem_->J_ >= J_im1_ && i>0) )
    {
      alpha_x_ = Zoom(alpha_im1_, alpha_i_, J_0_, J_im1_, normk_0_, false, true, true);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 1 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 1 not met"<<std::endl;

    // evaluate phi'(alpha_i)
    //{
    // Solve the adjoint problem
    invproblem_->SolveAdjointAcouProblem();
    invproblem_->SolveAdjointOptiProblem();

    // Build the gradient
    invproblem_->CalculateGradient();

    INVANA::MVDotProduct(*invproblem_->k_objgrad_,*invproblem_->k_d_,*invproblem_->scatra_discret_->ElementColMap(),&normk_);
    //}

    if(std::abs(normk_) <= -c2_*normk_0_)
    {
      alpha_x_ = alpha_i_;
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 2 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 2 not met: norm_ "<<normk_<<", norm_0_ "<<normk_0_<<std::endl;

    if(normk_ >= 0)
    {
      alpha_x_ = Zoom(alpha_i_,alpha_im1_, J_0_, invproblem_->J_, normk_0_, true, true, true);
      if(!invproblem_->myrank_)
        std::cout<<"*************** condition 3 met"<<std::endl;
      break;
    }
    else if(!invproblem_->myrank_)
      std::cout<<"*************** condition 3 not met: norm_ "<<normk_<<", norm_0_ "<<normk_0_<<std::endl;

    alpha_im1_ = alpha_i_;
    alpha_i_ *= 2.0;
    //alpha_i_ = (alpha_i_+alpha_max_)/2.0;
  }

  if(alpha_x_ != 0.0)
  {
    if(!invproblem_->myrank_)
      std::cout<<"*************** segmentation line search succeeded"<<std::endl;
    return true;
  }
  else
  {
    if(!invproblem_->myrank_)
      std::cout<<"*************** segmentation line search failed"<<std::endl;
    return false;
  }

}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::Zoom(double alpha_lo, double alpha_hi, double J_0, double J_alpha_lo, double normphi_0, bool turn, bool sca, bool acou)
{
  double normphi_j = 0.0;
  double alpha_j = 0.0;

  for(int j=0; j<itermax_; ++j)
  {
    if(turn)
      alpha_j = alpha_hi+(alpha_lo-alpha_hi)/3.0;
    else
    {
      alpha_j = (alpha_lo + alpha_hi) / 2.0;
      //if(J_0_/invproblem_->J_>1.0e-3)
      //  alpha_j = alpha_lo + (alpha_hi-alpha_lo)*J_0/invproblem_->J_/2.0;
      //else
      //  alpha_j = alpha_lo+(alpha_hi-alpha_lo)/10.0;
    }
    if(!invproblem_->myrank_)
      std::cout<<"*************** zoom "<<j<<": alpha_lo "<<alpha_lo<<" alpha_hi "<<alpha_hi<<" alpha_j "<<alpha_j<<std::endl;

    if(segmentation_)
    {
      invproblem_->k_->Update(1.0,*k_o_,0.0);
      invproblem_->k_->Update(alpha_j,*invproblem_->k_d_,1.0);

      invproblem_->segmentationsolver_->UpdateParameters(invproblem_->k_,invproblem_->scatra_matman_,invproblem_->acou_matman_,invproblem_->scatra_step_,invproblem_->acou_step_);
      invproblem_->scatra_matman_->ReplaceParams(*invproblem_->scatra_step_);
      invproblem_->acou_matman_->ReplaceParams(*invproblem_->acou_step_);
    }
    else
    {
      if(sca)
      {
        invproblem_->scatra_matman_->ResetParams();
        invproblem_->scatra_step_->Update(alpha_j,*invproblem_->scatra_d_,0.0);
        invproblem_->scatra_matman_->UpdateParams(invproblem_->scatra_step_);
      }
      if(acou && acouopt_)
      {
        invproblem_->acou_matman_->ResetParams();
        invproblem_->acou_step_->Update(alpha_j,*invproblem_->acou_d_,0.0);
        invproblem_->acou_matman_->UpdateParams(invproblem_->acou_step_);
      }
    }
    // do the forward problem to be able to evaluate the objective function
    invproblem_->SolveStandardProblemScatra();
    invproblem_->SolveStandardProblemAcou();

    // evaluate the objective function
    invproblem_->CalculateObjectiveFunctionValue(false);

    std::cout<<"bedingung: J_ "<<invproblem_->J_<<" J_0 "<<J_0<<" normphi_0 "<<normphi_0<<" J_alpha_lo "<<J_alpha_lo<<std::endl;
    if( invproblem_->J_ > J_0 + c1_ * alpha_j * normphi_0 || invproblem_->J_ >= J_alpha_lo )
      alpha_hi = alpha_j;
    else
    {
      // evaluate phi'(alpha_j)
      //{
      // Solve the adjoint problem
      invproblem_->SolveAdjointAcouProblem();
      invproblem_->SolveAdjointOptiProblem();

      // Build the gradient
      invproblem_->CalculateGradient();

      INVANA::MVDotProduct(*invproblem_->scatra_objgrad_,*invproblem_->scatra_d_,*invproblem_->scatra_matman_->ParamLayoutMapUnique(),&normsca_);
      if(acouopt_ && acou)
        INVANA::MVDotProduct(*invproblem_->acou_objgrad_,*invproblem_->acou_d_,*invproblem_->acou_matman_->ParamLayoutMapUnique(),&normacou_);
      if(segmentation_)
        INVANA::MVDotProduct(*invproblem_->k_objgrad_,*invproblem_->k_d_,*invproblem_->scatra_discret_->ElementColMap(),&normk_);

      //}
      if(segmentation_)
        normphi_j = normk_;
      else if(sca==true && acou==false)
        normphi_j = normsca_;
      else if(sca==false && acou==true)
        normphi_j = normacou_;
      else if(sca==true && acou==true)
        normphi_j = normsca_ + normacou_;

      std::cout<<"norm "<<normphi_j<<" "<<normphi_0<<std::endl;
      if(std::abs(normphi_j) <= -c2_*normphi_0)
          return alpha_j;

      if( normphi_j * (alpha_hi - alpha_lo) >= 0 )
        alpha_hi = alpha_lo;

      alpha_lo = alpha_j;
    }
  }

  return 0.0;
}

ACOU::InvAnalysis::Segmentation::Segmentation(std::string materialfilename, double penaltypara, Teuchos::RCP<DRT::Discretization> scatradiscretization, Teuchos::RCP<DRT::DiscretizationHDG> acoudiscretization)
{
  penalty = penaltypara;
  scadisc = scatradiscretization;
  acoudisc = acoudiscretization;

  // read from the given file
  if (materialfilename=="none.material") dserror("No material file provided");

  // insert path to monitor file if necessary
  if (materialfilename[0]!='/')
  {
    std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::string::size_type pos = filename.rfind('/');
    if (pos!=std::string::npos)
    {
      std::string path = filename.substr(0,pos+1);
      materialfilename.insert(materialfilename.begin(), path.begin(), path.end());
    }
  }

  // open the file
  FILE* file = fopen(materialfilename.c_str(),"rb");
  if (file==NULL) dserror("Could not open material file %s",materialfilename.c_str());

  // prepare read in quantities
  char buffer[150000];
  fgets(buffer,150000,file);
  char* foundit = NULL;

  // read number of materials
  foundit = strstr(buffer,"nummats");
  foundit += strlen("nummats");
  nummats = strtol(foundit,&foundit,10);

  // prepare the materials
  materialtable.resize(nummats);
  for(unsigned i=0; i<nummats; ++i)
    materialtable[i].resize(4); // mu_a, D, c, rho

  // read the materials
  foundit = buffer;
  fgets(buffer,150000,file);
  for(unsigned int i=0; i<nummats; ++i)
  {
    for(int j=0; j<4; ++j)
    {
      materialtable[i][j] = strtod(foundit,&foundit);
    }
    fgets(buffer,150000,file);
    foundit = buffer;
  }
  /*
  {
    std::cout<<std::endl;
    std::cout<<"You specified "<<nummats<<" materials for segmentation with the following values:"<<std::endl;
    for(unsigned int i=0; i<nummats; ++i)
        std::cout<<"mu_a "<<materialtable[i][0]<<" D "<<materialtable[i][1]<<" c "<<materialtable[i][2]<<" rho "<<materialtable[i][3]<<std::endl;
  }
  */
}

void ACOU::InvAnalysis::Segmentation::InitParameters(
    Teuchos::RCP<Epetra_MultiVector> concentrations,
    Teuchos::RCP<INVANA::MatParManager> scatramatman,
    Teuchos::RCP<INVANA::MatParManager> acoumatman,
    Teuchos::RCP<Epetra_MultiVector> scatraparams,
    Teuchos::RCP<Epetra_MultiVector> acouparams)
{
  // all materials have the same initial concentration
  /*
  double c = tan(PI*(1./double(nummats+0.2)-0.5));
  for(unsigned int m=0; m<nummats; ++m)
    concentrations->operator ()(m)->PutScalar(c);
  */

  // the first material is assumed to be the default material and thus has highest concentration (default concentration defcon 0.9)
  double defcon = 0.9;
  double c = tan(PI*(defcon-0.5));
  concentrations->operator ()(0)->PutScalar(c);

  // the others make up the rest with equal parts
  c = tan(PI*((1.0-defcon)/(double(nummats)-1.)-0.5));
  for(unsigned int m=1; m<nummats; ++m)
    concentrations->operator ()(m)->PutScalar(c);


  // now, bring the concentration distribution to the parameter vectors
  UpdateParameters(concentrations,scatramatman,acoumatman,scatraparams,acouparams);

  // bring the parameters to the materials
  scatramatman->ReplaceParams(*scatraparams);
  acoumatman->ReplaceParams(*acouparams);

  return;
}

void ACOU::InvAnalysis::Segmentation::CalculateGradient(
    Teuchos::RCP<Epetra_MultiVector> sca_params,
    Teuchos::RCP<Epetra_MultiVector> acou_params,
    Teuchos::RCP<Epetra_MultiVector> sca_grad,
    Teuchos::RCP<Epetra_MultiVector> acou_grad,
    Teuchos::RCP<Epetra_MultiVector> k_params,
    Teuchos::RCP<Epetra_MultiVector> k_grad)
{
  // loop all entries
  // we make use of the knowledge, that there are two parameters per material:
  int numeleentries = sca_params->MyLength() / 2;
  double k_grad_val = 0.0;
  for(int i=0; i<numeleentries; ++i)
  {
    // for every entry, the chain rule  has to be applied forward and backward
    // for example the direction of optimization for the absorption coefficient (here abbreviated with a) is
    // - (dJ/da da/dk + dJ/dD dD/dk + dJ/dc dc/dk + dJ/dr dr/dk) dk/da
    // first the gradient for k is calculated, then the update for the absorption coefficient!

    double equ = -1.0;
    // sum constraint violation
    for(unsigned int m=0; m<nummats; ++m)
      equ += 1./PI*atan(k_params->operator ()(m)->operator [](i))+0.5;

    for(unsigned int m=0; m<nummats; ++m)
    {
      // calculate dJ/dk = dJ/da da/dk + dJ/dD dD/dk + dJ/dc dc/dk + dJ/dr dr/dk
      k_grad_val = sca_grad->operator ()(0)->operator [](i) * materialtable[m][0]
                 + sca_grad->operator ()(0)->operator [](i+numeleentries) * materialtable[m][1]
                 + acou_grad->operator ()(0)->operator [](i) * materialtable[m][2]
                 + acou_grad->operator ()(0)->operator [](i+numeleentries) * materialtable[m][3];

      // contribution from penalty term
      k_grad_val += penalty * equ;

      // chain rule for meta parametrization
      double k = k_params->operator ()(m)->operator [](i);
      k_grad_val *= 1./PI/(k*k+1.);

      // write it in the k gradient vector
      k_grad->operator ()(m)->operator [](i) = k_grad_val;
    }

  }

  return;
}

void ACOU::InvAnalysis::Segmentation::UpdateParameters(Teuchos::RCP<Epetra_MultiVector> k_p,
                                                       Teuchos::RCP<INVANA::MatParManager> scatramatman,
                                                       Teuchos::RCP<INVANA::MatParManager> acoumatman,
                                                       Teuchos::RCP<Epetra_MultiVector> sca_p,
                                                       Teuchos::RCP<Epetra_MultiVector> acou_p)
{
  sca_p->PutScalar(0.0);
  acou_p->PutScalar(0.0);

  // loop the global scatra elements
  for(int i=scadisc->ElementRowMap()->MinAllGID(); i<=scadisc->ElementRowMap()->MaxAllGID(); ++i)
  {
    // lid of the element
    int lid = scadisc->ElementRowMap()->LID(i);

    // determine local material values
    double loc_mat[4] = {0.0};
    if(lid>=0)
    {
      // determine material value as sum of concentration and correspondent value
      for(unsigned int m=0; m<nummats; ++m)
      {
        double meta_k = k_p->operator ()(m)->operator [](lid);
        // calculate from meta to normal
        double k = 1./PI*atan(meta_k)+0.5;
        loc_mat[0] += k * materialtable[m][0];
        loc_mat[1] += k * materialtable[m][1];
        loc_mat[2] += k * materialtable[m][2];
        loc_mat[3] += k * materialtable[m][3];
      }
    }

    // communicate material values to all processors
    double glo_mat[4] = {0.0};
    scadisc->Comm().SumAll(&loc_mat[0],&glo_mat[0],4);

    // write scatra values to scatra vector
    if(lid>=0)
    {
      int elematid = scadisc->lRowElement(lid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatramatman)->WriteValuesToVector(elematid,i,glo_mat[1],glo_mat[0],sca_p);
    }

    // write acou values to acou vector
    int agid = i - scadisc->ElementRowMap()->MinAllGID() + acoudisc->ElementRowMap()->MinAllGID();
    int alid = acoudisc->ElementRowMap()->LID(agid);
    if(alid >=0)
    {
      int elematid = acoudisc->lRowElement(alid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::AcouMatParManagerPerElement>(acoumatman)->WriteValuesToVector(elematid,agid,glo_mat[3],glo_mat[2],acou_p);
    }
  }

  return;
}
