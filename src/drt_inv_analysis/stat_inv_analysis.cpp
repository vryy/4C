/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/


#include "stat_inv_analysis.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"

#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimint.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseMatrix.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../linalg/linalg_solver.H"

#include "../drt_structure/stru_resulttest.H"

// needed to deal with materials
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::StatInvAnalysis::StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    output_(output)
{
  // input parameters inverse analysis
  //const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  // Store displacement
  displacement_ = Teuchos::rcp(new Epetra_Vector(*discret_->DofRowMap(),true));
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::StatInvAnalysis::Sample()
{

  SetParameters();
  SolveForwardProblem();
  //double error = CalcErrorMonitorFile();
  //IO::cout << "And The error is " << error << " !" << IO::endl;
  IO::cout << "And the displacement is  " << *displacement_ << IO::endl;

  return;
}


void STR::StatInvAnalysis::SolveForwardProblem()
{

  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

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
      ADAPTER::StructureBaseAlgorithm adapterbase(DRT::Problem::Instance()->StructuralDynamicParams(), structdis);
      ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

      // do restart
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        structadaptor.ReadRestart(restart);
      }
      structadaptor.Integrate();

      DRT::Problem::Instance()->AddFieldTest(structadaptor.CreateFieldTest());
      DRT::Problem::Instance()->TestAll(structadaptor.DofRowMap()->Comm());

      // get displacement of discretization at the end
      displacement_= structadaptor.Dispn();
      Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(structadaptor.DofRowMap()->Comm());
      Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);

      // time to go home...
    }
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  return;
}


//---------------------------------------------------------------------------------------------
void STR::StatInvAnalysis::CreateNewSample()
{
  return;
}


//--------------------------------------------------------------------------------------
// Setup Material Parameters in each element based on parametrization
void STR::StatInvAnalysis::SetParameters()
{
  // Variables for Random field
  double mat_par;
  // element center
  std::vector<double> ele_c_location;

  // flag have init stochmat??
  int stochmat_flag=0;

  // Get parameters from stochastic matlaw
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if (myrank == 0)
    IO::cout << "No. material laws considered: " << (int)mats.size() << IO::endl;
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
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
      // get material
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->lColElement(i)->Material().get());
      // get location of element
      std::vector<double> ele_center;
      ele_center = discret_->lColElement(i)->ElementCenterRefeCoords();

      // Eval Parametrization at Elementlocation
      //TODO implement class that parametrizes spatial field of parameters
     // mat_par = parametrization_->EvalAtLocation(ele_center);
      // HACK set beta to 3.7 for now until we have a parametrizatin class
      mat_par= 3.7E6;
      aaa_stopro->Init(mat_par,"beta");
    }
  } // EOF loop elements

}
double STR::StatInvAnalysis::CalcErrorMonitorFile()
{
  // input parameters inverse analysis
   const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
   int myrank = discret_->Comm().MyPID();

  int ndofs_ = 0;
  int nsteps_ =0;
  int nnodes_ =0;
  std::vector<int>          nodes_;         // list of node gids observed
  std::vector<std::vector<int> > dofs_;          // list of dofs on each node that are monitored
  Epetra_SerialDenseVector mcurve_; // measured displacement of the experiments (target value)
  std::vector<double> timesteps_;   // time values where displacements are measured

    char* foundit = NULL;
    string filename = statinvp.get<string>("MONITORFILE");
    if (filename=="none.monitor") dserror("No monitor file provided");
    FILE* file = fopen(filename.c_str(),"rb");
    if (file==NULL) dserror("Could not open monitor file %s",filename.c_str());
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
    dofs_.resize(nnodes_);
    for (int i=0; i<nnodes_; ++i)
    {
      fgets(buffer,150000,file);
      foundit = buffer;
      nodes_[i] = strtol(foundit,&foundit,10);
      int ndofs = strtol(foundit,&foundit,10);
      ndofs_ += ndofs;
      dofs_[i].resize(ndofs,-1);
      if (!myrank) printf("Monitored node %d ndofs %d dofs ",nodes_[i],(int)dofs_[i].size());
      for (int j=0; j<ndofs; ++j)
      {
        dofs_[i][j] = strtol(foundit,&foundit,10);
        if (!myrank) printf("%d ",dofs_[i][j]);
      }
      if (!myrank) printf("\n");
    }

    // total number of measured dofs
    int nmp_    = ndofs_*nsteps_;


    // read in measured curve

    {
      mcurve_ = Epetra_SerialDenseVector(nmp_);

      if (!myrank) printf("nsteps %d ndofs %d\n",nsteps_,ndofs_);

      // read comment lines
      foundit = buffer;
      fgets(buffer,150000,file);
      while(strstr(buffer,"#"))
        fgets(buffer,150000,file);

      // read in the values for each node in dirs directions
      int count=0;
      for (int i=0; i<nsteps_; ++i)
      {
        // read the time step
        timesteps_[i] = strtod(foundit,&foundit);
        for (int j=0; j<ndofs_; ++j)
          mcurve_[count++] = strtod(foundit,&foundit);
        fgets(buffer,150000,file);
        foundit = buffer;
      }
      if (count != nmp_) dserror("Number of measured disps wrong on input");
    }

  // compute nonlinear problem and obtain computed displacements
  // output at the last step
  Epetra_SerialDenseVector cvector;
  cvector = SolveForwardProblemOldStyle(false,nmp_,timesteps_,ndofs_,nnodes_,dofs_,nodes_);
  Epetra_SerialDenseVector rcurve(nmp_);
  for (int i=0; i<nmp_; i++)
     rcurve[i] = mcurve_[i] - cvector[i];
  double error   = rcurve.Norm2()/sqrt(nmp_);
  return error;

}

Epetra_SerialDenseVector STR::StatInvAnalysis::SolveForwardProblemOldStyle(bool outputtofile, int nmp_ ,std::vector<double>  timesteps_ ,int ndofs_,int nnodes_,std::vector<std::vector<int> > dofs_,std::vector<int>   nodes_)
{
  // input parameters for structural dynamics
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  // create a solver
  // get the solver number used for structural solver
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                      discret_->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));

  discret_->ComputeNullSpaceIfNecessary(solver->Params());
    ///  time integrator for structural dynamics
    Teuchos::RCP<STR::TimInt> sti_ = Teuchos::null;
    int myrank = discret_->Comm().MyPID();
    // get input parameter lists
    const Teuchos::ParameterList& ioflags
      = DRT::Problem::Instance()->IOParams();

    Teuchos::ParameterList xparams;
    xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
    xparams.set<int>("REDUCED_OUTPUT",0);
    // create time integrator
    sti_ = TimIntCreate(ioflags, sdyn, xparams, discret_, solver, solver, output_);
    if (sti_ == Teuchos::null) dserror("Failed in creating integrator.");
    // initialize time loop / Attention the Functions give back the
    // time and the step not timen and stepn value that is why we have
    // to use < instead of <= for the while loop
    double time = sti_->GetTime();
    const double timemax = sti_->GetTimeEnd();
    int step = sti_->GetStep();
    int writestep=0;
    const int stepmax = sti_->GetTimeNumStep();
    Epetra_SerialDenseVector cvector(nmp_);

    // time loop
    while ( (time < timemax) && (step < stepmax) )
    {
      // integrate time step
      // after this step we hold disn_, etc
      sti_->IntegrateStep();

      // calculate stresses, strains, energies
      sti_->PrepareOutput();

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      sti_->UpdateStepState();

      // update time and step
      sti_->UpdateStepTime();

      // Update Element
      sti_->UpdateStepElement();

      // print info about finished time step
      sti_->PrintStep();

      // write output
      if (outputtofile) sti_->OutputStep();

      // get current time ...
      time = sti_->GetTime();
      // ... and step
      step = sti_->GetStep();

      // get the displacements of the monitored timesteps
      {
        if (abs(time-timesteps_[writestep]) < 1.0e-7)
        {
          Epetra_SerialDenseVector cvector_arg = GetCalculatedCurve(*(sti_->DisNew()),ndofs_, nnodes_, dofs_, nodes_);
          if (!myrank)
            for (int j=0; j<ndofs_; ++j)
              cvector[writestep*ndofs_+j] = cvector_arg[j];
          writestep+=1;
        }

        // check if timestepsize is smaller than the tolerance above
        double deltat = sti_->GetTimeStepSize();
        if (deltat < 1.0e-7)
          dserror("your time step size is too small, you will have problems with the monitored steps, thus adapt the tolerance");
      }
    }

    if (!(writestep*ndofs_==nmp_)) dserror("# of monitored timesteps does not match # of timesteps extracted from the simulation ");

    return cvector;
  }

/*----------------------------------------------------------------------*/
/* */
Epetra_SerialDenseVector STR::StatInvAnalysis::GetCalculatedCurve(Epetra_Vector& disp, int ndofs_, int nnodes_,std::vector<std::vector<int> > dofs_,std::vector<int>   nodes_ )
{
  int myrank = discret_->Comm().MyPID();

  // values observed at current time step
  Epetra_SerialDenseVector lcvector_arg(ndofs_);
  Epetra_SerialDenseVector gcvector_arg(ndofs_);
  int count=0;
  for (int i=0; i<nnodes_; ++i)
  {
    int gnode = nodes_[i];
    if (!discret_->NodeRowMap()->MyGID(gnode))
    {
      count += (int)dofs_[i].size();
      continue;
    }
    DRT::Node* node = discret_->gNode(gnode);
    if (!node) dserror("Cannot find node on this proc");
    if (myrank != node->Owner())
    {
      count += (int)dofs_[i].size();
      continue;
    }
    for (int j=0; j<(int)dofs_[i].size(); ++j)
    {
      int ldof = dofs_[i][j];
      int gdof = discret_->Dof(node,ldof);
      if (!disp.Map().MyGID(gdof)) dserror("Cannot find dof on this proc");
      lcvector_arg[count+j] = disp[disp.Map().LID(gdof)];
    }
    count += (int)dofs_[i].size();
  }
  discret_->Comm().SumAll(&lcvector_arg[0],&gcvector_arg[0],ndofs_);

  return gcvector_arg;
}






