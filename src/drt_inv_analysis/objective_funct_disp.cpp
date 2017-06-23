/*----------------------------------------------------------------------*/
/*!
\file objective_funct_disp.cpp

\brief Displacement based objective function

\level 3

\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
*/
/*----------------------------------------------------------------------*/

#include "objective_funct_disp.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
/* standard constructor                                                 */
/*----------------------------------------------------------------------*/
INVANA::ObjectiveFunctDisp::ObjectiveFunctDisp(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
dofrowmap_(discret->DofRowMap()),
scalefac_(1.0)
{
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();

  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");

  ReadMonitor(invap.get<std::string>("MONITORFILE"));

  scaling_ = DRT::INPUT::IntegralValue<bool>(invap, "OBJECTIVEFUNCTSCAL");

  if (scaling_)
    scalefac_ = 1.0/sqrt(mstate_->GlobalLength());

  // estimation of the variance of the measurement noise
  var_estim_ = invap.get<double>("MEASVARESTIM");
}

/*----------------------------------------------------------------------*/
/* read monitor file                                        keh 10/13   */
/*----------------------------------------------------------------------*/
void INVANA::ObjectiveFunctDisp::ReadMonitor(std::string monitorfilename)
{
  int myrank = discret_->Comm().MyPID(); //only for print out

  int ndofs = 0;
  int nsteps = 0;
  int nnodes = 0;

  // list of node gids observed
  std::vector<int> nodes;
  // list of dofs (local numbering, i.e 0, 1, 2, ...) on each node that are monitored
  std::vector<std::vector<int> > dofs;
  // dof gids of the measured dofs
  std::vector<int> mdofs;
  // time values where displacements are measured
  double timestep = 0.0;

  // open monitor file
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

  // read number of steps
  foundit = strstr(buffer,"steps"); foundit += strlen("steps");
  nsteps = strtol(foundit,&foundit,10);

  // read number of nnodes
  foundit = strstr(buffer,"nnodes");
  foundit += strlen("nnodes");
  nnodes = strtol(foundit,&foundit,10);


  // read nodes
  nodes.resize(nnodes);
  dofs.resize(nnodes);
  for (int i=0; i<nnodes; ++i)
  {
    fgets(buffer,150000,file);
    foundit = buffer;
    nodes[i] = strtol(foundit,&foundit,10);
    int ndofs_act = strtol(foundit,&foundit,10);
    ndofs += ndofs_act;
    dofs[i].resize(ndofs,-1);

    if (!myrank) printf("Monitored node %d ndofs %d dofs ",nodes[i],(int)dofs[i].size());
    for (int j=0; j<ndofs; ++j)
    {
      dofs[i][j] = strtol(foundit,&foundit,10);
      if (discret_->HaveGlobalNode(nodes[i]))
      {
        DRT::Node* actnode = discret_->gNode(nodes[i]);
        std::vector<int> actdofs = discret_->Dof(actnode);
        if (dofrowmap_->LID(actdofs[dofs[i][j]]) != -1)
          mdofs.push_back(actdofs[dofs[i][j]]);
      }
      if (!myrank) printf("%d ",dofs[i][j]);
    }
    if (!myrank) printf("\n");
    ndofs = 0;
  }

  // setup map of measured dofs
  mdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1,mdofs.size(),mdofs.data(),0,discret_->Comm()));
  mextractor_ = Teuchos::rcp(new LINALG::MapExtractor(*dofrowmap_, mdofrowmap_));

  //initialize vector of measured values
  mstate_ = Teuchos::rcp(new Epetra_MultiVector(*mdofrowmap_,nsteps,true));

  // read comment lines
  foundit = buffer;
  fgets(buffer,150000,file);
  while(strstr(buffer,"#"))
    fgets(buffer,150000,file);


  // read in the values for each monitored time step
  for (int i=0; i<nsteps; i++)
  {
    timestep = strtod(foundit,&foundit);
    timesteps_.push_back(timestep);

    for (int j=0; j<nnodes; j++)
    {
      // if this processor does not own this node, skip it
      if (!(discret_->HaveGlobalNode(nodes[j])))
      {
        for (int k=0; k<(int)dofs[j].size(); k++)
        {
          strtod(foundit,&foundit);
        }
        continue;
      }

      DRT::Node* actnode = discret_->gNode(nodes[j]);
      std::vector<int> actdofs = discret_->Dof(actnode);

      for (int k=0; k<(int)dofs[j].size(); k++)
      {
        if (dofrowmap_->LID(actdofs[dofs[j][k]]) != -1)
        {
          int err=mstate_->ReplaceGlobalValue(actdofs[dofs[j][k]],i,strtod(foundit,&foundit));
          if (err==1) dserror("row not on this proc");
        }
        else
          strtod(foundit,&foundit);
      }
    }

    fgets(buffer,150000,file);
    foundit = buffer;
  }
}

/*----------------------------------------------------------------------*/
/* find step of measurement according to given time          keh 10/14  */
/*----------------------------------------------------------------------*/
int INVANA::ObjectiveFunctDisp::FindStep(double time)
{
  // find step of the evaluation according to time:
  int step=-1;
  for (int i=0; i<(int)timesteps_.size(); i++)
  {
    double dt=abs(timesteps_[i]-time);
    if (dt<1.0e-10)
      step=i;
  }

  return step;
}

/*----------------------------------------------------------------------*/
/* Evaluate value of the objective function                  keh 10/13  */
/*----------------------------------------------------------------------*/
void INVANA::ObjectiveFunctDisp::Evaluate(Teuchos::RCP<Epetra_Vector> state,
                                               double time,
                                               double& val)
{
  int step=FindStep(time);
  if (step == -1)
    dserror("no measurements for the requested time step");

  Teuchos::RCP<Epetra_Vector> sim = mextractor_->ExtractCondVector(state);

  // u_sim - u_meas;
  sim->Update(-1.0,*(*mstate_)(step),1.0);

  // sum up
  double norm;
  sim->Norm2(&norm);

  val = norm*norm/(2.0*var_estim_);

}

/*----------------------------------------------------------------------*/
/* Evaluate the gradient of the objective function                      */
/* w.r.t the displacements                                   keh 10/13  */
/*----------------------------------------------------------------------*/
void INVANA::ObjectiveFunctDisp::EvaluateGradient(Teuchos::RCP<Epetra_Vector> state,
                                                       double time,
                                                       Teuchos::RCP<Epetra_Vector> gradient)
{
  int step=FindStep(time);
  if (step == -1)
    dserror("no measurements for the requested time step");

  Teuchos::RCP<Epetra_Vector> sim = mextractor_->ExtractCondVector(state);

  // u_sim - u_meas;
  sim->Update(-1.0,*(*mstate_)(step),1.0);

  gradient->Update(1.0/var_estim_,*mextractor_->InsertCondVector(sim),0.0);

  return;
}
