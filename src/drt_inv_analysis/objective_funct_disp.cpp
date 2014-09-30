/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "objective_funct_disp.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*/
/* standard constructor                                                 */
/*----------------------------------------------------------------------*/
STR::INVANA::ObjectiveFunctDisp::ObjectiveFunctDisp(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
timesteps_(Teuchos::null),
msteps_(0)
{
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  // this is supposed to be the number of simulation steps in the primal AND the dual problem
  msteps_ = sdyn.get<int>("NUMSTEP");
  double timestep = sdyn.get<double>("TIMESTEP");

  // initialize the vector of time steps according to the structural dynamic params
  timesteps_ = Teuchos::rcp(new std::vector<double>(msteps_,0.0));
  for (int i=0; i<=msteps_-1; i++)
  {
    (*timesteps_)[i] = (i+1)*timestep;
  }

  // initialize vectors
  mdisp_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  //mask_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  mask_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));

  ReadMonitor(invap.get<std::string>("MONITORFILE"));

}

/*----------------------------------------------------------------------*/
/* read monitor file                                        keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctDisp::ReadMonitor(std::string monitorfilename)
{
  int myrank = discret_->Comm().MyPID();

  int ndofs = 0;
  int nsteps = 0;
  int nnodes = 0;

  // list of node gids observed
  std::vector<int> nodes;
  // list of dofs on each node that are monitored
  std::vector<std::vector<int> > dofs;
  // measured displacement of the experiments (target value)
  Epetra_SerialDenseVector mcurve_;
  // time values where displacements are measured
  double timestep = 0.0;

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
  if ( nsteps > msteps_ )
    dserror("number of measured time steps greater than simulated time steps");

  // read nnodes
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
      dofs[i][j] = strtol(foundit,&foundit,10); //actdofs[strtol(foundit,&foundit,10)];
      if (!myrank) printf("%d ",dofs[i][j]);
    }
    if (!myrank) printf("\n");
    ndofs = 0;
  }

  // read comment lines
    foundit = buffer;
    fgets(buffer,150000,file);
    while(strstr(buffer,"#"))
      fgets(buffer,150000,file);

    // read in the values for each node in dirs directions
    bool readnext = true;
    int count = 0;
    for (int i=0; i<msteps_; i++)
    {
      // read the next time step in the monitor file
      if (readnext)
        timestep = strtod(foundit,&foundit);

      //check whether this timestep was simulated:
      if ( (timestep-(*timesteps_)[i]) <= 1.0e-10 )
      {
        count += 1;
        readnext = true;
      }
      else
      {
        readnext = false;
        continue;
      }

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
          mdisp_->ReplaceGlobalValue(actdofs[dofs[j][k]],i,strtod(foundit,&foundit));
          mask_->ReplaceGlobalValue(actdofs[dofs[j][k]],i,1.0);
        }
      }

      fgets(buffer,150000,file);
      foundit = buffer;
    }

    // check whether reading was successful
    if ( nsteps != count )
      dserror("check your monitor file for consistency");

}

/*----------------------------------------------------------------------*/
/* Evaluate value of the objective function                  keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctDisp::Evaluate(Teuchos::RCP<Epetra_MultiVector> disp,
                                               double& val)
{
  Epetra_SerialDenseVector normvec(disp->NumVectors());

  Epetra_MultiVector tmpvec = Epetra_MultiVector(*dofrowmap_,msteps_,true);

  // tmpvec = u_sim - u_meas;
  tmpvec.Update(1.0,*mdisp_,0.0);
  tmpvec.Multiply(1.0,*mask_,*disp,-1.0);
  // (u_sim - u_meas)^2 for every element in the vector
  tmpvec.Multiply(1.0,tmpvec,tmpvec,0.0);
  // sum over every vector in the multivector
  tmpvec.Norm1(normvec.Values());
  // sum every entry
  val = 0.5*normvec.Norm1();

}

/*----------------------------------------------------------------------*/
/* Evaluate the gradient of the objective function                      */
/* w.r.t the displacements                                   keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctDisp::EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp,
                                                       Teuchos::RCP<Epetra_MultiVector> gradient)
{
  Epetra_MultiVector tmpvec = Epetra_MultiVector(*dofrowmap_,msteps_,true);
  tmpvec.Update(1.0,*mdisp_,0.0);

  tmpvec.Multiply(1.0,*mask_,*disp,-1.0);
  gradient->Update(1.0,tmpvec,0.0);

  return;
}
