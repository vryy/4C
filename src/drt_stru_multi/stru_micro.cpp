#if 0

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "stru_micro.H"
#include "../drt_lib/strugenalpha.H"

#include <vector>

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



StructureMulti::Micro::Micro(Teuchos::RefCountPtr<ParameterList> params,
                             Teuchos::RefCountPtr<DRT::Discretization> dis,
                             Teuchos::RefCountPtr<LINALG::Solver> solver,
                             Teuchos::RefCountPtr<IO::DiscretizationWriter> output)

  : StruGenAlpha(*params, *dis, *solver, *output),
    params_(params),
    solver_(solver),
    output_(output)
{
  // Determine dirichtoggle_ and its inverse since boundary conditions for
  // microscale simulations are due to the MicroBoundary condition and are
  // therefore not taken into account in the Constructor of StruGenAlpha

  StructureMulti::Micro::DetermineToggle();
}


void StructureMulti::Micro::DetermineToggle()
{
  RefCountPtr<DRT::Discretization> dis =
    DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  vector<DRT::Condition*> conds;
  dis->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!dis->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = dis->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      vector<int> dofs = dis->Dof(actnode);
      const unsigned numdf = dofs.size();

      for (unsigned j=0; j<numdf; ++j)
      {
        const int gid = dofs[j];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*dirichtoggle_)[lid] = 1.0;
      }
    }
  }
}

void EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd)
{
   RefCountPtr<Epetra_Vector> systemvector = disn_;

   RefCountPtr<DRT::Discretization> dis =
    DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  vector<DRT::Condition*> conds;
  dis->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!dis->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = dis->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);

      // nodal coordinates
      double x[3] = actnode->X();

      // boundary displacements are prescribed via the macroscopic
      // deformation gradient
      double disn_prescribed[3];
      for (int i=0; i<3;i++)
      {
        for (int j=0;j<3;j++)
        {
          disn_prescribed[i]=defgrd[i][j]*x[j]-x[j];
        }
      }

      vector<int> dofs = dis->Dof(actnode);

      for (int k=0; k<3; ++k)
      {
        const int gid = dofs[k];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*disn_)[lid] = disn_prescribed[k];
      }
    }
  }

}

void SetOldDisp(RefCountPtr<Epetra_Vector> disp) { dis_ = disp; }

ReturnNewDisp() { return rcp(new Epetra_Vector(*disn_)); }

void ClearDisp()
{
  dis_ = null;
}


#endif
#endif

#endif
