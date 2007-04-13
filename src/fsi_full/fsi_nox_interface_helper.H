#ifndef FSI_NOX_INTERFACE_HELPER_H
#define FSI_NOX_INTERFACE_HELPER_H


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
//! \brief loop the interface and call the operator for all dofs of
//         this processor.
/*----------------------------------------------------------------------*/
template <class Op>
void loop_interface(FIELD* structfield, Op& op, const Epetra_Vector& x)
{
  int disnum = 0;
  int numnp_total = structfield->dis[disnum].numnp;
  int numaf       = genprob.numaf;

  // We don't have a direct mapping, so loop all structural nodes and
  // select the dofs that belong to the local part of the interface.

  for (int i=0;i<numnp_total;i++)
  {
    NODE* actsnode  = &(structfield->dis[disnum].node[i]);
    GNODE* actsgnode = actsnode->gnode;

    NODE* actanode  = actsgnode->mfcpnode[numaf];
    if (actanode == NULL) continue;
    GNODE* actagnode = actanode->gnode;

    /* check for coupling nodes */
    if (actagnode->dirich == NULL)
      dserror("no dirich condition for coupled ALE node #%d",actanode->Id);

    if (actagnode->dirich->dirich_type != DIRICH_CONDITION::dirich_FSI) continue;

    for (int j=0;j<actanode->numdf;++j)
    {
      if (x.Map().MyGID(actsnode->dof[j]))
      {
        op(actsnode,j,x.Map().LID(actsnode->dof[j]));
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeDisplacements
{
  const Epetra_Vector& x_;
  int pos_;

public:
  DistributeDisplacements(const Epetra_Vector& x,int pos)
    : x_(x), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { actsnode->sol_mf.a.da[pos_][dof] = x_[LID]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeForces
{
  const Epetra_Vector& x_;
  int pos_;

public:
  DistributeForces(const Epetra_Vector& x,int pos)
    : x_(x), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    {
      NODE* actfnode = actsnode->gnode->mfcpnode[genprob.numff];
      actfnode->sol_mf.a.da[pos_][dof] = x_[LID];
    }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeDelta
{
  const Epetra_Vector& x_;
  const Epetra_Vector& oldx_;
  int pos_;

public:
  DistributeDelta(const Epetra_Vector& x,const Epetra_Vector& oldx,int pos)
    : x_(x), oldx_(oldx), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { actsnode->sol_mf.a.da[pos_][dof] = x_[LID]-oldx_[LID]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class DistributeSolution
{
  const Epetra_Vector& x_;
  const Epetra_Vector& F_;
  int pos_;

public:
  DistributeSolution(const Epetra_Vector& x,const Epetra_Vector& F,int pos)
    : x_(x), F_(F), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { actsnode->sol_mf.a.da[pos_][dof] = F_[LID]+x_[LID]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherResiduum
{
  Epetra_Vector& F_;
  int dispnp_;
  int dispn_;

public:
  GatherResiduum(Epetra_Vector& F,int dispnp,int dispn)
    : F_(F), dispnp_(dispnp), dispn_(dispn) {}

  void operator()(NODE* actsnode, int dof, int LID)
    {
      F_[LID] = actsnode->sol_mf.a.da[dispnp_][dof] -
                actsnode->sol_mf.a.da[dispn_ ][dof];
    }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherRelaxation
{
  Epetra_Vector& F_;
  int pos_;

public:
  GatherRelaxation(Epetra_Vector& F,int pos)
    : F_(F), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    {
      F_[LID] = actsnode->sol.a.da[pos_][dof];
    }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherDisplacements
{
  Epetra_Vector& soln_;
  int pos_;

public:
  GatherDisplacements(Epetra_Vector& soln,int pos)
    : soln_(soln), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    { soln_[LID] = actsnode->sol_mf.a.da[pos_][dof]; }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class GatherForces
{
  Epetra_Vector& soln_;
  int pos_;

public:
  GatherForces(Epetra_Vector& soln,int pos)
    : soln_(soln), pos_(pos) {}

  void operator()(NODE* actsnode, int dof, int LID)
    {
      NODE* actfnode = actsnode->gnode->mfcpnode[genprob.numff];
      soln_[LID] = actfnode->sol_mf.a.da[pos_][dof];
    }
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class OutputDisplacements
{
  std::ostream& out_;
  int dispnp_;
  int dispn_;

public:
  OutputDisplacements(std::ostream& out,int dispnp,int dispn)
    : out_(out), dispnp_(dispnp), dispn_(dispn) {}

  // a hack
  void operator()(NODE* actsnode, int dof, int LID)
    { if (dof==0)
        out_ << LID << " "
             << actsnode->sol_mf.a.da[dispnp_][0] << " "
             << actsnode->sol_mf.a.da[dispnp_][1] << " "
             << actsnode->sol_mf.a.da[dispn_][0] << " "
             << actsnode->sol_mf.a.da[dispn_][1] << " "
             << actsnode->sol_mf.a.da[dispnp_][0]-actsnode->sol_mf.a.da[dispn_][0] << " "
             << actsnode->sol_mf.a.da[dispnp_][1]-actsnode->sol_mf.a.da[dispn_][1] << "\n"; }
};

#endif
