/*!----------------------------------------------------------------------
\file drt_discret_utils.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "linalg_utils.H"
#ifdef D_SHELL8
#include "../drt_s8/shell8.H"
#endif


/*----------------------------------------------------------------------*
 |  compute nullspace of system (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::ComputeNullSpaceIfNecessary(
                                              ParameterList& solveparams,
                                              bool recompute)
{
  // see whether we have a sublist indicating usage of Trilinos::ML
  if (!solveparams.isSublist("ML Parameters")) return;

  // see whether we have previously computed the nullspace
  // and recomputation is enforced
  ParameterList& mllist = solveparams.sublist("ML Parameters");
  RefCountPtr<vector<double> > ns =
             mllist.get<RefCountPtr<vector<double> > >("nullspace",null);
  if (ns != null && !recompute) return;

  // do the usual tests
  if (!Filled()) dserror("FillComplete was not called on discretization");
  if (!HaveDofs()) dserror("Discretization has no dofs assigned");

  // no, we have not previously computed the nullspace
  // or want to recompute it anyway
  // -> compute nullspace
  ns = null;
  mllist.set<RefCountPtr<vector<double> > >("nullspace",null);
  const Epetra_Map* rowmap = DofRowMap();
  int numdf = 1; // default value for no. of degrees of freedom
  int dimns = 1; // default value for sirze of nullspace

  // get the first element of the discretization
  // Note that a processor might not have any elements
  // We assume that every proc has an element and they are of equal type
  if (!NumMyRowElements()) dserror("Proc does not have any elements");
  DRT::Element* ele = lRowElement(0);
  switch (ele->Type())
  {
    case DRT::Element::element_shell8:
      numdf = 6;
      dimns = 6;
    break;
    case DRT::Element::element_wall1:
      numdf = 2;
      dimns = 3;
    break;
    case DRT::Element::element_fluid3:
      numdf = 4;
      dimns = 4;
    break;
    case DRT::Element::element_fluid2:
      numdf = 3;
      dimns = 3;
    break;
    case DRT::Element::element_condif2:
      numdf = 1;
      dimns = 1;
    break;
    case DRT::Element::element_none:
    default:
      cout << "WARNING: Computation of nullspace not yet implemented for your discretization\n";
      fflush(stdout);
      numdf = 1;
      dimns = 1;
    break;
  }
  const int numproc = Comm().NumProc();
  int sumnumdf;
  Comm().SumAll(&numdf,&sumnumdf,1);
  if (sumnumdf != numdf*numproc) dserror("numdf not consistent among procs");
  int sumdimns;
  Comm().SumAll(&dimns,&sumdimns,1);
  if (sumdimns != dimns*numproc) dserror("dimns not consistent among procs");

  if (dimns>6) dserror("Nullspace size only upto 6 supported");
  mllist.set("PDE equations",numdf);
  mllist.set("null space: dimension",dimns);
  mllist.set("null space: type","pre-computed");
  mllist.set("null space: add default vectors",false);
  // allocate dimns times the local length of the rowmap
  const int lrows = rowmap->NumMyElements();
  ns = rcp(new vector<double>(dimns*lrows));
  double* nullsp = &((*ns)[0]);
  mllist.set<RefCountPtr<vector<double> > >("nullspace",ns);
  mllist.set("null space: vectors",nullsp);

  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(nullsp[i*lrows]);

  if (dimns==1 && numdf==1)
  {
    for (int i=0; i<lrows; ++i) nullsp[i] = 1.0;
    return;
  }

  // nodal center of the discretization
  double x0send[3] = {0.0,0.0,0.0};
  double x0[3];
  for (int i=0; i<NumMyRowNodes(); ++i)
    for (int j=0; j<3; ++j) x0send[j] += lRowNode(i)->X()[j];
  Comm().SumAll(x0send,x0,3);
  for (int i=0; i<3; ++i) x0[i] /= NumGlobalNodes();

#ifdef D_SHELL8
  // special for shell8
  Epetra_SerialDenseMatrix dir;
  if (ele->Type()==DRT::Element::element_shell8)
  {
    dir.Shape(NumMyRowNodes(),3);
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      DRT::Elements::Shell8* s8 =
        dynamic_cast<DRT::Elements::Shell8*>(actnode->Elements()[0]);
      if (!s8) dserror("Cannot cast to Shell8");
      int j;
      for (j=0; j<s8->NumNode(); ++j)
        if (s8->Nodes()[j]->Id() == actnode->Id()) break;
      if (j==s8->NumNode()) dserror("Can't find matching node - weird!");
      double h2 = (*s8->GetThickness())[j]/2.0;
      // get director
      const Epetra_SerialDenseMatrix* a3ref = s8->GetDirectors();
      dir(i,0) = (*a3ref)(0,j)*h2;
      dir(i,1) = (*a3ref)(1,j)*h2;
      dir(i,2) = (*a3ref)(2,j)*h2;
    }
  }
#endif

  /* the rigid body modes for structures are:
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       0          a3        -a2
  dy  |    0       0       0      -a3         0          a1
  dz  |    0       0       0       a2        -a1         0
  */

  // works straight for bricks as well
  if (ele->Type() == DRT::Element::element_shell8)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
#ifdef D_SHELL8
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = dir(i,2);
          mode[5][lid] = -dir(i,1);
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -dir(i,2);
          mode[4][lid] = 0.0;
          mode[5][lid] = dir(i,0);
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = dir(i,1);
          mode[4][lid] = -dir(i,0);
          mode[5][lid] = 0.0;
        break;
#endif
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // if (ele->Type() == DRT::Element::element_shell8)


  else if (ele->Type() == DRT::Element::element_wall1)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      const double* x = actnode->X();
      vector<int> dofs = Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
        break;
        default:
          dserror("Only dofs 0 - 1 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_wall1)


  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
  ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1
  */

  else if (ele->Type() == DRT::Element::element_fluid3)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      vector<int> dofs = Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 3 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid3)

  else if (ele->Type() == DRT::Element::element_fluid2)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      vector<int> dofs = Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = 0.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid2)

  else if (ele->Type() == DRT::Element::element_condif2)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      vector<int> dofs = Dof(actnode);
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        switch (j) // j is degree of freedom
        {
        case 0:
          mode[0][lid] = 1.0;
        break;
        default:
          dserror("Only dof 0 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_condif2)

  else ; // do nothing

  return;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
