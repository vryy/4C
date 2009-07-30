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

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "linalg_utils.H"
#include "drt_globalproblem.H"
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
  // see whether we have an aztec list
  if (!solveparams.isSublist("Aztec Parameters")) return;

  // get the aztec list and see whether we use downwinding
  ParameterList& azlist = solveparams.sublist("Aztec Parameters");
  // downwinding needs nodal block information, compute it
  if (!NumMyRowElements()) dserror("Proc does not have any elements");
  DRT::Element* dwele = lRowElement(0);
  int nv=0; // number of velocity dofs
  int np=0; // number of pressure dofs
  switch (dwele->Type())
  {
    case DRT::Element::element_beam2:
      nv = 3;
    break;
    case DRT::Element::element_beam2r:
      nv = 3;
    break;
    case DRT::Element::element_beam3:
      nv = 6;
    break;
    case DRT::Element::element_shell8:
      nv = 6;
    break;
    case DRT::Element::element_wall1:
      nv = 2;
    break;
    case DRT::Element::element_sosh8:
    case DRT::Element::element_so_hex8:
    case DRT::Element::element_so_hex27:
    case DRT::Element::element_so_tet4:
    case DRT::Element::element_ptet:
    case DRT::Element::element_so_tet10:
    case DRT::Element::element_so_weg6:
    case DRT::Element::element_sodisp:
    case DRT::Element::element_so_shw6:
    case DRT::Element::element_truss3:
      nv = 3;
    break;
    case DRT::Element::element_fluid3:
      nv = 3;
      np = 1;
    break;
    case DRT::Element::element_xfluid3:
      nv = 3;
      np = 1;
    break;
    case DRT::Element::element_xdiff3:
      nv = 1;
    break;
    case DRT::Element::element_combust3:
      nv = 3;
      np = 1;
    break;
    case DRT::Element::element_fluid2:
      nv = 2;
      np = 1;
    break;
    case DRT::Element::element_transport:
      nv = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
      if (DRT::Problem::Instance(0)->ProblemType() == "elch")
      {
        nv -= 1;
        np = 1;
      }
    break;
    case DRT::Element::element_ale2:
      nv = 2;
    break;
    case DRT::Element::element_ale3:
      nv = 3;
    break;
    case DRT::Element::element_none:
    default:
      dserror("Element type not supported by ML");
    break;
  }
  if (!(nv+np)) dserror("Cannot determine nodal block size");
  azlist.set<int>("downwinding nv",nv);
  azlist.set<int>("downwinding np",np);

  // see whether we have a sublist indicating usage of Trilinos::ML
  if (!solveparams.isSublist("ML Parameters")) return;

  // see whether we have previously computed the nullspace
  // and recomputation is enforced
  ParameterList& mllist = solveparams.sublist("ML Parameters");
  RCP<vector<double> > ns =
             mllist.get<RCP<vector<double> > >("nullspace",null);
  if (ns != null && !recompute) return;

  // do the usual tests
  if (!Filled()) dserror("FillComplete was not called on discretization");
  if (!HaveDofs()) dserror("Discretization has no dofs assigned");

  // no, we have not previously computed the nullspace
  // or want to recompute it anyway
  // -> compute nullspace
  ns = null;
  mllist.set<RCP<vector<double> > >("nullspace",null);
  // ML would not tolerate this rcp-ptr in its list otherwise
  mllist.set<bool>("ML validate parameter list",false);
  const Epetra_Map* rowmap = DofRowMap();
  int numdf = 1; // default value for no. of degrees of freedom
  int dimns = 1; // default value for size of nullspace

  // get the first element of the discretization
  // Note that a processor might not have any elements
  // We assume that every proc has an element and they are of equal type
  if (!NumMyRowElements()) dserror("Proc does not have any elements");
  DRT::Element* ele = lRowElement(0);
  switch (ele->Type())
  {
    case DRT::Element::element_beam2:
    case DRT::Element::element_beam2r:
      numdf = 3;
      dimns = 3;
    break;
    case DRT::Element::element_shell8:
    case DRT::Element::element_beam3:
      numdf = 6;
      dimns = 6;
    break;
    case DRT::Element::element_wall1:
    case DRT::Element::element_truss2:
    case DRT::Element::element_torsion2:
      numdf = 2;
      dimns = 3;
    break;
    case DRT::Element::element_sosh8:
    case DRT::Element::element_so_hex8:
    case DRT::Element::element_so_hex27:
    case DRT::Element::element_so_tet4:
    case DRT::Element::element_ptet:
    case DRT::Element::element_so_tet10:
    case DRT::Element::element_so_weg6:
    case DRT::Element::element_sodisp:
    case DRT::Element::element_so_shw6:
    case DRT::Element::element_truss3:
    case DRT::Element::element_torsion3:
      numdf = 3;
      dimns = 6;
    break;
    case DRT::Element::element_fluid3:
    case DRT::Element::element_xfluid3:
    case DRT::Element::element_combust3:
      numdf = 4;
      dimns = 4;
    break;
    case DRT::Element::element_xdiff3:
      numdf = 1;
      dimns = 1;
    break;
    case DRT::Element::element_fluid2:
      numdf = 3;
      dimns = 3;
    break;
    case DRT::Element::element_transport:
      numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
      dimns = numdf;
    break;
    case DRT::Element::element_ale2:
      numdf = 2;
      dimns = 3;
    break;
    case DRT::Element::element_ale3:
      numdf = 3;
      dimns = 6;
    break;
    case DRT::Element::element_none:
    default:
      dserror("Element type not supported by ML");
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
  mllist.set<RCP<vector<double> > >("nullspace",ns);
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
      DRT::ELEMENTS::Shell8* s8 =
        dynamic_cast<DRT::ELEMENTS::Shell8*>(actnode->Elements()[0]);
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
  if (ele->Type() == DRT::Element::element_shell8 ||
      ele->Type() == DRT::Element::element_ale3 ||
      ele->Type() == DRT::Element::element_so_hex8 ||
      ele->Type() == DRT::Element::element_sosh8 ||
      ele->Type() == DRT::Element::element_so_tet4 ||
      ele->Type() == DRT::Element::element_so_tet10 ||
      ele->Type() == DRT::Element::element_so_weg6 ||
      ele->Type() == DRT::Element::element_sodisp ||
      ele->Type() == DRT::Element::element_so_shw6 ||
      ele->Type() == DRT::Element::element_truss3 ||
      ele->Type() == DRT::Element::element_torsion3)
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


  else if (ele->Type() == DRT::Element::element_wall1 ||
           ele->Type() == DRT::Element::element_ale2 ||
           ele->Type() == DRT::Element::element_torsion2)
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


  /* the variable mode[i][j] describes how the position of a
   * node changes with respect to the j-th degree of freedom
   * in case that the i-th rigid body mode is applied to the
   * structure; the structure altogether always has 3 rigid body
   * modes in R^2 and 6 in R^3; these modes are translation in
   * each coordinate direction, respectively, and translation
   * around each axis, respectively. This is summed up in the
   * following table where in the left column x,y,z denote
   * translations in x-, y- and z-direction of a node due to
   * the application of a rigid body mode, whereas dx,dy,dz
   * denote increments of the node's rotational degrees of
   * freedom, which relate to a rotation around the x-,y-
   * and z-axis.
   *
        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
  -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       1          0          0
  dy  |    0       0       0       0          1          0
  dz  |    0       0       0       0          0          1

  for example the first line means: a translation of a node in
  x-direction may be caused either by a x-translation of the whole
  structure (which is rigid body mode 0) or by a rotation either
  around the y-axis or the z-axis. In case of a rotation dtheta around the
  y-axis for example the resulting x-translation is dtheta times the
  lever arm z - z0. Here z0 represents the z-coordinate of the point
  around which the structure is rotated. This point may be chosen
  arbitrarily and by the algorithms underlying to this method it is
  chosen automatically according to some mathematical considerations.
  Note that this holds true for infinitesimal rotations dtehta, only, of
  course.
  On the other hand e.g. the fourth column means that a rigid body
  rotation dtheta around the x-axis entails translations (-z+z0)*dtheta
  in y-direction and (y-y0)*dtheta in z-direction and a rotation
  increment 1*dtheta of the rotational degree of freedom related to the
  x-axis.
  */

  /* for beam2 elements the above table reduces to
   *
        xtrans   ytrans    zrot
        mode[0]  mode[1]   mode[2]
  -----------------------------------------------------------
  x   |    1       0       -y+y0
  y   |    0       1       x-x0
  dz  |    0       0       1
  note: for the here employed Timoshenko beam elements a rigid body
  rotation entails also an increment of the rotation degree of freedom
  dz which makes the director of the beam move accordingly; only then
  a rotation does not lead to any shear stress and is truely a rigid
  body rotation
  */

  /*two dimensional beam beam elements, where each node has 2 translational
   * degrees of freedom and one rotational degree of freedom*/
  else if (ele->Type() == DRT::Element::element_beam2 ||
           ele->Type() == DRT::Element::element_beam2r)
  {
    //looping through all nodes
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      vector<int> dofs = Dof(actnode);

      //looping through all degrees of freedom of a node
      for (unsigned j=0; j<dofs.size(); ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");
        // j is degree of freedom; each case refers to one line in the above table
        switch (j)
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
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_beam2 || ele->Type() == DRT::Element::element_beam2r)

  /* for beam3 elements the relation between rigid body modes and
   * increments on the degrees of freedom is non-trivial since
   * rotational increments in 3D are non-additive in general. In
   * general this relation may require calling all the elements.
   * However, in opposition to the SHELL8 element it is not
   * sufficient to just call a director saved in the element.
   * Rather for 3D non-linear beam elements a more complicated
   * relation has to be employed. As no use of beam3 elements
   * together with Algebraic Multigrid is supposed to be required
   * this case has not been treated here and if support is called
   * then just a dserror is returned*/

 else if (ele->Type() == DRT::Element::element_beam3)
 {
   dserror("No Algebraic Multigrid support by beam3 element");
 } // else if (ele->Type() == DRT::Element::element_beam3)


  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
  ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1
  */

  else if (ele->Type() == DRT::Element::element_fluid3 or ele->Type() == DRT::Element::element_xfluid3 or ele->Type() == DRT::Element::element_combust3)
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
  } // else if (ele->Type() == DRT::Element::element_fluid3 or ele->Type() == DRT::Element::element_xfluid3)

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
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 2 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid2)

  else if (ele->Type() == DRT::Element::element_transport)
  {
    for (int i=0; i<NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = lRowNode(i);
      vector<int> dofs = Dof(actnode);
      const unsigned int ndof = dofs.size();
      for (unsigned j=0; j<ndof; ++j)
      {
        const int dof = dofs[j];
        const int lid = rowmap->LID(dof);
        if (lid<0) dserror("Cannot find dof");

        for (unsigned k=0; k<ndof; ++k)
        {
          if (k == j)
            mode[k][lid] = 1.0;
          else
            mode[k][lid] = 0.0;
        }

      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_transport)


  else ; // do nothing

  return;
}




#endif  // #ifdef CCADISCRET
