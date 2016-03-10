/*!----------------------------------------------------------------------
\file drt_utils_nullspace.cpp
\brief A collection of helper methods for namespace DRT
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_nullspace.H"
#include "drt_discret.H"
#include "../drt_s8/shell8.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

#ifdef D_SHELL8
  // special for shell8
  Epetra_SerialDenseMatrix dir;
  DRT::Element* dwele = dis.lRowElement(0);
  if (dwele->ElementType()==DRT::ELEMENTS::Shell8Type::Instance())
  {
    dir.Shape(dis.NumMyRowNodes(),3);
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
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
//   if (ele->Type() == DRT::Element::element_shell8 ||
//       ele->Type() == DRT::Element::element_ale3 ||
//       ele->Type() == DRT::Element::element_so_hex8 ||
//       ele->Type() == DRT::Element::element_so_hex20 ||
//       ele->Type() == DRT::Element::element_so_hex27 ||
//       ele->Type() == DRT::Element::element_sosh8 ||
//       ele->Type() == DRT::Element::element_so_tet4 ||
//       ele->Type() == DRT::Element::element_so_tet10 ||
//       ele->Type() == DRT::Element::element_so_weg6 ||
//       ele->Type() == DRT::Element::element_sodisp ||
//       ele->Type() == DRT::Element::element_so_shw6 ||
//       ele->Type() == DRT::Element::element_truss3 ||
//       ele->Type() == DRT::Element::element_torsion3)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      std::vector<int> dofs = dis.Dof(0,actnode);  // use current dofset
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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeStructure2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

//   else if (ele->Type() == DRT::Element::element_wall1 ||
//            ele->Type() == DRT::Element::element_ale2 ||
//            ele->Type() == DRT::Element::element_torsion2)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      const double* x = actnode->X();
      std::vector<int> dofs = dis.Dof(0,actnode);
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeBeam3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* for beam3 elements the relation between rigid body modes and
   * increments on the degrees of freedom is non-trivial since
   * rotational increments in 3D are non-additive in general. In
   * general this relation may require calling all the elements.
   * However, in opposition to the SHELL8 element it is not
   * sufficient to just call a director saved in the element.
   * Rather to calculate proper increments for the rotational
   * degrees of freedom due to a rigid body rotation of the
   * complete structure, the triad at each node is required in
   * order to transform non-additive increments into additive ones.
   * However, the beam3 element currently does not save the nodal
   * triads as a class variable, but only the triads at each Gauss
   * point. In the following a wrong (!!!) dummy version is implemneted
   * but commented out. In this dummy version the rotational degrees of
   * freedom are treated identically to the additive translational
   * degrees of freedom. Activating and using this part of the code
   * quickly reveals the problems of such a naive implemnetation.
   * Usually the equation solver simply does not work with this
   * dummy code, i.e. the iterative solution process does not converge.
   * If Algebraic Multigrid methods should be really used for beam3
   * elements, one first has to develop efficient special methods for
   * these elements. Currently trying to use Algebraic multigrid methods
   * for beam3 elements just amounts to an error as no properly working
   * implementation has been available so far*/

//   else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3r)
  {
    //looping through all nodes
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      //getting pointer at current node
      DRT::Node* actnode = dis.lRowNode(i);

      //getting coordinates of current node
      const double* x = actnode->X();

      //getting number of degrees of freedom of current node
      std::vector<int> dofs = dis.Dof(actnode);

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
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] =  x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
        break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] =  x[0] - x0[0];
        break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] =  x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
        break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 0.0;
        break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 1.0;
          mode[5][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = 0.0;
          mode[5][lid] = 1.0;
        break;
        default:
          dserror("Only dofs 0 - 5 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)

  } // else if (ele->Type() == DRT::Element::element_beam3 || ele->Type() == DRT::Element::element_beam3r)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeXFluidDNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i) mode[i] = &(ns[i*lrows]);

  /* the rigid body modes for fluids are:
        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
  ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1
  */

//   else if (ele->Type() == DRT::Element::element_xfluid3 ||
//            ele->Type() == DRT::Element::element_combust3 ||
//            ele->Type() == DRT::Element::element_smoothrod ||
//            ele->Type() == DRT::Element::element_sosh8p8)
  {
    for (int i=0; i<dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      std::vector<int> dofs = dis.Dof(actnode);
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
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 6:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        case 7:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
        break;
        default:
          dserror("Only dofs 0 - 7 supported");
        break;
        } // switch (j)
      } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
    } // for (int i=0; i<NumMyRowNodes(); ++i)
  } // else if (ele->Type() == DRT::Element::element_fluid3 or ele->Type() == DRT::Element::element_xfluid3)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ComputeFluidDNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i=0; i<dimns; ++i)
    mode[i] = &(ns[i*lrows]);

  for (int i=0; i<dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    std::vector<int> dofs = dis.Dof(0,actnode);
    const unsigned int ndof = dofs.size();

    if (numdf>6) dserror("Cannot define more than 6 modes");
    for (unsigned j=0; j<ndof; ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid<0) dserror("Cannot find dof");

      for (unsigned k=0; k<ndof; ++k)
      {
        if (k%numdf == j%numdf)
          mode[k%numdf][lid] = 1.0;
        else
          mode[k%numdf][lid] = 0.0;
      }

    } // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  } // for (int i=0; i<NumMyRowNodes(); ++i)
}
