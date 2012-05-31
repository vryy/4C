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

#include "drt_discret.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "drt_globalproblem.H"
#include "drt_elementtype.H"
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
  if (!solveparams.isSublist("Aztec Parameters") &&
      !solveparams.isSublist("Belos Parameters") &&
      !solveparams.isSublist("Stratimikos Parameters")) return;

  int numdf = 1; // default value for no. of degrees of freedom
  int dimns = 1; // default value for size of nullspace

  // downwinding needs nodal block information, compute it
  if (!NumMyRowElements()) dserror("Proc does not have any elements");
  DRT::Element* dwele = lRowElement(0);
  int nv=0; // number of velocity dofs
  int np=0; // number of pressure dofs
#if 1
  dwele->ElementType().NodalBlockInformation( dwele, numdf, dimns, nv, np );

#else
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
    case DRT::Element::element_beam3ii:
      nv = 6;
    break;
    case DRT::Element::element_beam3eb:
      nv = 6;
    break;
    case DRT::Element::element_smoothrod:
      nv = 4;
    break;
    case DRT::Element::element_shell8:
      nv = 6;
    break;
    case DRT::Element::element_wall1:
      nv = 2;
    break;
    case DRT::Element::element_sosh8:
    case DRT::Element::element_so_hex8:
    case DRT::Element::element_so_hex8fbar:
    case DRT::Element::element_so_hex20:
    case DRT::Element::element_so_hex27:
    case DRT::Element::element_so_nurbs27:
    case DRT::Element::element_so_tet4:
    case DRT::Element::element_ptet:
    case DRT::Element::element_so_tet10:
    case DRT::Element::element_so_weg6:
    case DRT::Element::element_sodisp:
    case DRT::Element::element_so_shw6:
    case DRT::Element::element_truss3:
      nv = 3;
    break;
    //TODO: Clean up after Fluid2 element died (ehrl)
    case DRT::Element::element_fluid3:
      nv = dwele->NumDofPerNode(*(dwele->Nodes()[0]))-1;
      np = 1;
    break;
    case DRT::Element::element_sosh8p8:
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
        if (nv > 1) // only when we have more than 1 dof per node!
        {
          nv -= 1;
          np = 1;
        }
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
#endif
  if (!(nv+np)) dserror("Cannot determine nodal block size");

  // store nv and np at unique location in solver parameter list
  solveparams.sublist("NodalBlockInformation").set("nv",nv); // TODO improve names
  solveparams.sublist("NodalBlockInformation").set("np",np);
  solveparams.sublist("NodalBlockInformation").set("numdf",numdf);
  solveparams.sublist("NodalBlockInformation").set("dimns",dimns);

  if(solveparams.isSublist("Aztec Parameters"))
  {

    // get the aztec list and see whether we use downwinding
    ParameterList& azlist = solveparams.sublist("Aztec Parameters");

    azlist.set<int>("downwinding nv",nv);
    azlist.set<int>("downwinding np",np);
  }
  else if(solveparams.isSublist("Belos Parameters"))
  {
    // get the belos list and see whether we use downwinding
    ParameterList& beloslist = solveparams.sublist("Belos Parameters");

    beloslist.set<int>("downwinding nv",nv);
    beloslist.set<int>("downwinding np",np);
  }
  else if(solveparams.isSublist("Stratimikos Parameters"))
  {
    // no up and downwinding supported within Stratimikos...
  }
  else
  {
    dserror("no Aztec and no Belos list");
  }

  // adapt ML settings (if ML preconditioner is used)
  // see whether we have a sublist indicating usage of Trilinos::ML
  if (!solveparams.isSublist("ML Parameters") &&
      !solveparams.isSublist("MueLu Parameters") &&
      !solveparams.isSublist("MueLu (Contact) Parameters") &&
	  !solveparams.isSublist("Stratimikos Parameters")) return;
  ParameterList* mllist_ptr = NULL;
  if (solveparams.isSublist("Stratimikos Parameters"))
  {
	// TODO: what about MueLu?
    if (solveparams.sublist("Stratimikos Parameters").get<string>("Preconditioner Type") != "ML")
        return;
    else
      mllist_ptr = &(solveparams.sublist("Stratimikos Parameters").sublist("Preconditioner Types").sublist("ML").sublist("ML Settings"));
  }
  else if (solveparams.isSublist("ML Parameters"))
    mllist_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
	mllist_ptr = &(solveparams.sublist("MueLu Parameters"));
  else if (solveparams.isSublist("MueLu (Contact) Parameters"))
        mllist_ptr = &(solveparams.sublist("MueLu (Contact) Parameters"));
  else return;

  // see whether we have previously computed the nullspace
  // and recomputation is enforced
  ParameterList& mllist = *mllist_ptr; //solveparams.sublist("ML Parameters");
  RCP<vector<double> > ns = mllist.get<RCP<vector<double> > >("nullspace",null);
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
  const Epetra_Map* rowmap = DofRowMap(0);

#if 0
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
    case DRT::Element::element_beam3ii:
    case DRT::Element::element_beam3eb:
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
    case DRT::Element::element_so_hex20:
    case DRT::Element::element_so_hex27:
    case DRT::Element::element_so_nurbs27:
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
    //case DRT::Element::element_fluid3:
    case DRT::Element::element_sosh8p8:
    case DRT::Element::element_xfluid3:
    case DRT::Element::element_combust3:
    case DRT::Element::element_smoothrod:
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
    case DRT::Element::element_fluid3:
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
#endif

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
  for (int i=0; i<3; ++i)
      x0[i] /= NumGlobalNodes();

  dwele->ElementType().ComputeNullSpace( *this, *ns, x0, numdf, dimns );
}

