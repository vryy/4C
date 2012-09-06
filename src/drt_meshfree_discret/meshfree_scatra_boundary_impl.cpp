/*--------------------------------------------------------------------------*/
/*!
  \file meshfree_scatra_boundary_impl.cpp

  \brief Internal implementation of meshfree scalar transport boundary cells

  <pre>
  Maintainer: Keijo Nissen
  nissen@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15253
  </pre>
*/
/*--------------------------------------------------------------------------*/

#include "../drt_lib/drt_node.H"            // to get position of node
#include "meshfree_scatra_boundary_impl.H"  // class declarations
#include "meshfree_scatra_cell.H"           // to get nen_ and knots
#include "drt_meshfree_cell_utils.H"        // to get Gauss points in real space
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H" // basis function evaluation
#include "../drt_lib/drt_globalproblem.H"   // in BodyForce(): DRT::Problem::Instance()

/*==========================================================================*
 * class MeshfreeScaTraBoundaryImplInterface                                *
 *==========================================================================*/

DRT::ELEMENTS::MeshfreeScaTraBoundaryImplInterface* DRT::ELEMENTS::MeshfreeScaTraBoundaryImplInterface::Impl(
    const DRT::Element* ele,
    const enum INPAR::SCATRA::ScaTraType scatratype)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));

  switch (ele->Shape())
  {
  case DRT::Element::quad4:{
    return MeshfreeScaTraBoundaryImpl<DRT::Element::quad4>::Instance(numdofpernode);
  }
  case DRT::Element::tri3: {
    return MeshfreeScaTraBoundaryImpl<DRT::Element::tri3 >::Instance(numdofpernode);
  }
  case DRT::Element::line2:{
    return MeshfreeScaTraBoundaryImpl<DRT::Element::line2>::Instance(numdofpernode);
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*==========================================================================*
 * class MeshfreeScaTraBoundaryImpl                                         *
 *==========================================================================*/

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraBoundaryImpl<distype>::MeshfreeScaTraBoundaryImpl
(int numdofpernode,int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    velint_(),
    edispnp_(true),
    funct_(),
    deriv_(),
    diffus_(numscal_,0),
    isale_(false),
    is_stationary_(false),
    is_genalpha_(false),
    is_incremental_(true)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  singleton access method                              (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraBoundaryImpl<distype> * DRT::ELEMENTS::MeshfreeScaTraBoundaryImpl<distype>::Instance(
  const int numdofpernode,
  const int numscal
  )
{
  static MeshfreeScaTraBoundaryImpl<distype> * instance;
  if (instance==NULL)
    instance = new MeshfreeScaTraBoundaryImpl<distype>(numdofpernode,numscal);
  return instance;
}


/*--------------------------------------------------------------------------*
 |  evaluate meshfree scatra boundary cell               (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeScaTraBoundaryImpl<distype>::Evaluate(
    DRT::ELEMENTS::MeshfreeTransportBoundary* ele,
    ParameterList&                    params,
    DRT::Discretization&              discretization,
    vector<int>&                      lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra
)
{
  dserror("No non-normal Neumann boundary conditions implemented.");
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeScaTraBoundaryImpl<distype>::EvaluateNeumann(
    DRT::ELEMENTS::MeshfreeTransportBoundary*   cell,
    ParameterList&                      params,
    DRT::Discretization&                discretization,
    DRT::Condition&                     condition,
    vector<int>&                        lm,
    Epetra_SerialDenseVector&           elevec1)
{
  //----------------------------------------------------------------------
  // get number of nodes
  //----------------------------------------------------------------------
  nen_ = cell->NumNode();

  //----------------------------------------------------------------------
  // get additional state vector for ALE case: grid displacement
  //----------------------------------------------------------------------

//  isale_ = params.get<bool>("isale");
//  if (isale_)
//  {
//    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",null);
//    if (dispnp==null) dserror("Cannot get state vector 'dispnp'");
//    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_);
//    // add nodal displacements to point coordinates
//    xyze_ += edispnp_;
//  }
//  else edispnp_.Clear();

  //----------------------------------------------------------------------
  // evaluate appropriate (time dependent?) boundary condition factor
  //----------------------------------------------------------------------

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  //----------------------------------------------------------------------
  // get switches, values, and spatial functions from the condition
  // (assumed to be constant on element boundary)
  // (all vectors of size numdofpernode_)
  //----------------------------------------------------------------------

  // determines whether dof is activated or not
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  // points to value of boundary condition at dof
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  // points to function number if condition is spatial function at dof
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPoints<distype>::Instance().GetGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad];        // read: current gauss xyz-coordinate
    fac = gw[iquad] * curvefac; // read: current gauss weight - scaled by time factor

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = DRT::MESHFREE::UTILS::maxent_basisfunction<nsd_>(funct_,deriv_,distng,cell->Params());
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    // factor given by spatial function
    double functfac;
    // id of spatial function of boundary condition
    int functnum = -1;

    for(int dof=0;dof<numdofpernode_;dof++)
    {
      if ((*onoff)[dof]) // is this dof activated?
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];
        {
          // evaluate function at current gauss point
          if (functnum>0)
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,cgxyz,0.0,NULL);
          else
            functfac = 1.0;
        }

        const double val_fac_functfac = (*val)[dof]*fac*functfac;

        for (int node=0;node<nen_;++node)
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points

  return 0;
}
