/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_calc.cpp

\brief evaluation of meshfree fluid terms at integration points

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell_boundary_calc.H"
#include "meshfree_fluid_cell.H"
#include "drt_meshfree_node.H"              // to get Gauss points in real space
#include "drt_meshfree_cell_utils.H"              // to get Gauss points in real space
#include "drt_meshfree_discret.H"              // to get Gauss points in real space

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_parameter_std.H"
#include "../drt_fluid_ele/fluid_ele_parameter_timint.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*
 |  Constructor                                       (public) nis Nov13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype>::MeshfreeFluidBoundaryCalc():
  kxyz_(bdrynsd_,bdrynek_),
  gxyz_(bdrynsd_,ngp_),
  gw_(ngp_),
  drs_(0.0),
  fac_(0.0),
  visc_(0.0),
  densaf_(1.0)
{
  // pointer to class FluidImplParameterTimInt for time integration
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
  // initialize also general parameter list, also it will be overwritten in derived subclasses
  fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate something at the Neumann boundary         (public) nis Nov13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype>::EvaluateAction(
  DRT::ELEMENTS::MeshfreeFluidBoundary* ele1,
  Teuchos::ParameterList&               params,
  DRT::Discretization&                  discretization,
  std::vector<int>&                     lm,
  Epetra_SerialDenseMatrix&             elemat1,
  Epetra_SerialDenseMatrix&             elemat2,
  Epetra_SerialDenseVector&             elevec1,
  Epetra_SerialDenseVector&             elevec2,
  Epetra_SerialDenseVector&             elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params,"action");

  // get status of Ale
  //const bool isale = ele1->ParentElement()->IsAle();

  switch(act)
  {
  default:
  {
    dserror("Unknown type of action for MeshfreeFluidBoundaryCalc!");
    break;
  }
  } // end of switch(act)

} // end MeshfreeFluidBoundaryCalc<distype>::EvaluateAction

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition    (public) nis Nov13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype>::EvaluateNeumann(
  DRT::ELEMENTS::MeshfreeFluidBoundary* cell,
  Teuchos::ParameterList&               params,
  DRT::Discretization&                  discretization,
  DRT::Condition&                       condition,
  std::vector<int>&                     lm,
  Epetra_SerialDenseVector&             elevec1_epetra,
  Epetra_SerialDenseMatrix*             elemat1_epetra)
{
  //----------------------------------------------------------------------
  // cast discretization pointer to derived meshfree type
  //----------------------------------------------------------------------
  discret_ = dynamic_cast<DRT::MESHFREE::MeshfreeDiscretization*>( &(discretization) );
  if (discret_==NULL)
    dserror("dynamic_cast of discretization to meshfree discretization failed!");

  // ---------------------------------------------------------------------
  // set size of all vectors of SerialDense element arrays
  // ---------------------------------------------------------------------

  // get number of nodes of boundary cell
  bdrynen_ = cell->NumNode();

  // resize matrices and vectors
  funct_.LightSize(bdrynen_);
  deriv_.LightShape(bdrynsd_,bdrynen_);
  nxyz_.LightShape(bdrynsd_,bdrynen_);

  //------------------------------------------------------------------------
  // compute timefac
  //------------------------------------------------------------------------

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = fldparatimint_->Time();
  if (time<0.0) usetime = false;

  // get time-curve factor/ n = - grad phi / |grad phi|
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  // get time factor for Neumann term
  const double timefac = fldparatimint_->TimeFacRhs();

  //------------------------------------------------------------------------
  // get local node coordinates
  //------------------------------------------------------------------------

  // matrix for nodes position in dimensions of parent cell
  LINALG::SerialDenseMatrix nxyz_parent(nsd_,bdrynen_,false);
  // get nodes position in dimensions of parent cell
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::SerialDenseMatrix>(cell,nxyz_parent);

  // add potential ALE displacements
  if (cell->ParentElement()->IsAle())
  {
    dserror("can't handel ALE for meshfree Neumann boundaries, yet!");
    // get nodal ALE-displacements
    Teuchos::RCP<const Epetra_Vector>  dispnp;
    std::vector<double>                mydispnp;
    dispnp = discret_->GetState("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    // add nodal ALE-displacements to nodal positions
    for (int inode=0;inode<bdrynen_;++inode)
      for(int idim=0;idim<nsd_;++idim)
        nxyz_parent(idim,inode) += mydispnp[numdofpernode_*inode+idim];
  }

  //------------------------------------------------------------------------
  // reduce dimension of nodal positions to dimension of boundary cell
  // (so far only plane boundaries aligned to coordinates axes can be
  // handeled)
  //------------------------------------------------------------------------

  // find coordinates in which boundary lies
  std::vector<int> dims;
  int ndims = 0;
  for(int idim=0; idim<nsd_; ++idim)
  {
    for(int inode=0; inode<(bdrynen_-1); ++inode)
    {
      if (!(abs(nxyz_parent(idim,inode)-nxyz_parent(idim,inode+1))<EPS12))
      {
        if (ndims<bdrynsd_)
        {
          dims.push_back(idim);
          ndims++;
          break;
        }
        else
          dserror("So far, meshfree Neumann boundary has to be aligned with a coordinate axis!");
      }
    }
  }
  if (!(ndims==bdrynsd_))
    dserror("Too few dimensions selected for meshfree boundary cell.");

  // copy node positions for reduced dimension
  for(int idim=0; idim<bdrynsd_; ++idim)
    for(int inode=0; inode<bdrynen_; ++inode)
      nxyz_(idim,inode) = nxyz_parent(dims[idim],inode);
  // fill matrix of cell knots position in xyz
  double const * kxyz_parent;
  for (int j=0; j<bdrynek_; j++){
    kxyz_parent =  cell->Knots()[j]->X();
    for (int k=0; k<bdrynsd_; k++){
      kxyz_(k,j) = kxyz_parent[dims[k]];
    }
  }

  //------------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //------------------------------------------------------------------------

  // get integrations points and weights in xyz-system
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(kxyz_, gxyz_, gw_);

  LINALG::SerialDenseMatrix distng(bdrynsd_,bdrynen_); // matrix for distance between node and Gauss point
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate

  //------------------------------------------------------------------------
  //  loop over integration points for current cell
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    //----------------------------------------------------------------------
    // get basis function values at current Gauss point
    //----------------------------------------------------------------------

    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz_[iquad]; // read: current gauss xyz-coordinate
    fac_ = gw_[iquad];    // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<bdrynen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nxyz_[i];
      for (int j=0; j<bdrynsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,bdrynsd_);
    if (error) dserror("Something went wrong when calculating the meshfree basis functions.");

    // get the required material information
    Teuchos::RCP<MAT::Material> material = cell->ParentElement()->Material();

    //----------------------------------------------------------------------
    // get density
    //----------------------------------------------------------------------

    // (evaluation always at integration point, in contrast to parent element)
    GetDensity(material);

    //------------------------------------------------------------------------
    // compute Neumann contibution of cell rhs
    //------------------------------------------------------------------------

    // aggregate all factors but factor given by spatial function
    const double fac_curve_time_dens = fac_*curvefac*timefac*densfac_;

    // factor given by spatial function
    double functfac = 1.0;

    // number of potential spatial function
    int functnum = -1;

    // loop over all dimensions
    for(int idim=0; idim<nsd_; ++idim)
    {
      // check whether this dof is activated
      if((*onoff)[idim])
      {
        // get function number
        if (func)
          functnum = (*func)[idim];

        // evaluate function at current gauss point
        if (functnum>0)
          functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,cgxyz,time,NULL);

        // aggregate all factors
        const double valfac = (*val)[idim]*fac_curve_time_dens*functfac;

        // loop over all nodes
        for(int inode=0; inode < bdrynen_; ++inode )
          elevec1_epetra[inode*numdofpernode_+idim] += funct_(inode)*valfac;

      }  // if (*onoff)
    } // for idim
  } // for iquad

  return 0;
}

/*----------------------------------------------------------------------*
 |  get density                                                vg 06/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<distype>::GetDensity(
  Teuchos::RCP<const MAT::Material>    material
)
{
  // initially set density and density factor for Neumann boundary conditions to 1.0
  // (the latter only changed for low-Mach-number flow/combustion problems)
  densaf_  = 1.0;
  densfac_ = 1.0;

  if (material->MaterialType() == INPAR::MAT::m_fluid)
  {
    // get material
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

    // incompressible flow
    densaf_ = actmat->Density();
  }
  else
    dserror("Material type is not supported for density evaluation for meshfree boundary cell!");

  // check whether there is zero or negative density
  if (densaf_ < EPS15)
    dserror("zero or negative density!");

  return;
} // MeshfreeFluidBoundaryCalc::GetDensity

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::MeshfreeFluidBoundaryCalc<DRT::Element::line2>;
