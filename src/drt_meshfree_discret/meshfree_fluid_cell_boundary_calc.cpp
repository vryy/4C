/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_calc.cpp

\brief evaluation of meshfree fluid terms at integration points

\maintainer Keijo Nissen

\level 3

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell_boundary_calc.H"
#include "meshfree_fluid_cell.H"
#include "drt_meshfree_node.H"
#include "drt_meshfree_cell_utils.H"
#include "drt_meshfree_discret.H"
#include "drt_meshfree_utils.H"

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
  nxyz_(Teuchos::rcp(new LINALG::SerialDenseMatrix())),
  pxyz_(Teuchos::rcp(new LINALG::SerialDenseMatrix(bdrynsd_,bdrynep_))),
  gxyz_(bdrynsd_,ngp_),
  gw_(ngp_),
  funct_(Teuchos::rcp(new LINALG::SerialDenseVector())),
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

  // through dserror for ALE
  if (cell->ParentElement()->IsAle())
    dserror("Can't handle ALE for meshfree Neumann boundaries, yet!");

  // ---------------------------------------------------------------------
  // set size of all vectors of SerialDense element arrays
  // ---------------------------------------------------------------------

  // get number of nodes of boundary cell
  bdrynen_ = cell->NumNode();

  // resize matrices and vectors
  funct_->LightSize(bdrynen_);
  nxyz_->LightShape(nsd_,bdrynen_);
  pxyz_->LightShape(nsd_,bdrynep_);

  //------------------------------------------------------------------------
  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  //------------------------------------------------------------------------

  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  //------------------------------------------------------------------------
  // compute factor of curve time and density
  //------------------------------------------------------------------------

  // find out whether we will use a time curve
  const double time = fldparatimint_->Time();

  //------------------------------------------------------------------------
  // get local node coordinates
  //------------------------------------------------------------------------

  // get nodes position in dimensions of parent cell
  double const * cnxyz;
  for (int j=0; j<bdrynen_; j++){
    cnxyz =  cell->Nodes()[j]->X();
    for (int k=0; k<nsd_; k++){
      (*nxyz_)(k,j) = cnxyz[k];
    }
  }

  // fill matrix of cell point positions in reduced dimensions
  double const * cpxyz;
  for (int j=0; j<bdrynep_; j++){
    cpxyz =  cell->Points()[j]->X();
    for (int k=0; k<nsd_; ++k){
      (*pxyz_)(k,j) = cpxyz[k];
    }
  }

  //------------------------------------------------------------------------
  // reduce dimension of nodal positions to dimension of boundary cell
  // (so far only plane boundaries aligned to coordinates axes can be
  // handeled)
  //------------------------------------------------------------------------

  // get matrix of cell node positions in reduced dimensions
  DRT::MESHFREE::MovePointOfReferenceToFace(nxyz_,pxyz_);
  DRT::MESHFREE::ReduceCloudDimensionOnFaces(nxyz_,pxyz_);
  if (nxyz_->M()!=bdrynsd_)
  {
    std::cout << "nxyz_: " << std::endl << *nxyz_ << std::endl;
    dserror("Dimension of boundary node cloud is %i instead of %i.",nxyz_->M(),bdrynsd_);
  }
  if (pxyz_->M()!=bdrynsd_)
  {
    std::cout << "pxyz_: " << std::endl << *pxyz_ << std::endl;
    dserror("Dimension of boundary point cloud is %i instead of %i.",pxyz_->M(),bdrynsd_);
  }

  //------------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //------------------------------------------------------------------------

  // get integrations points and weights in xyz-system
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(*pxyz_, gxyz_, gw_);
  if (ngp<0)
  {
    std::cout << "For cell " << cell->Id() << " with points";
    for(int i=0; i< bdrynep_; ++i)
      std::cout << " " << cell->Points()[i]->Id();
    std::cout << ":" << std::endl;
    dserror("Jacobi determinant is negative! Check elements in input file!");
  }

  LINALG::SerialDenseMatrix distng(bdrynsd_,bdrynen_); // matrix for distance between node and Gauss point
  double const * cgxyz; // read: current Gauss xyz-coordinate

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
      cnxyz = (*nxyz_)[i];
      for (int j=0; j<bdrynsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int err = discret_->GetSolutionApprox()->GetMeshfreeBasisFunction(bdrynsd_,Teuchos::rcpFromRef(distng),funct_);
    if (err>0) DRT::MESHFREE::OutputMeshfreeError(err);

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
    const double fac_time_dens = fac_ * densfac_ * fldparatimint_->TimeFacRhs();

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
          functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(idim,cgxyz,time);

        // aggregate all factors
        const double valfac = (*val)[idim]*fac_time_dens*functfac;

        // loop over all nodes
        for(int inode=0; inode < bdrynen_; ++inode )
          elevec1_epetra[inode*numdofpernode_+idim] += (*funct_)(inode)*valfac;

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
