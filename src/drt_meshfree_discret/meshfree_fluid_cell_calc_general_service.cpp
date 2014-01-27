/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_calc_general_service.cpp

\brief general service routines for calculation of meshfree fluid element

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*/
/*----------------------------------------------------------------------*/

#include "meshfree_fluid_cell.H"
#include "meshfree_fluid_cell_calc.H"

#include "drt_meshfree_discret.H"
#include "drt_meshfree_node.H"        //
#include "drt_meshfree_cell_utils.H"  // to get Gauss points in real space

#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"

#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 * Evaluate supporting methods of the element
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::EvaluateService(
  DRT::ELEMENTS::MeshfreeFluid* cell,
  Teuchos::ParameterList&       params,
  Teuchos::RCP<MAT::Material>&  mat,
  DRT::Discretization&          discretization,
  std::vector<int>&             lm,
  Epetra_SerialDenseMatrix&     elemat1,
  Epetra_SerialDenseMatrix&     elemat2,
  Epetra_SerialDenseVector&     elevec1,
  Epetra_SerialDenseVector&     elevec2,
  Epetra_SerialDenseVector&     elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  switch(act)
  {
    case FLD::calc_fluid_error:
    {
      // compute error for a known analytical solution
      return ComputeError(cell, params, mat, discretization, lm, elevec1);
      break;
    }
    case FLD::integrate_shape:
    {
      return IntegrateShapeFunction(cell, discretization, lm, elevec1);
      break;
    }
    default:
      dserror("Unknown type of action for Fluid");
    break;
  } // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 * Action type: Compute Error
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ComputeError(
  DRT::ELEMENTS::MeshfreeFluid*   cell,
  Teuchos::ParameterList&         params,
  Teuchos::RCP<MAT::Material>&    mat,
  DRT::Discretization&            discretization,
  std::vector<int>&               lm,
  Epetra_SerialDenseVector&       elevec1
  )
{
  dserror("ComputeError not yet implemented for max-ent, yet");
  return 0;
}

/*----------------------------------------------------------------------*
 * Action type: Compute Error
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ComputeError(
  DRT::ELEMENTS::MeshfreeFluid*   cell,
  Teuchos::ParameterList&         params,
  Teuchos::RCP<MAT::Material>&    mat,
  DRT::Discretization&            discretization,
  std::vector<int>&               lm,
  Epetra_SerialDenseVector&       elevec1,
  const DRT::UTILS::GaussIntegration& intpoints
  )
{
  dserror("ComputeError not yet implemented for max-ent, yet");
  return 0;
}

/*----------------------------------------------------------------------------*
 * Action type: Integrate shape function
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::IntegrateShapeFunction(
  DRT::ELEMENTS::MeshfreeFluid* cell,
  DRT::Discretization&          discretization,
  const std::vector<int>&       lm,
  Epetra_SerialDenseVector&     elevec1_epetra
  )
{
  //----------------------------------------------------------------------
  // cast discretization pointer to derived meshfree type
  //----------------------------------------------------------------------

  discret_ = dynamic_cast<DRT::MESHFREE::MeshfreeDiscretization*>( &(discretization) );
  if (discret_==NULL)
    dserror("dynamic_cast of discretization to meshfree discretization failed!");

  // ---------------------------------------------------------------------------
  // construct view on element vector
  // ---------------------------------------------------------------------------

  LINALG::SerialDenseVector elevec1(View,elevec1_epetra.Values(),elevec1_epetra.Length());

  //----------------------------------------------------------------------------
  // element geometry
  //----------------------------------------------------------------------------

  // get number of nodes
  nen_ = cell->NumNode();

  // resize matrices and vectors
  funct_.LightSize(nen_);
  deriv_.LightShape(0,0);

  // get global node coordinates
  nxyz_.LightShape(nsd_,nen_);
  double const * cnxyz;
  for (int j=0; j<nen_; j++){
    cnxyz =  cell->Nodes()[j]->X();
    for (int k=0; k<nsd_; k++){
      nxyz_(k,j) = cnxyz[k];
    }
  }

  // get global point coordinates
  double const * ckxyz;
  for (int j=0; j<nek_; j++){
    ckxyz =  cell->Points()[j]->X();
    for (int k=0; k<nsd_; k++){
      kxyz_(k,j) = ckxyz[k];
    }
  }

  // set element id
  eid_ = cell->Id();

  if (cell->IsAle())
  {
    dserror("Meshfree Fluid not tested for ALE, yet. Remove dserror at own risk.");
    LINALG::SerialDenseMatrix edispnp(nsd_,nen_,true);// need to be set to true??
    ExtractValuesFromGlobalVector(discretization,lm, &edispnp, NULL,"dispnp");

    // get new node and point positions for isale
    kxyz_ += edispnp;
    nxyz_ += edispnp;
  }

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(kxyz_, gxyz_, gw_);

  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
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
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nxyz_[i];
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->GetMeshfreeSolutionApprox()->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);

    if (error)
      dserror("Something went wrong when calculating the meshfree basis functions.");

    for (int ui=0; ui<nen_; ++ui) // loop rows  (test functions)
    {
      // integrated shape function is written into the pressure dof
      elevec1(numdofpernode_ * ui + nsd_) += fac_ * funct_(ui);
    }
  }

  return 0;
}


/*-------------------------------------------------------------------------------*
 | fill element matrix and vectors with the global values              nis Nov13 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ExtractValuesFromGlobalVector(
  const DRT::Discretization&   discretization, ///< discretization
  const std::vector<int>&      lm,             ///<
  LINALG::SerialDenseMatrix *  matrixtofill,   ///< vector field
  LINALG::SerialDenseVector *  vectortofill,   ///< scalar field
  const std::string            state)          ///< state of the global vector
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
  }
}

/*-------------------------------------------------------------------------------*
 | find elements of inflow section                                     nis Nov13 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::InflowElement(DRT::Element* ele)
{
  is_inflow_ele_ = false;

  std::vector<DRT::Condition*> myinflowcond;

  // check whether all nodes have a unique inflow condition
  DRT::UTILS::FindElementConditions(ele, "TurbulentInflowSection", myinflowcond);
  if (myinflowcond.size()>1)
    dserror("More than one inflow condition on one node!");

  if (myinflowcond.size()==1)
    is_inflow_ele_ = true;

  return;
}

// template classes
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::tri3>;
