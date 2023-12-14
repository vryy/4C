/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 2


 */
/*----------------------------------------------------------------------*/

#include "baci_scatra_ele_boundary_calc_poro.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_discretization_geometry_position_array.H"
#include "baci_fluid_rotsym_periodicbc.H"
#include "baci_lib_globalproblem.H"
#include "baci_nurbs_discret.H"
#include "baci_nurbs_discret_nurbs_utils.H"
#include "baci_scatra_ele.H"
#include "baci_scatra_ele_action.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_utils_singleton_owner.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Singleton access method                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcPoro<distype, probdim>>(
            new ScaTraEleBoundaryCalcPoro<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  Private constructor                                   hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::ScaTraEleBoundaryCalcPoro(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>::ScaTraEleBoundaryCalc(
          numdofpernode, numscal, disname),
      eporosity_(true),
      eprenp_(true),
      isnodalporosity_(false)
{
  return;
}

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    SCATRA::BoundaryAction action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_fps3i_surface_permeability:
    case SCATRA::BoundaryAction::calc_fs3i_surface_permeability:
    case SCATRA::BoundaryAction::calc_Neumann:
    case SCATRA::BoundaryAction::calc_Robin:
    case SCATRA::BoundaryAction::calc_normal_vectors:
    case SCATRA::BoundaryAction::integrate_shape_functions:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
    case SCATRA::BoundaryAction::add_convective_mass_flux:
    {
      // calculate integral of convective mass/heat flux
      // NOTE: since results are added to a global vector via normal assembly
      //       it would be wrong to suppress results for a ghosted boundary!

      // get actual values of transported scalars
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");

      // extract local values from the global vector
      std::vector<CORE::LINALG::Matrix<nen_, 1>> ephinp(
          my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phinp, ephinp, la[0].lm_);

      // get number of dofset associated with velocity related dofs
      const int ndsvel = my::scatraparams_->NdsVel();

      // get convective (velocity - mesh displacement) velocity at nodes
      Teuchos::RCP<const Epetra_Vector> convel =
          discretization.GetState(ndsvel, "convective velocity field");
      if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel(nsd_ele_ * nen_, -1);
      for (int inode = 0; inode < nen_; ++inode)
        for (int idim = 0; idim < nsd_ele_; ++idim)
          lmvel[inode * nsd_ele_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // we deal with a nsd_-dimensional flow field
      CORE::LINALG::Matrix<nsd_, nen_> econvel(true);

      // extract local values of convective velocity field from global state vector
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*convel, econvel, lmvel);

      // rotate the vector field in the case of rotationally symmetric boundary conditions
      my::rotsymmpbc_->RotateMyValuesIfNecessary(econvel);

      // construct location vector for pressure dofs
      std::vector<int> lmpre(nen_, -1);
      for (int inode = 0; inode < nen_; ++inode)
        lmpre[inode] = la[ndsvel].lm_[inode * numveldofpernode + nsd_ele_];

      // extract local values of pressure field from global state vector
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*convel, eprenp_, lmpre);

      // this is a hack. Check if the structure (assumed to be the dofset 1) has more DOFs than
      // dimension. If so, we assume that this is the porosity
      if (discretization.NumDof(1, ele->Nodes()[0]) == nsd_ele_ + 2)
      {
        isnodalporosity_ = true;

        // get number of dofset associated with velocity related dofs
        const int ndsdisp = my::scatraparams_->NdsDisp();

        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(ndsdisp, "dispnp");

        if (disp != Teuchos::null)
        {
          std::vector<double> mydisp(la[ndsdisp].lm_.size());
          DRT::UTILS::ExtractMyValues(*disp, mydisp, la[ndsdisp].lm_);

          for (int inode = 0; inode < nen_; ++inode)  // number of nodes
            eporosity_(inode, 0) = mydisp[nsd_ + (inode * (nsd_ele_ + 2))];
        }
        else
          dserror("Cannot get state vector displacement");
      }
      else
        isnodalporosity_ = false;

      // for the moment we ignore the return values of this method
      CalcConvectiveFlux(ele, ephinp, econvel, elevec1_epetra);
      // vector<double> locfluxintegral = CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
      // std::cout<<"locfluxintegral[0] = "<<locfluxintegral[0]<<std::endl;

      break;
    }
    default:
    {
      dserror("Invalid action parameter nr. %i!", action);
      break;
    }
  }  // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
 | calculate integral of convective flux across boundary    vuong 07/15 |
 | (overwrites method in ScaTraEleBoundaryCalc)                         |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
std::vector<double> DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::CalcConvectiveFlux(
    const DRT::FaceElement* ele, const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp,
    const CORE::LINALG::Matrix<nsd_, nen_>& evelnp, CORE::LINALG::SerialDenseVector& erhs)
{
  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  std::vector<double> integralflux(my::numscal_);

  // loop over all scalars
  for (int k = 0; k < my::numscal_; ++k)
  {
    integralflux[k] = 0.0;

    // loop over all integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndIntFac(intpoints, iquad, &(this->normal_));

      const double porosity = ComputePorosity(ele);

      // get velocity at integration point
      my::velint_.Multiply(evelnp, my::funct_);

      // normal velocity (note: normal_ is already a unit(!) normal)
      const double normvel = my::velint_.Dot(my::normal_);

      // scalar at integration point
      const double phi = my::funct_.Dot(ephinp[k]);

      const double val = porosity * phi * normvel * fac;
      integralflux[k] += val;
      // add contribution to provided vector (distribute over nodes using shape fct.)
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        erhs[fvi] += val * my::funct_(vi);
      }
    }
  }

  return integralflux;

}  // ScaTraEleBoundaryCalcPoro<distype>::ConvectiveFlux


/*----------------------------------------------------------------------*
 |  get the material constants  (protected)                  vuong 10/14|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::ComputePorosity(
    const DRT::FaceElement* ele  //!< the element we are dealing with
)
{
  double porosity = 0.0;

  if (isnodalporosity_)
  {
    porosity = eporosity_.Dot(my::funct_);
  }
  else
  {
    dserror("porosity calculation not yet implemented for non-node-based porosity!");
  }

  return porosity;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<CORE::FE::CellType::nurbs9, 3>;

BACI_NAMESPACE_CLOSE
