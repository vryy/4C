/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for low Mach number problems

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_boundary_calc_loma.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_fluid_rotsym_periodicbc.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_arrhenius_pv.hpp"
#include "4C_mat_arrhenius_temp.hpp"
#include "4C_mat_ferech_pv.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_mixfrac.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_mat_tempdepwater.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_mat_yoghurt.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcLoma<distype, probdim>>(
            new ScaTraEleBoundaryCalcLoma<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::ScaTraEleBoundaryCalcLoma(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::evaluate_action(
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
    case SCATRA::BoundaryAction::calc_loma_therm_press:
    {
      calc_loma_therm_press(ele, params, discretization, la);

      break;
    }

    default:
    {
      my::evaluate_action(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch (action)

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate loma therm pressure                              vg 03/09  |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::calc_loma_therm_press(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // get location vector associated with primary dofset
  std::vector<int>& lm = la[0].lm_;

  DRT::Element* parentele = ele->parent_element();
  // we dont know the parent element's lm vector; so we have to build it here
  const int nenparent = parentele->num_node();
  std::vector<int> lmparent(nenparent);
  std::vector<int> lmparentowner;
  std::vector<int> lmparentstride;
  parentele->LocationVector(discretization, lmparent, lmparentowner, lmparentstride);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = my::scatraparams_->NdsVel();

  // get velocity values at nodes
  const Teuchos::RCP<const Epetra_Vector> convel =
      discretization.GetState(ndsvel, "convective velocity field");

  // safety check
  if (convel == Teuchos::null) FOUR_C_THROW("Cannot get state vector convective velocity");

  // get values of velocity field from secondary dof-set
  const std::vector<int>& lmvel = la[ndsvel].lm_;
  std::vector<double> myconvel(lmvel.size());

  // extract local values of the global vectors
  CORE::FE::ExtractMyValues(*convel, myconvel, lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  my::rotsymmpbc_->rotate_my_values_if_necessary(myconvel);

  // define vector for normal diffusive and velocity fluxes
  std::vector<double> mynormdiffflux(lm.size());
  std::vector<double> mynormvel(lm.size());

  // determine constant outer normal to this element
  my::normal_ = my::get_const_normal(my::xyze_);

  // extract temperature flux vector for each node of the parent element
  CORE::LINALG::SerialDenseMatrix eflux(3, nenparent);
  DRT::Element* peleptr = (DRT::Element*)parentele;
  int k = my::numscal_ - 1;  // temperature is always last degree of freedom!!
  std::ostringstream temp;
  temp << k;
  std::string name = "flux_phi_" + temp.str();
  // try to get the pointer to the entry (and check if type is Teuchos::RCP<Epetra_MultiVector>)
  Teuchos::RCP<Epetra_MultiVector>* f = params.getPtr<Teuchos::RCP<Epetra_MultiVector>>(name);
  // check: field has been set and is not of type Teuchos::null
  if (f != nullptr)
    CORE::FE::ExtractMyNodeBasedValues(peleptr, eflux, *f, 3);
  else
    FOUR_C_THROW("MultiVector %s has not been found!", name.c_str());

  // calculate normal diffusive and velocity flux at each node of the
  // present boundary element
  for (int i = 0; i < nen_; ++i)
  {
    for (int j = 0; j < nenparent; ++j)
    {
      mynormdiffflux[i] = 0.0;
      mynormvel[i] = 0.0;
      for (int l = 0; l < nsd_; l++)
      {
        mynormdiffflux[i] += eflux(l, j) * my::normal_(l);
        mynormvel[i] += myconvel[i * nsd_ + l] * my::normal_(l);
      }
    }
  }

  // calculate integral of normal diffusive and velocity flux
  // NOTE: add integral value only for elements which are NOT ghosted!
  if (ele->Owner() == discretization.Comm().MyPID())
  {
    norm_diff_flux_and_vel_integral(ele, params, mynormdiffflux, mynormvel);
  }
}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions              fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::neumann_inflow(
    const DRT::FaceElement* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs)
{
  // set thermodynamic pressure
  thermpress_ = params.get<double>("thermodynamic pressure");

  // call base class routine
  my::neumann_inflow(ele, params, discretization, la, emat, erhs);

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::neumann_inflow


/*----------------------------------------------------------------------*
 | get density at integration point                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::get_density(
    Teuchos::RCP<const CORE::MAT::Material> material,
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephinp, const int k)
{
  // initialization
  double density(0.);

  // get density depending on material
  switch (material->MaterialType())
  {
    case CORE::Materials::m_matlist:
    {
      const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

      const int lastmatid = actmat->NumMat() - 1;

      if (actmat->MaterialById(lastmatid)->MaterialType() == CORE::Materials::m_arrhenius_temp)
      {
        // compute temperature and check whether it is positive
        const double temp = my::funct_.Dot(ephinp[my::numscal_ - 1]);
        if (temp < 0.0)
          FOUR_C_THROW(
              "Negative temperature in ScaTra Arrhenius temperature density evaluation on "
              "boundary!");

        // compute density based on temperature and thermodynamic pressure
        density = static_cast<const MAT::ArrheniusTemp*>(actmat->MaterialById(lastmatid).get())
                      ->ComputeDensity(temp, thermpress_);
      }
      else
        FOUR_C_THROW(
            "Type of material found in material list not supported, should be Arrhenius-type "
            "temperature!");

      break;
    }

    case CORE::Materials::m_mixfrac:
    {
      // compute density based on mixture fraction
      density = static_cast<const MAT::MixFrac*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case CORE::Materials::m_sutherland:
    {
      // compute temperature and check whether it is positive
      const double temp = my::funct_.Dot(ephinp[k]);
      if (temp < 0.0)
        FOUR_C_THROW("Negative temperature in ScaTra Sutherland density evaluation on boundary!");

      // compute density based on temperature and thermodynamic pressure
      density =
          static_cast<const MAT::Sutherland*>(material.get())->ComputeDensity(temp, thermpress_);

      break;
    }

    case CORE::Materials::m_tempdepwater:
    {
      // compute temperature and check whether it is positive
      const double temp = my::funct_.Dot(ephinp[k]);
      if (temp < 0.0)
        FOUR_C_THROW(
            "Negative temperature in ScaTra temperature-dependent water density evaluation on "
            "boundary!");

      // compute density based on temperature
      density = static_cast<const MAT::TempDepWater*>(material.get())->ComputeDensity(temp);

      break;
    }

    case CORE::Materials::m_arrhenius_pv:
    {
      // compute density based on progress variable
      density = static_cast<const MAT::ArrheniusPV*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case CORE::Materials::m_ferech_pv:
    {
      // compute density based on progress variable
      density = static_cast<const MAT::FerEchPV*>(material.get())
                    ->ComputeDensity(my::funct_.Dot(ephinp[k]));

      break;
    }

    case CORE::Materials::m_thermostvenant:
    {
      // get constant density
      density = static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get())->Density();

      break;
    }

    case CORE::Materials::m_yoghurt:
    {
      // get constant density
      density = static_cast<const MAT::Yoghurt*>(material.get())->Density();

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid material type!");
      break;
    }
  }

  return density;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::GetDensity


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux and velocity     vg 09/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::norm_diff_flux_and_vel_integral(
    const DRT::Element* ele, Teuchos::ParameterList& params,
    const std::vector<double>& enormdiffflux, const std::vector<double>& enormvel)
{
  // get variables for integrals of normal diffusive flux and velocity
  double normdifffluxint = params.get<double>("normal diffusive flux integral");
  double normvelint = params.get<double>("normal velocity integral");

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
  {
    const double fac = my::eval_shape_func_and_int_fac(intpoints, gpid);

    // compute integral of normal flux
    for (int node = 0; node < nen_; ++node)
    {
      normdifffluxint += my::funct_(node) * enormdiffflux[node] * fac;
      normvelint += my::funct_(node) * enormvel[node] * fac;
    }
  }  // loop over integration points

  // add contributions to the global values
  params.set<double>("normal diffusive flux integral", normdifffluxint);
  params.set<double>("normal velocity integral", normvelint);

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::norm_diff_flux_and_vel_integral


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<CORE::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
