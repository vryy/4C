/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_fluid_ele_boundary_calc.hpp"

#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_nurbs_discretization_utils.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_parameter_std.hpp"
#include "4C_fluid_ele_parameter_timint.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_material.hpp"
#include "4C_mat_carreauyasuda.hpp"
#include "4C_mat_fluid_linear_density_viscosity.hpp"
#include "4C_mat_fluid_murnaghantait.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_herschelbulkley.hpp"
#include "4C_mat_modpowerlaw.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidBoundaryImpl<distype>::FluidBoundaryImpl()
    :  // Discret::ELEMENTS::FluidBoundaryInterface(),
      xyze_(true),
      funct_(true),
      deriv_(true),
      unitnormal_(true),
      velint_(true),
      drs_(0.0),
      fac_(0.0),
      visc_(0.0),
      densaf_(1.0)
{
  // pointer to class FluidImplParameterTimInt for time integration
  fldparatimint_ = Discret::ELEMENTS::FluidEleParameterTimInt::instance();
  // initialize also general parameter list, also it will be overwritten in derived subclasses
  fldpara_ = Discret::ELEMENTS::FluidEleParameterStd::instance();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::evaluate_action(
    Discret::ELEMENTS::FluidBoundary* ele1, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = Core::UTILS::get_as_enum<FLD::BoundaryAction>(params, "action");

  // get status of Ale
  const bool isale = ele1->parent_element()->is_ale();

  switch (act)
  {
    case FLD::integrate_Shapefunction:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.get_state("dispnp");
        if (dispnp != Teuchos::null)
        {
          mydispnp.resize(lm.size());
          Core::FE::extract_my_values(*dispnp, mydispnp, lm);
        }
      }

      integrate_shape_function(ele1, params, discretization, lm, elevec1, mydispnp);
      break;
    }
    case FLD::calc_area:
    {
      if (ele1->owner() == discretization.get_comm().MyPID())
        area_calculation(ele1, params, discretization, lm);
      break;
    }
    case FLD::calc_pressure_bou_int:
    {
      if (ele1->owner() == discretization.get_comm().MyPID())
        pressure_boundary_integral(ele1, params, discretization, lm);
      break;
    }
    // general action to calculate the flow rate
    case FLD::calc_flowrate:
    {
      compute_flow_rate(ele1, params, discretization, lm, elevec1);
      break;
    }
    case FLD::flowratederiv:
    {
      flow_rate_deriv(
          ele1, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::Outletimpedance:
    {
      impedance_integration(ele1, params, discretization, lm, elevec1);
      break;
    }

    case FLD::dQdu:
    {
      d_qdu(ele1, params, discretization, lm, elevec1);
      break;
    }
    case FLD::ba_calc_node_normal:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.get_state("dispnp");
        if (dispnp != Teuchos::null)
        {
          mydispnp.resize(lm.size());
          Core::FE::extract_my_values(*dispnp, mydispnp, lm);
        }
      }

      Core::FE::element_node_normal<distype>(funct_, deriv_, fac_, unitnormal_, drs_, xsi_, xyze_,
          ele1, discretization, elevec1, mydispnp, Core::FE::is_nurbs<distype>,
          ele1->parent_element()->is_ale());
      break;
    }
    case FLD::calc_node_curvature:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      if (isale)
      {
        dispnp = discretization.get_state("dispnp");
        if (dispnp != Teuchos::null)
        {
          mydispnp.resize(lm.size());
          Core::FE::extract_my_values(*dispnp, mydispnp, lm);
        }
      }

      Teuchos::RCP<const Epetra_Vector> normals;
      std::vector<double> mynormals;

      normals = discretization.get_state("normals");
      if (normals != Teuchos::null)
      {
        mynormals.resize(lm.size());
        Core::FE::extract_my_values(*normals, mynormals, lm);
      }

      // what happens, if the mynormals vector is empty? (ehrl)
      FOUR_C_THROW(
          "the action calc_node_curvature has not been called by now. What happens, if the "
          "mynormal vector is empty");

      element_mean_curvature(ele1, params, discretization, lm, elevec1, mydispnp, mynormals);
      break;
    }
    case FLD::calc_Neumann_inflow:
    {
      neumann_inflow(ele1, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::calc_surface_tension:
    {
      // employs the divergence theorem acc. to Saksono eq. (24) and does not
      // require second derivatives.

      Teuchos::RCP<const Epetra_Vector> dispnp;
      std::vector<double> mydispnp;

      dispnp = discretization.get_state("dispnp");
      if (dispnp != Teuchos::null)
      {
        mydispnp.resize(lm.size());
        Core::FE::extract_my_values(*dispnp, mydispnp, lm);
      }

      // mynormals and mycurvature are not used in the function
      std::vector<double> mynormals;
      std::vector<double> mycurvature;

      element_surface_tension(
          ele1, params, discretization, lm, elevec1, mydispnp, mynormals, mycurvature);
      break;
    }
    case FLD::center_of_mass_calc:
    {
      // evaluate center of mass
      if (ele1->owner() == discretization.get_comm().MyPID())
        center_of_mass_calculation(ele1, params, discretization, lm);
      break;
    }
    case FLD::traction_velocity_component:
    {
      calc_traction_velocity_component(ele1, params, discretization, lm, elevec1);
      break;
    }
    case FLD::traction_Uv_integral_component:
    {
      compute_neumann_uv_integral(ele1, params, discretization, lm, elevec1);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of action for FluidBoundaryImpl!");
      break;
    }
  }  // end of switch(act)
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::FluidBoundaryImpl<distype>::evaluate_neumann(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseMatrix* elemat1_epetra)
{
  // find out whether we will use a time curve
  const double time = fldparatimint_->time();

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const auto* onoff = &condition.parameters().get<std::vector<int>>("ONOFF");
  const auto* val = &condition.parameters().get<std::vector<double>>("VAL");
  const auto* func = &condition.parameters().get<std::vector<int>>("FUNCT");
  const std::string* type = &condition.parameters().get<std::string>("TYPE");

  // get time factor for Neumann term
  const double timefac = fldparatimint_->time_fac_rhs();
  const double timefacn = (1.0 - fldparatimint_->theta()) * fldparatimint_->dt();

  // get Gaussrule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);


  // THESE ARE NEEDED IF DENSITY OR A SCALAR IS USED in the BC (which normally is NOT the case)
  //========================================================
  // get scalar vector
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.get_state("scaaf");
  if (scaaf == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'scaaf'");

  // extract local values from global vector
  std::vector<double> myscaaf(lm.size());
  Core::FE::extract_my_values(*scaaf, myscaaf, lm);

  Core::LinAlg::Matrix<bdrynen_, 1> escaaf(true);

  // insert scalar into element array
  // the scalar is stored to the pressure dof
  for (int inode = 0; inode < bdrynen_; ++inode)
  {
    escaaf(inode) = myscaaf[(nsd_) + (inode * numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermpress at n+alpha_F/n+1", 0.0);

  // extract pressure values from global velocity/pressure vector
  // (needed for weakly_compressible fluids)
  Core::LinAlg::Matrix<bdrynen_, 1> epreaf(true);
  Discret::ELEMENTS::FluidEleParameter* fldpara =
      Discret::ELEMENTS::FluidEleParameterStd::instance();
  if (fldpara->physical_type() == Inpar::FLUID::weakly_compressible or
      fldpara->physical_type() == Inpar::FLUID::weakly_compressible_stokes)
  {
    Teuchos::RCP<const Epetra_Vector> velaf = discretization.get_state("velaf");
    if (velaf == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'velaf'");

    // extract local values from global vector
    std::vector<double> myvelaf(lm.size());
    Core::FE::extract_my_values(*velaf, myvelaf, lm);

    // insert pressure into element array
    for (int inode = 0; inode < bdrynen_; inode++)
    {
      epreaf(inode) = myvelaf[nsd_ + inode * numdofpernode_];
    }
  }
  else
  {
    // insert pressure into element array
    for (int inode = 0; inode < bdrynen_; inode++)
    {
      epreaf(inode) = 0.0;
    }
  }
  //========================================================


  // add potential ALE displacements
  if (ele->parent_element()->is_ale())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < (nsd_); ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  // --------------------------------------------------

  // In the case of nurbs the normal vector is multiplied with normalfac
  double normalfac = 0.0;
  std::vector<Core::LinAlg::SerialDenseVector> myknots(bdrynsd_);
  Core::LinAlg::SerialDenseVector weights(bdrynen_);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if (Core::FE::is_nurbs<distype>)
  {
    std::vector<Core::LinAlg::SerialDenseVector> mypknots(nsd_);

    bool zero_size =
        Core::FE::Nurbs::get_knot_vector_and_weights_for_nurbs_boundary(ele, ele->surface_number(),
            ele->parent_element()->id(), discretization, mypknots, myknots, weights, normalfac);
    if (zero_size)
    {
      return 0;
    }
  }
  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < intpoints.ip().nquad; ++gpid)
  {
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, &myknots, &weights, Core::FE::is_nurbs<distype>);

    // get the required material information
    Teuchos::RCP<Core::Mat::Material> material = ele->parent_element()->material();

    // get density
    // (evaluation always at integration point, in contrast to parent element)
    get_density(material, escaaf, thermpressaf, epreaf);

    const double tol = 1e-8;
    if (densfac_ < (1.0 - tol) or densfac_ > (1.0 + tol))
    {
      // If you got this warning have a look in GetDensity(...) for clarification!
      std::cout << "                                                                    "
                << std::endl;
      std::cout << "                                                                    "
                << std::endl;
      std::cout << "          WARNING:                                                  "
                << std::endl;
      std::cout << "                  4C scales your NEUMANN BC with the DENSITY!!!   "
                << std::endl;
      std::cout << "                                                                    "
                << std::endl;
      std::cout << "                  Do you really want this?                          "
                << std::endl;
      std::cout << "                  Like, super sure about this?                      "
                << std::endl;
      std::cout << "                  Might want to think this through one more time...."
                << std::endl;
      std::cout << "                                                                    "
                << std::endl;
      std::cout << "                                                                    "
                << std::endl;
      std::cout << "        Fine... Do what you want...                                 "
                << std::endl;
      std::cout << "        densfac_=               " << densfac_ << std::endl;
    }

    const double fac_time_dens = fac_ * timefac * densfac_;
    const double fac_time_densn = fac_ * timefacn;

    // factor given by spatial function
    double functfac = 1.0;
    double functfacn = 1.0;

    // global coordinates of gausspoint
    Core::LinAlg::Matrix<(nsd_), 1> coordgp(0.0);

    // determine coordinates of current Gauss point
    coordgp.multiply(xyze_, funct_);

    // we need a 3D position vector for function evaluation!
    double coordgp3D[3];
    coordgp3D[0] = 0.0;
    coordgp3D[1] = 0.0;
    coordgp3D[2] = 0.0;
    for (int i = 0; i < nsd_; i++) coordgp3D[i] = coordgp(i);

    int functnum = -1;
    const double* coordgpref = coordgp3D;  // needed for function evaluation

    for (int idim = 0; idim < (nsd_); ++idim)
    {
      if (*type == "neum_live")
      {
        if ((*onoff)[idim])  // Is this dof activated
        {
          if (func) functnum = (*func)[idim];
          if (functnum > 0)
          {
            // evaluate function at current gauss point
            functfac = Global::Problem::instance()
                           ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                           .evaluate(coordgpref, time, idim);
            if (fldparatimint_->is_new_ost_implementation())
              functfacn = Global::Problem::instance()
                              ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                              .evaluate(coordgpref, time - fldparatimint_->dt(), idim);
          }
          else
          {
            functfac = 1.0;
            functfacn = 1.0;
          }

          const double valfac = (*val)[idim] * fac_time_dens * functfac;
          for (int inode = 0; inode < bdrynen_; ++inode)
          {
            elevec1_epetra[inode * numdofpernode_ + idim] += funct_(inode) * valfac;
            if (fldparatimint_->is_new_ost_implementation())
            {
              const double valfacn = (*val)[idim] * fac_time_densn * functfacn;
              elevec1_epetra[inode * numdofpernode_ + idim] += funct_(inode) * valfacn;
            }
          }  // end is_new_ost_implementation
        }    // if (*onoff)
      }
      else if (*type == "neum_pseudo_orthopressure")
      {
        if (idim != 0 and (*onoff)[idim])
          FOUR_C_THROW(
              "If you apply a pseudo_orthopressure load on the fluid, only a load in\n"
              "the first component (which corresponds to the normal direction) is valid!");

        if ((*onoff)[0])  // Do we have a load in normal direction?
        {
          if (func) functnum = (*func)[0];
          if (functnum > 0)
          {
            // evaluate function at current gauss point
            functfac = Global::Problem::instance()
                           ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                           .evaluate(coordgpref, time, idim);
            if (fldparatimint_->is_new_ost_implementation())
              functfacn = Global::Problem::instance()
                              ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                              .evaluate(coordgpref, time - fldparatimint_->dt(), idim);
          }
          else
          {
            functfac = 1.0;
            functfacn = 1.0;
          }

          const double valfac = (*val)[0] * fac_time_dens * functfac;
          for (int inode = 0; inode < bdrynen_; ++inode)
          {
            elevec1_epetra[inode * numdofpernode_ + idim] +=
                funct_(inode) * valfac * (-unitnormal_(idim));

            if (fldparatimint_->is_new_ost_implementation())
            {
              const double valfacn = (*val)[0] * fac_time_densn * functfacn;
              elevec1_epetra[inode * numdofpernode_ + idim] +=
                  funct_(inode) * valfacn * (-unitnormal_(idim));
            }
          }  // end is_new_ost_implementation
        }    // if (*onoff)
      }
      else
        FOUR_C_THROW(
            "The type '%s' is not supported in the fluid neumann condition!", type->c_str());

    }  // for(int idim=0; idim<(nsd_); ++idim)
  }

  return 0;
}



/*----------------------------------------------------------------------*
 | compute additional term at Neumann inflow boundary          vg 01/11 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::neumann_inflow(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseVector& elevec1)
{
  if (fldparatimint_->is_new_ost_implementation())
  {
    FOUR_C_THROW("NEUMANN INFLOW IS NOT IMPLEMENTED FOR NEW OST AS OF YET!");
  }  // end is_new_ost_implementation

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldparatimint_->time_fac();

  // get timefactor for right-hand side
  // One-step-Theta:            timefacrhs = theta*dt
  // BDF2:                      timefacrhs = 2/3 * dt
  // af-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // np-genalpha:               timefacrhs = (1/alpha_M) * gamma * dt
  // genalpha:                  timefacrhs = 1.0
  double timefacrhs = fldparatimint_->time_fac_rhs();

  // check ALE status
  const bool isale = ele->parent_element()->is_ale();

  // set flag for type of linearization to default value (fixed-point-like)
  bool is_newton = fldpara_->is_newton();

  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get global node coordinates for nsd_-dimensional domain
  // (nsd_: number of spatial dimensions of FluidBoundary element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;
  if (isale)
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < (nsd_); ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // get velocity and scalar vector at time n+alpha_F/n+1
  Teuchos::RCP<const Epetra_Vector> velaf = discretization.get_state("velaf");
  Teuchos::RCP<const Epetra_Vector> scaaf = discretization.get_state("scaaf");
  if (velaf == Teuchos::null or scaaf == Teuchos::null)
    FOUR_C_THROW("Cannot get state vector 'velaf' and/or 'scaaf'");

  // extract local values from global vector
  std::vector<double> myvelaf(lm.size());
  std::vector<double> myscaaf(lm.size());
  Core::FE::extract_my_values(*velaf, myvelaf, lm);
  Core::FE::extract_my_values(*scaaf, myscaaf, lm);

  // create Epetra objects for scalar array and velocities
  Core::LinAlg::Matrix<nsd_, bdrynen_> evelaf(true);
  Core::LinAlg::Matrix<bdrynen_, 1> epreaf(true);
  Core::LinAlg::Matrix<bdrynen_, 1> escaaf(true);

  // insert velocity, pressure and scalar into element array
  for (int inode = 0; inode < bdrynen_; ++inode)
  {
    for (int idim = 0; idim < (nsd_); ++idim)
    {
      evelaf(idim, inode) = myvelaf[idim + (inode * numdofpernode_)];
    }
    epreaf(inode) = myvelaf[nsd_ + inode * numdofpernode_];
    escaaf(inode) = myscaaf[(nsd_) + (inode * numdofpernode_)];
  }

  // get thermodynamic pressure at n+1/n+alpha_F
  const double thermpressaf = params.get<double>("thermpress at n+alpha_F/n+1", 1.0);

  // --------------------------------------------------
  // nurbs-specific stuff
  // --------------------------------------------------
  // normal vector multiplied by normalfac for nurbs
  double normalfac = 0.0;
  std::vector<Core::LinAlg::SerialDenseVector> mypknots(nsd_);
  std::vector<Core::LinAlg::SerialDenseVector> myknots(bdrynsd_);
  Core::LinAlg::SerialDenseVector weights(bdrynen_);

  // get knotvectors for parent element and surface element as well as weights
  // for isogeometric elements
  if (Core::FE::is_nurbs<distype>)
  {
    bool zero_size =
        Core::FE::Nurbs::get_knot_vector_and_weights_for_nurbs_boundary(ele, ele->surface_number(),
            ele->parent_element()->id(), discretization, mypknots, myknots, weights, normalfac);
    if (zero_size)
    {
      return;
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, &myknots, &weights, Core::FE::is_nurbs<distype>);

    // normal vector scaled by special factor in case of nurbs
    if (Core::FE::is_nurbs<distype>) unitnormal_.scale(normalfac);

    // compute velocity vector and normal velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double normvel = 0.0;
    velint_.multiply(evelaf, funct_);
    normvel = velint_.dot(unitnormal_);

    // check normal velocity -> further computation only required for
    // negative normal velocity, that is, inflow at this Neumann boundary
    if (normvel < -0.0001)
    {
      // get the required material information
      Teuchos::RCP<Core::Mat::Material> material = ele->parent_element()->material();

      // get density
      // (evaluation always at integration point, in contrast to parent element)
      get_density(material, escaaf, thermpressaf, epreaf);

      // extended integration factors for left- and right-hand side, respectively
      const double lhsfac = densaf_ * normvel * timefac * fac_;
      const double rhsfac = densaf_ * normvel * timefacrhs * fac_;

      // compute matrix contribution (fill diagonal elements)
      /*
              /                        \
             |                          |
           - |  v , rho * Du ( u o n )  |
             |                          |
              \                        /
      */
      for (int idim = 0; idim < nsd_; ++idim)  // loop over dimensions
      {
        for (int vi = 0; vi < bdrynen_; ++vi)  // loop over rows
        {
          const double vlhs = lhsfac * funct_(vi);

          const int fvi = numdofpernode_ * vi + idim;

          for (int ui = 0; ui < bdrynen_; ++ui)  // loop over columns
          {
            const int fui = numdofpernode_ * ui + idim;

            elemat1(fvi, fui) -= vlhs * funct_(ui);
          }  // end loop over columns
        }    // end loop over rows
      }      // end loop over dimensions

      // compute additional matrix contribution for Newton linearization
      if (is_newton)
      {
        // integration factor
        const double lhsnewtonfac = densaf_ * timefac * fac_;

        // dyadic product of unit normal vector and velocity vector
        Core::LinAlg::Matrix<nsd_, nsd_> n_x_u(true);
        n_x_u.multiply_nt(velint_, unitnormal_);

        /*
                /                        \
               |                          |
             - |  v , rho * u ( Du o n )  |
               |                          |
                \                        /

               rho * v_i * u_i * Du_j * n_j

        */
        for (int vi = 0; vi < bdrynen_; ++vi)  // loop rows
        {
          const double dens_dt_v = lhsnewtonfac * funct_(vi);

          for (int idimrow = 0; idimrow < nsd_; ++idimrow)  // loop row dim.
          {
            const int fvi = numdofpernode_ * vi + idimrow;

            for (int ui = 0; ui < bdrynen_; ++ui)  // loop columns
            {
              const double dens_dt_v_Du = dens_dt_v * funct_(ui);

              for (int idimcol = 0; idimcol < nsd_; ++idimcol)  // loop column dim.
              {
                const int fui = numdofpernode_ * ui + idimcol;

                elemat1(fvi, fui) -= dens_dt_v_Du * n_x_u(idimrow, idimcol);
              }  // end loop row dimensions
            }    // end loop rows
          }      // end loop column dimensions
        }        // end loop columns
      }          // end of Newton loop

      // compute rhs contribution
      Core::LinAlg::Matrix<nsd_, 1> vrhs(velint_, false);
      vrhs.scale(rhsfac);

      for (int vi = 0; vi < bdrynen_; ++vi)  // loop over rows
      {
        for (int idim = 0; idim < nsd_; ++idim)  // loop over dimensions
        {
          const int fvi = numdofpernode_ * vi + idim;

          elevec1(fvi) += funct_(vi) * vrhs(idim);
        }  // end loop over dimensions
      }    // end loop over rows
    }
  }

  return;
}  // Discret::ELEMENTS::FluidSurface::neumann_inflow


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)            gjb 07/07|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::integrate_shape_function(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, const std::vector<double>& edispnp)
{
  // get status of Ale
  const bool isale = ele->parent_element()->is_ale();

  // get Gaussrule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  if (isale)
  {
    FOUR_C_ASSERT(edispnp.size() != 0, "paranoid");

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < (nsd_); ++idim)
      {
        xyze_(idim, inode) += edispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // is not activated here Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < (nsd_); idim++)
      {
        elevec1(inode * numdofpernode_ + idim) += funct_(inode) * fac_;
      }
    }

  } /* end of loop over integration points gpid */


  return;
}  // Discret::ELEMENTS::FluidSurface::integrate_shape_function


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::element_mean_curvature(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, const std::vector<double>& edispnp,
    std::vector<double>& enormals)
{
  // get status of Ale
  const bool isale = ele->parent_element()->is_ale();

  // get Gauss rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // node normals &
  Core::LinAlg::Matrix<nsd_, bdrynen_> norm_elem(true);
  Core::LinAlg::Matrix<bdrynsd_, nsd_> dxyzdrs(true);

  // coordinates of current node in reference coordinates
  Core::LinAlg::Matrix<bdrynsd_, 1> xsi_node(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  if (isale)
  {
    FOUR_C_ASSERT(edispnp.size() != 0, "paranoid");

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += edispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // set normal vectors to length = 1.0
  // normal vector is coming from outside
  for (int inode = 0; inode < bdrynen_; ++inode)
  {
    // double length = 0.0;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      norm_elem(idim, inode) = enormals[numdofpernode_ * inode + idim];
    }
  }
  // compute normalized normal vector
  norm_elem.scale(1 / norm_elem.norm2());

  // get local node coordinates of the element
  // function gives back a matrix with the local node coordinates of the element (nsd_,bdrynen_)
  // the function gives back an Core::LinAlg::SerialDenseMatrix!!!
  Core::LinAlg::SerialDenseMatrix xsi_ele =
      Core::FE::get_ele_node_numbering_nodes_paramspace(distype);

  // ============================== loop over nodes ==========================
  for (int inode = 0; inode < bdrynen_; ++inode)
  {
    // the local node coordinates matrix is split to a vector containing the local coordinates of
    // the actual node
    for (int idim = 0; idim < bdrynsd_; idim++)
    {
      xsi_node(idim) = xsi_ele(idim, inode);
    }

    // get shape derivatives at this node
    // shape_function_2D_deriv1(deriv_, e0, e1, distype);
    Core::FE::shape_function<distype>(xsi_node, funct_);

    // the metric tensor and its determinant
    // Core::LinAlg::SerialDenseMatrix      metrictensor(nsd_,nsd_);
    Core::LinAlg::Matrix<bdrynsd_, bdrynsd_> metrictensor(true);

    // Addionally, compute metric tensor
    Core::FE::compute_metric_tensor_for_boundary_ele<distype>(xyze_, deriv_, metrictensor, drs_);

    dxyzdrs.multiply_nt(deriv_, xyze_);

    // calculate mean curvature H at node.
    double H = 0.0;
    Core::LinAlg::Matrix<bdrynsd_, nsd_> dn123drs(0.0);

    dn123drs.multiply_nt(deriv_, norm_elem);

    // Acc. to Bronstein ..."mittlere Kruemmung":
    // calculation of the mean curvature for a surface element
    if (bdrynsd_ == 2)
    {
      double L = 0.0, twoM = 0.0, N = 0.0;
      for (int i = 0; i < 3; i++)
      {
        L += (-1.0) * dxyzdrs(0, i) * dn123drs(0, i);
        twoM += (-1.0) * dxyzdrs(0, i) * dn123drs(1, i) - dxyzdrs(1, i) * dn123drs(0, i);
        N += (-1.0) * dxyzdrs(1, i) * dn123drs(1, i);
      }
      // mean curvature: H = 0.5*(k_1+k_2)
      H = 0.5 * (metrictensor(0, 0) * N - twoM * metrictensor(0, 1) + metrictensor(1, 1) * L) /
          (drs_ * drs_);
    }
    else
      FOUR_C_THROW(
          "Calcualtion of the mean curvature is only implemented for a 2D surface element");


    // get the number of elements adjacent to this node. Find out how many
    // will contribute to the interpolated mean curvature value.
    int contr_elements = 0;
    Core::Nodes::Node* thisNode = (ele->nodes())[inode];
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (thisNode == nullptr) FOUR_C_THROW("No node!\n");
#endif
    int NumElement = thisNode->num_element();
    Core::Elements::Element** ElementsPtr = thisNode->elements();

    // loop over adjacent Fluid elements
    for (int ele = 0; ele < NumElement; ele++)
    {
      Core::Elements::Element* Element = ElementsPtr[ele];

      // get surfaces
      std::vector<Teuchos::RCP<Core::Elements::Element>> surfaces = Element->surfaces();

      // loop over surfaces: how many free surfaces with this node on it?
      for (unsigned int surf = 0; surf < surfaces.size(); ++surf)
      {
        Teuchos::RCP<Core::Elements::Element> surface = surfaces[surf];
        Core::Nodes::Node** NodesPtr = surface->nodes();
        int numfsnodes = 0;
        bool hasthisnode = false;

        for (int surfnode = 0; surfnode < surface->num_node(); ++surfnode)
        {
          Core::Nodes::Node* checkNode = NodesPtr[surfnode];
          // check whether a free surface condition is active on this node
          if (checkNode->get_condition("FREESURFCoupling") != nullptr)
          {
            numfsnodes++;
          }
          if (checkNode->id() == thisNode->id())
          {
            hasthisnode = true;
          }
        }

        if (numfsnodes == surface->num_node() and hasthisnode)
        {
          // this is a free surface adjacent to this node.
          contr_elements++;
        }
      }
    }
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (!contr_elements) FOUR_C_THROW("No contributing elements found!\n");
#endif

    for (int idim = 0; idim < nsd_; ++idim)
    {
      elevec1[inode * numdofpernode_ + idim] = H / contr_elements;
    }
    elevec1[inode * numdofpernode_ + (numdofpernode_ - 1)] = 0.0;
  }  // END: loop over nodes

}  // Discret::ELEMENTS::FluidSurface::element_mean_curvature



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::element_surface_tension(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, const std::vector<double>& edispnp,
    std::vector<double>& enormals, std::vector<double>& ecurvature)
// Attention: mynormals and mycurvature are not used in the function
{
  // get status of Ale
  const bool isale = ele->parent_element()->is_ale();

  // get timefactor for left-hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // af-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // np-genalpha: timefac = (alpha_F/alpha_M) * gamma * dt
  // genalpha:    timefac =  alpha_F * gamma * dt
  const double timefac = fldparatimint_->time_fac();

  // isotropic and isothermal surface tension coefficient
  double SFgamma = 0.0;
  // get material data
  Teuchos::RCP<Core::Mat::Material> mat = ele->parent_element()->material();
  if (mat == Teuchos::null)
    FOUR_C_THROW("no mat from parent!");
  else if (mat->material_type() == Core::Materials::m_fluid)
  {
    const Mat::NewtonianFluid* actmat = static_cast<const Mat::NewtonianFluid*>(mat.get());
    SFgamma = actmat->gamma();
  }
  else
    FOUR_C_THROW("Newtonian fluid material expected but got type %d", mat->material_type());

  // get Gauss rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  if (isale)
  {
    FOUR_C_ASSERT(edispnp.size() != 0, "paranoid");

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += edispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // fac multiplied by the timefac
    const double fac_timefac = fac_ * timefac;

    // Compute dxyzdrs
    Core::LinAlg::Matrix<bdrynsd_, nsd_> dxyzdrs(true);
    dxyzdrs.multiply_nt(deriv_, xyze_);

    if (bdrynsd_ == 2)
    {
      double abs_dxyzdr = 0.0;
      double abs_dxyzds = 0.0;
      double pointproduct_rs = 0.0;

      for (int dim = 0; dim < 3; dim++)
      {
        abs_dxyzdr += dxyzdrs(0, dim) * dxyzdrs(0, dim);
        abs_dxyzds += dxyzdrs(1, dim) * dxyzdrs(1, dim);
        pointproduct_rs += dxyzdrs(0, dim) * dxyzdrs(1, dim);
      }
      abs_dxyzdr = sqrt(abs_dxyzdr);
      abs_dxyzds = sqrt(abs_dxyzds);

      for (int node = 0; node < bdrynen_; ++node)
      {
        for (int dim = 0; dim < 3; dim++)
        {
          // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
          // function) on Gamma_FS
          // See Saksono eq. (26)
          // discretized as surface gradient * ( Shapefunction-Matrix
          // transformed )

          // This uses a surface_gradient extracted from gauss general
          // formula for 2H...
          // this gives convincing results with TET elements, but HEX
          // elements seem more difficult -> due to edge problems?
          // too many nonlinear iterations
          elevec1[node * numdofpernode_ + dim] +=
              SFgamma * (-1.0) /
              (drs_ * drs_  //= abs_dxyzdr * abs_dxyzdr * abs_dxyzds * abs_dxyzds - pointproduct_rs
                            //* pointproduct_rs
                  ) *
              (abs_dxyzds * abs_dxyzds * deriv_(0, node) * dxyzdrs(0, dim) -
                  pointproduct_rs * deriv_(0, node) * dxyzdrs(1, dim) -
                  pointproduct_rs * deriv_(1, node) * dxyzdrs(0, dim) +
                  abs_dxyzdr * abs_dxyzdr * deriv_(1, node) * dxyzdrs(1, dim)) *
              fac_timefac;
        }
        elevec1[node * numdofpernode_ + 3] = 0.0;
      }
    }  // end if (nsd_==2)
    else if (bdrynsd_ == 1)
    {
      for (int inode = 0; inode < bdrynen_; ++inode)
      {
        for (int idim = 0; idim < 2; idim++)
        {
          // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
          // function) on Gamma_FS
          // See Saksono eq. (26)
          // discretized as surface gradient * ( Shapefunction-Matrix
          // transformed )
          // 2D: See Slikkerveer ep. (17)
          elevec1[inode * numdofpernode_ + idim] +=
              SFgamma / drs_ / drs_ * (-1.0) * deriv_(0, inode) * dxyzdrs(0, idim) * fac_timefac;
        }
      }
    }  // end if else (nsd_=1)
    else
      FOUR_C_THROW("There are no 3D boundary elements implemented");
  } /* end of loop over integration points gpid */
}  // Discret::ELEMENTS::FluidSurface::element_surface_tension

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::area_calculation(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm)
{
  //------------------------------------------------------------------
  // get and set density and viscosity (still required for following routines:
  // FluidImpedanceBc/FluidVolumetricSurfaceFlowBc/FluidCouplingBc::Area)
  //------------------------------------------------------------------
  Teuchos::RCP<Core::Mat::Material> mat = ele->parent_element()->material();
  if (mat->material_type() == Core::Materials::m_fluid)
  {
    const Mat::NewtonianFluid* actmat = static_cast<const Mat::NewtonianFluid*>(mat.get());
    densaf_ = actmat->density();
    visc_ = actmat->viscosity();
  }
  else if (mat->material_type() == Core::Materials::m_carreauyasuda)
  {
    const Mat::CarreauYasuda* actmat = static_cast<const Mat::CarreauYasuda*>(mat.get());
    densaf_ = actmat->density();

    const double nu_inf = actmat->nu_inf();  // parameter for infinite-shear viscosity

    // dynamic viscosity = kinematic viscosity * density
    visc_ = nu_inf * densaf_;
  }

  params.set<double>("density", densaf_);
  params.set<double>("viscosity", visc_);
  //------------------------------------------------------------------
  // end of get and set density and viscosity
  //------------------------------------------------------------------

  //------------------------------------------------------------------
  // start of actual area calculation
  //------------------------------------------------------------------
  // get node coordinates (nsd_: dimension of boundary element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;
  if (ele->parent_element()->is_ale())
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // get initial value for area
  double area = params.get<double>("area");

  // get Gauss rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // add to area integral
    area += fac_;
  }

  // set final value for area
  params.set<double>("area", area);
  //------------------------------------------------------------------
  // end of actual area calculation
  //------------------------------------------------------------------

}  // Discret::ELEMENTS::FluidSurface::area_calculation


/*----------------------------------------------------------------------*
 |                                                       ismail 04/2010 |
 |                                                           vg 06/2013 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::pressure_boundary_integral(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm)
{
  // extract pressure values from global velocity/pressure vector
  // renamed to "velaf" to be consistent in fluidimplicitintegration.cpp (krank 12/13)
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.get_state("velaf");
  if (velnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
  Core::FE::extract_my_values(*velnp, myvelnp, lm);

  Core::LinAlg::Matrix<1, bdrynen_> eprenp(true);
  for (int inode = 0; inode < bdrynen_; inode++)
  {
    eprenp(inode) = myvelnp[nsd_ + inode * numdofpernode_];
  }

  // get node coordinates (nsd_: dimension of boundary element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // add potential ALE displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;
  if (ele->parent_element()->is_ale())
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // get initial value for pressure boundary integral
  double press_int = params.get<double>("pressure boundary integral");

  // get Gauss rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // add to pressure boundary integral
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      press_int += funct_(inode) * eprenp(inode) * fac_;
    }
  }

  // set final value for pressure boundary integral
  params.set<double>("pressure boundary integral", press_int);

}  // Discret::ELEMENTS::FluidSurface::pressure_boundary_integral


/*----------------------------------------------------------------------*
 |                                                        ismail 10/2010|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::center_of_mass_calculation(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm)
{
  //------------------------------------------------------------------
  // This calculates the integrated the pressure from the
  // the actual pressure values
  //------------------------------------------------------------------

  // get integration rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  // Core::Geo::fill_initial_position_array<distype,nsd_,Core::LinAlg::SerialDenseMatrix>(ele,xyze_);
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->parent_element()->is_ale())
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // first evaluate the area of the surface element
  params.set<double>("area", 0.0);
  this->area_calculation(ele, params, discretization, lm);

  // get the surface element area
  const double elem_area = params.get<double>("area");

  Core::LinAlg::Matrix<(nsd_), 1> xyzGe(true);

  for (int i = 0; i < nsd_; i++)
  {
    // const IntegrationPoints2D  intpoints(gaussrule);
    for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
    {
      // Computation of the integration factor & shape function at the Gauss point & derivative of
      // the shape function at the Gauss point Computation of the unit normal vector at the Gauss
      // points Computation of nurb specific stuff is not activated here
      Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
          xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

      // global coordinates of gausspoint
      Core::LinAlg::Matrix<(nsd_), 1> coordgp(true);

      // determine coordinates of current Gauss point
      coordgp.multiply(xyze_, funct_);

      // Compute elment center of gravity
      xyzGe(i) += intpoints.ip().qwgt[gpid] * coordgp(i) * drs_;

    }  // end Gauss loop
    xyzGe(i) /= elem_area;
  }

  // Get the center of mass of the already calculate surface elements
  Teuchos::RCP<std::vector<double>> xyzG =
      params.get<Teuchos::RCP<std::vector<double>>>("center of mass");

  Teuchos::RCP<std::vector<double>> normal =
      params.get<Teuchos::RCP<std::vector<double>>>("normal");

  // Get the area of the of the already calculate surface elements
  double area = params.get<double>("total area");

  for (int i = 0; i < nsd_; i++)
  {
    (*xyzG)[i] = ((*xyzG)[i] * area + xyzGe(i) * elem_area) / (area + elem_area);
    (*normal)[i] = ((*normal)[i] * area + unitnormal_(i) * elem_area) / (area + elem_area);
  }

  // set new center of mass
  params.set("total area", area + elem_area);

}  // Discret::ELEMENTS::FluidSurface::center_of_mass_calculation



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::compute_flow_rate(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1)
{
  // get integration rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // extract local values from the global vectors
  // renamed to "velaf" to be consistent in fluidimplicitintegration.cpp (krank 12/13)
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.get_state("velaf");

  if (velnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
  Core::FE::extract_my_values(*velnp, myvelnp, lm);

  // allocate velocity vector
  Core::LinAlg::Matrix<nsd_, bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < bdrynen_; inode++)
  {
    for (int idim = 0; idim < nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * numdofpernode_)];
    }
  }

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  // Core::Geo::fill_initial_position_array<distype,nsd_,Core::LinAlg::SerialDenseMatrix>(ele,xyze_);
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->parent_element()->is_ale())
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // compute flowrate at gauss point
    velint_.multiply(evelnp, funct_);

    // flowrate = uint o normal
    const double flowrate = velint_.dot(unitnormal_);

    // store flowrate at first dof of each node
    // use negative value so that inflow is positiv
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      // see "A better consistency for low order stabilized finite element methods"
      // Jansen, Collis, Whiting, Shakib
      //
      // Here the principle is used to bring the flow rate to the outside world!!
      //
      // funct_ *  velint * n * fac
      //   |      |________________|
      //   |              |
      //   |         flow rate * fac  -> integral over Gamma
      //   |
      // flow rate is distributed to the single nodes of the element
      // = flow rate per node
      //
      // adding up all nodes (ghost elements are handled by the assembling strategy)
      // -> total flow rate at the desired boundary
      //
      // it can be interpreted as a rhs term
      //
      //  ( v , u o n)
      //               Gamma
      //
      elevec1[inode * numdofpernode_] += funct_(inode) * fac_ * flowrate;

      // alternative way:
      //
      //  velint * n * fac
      // |________________|
      //         |
      //    flow rate * fac  -> integral over Gamma
      //     = flow rate per element
      //
      //  adding up all elements (be aware of ghost elements!!)
      //  -> total flow rate at the desired boundary
      //     (is identical to the total flow rate computed above)
    }
  }
}  // Discret::ELEMENTS::FluidSurface::compute_flow_rate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::flow_rate_deriv(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // This function is only implemented for 3D
  if (bdrynsd_ != 2) FOUR_C_THROW("flow_rate_deriv is only implemented for 3D!");

  // get status of Ale
  const bool isale = ele->parent_element()->is_ale();

  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> edispnp;

  if (isale)
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'dispnp'");
    edispnp.resize(lm.size());
    Core::FE::extract_my_values(*dispnp, edispnp, lm);
  }

  // get integration rule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // order of accuracy of grid velocity determination
  const Teuchos::ParameterList& fdyn = Global::Problem::instance()->fluid_dynamic_params();
  const int gridvel = Core::UTILS::integral_value<Inpar::FLUID::Gridvel>(fdyn, "GRIDVEL");

  // normal vector
  Core::LinAlg::Matrix<nsd_, 1> normal(true);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  if (isale)
  {
    FOUR_C_ASSERT(edispnp.size() != 0, "paranoid");

    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += edispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // get nodal velocities and pressures
  Teuchos::RCP<const Epetra_Vector> convelnp = discretization.get_state("convectivevel");

  if (convelnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'convectivevel'");

  // extract local values from the global vectors
  std::vector<double> myconvelnp(lm.size());
  Core::FE::extract_my_values(*convelnp, myconvelnp, lm);

  // allocate velocities vector
  Core::LinAlg::Matrix<nsd_, bdrynen_> evelnp(true);

  for (int inode = 0; inode < bdrynen_; ++inode)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      evelnp(idim, inode) = myconvelnp[(numdofpernode_ * inode) + idim];
    }
  }


  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
    *----------------------------------------------------------------------*/
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // is not activated here Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // The integration factor is not multiplied with drs
    // since it is the same as the scaling factor for the unit normal
    // Therefore it cancels out!!
    const double fac = intpoints.ip().qwgt[gpid];

    // dxyzdrs vector -> normal which is not normalized
    Core::LinAlg::Matrix<bdrynsd_, nsd_> dxyzdrs(0.0);
    dxyzdrs.multiply_nt(deriv_, xyze_);
    normal(0, 0) = dxyzdrs(0, 1) * dxyzdrs(1, 2) - dxyzdrs(0, 2) * dxyzdrs(1, 1);
    normal(1, 0) = dxyzdrs(0, 2) * dxyzdrs(1, 0) - dxyzdrs(0, 0) * dxyzdrs(1, 2);
    normal(2, 0) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(0, 1) * dxyzdrs(1, 0);

    //-------------------------------------------------------------------
    //  Q
    Core::LinAlg::Matrix<3, 1> u(true);
    for (int dim = 0; dim < 3; ++dim)
      for (int node = 0; node < bdrynen_; ++node) u(dim) += funct_(node) * evelnp(dim, node);

    for (int dim = 0; dim < 3; ++dim) elevec3[0] += u(dim) * normal(dim, 0) * fac;

    if (params.get<bool>("flowrateonly", false) == false)
    {
      //-------------------------------------------------------------------
      // dQ/du
      for (int node = 0; node < bdrynen_; ++node)
      {
        for (int dim = 0; dim < 3; ++dim)
          elevec1[node * numdofpernode_ + dim] += funct_(node) * normal(dim, 0) * fac;
        elevec1[node * numdofpernode_ + 3] = 0.0;
      }

      //-------------------------------------------------------------------
      // dQ/dd

      // determine derivatives of surface normals wrt mesh displacements
      Core::LinAlg::Matrix<3, bdrynen_ * 3> normalderiv(true);

      for (int node = 0; node < bdrynen_; ++node)
      {
        normalderiv(0, 3 * node) = 0.;
        normalderiv(0, 3 * node + 1) =
            deriv_(0, node) * dxyzdrs(1, 2) - deriv_(1, node) * dxyzdrs(0, 2);
        normalderiv(0, 3 * node + 2) =
            deriv_(1, node) * dxyzdrs(0, 1) - deriv_(0, node) * dxyzdrs(1, 1);

        normalderiv(1, 3 * node) =
            deriv_(1, node) * dxyzdrs(0, 2) - deriv_(0, node) * dxyzdrs(1, 2);
        normalderiv(1, 3 * node + 1) = 0.;
        normalderiv(1, 3 * node + 2) =
            deriv_(0, node) * dxyzdrs(1, 0) - deriv_(1, node) * dxyzdrs(0, 0);

        normalderiv(2, 3 * node) =
            deriv_(0, node) * dxyzdrs(1, 1) - deriv_(1, node) * dxyzdrs(0, 1);
        normalderiv(2, 3 * node + 1) =
            deriv_(1, node) * dxyzdrs(0, 0) - deriv_(0, node) * dxyzdrs(1, 0);
        normalderiv(2, 3 * node + 2) = 0.;
      }

      for (int node = 0; node < bdrynen_; ++node)
      {
        for (int dim = 0; dim < 3; ++dim)
          for (int iterdim = 0; iterdim < 3; ++iterdim)
            elevec2[node * numdofpernode_ + dim] +=
                u(iterdim) * normalderiv(iterdim, 3 * node + dim) * fac;
        elevec2[node * numdofpernode_ + 3] = 0.0;
      }

      // consideration of grid velocity
      if (isale)
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) FOUR_C_THROW("invalid time step size");

        if (gridvel == Inpar::FLUID::BE)  // BE time discretization
        {
          for (int node = 0; node < bdrynen_; ++node)
          {
            for (int dim = 0; dim < 3; ++dim)
              elevec2[node * numdofpernode_ + dim] -=
                  1.0 / dt * funct_(node) * normal(dim, 0) * fac;
          }
        }
        else
          FOUR_C_THROW(
              "flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(du dd)

      for (int unode = 0; unode < bdrynen_; ++unode)
      {
        for (int udim = 0; udim < numdofpernode_; ++udim)
        {
          for (int nnode = 0; nnode < bdrynen_; ++nnode)
          {
            for (int ndim = 0; ndim < numdofpernode_; ++ndim)
            {
              if (udim == 3 or ndim == 3)
                elemat1(unode * numdofpernode_ + udim, nnode * numdofpernode_ + ndim) = 0.0;
              else
                elemat1(unode * numdofpernode_ + udim, nnode * numdofpernode_ + ndim) =
                    funct_(unode) * normalderiv(udim, 3 * nnode + ndim) * fac;
            }
          }
        }
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(dd)^2

      // determine second derivatives of surface normals wrt mesh displacements
      std::vector<Core::LinAlg::Matrix<bdrynen_ * 3, bdrynen_ * 3>> normalderiv2(3);

      for (int node1 = 0; node1 < bdrynen_; ++node1)
      {
        for (int node2 = 0; node2 < bdrynen_; ++node2)
        {
          double temp = deriv_(0, node1) * deriv_(1, node2) - deriv_(1, node1) * deriv_(0, node2);

          normalderiv2[0](node1 * 3 + 1, node2 * 3 + 2) = temp;
          normalderiv2[0](node1 * 3 + 2, node2 * 3 + 1) = -temp;

          normalderiv2[1](node1 * 3, node2 * 3 + 2) = -temp;
          normalderiv2[1](node1 * 3 + 2, node2 * 3) = temp;

          normalderiv2[2](node1 * 3, node2 * 3 + 1) = temp;
          normalderiv2[2](node1 * 3 + 1, node2 * 3) = -temp;
        }
      }

      for (int node1 = 0; node1 < bdrynen_; ++node1)
      {
        for (int dim1 = 0; dim1 < numdofpernode_; ++dim1)
        {
          for (int node2 = 0; node2 < bdrynen_; ++node2)
          {
            for (int dim2 = 0; dim2 < numdofpernode_; ++dim2)
            {
              if (dim1 == 3 or dim2 == 3)
                elemat2(node1 * numdofpernode_ + dim1, node2 * numdofpernode_ + dim2) = 0.0;
              else
              {
                for (int iterdim = 0; iterdim < 3; ++iterdim)
                  elemat2(node1 * numdofpernode_ + dim1, node2 * numdofpernode_ + dim2) +=
                      u(iterdim) * normalderiv2[iterdim](node1 * 3 + dim1, node2 * 3 + dim2) * fac;
              }
            }
          }
        }
      }

      // consideration of grid velocity
      if (isale)
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) FOUR_C_THROW("invalid time step size");

        if (gridvel == Inpar::FLUID::BE)
        {
          for (int node1 = 0; node1 < bdrynen_; ++node1)
          {
            for (int dim1 = 0; dim1 < 3; ++dim1)
            {
              for (int node2 = 0; node2 < bdrynen_; ++node2)
              {
                for (int dim2 = 0; dim2 < 3; ++dim2)
                {
                  elemat2(node1 * numdofpernode_ + dim1, node2 * numdofpernode_ + dim2) -=
                      (1.0 / dt * funct_(node1) * normalderiv(dim1, 3 * node2 + dim2) +
                          1.0 / dt * funct_(node2) * normalderiv(dim2, 3 * node1 + dim1)) *
                      fac;
                }
              }
            }
          }
        }
        else
          FOUR_C_THROW(
              "flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
    }
  }
}  // Discret::ELEMENTS::FluidSurface::flow_rate_deriv


/*----------------------------------------------------------------------*
 |  Impedance related parameters on boundary elements          AC 03/08  |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::impedance_integration(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1)
{
  const double tfacrhs = fldparatimint_->time_fac_rhs();

  const double pressure = params.get<double>("WindkesselPressure");

  // get Gaussrule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  if (ele->parent_element()->is_ale())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // add traction in the inward normal direction with norm 'pressure'
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    const double fac_facrhs_pres = fac_ * tfacrhs * pressure;

    for (int inode = 0; inode < bdrynen_; ++inode)
      for (int idim = 0; idim < nsd_; ++idim)
        elevec1[inode * numdofpernode_ + idim] +=
            fac_facrhs_pres * funct_(inode) * (-unitnormal_(idim));
  }

  return;
}  // Discret::ELEMENTS::FluidSurface::impedance_integration


/*---------------------------------------------------------------------------*
 |  linearization of flux w.r.t velocities on boundary elements  Thon 10/15  |
 *---------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::d_qdu(Discret::ELEMENTS::FluidBoundary* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1)
{
  // get Gaussrule
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  if (ele->parent_element()->is_ale())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  // compute dQ/du were (dQ/du)_i= (\phi_i \dot n)_Gamma
  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    for (int node = 0; node < bdrynen_; ++node)
    {
      for (int dim = 0; dim < nsd_; ++dim)
      {
        elevec1[node * numdofpernode_ + dim] += funct_(node) * unitnormal_(dim, 0) * fac_;
      }
      elevec1[node * numdofpernode_ + nsd_] = 0.0;
    }
  }

  params.set<double>("tfaclhs", fldparatimint_->time_fac());
}


/*----------------------------------------------------------------------*
 |  get density                                                vg 06/13 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::get_density(
    Teuchos::RCP<const Core::Mat::Material> material,
    const Core::LinAlg::Matrix<bdrynen_, 1>& escaaf, const double thermpressaf,
    const Core::LinAlg::Matrix<bdrynen_, 1>& epreaf)
{
  // initially set density and density factor for Neumann boundary conditions to 1.0
  // (the latter only changed for low-Mach-number flow/combustion problems)
  // (This is due to the nature of the Neumann BC's for these problems as they are usually given as
  // h_N=\rho*h.) (Thus in the input file for these problems, h is set which is then scaled by \rho
  // here.)
  //
  //          See Gravemeier, Wall 2011 IJNMF example 4.1.2
  //
  densaf_ = 1.0;
  densfac_ = 1.0;

  if (material->material_type() == Core::Materials::m_fluid)
  {
    const Mat::NewtonianFluid* actmat = static_cast<const Mat::NewtonianFluid*>(material.get());

    // varying density
    if (fldpara_->physical_type() == Inpar::FLUID::varying_density) densaf_ = funct_.dot(escaaf);
    // Boussinesq approximation: Calculation of delta rho
    else if (fldpara_->physical_type() == Inpar::FLUID::boussinesq)
      FOUR_C_THROW("Boussinesq approximation not yet supported for boundary terms!");
    else
      densaf_ = actmat->density();
  }
  else if (material->material_type() == Core::Materials::m_carreauyasuda)
  {
    const Mat::CarreauYasuda* actmat = static_cast<const Mat::CarreauYasuda*>(material.get());

    densaf_ = actmat->density();
  }
  else if (material->material_type() == Core::Materials::m_modpowerlaw)
  {
    const Mat::ModPowerLaw* actmat = static_cast<const Mat::ModPowerLaw*>(material.get());

    densaf_ = actmat->density();
  }
  else if (material->material_type() == Core::Materials::m_herschelbulkley)
  {
    const Mat::HerschelBulkley* actmat = static_cast<const Mat::HerschelBulkley*>(material.get());

    densaf_ = actmat->density();
  }
  else if (material->material_type() == Core::Materials::m_sutherland)
  {
    const Mat::Sutherland* actmat = static_cast<const Mat::Sutherland*>(material.get());

    // compute temperature at n+alpha_F or n+1
    const double tempaf = funct_.dot(escaaf);

    // compute density at n+alpha_F or n+1 based on temperature
    // and thermodynamic pressure
    densaf_ = actmat->compute_density(tempaf, thermpressaf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->material_type() == Core::Materials::m_fluid_linear_density_viscosity)
  {
    const Mat::LinearDensityViscosity* actmat =
        static_cast<const Mat::LinearDensityViscosity*>(material.get());

    // compute pressure at n+alpha_F or n+1
    const double preaf = funct_.dot(epreaf);

    // compute density at n+alpha_F or n+1 based on pressure
    densaf_ = actmat->compute_density(preaf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->material_type() == Core::Materials::m_fluid_murnaghantait)
  {
    const Mat::MurnaghanTaitFluid* actmat =
        static_cast<const Mat::MurnaghanTaitFluid*>(material.get());

    // compute pressure at n+alpha_F or n+1
    const double preaf = funct_.dot(epreaf);

    // compute density at n+alpha_F or n+1 based on pressure
    densaf_ = actmat->compute_density(preaf);

    // set density factor for Neumann boundary conditions to density for present material
    densfac_ = densaf_;
  }
  else if (material->material_type() == Core::Materials::m_fluidporo)
  {
    const Mat::FluidPoro* actmat = static_cast<const Mat::FluidPoro*>(material.get());

    densaf_ = actmat->density();
  }
  else if (material->material_type() == Core::Materials::m_matlist)
  {
    densaf_ = 1.0;
    densfac_ = 1.0;

    // This is only necessary if the BC is dependant on the density!

  }  // end else if m_matlist
  else
    FOUR_C_THROW("Material type is not supported for density evaluation for boundary element!");

  //  // check whether there is zero or negative density
  if (densaf_ < 1e-15) FOUR_C_THROW("zero or negative density!");



  return;
}  // FluidBoundaryImpl::GetDensity

/*----------------------------------------------------------------------*
 |  Evaluating the velocity component of the traction      ismail 05/11 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::calc_traction_velocity_component(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1)
{
  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> velnp = discretization.get_state("velaf");

  if (velnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'velaf'");

  std::vector<double> myvelnp(lm.size());
  Core::FE::extract_my_values(*velnp, myvelnp, lm);

  // allocate velocity vector
  Core::LinAlg::Matrix<nsd_, bdrynen_> evelnp(true);

  // split velocity and pressure, insert into element arrays
  for (int inode = 0; inode < bdrynen_; inode++)
  {
    for (int idim = 0; idim < nsd_; idim++)
    {
      evelnp(idim, inode) = myvelnp[idim + (inode * numdofpernode_)];
    }
  }


  Teuchos::RCP<Epetra_Vector> cond_velocities =
      params.get<Teuchos::RCP<Epetra_Vector>>("condition velocities");
  Teuchos::RCP<Epetra_Map> cond_dofrowmap =
      params.get<Teuchos::RCP<Epetra_Map>>("condition dofrowmap");

  double density = 0.0;  // inverse density of my parent element

  // get material of volume element this surface belongs to
  Teuchos::RCP<Core::Mat::Material> mat = ele->parent_element()->material();

  if (mat->material_type() != Core::Materials::m_carreauyasuda &&
      mat->material_type() != Core::Materials::m_modpowerlaw &&
      mat->material_type() != Core::Materials::m_herschelbulkley &&
      mat->material_type() != Core::Materials::m_fluid)
    FOUR_C_THROW("Material law is not a fluid");

  if (mat->material_type() == Core::Materials::m_fluid)
  {
    const Mat::NewtonianFluid* actmat = static_cast<const Mat::NewtonianFluid*>(mat.get());
    density = actmat->density();
  }
  else if (mat->material_type() == Core::Materials::m_carreauyasuda)
  {
    const Mat::CarreauYasuda* actmat = static_cast<const Mat::CarreauYasuda*>(mat.get());
    density = actmat->density();
  }
  else if (mat->material_type() == Core::Materials::m_modpowerlaw)
  {
    const Mat::ModPowerLaw* actmat = static_cast<const Mat::ModPowerLaw*>(mat.get());
    density = actmat->density();
  }
  else if (mat->material_type() == Core::Materials::m_herschelbulkley)
  {
    const Mat::HerschelBulkley* actmat = static_cast<const Mat::HerschelBulkley*>(mat.get());
    density = actmat->density();
  }
  else
    FOUR_C_THROW("Fluid material expected but got type %d", mat->material_type());

  //-------------------------------------------------------------------
  // get the tractions velocity component
  //-------------------------------------------------------------------

  // get Gaussrule
  //  const Core::FE::IntPointsAndWeights<bdrynsd_>
  //  intpoints(Discret::ELEMENTS::DisTypeToGaussRuleForExactSol<distype>::rule);
  const Core::FE::IntPointsAndWeights<bdrynsd_> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of FluidBoundary
  // element!)
  Core::Geo::fill_initial_position_array<distype, nsd_, Core::LinAlg::Matrix<nsd_, bdrynen_>>(
      ele, xyze_);

  // Add the deformation of the ALE mesh to the nodes coordinates
  // displacements
  Teuchos::RCP<const Epetra_Vector> dispnp;
  std::vector<double> mydispnp;

  if (ele->parent_element()->is_ale())
  {
    dispnp = discretization.get_state("dispnp");
    if (dispnp != Teuchos::null)
    {
      mydispnp.resize(lm.size());
      Core::FE::extract_my_values(*dispnp, mydispnp, lm);
    }
    FOUR_C_ASSERT(mydispnp.size() != 0, "paranoid");
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        xyze_(idim, inode) += mydispnp[numdofpernode_ * inode + idim];
      }
    }
  }

  const double timefac = fldparatimint_->time_fac_rhs();

  for (int gpid = 0; gpid < intpoints.ip().nquad; gpid++)
  {
    // Computation of the integration factor & shape function at the Gauss point & derivative of the
    // shape function at the Gauss point Computation of the unit normal vector at the Gauss points
    // Computation of nurb specific stuff is not activated here
    Core::FE::eval_shape_func_at_bou_int_point<distype>(funct_, deriv_, fac_, unitnormal_, drs_,
        xsi_, xyze_, intpoints, gpid, nullptr, nullptr, Core::FE::is_nurbs<distype>);

    // Get the velocity value at the corresponding Gauss point.
    std::vector<double> vel_gps(nsd_, 0.0);
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        vel_gps[idim] += myvelnp[inode * numdofpernode_ + idim] * funct_(inode);
      }
    }

    // Evaluate the normal velocity at the corresponding Gauss point
    double n_vel = 0.0;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      n_vel += vel_gps[idim] * (unitnormal_(idim));
    }
    // loop over all node and add the corresponding effect of the Neumann-Inflow condition
    for (int inode = 0; inode < bdrynen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        // evaluate the value of the Un.U at the corresponding Gauss point
        const double uV = n_vel * vel_gps[idim] * density;
        const double fac_thsl_pres_inve = fac_ * timefac * uV;

        // remove the Neumann-inflow contribution only if the normal velocity is an inflow velocity
        // i.e n_vel < 0
        if (n_vel < 0.0)
        {
          elevec1[inode * numdofpernode_ + idim] -= fac_thsl_pres_inve * funct_(inode);
        }
      }
      //      double radius = sqrt(pow(xyze_(0,inode),2.0)+pow(xyze_(1,inode),2.0));
      //      cout<<"n_vel("<<n_vel<<") vel: "<<n_vel<<" rad: "<<radius<<endl;
    }
  }
  return;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::FluidBoundaryImpl<distype>::compute_neumann_uv_integral(
    Discret::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1)
{
}
// Discret::ELEMENTS::FluidSurface::compute_neumann_uv_integral

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::line3>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::nurbs2>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::nurbs3>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::nurbs4>;
template class Discret::ELEMENTS::FluidBoundaryImpl<Core::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
