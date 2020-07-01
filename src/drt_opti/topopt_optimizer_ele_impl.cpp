/*---------------------------------------------------------------------*/
/*! \file

\brief element routines of the topology optimization element

\level 3


*/
/*---------------------------------------------------------------------*/


#include "topopt_optimizer_ele_impl.H"
#include "topopt_optimizer_ele_parameter.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../headers/definitions.h"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptImplInterface* DRT::ELEMENTS::TopOptImplInterface::Impl(
    const DRT::Element* ele)
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
    {
      return TopOptImpl<DRT::Element::hex8>::Instance();
    }
    case DRT::Element::hex20:
    {
      return TopOptImpl<DRT::Element::hex20>::Instance();
    }
    case DRT::Element::hex27:
    {
      return TopOptImpl<DRT::Element::hex27>::Instance();
    }
    case DRT::Element::tet4:
    {
      return TopOptImpl<DRT::Element::tet4>::Instance();
    }
    case DRT::Element::tet10:
    {
      return TopOptImpl<DRT::Element::tet10>::Instance();
    }
    case DRT::Element::quad4:
    {
      return TopOptImpl<DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return TopOptImpl<DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return TopOptImpl<DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return TopOptImpl<DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return TopOptImpl<DRT::Element::tri6>::Instance();
    }
    default:
    {
      dserror("Element shape %s not activated. Just do it.",
          DRT::DistypeToString(ele->Shape()).c_str());
      break;
    }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype>* DRT::ELEMENTS::TopOptImpl<distype>::Instance(bool create)
{
  static TopOptImpl<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new TopOptImpl<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype>::TopOptImpl()
    : xyze_(true),
      funct_(true),
      deriv_(true),
      derxy_(true),
      eadjointpres_(true),
      fluidvelint_(true),
      adjointvelint_(true),
      fluidvelxy_(true),
      adjointvelxy_(true),
      adjointpresxy_(true),
      adjointconvvel_(true),
      fluidvelint_old_(true),
      tau_(true),
      poroint_(0.0),
      intpoints_(distype),
      xsi_(true),
      xjm_(true),
      xji_(true),
      fac_(0.0),
      is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat)
{
  return EvaluateValues(ele, params, optidis, mat, intpoints_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat,
    DRT::UTILS::GaussIntegration& intpoints)
{
  double& objective = *params.get<double*>("objective_value");
  Epetra_SerialDenseVector& constraints =
      *params.get<Epetra_SerialDenseVector*>("constraint_values");

  Teuchos::RCP<DRT::Discretization> fluiddis =
      params.get<Teuchos::RCP<DRT::Discretization>>("fluiddis");

  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>> fluidvels =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>>>("fluidvel");

  std::map<int, LINALG::Matrix<nsd_, nen_>> efluidvels;

  LINALG::Matrix<nsd_, nen_> efluidvel;

  // extract element data of all time steps from fluid solution
  std::vector<int> fluidlm;
  DRT::UTILS::DisBasedLocationVector(*fluiddis, *ele, fluidlm, nsd_ + 1);

  for (std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator i = fluidvels->begin();
       i != fluidvels->end(); i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis, fluidlm, &efluidvel, NULL, i->second);
    efluidvels.insert(std::pair<int, LINALG::Matrix<nsd_, nen_>>(i->first, efluidvel));
  }

  // evaluate nodal porosities
  LINALG::Matrix<nen_, 1> edens(true);
  if (optiparams_->DensType() == INPAR::TOPOPT::dens_node_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    for (int nn = 0; nn < nen_; ++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();
      edens(nn, 0) = (*topopt_density)[lid];
    }
  }
  else if (optiparams_->DensType() == INPAR::TOPOPT::dens_ele_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    int lid = ele->LID();
    for (int nn = 0; nn < nen_; ++nn)  // set all values equal to hack a constant element porosity
                                       // on element level -> inefficient, but not relevant
      edens(nn, 0) = (*topopt_density)[lid];
  }
  else
    dserror("not implemented type of density function");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, xyze_);

  Values(ele->Id(), efluidvels, edens, objective, constraints, mat, intpoints);

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Values(const int eid,
    std::map<int, LINALG::Matrix<nsd_, nen_>>& efluidvel, LINALG::Matrix<nen_, 1>& edens,
    double& objective, Epetra_SerialDenseVector& constraints, Teuchos::RCP<MAT::Material> mat,
    DRT::UTILS::GaussIntegration& intpoints)
{
  double& densint = constraints[0];  // integrated density

  LINALG::Matrix<nsd_, nen_> evel(true);

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    double value = 0.0;

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad, eid);

    EvalPorosityAtIntPoint(edens);

    // volume constraint
    densint += fac_ * (edens.Dot(funct_) - optiparams_->VolBd());

    switch (optiparams_->TimeIntScheme())
    {
      case INPAR::FLUID::timeint_stationary:
      case INPAR::FLUID::timeint_one_step_theta:  // handle these two cases together
      {
        for (int timestep = 0; timestep <= optiparams_->NumTimesteps(); timestep++)
        {
          if (optiparams_->IsStationary()) timestep = 1;  // just one stationary time step 1

          evel = efluidvel.find(timestep)->second;

          fluidvelint_.Multiply(evel, funct_);
          fluidvelxy_.MultiplyNT(evel, derxy_);

          // weighting of the timesteps
          // if stationary, we have the second case with theta = 1, so all is ok
          double objfac = 0.0;
          if (timestep == 0)  // first time step -> factor 1-theta (old sol at first time step)
            objfac = 1.0 - optiparams_->ThetaObj();
          else if (timestep == optiparams_->NumTimesteps())  // last time step -> factor theta (new
                                                             // sol at last time step)
            objfac = optiparams_->ThetaObj();
          else  // all other time steps -> factor 1-theta as old sol, factor theta as new sol ->
                // overall factor 1
            objfac = 1.0;

          if (optiparams_->OptiCase() == INPAR::TOPOPT::optitest_workflow_without_fluiddata)
            value -= objfac * edens.Dot(funct_) * edens.Dot(funct_);
          else  // standard case
          {
            if (optiparams_->ObjDissipationTerm() == INPAR::TOPOPT::obj_diss_yes)
              value += objfac * poroint_ * fluidvelint_.Dot(fluidvelint_);

            for (int idim = 0; idim < nsd_; idim++)
            {
              for (int jdim = 0; jdim < nsd_; jdim++)
              {
                value += objfac * optiparams_->Viscosity() * fluidvelxy_(idim, jdim) *
                         (fluidvelxy_(idim, jdim) + fluidvelxy_(jdim, idim));
              }
            }
          }
        }
      }
      break;
      default:
      {
        dserror("unknown time integration scheme while evaluating objective gradient");
        break;
      }
    }
    objective += optiparams_->Dt() * fac_ * optiparams_->ObjDissipationFac() *
                 value;  // zero if no dissipation in obj-fcn
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradients(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat)
{
  return EvaluateGradients(ele, params, optidis, mat, intpoints_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradients(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat,
    DRT::UTILS::GaussIntegration& intpoints)
{
  Epetra_SerialDenseVector obj_deriv;
  Epetra_SerialDenseVector constr_deriv;  // volume constraint

  if (optiparams_->DensType() == INPAR::TOPOPT::dens_node_based)
  {
    obj_deriv.Size(nen_);
    constr_deriv.Size(nen_);
  }
  else if (optiparams_->DensType() == INPAR::TOPOPT::dens_ele_based)
  {
    obj_deriv.Size(1);
    constr_deriv.Size(1);
  }
  else
    dserror("not implemented type of density function");

  Teuchos::RCP<DRT::Discretization> fluiddis =
      params.get<Teuchos::RCP<DRT::Discretization>>("fluiddis");

  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>> fluidvels =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>>>("fluidvel");
  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>> adjointvels =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>>>("adjointvel");

  std::map<int, LINALG::Matrix<nsd_, nen_>> efluidvels;
  std::map<int, LINALG::Matrix<nsd_, nen_>> eadjointvels;
  std::map<int, LINALG::Matrix<nen_, 1>> eadjointpress;

  LINALG::Matrix<nsd_, nen_> efluidvel;
  LINALG::Matrix<nsd_, nen_> eadjointvel;
  LINALG::Matrix<nen_, 1> eadjointpres;

  // extract element data of all time steps from fluid solution
  std::vector<int> fluidlm;
  DRT::UTILS::DisBasedLocationVector(*fluiddis, *ele, fluidlm, nsd_ + 1);

  for (std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator i = fluidvels->begin();
       i != fluidvels->end(); i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis, fluidlm, &efluidvel, NULL, i->second);
    efluidvels.insert(std::pair<int, LINALG::Matrix<nsd_, nen_>>(i->first, efluidvel));
  }

  for (std::map<int, Teuchos::RCP<Epetra_Vector>>::iterator i = adjointvels->begin();
       i != adjointvels->end(); i++)
  {
    ExtractValuesFromGlobalVector(*fluiddis, fluidlm, &eadjointvel, &eadjointpres, i->second);
    eadjointvels.insert(std::pair<int, LINALG::Matrix<nsd_, nen_>>(i->first, eadjointvel));
    eadjointpress.insert(std::pair<int, LINALG::Matrix<nen_, 1>>(i->first, eadjointpres));
  }

  // evaluate nodal porosities
  LINALG::Matrix<nen_, 1> edens(true);
  if (optiparams_->DensType() == INPAR::TOPOPT::dens_node_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    for (int nn = 0; nn < nen_; ++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();
      edens(nn, 0) = (*topopt_density)[lid];
    }
  }
  else if (optiparams_->DensType() == INPAR::TOPOPT::dens_ele_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    int lid = ele->LID();
    for (int nn = 0; nn < nen_; ++nn)  // set all values equal to hack a constant element porosity
                                       // on element level -> inefficient, but not relevant
      edens(nn, 0) = (*topopt_density)[lid];
  }
  else
    dserror("not implemented type of density function");

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, xyze_);


  Gradients(ele->Id(), efluidvels, eadjointvels, eadjointpress, edens, obj_deriv, constr_deriv, mat,
      intpoints);

  // In advance, the map of the constraints and objectives derivation is unknown.
  // Thus, the assembly of both vectors is handled manually here.
  Teuchos::RCP<Epetra_MultiVector> obj_der =
      params.get<Teuchos::RCP<Epetra_Vector>>("objective_derivations");
  Teuchos::RCP<Epetra_MultiVector> constr_der =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("constraints_derivations");

  if (optiparams_->DensType() == INPAR::TOPOPT::dens_node_based)
  {
    std::vector<int> lm;             // the same as lm
    std::vector<int> dummylmstride;  // not required
    std::vector<int> lmowner;        // owners of the lm-dofs
    ele->LocationVector(optidis, lm, lmowner, dummylmstride);

    LINALG::Assemble(*obj_der, 0, obj_deriv, lm, lmowner);
    LINALG::Assemble(*constr_der, 0, constr_deriv, lm, lmowner);
  }
  else if (optiparams_->DensType() == INPAR::TOPOPT::dens_ele_based)
  {
    std::vector<int> lm;
    lm.push_back(ele->Id());
    std::vector<int> lmowner;
    lmowner.push_back(ele->Owner());

    LINALG::Assemble(*obj_der, 0, obj_deriv, lm, lmowner);
    LINALG::Assemble(*constr_der, 0, constr_deriv, lm, lmowner);
  }
  else
    dserror("not implemented type of density function");

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Gradients(const int eid,
    std::map<int, LINALG::Matrix<nsd_, nen_>>& efluidvels,
    std::map<int, LINALG::Matrix<nsd_, nen_>>& eadjointvels,
    std::map<int, LINALG::Matrix<nen_, 1>>& eadjointpress, LINALG::Matrix<nen_, 1>& edens,
    Epetra_SerialDenseVector& obj_deriv, Epetra_SerialDenseVector& constr_deriv,
    Teuchos::RCP<MAT::Material> mat, DRT::UTILS::GaussIntegration& intpoints)
{
  LINALG::Matrix<nsd_, nen_> evel(true);

  const double timefac = optiparams_->Theta();
  const double timefac_old = 1.0 - optiparams_->Theta();

  double dissipation_fac = optiparams_->ObjDissipationFac();

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    double value = 0.0;

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad, eid);

    EvalPorosityAtIntPoint(edens);

    switch (optiparams_->TimeIntScheme())
    {
      case INPAR::FLUID::timeint_stationary:
      case INPAR::FLUID::timeint_one_step_theta:  // handle these two cases together
      {
        for (int timestep = 0; timestep <= optiparams_->NumTimesteps(); timestep++)
        {
          if (optiparams_->IsStationary()) timestep = 1;  // just one stationary time step 1


          //------------------------------------------------------------------------
          //  standard part due to objective function
          //------------------------------------------------------------------------
          double objfac = 0.0;
          if (timestep == 0)
            objfac = 1.0 - optiparams_->ThetaObj();
          else if (timestep == optiparams_->NumTimesteps())
            objfac = optiparams_->ThetaObj();
          else
            objfac = 1.0;

          if (optiparams_->OptiCase() == INPAR::TOPOPT::optitest_workflow_without_fluiddata)
            value -= objfac * dissipation_fac * 2 * edens.Dot(funct_);
          else
          {
            if (efluidvels.find(timestep) == efluidvels.end())
              dserror("Fluid solution of time step %i not found!", timestep);

            if ((optiparams_->IsStationary() == false) and (timestep > 0))
              fluidvelint_old_ = fluidvelint_;

            evel = efluidvels.find(timestep)->second;
            fluidvelint_.Multiply(evel, funct_);

            if (optiparams_->ObjDissipationTerm() == INPAR::TOPOPT::obj_diss_yes)
              value += objfac * dissipation_fac *
                       fluidvelint_.Dot(fluidvelint_);  // dissipation part, zero if not used


            //------------------------------------------------------------------------
            //  adjoint part (no influence of time step 0)
            //------------------------------------------------------------------------
            if (timestep > 0)
            {
              if (eadjointvels.find(timestep) == eadjointvels.end())
                dserror("Adjoint solution of time step %i not found!", timestep);

              evel = eadjointvels.find(timestep)->second;
              adjointvelint_.Multiply(evel, funct_);
              adjointvelxy_.MultiplyNT(evel, derxy_);

              CalcStabParameter(fac_);

              value += timefac * adjointvelint_.Dot(fluidvelint_);
              if (optiparams_->IsStationary() == false)
                value += timefac_old * adjointvelint_.Dot(fluidvelint_old_);

              if (optiparams_->PSPG())
              {
                eadjointpres_ = eadjointpress.find(timestep)->second;
                adjointpresxy_.Multiply(derxy_, eadjointpres_);

                value += timefac * (tau_(1) * fluidvelint_.Dot(adjointpresxy_));  // adjoint part
                if (optiparams_->IsStationary() == false)
                  value +=
                      timefac_old * (tau_(1) * fluidvelint_old_.Dot(
                                                   adjointpresxy_));  // instationary adjoint part
              }

              if (optiparams_->SUPG())
              {
                adjointconvvel_.Multiply(adjointvelxy_, fluidvelint_);
                value += timefac * (optiparams_->Density() * tau_(0) *
                                       fluidvelint_.Dot(adjointconvvel_));  // adjoint part
                if (optiparams_->IsStationary() == false)
                  value +=
                      timefac_old * (optiparams_->Density() * tau_(0) *
                                        fluidvelint_old_.Dot(adjointconvvel_));  // adjoint part
              }
            }
          }
        }
      }
      break;
      default:
      {
        dserror("unknown time integration scheme while evaluating objective gradient");
        break;
      }
    }

    if (optiparams_->OptiCase() == INPAR::TOPOPT::optitest_workflow_without_fluiddata)
      value *= fac_ * optiparams_->Dt();  // scale with integration factor and time step size
    else                                  // standard case
      value *= fac_ * optiparams_->Dt() *
               poroderdens_;  // scale with integration factor, time step size and poro-derivative

    if (optiparams_->DensType() == INPAR::TOPOPT::dens_node_based)
    {
      for (int vi = 0; vi < nen_; vi++)
      {
        obj_deriv(vi) += value * funct_(vi);
        constr_deriv(vi) += fac_ * funct_(vi);
      }
    }
    else if (optiparams_->DensType() == INPAR::TOPOPT::dens_ele_based)
    {
      obj_deriv(0) += value;
      constr_deriv(0) += fac_;
    }
    else
      dserror("not implemented type of density function");
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  return;
}



/*----------------------------------------------------------------------*
  | evaluate shape functions and derivatives at int. point     gjb 08/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::GaussIntegration::iterator& iquad,  ///< actual integration point
    const int eleid                                 ///< the element id
)
{
  // coordinates of the current integration point
  const double* gpcoord = iquad.Point();
  for (int idim = 0; idim < nsd_; idim++) xsi_(idim) = gpcoord[idim];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_, funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  xjm_.MultiplyNT(deriv_, xyze_);
  const double det = xji_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // compute integration factor
  fac_ = iquad.Weight() * det;

  // compute global derivatives
  derxy_.Multiply(xji_, deriv_);

}  // TopOptImpl::CalcSubgrVelocity



/*---------------------------------------------------------------------------------*
 | evaluate porosity at gauss point                               winklmaier 05/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::EvalPorosityAtIntPoint(LINALG::Matrix<nen_, 1> edens)
{
  // poro = poro_max + (poro_min - poro_max) * rho * (1+fac) / (rho+fac)
  double densint = edens.Dot(funct_);

  if (optiparams_->SmearFac() > -1.0e-15)  // >=0 as it should be -> standard case
  {
    poroint_ = optiparams_->MaxPoro() + (optiparams_->MinPoro() - optiparams_->MaxPoro()) *
                                            densint * (1 + optiparams_->SmearFac()) /
                                            (densint + optiparams_->SmearFac());

    // dporo/ddens = (poro_min - poro_max) * (fac+fac*fac) / (rho+fac)
    poroderdens_ = (optiparams_->MinPoro() - optiparams_->MaxPoro()) *
                   (optiparams_->SmearFac() + optiparams_->SmearFac() * optiparams_->SmearFac()) /
                   pow(densint + optiparams_->SmearFac(), 2);
  }
  else  // special cases
  {
    INPAR::TOPOPT::OptiCase testcase = (INPAR::TOPOPT::OptiCase)round(-optiparams_->SmearFac());

    switch (testcase)
    {
      case INPAR::TOPOPT::optitest_channel:
      case INPAR::TOPOPT::optitest_channel_with_step:
      case INPAR::TOPOPT::optitest_cornerflow:
      {
        dserror("Discontinuous setting here! Algorithm should stop earlier!");
        break;
      }
      case INPAR::TOPOPT::optitest_lin_poro:
      {
        double diff = optiparams_->MaxPoro() - optiparams_->MinPoro();
        double pmax = optiparams_->MaxPoro();

        poroint_ = -diff * densint + pmax;
        poroderdens_ = -diff;
        break;
      }
      case INPAR::TOPOPT::optitest_quad_poro:
      {
        double diff = optiparams_->MaxPoro() - optiparams_->MinPoro();
        double pmax = optiparams_->MaxPoro();

        double k = 0.1;

        poroint_ = (diff - k * pmax) * densint * densint + (-2 * diff + k * pmax) * densint + pmax;
        poroderdens_ = 2 * (diff - k * pmax) * densint + (-2 * diff + k * pmax);
        break;
      }
      case INPAR::TOPOPT::optitest_cub_poro:
      {
        double diff = optiparams_->MaxPoro() - optiparams_->MinPoro();
        double pmax = optiparams_->MaxPoro();

        double k1 = -50000.0;
        double k2 = 0.1;

        poroint_ = (2 * diff + k1 - k2 * pmax) * densint * densint * densint +
                   (-3 * diff - 2 * k1 + k2 * pmax) * densint * densint + k1 * densint + pmax;
        poroderdens_ = 3 * (2 * diff + k1 - k2 * pmax) * densint * densint +
                       2 * (-3 * diff - 2 * k1 + k2 * pmax) * densint + k1;
        break;
      }
      default:
      {
        dserror("you should not be here with a testcase not handled above");
        break;
      }
    }
  }
}



/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter             winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::CalcStabParameter(const double vol)
{
  //---------------------------------------------------------------------
  // preliminary definition of values which will already be computed for
  // tau_M and later be used for tau_C again by some of the subsequent
  // stabilization parameter definitions
  //---------------------------------------------------------------------
  double traceG = 0.0;
  double Gnormu = 0.0;
  double Gvisc = 0.0;

  double strle = 0.0;
  double hk = 0.0;
  double fluidvel_norm = 0.0;
  double re12 = 0.0;
  double c3 = 0.0;

  // material parameters
  const double dens = optiparams_->Density();
  const double visc = optiparams_->Viscosity();

  //---------------------------------------------------------------------
  // first step: computation of tau_M with the following options
  // (both with or without inclusion of dt-part):
  // A) definition according to Taylor et al. (1998)
  //    -> see also Gravemeier and Wall (2010) for version for
  //       variable-density flow at low Mach number
  // B) combined definition according to Franca and Valentin (2000) as
  //    well as Barrenechea and Valentin (2002)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // C) definition according to Shakib (1989) / Shakib and Hughes (1991)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // D) definition according to Codina (1998)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // E) definition according to Franca et al. (2005) as well as Badia
  //    and Codina (2010)
  //    -> only for Darcy or Darcy-Stokes/Brinkman flow, hence only
  //       tau_Mp for this definition
  //---------------------------------------------------------------------
  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // computation depending on which parameter definition is used
  switch (optiparams_->TauType())
  {
    case INPAR::FLUID::tau_taylor_hughes_zarins:
    case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
    {
      /*

      literature:
      1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
         of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
         (1998) 155-196.
      2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
         multigrid method for large-eddy simulation of turbulent variable-
         density flow at low Mach number, J. Comput. Phys. 229 (2010)
         6047-6070.
         -> version for variable-density low-Mach-number flow as implemented
            here, which corresponds to version for incompressible flow as
            given in the previous publications when density is constant

                                                                             1
                       +-                                               -+ - -
                       |        2                                        |   2
                       | c_1*rho                                  2      |
            tau  = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
               M       |     2                                           |
                       |   dt                                            |
                       +-                                               -+

            with the constants and covariant metric tensor defined as follows:

            C   = 1.0 (not explicitly defined here),
            c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
            c_2 = 1.0 (not explicitly defined here),
            c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)

                    +-           -+   +-           -+   +-           -+
                    |             |   |             |   |             |
                    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
              G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
               ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                    |    i     j  |   |    i     j  |   |    i     j  |
                    +-           -+   +-           -+   +-           -+

                    +----
                     \
            G : G =   +   G   * G
                     /     ij    ij
                    +----
                     i,j
                               +----
                               \
            rho*u*G*rho*u  =   +   rho*u * G  *rho*u
                               /        i   ij      j
                              +----
                                i,j
      */

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = poroint_;
      if (optiparams_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins or
          optiparams_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
          optiparams_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_scaled)
        sigma_tot += 1.0 / optiparams_->Dt();

      // definition of constants as described above
      const double c1 = 4.0;
      c3 = 12.0 / mk;

      // computation of various values derived from covariant metric tensor
      // (trace of covariant metric tensor required for computation of tau_C below)
      double G;
      double normG = 0.0;
      const double dens_sqr = dens * dens;
      for (int nn = 0; nn < nsd_; ++nn)
      {
        const double dens_sqr_velint_nn = dens_sqr * fluidvelint_(nn);
        for (int mm = 0; mm < nsd_; ++mm)
        {
          traceG += xji_(nn, mm) * xji_(nn, mm);
        }
        for (int rr = 0; rr < nsd_; ++rr)
        {
          G = xji_(nn, 0) * xji_(rr, 0);
          for (int mm = 1; mm < nsd_; ++mm)
          {
            G += xji_(nn, mm) * xji_(rr, mm);
          }
          normG += G * G;
          Gnormu += dens_sqr_velint_nn * G * fluidvelint_(rr);
        }
      }

      // compute viscous part
      Gvisc = c3 * visc * visc * normG;

      // computation of stabilization parameters tau_Mu and tau_Mp
      // -> identical for the present definitions
      tau_(0) = 1.0 / (sqrt(c1 * dens_sqr * DSQR(sigma_tot) + Gnormu + Gvisc));
      tau_(1) = tau_(0);
      break;
    }

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
    {
      /*

      literature:
      1) L.P. Franca, F. Valentin, On an improved unusual stabilized
         finite element method for the advective-reactive-diffusive
         equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.
      2) G.R. Barrenechea, F. Valentin, An unusual stabilized finite
         element method for a generalized Stokes problem, Numer. Math.
         92 (2002) 652-677.


                    xi1,xi2 ^
                            |      /
                            |     /
                            |    /
                          1 +---+
                            |
                            |
                            |
                            +--------------> re1,re2
                                1

      */
      // get velocity norm
      fluidvel_norm = fluidvelint_.Norm2();

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      const double sigma_tot =
          1.0 / optiparams_->Dt() +
          poroint_;  // originally timefac, not dt here -> sums up with timesteps

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel_norm, strle, hk);

      // various parameter computations for case with dt:
      // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
      const double re01 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(strle));
      const double re11 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(hk));

      // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
      const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
      re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

      // respective "switching" parameters
      const double xi01 = std::max(re01, 1.0);
      const double xi11 = std::max(re11, 1.0);
      const double xi02 = std::max(re02, 1.0);
      const double xi12 = std::max(re12, 1.0);

      tau_(0) = DSQR(strle) / (DSQR(strle) * dens * sigma_tot * xi01 + (4.0 * visc / mk) * xi02);
      tau_(1) = DSQR(hk) / (DSQR(hk) * dens * sigma_tot * xi11 + (4.0 * visc / mk) * xi12);
      break;
    }

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
    {
      /*

       stabilization parameter as above without inclusion of dt-part

      */
      // get velocity norm
      fluidvel_norm = fluidvelint_.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel_norm, strle, hk);

      // various parameter computations for case without dt:
      // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
      double re01 = 4.0 * visc / (mk * dens * poroint_ * DSQR(strle));
      double re11 = 4.0 * visc / (mk * dens * poroint_ * DSQR(hk));

      // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
      const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
      re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

      // respective "switching" parameters
      const double xi01 = std::max(re01, 1.0);
      const double xi11 = std::max(re11, 1.0);
      const double xi02 = std::max(re02, 1.0);
      const double xi12 = std::max(re12, 1.0);

      tau_(0) = DSQR(strle) / (DSQR(strle) * dens * poroint_ * xi01 + (4.0 * visc / mk) * xi02);
      tau_(1) = DSQR(hk) / (DSQR(hk) * dens * poroint_ * xi11 + (4.0 * visc / mk) * xi12);
      break;
    }

    case INPAR::FLUID::tau_shakib_hughes_codina:
    case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
    {
      /*

      literature:
      1) F. Shakib, Finite element analysis of the compressible Euler and
         Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
         Stanford University, Stanford, CA, USA, 1989.
      2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
         computational fluid dynamics: IX. Fourier analysis of space-time
         Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
         Engrg. 87 (1991) 35-58.
      3) R. Codina, Stabilized finite element approximation of transient
         incompressible flows using orthogonal subscales, Comput. Methods
         Appl. Mech. Engrg. 191 (2002) 4295-4321.

         constants defined as in Shakib (1989) / Shakib and Hughes (1991),
         merely slightly different with respect to c_3:

         c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
         c_2 = 4.0,
         c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

         Codina (2002) proposed present version without dt and explicit
         definition of constants
         (condition for constants as defined here: c_2 <= sqrt(c_3)).

      */
      // get velocity norm
      fluidvel_norm = fluidvelint_.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel_norm, strle, hk);

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = poroint_;
      if (optiparams_->TauType() == INPAR::FLUID::tau_shakib_hughes_codina)
        sigma_tot += 1.0 / optiparams_->Dt();

      // definition of constants as described above
      const double c1 = 4.0;
      const double c2 = 4.0;
      c3 = 4.0 / (mk * mk);
      // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

      tau_(0) = 1.0 / (sqrt(c1 * DSQR(dens) * DSQR(sigma_tot) +
                            c2 * DSQR(dens) * DSQR(fluidvel_norm) / DSQR(strle) +
                            c3 * DSQR(visc) / (DSQR(strle) * DSQR(strle))));
      tau_(1) = 1.0 / (sqrt(c1 * DSQR(dens) * DSQR(sigma_tot) +
                            c2 * DSQR(dens) * DSQR(fluidvel_norm) / DSQR(hk) +
                            c3 * DSQR(visc) / (DSQR(hk) * DSQR(hk))));
      break;
    }
    case INPAR::FLUID::tau_codina:
    case INPAR::FLUID::tau_codina_wo_dt:
    {
      /*

        literature:
           R. Codina, Comparison of some finite element methods for solving
           the diffusion-convection-reaction equation, Comput. Methods
           Appl. Mech. Engrg. 156 (1998) 185-210.

           constants:
           c_1 = 1.0 (for version with dt), 0.0 (for version without dt),
           c_2 = 2.0,
           c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

           Codina (1998) proposed present version without dt.

      */
      // get velocity norm
      fluidvel_norm = fluidvelint_.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel_norm, strle, hk);

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = poroint_;
      if (optiparams_->TauType() == INPAR::FLUID::tau_codina) sigma_tot += 1.0 / optiparams_->Dt();

      // definition of constants as described above
      const double c1 = 1.0;
      const double c2 = 2.0;
      c3 = 4.0 / mk;

      tau_(0) = 1.0 / (sqrt(c1 * dens * sigma_tot + c2 * dens * fluidvel_norm / strle +
                            c3 * visc / DSQR(strle)));
      tau_(1) =
          1.0 /
          (sqrt(c1 * dens * sigma_tot + c2 * dens * fluidvel_norm / hk + c3 * visc / DSQR(hk)));
      break;
    }
    default:
    {
      dserror("unknown definition for tau_M\n %i  ", optiparams_->TauType());
      break;
    }
  }  // end switch (fldAdPara_->whichtau_)


  //---------------------------------------------------------------------
  // second step: computation of tau_C with the following options:
  // A) definition according to Taylor et al. (1998)
  // B) definition according to Whiting (1999)/Whiting and Jansen (2001)
  // C) scaled version of definition according to Taylor et al. (1998)
  // D) definition according to Wall (1999)
  // E) definition according to Codina (2002)
  // F) definition according to Badia and Codina (2010)
  //    (only for Darcy or Darcy-Stokes/Brinkman flow)
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (optiparams_->TauType())
  {
    case INPAR::FLUID::tau_taylor_hughes_zarins:
    case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
    {
      /*

      literature:
         C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
         of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
         (1998) 155-196.

                                                1/2
                             (c_2*rho*u*G*rho*u)
                      tau  = -------------------
                         C       trace (G)


         -> see respective definitions for computation of tau_M above

      */

      tau_(2) = sqrt(Gnormu) / traceG;
    }
    break;

    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
    {
      /*

      literature:
      1) C.H. Whiting, Stabilized finite element methods for fluid dynamics
         using a hierarchical basis, PhD thesis, Rensselaer Polytechnic
         Institute, Troy, NY, USA, 1999.
      2) C.H. Whiting, K.E. Jansen, A stabilized finite element method for
         the incompressible Navier-Stokes equations using a hierarchical
         basis, Int. J. Numer. Meth. Fluids 35 (2001) 93-116.

                                    1.0
                      tau  = ------------------
                         C    tau  * trace (G)
                                 M

         -> see respective definitions for computation of tau_M above

      */

      tau_(2) = 1.0 / (tau_(0) * traceG);
    }
    break;

    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
    {
      /*

        Caution: This is an experimental version of a stabilization
                 parameter definition which scales the definition
                 for tau_C by Taylor et al. (1998) in a similar
                 way as proposed below by Franca and Frey (1992)
                 and Wall (1999) by appropriately defining an
                 element Reynolds number based on the covariant
                 metric tensor.

                    /                        1/2    \
                    |  /                    \       |                       1/2
                    | |  c_2*rho*u*G*rho*u  |       |    (c_2*rho*u*G*rho*u)
        tau  =  MIN | | ------------------- | | 1.0 | *  -------------------
           C        | |          2          |       |         trace (G)
                    | \    c_3*mu *G:G      /       |
                    \                               /
                      |                     |
                      -----------------------
                      element Reynolds number
                        based on covariant
                          metric tensor

         -> see respective definitions for computation of tau_M above

      */

      // element Reynolds number based on covariant metric tensor
      const double reG = std::sqrt(Gnormu / Gvisc);

      // "switching" parameter
      const double xi_tau_c = std::min(reG, 1.0);

      tau_(2) = xi_tau_c * sqrt(Gnormu) / traceG;
    }
    break;

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
    {
      /*

      literature:
      1) L.P. Franca, S.L. Frey, Stabilized finite element methods:
         II. The incompressible Navier-Stokes equations, Comput. Methods
         Appl. Mech. Engrg. 99 (1992) 209-293.
      2) W.A. Wall, Fluid-Struktur-Interaktion mit stabilisierten Finiten
         Elementen, Dissertation, Universitaet Stuttgart, 1999.

                   xi_tau_c ^
                            |
                          1 |   +-----------
                            |  /
                            | /
                            |/
                            +--------------> re12
                                1

         -> see respective definitions for computation of tau_M above

      */

      // "switching" parameter
      const double xi_tau_c = std::min(re12, 1.0);

      tau_(2) = 0.5 * dens * fluidvel_norm * hk * xi_tau_c;
    }
    break;

    case INPAR::FLUID::tau_shakib_hughes_codina:
    case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
    case INPAR::FLUID::tau_codina:
    case INPAR::FLUID::tau_codina_wo_dt:
    {
      /*

      literature:
         R. Codina, Stabilized finite element approximations of transient
         incompressible flows using orthogonal subscales, Comput. Methods
         Appl. Mech. Engrg. 191 (2002) 4295-4321.

         -> see respective definitions for computation of tau_M above

      */

      tau_(2) = DSQR(hk) / (sqrt(c3) * tau_(1));
    }
    break;
    default:
    {
      dserror("unknown definition for tau_C\n %i  ", optiparams_->TauType());
      break;
    }
  }  // end switch (fldAdPara_->whichtau_)

  return;
}



/*----------------------------------------------------------------------*
 |  calculation of characteristic element length       winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::CalcCharEleLength(
    const double vol, const double fluidvel_norm, double& strle, double& hk) const
{
  // cast dimension to a double varibale -> pow()
  const double dim = double(nsd_);

  //! direction of flow (normed velocity vector)
  LINALG::Matrix<nsd_, 1> fluidvelino;

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mu
  //---------------------------------------------------------------------
  // a) streamlength due to Tezduyar et al. (1992) -> default
  // normed velocity vector
  if (fluidvel_norm >= 1e-6)
    fluidvelino.Update(1.0 / fluidvel_norm, fluidvelint_);
  else
  {
    fluidvelino.Clear();
    fluidvelino(0, 0) = 1.0;
  }

  LINALG::Matrix<nen_, 1> tmp;
  tmp.MultiplyTN(derxy_, fluidvelino);
  const double val = tmp.Norm1();
  strle = 2.0 / val;

  // b) volume-equivalent diameter (warning: 3-D formula!)
  // strle = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // c) cubic/square root of element volume/area
  // strle = std::pow(vol,1/dim);

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mp
  //---------------------------------------------------------------------
  // a) volume-equivalent diameter -> default for 3-D computations
  if (nsd_ == 3) hk = std::pow((6. * vol / M_PI), (1.0 / 3.0)) / sqrt(3.0);

  // b) square root of element area -> default for 2-D computations,
  // may also alternatively be used for 3-D computations
  else if (nsd_ == 2)
    hk = std::pow(vol, 1 / dim);
  // check for potential 1-D computations
  else
    dserror("element length calculation not implemented for 1-D computation!");

  return;
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization& discretization,       ///< discretization
    const std::vector<int>& lm,                ///<
    LINALG::Matrix<nsd_, nen_>* matrixtofill,  ///< vector field
    LINALG::Matrix<nen_, 1>* vectortofill,     ///< scalar field
    Teuchos::RCP<Epetra_Vector>& globalvector  ///< global vector
    ) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector, mymatrix, lm);

  for (int inode = 0; inode < nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for (int idim = 0; idim < nsd_; ++idim)  // number of dimensions
      {
        (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * (nsd_ + 1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL) (*vectortofill)(inode, 0) = mymatrix[nsd_ + (inode * (nsd_ + 1))];
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundaryImplInterface* DRT::ELEMENTS::TopOptBoundaryImplInterface::Impl(
    const DRT::Element* ele)
{
  switch (ele->Shape())
  {
    case DRT::Element::quad4:
    {
      return TopOptBoundaryImpl<DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return TopOptBoundaryImpl<DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return TopOptBoundaryImpl<DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return TopOptBoundaryImpl<DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return TopOptBoundaryImpl<DRT::Element::tri6>::Instance();
    }
    case DRT::Element::line2:
    {
      return TopOptBoundaryImpl<DRT::Element::line2>::Instance();
    }
    case DRT::Element::line3:
    {
      return TopOptBoundaryImpl<DRT::Element::line3>::Instance();
    }
    default:
    {
      dserror(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
    }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype>* DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Instance(
    bool create)
{
  static TopOptBoundaryImpl<distype>* instance;
  if (create)  // create instance if not present
  {
    if (instance == NULL)
    {
      instance = new TopOptBoundaryImpl<distype>();
    }
  }
  else  // delete instance if present
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype>::TopOptBoundaryImpl()
    : intpoints_(distype), xsi_(true), fac_(0.0)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat,
    std::vector<int>& lm)
{
  // TODO coming...

  // work is done
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryGradients(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& optidis, Teuchos::RCP<MAT::Material> mat,
    std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra)
{
  // TODO coming...

  // work is done
  return 0;
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    DRT::Discretization& discretization,       ///< discretization
    const std::vector<int>& lm,                ///<
    LINALG::Matrix<nsd_, nen_>* matrixtofill,  ///< vector field
    LINALG::Matrix<nen_, 1>* vectortofill,     ///< scalar field
    Teuchos::RCP<Epetra_Vector>& globalvector  ///< global vector
    ) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector, mymatrix, lm);

  for (int inode = 0; inode < nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for (int idim = 0; idim < nsd_; ++idim)  // number of dimensions
      {
        (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * (nsd_ + 1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL) (*vectortofill)(inode, 0) = mymatrix[nsd_ + (inode * (nsd_ + 1))];
  }
}
