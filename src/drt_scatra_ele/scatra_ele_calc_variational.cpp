/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for ion-transport equation

\level 2

\maintainer Martin Kronbichler
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_variational.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_calc.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/scatra_mat_var_chemdiffusion.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
// To read if it is a semi-implicit scheme
#include "../drt_scatra/scatra_timint_var_chemdiffusion.H"
// To read initial conditions from input file
#include "../drt_lib/drt_globalproblem.H"



/*----------------------------------------------------------------------*
 | singleton destruction                                    deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Done()
{
  // delete instance
  Instance(0, 0, "", this);

  return;
}

/*----------------------------------------------------------------------*
 | singleton access method                                  deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalVariational<distype>*
DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalVariational* delete_me)
{
  static std::map<std::string, ScaTraEleCalVariational<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalVariational<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::string, ScaTraEleCalVariational<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}

/*----------------------------------------------------------------------*
 | protected constructor for singletons                    deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalVariational<distype>::ScaTraEleCalVariational(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      ephi0_(my::numdofpernode_, LINALG::Matrix<my::nen_, 1>(true))  // size of vector
      ,
      IsSemImplicitFunctional_(false)
{
  // replace standard internal variable manager by internal variable manager for variational
  // chemical diffusion
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerVar<my::nsd_, my::nen_>(my::numscal_));

  // replace standard scatra diffusion manager by variational chemical diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerVar<my::nsd_, my::nen_>(my::numscal_));

  return;
}

/*----------------------------------------------------------------------*
 | Action type: Evaluate                                    deanda 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractInitialNodeValues(ele, params, discretization, la);

  // Checks if the functional is evaluated ina fully implicit or semi-implicit fashion
  IsSemImplicitFunctional_ = params.get<bool>("Is_semImplicit_Functional");

  // call base class routine
  my::Evaluate(ele, params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
      elevec2_epetra, elevec3_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 | Action type: Evaluate                                    deanda 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::ExtractInitialNodeValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // extract local values from the global vectors initial conditions
  Teuchos::RCP<const Epetra_Vector> phi0 = discretization.GetState("phi0");
  if (phi0 == Teuchos::null) dserror("Cannot get the initial state vector 'phi0'");
  const std::vector<int>& lm = la[0].lm_;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phi0, ephi0_, lm);

  // extract local values from the global vectors at time t(n)
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phin == Teuchos::null) dserror("Cannot get state vector 'phin'");
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phin, my::ephin_, lm);
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                deanda 08/17 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Sysmat(
    DRT::Element* ele,                   ///< the element whose matrix is calculated
    Epetra_SerialDenseMatrix& emat,      ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs,      ///< element rhs to calculate
    Epetra_SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
)
{
  //-----------------------------------------------------------------------------------------------
  // calculate material and stabilization parameters (one per transported scalar) at element center
  //-----------------------------------------------------------------------------------------------
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(my::numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(my::numscal_, 1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(my::numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // stabilization variables
  std::vector<double> tau(my::numscal_, 0.);
  std::vector<LINALG::Matrix<my::nen_, 1>> tauderpot(
      my::numscal_, LINALG::Matrix<my::nen_, 1>(true));

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    // get material parameters (evaluation at integration point)
    if (my::scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // get material fields interpolated (needed for elements where the positive values of c can get
    // lost locally)
    const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(ele->Material());
    if ((distype == DRT::Element::line3 || distype == DRT::Element::tri6 ||
            distype == DRT::Element::quad4 || distype == DRT::Element::quad9) &&
        actmat->ModelType() == MAT::PAR::chemdiff_fickean)
      GetMaterialInterpolatedFields(ele, my::ephinp_);

    // calculate contributions to element matrix and right-hand side at integration point
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;
    const double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

    // loop all scalars
    for (int k = 0; k < my::numscal_; ++k)
    {
      const double taufac = tau[k] * fac;
      const double timetaufac = my::scatraparatimint_->TimeFac() * taufac;
      const double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;

      // TODO compute rhs containing bodyforce
      double rhsint(0.0);
      my::GetRhsInt(rhsint, densnp[k], k);

      // Compute element matrix and rhs
      CalcMatAndRhs(emat, erhs, k, intpoints.IP().qwgt[iquad], fac, timefacfac, rhsfac, taufac,
          timetaufac, rhstaufac, tauderpot[k], rhsint);
    }  // end loop over scalar
  }
  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs                         deanda  08/17|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalcMatAndRhs(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
    Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
    const int k,                     //!< index of current scalar
    const double weight,             //!< quadrature weight at integration point
    const double fac,                //!< domain-integration factor
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    const double taufac,      //!< tau times domain-integration factor
    const double timetaufac,  //!< domain-integration factor times tau times time-integration factor
    const double
        rhstaufac,  //!< time-integration factor for rhs times tau times domain-integration factor
    LINALG::Matrix<my::nen_, 1>&
        tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
    double& rhsint  //!< rhs at Gauss point
)
{
  // Sets the offset for the Chemical potential field in the Matrix indexes
  const int offsetChemPotential = my::numscal_;
  //----------------------------------------------------------------------
  // 1) RHS: Residuals
  //----------------------------------------------------------------------
  // Concentration R_{c}
  CalcRhsGeneric(erhs, k, fac,
      DiffManager()->GetInternalEnergy1deriv(k) - VarManager()->ConjugatePhinp(k), my::funct_);

  // Chemical potential R_{\mu} (local part)
  CalcRhsGeneric(erhs, k + offsetChemPotential, fac, VarManager()->Phin(k) - VarManager()->Phinp(k),
      my::funct_);

  // Chemical potential R_{\mu} (global part)
  CalcRhsGeneric(erhs, k + offsetChemPotential, my::scatraparatimint_->Dt() * fac,
      DiffManager()->GetDissipationPot1deriv(k), my::derxy_);

  // Chemical potential R_{\mu} (body force / external source)
  CalcRhsGeneric(
      erhs, k + offsetChemPotential, my::scatraparatimint_->Dt() * fac, rhsint, my::funct_);

  //----------------------------------------------------------------------
  // 1) MATRIX: Tangent concentration
  //----------------------------------------------------------------------
  // K_{c , c}
  CalcMatGeneric(
      emat, k, k, fac, DiffManager()->GetInternalEnergy2deriv(k), my::funct_, my::funct_);

  // K_{\mu , \mu}
  CalcMatGeneric(emat, k + offsetChemPotential, k + offsetChemPotential,
      -my::scatraparatimint_->Dt() * fac, DiffManager()->GetDissipationPot2deriv(k), my::derxy_,
      my::derxy_);

  // K_{c, \mu}
  CalcMatGeneric(emat, k, k + offsetChemPotential, fac, -1.0, my::funct_, my::funct_);

  // K_{\mu , c}  //TODO see if there is a better way to copy the off diagonal tahn to re-do it
  CalcMatGeneric(emat, k + offsetChemPotential, k, fac, -1.0, my::funct_, my::funct_);

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of generic element RHS                   deanda 08/17  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalcRhsGeneric(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double&
        discrFact,  //!< discrete factors(time &temporal discretization, weights, Jacobian, etc...)
    const double& kernel,                      //!< integral kernel (excluding shape functions)
    const LINALG::Matrix<my::nen_, 1>& sfunct  //!< first shape function
    ) const
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;  // fvi = vi * my::numdofpernode_+ k;
    erhs[fvi] -= discrFact * kernel * sfunct(vi);
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of generic element RHS                   deanda 08/17  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalcRhsGeneric(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double&
        discrFact,  //!< discrete factors(time &temporal discretization, weights, Jacobian, etc...)
    const LINALG::Matrix<my::nsd_, 1>& kernel,  //!< integral kernel (excluding shape functions)
    const LINALG::Matrix<my::nsd_, my::nen_>& sfunct  //!< first shape function
    ) const
{
  LINALG::Matrix<my::nen_, 1> temporalParam;
  temporalParam.MultiplyTN(sfunct, kernel);
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;  // fvi = vi * my::numdofpernode_+ k;
    erhs[fvi] -= discrFact * temporalParam(vi);
  }



  return;
}

/*------------------------------------------------------------------- *
 |  calculation of generic element matrix                deanda 08/17  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalcMatGeneric(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int& row,                  //!< starting row position index
    const int& col,                  //!< starting column position index
    const double&
        discrFact,  //!< discrete factors(time &temporal discretization, weights, Jacobian, etc...)
    const double& kernel,                       //!< integral kernel (excluding shape functions)
    const LINALG::Matrix<my::nen_, 1>& sfunct,  //!< first shape function
    const LINALG::Matrix<my::nen_, 1>& tfunct   //!< second shape function
    ) const
{
  const double densamfac = discrFact * kernel;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const double v = densamfac * tfunct(vi);
    const int fvi = vi * my::numdofpernode_ + row;  // fvi = vi * my::numdofpernode_ + row;
    // const int fvi = vi*numdof_ + row; (and row is 0 or 1)
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + col;  // fui = ui * my::numdofpernode_ + col;
      emat(fvi, fui) += v * sfunct(ui);
    }
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of generic element matrix                deanda 08/17  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalcMatGeneric(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int& row,                  //!< starting row position index
    const int& col,                  //!< starting column position index
    const double&
        discrFact,  //!< discrete factors(time &temporal discretization, weights, Jacobian, etc...)
    const LINALG::Matrix<my::nsd_, my::nsd_>&
        kernel,  //!< integral kernel (excluding shape functions)
    const LINALG::Matrix<my::nsd_, my::nen_>& sfunct,  //!< first shape function
    const LINALG::Matrix<my::nsd_, my::nen_>& tfunct   //!< second shape function
    ) const
{
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + row;  // vi * my::numdofpernode_ + row;

    LINALG::Matrix<my::nsd_, my::nen_>
        temporalParam1;  // TODO move the multiplications inside the for loop
    temporalParam1.MultiplyNN(kernel, sfunct);

    LINALG::Matrix<my::nen_, my::nen_> temporalParam2;
    temporalParam2.MultiplyTN(tfunct, temporalParam1);

    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + col;  // ui * my::numdofpernode_ + col;
      emat(fvi, fui) += discrFact * temporalParam2(vi, ui);
    }
  }

  // TODO Check if this is better (or something matrix operations alike)

  return;
}

/*----------------------------------------------------------------------*
 | get material parameters                                  deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::GetMaterialParams(
    const DRT::Element* ele,      //!< the element we are dealing with
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< id of current gauss point (default = -1)
)
{
  // get material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // evaluate substance material
  if (material->MaterialType() == INPAR::MAT::m_var_chemdiffusion)
  {
    for (int k = 0; k < my::numscal_; ++k)
    {
      Materials(material, k, densn[k], densnp[k], densam[k], visc, iquad);
    }
  }
  else
    dserror("Material type not supported!");

  return;
}  // DRT::ELEMENTS::ScaTraEleCalVariational<distype>::GetMaterialParams


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                   deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point

)
{
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_var_chemdiffusion:
    {
      my::MatScaTra(material, k, densn, densnp, densam, visc, iquad);
      MatScaTra_Var_ChemDiff(material, k, densn, densnp, densam, visc, iquad);
      break;
    }

    default:
    {
      dserror("Material type %i is not supported!", material->MaterialType());
      break;
    }
  }

  return;
}  // ScaTraEleCalVariational<distype>::Materials

/*----------------------------------------------------------------------*
 |  Material MatScaTra_Var_ChemDiff                         deanda 08/17 |
 *----------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::MatScaTra_Var_ChemDiff(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point (default = -1)
)
{
  const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(material);

  // get reference Chemical potential
  DiffManager()->SetRefMu(VarManager()->ConjugatePhi0(k), k);  // SetRefMu(actmat->RefMu(),k);

  // get reference concentration
  DiffManager()->SetRefC(VarManager()->Phi0(k), k);  // SetRefC(actmat->RefC(),k);

  // get RT factor
  DiffManager()->SetRT(actmat->RefTemp(), actmat->GasConstant());

  // get mobility
  DiffManager()->SetMobility(DiffManager()->GetIsotropicDiff(k), DiffManager()->GetRT(), k);

  if (IsSemImplicitFunctional_)  //(SCATRA::TimIntVarChemDiffusionOST::IsSemImplicitFunctional())
  {
    // get Dissipation Potential 1st derivative
    DiffManager()->SetDissipationPot1deriv(actmat, VarManager()->Phin(k), DiffManager()->GetRefC(k),
        VarManager()->ConjugateField(k), DiffManager()->GetMobility(k), k);

    // get Dissipation Potential 2nd derivative
    DiffManager()->SetDissipationPot2deriv(actmat, VarManager()->Phin(k), DiffManager()->GetRefC(k),
        VarManager()->ConjugateField(k), DiffManager()->GetMobility(k), k);
  }  // if: Semi-implicit evaluation
  else
  {
    // get Dissipation Potential 1st derivative
    DiffManager()->SetDissipationPot1deriv(actmat, VarManager()->Phinp(k),
        DiffManager()->GetRefC(k), VarManager()->ConjugateField(k), DiffManager()->GetMobility(k),
        k);

    // get Dissipation Potential 2nd derivative
    DiffManager()->SetDissipationPot2deriv(actmat, VarManager()->Phinp(k),
        DiffManager()->GetRefC(k), VarManager()->ConjugateField(k), DiffManager()->GetMobility(k),
        k);
  }  // if: Fully implicit evaluation

  if ((distype != DRT::Element::line3 && distype != DRT::Element::tri6 &&
          distype != DRT::Element::quad4 && distype != DRT::Element::quad9) ||
      actmat->ModelType() !=
          MAT::PAR::chemdiff_fickean)  // Needed when local negative concentrations may arise due to
                                       // interpolation
  {
    // get Internal energy 1st derivative
    if (VarManager()->Phinp(k) < 0)
      dserror("Negative concentration in the quadrature points, before computing Internal energy");
    DiffManager()->SetInternalEnergy1deriv(actmat, VarManager()->Phinp(k),
        DiffManager()->GetRefMu(k), DiffManager()->GetRefC(k), DiffManager()->GetRT(), k);

    // get Internal energy 2nd derivative
    DiffManager()->SetInternalEnergy2deriv(actmat, VarManager()->Phinp(k),
        DiffManager()->GetRefMu(k), DiffManager()->GetRefC(k), DiffManager()->GetRT(), k);
  }

  return;
}  // ScaTraEleCalVariational<distype>::MatScaTra_Var_ChemDiff

/*----------------------------------------------------------------------*
 | get material interpolated fields                         deanda 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::GetMaterialInterpolatedFields(
    const DRT::Element* ele,  //!< the element we are dealing with
    const std::vector<LINALG::Matrix<my::nen_, 1>>&
        ephinp  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
)
{
  // get actual material
  Teuchos::RCP<const MAT::Material> material = ele->Material();
  const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(material);

  //! scalar interpolated internal energy 1st derivative at t_(n+1)
  std::vector<double> InterpolatedIntEnergynp_1derv;
  //! scalar interpolated internal energy 2nd derivative at t_(n+1)
  std::vector<double> InterpolatedIntEnergynp_2derv;

  std::vector<LINALG::Matrix<my::nen_, 1>> NodalValue1(
      my::numscal_, LINALG::Matrix<my::nen_, 1>(true));
  std::vector<LINALG::Matrix<my::nen_, 1>> NodalValue2(
      my::numscal_, LINALG::Matrix<my::nen_, 1>(true));

  for (int k = 0; k < my::numscal_; ++k)
  {
    // calculate interpolated values for internal energy at t_(n+1)
    for (unsigned j = 0; j < my::nen_; ++j)
    {
      if (ephinp.at(k)(j) < 0.0)
        dserror(
            "Negative concentration at the nodes, before trying to interpolate internal energy");
      // Evaluate internal energy first derivative at the nodes
      NodalValue1.at(k)(j) = actmat->ComputeInternalEnergy(ephinp.at(k)(j),
          DiffManager()->GetRefMu(k), DiffManager()->GetRefC(k), DiffManager()->GetRT(), 1);
      // Evaluate internal energy second derivative at the nodes
      NodalValue2.at(k)(j) = actmat->ComputeInternalEnergy(ephinp.at(k)(j),
          DiffManager()->GetRefMu(k), DiffManager()->GetRefC(k), DiffManager()->GetRT(), 2);
    }
    // Replace internal energy first derivative by its interpolated version
    DiffManager()->ChangeInternalEnergy1deriv(my::funct_.Dot(NodalValue1.at(k)), k);
    // Replace internal energy second derivative by its interpolated version
    DiffManager()->ChangeInternalEnergy2deriv(my::funct_.Dot(NodalValue2.at(k)), k);
  }

  // TODO: For the dual energies the for loop looks like the one below
  //  for (int k = ConjugateOffsetIndex; k < 2*my::numscal_; ++k)

  return;
}  // DRT::ELEMENTS::ScaTraEleCalVariational<distype>::GetMaterialParams

/*------------------------------------------------------------------------------*
 | set internal variables for variaitonal chemical diffusion        deanda 08/17 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables derived from chemical diffusion under a variational formulation
  VarManager()->SetInternalVariablesVariationalChemDiff(
      my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_, ephi0_);
  return;
}  // DRT::ELEMENTS::ScaTraEleCalVariational<distype>::SetInternalVariablesForMatAndRHS()


/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution       deanda 08/17 |
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& errors)
{
  // set constants for analytical solution
  const double t = my::scatraparatimint_->Time();

  // safety checks
  if (DRT::INPUT::get<SCATRA::Action>(params, "action") != SCATRA::calc_error)
    dserror("This area is supposed to be to evaluate errors!");
  if (my::numscal_ != 1) dserror("Numscal_ != 1 The solution is only coded for 1 species so far.");

  if (t != 0.0)
  {
    switch (DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag"))
    {
      case INPAR::SCATRA::calcerror_AnalyticSeries:
      {
        //   References:
        //   John Crank
        //   "Mathematics of Diffusion"
        //    Plane sheet problem
        //   1975, Ed II, 22-24 Eq(2.54)

        if (my::nsd_ == 1)
        {
          // Set constants for analytical solution
          // Series parameters
          if (!params.isParameter("error function number"))
            dserror(
                "You forgot to specify the number of iterations to be taken by the series. Specify "
                "it using the parameter CALCERRORNO in the datfile");
          const int series_end = params.get<int>("error function number");

          // Gets Dirichlet condition
          double Cext = params.get<double>("Dirichlet_values");
          // Gets current simulation time
          const double t = my::scatraparatimint_->Time();
          // Gets length of the domain
          double L = params.get<double>("length_domain");
          // Access material for post-process solution of chemical potential
          Teuchos::RCP<const MAT::Material> material = ele->Material();
          const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(material);


          // integration points and weights
          // more GP than usual due to (possible) cos/exp fcts in analytical solutions
          const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
              SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

          // working arrays
          LINALG::Matrix<1, 1> ChemPot_int(true);
          LINALG::Matrix<1, 1> ChemPot_Analytic(true);
          LINALG::Matrix<1, 1> Con_int(true);
          LINALG::Matrix<1, 1> Con_Analytic(true);
          LINALG::Matrix<1, 1> delta_ChemPot(true);
          LINALG::Matrix<1, 1> delta_Con(true);
          LINALG::Matrix<my::nsd_, 1> xint(true);

          // start loop over integration points
          for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
          {
            const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

            // density at t_(n)
            std::vector<double> densn(my::numscal_, 1.0);
            // density at t_(n+1) or t_(n+alpha_F)
            std::vector<double> densnp(my::numscal_, 1.0);
            // density at t_(n+alpha_M)
            std::vector<double> densam(my::numscal_, 1.0);

            // get values of concentration at integration point
            Con_int(0) = my::funct_.Dot(my::ephinp_[0]);
            // get chemical potential solution at integration point
            ChemPot_int = my::funct_.Dot(my::ephinp_[1]);
            // get global coordinate of integration point
            xint.Multiply(my::xyze_, my::funct_);

            // compute various constants
            const double D = DiffManager()->GetIsotropicDiff(0);
            const double mu_0 = DiffManager()->GetRefMu(0);
            const double c_0 = DiffManager()->GetRefC(0);
            const double RT = DiffManager()->GetRT();
            double Analytic_con;

            // Series parameters
            double sum1 = 0;
            double sum2 = 0;

            for (int m = 0; m < series_end; ++m)
            {
              sum1 = sum1 + pow(-1, m) * erfc(((2 * m + 1) * L - xint(0)) / (2 * sqrt(D * t)));
              sum2 = sum2 + pow(-1, m) * erfc(((2 * m + 1) * L + xint(0)) / (2 * sqrt(D * t)));
            }  // Series loop

            Analytic_con = (Cext - c_0) * (sum1 + sum2) + c_0;
            Con_Analytic = Analytic_con;

            // compute analytical solution for chemical potential (post-process relation, hence,
            // independet of dimmension)
            ChemPot_Analytic = actmat->ComputeInternalEnergy(Analytic_con, mu_0, c_0, RT, 1);

            // compute differences between analytical solution and numerical solution
            // delta_ChemPot = ChemPot_int - ChemPot_Analytic;
            delta_ChemPot.Update(1.0, ChemPot_int, -1.0, ChemPot_Analytic);
            delta_Con.Update(1.0, Con_int, -1.0, Con_Analytic);

            // add square to L2 error
            errors[0] += delta_Con(0) * delta_Con(0) * fac;  // Difference in concentration
            errors[1] +=
                delta_ChemPot(0) * delta_ChemPot(0) * fac;  // Difference in chemical potential
            errors[2] += Con_Analytic(0) * Con_Analytic(0) * fac;  // L2 Analytic concentration
            errors[3] += ChemPot_Analytic(0) * ChemPot_Analytic(0) * fac;  // L2 chemical potential

          }  // end of loop over integration points
        }    // if Analytic series dim ==1
        else
          dserror("Analytic solution not coded for this spatial dimension");
      }
      break;

      case INPAR::SCATRA::calcerror_byfunction:
      {
        const int errorfunctno = params.get<int>("error function number");

        // analytical solution
        double phi_exact(0.0);
        double deltaphi(0.0);
        //! spatial gradient of current scalar value
        LINALG::Matrix<my::nsd_, 1> gradphi(true);
        LINALG::Matrix<my::nsd_, 1> gradphi_exact(true);
        LINALG::Matrix<my::nsd_, 1> deltagradphi(true);
        LINALG::Matrix<1, 1> ChemPot_int(true);
        LINALG::Matrix<1, 1> ChemPot_Analytic(true);
        // LINALG::Matrix<1,1>         Con_int(true);
        LINALG::Matrix<1, 1> Con_Analytic(true);
        LINALG::Matrix<1, 1> delta_ChemPot(true);
        // LINALG::Matrix<1,1>         delta_Con(true);

        // more GP than usual due to (possible) cos/exp fcts in analytical solutions
        const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
            SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);
        // start loop over integration points
        for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
        {
          const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

          // get coordinates at integration point
          // gp reference coordinates
          LINALG::Matrix<my::nsd_, 1> xyzint(true);
          xyzint.Multiply(my::xyze_, my::funct_);

          // function evaluation requires a 3D position vector!!
          double position[3] = {0.0, 0.0, 0.0};

          for (unsigned dim = 0; dim < my::nsd_; ++dim) position[dim] = xyzint(dim);

          // Set constants for solution of dual variable
          // Access material for post-process solution of chemical potential
          Teuchos::RCP<const MAT::Material> material = ele->Material();
          const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
              Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(material);
          // compute various constants
          // const double D = DiffManager()->GetIsotropicDiff(0);
          const double mu_0 = DiffManager()->GetRefMu(0);
          const double c_0 = DiffManager()->GetRefC(0);
          const double RT = DiffManager()->GetRT();

          for (int k = 0; k < my::numscal_; k += 2)
          {
            // scalar at integration point at time step n+1
            const double phinp = my::funct_.Dot(my::ephinp_[k]);
            // get chemical potential solution at integration point
            ChemPot_int = my::funct_.Dot(my::ephinp_[k + 1]);
            // spatial gradient of current scalar value
            gradphi.Multiply(my::derxy_, my::ephinp_[k]);

            phi_exact = DRT::Problem::Instance()->Funct(errorfunctno - 1).Evaluate(k, position, t);
            Con_Analytic = phi_exact;

            // compute analytical solution for chemical potential (post-process relation, hence,
            // independet of dimmension)
            ChemPot_Analytic = actmat->ComputeInternalEnergy(phi_exact, mu_0, c_0, RT, 1);

            // compute differences between analytical solution and numerical solution
            // delta_ChemPot = ChemPot_int - ChemPot_Analytic;
            delta_ChemPot.Update(1.0, ChemPot_int, -1.0, ChemPot_Analytic);
            // delta_Con.Update(1.0,Con_int,-1.0,Con_Analytic);
            deltaphi = phinp - phi_exact;

            // add square to L2 error
            errors[k + 0] += deltaphi * deltaphi * fac;  // Difference in concentration
            // errors[k+0] += delta_Con(0)*delta_Con(0)*fac;       // Difference in concentration
            errors[k + 1] +=
                delta_ChemPot(0) * delta_ChemPot(0) * fac;  // Difference in chemical potential
            errors[k + 2] += Con_Analytic(0) * Con_Analytic(0) * fac;  // L2 Analytic concentration
            errors[k + 3] +=
                ChemPot_Analytic(0) * ChemPot_Analytic(0) * fac;  // L2 chemical potential
          }
        }  // loop over integration points
        break;
      }
      default:
      {
        my::CalErrorComparedToAnalytSolution(ele, params, errors);
        break;
      }
    }  // switch(errortype)
  }
  return;
}  // CalErrorComparedToAnalytSolution


/*----------------------------------------------------------------------------------------*
 | finite difference check on element level (for debugging only) (protected)  deanda 08/17 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalVariational<distype>::FDCheck(DRT::Element* ele,
    Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs,
    Epetra_SerialDenseVector& subgrdiff)
{
  // screen output
  std::cout << "FINITE DIFFERENCE CHECK FOR ELEMENT " << ele->Id();

  // make a copy of state variables to undo perturbations later
  std::vector<LINALG::Matrix<my::nen_, 1>> ephinp_original(my::numdofpernode_);
  for (int k = 0; k < my::numdofpernode_; ++k)
    for (unsigned i = 0; i < my::nen_; ++i) ephinp_original[k](i, 0) = my::ephinp_[k](i, 0);

  // generalized-alpha time integration requires a copy of history variables as well
  std::vector<LINALG::Matrix<my::nen_, 1>> ehist_original(my::numscal_);
  if (my::scatraparatimint_->IsGenAlpha())
  {
    for (int k = 0; k < my::numscal_; ++k)
      for (unsigned i = 0; i < my::nen_; ++i) ehist_original[k](i, 0) = my::ehist_[k](i, 0);
  }

  // initialize element matrix and vectors for perturbed state
  Epetra_SerialDenseMatrix emat_dummy(emat);
  Epetra_SerialDenseVector erhs_perturbed(erhs);
  Epetra_SerialDenseVector subgrdiff_dummy(subgrdiff);

  // initialize counter for failed finite difference checks
  unsigned counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over columns of element matrix by first looping over nodes and then over dofs at each node
  for (unsigned inode = 0; inode < my::nen_; ++inode)
  {
    for (int idof = 0; idof < my::numdofpernode_; ++idof)
    {
      // number of current column of element matrix
      unsigned col = inode * my::numdofpernode_ + idof;
      // clear element matrix and vectors for perturbed state
      emat_dummy.Scale(0.0);
      erhs_perturbed.Scale(0.0);
      subgrdiff_dummy.Scale(0.0);

      // fill state vectors with original state variables
      for (int k = 0; k < my::numdofpernode_; ++k)
        for (unsigned i = 0; i < my::nen_; ++i) my::ephinp_[k](i, 0) = ephinp_original[k](i, 0);
      if (my::scatraparatimint_->IsGenAlpha())
        for (int k = 0; k < my::numscal_; ++k)
          for (unsigned i = 0; i < my::nen_; ++i) my::ehist_[k](i, 0) = ehist_original[k](i, 0);

      // impose perturbation
      if (my::scatraparatimint_->IsGenAlpha())
      {
        // perturbation of phi(n+alphaF), not of phi(n+1) => scale epsilon by factor alphaF
        my::ephinp_[idof](inode, 0) +=
            my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();

        // perturbation of phi(n+alphaF) by alphaF*epsilon corresponds to perturbation of phidtam
        // (stored in ehist_) by alphaM*epsilon/(gamma*dt); note: alphaF/timefac = alphaM/(gamma*dt)
        if (idof < my::numscal_)
          my::ehist_[idof](inode, 0) += my::scatraparatimint_->AlphaF() /
                                        my::scatraparatimint_->TimeFac() *
                                        my::scatrapara_->FDCheckEps();
      }
      else
        my::ephinp_[idof](inode, 0) += my::scatrapara_->FDCheckEps();

      // calculate element right-hand side vector for perturbed state
      Sysmat(ele, emat_dummy, erhs_perturbed, subgrdiff_dummy);

      // Now we compare the difference between the current entries in the element matrix
      // and their finite difference approximations according to
      // entries ?= (-erhs_perturbed + erhs_original) / epsilon

      // Note that the element right-hand side equals the negative element residual.
      // To account for errors due to numerical cancellation, we additionally consider
      // entries - erhs_original / epsilon ?= -erhs_perturbed / epsilon

      // Note that we still need to evaluate the first comparison as well. For small entries in the
      // element matrix, the second comparison might yield good agreement in spite of the entries
      // being wrong!
      for (unsigned row = 0; row < static_cast<unsigned>(my::numdofpernode_ * my::nen_); ++row)
      {
        // get current entry in original element matrix
        const double entry = emat(row, col);

        // finite difference suggestion (first divide by epsilon and then subtract for better
        // conditioning)
        const double fdval = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps() +
                             erhs(row) / my::scatrapara_->FDCheckEps();

        // confirm accuracy of first comparison
        if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in first comparison
        const double abserr1 = entry - fdval;
        if (abs(abserr1) > abs(maxabserr)) maxabserr = abserr1;
        double relerr1(0.);
        if (abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if (abs(fdval) > 1.e-17)
          relerr1 = abserr1 / abs(fdval);
        if (abs(relerr1) > abs(maxrelerr)) maxrelerr = relerr1;

        // evaluate first comparison
        if (abs(relerr1) > my::scatrapara_->FDCheckTol())
        {
          if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
          std::cout << "emat[" << row << "," << col << "]:  " << entry << "   ";
          std::cout << "finite difference suggestion:  " << fdval << "   ";
          std::cout << "absolute error:  " << abserr1 << "   ";
          std::cout << "relative error:  " << relerr1 << std::endl;
          std::cout << "numdofpernode_:  " << my::numdofpernode_ << std::endl;
          counter++;
        }

        // first comparison OK
        else
        {
          // left-hand side in second comparison
          const double left = entry - erhs(row) / my::scatrapara_->FDCheckEps();

          // right-hand side in second comparison
          const double right = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps();

          // confirm accuracy of second comparison
          if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
            dserror("Finite difference check involves values too close to numerical zero!");

          // absolute and relative errors in second comparison
          const double abserr2 = left - right;
          if (abs(abserr2) > abs(maxabserr)) maxabserr = abserr2;
          double relerr2(0.);
          if (abs(left) > 1.e-17)
            relerr2 = abserr2 / abs(left);
          else if (abs(right) > 1.e-17)
            relerr2 = abserr2 / abs(right);
          if (abs(relerr2) > abs(maxrelerr)) maxrelerr = relerr2;

          // evaluate second comparison
          if (abs(relerr2) > my::scatrapara_->FDCheckTol())
          {
            if (!counter) std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
            std::cout << "emat[" << row << "," << col << "]-erhs[" << row << "]/eps:  " << left
                      << "   ";
            std::cout << "-erhs_perturbed[" << row << "]/eps:  " << right << "   ";
            std::cout << "absolute error:  " << abserr2 << "   ";
            std::cout << "relative error:  " << relerr2 << std::endl;

            counter++;
          }
        }
      }
    }
  }

  // screen output in case finite difference check is passed
  if (!counter)
    std::cout << " --> PASSED WITH MAXIMUM ABSOLUTE ERROR " << maxabserr
              << " AND MAXIMUM RELATIVE ERROR " << maxrelerr << std::endl;

  // undo perturbations of state variables
  for (int k = 0; k < my::numdofpernode_; ++k)
    for (unsigned i = 0; i < my::nen_; ++i) my::ephinp_[k](i, 0) = ephinp_original[k](i, 0);
  if (my::scatraparatimint_->IsGenAlpha())
    for (int k = 0; k < my::numscal_; ++k)
      for (unsigned i = 0; i < my::nen_; ++i) my::ehist_[k](i, 0) = ehist_original[k](i, 0);

  return;
}

// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalVariational<DRT::Element::nurbs27>;
