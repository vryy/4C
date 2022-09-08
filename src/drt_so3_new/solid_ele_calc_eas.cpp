/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       with EAS element technology
\level 1
*/
/*----------------------------------------------------------------------*/

#include "solid_utils.H"
#include "solid_ele_calc.H"
#include "solid_ele_calc_eas.H"
#include "solid_ele_interface.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "solid_ele.H"
#include "../drt_mat/so3_material.H"
#include "solid_ele_eas_utils.H"

template <DRT::Element::DiscretizationType distype, int neas>
DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::SolidEleCalcEas()
    : DRT::ELEMENTS::SolidEleInterface::SolidEleInterface(),
      DRT::ELEMENTS::SolidEleCalc<distype>::SolidEleCalc()
{
}

template <DRT::Element::DiscretizationType distype, int neas>
int DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::nln_force_stiff_mass(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* stiff_ep, Epetra_SerialDenseMatrix* mass_ep,
    Epetra_SerialDenseVector* force_ep, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra, DRT::UTILS::GaussIntegration* intpoints)
{
  std::cout << std::scientific;
  std::cout.precision(8);

  DRT::UTILS::GaussIntegration* my_intrule = my::default_integration_.getRawPtr();
  if (intpoints) my_intrule = intpoints;

  // get view
  Teuchos::RCP<LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_* my::nen_>> stiff = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_* my::nen_>> mass = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<my::nsd_ * my::nen_, 1>> force;
  if (stiff_ep)
    stiff =
        Teuchos::rcp(new LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>(*stiff_ep, true));
  if (mass_ep)
    mass =
        Teuchos::rcp(new LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>(*mass_ep, true));
  if (force_ep) force = Teuchos::rcp(new LINALG::Matrix<my::nsd_ * my::nen_, 1>(*force_ep, true));

  my::FillPositionArray(ele, discretization, lm);

  PrepareEAS(ele);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < my_intrule->NumPoints(); ++gp)
  {
    my::EvaluateShapeDeriv(*my_intrule, gp);

    my::JacobianMapping(*my_intrule, gp);

    EvaluateEASshape(ele->GetEAStype());

    my::Strains();
    my::StrainGradient();
    EnhanceEASstrains(ele);



    my::Material(&defgrd_enh_, &(my::gl_strain_), ele, params, gp);


    if (force != Teuchos::null) my::AddInternalForce(*force, my::fac_);
    //    if (!ele->Id()/*&&!gp*/) std::cout << "force "<<*force<< std::endl;

    if (stiff != Teuchos::null)
    {
      my::AddElasticStiff(*stiff, my::fac_);
      my::AddGeometricStiff(*stiff, my::fac_);
    }

    if (mass != Teuchos::null) my::AddMass(*mass, ele->SolidMaterial()->Density() * my::fac_);

    IntegrateEAS(ele);

  }  // gp loop


  if (stiff_ep != nullptr || force_ep != nullptr)
  {
    LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_> m(true);
    LINALG::Matrix<my::nen_ * my::nsd_, 1> f(true);
    //    if (!ele->Id())
    //      std::cout << "Kbb: " <<
    //      LINALG::Matrix<21,21>(ele->GetElementCondensator()->GetKbbInv(),true);
    ele->GetElementCondensator()->Invert();
    ele->GetElementCondensator()->Condense(
        m.A(), f.A(), STR::ELEMENTS::CondensatorBase::block_disp);
    if (stiff != Teuchos::null) stiff->Update(1., m, 1.);
    if (force != Teuchos::null) force->Update(1., f, 1.);
    //    if (!ele->Id())std::cout<<__FILE__<<__LINE__<<std::endl;
  }

  //  if(!ele->Id() && force!=Teuchos::null)
  //    std::cout << *force << std::endl;
  return 0;
}

template <DRT::Element::DiscretizationType distype, int neas>
int DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::nln_force_stiff_mass_gemm(
    DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
    Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
    Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra,
    DRT::UTILS::GaussIntegration* intpoints)
{
  dserror("not implemented");
  return 0;
}

template <DRT::Element::DiscretizationType distype, int neas>
int DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::UpdateElement(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
    Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra)
{
  return 0;
}

template <DRT::Element::DiscretizationType distype, int neas>
void DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::Setup(
    DRT::ELEMENTS::Solid* ele, DRT::INPUT::LineDefinition* linedef)
{
  if (linedef) DRT::ELEMENTS::SolidEleCalc<distype>::Setup(ele, linedef);

  ele->SetElementCondensator(Teuchos::rcp(new STR::ELEMENTS::Condensator<neas>));
  ele->GetElementCondensator()->AddBlock(
      STR::ELEMENTS::CondensatorBase::block_disp, my::nen_ * my::nsd_);
}

template <DRT::Element::DiscretizationType distype, int neas>
void DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::PrepareEAS(DRT::ELEMENTS::Solid* ele)
{
  ele->GetElementCondensator()->Clear();
  my::xi_.Clear();
  DRT::UTILS::shape_function<distype>(my::xi_, my::shapefunction_);
  DRT::UTILS::shape_function_deriv1<distype>(my::xi_, my::deriv_);

  // for now, this is not the inverse
  my::invJ_.Multiply(my::deriv_, my::xrefe_);

  // first, build T0^T transformation matrix which maps the M-matrix
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  // todo: make this work for 2D
  T0invT_(0, 0) = my::invJ_(0, 0) * my::invJ_(0, 0);
  T0invT_(1, 0) = my::invJ_(1, 0) * my::invJ_(1, 0);
  T0invT_(2, 0) = my::invJ_(2, 0) * my::invJ_(2, 0);
  T0invT_(3, 0) = 2 * my::invJ_(0, 0) * my::invJ_(1, 0);
  T0invT_(4, 0) = 2 * my::invJ_(1, 0) * my::invJ_(2, 0);
  T0invT_(5, 0) = 2 * my::invJ_(0, 0) * my::invJ_(2, 0);

  T0invT_(0, 1) = my::invJ_(0, 1) * my::invJ_(0, 1);
  T0invT_(1, 1) = my::invJ_(1, 1) * my::invJ_(1, 1);
  T0invT_(2, 1) = my::invJ_(2, 1) * my::invJ_(2, 1);
  T0invT_(3, 1) = 2 * my::invJ_(0, 1) * my::invJ_(1, 1);
  T0invT_(4, 1) = 2 * my::invJ_(1, 1) * my::invJ_(2, 1);
  T0invT_(5, 1) = 2 * my::invJ_(0, 1) * my::invJ_(2, 1);

  T0invT_(0, 2) = my::invJ_(0, 2) * my::invJ_(0, 2);
  T0invT_(1, 2) = my::invJ_(1, 2) * my::invJ_(1, 2);
  T0invT_(2, 2) = my::invJ_(2, 2) * my::invJ_(2, 2);
  T0invT_(3, 2) = 2 * my::invJ_(0, 2) * my::invJ_(1, 2);
  T0invT_(4, 2) = 2 * my::invJ_(1, 2) * my::invJ_(2, 2);
  T0invT_(5, 2) = 2 * my::invJ_(0, 2) * my::invJ_(2, 2);

  T0invT_(0, 3) = my::invJ_(0, 0) * my::invJ_(0, 1);
  T0invT_(1, 3) = my::invJ_(1, 0) * my::invJ_(1, 1);
  T0invT_(2, 3) = my::invJ_(2, 0) * my::invJ_(2, 1);
  T0invT_(3, 3) = my::invJ_(0, 0) * my::invJ_(1, 1) + my::invJ_(1, 0) * my::invJ_(0, 1);
  T0invT_(4, 3) = my::invJ_(1, 0) * my::invJ_(2, 1) + my::invJ_(2, 0) * my::invJ_(1, 1);
  T0invT_(5, 3) = my::invJ_(0, 0) * my::invJ_(2, 1) + my::invJ_(2, 0) * my::invJ_(0, 1);


  T0invT_(0, 4) = my::invJ_(0, 1) * my::invJ_(0, 2);
  T0invT_(1, 4) = my::invJ_(1, 1) * my::invJ_(1, 2);
  T0invT_(2, 4) = my::invJ_(2, 1) * my::invJ_(2, 2);
  T0invT_(3, 4) = my::invJ_(0, 1) * my::invJ_(1, 2) + my::invJ_(1, 1) * my::invJ_(0, 2);
  T0invT_(4, 4) = my::invJ_(1, 1) * my::invJ_(2, 2) + my::invJ_(2, 1) * my::invJ_(1, 2);
  T0invT_(5, 4) = my::invJ_(0, 1) * my::invJ_(2, 2) + my::invJ_(2, 1) * my::invJ_(0, 2);

  T0invT_(0, 5) = my::invJ_(0, 0) * my::invJ_(0, 2);
  T0invT_(1, 5) = my::invJ_(1, 0) * my::invJ_(1, 2);
  T0invT_(2, 5) = my::invJ_(2, 0) * my::invJ_(2, 2);
  T0invT_(3, 5) = my::invJ_(0, 0) * my::invJ_(1, 2) + my::invJ_(1, 0) * my::invJ_(0, 2);
  T0invT_(4, 5) = my::invJ_(1, 0) * my::invJ_(2, 2) + my::invJ_(2, 0) * my::invJ_(1, 2);
  T0invT_(5, 5) = my::invJ_(0, 0) * my::invJ_(2, 2) + my::invJ_(2, 0) * my::invJ_(0, 2);

  // now evaluate T0^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<my::numstr_, my::numstr_, 1> solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(T0invT_);
  int err2 = solve_for_inverseT0.Factor();
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) || (err2 != 0)) dserror("Inversion of T0inv (Jacobian0) failed");

  T0invT_.Scale(my::invJ_.Determinant());
}


template <DRT::Element::DiscretizationType distype, int neas>
void DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::EnhanceEASstrains(DRT::ELEMENTS::Solid* ele)
{
  my::gl_strain_.Multiply(
      1., shape_eas_, LINALG::Matrix<neas, 1>(ele->GetElementCondensator()->GetB(), true), 1.);
  if (ele->SolidMaterial()->NeedsDefgrd())
    STR::UTILS::CalcConsistentDefgrd<my::nsd_>(my::defgrd_, my::gl_strain_, defgrd_enh_);
  else
    defgrd_enh_ = my::defgrd_;
}


template <DRT::Element::DiscretizationType distype, int neas>
void DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::IntegrateEAS(DRT::ELEMENTS::Solid* ele)
{
  static LINALG::Matrix<my::numstr_, neas> cM;
  cM.Clear();
  cM.Multiply(1., my::cmat_, shape_eas_, 0.);
  LINALG::Matrix<neas, neas> kbb(ele->GetElementCondensator()->GetKbbInv(), true);
  kbb.MultiplyTN(my::fac_, shape_eas_, cM, 1.);

  LINALG::Matrix<neas, 1> rb(ele->GetElementCondensator()->GetRb(), true);
  rb.MultiplyTN(my::fac_, shape_eas_, my::pk2_, 1.);

  static LINALG::Matrix<neas, my::nen_ * my::nsd_> tmp;
  tmp.Clear();
  tmp.MultiplyTN(my::fac_, shape_eas_, my::cb_, 0.);
  LINALG::Matrix<neas, my::nen_ * my::nsd_>(
      ele->GetElementCondensator()->GetKba(STR::ELEMENTS::CondensatorBase::block_disp), true)
      .Update(1., tmp, 1.);
  LINALG::Matrix<my::nen_ * my::nsd_, neas>(
      ele->GetElementCondensator()->GetKab(STR::ELEMENTS::CondensatorBase::block_disp), true)
      .UpdateT(1., tmp, 1.);
}


template <DRT::Element::DiscretizationType distype, int neas>
void DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::EvaluateEASshape(
    const ::STR::ELEMENTS::EASType eastype)
{
  static LINALG::Matrix<my::numstr_, neas> M_gp(false);
  M_gp.Clear();
  switch (eastype)
  {
    case ::STR::ELEMENTS::EASType::eastype_h8_9:
    {
      M_gp(0, 0) = my::xi_(0);
      M_gp(1, 1) = my::xi_(1);
      M_gp(2, 2) = my::xi_(2);

      M_gp(3, 3) = my::xi_(0);
      M_gp(3, 4) = my::xi_(1);
      M_gp(4, 5) = my::xi_(1);
      M_gp(4, 6) = my::xi_(2);
      M_gp(5, 7) = my::xi_(0);
      M_gp(5, 8) = my::xi_(2);

      break;
    }
    case ::STR::ELEMENTS::EASType::eastype_h8_21:
    {
      M_gp(0, 0) = my::xi_(0);
      M_gp(0, 15) = my::xi_(0) * my::xi_(1);
      M_gp(0, 16) = my::xi_(0) * my::xi_(2);
      M_gp(1, 1) = my::xi_(1);
      M_gp(1, 17) = my::xi_(0) * my::xi_(1);
      M_gp(1, 18) = my::xi_(1) * my::xi_(2);
      M_gp(2, 2) = my::xi_(2);
      M_gp(2, 19) = my::xi_(0) * my::xi_(2);
      M_gp(2, 20) = my::xi_(1) * my::xi_(2);

      M_gp(3, 3) = my::xi_(0);
      M_gp(3, 4) = my::xi_(1);
      M_gp(3, 9) = my::xi_(0) * my::xi_(2);
      M_gp(3, 10) = my::xi_(1) * my::xi_(2);
      M_gp(4, 5) = my::xi_(1);
      M_gp(4, 6) = my::xi_(2);
      M_gp(4, 11) = my::xi_(0) * my::xi_(1);
      M_gp(4, 12) = my::xi_(0) * my::xi_(2);
      M_gp(5, 7) = my::xi_(0);
      M_gp(5, 8) = my::xi_(2);
      M_gp(5, 13) = my::xi_(0) * my::xi_(1);
      M_gp(5, 14) = my::xi_(1) * my::xi_(2);

      break;
    }
    default:
      dserror("unknown EAS type");
      break;
  }
  shape_eas_.Multiply(1. / my::detJ_, T0invT_, M_gp, 0.);
}



template <DRT::Element::DiscretizationType distype, int neas>
int DRT::ELEMENTS::SolidEleCalcEas<distype, neas>::RecoverCondensed(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
    Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra)
{
  // get access to the interface parameters
  const double step_length = ele->ParamsInterface().GetStepLength();

  // first, store the eas state of the previous accepted Newton step
  ele->ParamsInterface().SumIntoMyPreviousSolNorm(
      NOX::NLN::StatusTest::quantity_eas, neas, ele->GetElementCondensator()->GetB(), ele->Owner());

  if (ele->ParamsInterface().IsDefaultStep())
  {
    Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
    if (res == Teuchos::null) dserror("where's my residual displacement");
    std::vector<double> myres(lm.size());
    DRT::UTILS::ExtractMyValues(*res, myres, lm);
    std::map<STR::ELEMENTS::CondensatorBase::BlockName, double*> incr;
    incr[STR::ELEMENTS::CondensatorBase::block_disp] = &(myres[0]);
    ele->GetElementCondensator()->Recover(incr, step_length);
  }
  else
  {
    ele->GetElementCondensator()->ReduceStepLength(step_length);
  }

  ele->ParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, neas,
      ele->GetElementCondensator()->GetBincr(), ele->GetElementCondensator()->GetB(), step_length,
      ele->Owner());

  //  if (!ele->Id())std::cout <<
  //  LINALG::Matrix<21,1>(ele->GetElementCondensator()->GetBincr(),true) << std::endl;
  return 0;
}

// template classes
template class DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 9>;
template class DRT::ELEMENTS::SolidEleCalcEas<DRT::Element::hex8, 21>;
