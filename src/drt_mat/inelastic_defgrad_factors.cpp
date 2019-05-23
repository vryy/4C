/*----------------------------------------------------------------------*/
/*!
\brief evaluation of inelastic deformation gradients and their derivatives

\level 3

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/

/* headers */
#include "inelastic_defgrad_factors.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_structure.H"


/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors() : params_(NULL) {}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradFactors::InelasticDefgradFactors(MAT::PAR::Parameter* params) : params_(params)
{
}


/*--------------------------------------------------------------------*
 | factory                                              schmidt 03/18 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::InelasticDefgradFactors> MAT::InelasticDefgradFactors::Factory(int matnum)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");

  // another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case INPAR::MAT::mfi_lin_scalar_iso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarIso(curmat));
      MAT::PAR::InelasticDefgradLinScalarIso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradLinScalarIso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradLinScalarIso(params));
    }
    case INPAR::MAT::mfi_lin_scalar_aniso:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::InelasticDefgradLinScalarAniso(curmat));
      MAT::PAR::InelasticDefgradLinScalarAniso* params =
          dynamic_cast<MAT::PAR::InelasticDefgradLinScalarAniso*>(curmat->Parameter());
      return Teuchos::rcp(new InelasticDefgradLinScalarAniso(params));
    }
    default:
      dserror("cannot deal with type %d", curmat->Type());
  }
  // dummy return
  return Teuchos::null;
}

/*--------------------------------------------------------------------*
 | constructor                                          schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      Scalar1_(matdata->GetInt("SCALAR1")),
      Scalar1growthfac_(matdata->GetDouble("SCALAR1_GrowthFac")),
      Scalar1refconc_(matdata->GetDouble("SCALAR1_RefConc"))
{
  // safety checks
  // in case not all scatra dofs are transported scalars, the last scatra dof is a potential and can
  // not be treated as a concentration but it is treated like that in so3_scatra_evaluate.cpp in the
  // PreEvaluate method!
  if (Scalar1_ != 1) dserror("At the moment it is only possible that SCALAR1 induces growth");
  if (Scalar1growthfac_ < 0.0)
    dserror("The influence of scalar field SCALAR1 to growth can't be negativ");
  if (Scalar1refconc_ < 0.0) dserror("The reference concentration of SCALAR1 can't be negative");

  // check correct masslin type
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sdyn, "MASSLIN") != INPAR::STR::ml_none)
    dserror(
        "If you use the material 'InelasticDefgradLinScalarIso' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");

  return;
}


/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso()
    : InelasticDefgradFactors(), gp_(-1)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarIso::InelasticDefgradLinScalarIso(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params), gp_(-1)
{
}


/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
    LINALG::Matrix<3, 3>& iFinM)          ///< Inverse inelastic deformation gradient
{
  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();

  // calculate growth
  const double detjacobian = defgrad->Determinant();
  const double deltagrowth =
      Sc1growthfac * detjacobian * (concentrations_->at(gp_).at(Sc1 - 1) - Sc1refconc);
  const double isoinelasticdefo = pow(1.0 + deltagrowth, (1.0 / 3.0));

  // calculate inverse inelastic deformation gradient (FinM is modeled, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0 / isoinelasticdefo;

  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
    const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                         ///< inelastic contribution as 3x3 matrix
    const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
    const LINALG::Matrix<6, 9>&
        dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                   ///< deformation gradient of current inelastic contribution
    LINALG::Matrix<6, 6>& cmatadd)  ///< Additional elasticity tensor
{
  // static variables
  static LINALG::Matrix<9, 6> diFinjdC(true);
  static LINALG::Matrix<9, 1> id9x1(true);

  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();

  // determinant of jacobian
  const double detjacobian = defgrad->Determinant();

  // calculate scalefac
  const double temp =
      Sc1growthfac * detjacobian * (concentrations_->at(gp_).at(Sc1 - 1) - Sc1refconc);
  const double scalefac = -temp / (6.0 * pow(1.0 + temp, 4.0 / 3.0));

  // diFinjdC = - growthfac*det(F)*(c - c_{ref})/(6*[1 + growthfac*det(F)*(c - c_{ref})]^(4/3)) * I
  diFinjdC.MultiplyNT(scalefac, id9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
    const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                         ///< inelastic contribution as 3x3 matrix
    const LINALG::Matrix<6, 9>&
        dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                   ///< deformation gradient of current inelastic contribution
    LINALG::Matrix<6, 1>&
        dstressdc)  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. concentration
{
  static LINALG::Matrix<9, 1> id9x1(true);
  // prepare id9x1 (identity matrix written as a 9x1 vector)
  for (int i = 0; i < 3; ++i) id9x1(i) = 1.0;

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();

  // determinant of jacobian
  const double detjacobian = defgrad->Determinant();

  // calculate scalefac
  const double temp = Sc1growthfac * detjacobian;
  const double scalefac =
      -temp /
      (3.0 * pow(1.0 + temp * (concentrations_->at(gp_).at(Sc1 - 1) - Sc1refconc), 4.0 / 3.0));

  // dstressdc = dSdiFinj : diFinjdc
  // diFinjdc = - growthfac*det(F)/(3*[1 + growthfac*det(F)*(c-c_{ref})]^(4/3)) I
  dstressdc.MultiplyNN(scalefac, dSdiFinj, id9x1, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarIso::PreEvaluate(
    Teuchos::ParameterList& params)  ///< parameter list as handed in from the element
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  return;
}


/*--------------------------------------------------------------------*
 | constructor                                          schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::PAR::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      Scalar1_(matdata->GetInt("SCALAR1")),
      Scalar1growthfac_(matdata->GetDouble("SCALAR1_GrowthFac")),
      Scalar1refconc_(matdata->GetDouble("SCALAR1_RefConc")),
      growthdirmat_(true)
{
  // safety checks
  // in the case when not all scatra dofs are transported scalars, the last scatra dof is a
  // potential and can not be treated as a concentration but it is treated like that in
  // so3_scatra_evaluate.cpp in the PreEvaluate!
  if (Scalar1_ != 1) dserror("At the moment it is only possible that SCALAR1 induces growth");
  if (Scalar1growthfac_ < 0.0)
    dserror("The influence of scalar field SCALAR1 to growth can't be negativ");
  if (Scalar1refconc_ < 0.0) dserror("The reference concentration of SCALAR1 can't be negative");

  // check correct masslin type
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  if (DRT::INPUT::IntegralValue<INPAR::STR::MassLin>(sdyn, "MASSLIN") != INPAR::STR::ml_none)
    dserror(
        "If you use the material 'InelasticDefgradLinScalarAniso' please set 'MASSLIN' in the "
        "STRUCTURAL DYNAMIC Section to 'None', or feel free to implement other possibility!");

  // get the growth direction vector from the input file and check length
  std::vector<double> growthdirvec = *(matdata->Get<std::vector<double>>("GrowthDirection"));
  if (growthdirvec.size() != 3)
    dserror(
        "Since we have a 3D problem here, vector that defines the growth direction also needs to "
        "have the size 3!");

  // fill matrix that determines the growth direction
  const double growthdirvecnorm =
      std::sqrt(pow(growthdirvec[0], 2.0) + pow(growthdirvec[1], 2.0) + pow(growthdirvec[2], 2.0));
  const double invquadrgrowthdirvecnorm = 1.0 / (growthdirvecnorm * growthdirvecnorm);

  // loop over all rows and colomns to fill the matrix and scale it correctly on the fly
  for (unsigned i = 0; i < 3; ++i)
  {
    for (unsigned j = 0; j < 3; ++j)
    {
      growthdirmat_(i, j) = invquadrgrowthdirvecnorm * growthdirvec[i] * growthdirvec[j];
    }
  }

  return;
}


/*--------------------------------------------------------------------*
 | construct empty material                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso()
    : InelasticDefgradFactors(), gp_(-1)
{
}


/*--------------------------------------------------------------------*
 | construct material with specific material params     schmidt 03/18 |
 *--------------------------------------------------------------------*/
MAT::InelasticDefgradLinScalarAniso::InelasticDefgradLinScalarAniso(MAT::PAR::Parameter* params)
    : InelasticDefgradFactors(params), gp_(-1)
{
}


/*--------------------------------------------------------------------*
 | evaluate the inverse of the inelastic deformation                  |
 | gradient                                             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateInverseInelasticDefGrad(
    const LINALG::Matrix<3, 3>* defgrad,  ///< Deformation gradient
    LINALG::Matrix<3, 3>& iFinM)          ///< Inverse inelastic deformation gradient
{
  // init and clear variable
  static LINALG::Matrix<3, 3> FinM(true);
  FinM.Clear();

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();
  const double detjacobian = defgrad->Determinant();

  // calculate growth
  const double deltagrowth =
      Sc1growthfac * detjacobian * (concentrations_->at(gp_).at(Sc1 - 1) - Sc1refconc);

  // calculate inelastic deformation gradient and its inverse
  for (int i = 0; i < 3; ++i) FinM(i, i) = 1.0;

  // finalize inelastic deformation gradient matrix (FinM is calculated, such that the volume change
  // is a linear function of the scalar (mapped to reference frame) that causes it)
  FinM.Update(deltagrowth, Parameter()->Growthdirmat(), 1.0);
  // calculate inverse of inelastic deformation gradient matrix
  iFinM.Invert(FinM);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate additional contribution to cmat             schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateAdditionalCmat(
    const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
    const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                         ///< inelastic contribution as 3x3 matrix
    const LINALG::Matrix<6, 1>& iCV,     ///< Inverse right Cauchy-Green tensor
    const LINALG::Matrix<6, 9>&
        dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                   ///< deformation gradient of current inelastic contribution
    LINALG::Matrix<6, 6>& cmatadd)  ///< Additional elasticity tensor
{
  static LINALG::Matrix<3, 3> temp(true);
  static LINALG::Matrix<3, 3> iFinjGiFinj(true);
  static LINALG::Matrix<9, 1> diFinjddetJ9x1(true);
  static LINALG::Matrix<9, 6> diFinjdC(true);

  // get parameters
  const double Sc1 = Parameter()->Scalar1();
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  const double Sc1refconc = Parameter()->Scalar1refconc();

  // determinant of jacobian
  const double detjacobian = defgrad->Determinant();

  const double scalefac =
      -Sc1growthfac * detjacobian * (concentrations_->at(gp_).at(Sc1 - 1) - Sc1refconc) / 2.0;

  // calculate F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor of F_{in}
  temp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  iFinjGiFinj.MultiplyNN(1.0, temp, iFinjM, 0.0);
  Matrix3x3to9x1(iFinjGiFinj, diFinjddetJ9x1);

  // diFinjdC = - growthfac*det(F)*(c-c_{ref})/2.0 (F_{in,j}^{-1} . G . F_{in,j}^{-1}) \otimes
  // C^{-1}
  diFinjdC.MultiplyNT(scalefac, diFinjddetJ9x1, iCV, 0.0);

  // cmatadd = 2 dSdiFinj : diFinjdC
  cmatadd.MultiplyNN(2.0, dSdiFinj, diFinjdC, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate the off-diagonal contribution to the                      |
 | stiffness matrix (for monolithic calculation)        schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::EvaluateODStiffMat(
    const LINALG::Matrix<3, 3>* const defgrad,  ///< Deformation gradient
    const LINALG::Matrix<3, 3>& iFinjM,  ///< Inverse inelastic deformation gradient of current
                                         ///< inelastic contribution as 3x3 matrix
    const LINALG::Matrix<6, 9>&
        dSdiFinj,  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. the inverse inelastic
                   ///< deformation gradient of current inelastic contribution
    LINALG::Matrix<6, 1>&
        dstressdc)  ///< Derivative of 2nd Piola Kirchhoff stresses w.r.t. concentration
{
  // static variables
  static LINALG::Matrix<3, 3> tmp(true);
  static LINALG::Matrix<3, 3> diFinjdcM(true);
  static LINALG::Matrix<9, 1> diFinjdc9x1(true);

  // get parameters
  const double Sc1growthfac = Parameter()->Scalar1growthfac();
  // determinant of jacobian
  const double detjacobian = defgrad->Determinant();
  // scalefac
  const double scalefac = -Sc1growthfac * detjacobian;

  // diFinjdc = - growthfac*det(F) F_{in,j}^{-1} . G . F_{in,j}^{-1} with F_{in,j}, the j-th factor
  // of F_{in}
  tmp.MultiplyNN(1.0, iFinjM, Parameter()->Growthdirmat(), 0.0);
  diFinjdcM.MultiplyNN(scalefac, tmp, iFinjM, 0.0);
  Matrix3x3to9x1(diFinjdcM, diFinjdc9x1);

  // dstressdc = dSdiFinj : diFinjdc
  dstressdc.MultiplyNN(1.0, dSdiFinj, diFinjdc9x1, 1.0);

  return;
}


/*--------------------------------------------------------------------*
 | pre evaluate                                         schmidt 03/18 |
 *--------------------------------------------------------------------*/
void MAT::InelasticDefgradLinScalarAniso::PreEvaluate(
    Teuchos::ParameterList& params)  ///< parameter list as handed in from the element
{
  // get Gauss point number
  gp_ = params.get<int>("gp", -1);
  if (gp_ == -1) dserror("No Gauss point number provided in material.");

  // set pointer to vector of gp_conc, only if gp is 0, because this is the first gp
  if (gp_ == 0)
    concentrations_ = params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc",
        Teuchos::rcp(new std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0))));

  return;
}
