/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of spatial functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_functions.H"
#include "lib_globalproblem.H"
#include "mat_par_bundle.H"
#include "mat_newtonianfluid.H"
#include "mat_fluid_murnaghantait.H"
#include "mat_fluid_linear_density_viscosity.H"
#include "mat_fluid_weakly_compressible.H"
#include "mat_stvenantkirchhoff.H"
#include <Teuchos_RCP.hpp>

namespace
{
  /// returns Weakly Compressible Fluid quick access parameters from given material id
  const MAT::PAR::WeaklyCompressibleFluid& GetWeaklyCompressibleFluidMatPars(int mat_id)
  {
    Teuchos::RCP<MAT::PAR::Material> mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
    if (mat->Type() != INPAR::MAT::m_fluid_weakly_compressible)
      dserror("Material %d is not a weakly compressible fluid", mat_id);
    MAT::PAR::Parameter* params = mat->Parameter();
    auto* fparams = dynamic_cast<MAT::PAR::WeaklyCompressibleFluid*>(params);
    if (!fparams) dserror("Material does not cast to Weakly compressible fluid");
    return *fparams;
  }

  /// returns Newton Fluid quick access parameters from given material id
  const MAT::PAR::NewtonianFluid& GetNewtonianFluidMatPars(int mat_id)
  {
    Teuchos::RCP<MAT::PAR::Material> mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
    if (mat->Type() != INPAR::MAT::m_fluid) dserror("Material %d is not a fluid", mat_id);
    MAT::PAR::Parameter* params = mat->Parameter();
    auto* fparams = dynamic_cast<MAT::PAR::NewtonianFluid*>(params);
    if (!fparams) dserror("Material does not cast to Newtonian fluid");
    return *fparams;
  }

  /// returns St. Venant Kirchhof quick access parameters from given material id
  const MAT::PAR::StVenantKirchhoff& GetSVKMatPars(int mat_id)
  {
    Teuchos::RCP<MAT::PAR::Material> mat = DRT::Problem::Instance()->Materials()->ById(mat_id);
    if (mat->Type() != INPAR::MAT::m_stvenant)
      dserror("Material %d is not a St.Venant-Kirchhoff structure material", mat_id);
    MAT::PAR::Parameter* params = mat->Parameter();
    auto* fparams = dynamic_cast<MAT::PAR::StVenantKirchhoff*>(params);
    if (!fparams) dserror("Material does not cast to St.Venant-Kirchhoff structure material");
    return *fparams;
  }
}  // namespace

void FLD::AddValidFluidFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines)
{
  DRT::INPUT::LineDefinition beltrami;
  beltrami.AddTag("BELTRAMI").AddNamedDouble("c1");

  DRT::INPUT::LineDefinition channelweaklycompressible;
  channelweaklycompressible.AddTag("CHANNELWEAKLYCOMPRESSIBLE");

  DRT::INPUT::LineDefinition correctiontermchannelweaklycompressible;
  correctiontermchannelweaklycompressible.AddTag("CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLE");

  DRT::INPUT::LineDefinition weaklycompressiblepoiseuille;
  weaklycompressiblepoiseuille.AddTag("WEAKLYCOMPRESSIBLE_POISEUILLE")
      .AddNamedInt("MAT")
      .AddNamedDouble("L")
      .AddNamedDouble("R")
      .AddNamedDouble("U");

  DRT::INPUT::LineDefinition weaklycompressiblepoiseuilleforce;
  weaklycompressiblepoiseuilleforce.AddTag("WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE")
      .AddNamedInt("MAT")
      .AddNamedDouble("L")
      .AddNamedDouble("R")
      .AddNamedDouble("U");

  DRT::INPUT::LineDefinition weaklycompressiblemanufacturedflow;
  weaklycompressiblemanufacturedflow.AddTag("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressiblemanufacturedflowforce;
  weaklycompressiblemanufacturedflowforce.AddTag("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfd;
  weaklycompressibleetiennecfd.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD").AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfdforce;
  weaklycompressibleetiennecfdforce.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfdviscosity;
  weaklycompressibleetiennecfdviscosity.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluid;
  weaklycompressibleetiennefsifluid.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluidforce;
  weaklycompressibleetiennefsifluidforce.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluidviscosity;
  weaklycompressibleetiennefsifluidviscosity
      .AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition beltramiup;
  beltramiup.AddTag("BELTRAMI-UP").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition beltramigradu;
  beltramigradu.AddTag("BELTRAMI-GRADU").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition beltramirhs;
  beltramirhs.AddTag("BELTRAMI-RHS")
      .AddNamedInt("MAT")
      .AddNamedInt("ISSTAT")
      .AddNamedInt("ISSTOKES");

  DRT::INPUT::LineDefinition kimmoinup;
  kimmoinup.AddTag("KIMMOIN-UP").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition kimmoingradu;
  kimmoingradu.AddTag("KIMMOIN-GRADU").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition kimmoinrhs;
  kimmoinrhs.AddTag("KIMMOIN-RHS").AddNamedInt("MAT").AddNamedInt("ISSTAT").AddNamedInt("ISSTOKES");

  DRT::INPUT::LineDefinition kimmoinstress;
  kimmoinstress.AddTag("KIMMOIN-STRESS")
      .AddNamedInt("MAT")
      .AddNamedInt("ISSTAT")
      .AddNamedDouble("AMPLITUDE");

  lines->Add(beltrami);
  lines->Add(channelweaklycompressible);
  lines->Add(correctiontermchannelweaklycompressible);
  lines->Add(weaklycompressiblepoiseuille);
  lines->Add(weaklycompressiblepoiseuilleforce);
  lines->Add(weaklycompressiblemanufacturedflow);
  lines->Add(weaklycompressiblemanufacturedflowforce);
  lines->Add(weaklycompressibleetiennecfd);
  lines->Add(weaklycompressibleetiennecfdforce);
  lines->Add(weaklycompressibleetiennecfdviscosity);
  lines->Add(weaklycompressibleetiennefsifluid);
  lines->Add(weaklycompressibleetiennefsifluidforce);
  lines->Add(weaklycompressibleetiennefsifluidviscosity);
  lines->Add(beltramiup);
  lines->Add(beltramigradu);
  lines->Add(beltramirhs);
  lines->Add(kimmoinup);
  lines->Add(kimmoingradu);
  lines->Add(kimmoinrhs);
  lines->Add(kimmoinstress);
}

Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> FLD::TryCreateFluidFunction(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs)
{
  if (function_line_defs.size() != 1) return Teuchos::null;

  const auto& function_lin_def = function_line_defs.front();

  if (function_lin_def->HaveNamed("BELTRAMI"))
  {
    double c1;
    function_lin_def->ExtractDouble("c1", c1);

    return Teuchos::rcp(new FLD::BeltramiFunction(c1));
  }
  else if (function_lin_def->HaveNamed("CHANNELWEAKLYCOMPRESSIBLE"))
  {
    return Teuchos::rcp(new FLD::ChannelWeaklyCompressibleFunction());
  }
  else if (function_lin_def->HaveNamed("CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLE"))
  {
    return Teuchos::rcp(new FLD::CorrectionTermChannelWeaklyCompressibleFunction());
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_POISEUILLE"))
  {
    // read data
    int mat_id = -1;
    double L = 0.0;
    double R = 0.0;
    double U = 0.0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractDouble("L", L);
    function_lin_def->ExtractDouble("R", R);
    function_lin_def->ExtractDouble("U", U);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_POISEUILLE");
    if (L <= 0) dserror("Please give a (reasonable) 'L' in WEAKLYCOMPRESSIBLE_POISEUILLE");
    if (R <= 0) dserror("Please give a (reasonable) 'R' in WEAKLYCOMPRESSIBLE_POISEUILLE");
    if (U <= 0) dserror("Please give a (reasonable) 'U' in WEAKLYCOMPRESSIBLE_POISEUILLE");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressiblePoiseuilleFunction(fparams, L, R, U));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE"))
  {
    // read data
    int mat_id = -1;
    double L = 0.0;
    double R = 0.0;
    double U = 0.0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractDouble("L", L);
    function_lin_def->ExtractDouble("R", R);
    function_lin_def->ExtractDouble("U", U);

    if (mat_id <= 0)
      dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
    if (L <= 0) dserror("Please give a (reasonable) 'L' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
    if (R <= 0) dserror("Please give a (reasonable) 'R' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
    if (U <= 0) dserror("Please give a (reasonable) 'U' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressiblePoiseuilleForceFunction(fparams, L, R, U));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW"))
  {
    // read data
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0)
      dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressibleManufacturedFlowFunction(fparams));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE"))
  {
    // read data
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0)
      dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressibleManufacturedFlowForceFunction(fparams));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD"))
  {
    // read data
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDFunction(fparams));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE"))
  {
    // read data
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0)
      dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDForceFunction(fparams));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY"))
  {
    // read data
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0)
      dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY");

    // get materials
    auto fparams = GetWeaklyCompressibleFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDViscosityFunction(fparams));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID"))
  {
    // read data
    int mat_id_fluid = -1;
    int mat_id_struc = -1;

    function_lin_def->ExtractInt("MAT_FLUID", mat_id_fluid);
    function_lin_def->ExtractInt("MAT_STRUC", mat_id_struc);

    if (mat_id_fluid <= 0)
      dserror("Please give a (reasonable) 'MAT_FLUID' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID");
    if (mat_id_struc <= 0)
      dserror("Please give a (reasonable) 'MAT_STRUC' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID");

    // get materials
    auto fparams_fluid = GetWeaklyCompressibleFluidMatPars(mat_id_fluid);
    auto fparams_struc = GetSVKMatPars(mat_id_struc);

    return Teuchos::rcp(
        new FLD::WeaklyCompressibleEtienneFSIFluidFunction(fparams_fluid, fparams_struc));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE"))
  {
    // read data
    int mat_id_fluid = -1;
    int mat_id_struc = -1;

    function_lin_def->ExtractInt("MAT_FLUID", mat_id_fluid);
    function_lin_def->ExtractInt("MAT_STRUC", mat_id_struc);

    if (mat_id_fluid <= 0)
    {
      dserror(
          "Please give a (reasonable) 'MAT_FLUID' in "
          "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE");
    }
    if (mat_id_struc <= 0)
    {
      dserror(
          "Please give a (reasonable) 'MAT_STRUC' in "
          "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE");
    }

    // get materials
    auto fparams_fluid = GetWeaklyCompressibleFluidMatPars(mat_id_fluid);
    auto fparams_struc = GetSVKMatPars(mat_id_struc);

    return Teuchos::rcp(
        new FLD::WeaklyCompressibleEtienneFSIFluidForceFunction(fparams_fluid, fparams_struc));
  }
  else if (function_lin_def->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY"))
  {
    // read data
    int mat_id_fluid = -1;
    int mat_id_struc = -1;

    function_lin_def->ExtractInt("MAT_FLUID", mat_id_fluid);
    function_lin_def->ExtractInt("MAT_STRUC", mat_id_struc);

    if (mat_id_fluid <= 0)
    {
      dserror(
          "Please give a (reasonable) 'MAT_FLUID' in "
          "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY");
    }
    if (mat_id_struc <= 0)
    {
      dserror(
          "Please give a (reasonable) 'MAT_STRUC' in "
          "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY");
    }

    // get materials
    auto fparams_fluid = GetWeaklyCompressibleFluidMatPars(mat_id_fluid);
    auto fparams_struc = GetSVKMatPars(mat_id_struc);

    return Teuchos::rcp(
        new FLD::WeaklyCompressibleEtienneFSIFluidViscosityFunction(fparams_fluid, fparams_struc));
  }
  else if (function_lin_def->HaveNamed("BELTRAMI-UP"))
  {
    // read material
    int mat_id = -1;

    function_lin_def->ExtractInt("MAT", mat_id);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-UP");

    // get material
    auto fparams = GetNewtonianFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::BeltramiUP(fparams));
  }
  else if (function_lin_def->HaveNamed("BELTRAMI-RHS"))
  {
    // read material
    int mat_id = -1;
    int is_stokes = 0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractInt("ISSTOKES", is_stokes);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-RHS");

    // get material
    auto fparams = GetNewtonianFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::BeltramiRHS(fparams, (bool)is_stokes));
  }
  else if (function_lin_def->HaveNamed("KIMMOIN-UP"))
  {
    // read material
    int mat_id = -1;
    int is_stationary = 0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractInt("ISSTAT", is_stationary);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-UP");

    // get material
    auto fparams = GetNewtonianFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::KimMoinUP(fparams, (bool)is_stationary));
  }
  else if (function_lin_def->HaveNamed("KIMMOIN-RHS"))
  {
    // read material
    int mat_id = -1;
    int is_stationary = 0;
    int is_stokes = 0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractInt("ISSTAT", is_stationary);
    function_lin_def->ExtractInt("ISSTOKES", is_stokes);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-RHS");

    // get material
    auto fparams = GetNewtonianFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::KimMoinRHS(fparams, (bool)is_stationary, (bool)is_stokes));
  }
  else if (function_lin_def->HaveNamed("KIMMOIN-STRESS"))
  {
    // read material
    int mat_id = -1;
    int is_stationary = 0;
    double amplitude = 1.0;

    function_lin_def->ExtractInt("MAT", mat_id);
    function_lin_def->ExtractInt("ISSTAT", is_stationary);
    function_lin_def->ExtractDouble("AMPLITUDE", amplitude);

    if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-STRESS");

    // get material
    auto fparams = GetNewtonianFluidMatPars(mat_id);

    return Teuchos::rcp(new FLD::KimMoinStress(fparams, (bool)is_stationary, amplitude));
  }
  else
  {
    return Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime>(NULL);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiFunction::BeltramiFunction(double c1) : c1_(c1) {}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  double a = M_PI / 4.0;
  double d = M_PI / 2.0;

  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id == -1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const auto* actmat = dynamic_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double dynvisc = actmat->viscosity_;
  double kinvisc = dynvisc / dens;
  double tempfac = exp(-c1_ * kinvisc * d * d * t);

  switch (component)
  {
    case 0:
      return -a *
             (exp(a * xp[0]) * sin(a * xp[1] + d * xp[2]) +
                 exp(a * xp[2]) * cos(a * xp[0] + d * xp[1])) *
             tempfac;
    case 1:
      return -a *
             (exp(a * xp[1]) * sin(a * xp[2] + d * xp[0]) +
                 exp(a * xp[0]) * cos(a * xp[1] + d * xp[2])) *
             tempfac;
    case 2:
      return -a *
             (exp(a * xp[2]) * sin(a * xp[0] + d * xp[1]) +
                 exp(a * xp[1]) * cos(a * xp[2] + d * xp[0])) *
             tempfac;
    case 3:
      return -a * a / 2 * dens *
             (exp(2 * a * xp[0]) + exp(2 * a * xp[1]) + exp(2 * a * xp[2]) +
                 2 * sin(a * xp[0] + d * xp[1]) * cos(a * xp[2] + d * xp[0]) *
                     exp(a * (xp[1] + xp[2])) +
                 2 * sin(a * xp[1] + d * xp[2]) * cos(a * xp[0] + d * xp[1]) *
                     exp(a * (xp[2] + xp[0])) +
                 2 * sin(a * xp[2] + d * xp[0]) * cos(a * xp[1] + d * xp[2]) *
                     exp(a * (xp[0] + xp[1]))) *
             tempfac;
    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  double a = M_PI / 4.0;
  double d = M_PI / 2.0;

  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
  if (id == -1) dserror("Newtonian fluid material could not be found");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const auto* actmat = dynamic_cast<const MAT::PAR::NewtonianFluid*>(mat);
  double dens = actmat->density_;
  double dynvisc = actmat->viscosity_;
  double kinvisc = dynvisc / dens;
  double der1tempfac = (-c1_ * kinvisc * d * d) * exp(-c1_ * kinvisc * d * d * t);
  double der2tempfac =
      (-c1_ * kinvisc * d * d) * (-c1_ * kinvisc * d * d) * exp(-c1_ * kinvisc * d * d * t);

  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (component)
    {
      case 0:
        res[1] = -a *
                 (exp(a * xp[0]) * sin(a * xp[1] + d * xp[2]) +
                     exp(a * xp[2]) * cos(a * xp[0] + d * xp[1])) *
                 der1tempfac;
      case 1:
        res[1] = -a *
                 (exp(a * xp[1]) * sin(a * xp[2] + d * xp[0]) +
                     exp(a * xp[0]) * cos(a * xp[1] + d * xp[2])) *
                 der1tempfac;
      case 2:
        res[1] = -a *
                 (exp(a * xp[2]) * sin(a * xp[0] + d * xp[1]) +
                     exp(a * xp[1]) * cos(a * xp[2] + d * xp[0])) *
                 der1tempfac;
      case 3:
        res[1] = -a * a / 2 * dens *
                 (exp(2 * a * xp[0]) + exp(2 * a * xp[1]) + exp(2 * a * xp[2]) +
                     2 * sin(a * xp[0] + d * xp[1]) * cos(a * xp[2] + d * xp[0]) *
                         exp(a * (xp[1] + xp[2])) +
                     2 * sin(a * xp[1] + d * xp[2]) * cos(a * xp[0] + d * xp[1]) *
                         exp(a * (xp[2] + xp[0])) +
                     2 * sin(a * xp[2] + d * xp[0]) * cos(a * xp[1] + d * xp[2]) *
                         exp(a * (xp[0] + xp[1]))) *
                 der1tempfac;
      default:
        res[1] = 1.0;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (component)
    {
      case 0:
        res[2] = -a *
                 (exp(a * xp[0]) * sin(a * xp[1] + d * xp[2]) +
                     exp(a * xp[2]) * cos(a * xp[0] + d * xp[1])) *
                 der2tempfac;
      case 1:
        res[2] = -a *
                 (exp(a * xp[1]) * sin(a * xp[2] + d * xp[0]) +
                     exp(a * xp[0]) * cos(a * xp[1] + d * xp[2])) *
                 der2tempfac;
      case 2:
        res[2] = -a *
                 (exp(a * xp[2]) * sin(a * xp[0] + d * xp[1]) +
                     exp(a * xp[1]) * cos(a * xp[2] + d * xp[0])) *
                 der2tempfac;
      case 3:
        res[2] = -a * a / 2 * dens *
                 (exp(2 * a * xp[0]) + exp(2 * a * xp[1]) + exp(2 * a * xp[2]) +
                     2 * sin(a * xp[0] + d * xp[1]) * cos(a * xp[2] + d * xp[0]) *
                         exp(a * (xp[1] + xp[2])) +
                     2 * sin(a * xp[1] + d * xp[2]) * cos(a * xp[0] + d * xp[1]) *
                         exp(a * (xp[2] + xp[0])) +
                     2 * sin(a * xp[2] + d * xp[0]) * cos(a * xp[1] + d * xp[2]) *
                         exp(a * (xp[0] + xp[1]))) *
                 der2tempfac;
      default:
        res[2] = 1.0;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::ChannelWeaklyCompressibleFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid_murnaghantait);
  if (id == -1)
  {
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(
        INPAR::MAT::m_fluid_linear_density_viscosity);
    if (id == -1)
      dserror(
          "Fluid with Murnaghan-Tait equation of state or with "
          "linear law (pressure-dependent) for the density and the viscosity could not be found");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const auto* actmat = dynamic_cast<const MAT::PAR::LinearDensityViscosity*>(mat);

    double x = xp[0];
    double y = xp[1];

    double length = 10.0;
    double radius = 1.0;
    double aspect_ratio = radius / length;
    double mean_velocity_channel_exit = 1.0;
    double reference_viscosity = actmat->refviscosity_;
    double reference_pressure = actmat->refpressure_;
    double coefficient_density = actmat->coeffdensity_;
    double coefficient_viscosity = actmat->coeffviscosity_;
    double coefficient_density_adim =
        (3.0 * coefficient_density * reference_viscosity * length * mean_velocity_channel_exit) /
        std::pow(radius, 2.0);
    double coefficient_viscosity_adim =
        (3.0 * coefficient_viscosity * reference_viscosity * length * mean_velocity_channel_exit) /
        std::pow(radius, 2.0);

    // parameters according with the paper
    double z = x / length;
    double r = y / radius;
    double alfa = aspect_ratio;
    double beta = coefficient_viscosity_adim;
    double epsilon = coefficient_density_adim;
    double a = aspect_ratio;
    double B = alfa * beta;
    double lambda = 1.0 + (1.0 / 5.0) * std::pow(B, 2.0) + (11.0 / 175.0) * std::pow(B, 4.0) +
                    (533.0 / 23625.0) * std::pow(B, 6.0) + (5231.0 / 606375.0) * std::pow(B, 8.0);
    double p_0_hat = std::cosh(alfa * beta * lambda * r) / std::cosh(alfa * beta * lambda);
    double u_r1_hat =
        -(11.0 * r * std::pow(1.0 - std::pow(r, 2.0), 2.0)) / 40.0 * std::pow(B, 2.0) *
        (1.0 + ((173.0 - 85.0 * std::pow(r, 2.0)) / (770.0)) * std::pow(B, 2.0) +
            ((5793.0 - 7190.0 * std::pow(r, 2.0) + 3965.0 * std::pow(r, 4.0)) / (83160.0)) *
                std::pow(B, 4.0) +
            ((7435723.0 - 16839665.0 * std::pow(r, 2.0) + 16836225.0 * std::pow(r, 4.0) -
                 5021275.0 * std::pow(r, 6.0)) /
                (320166000.0)) *
                std::pow(B, 6.0));
    double u_r1_hat_first =
        (11.0 * std::pow(B, 2.0) * std::pow(std::pow(r, 2.0) - 1.0, 2.0) *
            (((4099.0 * std::pow(r, 6.0)) / 261360.0 - (32069.0 * std::pow(r, 4.0)) / 609840.0 +
                 (3367933.0 * std::pow(r, 2.0)) / 64033200.0 - 7435723.0 / 320166000.0) *
                    std::pow(B, 6.0) +
                (-(793.0 * std::pow(r, 4.0)) / 16632.0 + (719.0 * std::pow(r, 2.0)) / 8316.0 -
                    1931.0 / 27720.0) *
                    std::pow(B, 4.0) +
                ((17.0 * std::pow(r, 2.0)) / 154.0 - 173.0 / 770.0) * std::pow(B, 2.0) - 1.0)) /
            40.0 +
        (11.0 * std::pow(B, 2.0) * std::pow(r, 2.0) * (std::pow(r, 2.0) - 1.0) *
            (((4099.0 * std::pow(r, 6.0)) / 261360.0 - (32069.0 * std::pow(r, 4.0)) / 609840.0 +
                 (3367933.0 * std::pow(r, 2.0)) / 64033200.0 - 7435723.0 / 320166000.0) *
                    std::pow(B, 6.0) +
                (-(793.0 * std::pow(r, 4.0)) / 16632.0 + (719.0 * std::pow(r, 2.0)) / 8316.0 -
                    1931.0 / 27720.0) *
                    std::pow(B, 4.0) +
                ((17.0 * std::pow(r, 2.0)) / 154.0 - 173.0 / 770.0) * std::pow(B, 2.0) - 1.0)) /
            10.0 +
        (11.0 * std::pow(B, 2.0) * r * std::pow(std::pow(r, 2.0) - 1.0, 2.0) *
            (((4099.0 * std::pow(r, 5.0)) / 43560.0 - (32069.0 * std::pow(r, 3.0)) / 152460.0 +
                 (3367933.0 * r) / 32016600.0) *
                    std::pow(B, 6.0) +
                ((719.0 * r) / 4158.0 - (793.0 * std::pow(r, 3.0)) / 4158.0) * std::pow(B, 4.0) +
                (17.0 * r * std::pow(B, 2.0)) / 77.0)) /
            40.0;
    double h = 1.0 / std::pow(beta, 2.0) *
               (-1.0 + ((11.0 - 10.0 * std::pow(r, 2.0)) / (15.0)) * std::pow(B, 2.0) +
                   ((359.0 - 126.0 * std::pow(r, 2.0) + 35.0 * std::pow(r, 4.0)) / (1260.0)) *
                       std::pow(B, 4.0) +
                   ((13761.0 - 17790.0 * std::pow(r, 2.0) + 34125.0 * std::pow(r, 4.0) -
                        17500.0 * std::pow(r, 6.0)) /
                       (94500.0)) *
                       std::pow(B, 6.0) +
                   ((225311.0 - 614515.0 * std::pow(r, 2.0) + 1492755.0 * std::pow(r, 4.0) -
                        1324785.0 * std::pow(r, 6.0) + 394350.0 * std::pow(r, 8.0)) /
                       (3118500.0)) *
                       std::pow(B, 8.0));
    double h_1 = 1.0 / std::pow(beta, 2.0) *
                 (-1.0 + ((11.0 - 10.0) / (15.0)) * std::pow(B, 2.0) +
                     ((359.0 - 126.0 + 35.0) / (1260.0)) * std::pow(B, 4.0) +
                     ((13761.0 - 17790.0 + 34125.0 - 17500.0) / (94500.0)) * std::pow(B, 6.0) +
                     ((225311.0 - 614515.0 + 1492755.0 - 1324785.0 + 394350.0) / (3118500.0)) *
                         std::pow(B, 8.0));

    CORE::LINALG::Matrix<2, 1> u;
    double p;

    u(0) = -(3.0 * std::log(p_0_hat)) / (std::pow(B, 2.0) * lambda) +
           epsilon * ((3 * (std::tanh(B * lambda) - r * std::tanh(B * lambda * r)) +
                          std::log(std::pow(p_0_hat, 3.0)) / (B * lambda)) /
                             (beta * (3 * std::tanh(B * lambda) - 2 * B)) +
                         (std::exp(lambda * beta * (1.0 - z))) / (lambda * beta) *
                             ((p_0_hat * std::log(std::pow(p_0_hat, 3.0))) / (std::pow(B, 2.0)) +
                                 u_r1_hat_first));
    u(1) = epsilon * u_r1_hat * std::exp(lambda * beta * (1.0 - z));
    p = (p_0_hat * std::exp(lambda * beta * (1.0 - z)) - 1.0) / beta +
        epsilon * p_0_hat * std::exp(lambda * beta * (1.0 - z)) *
            ((lambda * a *
                 (1.0 - z + a * (r * std::tanh(B * lambda * r) - std::tanh(B * lambda)))) /
                    (3.0 * std::tanh(B * lambda) - 2.0 * B) +
                p_0_hat * h * std::exp(lambda * beta * (1.0 - z)) - h_1);

    switch (component)
    {
      case 0:
        return u(0) * mean_velocity_channel_exit;
      case 1:
        return u(1) * mean_velocity_channel_exit * radius / length;
      case 2:
        return p * (3.0 * reference_viscosity * length * mean_velocity_channel_exit /
                       std::pow(radius, 2.0)) +
               reference_pressure;

      default:
        return 1.0;
    }
  }
  else
  {
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const auto* actmat = dynamic_cast<const MAT::PAR::MurnaghanTaitFluid*>(mat);

    if (actmat->matparameter_ != 1.0)
    {
      dserror("The analytical solution is only valid for material parameter = 1");
    }

    double x = xp[0];
    double y = xp[1];

    double length = 10.0;
    double radius = 1.0;
    double mean_velocity_channel_exit = 1.0;
    double aspect_ratio = radius / length;
    double viscosity = actmat->viscosity_;
    double reference_pressure = actmat->refpressure_;
    double reference_bulk_modulus = actmat->refbulkmodulus_;
    double linear_coefficient_density =
        (3.0 * (1.0 / reference_bulk_modulus) * viscosity * length * mean_velocity_channel_exit) /
        std::pow(radius, 2.0);

    CORE::LINALG::Matrix<2, 1> u;
    double p;
    CORE::LINALG::Matrix<2, 2> dervel;

    u(0) = 3.0 / 2.0 * (1.0 - std::pow(y / radius, 2.0)) *
           (1.0 + linear_coefficient_density * (x / length - 1.0));
    u(1) = 0.0;
    p = 1.0 - x / length -
        linear_coefficient_density *
            (1.0 / 6.0 * std::pow(aspect_ratio, 2.0) * (std::pow(y / radius, 2.0) - 1.0) +
                1.0 / 2.0 * std::pow(1.0 - x / length, 2.0));

    dervel(0, 0) =
        3.0 / 2.0 * (1.0 - std::pow(y / radius, 2.0)) * linear_coefficient_density / length;
    dervel(0, 1) =
        -3.0 / std::pow(radius, 2.0) * y * (1.0 + linear_coefficient_density * (x / length - 1.0));
    dervel(1, 0) = 0.0;
    dervel(1, 1) = 0.0;

    // scaling correctly the variables
    u(0) = u(0) * mean_velocity_channel_exit;
    u(1) = u(1) * mean_velocity_channel_exit * radius / length;
    p = p * (3.0 * viscosity * length * mean_velocity_channel_exit / std::pow(radius, 2.0)) +
        reference_pressure;
    dervel(0, 0) = dervel(0, 0) * mean_velocity_channel_exit;
    dervel(0, 1) = dervel(0, 1) * mean_velocity_channel_exit;
    dervel(1, 0) = dervel(1, 0) * mean_velocity_channel_exit * radius / length;
    dervel(1, 1) = dervel(1, 1) * mean_velocity_channel_exit * radius / length;

    switch (component)
    {
      case 0:
        return u(0);
      case 1:
        return u(1);
      case 2:
        return p;
      case 3:
        return dervel(0, 0);
      case 4:
        return dervel(0, 1);
      case 5:
        return dervel(1, 0);
      case 6:
        return dervel(1, 1);

      default:
        return 1.0;
    }
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::ChannelWeaklyCompressibleFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::CorrectionTermChannelWeaklyCompressibleFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid_murnaghantait);
  if (id == -1)
  {
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(
        INPAR::MAT::m_fluid_linear_density_viscosity);
    if (id == -1)
      dserror(
          "Fluid with Murnaghan-Tait equation of state or with "
          "linear law (pressure-dependent) for the density and the viscosity could not be found");
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const auto* actmat = dynamic_cast<const MAT::PAR::LinearDensityViscosity*>(mat);

    double x = xp[0];
    double y = xp[1];

    if (actmat->coeffviscosity_ > 1.0e-5)
    {
      dserror("The correction term is only valid for viscosity coefficient -> 0");
    }

    double L = 10.0;
    double R = 1.0;
    double U = 1.0;
    double mu = actmat->refviscosity_;
    double K0 = 1.0 / actmat->coeffdensity_;

    double Corrterm =
        ((162.0 * std::pow(R, 4.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) * x -
             162.0 * L * std::pow(R, 4.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) +
             162.0 * L * std::pow(R, 2.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) *
                 std::pow(y, 2.0) -
             162.0 * std::pow(R, 2.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) * x *
                 std::pow(y, 2.0)) *
                K0 +
            243.0 * std::pow(L, 2.0) * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) -
            243.0 * std::pow(L, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 2.0) -
            486.0 * L * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * x +
            486.0 * L * std::pow(U, 4.0) * std::pow(mu, 3.0) * x * std::pow(y, 2.0) -
            27.0 * std::pow(R, 4.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) +
            243.0 * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(x, 2.0) +
            54.0 * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 2.0) -
            243.0 * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(x, 2.0) * std::pow(y, 2.0) -
            27.0 * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 4.0)) /
        (-4.0 * std::pow(R, 8.0) * std::pow(K0, 3.0) +
            (12.0 * std::pow(R, 6.0) * U * mu * x - 12.0 * L * std::pow(R, 6.0) * U * mu) *
                std::pow(K0, 2.0) +
            (18.0 * std::pow(L, 2.0) * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) -
                36.0 * L * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * x -
                6.0 * std::pow(R, 6.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) +
                18.0 * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * std::pow(x, 2.0) +
                6.0 * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * std::pow(y, 2.0)) *
                K0);

    switch (component)
    {
      case 0:
        return Corrterm;

      default:
        return 0.0;
    }
  }
  else
  {
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const auto* actmat = dynamic_cast<const MAT::PAR::MurnaghanTaitFluid*>(mat);

    if (actmat->matparameter_ != 1.0)
    {
      dserror("The correction term is only valid for material parameter = 1");
    }

    double x = xp[0];
    double y = xp[1];

    double L = 10.0;
    double R = 1.0;
    double U = 1.0;
    double mu = actmat->viscosity_;
    double K0 = actmat->refbulkmodulus_;

    double Corrterm =
        ((162.0 * std::pow(R, 4.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) * x -
             162.0 * L * std::pow(R, 4.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) +
             162.0 * L * std::pow(R, 2.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) *
                 std::pow(y, 2.0) -
             162.0 * std::pow(R, 2.0) * std::pow(U, 3.0) * std::pow(mu, 2.0) * x *
                 std::pow(y, 2.0)) *
                K0 +
            243.0 * std::pow(L, 2.0) * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) -
            243.0 * std::pow(L, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 2.0) -
            486.0 * L * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * x +
            486.0 * L * std::pow(U, 4.0) * std::pow(mu, 3.0) * x * std::pow(y, 2.0) -
            27.0 * std::pow(R, 4.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) +
            243.0 * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(x, 2.0) +
            54.0 * std::pow(R, 2.0) * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 2.0) -
            243.0 * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(x, 2.0) * std::pow(y, 2.0) -
            27.0 * std::pow(U, 4.0) * std::pow(mu, 3.0) * std::pow(y, 4.0)) /
        (-4.0 * std::pow(R, 8.0) * std::pow(K0, 3.0) +
            (12.0 * std::pow(R, 6.0) * U * mu * x - 12.0 * L * std::pow(R, 6.0) * U * mu) *
                std::pow(K0, 2.0) +
            (18.0 * std::pow(L, 2.0) * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) -
                36.0 * L * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * x -
                6.0 * std::pow(R, 6.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) +
                18.0 * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * std::pow(x, 2.0) +
                6.0 * std::pow(R, 4.0) * std::pow(U, 2.0) * std::pow(mu, 2.0) * std::pow(y, 2.0)) *
                K0);

    switch (component)
    {
      case 0:
        return Corrterm;

      default:
        return 0.0;
    }
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::CorrectionTermChannelWeaklyCompressibleFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressiblePoiseuilleFunction::WeaklyCompressiblePoiseuilleFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U)
    : length_(L),
      halfheight_(R),
      meanvelocityexit_(U),
      viscosity_(fparams.viscosity_),
      refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressiblePoiseuilleFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double L = length_;
  double R = halfheight_;
  double U = meanvelocityexit_;
  double mu = viscosity_;
  double r0 = refdensity_;
  double p0 = refpressure_;
  double epsilon = comprcoeff_;

  // initialize variables
  CORE::LINALG::Matrix<3, 1> L_ex;
  double p_ex;
  CORE::LINALG::Matrix<2, 1> v_ex;
  double r_ex;
  CORE::LINALG::Matrix<2, 1> w_ex;

  // evaluate variables
  L_ex(0) = -sqrt(2. * mu) / 2. *
            (+4. / 3. * (9. / 2. * mu * pow(U, 2.) / pow(R, 2.) * (1. - pow(y / R, 2.)) * epsilon));
  L_ex(1) = -sqrt(2. * mu) / 2. *
            (-2. / 3. * (9. / 2. * mu * pow(U, 2.) / pow(R, 2.) * (1. - pow(y / R, 2.)) * epsilon));
  L_ex(2) = -sqrt(1. * mu) * (-3. * U / pow(R, 2.) * y + 9. * mu * L * pow(U, 2.) / pow(R, 4.) *
                                                             (1. - x / L) * y * epsilon);
  p_ex = p0 + 3. * mu * L * U / pow(R, 2.) * (1. - x / L) -
         3. / 2. * pow(mu * U / R, 2.) *
             (3. * pow(L / R, 2.) * pow(1. - x / L, 2.) - (1. - pow(y / R, 2.))) * epsilon;
  v_ex(0) = 3. / 2. * U * (1. - pow(y / R, 2.)) - 9. / 2. * mu * L * pow(U, 2.) / pow(R, 2.) *
                                                      (1. - x / L) * (1. - pow(y / R, 2.)) *
                                                      epsilon;
  v_ex(1) = 0.;

  // density - momentum formulation
  r_ex = r0 + epsilon * (p_ex - p0);
  w_ex(0) = r_ex * v_ex(0);
  w_ex(1) = r_ex * v_ex(1);

  switch (component)
  {
    case 0:
      return r_ex;
    case 1:
      return w_ex(0);
    case 2:
      return w_ex(1);
    case 3:
      return L_ex(0);
    case 4:
      return L_ex(1);
    case 5:
      return L_ex(2);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressiblePoiseuilleFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressiblePoiseuilleForceFunction::WeaklyCompressiblePoiseuilleForceFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams, double L, double R, double U)
    : length_(L),
      halfheight_(R),
      meanvelocityexit_(U),
      viscosity_(fparams.viscosity_),
      refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressiblePoiseuilleForceFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double L = length_;
  double R = halfheight_;
  double U = meanvelocityexit_;
  double mu = viscosity_;
  double r0 = refdensity_;
  double epsilon = comprcoeff_;

  // initialize variables
  double f_r_ex;
  CORE::LINALG::Matrix<2, 1> f_w_ex;

  // evaluate variables
  f_r_ex =
      -(9. * pow(U, 2.) * epsilon * mu * (pow(R, 2.) - pow(y, 2.)) *
          (2. * pow(R, 4.) - 2. * pow(R, 4.) * r0 +
              27. * pow(L, 2.) * pow(U, 2.) * pow(epsilon, 2.) * pow(mu, 2.) -
              3. * pow(R, 2.) * pow(U, 2.) * pow(epsilon, 2.) * pow(mu, 2.) +
              27. * pow(U, 2.) * pow(epsilon, 2.) * pow(mu, 2.) * pow(x, 2.) +
              3. * pow(U, 2.) * pow(epsilon, 2.) * pow(mu, 2.) * pow(y, 2.) -
              54. * L * pow(U, 2.) * pow(epsilon, 2.) * pow(mu, 2.) * x -
              18. * L * pow(R, 2.) * U * epsilon * mu + 18. * pow(R, 2.) * U * epsilon * mu * x)) /
      (4. * pow(R, 8.));
  f_w_ex(0) = 0.0;
  f_w_ex(1) = 0.0;

  switch (component)
  {
    case 0:
      return f_r_ex;
    case 1:
      return f_w_ex(0);
    case 2:
      return f_w_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressiblePoiseuilleForceFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleManufacturedFlowFunction::WeaklyCompressibleManufacturedFlowFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams)
    : viscosity_(fparams.viscosity_),
      refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleManufacturedFlowFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double mu = viscosity_;
  double r0 = refdensity_;
  double p0 = refpressure_;
  double epsilon = comprcoeff_;

  // initialize variables
  CORE::LINALG::Matrix<3, 1> L_ex;
  double p_ex;
  CORE::LINALG::Matrix<2, 1> v_ex;
  double r_ex;
  CORE::LINALG::Matrix<2, 1> w_ex;

  // evaluate variables
  L_ex(0) =
      -sqrt(2. * mu) / 2. *
      (+4. / 3. *
              (M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y) -
                  epsilon *
                      ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0))) -
          2. / 3. *
              (-epsilon *
                      ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0)) -
                  M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)));
  L_ex(1) =
      -sqrt(2. * mu) / 2. *
      (-2. / 3. *
              (M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y) -
                  epsilon *
                      ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0))) +
          4. / 3. *
              (-epsilon *
                      ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0)) -
                  M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)));
  L_ex(2) =
      -sqrt(1. * mu) *
      ((M_PI * cos(M_PI * y) * sin(M_PI * t) * sin(M_PI * x) -
           epsilon *
               (((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * y) * sin(M_PI * x)) / (2. * r0) -
                   ((std::pow(M_PI, 3.)) * x * (std::pow((sin(M_PI * t)), 2.)) *
                       sin(2. * M_PI * y)) /
                       (2. * r0))) +
          (-epsilon * (((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * y) * sin(M_PI * x)) /
                              (2. * r0) -
                          ((std::pow(M_PI, 3.)) * y * (std::pow((sin(M_PI * t)), 2.)) *
                              sin(2. * M_PI * x)) /
                              (2. * r0)) -
              M_PI * cos(M_PI * y) * sin(M_PI * t) * sin(M_PI * x)));
  p_ex = M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y);
  v_ex(0) = sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y) -
            epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0));
  v_ex(1) = cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t) -
            epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                              (8. * r0) -
                          (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0));

  // density - momentum formulation
  r_ex = r0 + epsilon * (p_ex - p0);
  w_ex(0) = r_ex * v_ex(0);
  w_ex(1) = r_ex * v_ex(1);

  switch (component)
  {
    case 0:
      return r_ex;
    case 1:
      return w_ex(0);
    case 2:
      return w_ex(1);
    case 3:
      return L_ex(0);
    case 4:
      return L_ex(1);
    case 5:
      return L_ex(2);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleManufacturedFlowFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleManufacturedFlowForceFunction::
    WeaklyCompressibleManufacturedFlowForceFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams)
    : viscosity_(fparams.viscosity_),
      refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleManufacturedFlowForceFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double mu = viscosity_;
  double r0 = refdensity_;
  double p0 = refpressure_;
  double epsilon = comprcoeff_;

  // initialize variables
  double f_r_ex;
  CORE::LINALG::Matrix<2, 1> f_w_ex;

  // evaluate variables
  f_r_ex = (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y) -
           (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0)) -
               M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
               (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) -
           (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                           (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                              (8. * r0) +
                          ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                              (2. * r0)) +
               M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
               (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) +
           (std::pow(M_PI, 2.)) * epsilon * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y) *
               (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                               (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                                  (8. * r0) +
                              (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
                   sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) -
           (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t) *
               (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                               (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                                  (8. * r0) -
                              (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
                   cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t));
  f_w_ex(0) =
      (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                      (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                         (8. * r0) +
                     (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
          sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                             (2. * r0)) +
              M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) -
      mu * (epsilon *
                   (((std::pow(M_PI, 3.)) * (std::pow((sin(M_PI * t)), 2.)) * sin(2. * M_PI * x)) /
                           (2. * r0) +
                       ((std::pow(M_PI, 3.)) * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) /
                           (2. * r0)) +
               epsilon * (((std::pow(M_PI, 3.)) * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) /
                                 (2. * r0) +
                             ((std::pow(M_PI, 4.)) * x * cos(2. * M_PI * y) *
                                 (std::pow((sin(M_PI * t)), 2.))) /
                                 r0) +
               (2. * epsilon *
                   (((std::pow(M_PI, 3.)) * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) /
                           (2. * r0) +
                       ((std::pow(M_PI, 3.)) * cos(M_PI * x) * (std::pow((sin(M_PI * t)), 2.)) *
                           sin(M_PI * x)) /
                           r0)) /
                   3. -
               2. * (std::pow(M_PI, 2.)) * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) -
      (epsilon * (((std::pow(M_PI, 2.)) * cos(M_PI * t) * sin(M_PI * t) *
                      (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                         (4. * r0) -
                     ((std::pow(M_PI, 2.)) * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) /
                         (2. * r0)) -
          M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) +
      2. *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
              sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                             (2. * r0)) -
              M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) +
      (epsilon *
              (((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * y) * sin(M_PI * x)) / (2. * r0) -
                  ((std::pow(M_PI, 3.)) * x * (std::pow((sin(M_PI * t)), 2.)) *
                      sin(2. * M_PI * y)) /
                      (2. * r0)) -
          M_PI * cos(M_PI * y) * sin(M_PI * t) * sin(M_PI * x)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) -
      (std::pow(M_PI, 2.)) * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y) -
      (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
              sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) -
      (std::pow(M_PI, 2.)) * epsilon * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y) *
          (std::pow(
              (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                              (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                                 (8. * r0) +
                             (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
                  sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)),
              2.)) +
      (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
              sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t));
  f_w_ex(1) =
      mu * (epsilon * (((std::pow(M_PI, 3.)) * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) /
                              (2. * r0) -
                          ((std::pow(M_PI, 4.)) * y * cos(2. * M_PI * x) *
                              (std::pow((sin(M_PI * t)), 2.))) /
                              r0) -
               epsilon *
                   (((std::pow(M_PI, 3.)) * (std::pow((sin(M_PI * t)), 2.)) * sin(2. * M_PI * y)) /
                           (2. * r0) -
                       ((std::pow(M_PI, 3.)) * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) /
                           (2. * r0)) +
               (2. * epsilon *
                   (((std::pow(M_PI, 3.)) * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) /
                           (2. * r0) -
                       ((std::pow(M_PI, 3.)) * cos(M_PI * y) * (std::pow((sin(M_PI * t)), 2.)) *
                           sin(M_PI * y)) /
                           r0)) /
                   3. +
               2. * (std::pow(M_PI, 2.)) * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) -
      (epsilon *
              (((std::pow(M_PI, 2.)) * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) / (2. * r0) +
                  ((std::pow(M_PI, 2.)) * cos(M_PI * t) * sin(M_PI * t) *
                      (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                      (4. * r0)) -
          M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) +
      (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                      (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                         (8. * r0) +
                     (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
          sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (epsilon * (((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * y) * sin(M_PI * x)) /
                             (2. * r0) -
                         ((std::pow(M_PI, 3.)) * y * (std::pow((sin(M_PI * t)), 2.)) *
                             sin(2. * M_PI * x)) /
                             (2. * r0)) +
              M_PI * cos(M_PI * y) * sin(M_PI * t) * sin(M_PI * x)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) +
      2. *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                             (2. * r0)) +
              M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) +
      (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                      (2. * M_PI * cos(2. * M_PI * x) + 2. * M_PI * cos(2. * M_PI * y))) /
                         (8. * r0) +
                     ((std::pow(M_PI, 2.)) * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y)) /
                         (2. * r0)) -
          M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y)) *
          (r0 - epsilon * (p0 - M_PI * cos(M_PI * x) * sin(M_PI * t) * sin(M_PI * y))) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) +
      (std::pow(M_PI, 2.)) * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t) +
      (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t) *
          (std::pow(
              (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                              (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                                 (8. * r0) -
                             (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
                  cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)),
              2.)) -
      (std::pow(M_PI, 2.)) * epsilon * cos(M_PI * t) * cos(M_PI * x) * sin(M_PI * y) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t)) -
      (std::pow(M_PI, 2.)) * epsilon * sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * x) + 2. * M_PI * x * cos(2. * M_PI * y))) /
                             (8. * r0) +
                         (M_PI * cos(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) / (2. * r0)) -
              sin(M_PI * t) * sin(M_PI * x) * sin(M_PI * y)) *
          (epsilon * ((M_PI * (std::pow((sin(M_PI * t)), 2.)) *
                          (sin(2. * M_PI * y) + 2. * M_PI * y * cos(2. * M_PI * x))) /
                             (8. * r0) -
                         (M_PI * cos(M_PI * t) * cos(M_PI * x) * cos(M_PI * y)) / (2. * r0)) -
              cos(M_PI * x) * cos(M_PI * y) * sin(M_PI * t));

  switch (component)
  {
    case 0:
      return f_r_ex;
    case 1:
      return f_w_ex(0);
    case 2:
      return f_w_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleManufacturedFlowForceFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneCFDFunction::WeaklyCompressibleEtienneCFDFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams)
    : refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneCFDFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double r0 = refdensity_;
  double p0 = refpressure_;
  double epsilon = comprcoeff_;

  // initialize variables
  CORE::LINALG::Matrix<3, 1> L_ex;
  double p_ex;
  CORE::LINALG::Matrix<2, 1> v_ex;
  double r_ex;
  CORE::LINALG::Matrix<2, 1> w_ex;

  // evaluate variables
  L_ex(0) = +10. * x * (std::pow(y, 3.)) *
            (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
            (std::pow((-10. * (std::pow(x, 5.)) + 25. * (std::pow(x, 4.)) -
                          20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                4.)) *
            (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
                20. * (std::pow(x, 2.)) + 5. * y - 4.);
  L_ex(1) = -10. * x * (std::pow(y, 3.)) *
            (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
            (std::pow((-10. * (std::pow(x, 5.)) + 25. * (std::pow(x, 4.)) -
                          20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                4.)) *
            (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
                20. * (std::pow(x, 2.)) + 5. * y - 4.);
  L_ex(2) =
      20. * (std::pow(y, 3.)) *
          ((std::pow(
               (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 5.)) /
                  5. +
              (std::pow(y, 5.)) / 5.) +
      12. * (std::pow(y, 2.)) *
          ((std::pow(
               (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 6.)) /
                  6. -
              (std::pow(y, 6.)) / 6.) +
      (std::pow(y, 8.)) -
      (std::pow(y, 4.)) *
          (std::pow(
              (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 4.)) *
          (std::pow((10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) +
                        10. * (std::pow(x, 2.)) * (2. * x - 2.) * (x - 1. / 2.) +
                        20. * x * (std::pow((x - 1.), 2.)) * (x - 1. / 2.)),
              2.)) -
      (std::pow(y, 4.)) *
          (std::pow(
              (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 4.)) *
          (y + 10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.) *
          (40. * x * (std::pow((x - 1.), 2.)) + 20. * (std::pow(x, 2.)) * (x - 1. / 2.) +
              20. * (std::pow(x, 2.)) * (2. * x - 2.) +
              20. * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) +
              40. * x * (2. * x - 2.) * (x - 1. / 2.)) -
      4. * (std::pow(y, 4.)) *
          (std::pow(
              (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 3.)) *
          (y + 10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.) *
          (std::pow((10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) +
                        10. * (std::pow(x, 2.)) * (2. * x - 2.) * (x - 1. / 2.) +
                        20. * x * (std::pow((x - 1.), 2.)) * (x - 1. / 2.)),
              2.));
  p_ex = (std::pow(x, 2.)) + (std::pow(y, 2.));
  v_ex(0) =
      -5. * (std::pow(y, 4.)) *
          ((std::pow(
               (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 5.)) /
                  5. +
              (std::pow(y, 5.)) / 5.) -
      4. * (std::pow(y, 3.)) *
          ((std::pow(
               (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 6.)) /
                  6. -
              (std::pow(y, 6.)) / 6.);
  v_ex(1) =
      (std::pow(y, 4.)) *
      (std::pow((10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.), 4.)) *
      (y + 10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) * (x - 1. / 2.) - 1.) *
      (10. * (std::pow(x, 2.)) * (std::pow((x - 1.), 2.)) +
          10. * (std::pow(x, 2.)) * (2. * x - 2.) * (x - 1. / 2.) +
          20. * x * (std::pow((x - 1.), 2.)) * (x - 1. / 2.));

  // density - momentum formulation
  r_ex = r0 + epsilon * (p_ex - p0);
  w_ex(0) = r_ex * v_ex(0);
  w_ex(1) = r_ex * v_ex(1);

  switch (component)
  {
    case 0:
      return r_ex;
    case 1:
      return w_ex(0);
    case 2:
      return w_ex(1);
    case 3:
      return L_ex(0);
    case 4:
      return L_ex(1);
    case 5:
      return L_ex(2);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneCFDFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneCFDForceFunction::WeaklyCompressibleEtienneCFDForceFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams)
    : refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneCFDForceFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double r0 = refdensity_;

  // initialize variables
  double f_r_ex;
  CORE::LINALG::Matrix<2, 1> f_w_ex;

  // evaluate variables
  f_r_ex = 0;
  f_w_ex(0) =
      2. * x -
      (y *
          (20. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              12. * (std::pow(y, 2.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.) -
              (std::pow(y, 8.)) +
              (std::pow(y, 4.)) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.) +
              100. * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) -
              400. * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.))) /
          5. +
      ((std::pow(x, 2.)) / 10. + (std::pow(y, 2.)) / 10. + 1. / 10.) *
          (24. * y *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.) -
              60. * (std::pow(y, 2.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) +
              16. * (std::pow(y, 7.)) -
              (std::pow(y, 4.)) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) -
              4. * (std::pow(y, 3.)) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.) -
              400. * (std::pow(x, 2.)) * (std::pow(y, 3.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) +
              400. * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) +
              1600. * (std::pow(x, 2.)) * (std::pow(y, 3.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.)) -
      20. * (std::pow(y, 3.)) * ((std::pow(x, 2.)) / 10. + (std::pow(y, 2.)) / 10. + 1. / 10.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              3.)) *
          (48. * x + 5. * y - 60. * x * y + 375. * (std::pow(x, 2.)) * y -
              2900. * (std::pow(x, 3.)) * y + 13275. * (std::pow(x, 4.)) * y -
              31050. * (std::pow(x, 5.)) * y + 38350. * (std::pow(x, 6.)) * y -
              24000. * (std::pow(x, 7.)) * y + 6000. * (std::pow(x, 8.)) * y -
              360. * (std::pow(x, 2.)) + 3120. * (std::pow(x, 3.)) - 15620. * (std::pow(x, 4.)) +
              52080. * (std::pow(x, 5.)) - 166360. * (std::pow(x, 6.)) +
              504000. * (std::pow(x, 7.)) - 1141500. * (std::pow(x, 8.)) +
              1737200. * (std::pow(x, 9.)) - 1720400. * (std::pow(x, 10.)) +
              1066800. * (std::pow(x, 11.)) - 377000. * (std::pow(x, 12.)) +
              58000. * (std::pow(x, 13.)) - 4.) +
      4. * (std::pow(x, 2.)) * (std::pow(y, 3.)) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
              20. * (std::pow(x, 2.)) + 5. * y - 4.) +
      10. * r0 * x * (std::pow(y, 4.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) +
      40. * r0 * x * (std::pow(y, 3.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
              5. * (std::pow(x, 2.)) + y - 1.) -
      20. * r0 * x * (std::pow(y, 3.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
              20. * (std::pow(x, 2.)) + 5. * y - 4.) -
      10. * r0 * x * (std::pow(y, 4.)) *
          (12. * (std::pow(y, 2.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.) -
              20. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) +
              (std::pow(y, 8.))) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
              5. * (std::pow(x, 2.)) + y - 1.);
  f_w_ex(1) =
      2. * y -
      (x *
          (20. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              12. * (std::pow(y, 2.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.) -
              (std::pow(y, 8.)) +
              (std::pow(y, 4.)) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.) +
              100. * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) -
              400. * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.))) /
          5. -
      ((std::pow(x, 2.)) / 10. + (std::pow(y, 2.)) / 10. + 1. / 10.) *
          (120. * x * (std::pow(y, 2.)) *
                  (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      5.)) -
              8000. * (std::pow(x, 3.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 3.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) -
              200. * x * (std::pow(y, 3.)) *
                  (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) +
              (std::pow(y, 4.)) * (600. * (std::pow(x, 2.)) - 600. * x + 120.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.) +
              12000. * (std::pow(x, 3.)) * (std::pow(y, 4.)) *
                  (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 3.)) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      2.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.) +
              30. * x * (std::pow(y, 4.)) *
                  (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      4.)) -
              120. * x * (std::pow(y, 4.)) *
                  (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
                  (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
                  (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                      3.)) *
                  (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                      5. * (std::pow(x, 2.)) + y - 1.)) -
      4. * x * (std::pow(y, 4.)) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
              20. * (std::pow(x, 2.)) + 5. * y - 4.) +
      100. * r0 * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) +
      100. * r0 * (std::pow(x, 2.)) * (std::pow(y, 8.)) *
          (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              8.)) *
          (20. * (std::pow(x, 5.)) - 50. * (std::pow(x, 4.)) + 40. * (std::pow(x, 3.)) -
              10. * (std::pow(x, 2.)) + 2. * y - 2.) +
      800. * r0 * (std::pow(x, 2.)) * (std::pow(y, 7.)) *
          (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              8.)) *
          (std::pow((10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
                        5. * (std::pow(x, 2.)) + y - 1.),
              2.)) +
      r0 * (std::pow(y, 4.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (200. * (std::pow(x, 3.)) - 300. * (std::pow(x, 2.)) + 120. * x - 10.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
              5. * (std::pow(x, 2.)) + y - 1.) -
      80. * x * (std::pow(y, 2.)) * ((std::pow(x, 2.)) / 10. + (std::pow(y, 2.)) / 10. + 1. / 10.) *
          (5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              4.)) *
          (30. * (std::pow(x, 5.)) - 75. * (std::pow(x, 4.)) + 60. * (std::pow(x, 3.)) -
              15. * (std::pow(x, 2.)) + 5. * y - 3.) -
      400. * r0 * (std::pow(x, 2.)) * (std::pow(y, 4.)) *
          (5. * (std::pow(y, 4.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       5.)) /
                          5. -
                      (std::pow(y, 5.)) / 5.) -
              4. * (std::pow(y, 3.)) *
                  ((std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) -
                                 20. * (std::pow(x, 3.)) + 5. * (std::pow(x, 2.)) + 1.),
                       6.)) /
                          6. -
                      (std::pow(y, 6.)) / 6.)) *
          (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              3.)) *
          (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
              5. * (std::pow(x, 2.)) + y - 1.) -
      100. * r0 * (std::pow(x, 2.)) * (std::pow(y, 7.)) *
          (std::pow((5. * (std::pow(x, 3.)) - 10. * (std::pow(x, 2.)) + 6. * x - 1.), 2.)) *
          (std::pow((25. * (std::pow(x, 4.)) - 10. * (std::pow(x, 5.)) - 20. * (std::pow(x, 3.)) +
                        5. * (std::pow(x, 2.)) + 1.),
              8.)) *
          (40. * (std::pow(x, 5.)) - 100. * (std::pow(x, 4.)) + 80. * (std::pow(x, 3.)) -
              20. * (std::pow(x, 2.)) + 5. * y - 4.) *
          (10. * (std::pow(x, 5.)) - 25. * (std::pow(x, 4.)) + 20. * (std::pow(x, 3.)) -
              5. * (std::pow(x, 2.)) + y - 1.);

  switch (component)
  {
    case 0:
      return f_r_ex;
    case 1:
      return f_w_ex(0);
    case 2:
      return f_w_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneCFDForceFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneCFDViscosityFunction::WeaklyCompressibleEtienneCFDViscosityFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams)
    : refdensity_(fparams.refdensity_),
      refpressure_(fparams.refpressure_),
      comprcoeff_(fparams.comprcoeff_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneCFDViscosityFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];

  // initialize variables
  double mu_ex;

  // evaluate variables
  mu_ex = (1. + (std::pow(x, 2.)) + (std::pow(y, 2.))) / 10.;

  switch (component)
  {
    case 0:
      return mu_ex;

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneCFDViscosityFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneFSIFluidFunction::WeaklyCompressibleEtienneFSIFluidFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
    const MAT::PAR::StVenantKirchhoff& fparams_struc)
    : refdensity_(0.0),
      refpressure_(0.0),
      comprcoeff_(0.0),
      youngmodulus_(0.0),
      poissonratio_(0.0),
      strucdensity_(0.0)
{
  // get data
  refdensity_ = fparams_fluid.refdensity_;
  refpressure_ = fparams_fluid.refpressure_;
  comprcoeff_ = fparams_fluid.comprcoeff_;

  youngmodulus_ = fparams_struc.youngs_;
  poissonratio_ = fparams_struc.poissonratio_;
  strucdensity_ = fparams_struc.density_;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneFSIFluidFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double r0 = refdensity_;
  double p0 = refpressure_;
  double epsilon = comprcoeff_;
  double E = youngmodulus_;
  double v = poissonratio_;

  // initialize variables
  CORE::LINALG::Matrix<3, 1> L_ex;
  double p_ex;
  CORE::LINALG::Matrix<2, 1> v_ex;
  double r_ex;
  CORE::LINALG::Matrix<2, 1> w_ex;

  // evaluate variables
  L_ex(0) = (M_PI * sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                (cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                    cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) + 20.) *
                (40. * y - cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                    cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) - 20.)) /
                12000. -
            (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 10.) *
                (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                    20. * y + 10.)) /
                750.;
  L_ex(1) = (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 10.) *
                (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                    20. * y + 10.)) /
                1500. -
            (M_PI * sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                (cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                    cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) + 20.) *
                (40. * y - cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                    cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) - 20.)) /
                6000.;
  L_ex(2) =
      (2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
          (std::pow(
              (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.), 2.))) /
          25. +
      3. * (std::pow(y, 2.)) * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) -
      6. * (std::pow(y, 2.)) * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) +
      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) / 5. +
      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
          (cos(2. * M_PI * x) - 1.)) /
          5. -
      ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
          ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) / 20. -
              1.) *
          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
          (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
          100. -
      ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
          (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
          (y +
              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                  20. -
              1.)) /
          100. +
      ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
          ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) / 20. -
              1.) *
          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
          (y +
              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                  20. -
              1.)) /
          5.;
  p_ex =
      -((E * M_PI * cos(2. * M_PI * t) *
            (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
            (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
            ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                 (std::pow(
                     (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                     2.))) /
                    25. +
                (3. *
                    (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.),
                        2.)) *
                    ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                    400. -
                (3. *
                    (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.),
                        2.)) *
                    ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                    200. +
                ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                    5. +
                ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                    (cos(2. * M_PI * x) - 1.)) /
                    5. +
                ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                    (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.),
                        2.)) *
                    ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                    (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                        2.))) /
                    40000.)) /
              (6. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.)) +
          (E * M_PI * cos(2. * M_PI * t) *
              (std::pow(
                  (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.) -
                      20.),
                  2.)) *
              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
              (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                  6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                  3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
              (1200. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.))) /
      (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
           (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
           ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                (std::pow(
                    (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                    2.))) /
                   25. +
               (3. *
                   (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                     (cos(2. * M_PI * x) - 1.) -
                                 20.),
                       2.)) *
                   ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                   400. -
               (3. *
                   (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                     (cos(2. * M_PI * x) - 1.) -
                                 20.),
                       2.)) *
                   ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                   200. +
               ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                   5. +
               ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                   (cos(2. * M_PI * x) - 1.)) /
                   5. +
               ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                     (cos(2. * M_PI * x) - 1.) -
                                 20.),
                       2.)) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                       2.))) /
                   40000.)) /
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) -
          (100. *
              ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                   (std::pow(
                       (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                       2.))) /
                      25. +
                  (3. *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                      400. -
                  (3. *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                      200. +
                  ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                      (std::pow((sin(2. * M_PI * x)), 2.))) /
                      5. +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                      (cos(2. * M_PI * x) - 1.)) /
                      5. +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      40000.)) /
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) +
          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
              (std::pow(
                  (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.) -
                      20.),
                  2.)) *
              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
              (std::pow(
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
              (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.)) +
          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
              sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
              (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 10.) *
              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
              (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                  2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                  10.)) /
              (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                             (std::pow((cos(2. * M_PI * x) -
                                           2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                 2.)) +
                         100.)));
  v_ex(0) =
      (((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) *
          (std::pow(
              (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.), 3.))) /
          125. -
      2. * y *
          ((((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
               (std::pow(
                   (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                   2.))) /
                  25. -
              (std::pow(y, 2.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) -
      (std::pow(y, 3.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.);
  v_ex(1) =
      -(M_PI * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
          10. -
      y * (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
          ((M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) * sin(2. * M_PI * (t + 1. / 4.))) / 10. -
              (M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  10.) *
          ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) / 20. -
              1.) *
          (y +
              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                  20. -
              1.);

  // density - momentum formulation
  r_ex = r0 + epsilon * (p_ex - p0);
  w_ex(0) = r_ex * v_ex(0);
  w_ex(1) = r_ex * v_ex(1);

  switch (component)
  {
    case 0:
      return r_ex;
    case 1:
      return w_ex(0);
    case 2:
      return w_ex(1);
    case 3:
      return L_ex(0);
    case 4:
      return L_ex(1);
    case 5:
      return L_ex(2);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneFSIFluidFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneFSIFluidForceFunction::WeaklyCompressibleEtienneFSIFluidForceFunction(
    const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
    const MAT::PAR::StVenantKirchhoff& fparams_struc)
    : refdensity_(0.0),
      refpressure_(0.0),
      comprcoeff_(0.0),
      youngmodulus_(0.0),
      poissonratio_(0.0),
      strucdensity_(0.0)
{
  // get data
  refdensity_ = fparams_fluid.refdensity_;
  refpressure_ = fparams_fluid.refpressure_;
  comprcoeff_ = fparams_fluid.comprcoeff_;

  youngmodulus_ = fparams_struc.youngs_;
  poissonratio_ = fparams_struc.poissonratio_;
  strucdensity_ = fparams_struc.density_;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneFSIFluidForceFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double r0 = refdensity_;
  double E = youngmodulus_;
  double v = poissonratio_;

  // initialize variables
  double f_r_ex;
  CORE::LINALG::Matrix<2, 1> f_w_ex;

  // evaluate variables
  f_r_ex = 0.;
  f_w_ex(0) =
      (r0 * (2. * y *
                    ((14. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * t) *
                         (std::pow(
                             (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                 5.),
                             2.))) /
                            25. -
                        14. * M_PI * (std::pow(y, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * t) +
                        (4. * M_PI * cos(M_PI * x) * sin(2. * M_PI * t) *
                            (std::pow((sin(M_PI * x)), 3.)) *
                            ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.)) /
                            25.) +
                (28. * M_PI * (std::pow(y, 3.)) * cos(2. * M_PI * t) * sin(2. * M_PI * t)) / 3. -
                (28. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * t) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        3.))) /
                    375. -
                (6. * M_PI * cos(M_PI * x) * sin(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 3.)) *
                    ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        2.))) /
                    125.) +
          (r0 * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
              ((M_PI * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                   (cos(2. * M_PI * x) - 1.)) /
                      10. +
                  (M_PI * y * cos(2. * M_PI * t) *
                      ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                              20. -
                          1.) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (y +
                          (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                              (cos(2. * M_PI * x) - 1.)) /
                              20. -
                          1.)) /
                      10.) *
              ((std::pow((cos(2. * M_PI * t)), 2.)) * (std::pow((cos(M_PI * x)), 2.)) *
                      (std::pow((sin(M_PI * x)), 6.)) -
                  50. * (std::pow(y, 2.)) +
                  10. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                  25.)) /
              25. -
          (M_PI * r0 * cos(2. * M_PI * t) * (7. * cos(4. * M_PI * t) - 9.) *
              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
              (2. * y *
                      ((((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. -
                          (std::pow(y, 2.)) *
                              ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) -
                  (((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) *
                      (std::pow(
                          (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                              5.),
                          3.))) /
                      125. +
                  (std::pow(y, 3.)) *
                      ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.)) *
              (cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                  cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) + 20.) *
              (40. * y - cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                  cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) - 20.)) /
              8000. -
          (2. * M_PI * r0 * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 2.)) *
              (7. * (std::pow((cos(2. * M_PI * t)), 2.)) - 8.) *
              (4. * (std::pow((cos(M_PI * x)), 2.)) - 1.) *
              (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.) *
              (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) - 10. * y +
                  5.) *
              (2. * y *
                      ((((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. -
                          (std::pow(y, 2.)) *
                              ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) -
                  (((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) *
                      (std::pow(
                          (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                              5.),
                          3.))) /
                      125. +
                  (std::pow(y, 3.)) *
                      ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.))) /
              125.) +
      ((5. * E *
           (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
               6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
               3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
           (12. * y * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) -
               6. * y * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) +
               ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                       2.))) /
                   100. +
               ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                       2.))) /
                   100. +
               ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   (y +
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.)) /
                   100. -
               ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                   ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.)) /
                   5. -
               ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                   ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                   (y +
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.)) /
                   5.)) /
              (3. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) *
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.))) /
                               25. +
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                               400. -
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                               cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) *
                          sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)))) +
          (5. * E *
              ((M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                   (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                   (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                       (std::pow((sin(2. * M_PI * x)), 2.))) *
                   (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       20.) *
                   (40. * y - sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       20.)) /
                      6000. -
                  (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      375.) *
              (10. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                  6. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                  24. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * cos(2. * M_PI * x) *
                      sin(2. * M_PI * x))) /
              (3. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) *
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.))) /
                               25. +
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                               400. -
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                               cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) *
                          sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)))) +
          (5. * E *
              (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                  6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                  3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
              ((M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                   (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                   (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                       (std::pow((sin(2. * M_PI * x)), 2.))) *
                   (2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       2. * M_PI * (std::pow((cos(2. * M_PI * x)), 2.)) *
                           sin(2. * M_PI * (t + 1. / 4.)) +
                       2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                           sin(2. * M_PI * (t + 1. / 4.))) *
                   (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       20.)) /
                      6000. -
                  (M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                      (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                      (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                          (std::pow((sin(2. * M_PI * x)), 2.))) *
                      (2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                          2. * M_PI * (std::pow((cos(2. * M_PI * x)), 2.)) *
                              sin(2. * M_PI * (t + 1. / 4.)) +
                          2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                              sin(2. * M_PI * (t + 1. / 4.))) *
                      (40. * y - sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                          cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                          20.)) /
                      6000. +
                  (M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                      (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                          cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                          20.) *
                      (40. * y - sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                          cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                          20.)) /
                      6000. -
                  (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.))) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      375. -
                  (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.))) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.)) /
                      375. +
                  ((std::pow(M_PI, 2.)) * cos(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      375. +
                  ((std::pow(M_PI, 2.)) * cos(3. * M_PI * x) * sin(M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      125.)) /
              (3. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) *
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.))) /
                               25. +
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                               400. -
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                               cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) *
                          sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)))) +
          (5. * E *
              ((M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                   (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                   (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                       (std::pow((sin(2. * M_PI * x)), 2.))) *
                   (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       20.) *
                   (40. * y - sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       20.)) /
                      6000. -
                  (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      375.) *
              (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                  6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                  3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
              ((100. *
                   ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                            3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                                (std::pow((sin(M_PI * x)), 2.))) *
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                            5.)) /
                           25. +
                       (3. *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           200. -
                       (3. *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           100. +
                       (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                           cos(2. * M_PI * (t + 1. / 4.))) /
                           5. +
                       (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                           cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                           5. -
                       (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                           sin(2. * M_PI * x)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           20000. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (2. * M_PI * sin(2. * M_PI * x) -
                               8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                           20000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) *
                      ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                               3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                                   (std::pow((sin(M_PI * x)), 2.))) *
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.)) /
                              25. +
                          (3. *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              200. -
                          (3. *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              100. +
                          (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.))) /
                              5. +
                          (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. -
                          (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                              sin(2. * M_PI * x)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              20000. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (2. * M_PI * sin(2. * M_PI * x) -
                                  8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) +
                                  1.)) /
                              20000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (2. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (200. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      (std::pow(
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.),
                          2.)) +
                  (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          3.)) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      (std::pow(
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.),
                          2.)) +
                  ((std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          3.))) /
                      (100. *
                          (std::pow(
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.),
                              2.))) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                              sin(2. * M_PI * (t + 1. / 4.)) -
                          2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                              (cos(2. * M_PI * x) - 1.)) *
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                              (cos(2. * M_PI * x) - 1.) -
                          20.) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                      (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) -
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.))) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  ((std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  (3. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(3. * M_PI * x) *
                      sin(M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) -
                          2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 4.)) -
                          2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 2.)) +
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.)))) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) -
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 3.)) *
                      sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. *
                          (std::pow(
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.),
                              2.))))) /
              (3. * (v + 1.) *
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) *
                  (std::pow(
                      (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) *
                           ((2. *
                                (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                  (std::pow((sin(M_PI * x)), 3.)) +
                                              5.),
                                    2.)) *
                                ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                                   25. +
                               (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                                   (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                     (cos(2. * M_PI * x) - 1.) -
                                                 20.),
                                       2.))) /
                                   400. -
                               (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                                   (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                     (cos(2. * M_PI * x) - 1.) -
                                                 20.),
                                       2.))) /
                                   200. +
                               ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                   (std::pow((sin(2. * M_PI * x)), 2.))) /
                                   5. +
                               ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                                   cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                                   5. +
                               ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                   (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                     (cos(2. * M_PI * x) - 1.) -
                                                 20.),
                                       2.)) *
                                   (std::pow((cos(2. * M_PI * x) -
                                                 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                       2.))) /
                                   40000.)) /
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.) -
                          (100. *
                              ((2. *
                                   (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                     (std::pow((sin(M_PI * x)), 3.)) +
                                                 5.),
                                       2.)) *
                                   ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                                      25. +
                                  (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                                      (std::pow(
                                          (sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                  (cos(2. * M_PI * x) - 1.) -
                                              20.),
                                          2.))) /
                                      400. -
                                  (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                                      (std::pow(
                                          (sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                  (cos(2. * M_PI * x) - 1.) -
                                              20.),
                                          2.))) /
                                      200. +
                                  ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                      (std::pow((sin(2. * M_PI * x)), 2.))) /
                                      5. +
                                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                                      cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                                      5. +
                                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                      (std::pow(
                                          (sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                  (cos(2. * M_PI * x) - 1.) -
                                              20.),
                                          2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.))) /
                                      40000.)) /
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.) +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              (200. *
                                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                          (std::pow(
                                              (cos(2. * M_PI * x) -
                                                  2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                              2.)) +
                                      100.)) +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * (t + 1. / 4.)) *
                              cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                              (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                      (std::pow((sin(M_PI * x)), 3.)) +
                                  10.) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) +
                                  1.) *
                              (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                                  2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                      sin(M_PI * x) +
                                  10.)) /
                              (25. *
                                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                          (std::pow(
                                              (cos(2. * M_PI * x) -
                                                  2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                              2.)) +
                                      100.))),
                      2.))) -
          (10. * E * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
              ((M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
                   (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
                   (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                       (std::pow((sin(2. * M_PI * x)), 2.))) *
                   (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       20.) *
                   (40. * y - sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                       cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                       20.)) /
                      6000. -
                  (M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) -
                          20. * y + 10.)) /
                      375.) *
              (2. * M_PI * sin(2. * M_PI * x) -
                  8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
              (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                  6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                  3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
              (3. * (v + 1.) *
                  (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                    (std::pow((cos(2. * M_PI * x) -
                                                  2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                        2.)) +
                                100.),
                      2.)) *
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.))) /
                               25. +
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                               400. -
                           (3. *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                               cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) *
                          sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.))))) +
      (((E * M_PI * cos(2. * M_PI * t) *
            (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
            (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
            ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                 (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                     3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                         (std::pow((sin(M_PI * x)), 2.))) *
                 (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.)) /
                    25. +
                (3. *
                    (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                            sin(2. * M_PI * (t + 1. / 4.)) -
                        2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.)) *
                    (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.) -
                        20.) *
                    ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                    200. -
                (3. *
                    (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                            sin(2. * M_PI * (t + 1. / 4.)) -
                        2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.)) *
                    (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.) -
                        20.) *
                    ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                    100. +
                (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                    cos(2. * M_PI * (t + 1. / 4.))) /
                    5. +
                (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                    (cos(2. * M_PI * x) - 1.)) /
                    5. -
                (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                    sin(2. * M_PI * x)) /
                    5. +
                ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                    (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                            sin(2. * M_PI * (t + 1. / 4.)) -
                        2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.)) *
                    (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                            (cos(2. * M_PI * x) - 1.) -
                        20.) *
                    ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                    (std::pow((cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                        2.))) /
                    20000. +
                ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                    (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.),
                        2.)) *
                    ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                    (2. * M_PI * sin(2. * M_PI * x) -
                        8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                    (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                    20000.)) /
               (6. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.)) +
           (E * M_PI * cos(2. * M_PI * t) *
               (10. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                   6. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                   24. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * cos(2. * M_PI * x) *
                       sin(2. * M_PI * x)) *
               (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
               ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        2.))) /
                       25. +
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                       400. -
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                       200. +
                   ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                       (std::pow((sin(2. * M_PI * x)), 2.))) /
                       5. +
                   ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                       5. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.))) /
                       40000.)) /
               (6. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.)) +
           (E * M_PI * cos(2. * M_PI * t) *
               (2. * M_PI * sin(2. * M_PI * x) -
                   8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
               (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                   6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                   3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
               ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        2.))) /
                       25. +
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                       400. -
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                       200. +
                   ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                       (std::pow((sin(2. * M_PI * x)), 2.))) /
                       5. +
                   ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                       5. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.))) /
                       40000.)) /
               (6. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.)) +
           (E * M_PI * cos(2. * M_PI * t) *
               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                 (cos(2. * M_PI * x) - 1.) -
                             20.),
                   2.)) *
               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
               (10. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * x) -
                   6. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * x) +
                   24. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * cos(2. * M_PI * x) *
                       sin(2. * M_PI * x)) *
               (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
               (1200. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.)) -
           (E * (std::pow(M_PI, 3.)) * (std::pow((cos(2. * M_PI * t)), 3.)) *
               (2. * M_PI * sin(2. * M_PI * x) -
                   8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
               (std::pow(
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
               (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                   6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                   3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
               ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        2.))) /
                       25. +
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                       400. -
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                       200. +
                   ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                       (std::pow((sin(2. * M_PI * x)), 2.))) /
                       5. +
                   ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                       5. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.))) /
                       40000.)) /
               (3. * (v + 1.) *
                   (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.),
                       2.))) +
           (E * M_PI * cos(2. * M_PI * t) *
               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                 (cos(2. * M_PI * x) - 1.) -
                             20.),
                   2.)) *
               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
               (2. * M_PI * sin(2. * M_PI * x) -
                   8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
               (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                   6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                   3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
               (1200. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.)) -
           (E * (std::pow(M_PI, 3.)) * (std::pow((cos(2. * M_PI * t)), 3.)) *
               (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                 (cos(2. * M_PI * x) - 1.) -
                             20.),
                   2.)) *
               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
               (2. * M_PI * sin(2. * M_PI * x) -
                   8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
               (std::pow(
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
               (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                   6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                   3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
               (600. * (v + 1.) *
                   (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.),
                       2.))) +
           (E * M_PI * cos(2. * M_PI * t) *
               (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) * sin(2. * M_PI * (t + 1. / 4.)) -
                   2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) *
               (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.) -
                   20.) *
               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
               (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
               (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                   6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                   3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
               (600. * (v + 1.) *
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.)) +
                       100.))) /
              (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (std::pow(
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.),
                            2.))) /
                           25. +
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           400. -
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           200. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                           (std::pow((sin(2. * M_PI * x)), 2.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.))) +
          (((E * M_PI * cos(2. * M_PI * t) *
                (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                    6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                    3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
                ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                     (std::pow(
                         (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                             5.),
                         2.))) /
                        25. +
                    (3. *
                        (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                          (cos(2. * M_PI * x) - 1.) -
                                      20.),
                            2.)) *
                        ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                        400. -
                    (3. *
                        (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                          (cos(2. * M_PI * x) - 1.) -
                                      20.),
                            2.)) *
                        ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                        200. +
                    ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                        (std::pow((sin(2. * M_PI * x)), 2.))) /
                        5. +
                    ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                        5. +
                    ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                        (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                          (cos(2. * M_PI * x) - 1.) -
                                      20.),
                            2.)) *
                        ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                        (std::pow(
                            (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                            2.))) /
                        40000.)) /
                   (6. * (v + 1.) *
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.)) +
                           100.)) +
               (E * M_PI * cos(2. * M_PI * t) *
                   (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                     (cos(2. * M_PI * x) - 1.) -
                                 20.),
                       2.)) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                   (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
                       6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
                       3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
                   (1200. * (v + 1.) *
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.)) +
                           100.))) *
              ((100. *
                   ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                            3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                                (std::pow((sin(M_PI * x)), 2.))) *
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                            5.)) /
                           25. +
                       (3. *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           200. -
                       (3. *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           100. +
                       (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                           cos(2. * M_PI * (t + 1. / 4.))) /
                           5. +
                       (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                           cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                           5. -
                       (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                           sin(2. * M_PI * x)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                   sin(2. * M_PI * (t + 1. / 4.)) -
                               2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.)) *
                           (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                   (cos(2. * M_PI * x) - 1.) -
                               20.) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           20000. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (2. * M_PI * sin(2. * M_PI * x) -
                               8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                           20000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) *
                      ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                               3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                                   (std::pow((sin(M_PI * x)), 2.))) *
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.)) /
                              25. +
                          (3. *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              200. -
                          (3. *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              100. +
                          (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.))) /
                              5. +
                          (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. -
                          (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                              sin(2. * M_PI * x)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                      sin(2. * M_PI * (t + 1. / 4.)) -
                                  2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.)) *
                              (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                      (cos(2. * M_PI * x) - 1.) -
                                  20.) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              20000. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (2. * M_PI * sin(2. * M_PI * x) -
                                  8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) +
                                  1.)) /
                              20000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (2. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (200. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      (std::pow(
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.),
                          2.)) +
                  (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          3.)) *
                      ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                      (std::pow(
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.),
                          2.)) +
                  ((std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          3.))) /
                      (100. *
                          (std::pow(
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.),
                              2.))) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                              sin(2. * M_PI * (t + 1. / 4.)) -
                          2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                              (cos(2. * M_PI * x) - 1.)) *
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                              (cos(2. * M_PI * x) - 1.) -
                          20.) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) -
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                      (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) -
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.))) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  ((std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  (3. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(3. * M_PI * x) *
                      sin(M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) -
                          2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 4.)) -
                          2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 2.)) +
                          6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                              (std::pow((sin(M_PI * x)), 2.)))) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) -
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)) +
                  (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 3.)) *
                      sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * M_PI * sin(2. * M_PI * x) -
                          8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. *
                          (std::pow(
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.),
                              2.))))) /
              (std::pow(
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.)) *
                            ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                               25. +
                           (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.))) /
                               400. -
                           (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.))) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                               cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.)) *
                               ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                                  25. +
                              (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.))) /
                                  400. -
                              (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.))) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                                  cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * (t + 1. / 4.)) * cos(2. * M_PI * t) *
                          sin(M_PI * x) * sin(3. * M_PI * x) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.))),
                  2.)));
  f_w_ex(1) =
      (5. * E *
          ((M_PI * sin(2. * M_PI * (t + 1. / 4.)) *
               (7. * (std::pow((sin(2. * M_PI * (t + 1. / 4.))), 2.)) - 8.) *
               (cos(2. * M_PI * x) - (std::pow((cos(2. * M_PI * x)), 2.)) +
                   (std::pow((sin(2. * M_PI * x)), 2.))) *
               (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) -
                   cos(2. * M_PI * x) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) +
                   20.)) /
                  75. +
              (2. * M_PI * sin(M_PI * x) * sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.)) /
                  75.) *
          (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
              6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
              3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
          (3. * (v + 1.) *
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) *
              (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (std::pow(
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.),
                            2.))) /
                           25. +
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           400. -
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           200. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                           (std::pow((sin(2. * M_PI * x)), 2.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)))) -
      r0 *
          (2. * y *
                  ((((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                       (std::pow(
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.),
                           2.))) /
                          25. -
                      (std::pow(y, 2.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) -
              (((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) *
                  (std::pow(
                      (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                      3.))) /
                  125. +
              (std::pow(y, 3.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.)) *
          (((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
               ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                       20. -
                   1.) *
               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
               (std::pow(
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  100. -
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                  5. +
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  100. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (sin(2. * M_PI * x) - 2. * sin(4. * M_PI * x)) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5.) -
      r0 * ((14. * (std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * cos(2. * M_PI * (t + 1. / 4.)) *
                sin(2. * M_PI * (t + 1. / 4.)) *
                ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                        20. -
                    1.) *
                (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                (y +
                    (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                        20. -
                    1.)) /
                   5. -
               ((std::pow(M_PI, 2.)) * y * sin(2. * M_PI * t) *
                   ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                   (y +
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.)) /
                   5. -
               ((std::pow(M_PI, 2.)) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                   (cos(2. * M_PI * x) - 1.)) /
                   5. +
               ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                   cos(2. * M_PI * (t + 1. / 4.)) *
                   ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                        (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.) *
                   ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (cos(2. * M_PI * x) - 1.) *
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                   100. +
               ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                   cos(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                   (cos(2. * M_PI * x) - 1.) *
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                   (y +
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           20. -
                       1.)) /
                   100.) +
      (5. * E *
          (10. * M_PI * cos(2. * M_PI * t) * sin(2. * M_PI * x) -
              6. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(2. * M_PI * x) +
              24. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * cos(2. * M_PI * x) *
                  sin(2. * M_PI * x)) *
          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
               (std::pow(
                   (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                   2.))) /
                  25. -
              6. * (std::pow(y, 2.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) +
              3. * (std::pow(y, 2.)) *
                  ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) +
              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                  5. +
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  100. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  100. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5.)) /
          (3. * (v + 1.) *
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) *
              (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (std::pow(
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.),
                            2.))) /
                           25. +
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           400. -
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           200. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                           (std::pow((sin(2. * M_PI * x)), 2.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)))) +
      (5. * E *
          (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
              6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
              3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
          ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
               (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                   3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                       (std::pow((sin(M_PI * x)), 2.))) *
               (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.)) /
                  25. +
              (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                  cos(2. * M_PI * (t + 1. / 4.))) /
                  5. +
              (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                  sin(2. * M_PI * x)) /
                  5. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) * sin(2. * M_PI * (t + 1. / 4.))) /
                          10. -
                      (M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          10.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  50. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  50. +
              (8. * (std::pow(M_PI, 3.)) * y * cos(2. * M_PI * t) *
                  (std::pow((sin(2. * M_PI * x)), 2.)) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                  50. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                  ((M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) * sin(2. * M_PI * (t + 1. / 4.))) /
                          10. -
                      (M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          10.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5. -
              (2. * (std::pow(M_PI, 3.)) * y * cos(2. * M_PI * t) * cos(2. * M_PI * x) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                  ((M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) * sin(2. * M_PI * (t + 1. / 4.))) /
                          10. -
                      (M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          10.) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.)) /
                  5.)) /
          (3. * (v + 1.) *
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) *
              (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (std::pow(
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.),
                            2.))) /
                           25. +
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           400. -
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           200. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                           (std::pow((sin(2. * M_PI * x)), 2.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.)))) +
      (5. * E *
          (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
              6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
              3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
               (std::pow(
                   (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                   2.))) /
                  25. -
              6. * (std::pow(y, 2.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) +
              3. * (std::pow(y, 2.)) *
                  ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) +
              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                  5. +
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  100. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  100. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5.) *
          ((100. *
               ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                    (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                        3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                            (std::pow((sin(M_PI * x)), 2.))) *
                    (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.)) /
                       25. +
                   (3. *
                       (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                               sin(2. * M_PI * (t + 1. / 4.)) -
                           2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.)) *
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.) -
                           20.) *
                       ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                       200. -
                   (3. *
                       (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                               sin(2. * M_PI * (t + 1. / 4.)) -
                           2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.)) *
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.) -
                           20.) *
                       ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                       100. +
                   (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                       cos(2. * M_PI * (t + 1. / 4.))) /
                       5. +
                   (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                       cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                       5. -
                   (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                       sin(2. * M_PI * x)) /
                       5. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                               sin(2. * M_PI * (t + 1. / 4.)) -
                           2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.)) *
                       (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                               (cos(2. * M_PI * x) - 1.) -
                           20.) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.))) /
                       20000. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (2. * M_PI * sin(2. * M_PI * x) -
                           8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                       20000.)) /
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) -
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  ((4. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                       (M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                           3. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                               (std::pow((sin(M_PI * x)), 2.))) *
                       (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                           5.)) /
                          25. +
                      (3. *
                          (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                  sin(2. * M_PI * (t + 1. / 4.)) -
                              2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.)) *
                          (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.) -
                              20.) *
                          ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                          200. -
                      (3. *
                          (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                  sin(2. * M_PI * (t + 1. / 4.)) -
                              2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.)) *
                          (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.) -
                              20.) *
                          ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                          100. +
                      (2. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                          cos(2. * M_PI * (t + 1. / 4.))) /
                          5. +
                      (2. * (std::pow(M_PI, 3.)) * sin(2. * M_PI * x) *
                          cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                          5. -
                      (4. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * x) * sin(2. * M_PI * t) *
                          sin(2. * M_PI * x)) /
                          5. +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                                  sin(2. * M_PI * (t + 1. / 4.)) -
                              2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.)) *
                          (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                  (cos(2. * M_PI * x) - 1.) -
                              20.) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          20000. +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (2. * M_PI * sin(2. * M_PI * x) -
                              8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                          20000.)) /
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) -
              (2. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                       (std::pow(
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.),
                           2.))) /
                          25. +
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                          400. -
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                          200. +
                      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                          (std::pow((sin(2. * M_PI * x)), 2.))) /
                          5. +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          5. +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          40000.)) /
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) -
              (200. * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                       (std::pow(
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.),
                           2.))) /
                          25. +
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                          400. -
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                          200. +
                      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                          (std::pow((sin(2. * M_PI * x)), 2.))) /
                          5. +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          5. +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          40000.)) /
                  (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                    (std::pow((cos(2. * M_PI * x) -
                                                  2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                        2.)) +
                                100.),
                      2.)) +
              (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 3.)) *
                  ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                       (std::pow(
                           (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                               5.),
                           2.))) /
                          25. +
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                          400. -
                      (3. *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                          200. +
                      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                          (std::pow((sin(2. * M_PI * x)), 2.))) /
                          5. +
                      ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          5. +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          40000.)) /
                  (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                    (std::pow((cos(2. * M_PI * x) -
                                                  2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                        2.)) +
                                100.),
                      2.)) +
              ((std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 4.)) *
                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                    (cos(2. * M_PI * x) - 1.) -
                                20.),
                      2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 3.))) /
                  (100. * (std::pow(
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.),
                              2.))) -
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (2. * M_PI * (std::pow((sin(2. * M_PI * x)), 2.)) *
                          sin(2. * M_PI * (t + 1. / 4.)) -
                      2. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) *
                  (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.) -
                      20.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.)) -
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                    (cos(2. * M_PI * x) - 1.) -
                                20.),
                      2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.)) /
                  (100. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.)) -
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 4.)) -
                      6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                          (std::pow((sin(M_PI * x)), 2.))) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.)) +
              ((std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(M_PI * x) * sin(3. * M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.)) +
              (3. * (std::pow(M_PI, 3.)) * cos(2. * M_PI * t) * cos(3. * M_PI * x) * sin(M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.)) +
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) -
                      2. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 4.)) -
                      2. * M_PI * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 2.)) +
                      6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 2.)) *
                          (std::pow((sin(M_PI * x)), 2.)))) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.)) -
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.)) +
              (2. * (std::pow(M_PI, 4.)) * (std::pow((cos(2. * M_PI * t)), 3.)) * sin(M_PI * x) *
                  sin(3. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * M_PI * sin(2. * M_PI * x) -
                      8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * (std::pow(
                             ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.),
                             2.))))) /
          (3. * (v + 1.) *
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.)) +
                  100.) *
              (std::pow(
                  (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.)) *
                       ((2. *
                            (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                              (std::pow((sin(M_PI * x)), 3.)) +
                                          5.),
                                2.)) *
                            ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                               25. +
                           (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.))) /
                               400. -
                           (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.))) /
                               200. +
                           ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                               (std::pow((sin(2. * M_PI * x)), 2.))) /
                               5. +
                           ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                               cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                               5. +
                           ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                               ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                               (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                 (cos(2. * M_PI * x) - 1.) -
                                             20.),
                                   2.)) *
                               (std::pow((cos(2. * M_PI * x) -
                                             2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                   2.))) /
                               40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) -
                      (100. *
                          ((2. *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.)) *
                               ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.)) /
                                  25. +
                              (3. * ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.))) /
                                  400. -
                              (3. * ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.))) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * (t + 1. / 4.)) *
                                  cos(2. * M_PI * x) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.) +
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (std::pow((sin(2. * M_PI * (t + 1. / 4.)) * sin(2. * M_PI * x) *
                                            (cos(2. * M_PI * x) - 1.) -
                                        20.),
                              2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.))) /
                          (200. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                      ((std::pow(M_PI, 2.)) * sin(2. * M_PI * (t + 1. / 4.)) * cos(2. * M_PI * t) *
                          sin(M_PI * x) * sin(3. * M_PI * x) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) *
                                  (std::pow((sin(M_PI * x)), 3.)) +
                              10.) *
                          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                          (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                              2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                                  sin(M_PI * x) +
                              10.)) /
                          (25. *
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.))),
                  2.))) -
      (M_PI * r0 * sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
          ((M_PI * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
               (cos(2. * M_PI * x) - 1.)) /
                  10. +
              (M_PI * y * cos(2. * M_PI * t) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  10.) *
          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
          (cos(2. * M_PI * t) * sin(2. * M_PI * x) -
              cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) + 20.) *
          (40. * y - cos(2. * M_PI * t) * sin(2. * M_PI * x) +
              cos(2. * M_PI * t) * cos(2. * M_PI * x) * sin(2. * M_PI * x) - 20.)) /
          2000. -
      (M_PI * r0 * cos(2. * M_PI * t) * (std::pow((sin(M_PI * x)), 2.)) *
          ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (2. * cos(2. * M_PI * x) + 1.) *
          ((M_PI * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
               (cos(2. * M_PI * x) - 1.)) /
                  10. +
              (M_PI * y * cos(2. * M_PI * t) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  10.) *
          (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.) *
          (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) - 10. * y + 5.)) /
          125. -
      (10. * E * (std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
          (2. * M_PI * sin(2. * M_PI * x) - 8. * M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * x)) *
          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
          (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
              6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
              3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.) *
          ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
               (std::pow(
                   (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                   2.))) /
                  25. -
              6. * (std::pow(y, 2.)) * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) +
              3. * (std::pow(y, 2.)) *
                  ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 3. - 8. / 3.) +
              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) * (std::pow((sin(2. * M_PI * x)), 2.))) /
                  5. +
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  100. -
              ((std::pow(M_PI, 2.)) * y * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  100. +
              ((std::pow(M_PI, 2.)) * y * cos(2. * M_PI * t) * sin(2. * M_PI * x) *
                  ((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) * (4. * cos(2. * M_PI * x) - 1.) *
                  (y +
                      (sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                          (cos(2. * M_PI * x) - 1.)) /
                          20. -
                      1.)) /
                  5.)) /
          (3. * (v + 1.) *
              (std::pow(((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                (std::pow((cos(2. * M_PI * x) -
                                              2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                    2.)) +
                            100.),
                  2.)) *
              (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                   (std::pow(
                       (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
                   ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                        (std::pow(
                            (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                                5.),
                            2.))) /
                           25. +
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                           400. -
                       (3. *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                           200. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                           (std::pow((sin(2. * M_PI * x)), 2.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5. +
                       ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                           (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                             (cos(2. * M_PI * x) - 1.) -
                                         20.),
                               2.)) *
                           ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                           (std::pow((cos(2. * M_PI * x) -
                                         2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                               2.))) /
                           40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) -
                  (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                               (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                                 (std::pow((sin(M_PI * x)), 3.)) +
                                             5.),
                                   2.))) /
                                  25. +
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                                  400. -
                              (3. *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                                  200. +
                              ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                                  (std::pow((sin(2. * M_PI * x)), 2.))) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                                  cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                                  5. +
                              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                    (cos(2. * M_PI * x) - 1.) -
                                                20.),
                                      2.)) *
                                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.))) /
                                  40000.)) /
                      ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.)) +
                          100.) +
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                      (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                        (cos(2. * M_PI * x) - 1.) -
                                    20.),
                          2.)) *
                      ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (std::pow(
                          (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                          2.))) /
                      (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                      (std::pow((cos(2. * M_PI * x) -
                                                    2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                          2.)) +
                                  100.)) +
                  ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                      sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                          10.) *
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                      (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                          2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) *
                              sin(M_PI * x) +
                          10.)) /
                      (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                     (std::pow((cos(2. * M_PI * x) -
                                                   2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                         2.)) +
                                 100.))));

  switch (component)
  {
    case 0:
      return f_r_ex;
    case 1:
      return f_w_ex(0);
    case 2:
      return f_w_ex(1);

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneFSIFluidForceFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::WeaklyCompressibleEtienneFSIFluidViscosityFunction::
    WeaklyCompressibleEtienneFSIFluidViscosityFunction(
        const MAT::PAR::WeaklyCompressibleFluid& fparams_fluid,
        const MAT::PAR::StVenantKirchhoff& fparams_struc)
    : refdensity_(0.0),
      refpressure_(0.0),
      comprcoeff_(0.0),
      youngmodulus_(0.0),
      poissonratio_(0.0),
      strucdensity_(0.0)
{
  // get data
  refdensity_ = fparams_fluid.refdensity_;
  refpressure_ = fparams_fluid.refpressure_;
  comprcoeff_ = fparams_fluid.comprcoeff_;

  youngmodulus_ = fparams_struc.youngs_;
  poissonratio_ = fparams_struc.poissonratio_;
  strucdensity_ = fparams_struc.density_;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::WeaklyCompressibleEtienneFSIFluidViscosityFunction::Evaluate(
    const double* xp, const double t, const std::size_t component)
{
  // ease notation
  double x = xp[0];
  double E = youngmodulus_;
  double v = poissonratio_;

  // initialize variables
  double mu_ex;

  // evaluate variables
  mu_ex =
      -(5. * E *
          (5. * cos(2. * M_PI * t) * cos(2. * M_PI * x) - 3. * M_PI * cos(2. * M_PI * t) +
              6. * M_PI * cos(2. * M_PI * t) * (std::pow((cos(2. * M_PI * x)), 2.)) -
              3. * M_PI * cos(2. * M_PI * t) * cos(2. * M_PI * x) + 30.)) /
      (3. * (v + 1.) *
          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) +
              100.) *
          (((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
               (std::pow(
                   (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.)) *
               ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                    (std::pow(
                        (cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) + 5.),
                        2.))) /
                       25. +
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                       400. -
                   (3. *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                       200. +
                   ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                       (std::pow((sin(2. * M_PI * x)), 2.))) /
                       5. +
                   ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) *
                       (cos(2. * M_PI * x) - 1.)) /
                       5. +
                   ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                       (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                         (cos(2. * M_PI * x) - 1.) -
                                     20.),
                           2.)) *
                       ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                       (std::pow(
                           (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                           2.))) /
                       40000.)) /
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) -
              (100. * ((2. * ((7. * (std::pow((cos(2. * M_PI * t)), 2.))) / 2. - 4.) *
                           (std::pow((cos(2. * M_PI * t) * cos(M_PI * x) *
                                             (std::pow((sin(M_PI * x)), 3.)) +
                                         5.),
                               2.))) /
                              25. +
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 6. - 3. / 2.)) /
                              400. -
                          (3. *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 4. - 9. / 4.)) /
                              200. +
                          ((std::pow(M_PI, 2.)) * sin(2. * M_PI * t) *
                              (std::pow((sin(2. * M_PI * x)), 2.))) /
                              5. +
                          ((std::pow(M_PI, 2.)) * cos(2. * M_PI * x) *
                              cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
                              5. +
                          ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                              (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                                (cos(2. * M_PI * x) - 1.) -
                                            20.),
                                  2.)) *
                              ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                              (std::pow((cos(2. * M_PI * x) -
                                            2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                  2.))) /
                              40000.)) /
                  ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                          (std::pow(
                              (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                              2.)) +
                      100.) +
              ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                  (std::pow((sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                                    (cos(2. * M_PI * x) - 1.) -
                                20.),
                      2.)) *
                  ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (std::pow(
                      (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.), 2.))) /
                  (200. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                  (std::pow((cos(2. * M_PI * x) -
                                                2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                      2.)) +
                              100.)) +
              ((std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * sin(M_PI * x) * sin(3. * M_PI * x) *
                  sin(2. * M_PI * (t + 1. / 4.)) * ((7. * cos(4. * M_PI * t)) / 2. - 9. / 2.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * (std::pow((sin(M_PI * x)), 3.)) +
                      10.) *
                  (cos(2. * M_PI * x) - 2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.) *
                  (2. * cos(2. * M_PI * t) * cos(M_PI * x) * sin(M_PI * x) -
                      2. * cos(2. * M_PI * t) * (std::pow((cos(M_PI * x)), 3.)) * sin(M_PI * x) +
                      10.)) /
                  (25. * ((std::pow(M_PI, 2.)) * (std::pow((cos(2. * M_PI * t)), 2.)) *
                                 (std::pow((cos(2. * M_PI * x) -
                                               2. * (std::pow((cos(2. * M_PI * x)), 2.)) + 1.),
                                     2.)) +
                             100.))));

  switch (component)
  {
    case 0:
      return mu_ex;

    default:
      return 1.0;
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::WeaklyCompressibleEtienneFSIFluidViscosityFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    res[1] = 0.0;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    res[2] = 0.0;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiUP::BeltramiUP(const MAT::PAR::NewtonianFluid& fparams)
    : density_(fparams.density_), kinviscosity_(fparams.viscosity_ / fparams.density_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiUP::Evaluate(const double* xp, const double t, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = M_PI / 4.0;
  double b = M_PI / 4.0;
  double c = a * a + b * b + a * b;

  double K1 = exp(a * (x - z) + b * (y - z));  // =K4
  double K2 = exp(a * (z - y) + b * (x - y));  // =K5
  double K3 = exp(a * (y - x) + b * (z - x));  // =K6


  switch (component)
  {
    case 0:
      return b * K1 - a * K2;
    case 1:
      return b * K3 - a * K1;
    case 2:
      return b * K2 - a * K3;
    case 3:
      return c * (1.0 / K3 + 1.0 / K2 + 1.0 / K1) * density_;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiUP::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (component)
    {
      case 0:
        res[1] = 0;
      case 1:
        res[1] = 0;
      case 2:
        res[1] = 0;
      case 3:
        res[1] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (component)
    {
      case 0:
        res[2] = 0;
      case 1:
        res[2] = 0;
      case 2:
        res[2] = 0;
      case 3:
        res[2] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiGradU::BeltramiGradU(const MAT::PAR::NewtonianFluid& fparams)
    : kinviscosity_(fparams.viscosity_ / fparams.density_)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiGradU::Evaluate(const double* xp, const double t, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = M_PI / 4.0;
  double b = M_PI / 4.0;

  double K1 = exp(a * (x - z) + b * (y - z));  // =K4
  double K2 = exp(a * (z - y) + b * (x - y));  // =K5
  double K3 = exp(a * (y - x) + b * (z - x));  // =K6


  switch (component)
  {
    case 0:  // u,x
      return a * b * (K1 - K2);
    case 1:  // u,y
      return b * b * K1 + a * (a + b) * K2;
    case 2:  // u,z
      return -b * (a + b) * K1 - a * a * K2;
    case 3:  // v,x
      return -b * (a + b) * K3 - a * a * K1;
    case 4:  // v,y
      return a * b * K3 - a * b * K1;
    case 5:  // v,z
      return b * b * K3 + a * (a + b) * K1;
    case 6:  // w,x
      return b * b * K2 + a * (a + b) * K3;
    case 7:  // w,y
      return -b * (a + b) * K2 - a * a * K3;
    case 8:  // w,z
      return a * b * (K2 - K3);
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiGradU::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (component)
    {
      case 0:
        res[1] = 0;
      case 1:
        res[1] = 0;
      case 2:
        res[1] = 0;
      case 3:
        res[1] = 0;
      case 4:
        res[1] = 0;
      case 5:
        res[1] = 0;
      case 6:
        res[1] = 0;
      case 7:
        res[1] = 0;
      case 8:
        res[1] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (component)
    {
      case 0:
        res[2] = 0;
      case 1:
        res[2] = 0;
      case 2:
        res[2] = 0;
      case 3:
        res[2] = 0;
      case 4:
        res[2] = 0;
      case 5:
        res[2] = 0;
      case 6:
        res[2] = 0;
      case 7:
        res[2] = 0;
      case 8:
        res[2] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::BeltramiRHS::BeltramiRHS(const MAT::PAR::NewtonianFluid& fparams, bool is_stokes)
    : kinviscosity_(fparams.viscosity_ / fparams.density_), is_stokes_(is_stokes)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::BeltramiRHS::Evaluate(const double* xp, const double t, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  double z = xp[2];

  double a = M_PI / 4.0;
  double b = M_PI / 4.0;
  double c = a * a + b * b + a * b;

  double K1 = exp(a * (x - z) + b * (y - z));  // =K4
  double K2 = exp(a * (z - y) + b * (x - y));  // =K5
  double K3 = exp(a * (y - x) + b * (z - x));  // =K6

  double t1 = b * K1 - a * K2;
  double t2 = b * K3 - a * K1;
  double t3 = b * K2 - a * K3;

  double conv_x = 0.0;
  double conv_y = 0.0;
  double conv_z = 0.0;

  if (!is_stokes_)
  {
    conv_x = t1 * (a * b * K1 - a * b * K2) + t2 * (b * b * K1 + a * (a + b) * K2) +
             t3 * (-b * (a + b) * K1 - a * a * K2);
    conv_y = t1 * (-b * (a + b) * K3 - a * a * K1) + t2 * (a * b * K3 - a * b * K1) +
             t3 * (b * b * K3 + a * (a + b) * K1);
    conv_z = t1 * (b * b * K2 + a * (a + b) * K3) + t2 * (-b * (a + b) * K2 - a * a * K3) +
             t3 * (a * b * K2 - a * b * K3);
  }

  switch (component)
  {
    case 0:
      return c * ((a + b) / K3 - b / K2 - a / K1 - 2. * kinviscosity_ * (b * K1 - a * K2)) + conv_x;
    case 1:
      return c * (-a / K3 + (a + b) / K2 - b / K1 - 2. * kinviscosity_ * (b * K3 - a * K1)) +
             conv_y;
    case 2:
      return c * (-b / K3 - a / K2 + (a + b) / K1 - 2. * kinviscosity_ * (b * K2 - a * K3)) +
             conv_z;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::BeltramiRHS::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    switch (component)
    {
      case 0:
        res[1] = 0;
      case 1:
        res[1] = 0;
      case 2:
        res[1] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    switch (component)
    {
      case 0:
        res[2] = 0;
      case 1:
        res[2] = 0;
      case 2:
        res[2] = 0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinUP::KimMoinUP(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary)
    : density_(fparams.density_),
      kinviscosity_(fparams.viscosity_ / fparams.density_),
      is_stationary_(is_stationary)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinUP::Evaluate(const double* xp, const double t, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if (!is_stationary_)
  {
    gu = exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    gp = exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
  }

  switch (component)
  {
    case 0:
      return -cos(a_pi_x) * sin(a_pi_y) * gu;
    case 1:
      return sin(a_pi_x) * cos(a_pi_y) * gu;
    case 2:
      return 0.0;
    case 3:
      return -1. / 4. * (cos(2.0 * a_pi_x) + cos(2.0 * a_pi_y)) * gp * density_;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinUP::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if (!is_stationary_)
    {
      gu = (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
      gp = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    switch (component)
    {
      case 0:
        res[1] = -cos(a_pi_x) * sin(a_pi_y) * gu;
      case 1:
        res[1] = sin(a_pi_x) * cos(a_pi_y) * gu;
      case 2:
        res[1] = 0.0;
      case 3:
        res[1] = -1. / 4. * (cos(2.0 * a_pi_x) + cos(2.0 * a_pi_y)) * gp * density_;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if (!is_stationary_)
    {
      gu = (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
      gp = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    switch (component)
    {
      case 0:
        res[2] = -cos(a_pi_x) * sin(a_pi_y) * gu;
      case 1:
        res[2] = sin(a_pi_x) * cos(a_pi_y) * gu;
      case 2:
        res[2] = 0.0;
      case 3:
        res[2] = -1. / 4. * (cos(2.0 * a_pi_x) + cos(2.0 * a_pi_y)) * gp * density_;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinGradU::KimMoinGradU(const MAT::PAR::NewtonianFluid& fparams, bool is_stationary)
    : kinviscosity_(fparams.viscosity_ / fparams.density_), is_stationary_(is_stationary)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinGradU::Evaluate(const double* xp, const double t, const std::size_t component)
{
  // double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  // time factors
  double gu = 1.0;

  if (!is_stationary_)
  {
    gu = exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
  }


  switch (component)
  {
    case 0:  // u,x
      return sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
    case 1:  // u,y
      return -cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
    case 2:  // u,z
      return 0.0;
    case 3:  // v,x
      return cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
    case 4:  // v,y
      return -sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
    case 5:  // v,z
      return 0.0;
    case 6:  // w,x
      return 0.0;
    case 7:  // w,y
      return 0.0;
    case 8:  // w,z
      return 0.0;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinGradU::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  // double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;

    if (!is_stationary_)
    {
      gu = (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    switch (component)
    {
      case 0:  // u,x
        res[1] = sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
      case 1:  // u,y
        res[1] = -cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
      case 2:  // u,z
        res[1] = 0.0;
      case 3:  // v,x
        res[1] = cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
      case 4:  // v,y
        res[1] = -sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
      case 5:  // v,z
        res[1] = 0.0;
      case 6:  // w,x
        res[1] = 0.0;
      case 7:  // w,y
        res[1] = 0.0;
      case 8:  // w,z
        res[1] = 0.0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;

    if (!is_stationary_)
    {
      gu = (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           (-2.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    switch (component)
    {
      case 0:  // u,x
        res[2] = sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
      case 1:  // u,y
        res[2] = -cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
      case 2:  // u,z
        res[2] = 0.0;
      case 3:  // v,x
        res[2] = cos(a_pi_x) * cos(a_pi_y) * a * M_PI * gu;
      case 4:  // v,y
        res[2] = -sin(a_pi_x) * sin(a_pi_y) * a * M_PI * gu;
      case 5:  // v,z
        res[2] = 0.0;
      case 6:  // w,x
        res[2] = 0.0;
      case 7:  // w,y
        res[2] = 0.0;
      case 8:  // w,z
        res[2] = 0.0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinRHS::KimMoinRHS(
    const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, bool is_stokes)
    : kinviscosity_(fparams.viscosity_ / fparams.density_),
      is_stationary_(is_stationary),
      is_stokes_(is_stokes)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinRHS::Evaluate(const double* xp, const double t, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if (!is_stationary_)
  {
    gu = exp(-2.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    gp = exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
  }


  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  double conv_x = 0.0;
  double conv_y = 0.0;
  // double conv_z = 0.0;

  if (!is_stokes_)
  {
    conv_x = -a * M_PI * sin(a_pi_x) * cos(a_pi_x);
    conv_y = -a * M_PI * sin(a_pi_y) * cos(a_pi_y);
  }

  double visc_x = 0.0;
  double visc_y = 0.0;

  if (is_stationary_)
  {
    visc_x = kinviscosity_ * 2. * a * a * M_PI * M_PI *
             (-cos(a_pi_x) * sin(a_pi_y));  // * (gu = 1) for stationary
    visc_y = kinviscosity_ * 2. * a * a * M_PI * M_PI *
             (sin(a_pi_x) * cos(a_pi_y));  // * (gu = 1) for stationary
  }

  // in case of instationary: du/dt - \nu \laplacian(u) = 0

  switch (component)
  {
    case 0:
      return 0.5 * a * M_PI * sin(2. * a_pi_x) * gp + visc_x + conv_x * gu * gu;
    case 1:
      return 0.5 * a * M_PI * sin(2. * a_pi_y) * gp + visc_y + conv_y * gu * gu;
    case 2:
      return 0.0;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> FLD::KimMoinRHS::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component)
{
  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  double a = 2.0;

  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if (!is_stationary_)
    {
      gu = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
      gp = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    double a_pi_x = a * M_PI * x;
    double a_pi_y = a * M_PI * y;

    double conv_x = 0.0;
    double conv_y = 0.0;
    // double conv_z = 0.0;

    if (!is_stokes_)
    {
      conv_x = -a * M_PI * sin(a_pi_x) * cos(a_pi_x);
      conv_y = -a * M_PI * sin(a_pi_y) * cos(a_pi_y);
    }

    // in case of instationary: du/dt - \nu \laplacian(u) = 0

    switch (component)
    {
      case 0:
        res[1] = 0.5 * a * M_PI * sin(2. * a_pi_x) * gp + conv_x * gu;
      case 1:
        res[1] = 0.5 * a * M_PI * sin(2. * a_pi_y) * gp + conv_y * gu;
      case 2:
        res[1] = 0.0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // time factors
    double gu = 1.0;
    double gp = 1.0;

    if (!is_stationary_)
    {
      gu = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
      gp = (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           (-4.0 * a * a * M_PI * M_PI * kinviscosity_) *
           exp(-4.0 * a * a * M_PI * M_PI * t * kinviscosity_);
    }

    double a_pi_x = a * M_PI * x;
    double a_pi_y = a * M_PI * y;

    double conv_x = 0.0;
    double conv_y = 0.0;
    // double conv_z = 0.0;

    if (!is_stokes_)
    {
      conv_x = -a * M_PI * sin(a_pi_x) * cos(a_pi_x);
      conv_y = -a * M_PI * sin(a_pi_y) * cos(a_pi_y);
    }

    // in case of instationary: du/dt - \nu \laplacian(u) = 0

    switch (component)
    {
      case 0:
        res[2] = 0.5 * a * M_PI * sin(2. * a_pi_x) * gp + conv_x * gu;
      case 1:
        res[2] = 0.5 * a * M_PI * sin(2. * a_pi_y) * gp + conv_y * gu;
      case 2:
        res[2] = 0.0;
      default:
        dserror("wrong component %d", component);
        break;
    }
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::KimMoinStress::KimMoinStress(
    const MAT::PAR::NewtonianFluid& fparams, bool is_stationary, double amplitude)
    : kinviscosity_(fparams.viscosity_ / fparams.density_),
      density_(fparams.density_),
      is_stationary_(is_stationary),
      amplitude_(amplitude)
{
  if (amplitude_ != 1.0)
    dserror(
        "At the moment other implementation of Kimmoin Functions do not include the amplitude "
        "functionality!");
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FLD::KimMoinStress::Evaluate(const double* xp, const double t, const std::size_t component)
{
  // double visc = 1.0;

  double x = xp[0];
  double y = xp[1];
  // double z = xp[2];

  const double a = 2;  // in general does not need to be equal ... PARAM_A;
  const double b = 2;  // in general does not need to be equal ... PARAM_B;

  double a_pi_x = a * M_PI * x;
  double a_pi_y = a * M_PI * y;

  // time factors
  double gu = 1.0;
  double gp = 1.0;

  if (!is_stationary_)
  {
    gu = exp(-2.0 * b * b * M_PI * M_PI * t * kinviscosity_);
    gp = exp(-4.0 * b * b * M_PI * M_PI * t * kinviscosity_);
  }

  double fac =
      2 * kinviscosity_ * density_ * gu * M_PI * a * sin(a_pi_x) * sin(a_pi_y) * amplitude_;
  double p = -1. / 4. * (cos(2.0 * a_pi_x) + cos(2.0 * a_pi_y)) * gp * density_;


  switch (component)
  {
    case 0:  // sigma_xx
      return (fac - p);
    case 1:  // sigma_yy
      return (-fac - p);
    case 2:  // sigma_zz
      return (-p);
    case 3:  // sigma_xy
      return 0.0;
    case 4:  // sigma_yz
      return 0.0;
    case 5:  // sigma_zx
      return 0.0;
    default:
      dserror("wrong component %d", component);
      break;
  }

  return 1.0;
}