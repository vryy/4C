/*-----------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of functions for structure problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_functions.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_stvenantkirchhoff.H"
#include "baci_utils_function_manager.H"

BACI_NAMESPACE_OPEN

namespace
{
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

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  Teuchos::RCP<CORE::UTILS::FunctionOfSpaceTime> CreateStructureFunction(
      const std::vector<INPUT::LineDefinition>& function_line_defs)
  {
    if (function_line_defs.size() != 1) return Teuchos::null;

    const auto& function_lin_def = function_line_defs.front();

    if (function_lin_def.HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE"))
    {
      // read data
      int mat_id_struc = -1;

      function_lin_def.ExtractInt("MAT_STRUC", mat_id_struc);

      if (mat_id_struc <= 0)
        dserror(
            "Please give a (reasonable) 'MAT_STRUC' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE");

      // get materials
      auto fparams = GetSVKMatPars(mat_id_struc);

      return Teuchos::rcp(new STR::WeaklyCompressibleEtienneFSIStructureFunction(fparams));
    }
    else if (function_lin_def.HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE"))
    {
      // read data
      int mat_id_struc = -1;

      function_lin_def.ExtractInt("MAT_STRUC", mat_id_struc);

      if (mat_id_struc <= 0)
      {
        dserror(
            "Please give a (reasonable) 'MAT_STRUC' in "
            "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE");
      }

      // get materials
      auto fparams = GetSVKMatPars(mat_id_struc);

      return Teuchos::rcp(new STR::WeaklyCompressibleEtienneFSIStructureForceFunction(fparams));
    }
    else
    {
      return Teuchos::RCP<CORE::UTILS::FunctionOfSpaceTime>(nullptr);
    }
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::AddValidStructureFunctions(CORE::UTILS::FunctionManager& function_manager)
{
  std::vector<INPUT::LineDefinition> lines;
  lines.emplace_back(INPUT::LineDefinition::Builder()
                         .AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE")
                         .AddNamedInt("MAT_STRUC")
                         .Build());
  lines.emplace_back(INPUT::LineDefinition::Builder()
                         .AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE")
                         .AddNamedInt("MAT_STRUC")
                         .Build());

  function_manager.AddFunctionDefinition(std::move(lines), CreateStructureFunction);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STR::WeaklyCompressibleEtienneFSIStructureFunction::WeaklyCompressibleEtienneFSIStructureFunction(
    const MAT::PAR::StVenantKirchhoff& fparams)
    : poissonratio_(fparams.poissonratio_)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double STR::WeaklyCompressibleEtienneFSIStructureFunction::Evaluate(
    const double* xp, const double t, const std::size_t component) const
{
  // ease notation
  double x = xp[0];
  double y = xp[1];

  // initialize variables
  CORE::LINALG::Matrix<2, 1> u_ex;

  // evaluate variables
  u_ex(0) = -((cos(2. * M_PI * t) * cos(2. * M_PI * x)) / 6. + 1.) * (y - 1.);
  u_ex(1) =
      -(sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) / 20.;

  switch (component)
  {
    case 0:
      return u_ex(0);
    case 1:
      return u_ex(1);

    default:
      return 1.0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> STR::WeaklyCompressibleEtienneFSIStructureFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component) const
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(xp, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // ease notation
    double x = xp[0];
    double y = xp[1];

    // initialize variables
    CORE::LINALG::Matrix<2, 1> dudt_ex;

    // evaluate variables
    dudt_ex(0) = (M_PI * cos(2. * M_PI * x) * sin(2. * M_PI * t) * (y - 1.)) / 3.;
    dudt_ex(1) =
        -(M_PI * sin(2. * M_PI * x) * cos(2. * M_PI * (t + 1. / 4.)) * (cos(2. * M_PI * x) - 1.)) /
        10.;

    switch (component)
    {
      case 0:
        res[1] = dudt_ex(0);
        break;
      case 1:
        res[1] = dudt_ex(1);
        break;

      default:
        res[1] = 0.0;
    }
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
STR::WeaklyCompressibleEtienneFSIStructureForceFunction::
    WeaklyCompressibleEtienneFSIStructureForceFunction(const MAT::PAR::StVenantKirchhoff& fparams)
    : youngmodulus_(fparams.youngs_),
      poissonratio_(fparams.poissonratio_),
      strucdensity_(fparams.density_)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double STR::WeaklyCompressibleEtienneFSIStructureForceFunction::Evaluate(
    const double* xp, const double t, const std::size_t component) const
{
  // ease notation
  double x = xp[0];
  double y = xp[1];
  double E = youngmodulus_;
  double v = poissonratio_;
  double r = strucdensity_;

  // initialize variables
  CORE::LINALG::Matrix<2, 1> f_u_ex;

  // evaluate variables
  f_u_ex(0) = (2. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * t) * cos(2. * M_PI * x) * (y - 1.) *
                  (E - r - E * v + r * v + 2. * r * (std::pow(v, 2.)))) /
              (3. * (2. * (std::pow(v, 2.)) + v - 1.));
  f_u_ex(1) = ((std::pow(M_PI, 2.)) * r * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                  (cos(2. * M_PI * x) - 1.)) /
                  5. -
              (E * ((M_PI * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.))) / 3. +
                       (3. * (std::pow(M_PI, 2.)) * cos(2. * M_PI * x) * sin(2. * M_PI * x) *
                           sin(2. * M_PI * (t + 1. / 4.))) /
                           5. +
                       ((std::pow(M_PI, 2.)) * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.)) *
                           (cos(2. * M_PI * x) - 1.)) /
                           5.)) /
                  (2. * v + 2.) +
              (E * M_PI * v * sin(2. * M_PI * x) * sin(2. * M_PI * (t + 1. / 4.))) /
                  (3. * (2. * v - 1.) * (v + 1.));

  switch (component)
  {
    case 0:
      return f_u_ex(0);
    case 1:
      return f_u_ex(1);

    default:
      return 1.0;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> STR::WeaklyCompressibleEtienneFSIStructureForceFunction::EvaluateTimeDerivative(
    const double* xp, const double t, const unsigned deg, const std::size_t component) const
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

BACI_NAMESPACE_CLOSE
