/*----------------------------------------------------------------------*/
/*! \file
\brief auxiliary material for macro-scale elements in multi-scale simulations of scalar transport
problems. This material handles the communication between micro and macro materials

The functions implemented in this file have to be separated from the remainder
of the ScatraMicroMacroCoupling class in the files
baci_mat_scatra_micro_macro_coupling.hpp.{H;cpp}.
The reason is that the files scatra_mat_multiscale_gp.{H;cpp} are not compiled
when building the preprocessing and postprocessing filters (cf. CMakeLists.txt),
and hence the ScatraMultiScaleGP class is not available. However, that class
is used by all functions in this file, and therefore they are reimplemented as
empty dummy functions throwing dserrors in filter_commmon/filter_evaluation.cpp.
When building the filters, that file is compiled instead of this one, such that
the ScatraMultiScaleGP class is not required and the linker is happy. The
dummy functions are never called by the filters anyway, and if they will be one
day, the corresponding dserrors will indicate that the current code structure
needs to be reworked. If functions are added, removed, or modified in this file,
the file filter_commmon/filter_evaluation.cpp needs to be adapted accordingly.

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra_micro_macro_coupling.hpp"
#include "4C_mat_scatra_multiscale_gp.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::Initialize(const int ele_id, const int gp_id, const bool is_ale)
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // initialize multi-scale scalar transport material
  if (matgp_.find(gp_id) == matgp_.end())
  {
    // instantiate and initialize multi-scale scalar transport submaterial at macro-scale Gauss
    // point
    matgp_[gp_id] = Teuchos::rcp(new ScatraMultiScaleGP(ele_id, gp_id, MicroDisNum(), is_ale));
    matgp_[gp_id]->Init();
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::PrepareTimeStep(
    const int gp_id, const std::vector<double>& phinp_macro) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // prepare time step on micro scale
  matgp_.at(gp_id)->PrepareTimeStep(phinp_macro);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::Evaluate(const int gp_id,
    const std::vector<double>& phinp_macro, double& q_micro, std::vector<double>& dq_dphi_micro,
    const double detF, const bool solve) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // evaluate multi-scale scalar transport sub material at macro-scale Gauss point
  matgp_.at(gp_id)->Evaluate(phinp_macro, q_micro, dq_dphi_micro, detF, solve);
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::ScatraMicroMacroCoupling::EvaluateMeanConcentration(const int gp_id) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // evaluate mean concentration on micro scale
  return matgp_.at(gp_id)->EvaluateMeanConcentration();
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
double MAT::ScatraMicroMacroCoupling::EvaluateMeanConcentrationTimeDerivative(const int gp_id) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // evaluate mean concentration time derivative on micro scale
  return matgp_.at(gp_id)->EvaluateMeanConcentrationTimeDerivative();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::Update(const int gp_id) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // update multi-scale scalar transport submaterial at macro-scale Gauss point
  matgp_.at(gp_id)->Update();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::Output(const int gp_id) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // create output on micro scale
  matgp_.at(gp_id)->Output();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::ReadRestart(const int gp_id) const
{
  // safety check
  if (gp_id < 0) FOUR_C_THROW("Invalid macro-scale Gauss point ID!");

  // read restart on micro scale
  matgp_.at(gp_id)->ReadRestart();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMicroMacroCoupling::SetTimeStepping(
    const int gp_id, const double dt, const double time, const int step)
{
  matgp_.at(gp_id)->SetTimeStepping(dt, time, step);
}

FOUR_C_NAMESPACE_CLOSE
