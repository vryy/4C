/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation of multi-scale scalar transport material

The functions implemented in this file have to be separated from the remainder
of the ScatraMultiScale class in the files scatra_mat_multiscale.{H;cpp}.
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

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#include "scatra_multiscale.H"

#include "scatra_multiscale_gp.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

/*--------------------------------------------------------------------*
 | initialize multi-scale scalar transport material        fang 02/16 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Initialize(const int ele_id,  //!< macro-scale element ID
    const int gp_id                                       //!< macro-scale Gauss point ID
)
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // initialize multi-scale scalar transport material
  if (matgp_.find(gp_id) == matgp_.end())
  {
    // instantiate and initialize multi-scale scalar transport submaterial at macro-scale Gauss
    // point
    matgp_[gp_id] = Teuchos::rcp(new ScatraMultiScaleGP(ele_id, gp_id, MicroDisNum()));
    matgp_[gp_id]->Init();
  }

  return;
}


/*--------------------------------------------------------------------*
 | prepare time step on micro scale                        fang 02/16 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::PrepareTimeStep(const int gp_id,  //!< macro-scale Gauss point ID
    const std::vector<double>& phinp_macro                    //!< macro-scale state variables
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // prepare time step on micro scale
  matgp_.at(gp_id)->PrepareTimeStep(phinp_macro);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate multi-scale scalar transport material          fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Evaluate(const int gp_id,  //!< macro-scale Gauss point ID
    const std::vector<double>& phinp_macro,            //!< macro-scale state variables
    double& q_micro,                                   //!< micro-scale flux
    std::vector<double>&
        dq_dphi_micro,  //!< derivatives of micro-scale flux w.r.t. macro-scale state variables
    const bool solve    //!< flag indicating whether micro-scale problem should be solved
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // evaluate multi-scale scalar transport submaterial at macro-scale Gauss point
  matgp_.at(gp_id)->Evaluate(phinp_macro, q_micro, dq_dphi_micro, solve);

  return;
}


/*--------------------------------------------------------------------*
 | evaluate mean concentration on micro scale              fang 08/17 |
 *--------------------------------------------------------------------*/
double MAT::ScatraMultiScale::EvaluateMeanConcentration(
    const int gp_id  //!< macro-scale Gauss point ID
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // evaluate mean concentration on micro scale
  return matgp_.at(gp_id)->EvaluateMeanConcentration();
}


/*-------------------------------------------------------------------------*
 | evaluate mean concentration time derivative on micro scale   fang 03/18 |
 *-------------------------------------------------------------------------*/
double MAT::ScatraMultiScale::EvaluateMeanConcentrationTimeDerivative(
    const int gp_id  //!< macro-scale Gauss point ID
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // evaluate mean concentration time derivative on micro scale
  return matgp_.at(gp_id)->EvaluateMeanConcentrationTimeDerivative();
}


/*--------------------------------------------------------------------*
 | update multi-scale scalar transport material            fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Update(const int gp_id  //!< macro-scale Gauss point ID
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // update multi-scale scalar transport submaterial at macro-scale Gauss point
  matgp_.at(gp_id)->Update();

  return;
}


/*--------------------------------------------------------------------*
 | create output on micro scale                            fang 02/16 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Output(const int gp_id  //!< macro-scale Gauss point ID
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // create output on micro scale
  matgp_.at(gp_id)->Output();

  return;
}


/*--------------------------------------------------------------------*
 | read restart on micro scale                             fang 03/16 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::ReadRestart(const int gp_id  //!< macro-scale Gauss point ID
    ) const
{
  // safety check
  if (gp_id < 0) dserror("Invalid macro-scale Gauss point ID!");

  // read restart on micro scale
  matgp_.at(gp_id)->ReadRestart();

  return;
}
