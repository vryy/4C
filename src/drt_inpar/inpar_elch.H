/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electrochemistry

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_ELCH_H
#define INPAR_ELCH_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

namespace INPAR
{
  namespace ELCH
  {
    /// possible types of moving boundary simulation
    enum ElchMovingBoundary
    {
      elch_mov_bndry_no,
      elch_mov_bndry_pseudo_transient,
      elch_mov_bndry_fully_transient,
    };

    /// type of closing equation for electric potential
    enum EquPot
    {
      equpot_undefined,
      equpot_enc,
      equpot_enc_pde,
      equpot_enc_pde_elim,
      equpot_poisson,
      equpot_laplace,
      equpot_divi
    };

    /// type of electrode kinetics
    enum ElectrodeKinetics
    {
      butler_volmer,
      butler_volmer_yang1997,
      tafel,
      linear,
      butler_volmer_newman,
      butler_volmer_bard,
      nernst,
      zero
    };

    enum DiffCondMat
    {
      diffcondmat_newman,
      diffcondmat_ion,
      diffcondmat_undefined
    };

    enum ApproxElectResist
    {
      approxelctresist_relpotcur,
      approxelctresist_effleninitcond,
      approxelctresist_efflenintegcond
    };

    enum class CCCVHalfCyclePhase
    {
      undefined,           //!< undefined mode
      constant_current,    //!< constant-current (CC) mode
      constant_voltage,    //!< constant-voltage (CV) mode
      relaxation,          //!< relaxation (RX) mode
      initital_relaxation  //!< initial relaxation mode
    };

    // permittivität (not a general definition, but only once used)
    const double epsilon_const = 1.e-4;

    /// set the elch parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific elch conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace ELCH
}  // namespace INPAR
#endif  // #ifndef INPAR_ELCH_H
