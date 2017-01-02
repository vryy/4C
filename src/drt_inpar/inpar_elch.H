/*----------------------------------------------------------------------*/
/*!
\file inpar_elch.H

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
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
}

namespace INPAR
{
  /// Input parameters for ELCH
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

    // faraday constant (general definition)
    const double faraday_const = 96485.3399;

    // gas constant (general definition)
    const double gas_const = 8.314472;

    // permittivität (not a general definition, but only once used)
    const double epsilon_const = 1.e-4;

    /// set the elch parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific elch conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);
  }  // namespace ELCH
}  // namespace INPAR
#endif  // #ifndef INPAR_ELCH_H
