/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static diffusion-conduction parameters required for element
evaluation

This singleton class holds all static diffusion-conduction parameters required for element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with more general parameter classes holding
additional static parameters required for scalar transport element evaluation.


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H
#define BACI_SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H

#include "baci_scatra_ele_parameter_base.H"

#include <string>

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterStd;

    // class implementation
    class ScaTraEleParameterElchDiffCond : public ScaTraEleParameterBase
    {
     public:
      //! singleton access method
      static ScaTraEleParameterElchDiffCond* Instance(
          const std::string& disname  //!< name of discretization
      );

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) override;

      //! return flag for current as solution variable
      bool CurSolVar() const { return cursolvar_; }

      //! return flag for diffusion potential
      bool DiffusionCoeffBased() const { return diffusioncoefbased_; }

      //! return Newman constants
      double NewmanConstA() const { return newmanconsta_; }
      double NewmanConstB() const { return newmanconstb_; }
      double NewmanConstC() const { return newmanconstc_; }
      double Epsilon0() const { return epsilon_0_; }

     private:
      //! private constructor for singletons
      ScaTraEleParameterElchDiffCond(const std::string& disname  //!< name of discretization
      );

      //! flag if current is used as a solution variable
      bool cursolvar_;

      //! mat_diffcond: flag if diffusion potential is based on diffusion coefficients or
      //! transference number
      bool diffusioncoefbased_;

      //! switch for dilute and concentrated solution theory (diffusion potential in current
      //! equation):
      //!    A          B
      //!   |--|  |----------|
      //!   z_1 + (z_2 - z_1) t_1
      //! ------------------------ (RT/F kappa 1/c_k grad c_k)
      //!      z_1 z_2
      //!     |________|
      //!         C
      double newmanconsta_;
      double newmanconstb_;
      double newmanconstc_;
      double epsilon_0_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
