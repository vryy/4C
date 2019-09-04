/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static diffusion-conduction parameters required for element
evaluation

This singleton class holds all static diffusion-conduction parameters required for element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with more general parameter classes holding
additional static parameters required for scalar transport element evaluation.

\maintainer Christoph Schmidt

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H
#define SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H

#include "scatra_ele_parameter_base.H"

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
      //! destructor
      virtual ~ScaTraEleParameterElchDiffCond() { return; };

      //! singleton access method
      static ScaTraEleParameterElchDiffCond* Instance(
          const std::string& disname,  //!< name of discretization
          const ScaTraEleParameterElchDiffCond* delete_me =
              NULL  //!< creation/destruction indication
      );

      //! singleton destruction
      virtual void Done();

      //! set parameters
      void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
      );

      //! return flag for current as solution variable
      bool CurSolVar() const { return cursolvar_; };

      //! return flag for diffusion potential
      bool DiffusionCoeffBased() const { return diffusioncoefbased_; };

      //! return Newman constants
      double NewmanConstA() const { return newmanconsta_; };
      double NewmanConstB() const { return newmanconstb_; };
      double NewmanConstC() const { return newmanconstc_; };

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
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
