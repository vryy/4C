/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch_diffcond.H

\brief element parameter class for electrochemistry problems governed by diffusion-conduction formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H
#define SCATRA_ELE_PARAMETER_ELCH_DIFFCOND_H

#include "scatra_ele_parameter_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleParameterElchDiffCond : public ScaTraEleParameterElch
    {
      typedef DRT::ELEMENTS::ScaTraEleParameterElch myelch;

    public:

      //! destructor
      virtual ~ScaTraEleParameterElchDiffCond(){return;};

      //! singleton access method
      static ScaTraEleParameterElchDiffCond* Instance(const std::string& disname,bool create=true);

      //! singleton destruction
      virtual void Done();

      //! set parameters for electrochemistry problems with diffusion-conduction formulation
      void SetElementGeneralParameters(Teuchos::ParameterList& params);

      //! return flag for current as solution variable
      const bool CurSolVar() const {return cursolvar_;};

      //! return flag for diffusion potential
      const bool DiffusionCoeffBased() const {return diffusioncoefbased_;};

      //! return Newman constants
      const double NewmanConstA() const {return newmanconsta_;};
      const double NewmanConstB() const {return newmanconstb_;};
      const double NewmanConstC() const {return newmanconstc_;};

    private:

      //! private constructor for singletons
      ScaTraEleParameterElchDiffCond(const std::string& disname);

      //! flag if current is used as a solution variable
      bool cursolvar_;

      //! mat_diffcond: flag if diffusion potential is based on diffusion coefficients or transference number
      bool diffusioncoefbased_;

      //! switch for dilute and concentrated solution theory (diffusion potential in current equation):
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
  } // namespace ELEMENTS
} // namespace DRT
#endif
