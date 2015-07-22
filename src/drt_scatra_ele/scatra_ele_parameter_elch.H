/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.H

\brief element parameter class for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_ELCH_H
#define SCATRA_ELE_PARAMETER_ELCH_H

#include "../drt_inpar/inpar_elch.H"

#include "scatra_ele_parameter.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    class ScaTraEleParameterElch : public ScaTraEleParameter
    {
      typedef DRT::ELEMENTS::ScaTraEleParameter my;

    public:

      //! destructor
      virtual ~ScaTraEleParameterElch(){return;};

      //! singleton access method
      static ScaTraEleParameterElch* Instance(const std::string& disname, bool create=true);

      //! singleton destruction
      virtual void Done();

      //! set parameters for electrochemistry problems
      virtual void SetElementGeneralParameters(Teuchos::ParameterList& params);

      //! return type of closing equation for electric potential
      const INPAR::ELCH::EquPot EquPot() const {return equpot_;};

      //! return Faraday constant
      const double Faraday() const {return faraday_;};

      //! return dielectric constant
      const double Epsilon() const {return epsilon_;};

      //! return constant F/RT
      const double FRT() const {return frt_;};

    protected:

      //! protected constructor for singletons
      ScaTraEleParameterElch(const std::string& disname);

      //! equation used for closing of the elch-system
      enum INPAR::ELCH::EquPot equpot_;

      //! Faraday constant
      const double faraday_;

      //! dielectric constant
      const double epsilon_;

      //! pre-calculation of regularly used constant F/RT
      //! (a division is much more expensive than a multiplication)
      double frt_;
    };
  } // namespace ELEMENTS
} // namespace DRT
#endif
