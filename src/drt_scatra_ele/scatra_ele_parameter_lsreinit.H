/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_lsreinit.H

\brief Setting of scatra parameter for element evaluation of reinitialization equation

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_LSREINIT_H
#define SCATRA_ELE_PARAMETER_LSREINIT_H

#include "scatra_ele_parameter.H"

#include "../drt_inpar/inpar_levelset.H"


namespace DRT
{
namespace ELEMENTS
{

class ScaTraEleParameterLsReinit : public ScaTraEleParameter
{
protected:

  /// (private) protected constructor, since we are a Singleton.
  /// this constructor is called from a derived class
  /// -> therefore, it has to be protected instead of private
  ScaTraEleParameterLsReinit();

public:

  virtual ~ScaTraEleParameterLsReinit() {}

  /// Singleton access method
  static ScaTraEleParameterLsReinit* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /*========================================================================*/
  //! @name set-routines
  /*========================================================================*/

  //! general fluid parameter are set
  void SetElementLsReinitScaTraParameter(Teuchos::ParameterList& params);   //> parameter list

  /*========================================================================*/
  //! @name access-routines
  /*========================================================================*/

  INPAR::SCATRA::SmoothedSignType SignType(){return signtype_;};
  INPAR::SCATRA::CharEleLengthReinit CharEleLengthReinit(){return charelelengthreinit_;};
  double InterfaceThicknessFac(){return interfacethicknessfac_;};
  bool UseProjectedVel(){return useprojectedreinitvel_;};
  INPAR::SCATRA::LinReinit LinForm(){return linform_;};
  INPAR::SCATRA::ArtDiff ArtDiff(){return artdiff_;};

private:

  // sign function for phi
  INPAR::SCATRA::SmoothedSignType signtype_;
  // element length for smoothing of sign function
  INPAR::SCATRA::CharEleLengthReinit charelelengthreinit_;
  // interface thickness factor (multiple of characteristic element length)
  double interfacethicknessfac_;
  // from of velocity evaluation
  bool useprojectedreinitvel_;
  // form of linearization of nonlinear terms
  INPAR::SCATRA::LinReinit linform_;
  // form of artificial diffusion
  INPAR::SCATRA::ArtDiff artdiff_;

};

}
}

#endif
