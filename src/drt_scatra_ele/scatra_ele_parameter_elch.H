/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_elch.H

\brief Setting of elch scatra parameter for element evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_ELCH_H
#define SCATRA_ELE_PARAMETER_ELCH_H

#include "scatra_ele_parameter.H"

#include "../drt_inpar/inpar_elch.H"


namespace DRT
{
namespace ELEMENTS
{

class ScaTraEleParameterElch : public ScaTraEleParameter
{

public:

  virtual ~ScaTraEleParameterElch() {}

  /// Singleton access method
  static ScaTraEleParameterElch* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /*========================================================================*/
  //! @name set-routines
  /*========================================================================*/

  //! elch-specific parameters are set
  void SetElementElchScaTraParameter(Teuchos::ParameterList& params,   //> parameter list
                                     int myrank);                      //> proc id

  //! elch diffusion-conduction formulation specific parameters are set
  void SetElementElchDiffCondScaTraParameter(
                                     Teuchos::ParameterList& params,   //> parameter list
                                     int myrank);                      //> proc id

  /*========================================================================*/
  //! @name access-routines
  /*========================================================================*/

  /*--------------------------------------------------------*/
  //! @name general elch parameters for Nernst-Planck equation
  /*--------------------------------------------------------*/

  INPAR::ELCH::ElchType ElchType(){return elchtype_;};

  bool NernstPlanck(){return nernstplanck_;};

  double FRT(){return frt_;};

  /*--------------------------------------------------------*/
  //! @name elch parameters for Nernst-Planck equation
  /*--------------------------------------------------------*/

  /*--------------------------------------------------------*/
  //! @name elch parameters for diffusion-conduction equation
  /*--------------------------------------------------------*/

  bool CurSolVar(){return cursolvar_;};

  bool DiffusionCoeffBased(){return diffusioncoefbased_;};

  INPAR::ELCH::EquPot EquPot(){return equpot_;};

  double NewmanConstA(){return newmanconsta_;};

  double NewmanConstB(){return newmanconstb_;};

  double NewmanConstC(){return newmanconstc_;};

protected:

  /// protected Constructor since we are a Singleton.
  ScaTraEleParameterElch();

  //! type of elch problem
  INPAR::ELCH::ElchType elchtype_;

  //! activation of Nernst-Planck or diffusion-conduction formulation
  bool nernstplanck_;

  //! factor frt
  double frt_;

  bool cursolvar_;

  bool diffusioncoefbased_;

  enum INPAR::ELCH::EquPot equpot_;


  /// dilute solution theory (diffusion potential in current equation):
  ///    A          B
  ///   |--|  |----------|
  ///   z_1 + (z_2 - z_1) t_1
  /// ------------------------ (RT/F kappa 1/c_k grad c_k)
  ///      z_1 z_2
  ///     |________|
  ///         C
  double newmanconsta_;

  double newmanconstb_;

  double newmanconstc_;

private:

};

}
}

#endif
