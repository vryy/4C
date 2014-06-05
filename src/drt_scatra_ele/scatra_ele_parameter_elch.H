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
protected:

  /// (private) protected constructor, since we are a Singleton.
  /// this constructor is called from a derived class
  /// -> therefore, it has to be protected instead of private
  ScaTraEleParameterElch();

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
  //! @name general elch parameters
  /*--------------------------------------------------------*/

  INPAR::ELCH::ElchType ElchType(){return elchtype_;};

  INPAR::ELCH::EquPot EquPot(){return equpot_;};

  double FRT(){return frt_;};

  /*--------------------------------------------------------------*/
  //! @name elch parameters for Nernst-Planck formulation
  /*--------------------------------------------------------------*/

  /*--------------------------------------------------------------*/
  //! @name elch parameters for diffusion-conduction formulation
  /*--------------------------------------------------------------*/

  bool CurSolVar(){return cursolvar_;};

  bool DiffusionCoeffBased(){return diffusioncoefbased_;};

  const double NewmanConstA(){return newmanconsta_;};

  const double NewmanConstB(){return newmanconstb_;};

  const double NewmanConstC(){return newmanconstc_;};

private:

  //! type of elch problem
  enum INPAR::ELCH::ElchType elchtype_;

  //! equation used for closing of the elch-system
  enum INPAR::ELCH::EquPot equpot_;

  //! factor frt
  double frt_;

  //! flag if current is used as a solution variable
  bool cursolvar_;

  // mat_diffcond: flag if diffusion potential is based on diffusion coefficients or transference number
  bool diffusioncoefbased_;

  /// switch for dilute and concentrated solution theory (diffusion potential in current equation):
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
};

} // namespace ELEMENTS
} // namespace DRT

#endif
