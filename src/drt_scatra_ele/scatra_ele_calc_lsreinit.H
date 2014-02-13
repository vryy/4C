/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_lsreinit.H

\brief evaluation of scatra elements for reinitialization equation

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_LSREINIT_H
#define SCATRA_ELE_CALC_LSREINIT_H

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcLsReinit : public ScaTraEleCalc<distype>
{
private:
    //! private constructor, since we are a Singleton.
    ScaTraEleCalcLsReinit(const int numdofpernode,const int numscal);

    typedef ScaTraEleCalc<distype> my;

public:

  virtual ~ScaTraEleCalcLsReinit() {}

  //! Singleton access method
  static ScaTraEleCalcLsReinit<distype> * Instance(
      const int numdofpernode,
      const int numscal,
      bool create=true );

  //! called upon destruction
  virtual void Done();

  /*========================================================================*/
  //! @name access routines
  /*========================================================================*/

  //! evaluate the element
  int Evaluate(DRT::ELEMENTS::Transport*     ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               const std::vector<int> &      lm,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                      Teuchos::ParameterList&       params,
                      DRT::Discretization &         discretization,
                      const std::vector<int> &      lm,
                      Epetra_SerialDenseMatrix&     elemat1_epetra,
                      Epetra_SerialDenseMatrix&     elemat2_epetra,
                      Epetra_SerialDenseVector&     elevec1_epetra,
                      Epetra_SerialDenseVector&     elevec2_epetra,
                      Epetra_SerialDenseVector&     elevec3_epetra);

protected:

  //! calculate matrix and rhs. Here the whole thing is hidden.
  void Sysmat(
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
  );

private:

  /*========================================================================*/
  //! @name overloaded methods for evaluation of individual terms
  /*========================================================================*/

  //! calculation of diffusive element matrix
  virtual void CalcMatDiff(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac,  //!< domain-integration factor times time-integration factor
    Teuchos::RCP<ScaTraEleDiffManager>  diffmanager  //!< diffusion manager
  );

  //! standard Galerkin diffusive term on right hand side
  void CalcRHSDiff(
    Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
    const int                     k,           //!< index of current scalar
    const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManager>  diffmanager,  //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1>&   gradphi       //!< scalar gradient at Gauss point
  );

  /*========================================================================*/
  //! @name additional method required for reinitialization
  /*========================================================================*/

  //! sign function
  void SignFunction(
    double&                           sign_phi,       //!< sign of phi
    const double                      charelelength,  //!< characteristic element length
    const double                      phizero,        //!< initial phi
    const LINALG::Matrix<my::nsd_,1>& gradphizero,    //!< gradient of initial phi
    const double                      phi,            //!< phi at time n+1
    const LINALG::Matrix<my::nsd_,1>& gradphi         //!< gradient of phi at time n+1
  );

  //! derivative of sign function
  void DerivSignFunction(
    double&                           deriv_sign,     //!< sign of phi
    const double                      charelelength,  //!< characteristic element length
    const double                      phizero         //!< initial phi
  );

  //! calculation of characteristic element length, i.e., interface thickness
  double CalcCharEleLengthReinit(
    const double                      vol,           //!< element volume
    const LINALG::Matrix<my::nsd_,1>& gradphizero   //!< gradient of initial phi
  );

  /*========================================================================*/
  //! @name penalty methods for reinitialization and related
  /*========================================================================*/

  //! calculation of element-wise denominator of penalty parameter
  void CalcElePenaltyParameter(
    double& penalty  //!< penalty parameter
  );

  //! calculate system matrix and rhs for correction step
  void SysmatCorrection(
    const double                          penalty, ///!< element penalty parameter
    Epetra_SerialDenseMatrix&             emat,    ///!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs     ///!< element rhs to calculate
  );

  //! calculation of penalty term on rhs
  void CalcRHSPenalty(
    Epetra_SerialDenseVector&     erhs,             //!< rhs vector
    const double                  fac,              //!< domain integration factor
    const double                  penalty,          //!< penalty parameter
    const double                  deriv_sign,       //!< derivative of sign function
    const double                  norm_gradphizero  //!< norm of gradient of initial phi
  );

  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  //! calculate system matrix and rhs for velocity projection
  void SysmatNodalVel(
  const int                             dir, ///< current spatial direction
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs ///< element rhs to calculate
  );

  /*========================================================================*/
  //! @name additional members
  /*========================================================================*/

  //! scalar at before reinitialization
  std::vector<LINALG::Matrix<my::nen_,1> > ephizero_;

};



/// Scatra diffusion manager
/*!

*/
template<int NSD>
class ScaTraEleDiffManagerLsReinit : public ScaTraEleDiffManager
{
  public:

    ScaTraEleDiffManagerLsReinit(int numscal):
    ScaTraEleDiffManager(numscal),
    diffdirectiontensor_(true),
    have_cross_wind_diff_(false)
    {return;}

    ~ScaTraEleDiffManagerLsReinit(){};

    void SetVelocityForCrossWindDiff(const LINALG::Matrix<NSD,1> velocity)
    {
      if (NSD != 3) dserror("Currently only 3d problems supported for crosswind diffusion");

      // compute tensor for anisotropic artificial diffusion
      // i.e., crosswind diffusion

      // get norm of velocity
      const double vel_norm_sq = velocity.Norm2() * velocity.Norm2();

      // compute tensor
      if (vel_norm_sq > 1.0e-8)
      {
        diffdirectiontensor_(0,0) = 1.0 - velocity(0,0) * velocity(0,0) / vel_norm_sq;
        diffdirectiontensor_(0,1) = - velocity(0,0) * velocity(1,0) / vel_norm_sq;
        diffdirectiontensor_(0,2) = - velocity(0,0) * velocity(2,0) / vel_norm_sq;
        diffdirectiontensor_(1,0) = diffdirectiontensor_(0,1);
        diffdirectiontensor_(1,1) = 1.0 - velocity(1,0) * velocity(1,0) / vel_norm_sq;
        diffdirectiontensor_(1,2) = - velocity(1,0) * velocity(2,0) / vel_norm_sq;
        diffdirectiontensor_(2,0) = diffdirectiontensor_(0,2);
        diffdirectiontensor_(2,1) = diffdirectiontensor_(1,2);
        diffdirectiontensor_(2,2) = 1.0 - velocity(2,0) * velocity(2,0) / vel_norm_sq;
      }
      else
       diffdirectiontensor_.Clear();

      // indicate that crosswind diffusion has to be used
      have_cross_wind_diff_ = true;
    }

    LINALG::Matrix<NSD,NSD> GetCrosswindTensor(){return diffdirectiontensor_;}

    bool HaveCrossWindDiff(){return have_cross_wind_diff_;}

  private:

    //! velocity for crosswind diffusion
    LINALG::Matrix<NSD,NSD> diffdirectiontensor_;

    //! flag for crosswind diffusion
    bool have_cross_wind_diff_;

};

}
}

#endif

