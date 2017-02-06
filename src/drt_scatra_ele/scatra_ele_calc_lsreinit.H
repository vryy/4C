/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_lsreinit.H

\brief evaluation of scatra elements for reinitialization equation

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_LSREINIT_H
#define SCATRA_ELE_CALC_LSREINIT_H

#include "scatra_ele_calc.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/geo_utils.H"

namespace DRT
{
namespace ELEMENTS
{
// forward declarations
template<int NSD> class ScaTraEleDiffManagerLsReinit;
template<int NSD, int NEN> class ScaTraEleInternalVariableManagerLsReinit;
class ScaTraEleParameterLsReinit;

template<DRT::Element::DiscretizationType distype, unsigned probDim>
class ScaTraEleCalcLsReinit : public ScaTraEleCalc<distype,probDim>
{
private:
    //! private constructor for singletons
    ScaTraEleCalcLsReinit(
        const int numdofpernode,
        const int numscal,
        const std::string& disname
        );

    typedef ScaTraEleCalc<distype,probDim> my;

public:

  virtual ~ScaTraEleCalcLsReinit() {}

  //! Singleton access method
  static ScaTraEleCalcLsReinit<distype, probDim> * Instance(
      const int numdofpernode,
      const int numscal,
      const std::string& disname,
      const ScaTraEleCalcLsReinit* delete_me = NULL );

  //! called upon destruction
  virtual void Done();

  /*========================================================================*/
  //! @name access routines
  /*========================================================================*/

  //! setup element evaluation
  int SetupCalc(
      DRT::Element*               ele,
      DRT::Discretization&        discretization
      );

  //! evaluate the element
  int Evaluate(DRT::Element*                 ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               DRT::Element::LocationArray&  la,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  int EvaluateAction(DRT::Element*                 ele,
                     Teuchos::ParameterList&       params,
                     DRT::Discretization &         discretization,
                     const SCATRA::Action&         action,
                     DRT::Element::LocationArray&  la,
                     Epetra_SerialDenseMatrix&     elemat1_epetra,
                     Epetra_SerialDenseMatrix&     elemat2_epetra,
                     Epetra_SerialDenseVector&     elevec1_epetra,
                     Epetra_SerialDenseVector&     elevec2_epetra,
                     Epetra_SerialDenseVector&     elevec3_epetra);

protected:

  //! calculate matrix and rhs. Here the whole thing is hidden. Hyperbolic reinit.
  void SysmatHyperbolic(
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
  );

  //! calculate matrix and rhs. Here the whole thing is hidden. Ellipitic reinit.
  void SysmatElliptic(
    Epetra_SerialDenseMatrix&         emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&         erhs,      //!< element rhs to calculate
    const GEO::BoundaryIntCellPtrs &  bcell    //!< interface for penalty term
  );

  /** \brief Evaluate the system matrix and right-hand-side for the ellipitic
   *  reinitialization
   *
   *  Build the system for a Newton Raphson scheme.
   *
   *  \param emat  (out) : element matrix (part of the tangential matrix)
   *  \param erhs  (out) : element right-hand-side vector (part of the residual vector)
   *  \param bcell (in)  : boundary integration cell (necessary for the penalty term)
   *
   *  \author hiermeier \date 12/16 */
  void EllipticNewtonSystem(
      Epetra_SerialDenseMatrix *             emat,
      Epetra_SerialDenseVector *             erhs,
      const LINALG::Matrix<my::nen_,1> &     el2sysmat_diag_inv,
      const GEO::BoundaryIntCellPtrs &       bcell );

private:

  /*========================================================================*/
  //! @name general evaluation methods
  /*========================================================================*/

  void EvalReinitialization(
      const Epetra_Vector &          phinp,
      const std::vector<int> &       lm,
      DRT::Element *                 ele,
      Teuchos::ParameterList &       params,
      DRT::Discretization &          discretization,
      Epetra_SerialDenseMatrix &     elemat1_epetra,
      Epetra_SerialDenseVector &     elevec1_epetra);

  void EvalReinitializationEmbedded(
      const std::vector<int> &       lm,
      DRT::Element *                 ele,
      Teuchos::ParameterList &       params,
      DRT::Discretization &          discretization,
      Epetra_SerialDenseMatrix &     elemat1_epetra,
      Epetra_SerialDenseVector &     elevec1_epetra);

  void EvalReinitializationStd(
      const Epetra_Vector &          phinp,
      const std::vector<int> &       lm,
      DRT::Element *                 ele,
      Teuchos::ParameterList &       params,
      DRT::Discretization &          discretization,
      Epetra_SerialDenseMatrix &     elemat1_epetra,
      Epetra_SerialDenseVector &     elevec1_epetra);

  /*========================================================================*/
  //! @name overloaded methods for evaluation of individual terms
  /*========================================================================*/

  //! calculation of diffusive element matrix
  virtual void CalcMatDiff(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac   //!< domain-integration factor times time-integration factor
  );

  //! standard Galerkin diffusive term on right hand side
  void CalcRHSDiff(
    Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
    const int                     k,           //!< index of current scalar
    const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    const LINALG::Matrix<my::nsd_,1>&   gradphi       //!< scalar gradient at Gauss point
  );

  /*========================================================================*/
  //! @name additional methods required for reinitialization
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

  //! calculation of interface penalty term for elliptic reinitialization
  void EvaluateInterfaceTerm(
    Epetra_SerialDenseMatrix *             emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector *             erhs,     //!< element vector to calculate
    const GEO::BoundaryIntCellPtrs &       boucell   //!< interface for penalty term
  );

  //! calculation of interface penalty term for elliptic reinitialization (gauss loop)
  template <DRT::Element::DiscretizationType celldistype>
  void CalcPenaltyTerm(
    Epetra_SerialDenseMatrix &             emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector &             erhs,     //!< element vector to calculate
    const GEO::BoundaryIntCell &           cell      //!< interface cell
  );


  /** \brief calculation of interface penalty term for elliptic reinitialization
   *         (special variant for the 0-D case)
   *
   *  \param emat (out) : element matrix to calculate
   *  \param erhs (out) : element vector to calculate
   *  \param cell (in)  : interface boundary integration cell
   *
   *  \author hiermeier \date 11/16 */
  void CalcPenaltyTerm_0D(
      Epetra_SerialDenseMatrix *            emat,
      Epetra_SerialDenseVector *            erhs,
      const GEO::BoundaryIntCell &          cell
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

  //! get diffusion manager for reinitialization
  Teuchos::RCP<ScaTraEleDiffManagerLsReinit<my::nsd_> > DiffManager(){return Teuchos::rcp_static_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(my::diffmanager_);};

  //! get internal variable manager for reinitialization
  Teuchos::RCP<ScaTraEleInternalVariableManagerLsReinit<my::nsd_,my::nen_> > VarManager(){return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerLsReinit<my::nsd_,my::nen_> >(my::scatravarmanager_);};

  /*========================================================================*/
  //! @name additional members
  /*========================================================================*/

  //! scalar at before reinitialization
  std::vector<LINALG::Matrix<my::nen_,1> > ephizero_;

  //! parameter class for reinitialization
  const ScaTraEleParameterLsReinit* lsreinitparams_;
};


/// Scatra diffusion manager
/*!
  advanced diffusion manager for reinitialization
  - enables crosswind diffusion
  - allows for negative diffusivity
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

    //! set the isotropic diffusion coefficient, which may be negative for elliptic reinitialization
    void SetIsotropicDiff(const double diff, const int k)
    {
      diff_[k] = diff;
      return;
    }

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

      return;
    }

    void Reset()
    {
      for (std::size_t kk=0; kk<diff_.size(); kk++)
      {
        diff_[kk]=0.0;
        sgdiff_[kk]=0.0;
      }
      diffdirectiontensor_.Clear();
      have_cross_wind_diff_ = false;

      return;
    }

    LINALG::Matrix<NSD,NSD> GetCrosswindTensor(){return diffdirectiontensor_;}

    bool HaveCrossWindDiff(){return have_cross_wind_diff_;}

  private:

    //! velocity for crosswind diffusion
    LINALG::Matrix<NSD,NSD> diffdirectiontensor_;

    //! flag for crosswind diffusion
    bool have_cross_wind_diff_;

};


/// ScaTraEleInternalVariableManager implementation
/*!
  advanced form for reinitialization: does not allow for setting all values at once and,
  therefore, provides set functions
*/
template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerLsReinit : public ScaTraEleInternalVariableManager<NSD, NEN>
{
  typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

  public:

  ScaTraEleInternalVariableManagerLsReinit(int numscal):
  ScaTraEleInternalVariableManager<NSD,NEN>(numscal)
  {return;}

  virtual ~ScaTraEleInternalVariableManagerLsReinit(){};

  /** \brief compute and set internal variables
   *
   * \param funct      (in) : array for shape functions
   * \param derxy      (in) : global derivatives of shape functions w.r.t x,y,z
   * \param ephinp     (in) : scalar at t_(n+1) or t_(n+alpha_F)
   * \param ephin      (in) : scalar at t_(n)
   * \param econvelnp  (in) : nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
   * \param ehist      (in) : history vector of transported scalars */
  void SetInternalVariables(
      const LINALG::Matrix<NEN,1>&                 funct,
      const LINALG::Matrix<NSD,NEN>&               derxy,
      const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,
      const std::vector<LINALG::Matrix<NEN,1> >&   ephin,
      const LINALG::Matrix<NSD,NEN>&               econvelnp,
      const std::vector<LINALG::Matrix<NEN,1> >&   ehist
      )
  {
    dserror("Setting all members at once is not allowed for reinitialization!");
    return;
  };

  /*========================================================================*/
  //! @name set methods for internal variables
  /*========================================================================*/
  /* Here, it is explicitly required to set the variables of the manager,
   * since, for instance, gradients my be obtained by using the derivative of
   * the shape functions or by interpolation of nodal gradient values computed
   * based on projections. Likewise the velocity of the reinitialization
   * equation according to Sussman requires special care. */

  //! set current scalar value
  void SetPhinp(const int k,double phinp) {my::phinp_[k] = phinp; return;};
  //! set previous scalar value
  void SetPhin(const int k,double phin) {my::phin_[k] = phin; return;};
  //! set spatial gradient of current scalar value
  void SetGradPhi(const int k,LINALG::Matrix<NSD,1>& gradphi) {my::gradphi_[k] = gradphi; return;};
  //! set convective term of current scalar value
  void SetConvPhi(const int k,double conv_phi) {my::conv_phi_[k] = conv_phi; return;};
  //! set convective velocity
  void SetConVel(LINALG::Matrix<NSD,1>& convel) {my::convelint_ = convel; return;};
  //! set history term of current scalar value
  void SetHist(const int k,double hist) {my::hist_[k] = hist; return;};
  //! set convective part in convective form
  virtual void SetConv(LINALG::Matrix<NEN,1>& conv) {my::conv_ = conv; return;};

  /*========================================================================*/
  //! @name manipulation methods for internal variables
  /*========================================================================*/

  //! set convective term of current scalar value
  void AddToConvPhi(const int k,double conv_phi) {dserror("Currently unused!"); return;};
  //{my::conv_phi_[k] += conv_phi;};
  //! set convective term of current scalar value
  void ScaleConvPhi(const int k,double scale) {dserror("Currently unused!"); return;};
  //{my::conv_phi_[k] *= scale;};

  /*========================================================================*/
  //! @name reset default values
  /*========================================================================*/

  void Reset()
  {
    for (int kk=0; kk<my::numscal_; kk++)
    {
      my::phinp_[kk]=0.0;
      my::phin_[kk]=0.0;
      (my::gradphi_[kk]).Clear();
      my::conv_phi_[kk]=0.0;
      my::hist_[kk]=0.0;
    }

    my::convelint_.Clear();
    my::conv_.Clear();

    return;
  }

};


}
}

namespace SCATRA
{
template<DRT::Element::DiscretizationType CELLDISTYPE>
struct CellTypeToOptGaussRule {};
template<> struct CellTypeToOptGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_36point;};
template<> struct CellTypeToOptGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_37point; };
}

#endif
