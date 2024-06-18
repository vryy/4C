/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for reinitialization equation

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_LSREINIT_HPP
#define FOUR_C_SCATRA_ELE_CALC_LSREINIT_HPP

#include "4C_config.hpp"

#include "4C_fem_geometry_geo_utils.hpp"
#include "4C_fem_geometry_integrationcell.hpp"
#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    template <int NSD>
    class ScaTraEleDiffManagerLsReinit;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerLsReinit;
    class ScaTraEleParameterLsReinit;

    template <Core::FE::CellType distype, unsigned probDim>
    class ScaTraEleCalcLsReinit : public ScaTraEleCalc<distype, probDim>
    {
     private:
      //! private constructor for singletons
      ScaTraEleCalcLsReinit(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype, probDim> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      //! Singleton access method
      static ScaTraEleCalcLsReinit<distype, probDim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /*========================================================================*/
      //! @name access routines
      /*========================================================================*/

      //! setup element evaluation
      int SetupCalc(
          Core::Elements::Element* ele, Core::FE::Discretization& discretization) override;

      //! evaluate the element
      int evaluate(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

     protected:
      //! calculate matrix and rhs. Here the whole thing is hidden. Hyperbolic reinit.
      void sysmat_hyperbolic(
          Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs   //!< element rhs to calculate
      );

      //! calculate matrix and rhs. Here the whole thing is hidden. Ellipitic reinit.
      void sysmat_elliptic(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs,                   //!< element rhs to calculate
          const Core::Geo::BoundaryIntCellPtrs& bcell              //!< interface for penalty term
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
      void elliptic_newton_system(Core::LinAlg::SerialDenseMatrix* emat,
          Core::LinAlg::SerialDenseVector* erhs,
          const Core::LinAlg::Matrix<nen_, 1>& el2sysmat_diag_inv,
          const Core::Geo::BoundaryIntCellPtrs& bcell);

     private:
      /*========================================================================*/
      //! @name general evaluation methods
      /*========================================================================*/

      void eval_reinitialization(const Epetra_Vector& phinp, const std::vector<int>& lm,
          Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra);

      void eval_reinitialization_embedded(const std::vector<int>& lm, Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra);

      void eval_reinitialization_std(const Epetra_Vector& phinp, const std::vector<int>& lm,
          Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra);

      /*========================================================================*/
      //! @name overloaded methods for evaluation of individual terms
      /*========================================================================*/

      //! calculation of diffusive element matrix
      void calc_mat_diff(Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                           //!< index of current scalar
          const double timefacfac  //!< domain-integration factor times time-integration factor
          ) override;

      //! standard Galerkin diffusive term on right hand side
      void calc_rhs_diff(Core::LinAlg::SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                           //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const Core::LinAlg::Matrix<nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
      );

      /*========================================================================*/
      //! @name additional methods required for reinitialization
      /*========================================================================*/

      //! sign function
      void sign_function(double& sign_phi,                   //!< sign of phi
          const double charelelength,                        //!< characteristic element length
          const double phizero,                              //!< initial phi
          const Core::LinAlg::Matrix<nsd_, 1>& gradphizero,  //!< gradient of initial phi
          const double phi,                                  //!< phi at time n+1
          const Core::LinAlg::Matrix<nsd_, 1>& gradphi       //!< gradient of phi at time n+1
      );

      //! derivative of sign function
      void deriv_sign_function(double& deriv_sign,  //!< sign of phi
          const double charelelength,               //!< characteristic element length
          const double phizero                      //!< initial phi
      );

      //! calculation of characteristic element length, i.e., interface thickness
      double calc_char_ele_length_reinit(const double vol,  //!< element volume
          const Core::LinAlg::Matrix<nsd_, 1>& gradphizero  //!< gradient of initial phi
      );

      /*========================================================================*/
      //! @name penalty methods for reinitialization and related
      /*========================================================================*/

      //! calculation of element-wise denominator of penalty parameter
      void calc_ele_penalty_parameter(double& penalty  //!< penalty parameter
      );

      //! calculate system matrix and rhs for correction step
      void sysmat_correction(const double penalty,  ///!< element penalty parameter
          Core::LinAlg::SerialDenseMatrix& emat,    ///!< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs     ///!< element rhs to calculate
      );

      //! calculation of penalty term on rhs
      void calc_rhs_penalty(Core::LinAlg::SerialDenseVector& erhs,  //!< rhs vector
          const double fac,                                         //!< domain integration factor
          const double penalty,                                     //!< penalty parameter
          const double deriv_sign,                                  //!< derivative of sign function
          const double norm_gradphizero  //!< norm of gradient of initial phi
      );

      //! calculation of interface penalty term for elliptic reinitialization
      void evaluate_interface_term(
          Core::LinAlg::SerialDenseMatrix* emat,       //!< element matrix to calculate
          Core::LinAlg::SerialDenseVector* erhs,       //!< element vector to calculate
          const Core::Geo::BoundaryIntCellPtrs& bcell  //!< interface for penalty term
      );

      //! calculation of interface penalty term for elliptic reinitialization (gauss loop)
      template <Core::FE::CellType celldistype>
      void calc_penalty_term(
          Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs,  //!< element vector to calculate
          const Core::Geo::BoundaryIntCell& cell  //!< interface cell
      );


      /** \brief calculation of interface penalty term for elliptic reinitialization
       *         (special variant for the 0-D case)
       *
       *  \param emat (out) : element matrix to calculate
       *  \param erhs (out) : element vector to calculate
       *  \param cell (in)  : interface boundary integration cell
       *
       *  \author hiermeier \date 11/16 */
      void calc_penalty_term_0_d(Core::LinAlg::SerialDenseMatrix* emat,
          Core::LinAlg::SerialDenseVector* erhs, const Core::Geo::BoundaryIntCell& cell);

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! calculate system matrix and rhs for velocity projection
      void sysmat_nodal_vel(const int dir,        ///< current spatial direction
          Core::LinAlg::SerialDenseMatrix& emat,  ///< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs   ///< element rhs to calculate
      );

      //! get diffusion manager for reinitialization
      Teuchos::RCP<ScaTraEleDiffManagerLsReinit<nsd_>> diff_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerLsReinit<nsd_>>(my::diffmanager_);
      };

      //! get internal variable manager for reinitialization
      Teuchos::RCP<ScaTraEleInternalVariableManagerLsReinit<nsd_, nen_>> var_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerLsReinit<nsd_, nen_>>(
            my::scatravarmanager_);
      };

      /*========================================================================*/
      //! @name additional members
      /*========================================================================*/

      //! scalar at before reinitialization
      std::vector<Core::LinAlg::Matrix<nen_, 1>> ephizero_;

      //! parameter class for reinitialization
      const ScaTraEleParameterLsReinit* lsreinitparams_;
    };


    /// Scatra diffusion manager
    /*!
      advanced diffusion manager for reinitialization
      - enables crosswind diffusion
      - allows for negative diffusivity
    */
    template <int NSD>
    class ScaTraEleDiffManagerLsReinit : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerLsReinit(int numscal)
          : ScaTraEleDiffManager(numscal), diffdirectiontensor_(true), have_cross_wind_diff_(false)
      {
        return;
      }

      //! set the isotropic diffusion coefficient, which may be negative for elliptic
      //! reinitialization
      void SetIsotropicDiff(const double diff, const int k)
      {
        diff_[k] = diff;
        return;
      }

      void set_velocity_for_cross_wind_diff(const Core::LinAlg::Matrix<NSD, 1> velocity)
      {
        if (NSD != 3) FOUR_C_THROW("Currently only 3d problems supported for crosswind diffusion");

        // compute tensor for anisotropic artificial diffusion
        // i.e., crosswind diffusion

        // get norm of velocity
        const double vel_norm_sq = velocity.Norm2() * velocity.Norm2();

        // compute tensor
        if (vel_norm_sq > 1.0e-8)
        {
          diffdirectiontensor_(0, 0) = 1.0 - velocity(0, 0) * velocity(0, 0) / vel_norm_sq;
          diffdirectiontensor_(0, 1) = -velocity(0, 0) * velocity(1, 0) / vel_norm_sq;
          diffdirectiontensor_(0, 2) = -velocity(0, 0) * velocity(2, 0) / vel_norm_sq;
          diffdirectiontensor_(1, 0) = diffdirectiontensor_(0, 1);
          diffdirectiontensor_(1, 1) = 1.0 - velocity(1, 0) * velocity(1, 0) / vel_norm_sq;
          diffdirectiontensor_(1, 2) = -velocity(1, 0) * velocity(2, 0) / vel_norm_sq;
          diffdirectiontensor_(2, 0) = diffdirectiontensor_(0, 2);
          diffdirectiontensor_(2, 1) = diffdirectiontensor_(1, 2);
          diffdirectiontensor_(2, 2) = 1.0 - velocity(2, 0) * velocity(2, 0) / vel_norm_sq;
        }
        else
          diffdirectiontensor_.Clear();

        // indicate that crosswind diffusion has to be used
        have_cross_wind_diff_ = true;

        return;
      }

      void Reset()
      {
        for (std::size_t kk = 0; kk < diff_.size(); kk++)
        {
          diff_[kk] = 0.0;
          sgdiff_[kk] = 0.0;
        }
        diffdirectiontensor_.Clear();
        have_cross_wind_diff_ = false;

        return;
      }

      Core::LinAlg::Matrix<NSD, NSD> GetCrosswindTensor() { return diffdirectiontensor_; }

      bool HaveCrossWindDiff() { return have_cross_wind_diff_; }

     private:
      //! velocity for crosswind diffusion
      Core::LinAlg::Matrix<NSD, NSD> diffdirectiontensor_;

      //! flag for crosswind diffusion
      bool have_cross_wind_diff_;
    };


    /// ScaTraEleInternalVariableManager implementation
    /*!
      advanced form for reinitialization: does not allow for setting all values at once and,
      therefore, provides set functions
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerLsReinit
        : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerLsReinit(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal)
      {
        return;
      }

      /** \brief compute and set internal variables
       *
       * \param funct      (in) : array for shape functions
       * \param derxy      (in) : global derivatives of shape functions w.r.t x,y,z
       * \param ephinp     (in) : scalar at t_(n+1) or t_(n+alpha_F)
       * \param ephin      (in) : scalar at t_(n)
       * \param econvelnp  (in) : nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
       * \param ehist      (in) : history vector of transported scalars */
      void set_internal_variables(const Core::LinAlg::Matrix<NEN, 1>& funct,
          const Core::LinAlg::Matrix<NSD, NEN>& derxy,
          const std::vector<Core::LinAlg::Matrix<NEN, 1>>& ephinp,
          const std::vector<Core::LinAlg::Matrix<NEN, 1>>& ephin,
          const Core::LinAlg::Matrix<NSD, NEN>& econvelnp,
          const std::vector<Core::LinAlg::Matrix<NEN, 1>>& ehist)
      {
        FOUR_C_THROW("Setting all members at once is not allowed for reinitialization!");
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
      void SetPhinp(const int k, double phinp)
      {
        my::phinp_[k] = phinp;
        return;
      };
      //! set previous scalar value
      void SetPhin(const int k, double phin)
      {
        my::phin_[k] = phin;
        return;
      };
      //! set spatial gradient of current scalar value
      void SetGradPhi(const int k, Core::LinAlg::Matrix<NSD, 1>& gradphi) override
      {
        my::gradphi_[k] = gradphi;
        return;
      };
      //! set convective term of current scalar value
      void SetConvPhi(const int k, double conv_phi) override
      {
        my::conv_phi_[k] = conv_phi;
        return;
      };
      //! set convective velocity
      void SetConVel(const int k, Core::LinAlg::Matrix<NSD, 1>& convel)
      {
        my::convelint_[k] = convel;
      };
      //! set history term of current scalar value
      void SetHist(const int k, double hist)
      {
        my::hist_[k] = hist;
        return;
      };
      //! set convective part in convective form
      virtual void SetConv(const int k, Core::LinAlg::Matrix<NEN, 1>& conv)
      {
        my::conv_[k] = conv;
      };

      /*========================================================================*/
      //! @name manipulation methods for internal variables
      /*========================================================================*/

      //! set convective term of current scalar value
      void AddToConvPhi(const int k, double conv_phi) override
      {
        FOUR_C_THROW("Currently unused!");
        return;
      };
      //{my::conv_phi_[k] += conv_phi;};
      //! set convective term of current scalar value
      void ScaleConvPhi(const int k, double scale) override
      {
        FOUR_C_THROW("Currently unused!");
        return;
      };
      //{my::conv_phi_[k] *= scale;};

      /*========================================================================*/
      //! @name reset default values
      /*========================================================================*/

      void Reset()
      {
        for (int kk = 0; kk < my::numscal_; kk++)
        {
          my::phinp_[kk] = 0.0;
          my::phin_[kk] = 0.0;
          (my::gradphi_[kk]).Clear();
          my::conv_phi_[kk] = 0.0;
          my::hist_[kk] = 0.0;
          my::convelint_[kk].Clear();
          my::conv_[kk].Clear();
        }



        return;
      }
    };


  }  // namespace ELEMENTS
}  // namespace Discret

namespace ScaTra
{
  template <Core::FE::CellType CELLDISTYPE>
  struct CellTypeToOptGaussRule
  {
  };
  template <>
  struct CellTypeToOptGaussRule<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_36point;
  };
  template <>
  struct CellTypeToOptGaussRule<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_37point;
  };
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
