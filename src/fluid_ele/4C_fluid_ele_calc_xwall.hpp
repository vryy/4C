/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of fluid element with xfem wall modeling

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_CALC_XWALL_HPP
#define FOUR_C_FLUID_ELE_CALC_XWALL_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidXWall;

    template <Core::FE::CellType distype, Discret::ELEMENTS::Fluid::EnrichmentType enrtype>
    class FluidEleCalcXWall : public FluidEleCalc<distype, enrtype>
    {
      typedef Discret::ELEMENTS::FluidEleCalc<distype, enrtype> my;

      using my::nen_;
      using my::nsd_;
      using my::numderiv2_;

     public:
      /// Singleton access method
      static FluidEleCalcXWall<distype, enrtype>* instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      int evaluate_service(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /// Evaluate supporting methods of the element for xwall
      /*!
        Interface function for supporting methods of the element
       */
      virtual int evaluate_service_x_wall(Discret::ELEMENTS::Fluid* ele,
          Teuchos::ParameterList& params, Teuchos::RCP<Core::Mat::Material>& mat,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3);

      int evaluate(Discret::ELEMENTS::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      void sysmat(const Core::LinAlg::Matrix<nsd_, nen_>& ebofoaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& ebofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& eprescpgn,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveln,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& fsevelaf,
          const Core::LinAlg::Matrix<nen_, 1>& fsescaaf,
          const Core::LinAlg::Matrix<nsd_, nen_>& evel_hat,
          const Core::LinAlg::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf, const Core::LinAlg::Matrix<nen_, 1>& epream,
          const Core::LinAlg::Matrix<nen_, 1>& epren, const Core::LinAlg::Matrix<nen_, 1>& eprenp,
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,
          const Core::LinAlg::Matrix<nen_, 1>& escaaf, const Core::LinAlg::Matrix<nen_, 1>& escaam,
          const Core::LinAlg::Matrix<nen_, 1>& escadtam,
          const Core::LinAlg::Matrix<nsd_, nen_>& eveldtam,
          const Core::LinAlg::Matrix<nen_, 1>& epredtam,
          const Core::LinAlg::Matrix<nen_, 1>& escabofoaf,
          const Core::LinAlg::Matrix<nen_, 1>& escabofon,
          const Core::LinAlg::Matrix<nsd_, nen_>& emhist,
          const Core::LinAlg::Matrix<nsd_, nen_>& edispnp,
          const Core::LinAlg::Matrix<nsd_, nen_>& egridv,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          const Core::LinAlg::Matrix<nen_, 1>& eporo,
          const Core::LinAlg::Matrix<nsd_, 2 * nen_>& egradphi,
          const Core::LinAlg::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const Core::Mat::Material> material, double& Cs_delta_sq,
          double& Ci_delta_sq, double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const Core::FE::GaussIntegration& intpoints) override;

      //! get ALE grid displacements and grid velocity for element
      void get_grid_disp_ale(Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::Matrix<nsd_, nen_>& edispnp) override;

      //! linearisation in the case of mesh motion 3-D
      void lin_mesh_motion_3d(
          Core::LinAlg::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                             ///< pressure at integration point
          const double& timefac,                           ///< time factor
          const double& timefacfac                         ///< = timefac x fac
          ) override;

     private:
      const static int enren_ = Core::FE::num_nodes<distype>;

      /// private Constructor since we are a Singleton.
      FluidEleCalcXWall();

      //! brief evaluate shape functions and their derivatives at integration point
      void eval_shape_func_and_derivs_at_int_point(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
          ) override;

      //! brief evaluate shape functions and their derivatives at integration point
      virtual void eval_std_shape_func_and_derivs_at_int_point(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
      );

      //! prepare special Gauss rule for quadrature
      virtual void prepare_gauss_rule();

      //! get properties of xwall element
      virtual void get_ele_properties(Discret::ELEMENTS::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Teuchos::RCP<Core::Mat::Material>& mat);

      //! brief get enrichment function
      virtual double spaldings_law(double dist, double utau);

      //! calculate derivative of enrichment wrt y+
      virtual double der_spaldings_law(double dist, double utau, double psi);

      //! calculate second derivative of enrichment wrt y+
      virtual double der2_spaldings_law(double dist, double utau, double psi, double derpsi);

      //! evaluate enrichment (1)
      virtual void eval_enrichment();

      //! evaluate enrichment (2)
      virtual double enrichment_shape_der(
          Core::LinAlg::Matrix<nsd_, 1>& derpsigp, Core::LinAlg::Matrix<numderiv2_, 1>& der2psigp);

      //! go increment of tauw for projection matrix
      virtual void x_wall_tau_w_inc_forward();

      //! go increment of tauw backwards for projection matrix
      virtual void x_wall_tau_w_inc_back();

      /*! \brief Calculate wall shear stress via gradient for xwall
       *
       *  \author bk \date 06/2014
       */

      virtual int tau_w_via_gradient(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      //! Get MK
      double get_mk() override;

      /*! \brief Calculate statilization parameter mk entry routine
       *
       *  \author bk \date 06/2014
       */

      virtual int calc_mk(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2);

      /*! \brief Calculate statilization parameter mk
       *
       *  \author bk \date 06/2014
       */
      virtual double calc_mk();

      /*! \brief Calculate Projection on updated shape functions (matrix and rhs)
       *
       *  \author bk \date 06/2014
       */
      virtual int x_wall_projection(Discret::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2);


      //! nodal values of wall distance
      Core::LinAlg::Matrix<enren_, 1> ewdist_;

      //! nodal values of tauw
      Core::LinAlg::Matrix<enren_, 1> etauw_;

      //! nodal values of inctauw
      Core::LinAlg::Matrix<enren_, 1> einctauw_;

      //! nodal values of ramp function
      Core::LinAlg::Matrix<enren_, 1> eramp_;

      //! nodal values of toggle vector
      Core::LinAlg::Matrix<enren_, 1> etoggle_;

      //! nodal values of psi
      Core::LinAlg::Matrix<enren_, 1> epsi_;

      //! nodal values of new psi
      Core::LinAlg::Matrix<enren_, 1> epsinew_;

      //! nodal values of old psi
      Core::LinAlg::Matrix<enren_, 1> epsiold_;

      //! nodal values of inctauw
      Core::LinAlg::Matrix<enren_, 1> eincwdist_;

      //! kinematic viscosity
      double visc_;

      //! inverse of viscosity
      double viscinv_;

      //! density
      double dens_;

      //! inverse of density
      double densinv_;

      //! bool if ramp functions active
      bool is_blending_ele_;

      //! bool if old and new psi should be calculated
      bool calcoldandnewpsi_;
      bool evaluate_sysmat_;

      //! node coordinates
      Core::LinAlg::Matrix<nsd_, enren_> xyze_;

      //! array for enr shape functions
      Core::LinAlg::Matrix<enren_, 1> functenr_;
      //! array for enr shape functions
      Core::LinAlg::Matrix<enren_, 1> funct_;

      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, enren_> derxyenr_;
      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, enren_> derxy_;

      //! global second derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<numderiv2_, enren_> derxyenr2_;

      //! global second derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<numderiv2_, enren_> derxy2_;

      //! array for shape function derivatives w.r.t r,s,t
      Core::LinAlg::Matrix<nsd_, enren_> deriv_;

      //! array for second derivatives of shape function w.r.t r,s,t
      Core::LinAlg::Matrix<numderiv2_, enren_> deriv2_;


      // find a definition of k and B in:
      // R. B. Dean, Reynolds number dependence of skin friction and other bulk flow variables
      // in two-dimensional rectangular duct flow, J. Fluid Eng. 100, 215 (1978)
      // k_ = 1.0/2.44 = 0.409836066
      // B_ = 5.17
      // constants for enrichment function
      const double k_;
      const double b_;
      //! pre-calculated expression exp(-k_*B_)
      const double expmkmb_;

      // element parameter mk
      double mk_;

      int numgpnorm_;
      int numgpnormow_;
      int numgpplane_;

      //! object to construct gauss points in several dimensions
      Teuchos::RCP<Core::FE::CollectedGaussPoints> cgp_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
