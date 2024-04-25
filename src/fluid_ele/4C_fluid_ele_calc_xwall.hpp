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

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidXWall;

    template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
    class FluidEleCalcXWall : public FluidEleCalc<distype, enrtype>
    {
      typedef DRT::ELEMENTS::FluidEleCalc<distype, enrtype> my;

      using my::nen_;
      using my::nsd_;
      using my::numderiv2_;

     public:
      /// Singleton access method
      static FluidEleCalcXWall<distype, enrtype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      int EvaluateService(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      /// Evaluate supporting methods of the element for xwall
      /*!
        Interface function for supporting methods of the element
       */
      virtual int EvaluateServiceXWall(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3);

      int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      void Sysmat(const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& ebofon,
          const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgn,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelam,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveln,
          const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& fsevelaf,
          const CORE::LINALG::Matrix<nen_, 1>& fsescaaf,
          const CORE::LINALG::Matrix<nsd_, nen_>& evel_hat,
          const CORE::LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
          const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& epream,
          const CORE::LINALG::Matrix<nen_, 1>& epren, const CORE::LINALG::Matrix<nen_, 1>& eprenp,
          const CORE::LINALG::Matrix<nsd_, nen_>& eaccam,
          const CORE::LINALG::Matrix<nen_, 1>& escaaf, const CORE::LINALG::Matrix<nen_, 1>& escaam,
          const CORE::LINALG::Matrix<nen_, 1>& escadtam,
          const CORE::LINALG::Matrix<nsd_, nen_>& eveldtam,
          const CORE::LINALG::Matrix<nen_, 1>& epredtam,
          const CORE::LINALG::Matrix<nen_, 1>& escabofoaf,
          const CORE::LINALG::Matrix<nen_, 1>& escabofon,
          const CORE::LINALG::Matrix<nsd_, nen_>& emhist,
          const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
          const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce,
          const CORE::LINALG::Matrix<nen_, 1>& eporo,
          const CORE::LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
          const CORE::LINALG::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
          const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
          Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
          double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
          const CORE::FE::GaussIntegration& intpoints) override;

      //! get ALE grid displacements and grid velocity for element
      void GetGridDispALE(DRT::Discretization& discretization, const std::vector<int>& lm,
          CORE::LINALG::Matrix<nsd_, nen_>& edispnp) override;

      //! linearisation in the case of mesh motion 3-D
      void LinMeshMotion_3D(
          CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,  ///< mesh motion
          const CORE::LINALG::Matrix<nsd_, nen_>& evelaf,  ///< velocity at time n+alpha_f / n+1
          const double& press,                             ///< pressure at integration point
          const double& timefac,                           ///< time factor
          const double& timefacfac                         ///< = timefac x fac
          ) override;

     private:
      const static int enren_ = CORE::FE::num_nodes<distype>;

      /// private Constructor since we are a Singleton.
      FluidEleCalcXWall();

      //! brief evaluate shape functions and their derivatives at integration point
      void EvalShapeFuncAndDerivsAtIntPoint(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
          ) override;

      //! brief evaluate shape functions and their derivatives at integration point
      virtual void EvalStdShapeFuncAndDerivsAtIntPoint(
          const double* gpcoords,  ///< actual integration point (coords)
          double gpweight          ///< actual integration point (weight)
      );

      //! prepare special Gauss rule for quadrature
      virtual void PrepareGaussRule();

      //! get properties of xwall element
      virtual void GetEleProperties(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat);

      //! brief get enrichment function
      virtual double SpaldingsLaw(double dist, double utau);

      //! calculate derivative of enrichment wrt y+
      virtual double DerSpaldingsLaw(double dist, double utau, double psi);

      //! calculate second derivative of enrichment wrt y+
      virtual double Der2SpaldingsLaw(double dist, double utau, double psi, double derpsi);

      //! evaluate enrichment (1)
      virtual void EvalEnrichment();

      //! evaluate enrichment (2)
      virtual double EnrichmentShapeDer(
          CORE::LINALG::Matrix<nsd_, 1>& derpsigp, CORE::LINALG::Matrix<numderiv2_, 1>& der2psigp);

      //! go increment of tauw for projection matrix
      virtual void XWallTauWIncForward();

      //! go increment of tauw backwards for projection matrix
      virtual void XWallTauWIncBack();

      /*! \brief Calculate wall shear stress via gradient for xwall
       *
       *  \author bk \date 06/2014
       */

      virtual int TauWViaGradient(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2);

      //! Get MK
      double GetMK() override;

      /*! \brief Calculate statilization parameter mk entry routine
       *
       *  \author bk \date 06/2014
       */

      virtual int CalcMK(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2);

      /*! \brief Calculate statilization parameter mk
       *
       *  \author bk \date 06/2014
       */
      virtual double CalcMK();

      /*! \brief Calculate Projection on updated shape functions (matrix and rhs)
       *
       *  \author bk \date 06/2014
       */
      virtual int XWallProjection(DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2);


      //! nodal values of wall distance
      CORE::LINALG::Matrix<enren_, 1> ewdist_;

      //! nodal values of tauw
      CORE::LINALG::Matrix<enren_, 1> etauw_;

      //! nodal values of inctauw
      CORE::LINALG::Matrix<enren_, 1> einctauw_;

      //! nodal values of ramp function
      CORE::LINALG::Matrix<enren_, 1> eramp_;

      //! nodal values of toggle vector
      CORE::LINALG::Matrix<enren_, 1> etoggle_;

      //! nodal values of psi
      CORE::LINALG::Matrix<enren_, 1> epsi_;

      //! nodal values of new psi
      CORE::LINALG::Matrix<enren_, 1> epsinew_;

      //! nodal values of old psi
      CORE::LINALG::Matrix<enren_, 1> epsiold_;

      //! nodal values of inctauw
      CORE::LINALG::Matrix<enren_, 1> eincwdist_;

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
      CORE::LINALG::Matrix<nsd_, enren_> xyze_;

      //! array for enr shape functions
      CORE::LINALG::Matrix<enren_, 1> functenr_;
      //! array for enr shape functions
      CORE::LINALG::Matrix<enren_, 1> funct_;

      //! global derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<nsd_, enren_> derxyenr_;
      //! global derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<nsd_, enren_> derxy_;

      //! global second derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<numderiv2_, enren_> derxyenr2_;

      //! global second derivatives of shape functions w.r.t x,y,z
      CORE::LINALG::Matrix<numderiv2_, enren_> derxy2_;

      //! array for shape function derivatives w.r.t r,s,t
      CORE::LINALG::Matrix<nsd_, enren_> deriv_;

      //! array for second derivatives of shape function w.r.t r,s,t
      CORE::LINALG::Matrix<numderiv2_, enren_> deriv2_;


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
      Teuchos::RCP<CORE::FE::CollectedGaussPoints> cgp_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
