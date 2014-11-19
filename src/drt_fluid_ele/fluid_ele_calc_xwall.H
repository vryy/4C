/*!----------------------------------------------------------------------
\file fluid_ele_calc_xwall.H

\brief main file containing routines for calculation of fluid element with xfem wall modeling

<pre>
Maintainer: Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XWALL_H
#define FLUID_ELE_CALC_XWALL_H

#include "fluid_ele_interface.H"
#include "fluid_ele_calc.H"

// includes gauss points
#include "../drt_lib/drt_utils.H"

namespace DRT
{
namespace ELEMENTS
{

  template<DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
  class FluidEleCalcXWall : public FluidEleCalc<distype,enrtype>
  {

    typedef DRT::ELEMENTS::FluidEleCalc<distype,enrtype> my;

public:

  /// dtor
  virtual ~FluidEleCalcXWall() {}

  /// Singleton access method
  static FluidEleCalcXWall<distype,enrtype> * Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  virtual int EvaluateService(
      DRT::ELEMENTS::Fluid*     ele,
      Teuchos::ParameterList&   params,
      Teuchos::RCP<MAT::Material> & mat,
      DRT::Discretization&      discretization,
      std::vector<int>&         lm,
      Epetra_SerialDenseMatrix& elemat1,
      Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1,
      Epetra_SerialDenseVector& elevec2,
      Epetra_SerialDenseVector& elevec3);

  /// Evaluate supporting methods of the element for xwall
  /*!
    Interface function for supporting methods of the element
   */
  virtual int EvaluateServiceXWall(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    Teuchos::RCP<MAT::Material> & mat,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3
    );

  virtual int Evaluate(DRT::ELEMENTS::Fluid*        ele,
      DRT::Discretization &         discretization,
      const std::vector<int> &      lm,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material> & mat,
      Epetra_SerialDenseMatrix&     elemat1_epetra,
      Epetra_SerialDenseMatrix&     elemat2_epetra,
      Epetra_SerialDenseVector&     elevec1_epetra,
      Epetra_SerialDenseVector&     elevec2_epetra,
      Epetra_SerialDenseVector&     elevec3_epetra,
      bool                          offdiag = false);

private:

  const static int enren_=DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  /// private Constructor since we are a Singleton.
  FluidEleCalcXWall();

  //! brief evaluate shape functions and their derivatives at integration point
  virtual void EvalShapeFuncAndDerivsAtIntPoint(const double* gpcoords,      ///< actual integration point (coords)
                                        double gpweight     ///< actual integration point (weight)
                                        );

  //! prepare special Gauss rule for quadrature
  virtual void PrepareGaussRule(Teuchos::ParameterList&    params);

  //! get properties of xwall element
  virtual void GetEleProperties(DRT::ELEMENTS::Fluid*    ele,
                                                            DRT::Discretization & discretization,
                                                            const std::vector<int> & lm,
                                                            Teuchos::ParameterList&    params,
                                                            Teuchos::RCP<MAT::Material> & mat);

  //! brief get enrichment function
  virtual double SpaldingsLaw(double dist, double utau);

  //! calculate derivative of enrichment wrt y+
  virtual double DerSpaldingsLaw(double dist, double utau, double psi);

  //! calculate second derivative of enrichment wrt y+
  virtual double Der2SpaldingsLaw(double dist, double utau, double psi,double derpsi);

  //! evaluate enrichment (1)
  virtual void EvalEnrichment();

  //! evaluate enrichment (2)
  virtual double EnrichmentShapeDer(
      LINALG::Matrix<my::nsd_, 1> &         derpsigp,
      LINALG::Matrix<my::numderiv2_,1> &    der2psigp);

  //! go increment of tauw for projection matrix
  virtual void XWallTauWIncForward();

  //! go increment of tauw backwards for projection matrix
  virtual void XWallTauWIncBack();

  /*! \brief Calculate wall shear stress via gradient for xwall
   *
   *  \author bk \date 06/2014
   */

  virtual int TauWViaGradient(
      DRT::ELEMENTS::Fluid*                ele,
      Teuchos::ParameterList&              params,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseVector&            elevec1,
      Epetra_SerialDenseVector&            elevec2);

  //! Get MK
  virtual double GetMK();

  /*! \brief Calculate statilization parameter mk entry routine
   *
   *  \author bk \date 06/2014
   */

  virtual int CalcMK(
      DRT::ELEMENTS::Fluid*                ele,
      Teuchos::ParameterList&              params,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseVector&            elevec1,
      Epetra_SerialDenseVector&            elevec2);

  /*! \brief Calculate statilization parameter mk
   *
   *  \author bk \date 06/2014
   */
  virtual double CalcMK();

  /*! \brief Calculate Projection on updated shape functions (matrix and rhs)
   *
   *  \author bk \date 06/2014
   */
  virtual int XWallProjection(
      DRT::ELEMENTS::Fluid*                ele,
      Teuchos::ParameterList&              params,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseMatrix&            elemat1,
      Epetra_SerialDenseMatrix&            elemat2);

  /*! \brief Calculate Projection on updated shape functions (matrix and rhs) including continuity constraint
   *
   *  \author bk \date 07/2014
   */
  virtual int XWallProjectionWithContinuityConstraint(
      DRT::ELEMENTS::Fluid*                ele,
      Teuchos::ParameterList&              params,
      DRT::Discretization &                discretization,
      const std::vector<int> &             lm,
      Teuchos::RCP<MAT::Material> &        mat,
      Epetra_SerialDenseMatrix&            elemat1,
      Epetra_SerialDenseMatrix&            elemat2,
      Epetra_SerialDenseVector&            elevec1);

  //! nodal values of wall distance
  LINALG::Matrix<enren_,1>    ewdist_;

  //! nodal values of tauw
  LINALG::Matrix<enren_,1>    etauw_;

  //! nodal values of inctauw
  LINALG::Matrix<enren_,1>    einctauw_;

  //! nodal values of ramp function
  LINALG::Matrix<enren_,1>    eramp_;

  //! nodal values of toggle vector
  LINALG::Matrix<enren_,1>    etoggle_;

  //! nodal values of psi
  LINALG::Matrix<enren_,1>    epsi_;

  //! nodal values of new psi
  LINALG::Matrix<enren_,1>    epsinew_;

  //! nodal values of old psi
  LINALG::Matrix<enren_,1>    epsiold_;

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

  //! node coordinates
  LINALG::Matrix<my::nsd_,enren_> xyze_;

  //! array for enr shape functions
  LINALG::Matrix<enren_,1> functenr_;
  //! array for enr shape functions
  LINALG::Matrix<enren_,1> funct_;

  //! global derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<my::nsd_,enren_> derxyenr_;
  //! global derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<my::nsd_,enren_> derxy_;

  //! global second derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<my::numderiv2_,enren_> derxyenr2_;

  //! global second derivatives of shape functions w.r.t x,y,z
  LINALG::Matrix<my::numderiv2_,enren_> derxy2_;

  //! array for shape function derivatives w.r.t r,s,t
  LINALG::Matrix<my::nsd_,enren_> deriv_;

  //! array for second derivatives of shape function w.r.t r,s,t
  LINALG::Matrix<my::numderiv2_,enren_> deriv2_;


  // find a definition of k and B in:
  // R. B. Dean, Reynolds number dependence of skin friction and other bulk flow variables
  // in two-dimensional rectangular duct flow, J. Fluid Eng. 100, 215 (1978)
  // k_ = 1.0/2.44 = 0.409836066
  // B_ = 5.17
  // constants for enrichment function
  const double k_;
  const double B_;
  //! pre-calculated expression exp(-k_*B_)
  const double expmkmb_;

  // element parameter mk
  double mk_;

  //! object to construct gauss points in several dimensions
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp_;
  };
  }
}

#endif
