/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.H

\brief Internal implementation of XFluid element interface coupling

<pre>
Maintainer: Raffaela Kruse /Benedikt Schott
            kruse@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_H
#define FLUID_ELE_CALC_XFEM_H

#include "../drt_lib/drt_utils.H"

#include "fluid_ele_calc.H"
#include "fluid_ele_calc_xfem_coupling.H"

namespace DRT
{
namespace ELEMENTS
{
class FluidEleParameterXFEM;

/// Fluid element interface coupling implementation with XFEM
/*!
  This internal class keeps all the working arrays needed to
  calculate the interface stabilization for XFEM with fluid elements.
  The method ElementXfemInterfaceHybridLM() provides a clean and fast element implementation
  for Mixed/Stress/Hybrid interface coupling, using either Cauchy stress-based or viscous
  stress-based Lagrange multipliers.
  The method ElementXfemInterfaceNIT() provides a clean and fast element implementation
  for interface coupling using Nitsche's method (only the cut element sided mortaring
  for xfluid and xfluidfluid applications).
  The method ElementXfemInterfaceNIT2() provides a clean and fast element implementation
  for interface coupling using Nitsche's method also for non-"cut element sided mortaring"
  (can be used for embedded mortaring in case of Xfluidfluid).

  <h3>Purpose</h3>

  The fluid element will allocate exactly one object of this class for all
  fluid elements with the same number of nodes in the mesh. This
  allows us to use exactly matching working arrays (and keep them
  around.)

  The code is meant to be as clean as possible. This is the only way
  to keep it fast. The number of working arrays has to be reduced to
  a minimum so that the element fits into the cache. (There might be
  room for improvements.)

  <h3>Usability</h3>

  The calculations are done by the EvaluateXfemInterface...() methods.

  \author schott
  \date 04/12
 */
template<DRT::Element::DiscretizationType distype>
class FluidEleCalcXFEM : public FluidEleCalc<distype>
{

  /// private Constructor since we are a Singleton.
  FluidEleCalcXFEM();

  typedef FluidEleCalc<distype> my;

private:

  /// pointer to the cast object, fluid parameter list for XFEM
  DRT::ELEMENTS::FluidEleParameterXFEM* fldparaxfem_;


  /// number of stress-dof
  static const int numstressdof_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

  //! @name useful constants for DOF-index numbering
  //@{
  static const unsigned Velx = 0;
  static const unsigned Vely = 1;
  static const unsigned Velz = 2;
  static const unsigned Pres = 3;

  static const unsigned Sigmaxx = 0;
  static const unsigned Sigmaxy = 1;
  static const unsigned Sigmaxz = 2;
  static const unsigned Sigmayx = 1;
  static const unsigned Sigmayy = 3;
  static const unsigned Sigmayz = 4;
  static const unsigned Sigmazx = 2;
  static const unsigned Sigmazy = 4;
  static const unsigned Sigmazz = 5;
  //@}

  /// get stress dof-index
  unsigned stressIndex(unsigned xi, unsigned xj){ return (xi * xj) ? xi + xj + 1 : xi + xj;  }

public:
  virtual ~FluidEleCalcXFEM() {}

  /// Singleton access method
  static FluidEleCalcXFEM<distype> * Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /// evaluate the XFEM cut element
  virtual int EvaluateXFEM(
      DRT::ELEMENTS::Fluid*                             ele,
      DRT::Discretization &                             discretization,
      const std::vector<int> &                          lm,
      Teuchos::ParameterList&                           params,
      Teuchos::RCP<MAT::Material> &                     mat,
      Epetra_SerialDenseMatrix&                         elemat1_epetra,
      Epetra_SerialDenseMatrix&                         elemat2_epetra,
      Epetra_SerialDenseVector&                         elevec1_epetra,
      Epetra_SerialDenseVector&                         elevec2_epetra,
      Epetra_SerialDenseVector&                         elevec3_epetra,
      const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
      const GEO::CUT::plain_volumecell_set &            cells,
      bool                                              offdiag = false);

  /// evaluate the shape functions in the XFEM
  virtual int IntegrateShapeFunctionXFEM(
      DRT::ELEMENTS::Fluid*                             ele,
      DRT::Discretization &                             discretization,
      const std::vector<int> &                          lm,
      Epetra_SerialDenseVector&                         elevec1_epetra,
      const std::vector<DRT::UTILS::GaussIntegration> & intpoints,
      const GEO::CUT::plain_volumecell_set &            cells);

  /// error computation
  virtual int ComputeError(DRT::ELEMENTS::Fluid*        ele,
      Teuchos::ParameterList&                           params,
      Teuchos::RCP<MAT::Material>&                      mat,
      DRT::Discretization&                              discretization,
      std::vector<int>&                                 lm,
      Epetra_SerialDenseVector&                         ele_dom_norms);

  virtual int ComputeError(DRT::ELEMENTS::Fluid*        ele,
      Teuchos::ParameterList&                           params,
      Teuchos::RCP<MAT::Material>&                      mat,
      DRT::Discretization&                              discretization,
      std::vector<int>&                                 lm,
      Epetra_SerialDenseVector&                         ele_dom_norms,
      const DRT::UTILS::GaussIntegration &              intpoints);

  virtual int ComputeErrorInterface(
      DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
      DRT::Discretization &                                               dis,               ///< background discretization
      const std::vector<int> &                                            lm,                ///< element local map
      Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
      Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
      DRT::Discretization &                                               cutdis,            ///< cut discretization
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
      Teuchos::ParameterList&                                             params,            ///< parameter list
      const GEO::CUT::plain_volumecell_set&                               vcSet              ///< volumecell sets in this element
  );

  virtual int ComputeErrorInterfaceXFluidFluid(
      DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
      DRT::Discretization &                                               dis,               ///< background discretization
      const std::vector<int> &                                            lm,                ///< element local map
      Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
      Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
      DRT::Discretization &                                               cutdis,            ///< cut discretization
      DRT::Discretization &                                               embdis,            ///< embedded discretization
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
      Teuchos::ParameterList&                                             params,            ///< parameter list
      const GEO::CUT::plain_volumecell_set&                               vcSet,             ///< volumecell sets in this element
      std::map<int,int> &                                                 boundary_emb_gid_map
  );


  /// add terms from mixed/hybrid Lagrange multiplier coupling approach to element matrix and rhs
  virtual void ElementXfemInterfaceHybridLM(
      DRT::ELEMENTS::Fluid *                                            ele,                      ///< fluid element
      DRT::Discretization &                                             dis,                      ///< background discretization
      const std::vector<int> &                                          lm,                       ///< element local map
      const std::vector<DRT::UTILS::GaussIntegration> &                 intpoints,                ///< element gauss points
      DRT::Discretization &                                             cutdis,                   ///< boundary discretization
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> >&       bcells,                   ///< boundary cells
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & bintpoints,               ///< boundary integration points
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > &           side_coupling,            ///< side coupling matrices
      Teuchos::ParameterList&                                           params,                   ///< parameter list
      Epetra_SerialDenseMatrix&                                         elemat1_epetra,           ///< local system matrix of intersected element
      Epetra_SerialDenseVector&                                         elevec1_epetra,           ///< local element vector of intersected element
      Epetra_SerialDenseMatrix&                                         Cuiui,                    ///< coupling matrix of a side with itself
      const GEO::CUT::plain_volumecell_set &                            vcSet,                    ///< set of plain volume cells
      bool                                                              fluidfluidcoupling = false///< indicates fluid-fluid coupling context
  );


  /// add Nitsche (NIT) interface condition to element matrix and rhs
  virtual void ElementXfemInterfaceNIT(
      DRT::ELEMENTS::Fluid *                                              ele,                ///< fluid element
      DRT::Discretization &                                               dis,               ///< background discretization
      const std::vector<int> &                                            lm,                ///< element local map
      DRT::Discretization &                                               cutdis,            ///< cut discretization
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
      Teuchos::ParameterList&                                             params,            ///< parameter list
      Epetra_SerialDenseMatrix&                                           elemat1_epetra,    ///< element matrix
      Epetra_SerialDenseVector&                                           elevec1_epetra,    ///< element vector
      Epetra_SerialDenseMatrix&                                           Cuiui,             ///< ui-ui coupling matrix
      const GEO::CUT::plain_volumecell_set&                               vcSet,             ///< volumecell sets in this element
      bool                                                                fluidfluidcoupling = false ///< indicates fluid-fluid coupling context
  );

  /// add two-sided Nitsche (NIT2) interface condition to element matrix and rhs
  virtual void ElementXfemInterfaceNIT2(
      DRT::ELEMENTS::Fluid *                                              ele,
      DRT::Discretization &                                               dis,
      const std::vector<int> &                                            lm,
      DRT::Discretization &                                               cutdis,
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,
      Teuchos::ParameterList&                                             params,
      DRT::Discretization &                                               alediscret,
      std::map<int,int> &                                                 boundary_emb_gid_map,
      Epetra_SerialDenseMatrix&                                           elemat1_epetra,
      Epetra_SerialDenseVector&                                           elevec1_epetra,
      Epetra_SerialDenseMatrix&                                           Cuiui,
      const GEO::CUT::plain_volumecell_set&                               vcSet
  );

  ///
  virtual void CalculateContinuityXFEM(
      DRT::ELEMENTS::Fluid *               ele,            ///< fluid element
      DRT::Discretization &                dis,            ///< discretization
      const std::vector<int> &             lm,             ///< local map
      Epetra_SerialDenseVector&            elevec1_epetra, ///< element vector
      const DRT::UTILS::GaussIntegration & intpoints       ///< integration points
  );

  ///
  virtual void CalculateContinuityXFEM(
      DRT::ELEMENTS::Fluid *      ele,                ///< fluid element
      DRT::Discretization &       dis,                ///< discretization
      const std::vector<int> &    lm,                 ///< local map
      Epetra_SerialDenseVector&   elevec1_epetra      ///< element vector
  );

private:

  //! evaluate analytical reference solution
  void AnalyticalReference(
      const int                               calcerr,     ///< which reference solution
      LINALG::Matrix<my::nsd_,1> &            u,           ///< exact jump vector (coupled)
      LINALG::Matrix<my::nsd_,my::nsd_> &     grad_u,      ///< exact velocity gradient
      double &                                p,           ///< exact pressure
      LINALG::Matrix<my::nsd_,1> &            xyzint,      ///< xyz position of gaussian point
      const double &                          t,           ///< time
      Teuchos::RCP<MAT::Material>             mat = Teuchos::null
  );

  //! build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
  void NIT_BuildPatchCuiui(
      Epetra_SerialDenseMatrix &                              Cuiui,            ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
      std::map<int, std::vector<Epetra_SerialDenseMatrix> >&  Cuiui_coupling    ///< Cuiui matrices for all cutting sides
  );

  //! compute viscous part of Nitsche's penalty term scaling for Nitsche's method
  double NIT_ComputeNitscheStabfac(
      const DRT::Element::DiscretizationType ele_distype, ///< the discretization type of the element w.r.t which the stabilization factor is computed
      const double h_k,                                   ///< the characteristic element length
      const int sid = -1                                  ///< Id of the side element for embedded sided fluidfluid formulations
      );

  //! get the constant which satisfies the trace inequality depending on the spatial dimension and polynomial order of the element
  double NIT_getTraceEstimateConstant(const DRT::Element::DiscretizationType ele_distype);

  //! compute full Nitsche's penalty term scaling for Nitsche's method (xfluid/xfsi)
  void NIT_ComputeFullStabfac(
      double &                         NIT_full_stab_fac,     ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
      const double                     NIT_visc_stab_fac,     ///< Nitsche's viscous scaling part of penalty term
      INPAR::XFEM::ConvStabScaling     conv_stab_scaling,     ///< type of convective stabilization for xfluid
      const double                     veln_normal            ///< interface-normal velocity contribution
  );

  //! compute stabilization factors for Nitsche's method (fluid-fluid coupling)
  void NIT_ComputeFullStabfacFluidFluid(
      double &                         NIT_full_stab_fac,             ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
      double &                         avg_conv_stab_fac,             ///< to be filled: scaling of the convective average coupling term for fluidfluid problems
      double &                         velgrad_interface_fac,         ///< to be filled: stabilization parameter for velocity-gradients penalty term
      double &                         press_coupling_fac,            ///< to be filled: stabilization parameter for pressure-coupling penalty term
      INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling,         ///< type of convective stabilization for fluid-fluid problem
      const double                     veln_normal,                   ///< interface-normal velocity contribution
      const LINALG::Matrix<my::nsd_,1> velint,                        ///< interface velocity
      const double                     h_k,                           ///< characteristic element length
      const double                     NIT_visc_stab_fac,             ///< Nitsche's viscous scaling part of penalty term
      const double                     gamma_ghost_penalty = 0.0,     ///< ghost-penalty parameter
      const double                     gamma_press_coupling = 0.0     ///< factor of stabilization parameter for pressure-coupling penalty term
  );

  //! create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices
  void PatchLocationVector(
      std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
      DRT::Discretization &                                   cutdis,                  ///< cut discretization
      std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
      std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling           ///< coupling matrices
  );

  //! create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices (for embedded coupling)
  void PatchLocationVector(
      std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
      DRT::Discretization &                                   cutdis,                  ///< cut discretization
      std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
      std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
      std::map<int,int> &                                     boundary_emb_gid_map     ///< map between boundary sid and corresponding embedded element id
  );

  //! prepare coupling matrices, that include contributions from convective stabilization
  void HybridLM_CreateConvStabMatrices(
      std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
      DRT::Discretization &                                   cutdis,                  ///< cut discretization
      std::map<int, std::vector<Epetra_SerialDenseMatrix> > & conv_side_coupling       ///< contributions to coupling matrices from convective stabilizations
  );

  //! build traction vector w.r.t fluid domain
  void BuildTractionVector(
      LINALG::Matrix<my::nsd_,1> &  traction,   ///< traction vector
      double &                      press,      ///< pressure at gaussian point
      LINALG::Matrix<my::nsd_,1> &  normal      ///< normal vector
  );

  //! assemble side's interface force
  void AssembleInterfaceForce(
      Teuchos::RCP<Epetra_Vector>            iforcecol, ///< interface force column vector
      DRT::Discretization &                  cutdis,    ///< cut discretization
      std::vector<int> &                     lm,        ///< local dof map
      Epetra_SerialDenseVector &             iforce     ///< interface force vector
  );

  //! compute transformation factor for surface integration, normal, local and global gp coordinates
  void ComputeSurfaceTransformation(
      double &                    drs,         ///< surface transformation factor
      LINALG::Matrix<3,1> &       x_gp_lin,    ///< global coordiantes of gaussian point
      LINALG::Matrix<3,1> &       normal,      ///< normal vector on boundary cell
      GEO::CUT::BoundaryCell *    bc,          ///< boundary cell
      const LINALG::Matrix<2,1> & eta          ///< local coordinates of gaussian point w.r.t boundarycell
  );

  //! evaluate shape function and derivative at point with local coordinates rst
  void EvalFuncAndDeriv( LINALG::Matrix<3,1> &  rst );

  //! build matrices from volume-based terms for Cauchy & viscous stress-based mixed/hybrid LM-coupling
  //! \author kruse \date 06/14
  void HybridLM_Build_VolBased(
      const INPAR::XFEM::Hybrid_LM_L2_Proj                                                        hybrid_lm_l2_proj,
      const std::vector<DRT::UTILS::GaussIntegration> &                                           intpoints,
      const GEO::CUT::plain_volumecell_set&                                                       cells,
      const LINALG::Matrix<my::nsd_, my::nen_> &                                                  evelaf,            ///< element velocity
      const LINALG::Matrix<my::nen_,1>&                                                           epreaf,            ///< element pressure
      LINALG::Matrix<my::nen_,my::nen_>&                                                          bK_ss,             ///< block K_ss matrix
      LINALG::Matrix<my::nen_,my::nen_>&                                                          invbK_ss,          ///< inverse of block K_ss matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> &   K_su,              ///< K_su matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                           rhs_s,             ///< rhs_s vector
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> &   K_us,              ///< K_us matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                  K_uu,              ///< K_uu matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1> &                               rhs_uu,            ///< rhs_u(u) vector
      const bool                                                                                  is_MHVS = false,   ///< viscous (true) or Cauchy (false) stress-based LM
      const double                                                                                mhvs_param = 1.0   ///< stabilizing parameter for viscous stress-based LM
  );

  //! evaluate matrices from volume-based terms for Cauchy stress-based mixed/hybrid LM coupling at current Gauss-point
  void MHCS_Evaluate_VolBased(
      const LINALG::Matrix<my::nsd_,my::nen_>&                                                  evelaf,   ///< element velocity
      const LINALG::Matrix<my::nen_,1>&                                                         epreaf,   ///< element pressure
      LINALG::Matrix<my::nen_,my::nen_>&                                                        bK_ss,    ///< block K_ss matrix
      LINALG::Matrix<my::nen_,my::nen_>&                                                        invbK_ss, ///< inverse of block K_ss matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_>&  K_su,     ///< K_su matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s     ///< rhs_s vector
  );

  //! evaluate matrices from volume-based terms for viscous stress-based mixed/hybrid LM coupling at current Gauss-point
  //! \author kruse \date 06/14
  void MHVS_Evaluate_VolBased(
      const LINALG::Matrix<my::nsd_,my::nen_> &                                                 evelaf,   ///< element velocity
      LINALG::Matrix<my::nen_,my::nen_> &                                                       bK_ss,    ///< block K_ss matrix
      LINALG::Matrix<my::nen_,my::nen_> &                                                       invbK_ss, ///< inverse of block K_ss matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> & K_su,     ///< K_su matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s,    ///< rhs_s vector
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> & K_us,     ///< K_us matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                K_uu,     ///< K_uu matrix
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_, 1> &                             rhs_uu,   ///< rhs_uu vector
      const double &                                                                            mhvs_param///< stabilizing parameter
  );

  //! evaluate matrices from surface-based terms for Cauchy & viscous stress-based mixed/hybrid LM coupling at current Gauss-point
  //! \author kruse \date 06/14
  void HybridLM_Evaluate_SurfBased(
      Teuchos::RCP<DRT::ELEMENTS::XFLUID::HybridLMInterface<distype> > &                        si,
      const LINALG::Matrix<my::nen_,my::nen_> &                                                 bK_ss,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,numstressdof_,my::numdofpernode_> & K_su,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::numdofpernode_,numstressdof_> & K_us,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,numstressdof_,1> &                         rhs_s,
      const LINALG::Matrix<my::nen_,1> &                                                        epreaf,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,my::nsd_> &                K_uu,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_,1> &                              rhs_uu,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,my::nsd_,1> &                       G_up,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,1,my::nsd_> &                       G_pu,
      LINALG::BlockMatrix<LINALG::Matrix<my::nen_,1>,my::nsd_,1> &                              rhs_up,
      LINALG::Matrix<my::nen_,1> &                                                              rhs_pu,
      const LINALG::Matrix<3,1> &                                                               normal,
      const double &                                                                            timesurffac,
      const bool                                                                                eval_side_coupling,
      const bool                                                                                is_MHVS
  );

  //! compute characteristic element length h_k
  //! \author kruse \date 06/14
  double ComputeCharEleLength(
      DRT::Element *                                                        ele,                   ///< fluid element
      Epetra_SerialDenseMatrix &                                            ele_xyze,              ///< element coordinates
      const GEO::CUT::plain_volumecell_set &                                vcSet,                 ///< volumecell sets for volume integration
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &          bcells,                ///< bcells for boundary cell integration
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &     bintpoints,            ///< integration points for boundary cell integration
      Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<distype> >  emb = Teuchos::null,   ///< pointer to the embedded coupling implementation
      DRT::Element *                                                        face = NULL            ///< side element in 3D
  );

  //! compute the measure of the elements surface with given local id
  double ComputeMeasFace(
      DRT::Element *            ele,              ///< fluid element
      Epetra_SerialDenseMatrix& ele_xyze,         ///< element coordinates
      const int                 local_face_id     ///< the local id of the face w.r.t the fluid element
  );

  //! pre-compute the measure of the element's intersecting surface
  double ComputeMeasCutSurf(
      const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary cell integration points
      const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells             ///< boundary cells
  );

  //! compute volume-equivalent diameter
  inline double ComputeVolEqDiameter(double vol)
  {
    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

    return hk;
  }

};


} //namespace ELEMENTS
} //namespace DRT

#endif

