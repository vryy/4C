/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem.H

\brief Internal implementation of XFluid element interface coupling

<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_XFEM_H
#define FLUID_ELE_CALC_XFEM_H

#include "../drt_inpar/inpar_xfem.H"

#include "fluid_ele_calc.H"
#include "fluid_ele_calc_xfem_coupling.H"

namespace DRT
{
namespace ELEMENTS
{

/// Fluid element interface coupling implementation with XFEM
/*!
  This internal class keeps all the working arrays needed to
  calculate the interface stabilization for XFEM with fluid elements.
  The method ElementXfemInterfaceMSH() provides a clean and fast element implementation
  for Mixed/Stress/Hybrid interface coupling.
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

public:
    virtual ~FluidEleCalcXFEM() {}

    /// Singleton access method
    static FluidEleCalcXFEM<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();

    /// Evaluate the XFEM cut element
    virtual int EvaluateXFEM(DRT::ELEMENTS::Fluid*                             ele,
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
                             std::string&                                      VCellGaussPts,
                             const GEO::CUT::plain_volumecell_set &            cells,
                             bool                                              offdiag = false);

    /// error computation
    virtual int ComputeError(DRT::ELEMENTS::Fluid*         ele,
                             Teuchos::ParameterList&       params,
                             Teuchos::RCP<MAT::Material>&  mat,
                             DRT::Discretization&          discretization,
                             std::vector<int>&             lm,
                             Epetra_SerialDenseVector&     ele_dom_norms);

    virtual int ComputeError(DRT::ELEMENTS::Fluid*                ele,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material>&         mat,
                             DRT::Discretization&                 discretization,
                             std::vector<int>&                    lm,
                             Epetra_SerialDenseVector&            ele_dom_norms,
                             const DRT::UTILS::GaussIntegration & intpoints);

    virtual int ComputeErrorInterface(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
        Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
        Teuchos::ParameterList&                                             params,            ///< parameter list
        const GEO::CUT::plain_volumecell_set&                               vcSet              ///< volumecell sets in this element
    );

    virtual int ComputeErrorInterfacefluidfluidcoupling(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        Teuchos::RCP<MAT::Material>&                                        mat,               ///< material
        Epetra_SerialDenseVector&                                           ele_interf_norms,  /// squared element interface norms
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        DRT::Discretization &                                               embdis,            ///< embedded discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
        Teuchos::ParameterList&                                             params,            ///< parameter list
        const GEO::CUT::plain_volumecell_set&                               vcSet,             ///< volumecell sets in this element
        std::map<int,int> &                                                 boundary_emb_gid_map
      );

    /// add mixed/stress/hybrid (MSH) interface condition to element matrix and rhs
    virtual void ElementXfemInterfaceMSH(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
        DRT::Discretization &                                               dis,               ///< background discretization
        const std::vector<int> &                                            lm,                ///< element local map
        const std::vector<DRT::UTILS::GaussIntegration> &                   intpoints,         ///< background element integration points
        DRT::Discretization &                                               cutdis,            ///< cut discretization
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells,            ///< boundary cells
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary integration points
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > &             side_coupling,     ///< side coupling matrices
        Teuchos::ParameterList&                                             params,            ///< parameter list
        Epetra_SerialDenseMatrix&                                           elemat1_epetra,    ///< element matrix
        Epetra_SerialDenseVector&                                           elevec1_epetra,    ///< element vector
        Epetra_SerialDenseMatrix&                                           Cuiui,             ///< ui-ui coupling matrix
        std::string&                                                        VCellGaussPts,     ///< Method of volumecell gauss point generation
        const GEO::CUT::plain_volumecell_set&                               vcSet,              ///< Volumecells present in the considered set
        bool                                                                fluidfluidcoupling ///< Is this xfluidfluid problem?
    );

    /// add Nitsche (NIT) interface condition to element matrix and rhs
    virtual void ElementXfemInterfaceNIT(
        DRT::ELEMENTS::Fluid *                                              ele,               ///< fluid element
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
        const GEO::CUT::plain_volumecell_set&                               vcSet,						 ///< volumecell sets in this element
        bool                                                                fluidfluidcoupling ///< Is this xfluidfluid problem?
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
    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid *               ele,            ///< fluid element
                                         DRT::Discretization &                dis,            ///< discretization
                                         const std::vector<int> &             lm,             ///< local map
                                         Epetra_SerialDenseVector&            elevec1_epetra, ///< element vector
                                         const DRT::UTILS::GaussIntegration & intpoints       ///< integration points
      );

    ///
    virtual void CalculateContinuityXFEM(DRT::ELEMENTS::Fluid *      ele,                ///< fluid element
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
        const double &                          t            ///< time
        );

    //! build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
    void NIT_BuildPatchCuiui(
        Epetra_SerialDenseMatrix &                              Cuiui,            ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
        std::map<int, std::vector<Epetra_SerialDenseMatrix> >&  Cuiui_coupling    ///< Cuiui matrices for all cutting sides
    );

    //! compute stabilization factors for Nitsche's method
    void NIT_ComputeStabfac(
        bool                            fluidfluidcoupling,
        double &                        stabfac_visc,         ///< Nitsche stabilization factor
        double &                        stabfac_conv,         ///< stabilization factor 2
        const double                    NIT_stab_fac,
        INPAR::XFEM::ConvStabScaling    conv_stab_scaling,
        const double                    conv_stab_fac,
        const double                    veln_normal,
        const double                    gamma_ghost_penalty,
        double &                        velgrad_interface_fac,
        const double                    gamma_press_coupling,
        double &                        press_coupling_fac,
        const double                    h_k
    );

    //! pre-compute the measure of the element's intersecting surface
    double ComputeMeasSurf(
        const std::map<int, std::vector<DRT::UTILS::GaussIntegration> > &   bintpoints,        ///< boundary cell integration points
        const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > &        bcells             ///< boundary cells
    );

    //! create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices
    void PatchLocationVector(
        std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
        DRT::Discretization &                                   cutdis,                  ///< cut discretization
        std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
        std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
        string                                                  coupl_method             ///< coupling method
    );

    //! create location vector w.r.t patch of intersecting boundary elements and reshape coupling matrices (for embedded coupling)
    void PatchLocationVector(
        std::set<int> &                                         begids,                  ///< ids of intersecting boundary elements
        DRT::Discretization &                                   cutdis,                  ///< cut discretization
        std::vector<int> &                                      patchelementslmv,        ///< lm vector for patch of boundary elements
        std::vector<int> &                                      patchelementslmowner,    ///< lmowner vector for patch of boundary elements
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > & Cuiui_coupling,          ///< coupling matrices
        std::map<int,int> &                                     boundary_emb_gid_map,    ///< map between boundary sid and corresponding embedded element id
        string                                                  coupl_method             ///< coupling method
    );

    //! build traction vector
    void buildTractionVector(
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

    //! build volume-based K matrices for coupling stress fields
    void MSH_Build_K_Matrices(
        INPAR::XFEM::MSH_L2_Proj                                               msh_l2_proj,       ///< full or partial l2 projection for MSH method
        const std::vector<DRT::UTILS::GaussIntegration> &                      intpoints,         ///< background element integration points
        std::string &                                                          VCellGaussPts,     ///< volumecell gaussian points method
        const GEO::CUT::plain_volumecell_set&                                  cells,             ///< volumecells in the present set
        LINALG::Matrix<my::nsd_,my::nen_>&                                     evelaf,            ///< element velocity
        LINALG::Matrix<my::nen_,1>&                                            epreaf,            ///< element pressure
        LINALG::Matrix<my::nen_,my::nen_>&                                     bK_ss,             ///< block K_ss matrix
        LINALG::Matrix<my::nen_,my::nen_>&                                     invbK_ss,          ///< inverse of block K_ss matrix
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)>& K_su,              ///< K_su matrix
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>&                rhs                ///< rhs vector
        );

    //! evaluate volume-based K matrix terms at Gaussian point
    void MSH_EvaluateMatrices(
        LINALG::Matrix<my::nsd_,my::nen_>&    evelaf,                                ///< element velocity
        LINALG::Matrix<my::nen_,1>&           epreaf,                                ///< element pressure
        LINALG::Matrix<my::nen_,my::nen_>&    bK_ss,                                 ///< block K_ss matrix
        LINALG::Matrix<my::nen_,my::nen_>&    invbK_ss,                              ///< inverse of block K_ss matrix
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,my::nen_>,6,(my::nsd_+1)>& K_su, ///< K_su matrix
        LINALG::BlockMatrix<LINALG::Matrix<my::nen_,   1>,6,1>&                rhs   ///< rhs vector
        );

//    Teuchos::RCP<DRT::ELEMENTS::XFLUID::SideInterface<distype> > CreateSideImpl(
//        bool fluidfluidcoupling,
//        DRT::Element * side,
//        Epetra_SerialDenseMatrix & side_xyze,
//        std::vector<Epetra_SerialDenseMatrix> & side_matrices,
//        std::map<int, std::vector<Epetra_SerialDenseMatrix> > Cuiui_coupling
//    );

    inline double HK(double vol)
    {
      // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
      const double hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

      return hk;
    }


};


} //namespace ELEMENTS
} //namespace DRT

#endif

