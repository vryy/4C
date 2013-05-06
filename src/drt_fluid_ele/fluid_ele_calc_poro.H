#ifndef FLUID_ELE_CALC_PORO_H
#define FLUID_ELE_CALC_PORO_H

#include "fluid_ele_calc.H"

// forward declarations
namespace MAT
{
  class StructPoro;
}


namespace DRT
{
namespace ELEMENTS
{
  class So_Poro_Interface;

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcPoro : public FluidEleCalc<distype>
{
    typedef FluidEleCalc<distype> my;

protected:
    /// private Constructor since we are a Singleton.
    FluidEleCalcPoro();

public:
    virtual ~FluidEleCalcPoro() {}

    /// Singleton access method
    static FluidEleCalcPoro<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();

    /*!
    \brief calculate element matrix and rhs for porous flow

    \param eid              (i) element id
    \param discretization   (i) fluid discretization the element belongs to
    \param lm               (i) location matrix of element
    \param mat              (i) material
    \param elemat1_epetra   (o) element matrix to calculate
    \param elemat2_epetra   (o) element matrix to calculate
    \param elevec1_epetra   (o) element vector to calculate
    \param elevec2_epetra   (o) element vector to calculate
    \param elevec3_epetra   (o) element vector to calculate
    \param offdiag          (i) flag indicating wether diagonal or off diagonal blocks are to be calculated

    */
    virtual int Evaluate(DRT::ELEMENTS::Fluid*         ele,
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

    /// Evaluate the element at specified gauss points for porous flow
    virtual int Evaluate(DRT::ELEMENTS::Fluid*                    ele,
                             DRT::Discretization &                discretization,
                             const std::vector<int> &             lm,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material> &        mat,
                             Epetra_SerialDenseMatrix&            elemat1_epetra,
                             Epetra_SerialDenseMatrix&            elemat2_epetra,
                             Epetra_SerialDenseVector&            elevec1_epetra,
                             Epetra_SerialDenseVector&            elevec2_epetra,
                             Epetra_SerialDenseVector&            elevec3_epetra,
                             const DRT::UTILS::GaussIntegration & intpoints);

    /// Evaluate the element at specified gauss points for porous flow
    virtual int EvaluateOD(DRT::ELEMENTS::Fluid*                  ele,
                             DRT::Discretization &                discretization,
                             const std::vector<int> &             lm,
                             Teuchos::ParameterList&              params,
                             Teuchos::RCP<MAT::Material> &        mat,
                             Epetra_SerialDenseMatrix&            elemat1_epetra,
                             Epetra_SerialDenseMatrix&            elemat2_epetra,
                             Epetra_SerialDenseVector&            elevec1_epetra,
                             Epetra_SerialDenseVector&            elevec2_epetra,
                             Epetra_SerialDenseVector&            elevec3_epetra,
                             const DRT::UTILS::GaussIntegration & intpoints);


protected:

    /*!
      \brief evaluate function for Fluid element for porous flow

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.
     */
    int Evaluate( Teuchos::ParameterList&                                         params,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       ebofoaf,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat1,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,1> &                       elevec1,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelaf,
                  const LINALG::Matrix<my::nen_,1>    &                           epreaf,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelnp,
                  const LINALG::Matrix<my::nen_,1>    &                           eprenp,
                  const LINALG::Matrix<my::nen_,1>    &                           epren,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
                  const LINALG::Matrix<my::nen_,1>&                               echist,
                  const LINALG::Matrix<my::nen_,1>    &                           epressnp_timederiv,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispnp,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispn,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       egridv,
                  const LINALG::Matrix<my::nen_,1>&                               escaaf,
                  const LINALG::Matrix<my::nen_,1>*                               eporositynp,
                  const LINALG::Matrix<my::nen_,1>*                               eporositydot,
                  const LINALG::Matrix<my::nen_,1>*                               eporositydotn,
                  Teuchos::RCP<MAT::Material>                                     mat,
                  bool                                                            isale,
                  const DRT::UTILS::GaussIntegration &                            intpoints);

    /*!
      \brief evaluate function for Fluid element for porous flow

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.
     */
    int EvaluateOD(
        Teuchos::ParameterList&                                           params,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> &  elemat1,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
        const LINALG::Matrix<my::nen_, 1> &                               epreaf,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
        const LINALG::Matrix<my::nen_, 1> &                               eprenp,
        const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
        const LINALG::Matrix<my::nen_,1>&                                 escaaf,
        Teuchos::RCP<MAT::Material>                                       mat,
        bool                                                              isale,
        const DRT::UTILS::GaussIntegration &                              intpoints);

    /*!
    \brief calculate element matrix and rhs for porous flow

    \param params           (i) element parameter list
    \param ebofoaf          (i) body force at n+alpha_F/n+1
    \param evelaf           (i) nodal velocities at n+alpha_F/n+1
    \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf           (i) nodal pressure at n+alpha_F/n+1
    \param eprenp           (i) nodal pressure at n+alpha_F/n+1
    \param eaccam           (i) nodal accelerations at n+alpha_M
    \param emhist           (i) time rhs for momentum equation
    \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
    \param edispn           (i) nodal displacements at n (on moving mesh)
    \param egridv           (i) grid velocity at n+1
    \param estif            (o) element matrix to calculate
    \param emesh            (o) linearization wrt mesh motion
    \param eforce           (o) element rhs to calculate
    \param material         (i) fluid material
    \param isale            (i) ALE flag
    \param intpoints        (i) Gaussian integration points

    */
  void Sysmat(
      Teuchos::ParameterList&                                       params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      ebofoaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      evelaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      evelnp,
      const LINALG::Matrix<my::nen_,1>&                             epreaf,
      const LINALG::Matrix<my::nen_,1>&                             eprenp,
      const LINALG::Matrix<my::nen_,1>&                             epren,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      eaccam,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      emhist,
      const LINALG::Matrix<my::nen_,1>&                             echist,
      const LINALG::Matrix<my::nen_,1>    &                         epressnp_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispn,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      egridv,
      const LINALG::Matrix<my::nen_,1>&                             escaaf,
      const LINALG::Matrix<my::nen_,1>*                             eporositynp,
      const LINALG::Matrix<my::nen_,1>*                             eporositydot,
      const LINALG::Matrix<my::nen_,1>*                             eporositydotn,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&                      eforce,
      Teuchos::RCP<const MAT::Material>                             material,
      bool                                                          isale,
      const DRT::UTILS::GaussIntegration &                          intpoints);

  void SysmatOD(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_,my::nsd_ * my::nen_>&  ecoupl,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
      Teuchos::RCP<const MAT::Material>                               material,
      bool                                                            isale,
      const DRT::UTILS::GaussIntegration &                            intpoints);


  //! linearisation in the case of mesh motion 3-D for Poroelasticity
  void LinMeshMotion_3D_OD(
      LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>&         ecoupl_u,
      LINALG::Matrix< my::nen_, my::nsd_ * my::nen_>&                   ecoupl_p,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
      const LINALG::Matrix<my::nen_, 1>&                                epreaf,
      const double &                                                    press,
      const double &                                                    press_dot,
      const double &                                                    porosity,
      const LINALG::Matrix<my::nsd_, 1>&                                gradphi,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    J,
      const double &                                                    refporositydot,
      const double & timefac, const double &                            timefacfac
      );

  //! linearisation in the case of mesh motion 2-D for Poroelasticity
  void LinMeshMotion_2D_OD(
      LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>&         ecoupl_u,
      LINALG::Matrix< my::nen_, my::nsd_ * my::nen_>&                   ecoupl_p,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
      const LINALG::Matrix<my::nen_, 1>&                                epreaf,
      const double &                                                    press,
      const double &                                                    press_dot,
      const double &                                                    porosity,
      const LINALG::Matrix<my::nsd_, 1>&                                gradphi,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    J,
      const double &                                                    refporositydot,
      const double & timefac, const double &                            timefacfac
      );

  //! Compute element matrix entries: PSPG
  void PSPG(
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &  estif_q_u,      ///< block (weighting function q x u)
      LINALG::Matrix<my::nen_,my::nen_> &            ppmat,          ///< block (weighting function q x p)
      LINALG::Matrix<my::nen_,1> &                   preforce,       ///< rhs forces pressure
      LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &   lin_resM_Du,    ///< linearisation of the stabilization residual
      const double &                                 fac3,           ///< factor for residual in current subgrid velocities
      const double &                                 timefacfac,     ///< = timefac x fac
      const double &                                 timefacfacpre, ///< = timefacpre x fac
      const double &                                 rhsfac       ///< right-hand-side factor for residuals
      );

  void ComputeFDerivative(const LINALG::Matrix<my::nsd_, my::nen_>&         edispnp,
                          const LINALG::Matrix<my::nsd_,my::nsd_>&          defgrd_inv,
                                LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>& F_x,
                                LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>& F_X);

  void ComputeContiTimeRHS(Teuchos::ParameterList&                   params,
                           const LINALG::Matrix<my::nen_,1>&         echist,
                           const LINALG::Matrix<my::nen_, 1>&        epren,
                           const LINALG::Matrix<my::nsd_, my::nen_>& edispn,
                           LINALG::Matrix<my::nen_,1>&               preforce,
                           const double&                             rhsfac,
                           const double&                             scale
                          );

  void ComputeGradients(const double&                                      J,
                        const double&                                      dphidp,
                        const double&                                      dphidJ,
                        const LINALG::Matrix<my::nsd_*my::nsd_,1>&         defgrd_IT_vec,
                        const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_x,
                        const LINALG::Matrix<my::nen_,1>*                  eporositynp,
                        LINALG::Matrix<my::nsd_,1>&                        gradJ,
                        LINALG::Matrix<my::nsd_,1>&                        grad_porosity);

  virtual void ComputePorosityGradient(
                        const double&                                      dphidp,
                        const double&                                      dphidJ,
                        LINALG::Matrix<my::nsd_,1>&                        gradJ,
                        const LINALG::Matrix<my::nen_,1>*                  eporositynp,
                        LINALG::Matrix<my::nsd_,1>&                        grad_porosity);

  void ComputeLinearization(
                              const double&                                      dphi_dp,
                              const double&                                      dphi_dpp,
                              const double&                                      dphi_dJp,
                              const LINALG::Matrix<my::nsd_,1>&                  gradJ,
                              LINALG::Matrix<my::nsd_,my::nen_>&                 dgradphi_dp);

  void ComputeLinearizationOD(const double&                                      J,
                              const double&                                      dphi_dJ,
                              const double&                                      dphi_dJJ,
                              const double&                                      dphi_dJp,
                              const LINALG::Matrix<my::nsd_,my::nsd_>&           defgrd_inv,
                              const LINALG::Matrix<my::nsd_*my::nsd_,1>&         defgrd_IT_vec,
                              const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_x,
                              const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_X,
                              const LINALG::Matrix<my::nsd_,1>&                  gradJ,
                              LINALG::Matrix<1,my::nsd_*my::nen_>&               dJ_dus,
                              LINALG::Matrix<1,my::nsd_*my::nen_>&               dphi_dus,
                              LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&        dgradphi_dus);

  virtual void ComputePorosity( Teuchos::ParameterList&             params,
                             const double&                          press,
                             const double&                          J,
                             const int&                             gp,
                             const LINALG::Matrix<my::nen_,1>&      shapfct,
                             const LINALG::Matrix<my::nen_,1>*      myporosity,
                             double&                                porosity,
                             double*                                dphi_dp,
                             double*                                dphi_dJ,
                             double*                                dphi_dJdp,
                             double*                                dphi_dJJ,
                             double*                                dphi_dpp,
                             bool                                   save);

  void GaussPointLoopOD(
                          Teuchos::ParameterList&                                         params,
                          const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
                          const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
                          const LINALG::Matrix<my::nen_, 1>&                              epreaf,
                          const LINALG::Matrix<my::nen_, 1>&                              eprenp,
                          const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
                          const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
                          const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
                          const LINALG::Matrix<my::nen_,1>&                               escaaf,
                          const LINALG::Matrix<my::nen_,1>*                               eporositynp,
                          LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
                          LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&       ecoupl_u,
                          LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&                  ecoupl_p,
                          Teuchos::RCP<const MAT::Material>                               material,
                          const DRT::UTILS::GaussIntegration &                            intpoints);

  void FillMatrixMomentumOD(
      const double&                                               timefacfac,
      const double&                                               porosity,
      const LINALG::Matrix<my::nsd_,1>&                           gridvelint,
      const LINALG::Matrix<my::nsd_,1>&                           grad_porosity,
      const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
      LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&   ecoupl_u);

  void FillMatrixContiOD(
      const double&                                               timefacfacpre,
      const double&                                               porosity,
      const double&                                               J,
      const double&                                               dphi_dJ,
      const double&                                               dphi_dJJ,
      const double&                                               dphi_dJdp,
      const double&                                               press_dot,
      const double&                                               gridvdiv,
      const LINALG::Matrix<my::nsd_,1>&                           gridvelint,
      const LINALG::Matrix<my::nsd_,1>&                           grad_porosity,
      const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dJ_dus,
      LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&              ecoupl_p);


  void PreEvaluate
    (
      Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::ELEMENTS::Fluid*                 ele, ///< fluid element
      DRT::Discretization&      discretization  ///< pointer to discretization for de-assembly
    );

  /// computation of material derivatives
  double SetupMaterialDerivatives();

  /// access structure material of corresponding solid (poro) element
  void GetStructMaterial();

  ///effective viscosity (zero for Darcy flow )
  double visceff_;

  ///first derivatives of shape functions w.r.t. material coordinates
  LINALG::Matrix<my::nsd_,my::nen_> N_XYZ_;
  /*!
  \brief second derivatives of shape functions w.r.t. material coordinates XYZ
    3D: (N,XX ; N,YZ ; N,ZZ ; N,XY ; N,XZ ; N,YX ; N,YZ;  N,ZX ; N,ZY)
    2D: (N,XX ; N,YZ ; N,XY ; N,YX )
   */
  LINALG::Matrix<my::numderiv2_,my::nen_> N_XYZ2_;
  ///second derivatives of shape functions w.r.t. material coordinates XYZ (ordered in symmetric matrix)
  LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> N_XYZ2full_;
  /// material coordinates
  LINALG::Matrix<my::nsd_,my::nen_> xyze0_;

  //! vector containing all values from previous timelevel n for continuity equation
  double histcon_;

  //! material of corresponding solid element
  //Teuchos::RCP<MAT::StructPoro> structmat_;

  DRT::ELEMENTS::So_Poro_Interface* so_interface_;

};
}
}

#endif

