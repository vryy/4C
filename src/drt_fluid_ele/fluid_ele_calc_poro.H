#ifndef FLUID_ELE_CALC_PORO_H
#define FLUID_ELE_CALC_PORO_H

#include "fluid_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcPoro : public FluidEleCalc<distype>
{
    /// private Constructor since we are a Singleton.
    FluidEleCalcPoro();

    typedef FluidEleCalc<distype> my;

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


private:

    /*!
      \brief evaluate function for Fluid element for porous flow

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.
     */
    int Evaluate( Teuchos::ParameterList&                                         params,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       ebofoaf,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat1,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat2,
                  LINALG::Matrix<(my::nsd_+1)*my::nen_,1> &                       elevec1,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelaf,
                  const LINALG::Matrix<my::nen_,1>    &                           epreaf,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       evelnp,
                  const LINALG::Matrix<my::nen_,1>    &                           eprenp,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
                  const LINALG::Matrix<my::nen_,1>    &                           epren,
                  const LINALG::Matrix<my::nen_,1>    &                           epressn_timederiv,
                  const LINALG::Matrix<my::nen_,1>    &                           epressnp_timederiv,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispnp,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       edispn,
                  const LINALG::Matrix<my::nsd_,my::nen_> &                       egridv,
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
        const LINALG::Matrix<my::nsd_, my::nen_> &                        ebofoaf,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> &  elemat1,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
        const LINALG::Matrix<my::nen_, 1> &                               epreaf,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
        const LINALG::Matrix<my::nen_, 1> &                               eprenp,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        emhist,
        const LINALG::Matrix<my::nen_, 1> &                               epren,
        const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        eaccam,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        edispn,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
        Teuchos::RCP<MAT::Material>                                       mat,
        bool                                                              isale,
        const DRT::UTILS::GaussIntegration &                              intpoints);

    /*!
    \brief calculate element matrix and rhs for porous flow

    \param ebofoaf          (i) body force at n+alpha_F/n+1
    \param evelaf           (i) nodal velocities at n+alpha_F/n+1
    \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf           (i) nodal pressure at n+alpha_F/n+1
    \param eprenp           (i) nodal pressure at n+alpha_F/n+1
    \param eaccam           (i) nodal accelerations at n+alpha_M
    \param emhist           (i) time rhs for momentum equation
    \param edispnp          (i) nodal displacements (on moving mesh)
    \param egridv           (i) grid velocity (on moving mesh)
    \param estif            (o) element matrix to calculate
    \param emesh            (o) linearization wrt mesh motion
    \param eforce           (o) element rhs to calculate
    \param material         (i) fluid material
    \param isale            (i) ALE flag
    \param intpoints        (i) Gaussian integration points

    */
  void Sysmat(
      const LINALG::Matrix<my::nsd_,my::nen_>&                      ebofoaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      evelaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      evelnp,
      const LINALG::Matrix<my::nen_,1>&                             epreaf,
      const LINALG::Matrix<my::nen_,1>&                             eprenp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      eaccam,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      emhist,
      const LINALG::Matrix<my::nen_,1>    &                         epren,
      const LINALG::Matrix<my::nen_,1>    &                         epressn_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                         epressnp_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispn,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      egridv,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  emesh,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&                      eforce,
      Teuchos::RCP<const MAT::Material>                             material,
      bool                                                          isale,
      const DRT::UTILS::GaussIntegration &                          intpoints);

  void SysmatOD(
      const LINALG::Matrix<my::nsd_, my::nen_>&                       ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       emhist,
      const LINALG::Matrix<my::nen_, 1> &                             epren,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_,my::nsd_ * my::nen_>&  ecoupl,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
      Teuchos::RCP<const MAT::Material>                               material,
      bool                                                            isale,
      const DRT::UTILS::GaussIntegration &                            intpoints);


  //! linearisation in the case of mesh motion 3-D for Poroelasticity
  void LinMeshMotion_3D_OD(
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_) * my::nen_>& emesh,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
      const double &                                                    press,
      const double &                                                    press_dot,
      const double &                                                    porosity,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    J,
      LINALG::Matrix<1, my::nsd_>&                                      gradJ,
      const double & timefac, const double &                            timefacfac
      );

  //! linearisation in the case of mesh motion 2-D for Poroelasticity
  void LinMeshMotion_2D_OD(
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_) * my::nen_>& emesh,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                         egridv,
      const double &                                                    press,
      const double &                                                    press_dot,
      const double &                                                    porosity,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    J,
      LINALG::Matrix<1, my::nsd_>&                                      gradJ,
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
      const double &                                 timefacfacpre,
      const double &                                 rhsfac       ///< right-hand-side factor for residuals
      );

  void CalcAuxiliaryDerivatives(const LINALG::Matrix<my::nsd_, my::nen_>& edispnp,
                                const LINALG::Matrix<my::nsd_,my::nsd_>& defgrd,
                                const LINALG::Matrix<my::nsd_,my::nsd_>& defgrd_inv,
                                LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>& N_X_x,
                                LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>& F_x,
                                LINALG::Matrix<my::nsd_,my::nsd_*my::nen_>& Finv_N_X_x);

  double visceff_;

};
}
}

#endif

