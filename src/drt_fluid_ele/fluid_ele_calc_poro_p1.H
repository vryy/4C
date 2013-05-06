/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro_p1.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef FLUID_ELE_CALC_PORO_P1_H_
#define FLUID_ELE_CALC_PORO_P1_H_

#include "fluid_ele_calc_poro.H"

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
class FluidEleCalcPoroP1 : public FluidEleCalcPoro<distype>
{
public:
    /// private Constructor since we are a Singleton.
    FluidEleCalcPoroP1();

    typedef FluidEleCalcPoro<distype> my;

    virtual ~FluidEleCalcPoroP1() {}

    /// Singleton access method
    static FluidEleCalcPoroP1<distype> * Instance( bool create=true );

    /// called upon destruction
    virtual void Done();

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
                  const LINALG::Matrix<my::nen_,1>&                               eporositynp,
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
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, (my::nsd_ + 1) * my::nen_> &  elemat1,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
        const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
        const LINALG::Matrix<my::nen_, 1> &                               epreaf,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
        const LINALG::Matrix<my::nen_, 1> &                               eprenp,
        const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
        const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
        const LINALG::Matrix<my::nen_,1>&                                 escaaf,
        const LINALG::Matrix<my::nen_,1>*                                 eporositynp,
        Teuchos::RCP<MAT::Material>                                       mat,
        bool                                                              isale,
        const DRT::UTILS::GaussIntegration &                              intpoints);

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
        const LINALG::Matrix<my::nen_,1>*                               eporositynp,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_,(my::nsd_ + 1) * my::nen_>&  ecoupl,
        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
        Teuchos::RCP<const MAT::Material>                               material,
        bool                                                            isale,
        const DRT::UTILS::GaussIntegration &                            intpoints);

    void GaussPointLoopP1OD(
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
                            LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>&                  ecouplp1_u,
                            LINALG::Matrix<my::nen_, my::nen_>&                             ecouplp1_p,
                            Teuchos::RCP<const MAT::Material>                               material,
                            const DRT::UTILS::GaussIntegration &                            intpoints);

    virtual void ComputePorosity( Teuchos::ParameterList& params,
                               const double& press,
                               const double& J,
                               const int& gp,
                               const LINALG::Matrix<my::nen_,1>&       shapfct,
                               const LINALG::Matrix<my::nen_,1>*           myporosity,
                               double& porosity,
                               double* dphi_dp,
                               double* dphi_dJ,
                               double* dphi_dJdp,
                               double* dphi_dJJ,
                               double* dphi_dpp,
                               bool save);

    virtual void ComputePorosityGradient(
        const double&                                      dphidp,
        const double&                                      dphidJ,
        LINALG::Matrix<my::nsd_,1>&                        gradJ,
        const LINALG::Matrix<my::nen_,1>*                  eporositynp,
        LINALG::Matrix<my::nsd_,1>&                        grad_porosity);

    virtual void ComputeContiTimeRHS(Teuchos::ParameterList&                   params,
                             const LINALG::Matrix<my::nen_,1>&         echist,
                             const LINALG::Matrix<my::nen_, 1>&        epren,
                             const LINALG::Matrix<my::nsd_, my::nen_>& edispn,
                             LINALG::Matrix<my::nen_,1>&               preforce,
                             const double&                             rhsfac
                            );
};

}//ELEMENTS
}//DRT


#endif /* FLUID_ELE_CALC_PORO_P1_H_ */
