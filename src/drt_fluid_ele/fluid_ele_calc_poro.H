/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_calc_poro.H

 \brief Internal implementation of poro Fluid element (standard poro fluid)

 <pre>
   \maintainer Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>

 \level 2

 *----------------------------------------------------------------------*/

#ifndef FLUID_ELE_CALC_PORO_H
#define FLUID_ELE_CALC_PORO_H

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "../drt_inpar/inpar_structure.H"

#include "fluid_ele_calc.H"

/*---------------------------------------------------------------------*
 | forward declarations                                                 |
 *---------------------------------------------------------------------*/
namespace MAT
{
  class StructPoro;
}


namespace DRT
{
namespace ELEMENTS
{

/// Class for Evaluating boundary integrals for porous media problems
/*!
 This class is derived from the FluidEleCalc class, i.e. it is capable of
 evaluated all integrals implemented there. It will do so, if the evaluate action
 given by the control routine is not known (see method Evaluate).

 The main methods are the Evaluate and the EvaluateOD routines. Therein,
 the stiffness matrixes of a porous fluid problem are evaluated. OD means
 off diagonal, indicating linearizations with respect to structural degrees of freedom,
 that will be assembled into off diagonal entries in the global system matrix.
 The terms are eventually evaluated in the GaussPointLoop.. methods

 This a calculation class implemented as a singleton, like all calc classes in fluid
 (see comments on base classes for more details). In short this means that on instance
 exists for every discretization type of the boundary element (because of the template).

 \author vuong 10/14
 */

template<DRT::Element::DiscretizationType distype>
class FluidEleCalcPoro : public FluidEleCalc<distype>
{
  typedef DRT::ELEMENTS::FluidEleCalc<distype,DRT::ELEMENTS::Fluid::none> my;

protected:
  /// private Constructor since we are a Singleton.
  FluidEleCalcPoro();

public:
  virtual ~FluidEleCalcPoro() {}

  /// Singleton access method
  static FluidEleCalcPoro<distype>* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  /*!
    \brief calculate element matrix and rhs for porous flow

    \param eid              (i) element id
    \param discretization   (i) fluid discretization the element belongs to
    \param lm               (i) location matrix of element
    \param params           (i) element parameter list
    \param mat              (i) material
    \param elemat1_epetra   (o) element matrix to calculate
    \param elemat2_epetra   (o) element matrix to calculate
    \param elevec1_epetra   (o) element vector to calculate
    \param elevec2_epetra   (o) element vector to calculate
    \param elevec3_epetra   (o) element vector to calculate
    \param offdiag          (i) flag indicating wether diagonal or off diagonal blocks are to be calculated

   */
  virtual int Evaluate(
      DRT::ELEMENTS::Fluid*         ele,
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

  /*!
    \brief calculate element matrix and rhs for porous flow (2)

    \param eid              (i) element id
    \param discretization   (i) fluid discretization the element belongs to
    \param lm               (i) location matrix of element
    \param params           (i) element parameter list
    \param mat              (i) material
    \param elemat1_epetra   (o) element matrix to calculate
    \param elemat2_epetra   (o) element matrix to calculate
    \param elevec1_epetra   (o) element vector to calculate
    \param elevec2_epetra   (o) element vector to calculate
    \param elevec3_epetra   (o) element vector to calculate
    \param intpoints        (i) Gaussian integration points

   */
  virtual int Evaluate(
      DRT::ELEMENTS::Fluid*                ele,
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

  /*!
    \brief Evaluate coupling terms (off diagonal terms) of the element at specified gauss points for porous flow

    \param eid              (i) element id
    \param discretization   (i) fluid discretization the element belongs to
    \param lm               (i) location matrix of element
    \param params           (i) element parameter list
    \param mat              (i) material
    \param elemat1_epetra   (o) element matrix to calculate
    \param elemat2_epetra   (o) element matrix to calculate
    \param elevec1_epetra   (o) element vector to calculate
    \param elevec2_epetra   (o) element vector to calculate
    \param elevec3_epetra   (o) element vector to calculate
    \param intpoints        (i) Gaussian integration points

   */
  virtual int EvaluateOD(
      DRT::ELEMENTS::Fluid*                ele,
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


  /// Evaluate supporting methods of the element
  /*!
      Interface function for supporting methods of the element
   */
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
      Epetra_SerialDenseVector& elevec3
  );

protected:

  /*!
      \brief evaluate function for Fluid element for porous flow

      Specific evaluate function without any knowledge about DRT objects. This
      way the element evaluation is independent of the specific mesh storage.

      \param params           (i) element parameter list
      \param ebofoaf          (i) body force at time step n+alpha_F/n+1
      \param elemat1          (o) element matrix to be filled
      \param elevec1          (o) element rhs vector to be filled
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1
      \param epreaf           (i) nodal pressure at time step n+alpha_F/n+1
      \param evelnp           (i) nodal velocities time step at n+1
      \param eveln            (i) nodal velocities time step at n
      \param eprenp           (i) nodal pressure at time step n+alpha_F/n+1
      \param epren            (i) nodal pressure at time step n
      \param emhist           (i) time rhs for momentum equation
      \param echist           (i) time rhs for continuity equation
      \param epressnp_timederiv (i) nodal pressure time derivative at time step n+alpha_F/n+1
      \param eaccam           (i) nodal accelerations at time step n+alpha_M
      \param edispnp          (i) nodal displacements at time step n+1 (on moving mesh)
      \param edispn           (i) nodal displacements at time step n (on moving mesh)
      \param egridv           (i) grid velocity at time step n+1
      \param egridvn          (i) grid velocity at time step n
      \param escaaf           (i) nodal scalar at time step n+alpha_F/n+1
      \param material         (i) fluid material
      \param isale            (i) ALE flag
      \param intpoints        (i) Gaussian integration points
   */
  int Evaluate( Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       ebofoaf,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_> &   elemat1,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,1> &                       elevec1,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       evelaf,
      const LINALG::Matrix<my::nen_,1>    &                           epreaf,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       evelnp,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       eveln,
      const LINALG::Matrix<my::nen_,1>    &                           eprenp,
      const LINALG::Matrix<my::nen_,1>    &                           epren,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>    &                           epressnp_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                           epressam_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                           epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       eaccam,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       edispn,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       egridv,
      const LINALG::Matrix<my::nsd_,my::nen_> &                       egridvn,
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

          \param params           (i) element parameter list
          \param elemat1          (o) element matrix to be filled
          \param elevec1          (o) element rhs vector to be filled
          \param evelaf           (i) nodal velocities at n+alpha_F/n+1
          \param epreaf           (i) nodal pressure at n+alpha_F/n+1
          \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
          \param eveln            (i) nodal velocities at n (np_genalpha)
          \param eprenp           (i) nodal pressure at n+1
          \param epren            (i) nodal pressure at n
          \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
          \param eaccam           (i) nodal accelerations at n+alpha_M
          \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
          \param edispn           (i) nodal displacements at n (on moving mesh)
          \param egridv           (i) grid velocity at n+1
          \param egridvn          (i) grid velocity at n
          \param escaaf           (i) nodal scalar at n+alpha_F/n+1
          \param eporositynp      (i) nodal porosity at n+1
          \param mat              (i) fluid material
          \param isale            (i) ALE flag
          \param intpoints        (i) Gaussian integration points
   */
  int EvaluateOD(
      Teuchos::ParameterList&                                           params,
      const LINALG::Matrix<my::nsd_,my::nen_> &                         ebofoaf,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, my::nsd_ * my::nen_> &  elemat1,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1> &                    elevec1,
      const LINALG::Matrix<my::nsd_,my::nen_> &                         evelaf,
      const LINALG::Matrix<my::nen_, 1> &                               epreaf,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        eveln,
      const LINALG::Matrix<my::nen_, 1> &                               eprenp,
      const LINALG::Matrix<my::nen_, 1> &                               epren,
      const LINALG::Matrix<my::nen_, 1> &                               epressnp_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                             epressam_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                             epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                          eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        edispn,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        egridv,
      const LINALG::Matrix<my::nsd_, my::nen_> &                        egridvn,
      const LINALG::Matrix<my::nen_,1>&                                 escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                          emhist,
      const LINALG::Matrix<my::nen_,1>&                                 echist,
      const LINALG::Matrix<my::nen_,1>*                                 eporositynp,
      Teuchos::RCP<MAT::Material>                                       mat,
      bool                                                              isale,
      const DRT::UTILS::GaussIntegration &                              intpoints);

  /*!
    \brief calculate element matrix and rhs for porous flow

    \param params           (i) element parameter list
    \param ebofoaf          (i) body force at n+alpha_F/n+1
    \param evelaf           (i) nodal velocities at n+alpha_F/n+1
    \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
    \param eveln            (i) nodal velocities at n (np_genalpha)
    \param epreaf           (i) nodal pressure at n+alpha_F/n+1
    \param eprenp           (i) nodal pressure at n+alpha_F/n+1
    \param epren            (i) nodal pressure at n
    \param eaccam           (i) nodal accelerations at n+alpha_M
    \param emhist           (i) time rhs for momentum equation
    \param echist           (i) time rhs for continuity equation
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
    \param edispn           (i) nodal displacements at n (on moving mesh)
    \param egridv           (i) grid velocity at n+1
    \param egridvn          (i) grid velocity at n
    \param escaaf           (i) nodal scalar at n+alpha_F/n+1
    \param eporositynp      (i) nodal porosity at n+1
    \param eporositydot     (i) nodal porosity time derivative at n+1
    \param eporositydot     (i) nodal porosity time derivative at n
    \param estif            (o) element matrix to calculate
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
      const LINALG::Matrix<my::nsd_,my::nen_>&                      eveln,
      const LINALG::Matrix<my::nen_,1>&                             epreaf,
      const LINALG::Matrix<my::nen_,1>&                             eprenp,
      const LINALG::Matrix<my::nen_,1>&                             epren,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      eaccam,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      emhist,
      const LINALG::Matrix<my::nen_,1>&                             echist,
      const LINALG::Matrix<my::nen_,1>    &                         epressnp_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                         epressam_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                         epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      edispn,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      egridv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                      egridvn,
      const LINALG::Matrix<my::nen_,1>&                             escaaf,
      const LINALG::Matrix<my::nen_,1>*                             eporositynp,
      const LINALG::Matrix<my::nen_,1>*                             eporositydot,
      const LINALG::Matrix<my::nen_,1>*                             eporositydotn,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,(my::nsd_+1)*my::nen_>&  estif,
      LINALG::Matrix<(my::nsd_+1)*my::nen_,1>&                      eforce,
      Teuchos::RCP<const MAT::Material>                             material,
      bool                                                          isale,
      const DRT::UTILS::GaussIntegration &                          intpoints);

  /*!
    \brief calculate off diagonal element matrix and rhs for porous flow

    \param params           (i) element parameter list
    \param evelaf           (i) nodal velocities at n+alpha_F/n+1
    \param evelnp           (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf           (i) nodal pressure at n+alpha_F/n+1
    \param eprenp           (i) nodal pressure at n+alpha_F/n+1
    \param epren            (i) nodal pressure at n
    \param eaccam           (i) nodal accelerations at n+alpha_M
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param edispnp          (i) nodal displacements at n+1 (on moving mesh)
    \param edispn           (i) nodal displacements at n (on moving mesh)
    \param egridv           (i) grid velocity at n+1
    \param egridvn          (i) grid velocity at n
    \param escaaf           (i) nodal scalar at n+alpha_F/n+1
     param emhist           (i) time rhs for momentum equation
    \param echist           (i) time rhs for continuity equation
    \param eporositynp      (i) nodal porosity at n+1
    \param ecoupl           (o) element matrix to calculate
    \param eforce           (o) element rhs to calculate
    \param material         (i) fluid material
    \param isale            (i) ALE flag
    \param intpoints        (i) Gaussian integration points
   */
  void SysmatOD(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1>&                              epren,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                           epressam_timederiv,
      const LINALG::Matrix<my::nen_,1>    &                           epressn_timederiv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>*                               eporositynp,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_,my::nsd_ * my::nen_>&  ecoupl,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
      Teuchos::RCP<const MAT::Material>                               material,
      bool                                                            isale,
      const DRT::UTILS::GaussIntegration &                            intpoints);

  /*!
    \brief linearisation of momentum equation in the case of mesh motion 3-D for Poroelasticity

    \param ecoupl_u           (o) coupling element matrix for fluid velocity
    \param dphi_dp            (i) derivative of porosity w.r.t. fluid pressure
    \param dphi_dJ            (i) derivative of porosity w.r.t. jacobian determinant
    \param refporositydot     (i) time derivative of reference porosity
    \param timefac            (i) time factor
    \param timefacfac         (i) time factor * integration factor
   */
  void LinMeshMotion_3D_OD(
      LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>&         ecoupl_u,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    refporositydot,
      const double &                                                    timefac,
      const double &                                                    timefacfac
  );

  /*!
    \brief linearisation of momentum equation in the case of mesh motion 2-D for Poroelasticity

    \param ecoupl_u           (o) coupling element matrix for fluid velocity
    \param dphi_dp            (i) derivative of porosity w.r.t. fluid pressure
    \param dphi_dJ            (i) derivative of porosity w.r.t. jacobian determinant
    \param refporositydot     (i) time derivative of reference porosity
    \param timefac            (i) time factor
    \param timefacfac         (i) time factor * integration factor
   */
  void LinMeshMotion_2D_OD(
      LINALG::Matrix<my::nsd_ * my::nen_, my::nsd_ * my::nen_>&         ecoupl_u,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    refporositydot,
      const double &                                                    timefac,
      const double &                                                    timefacfac
  );

  /*!
    \brief linearisation of continuity equation in the case of mesh motion 2-D for Poroelasticity

    \param ecoupl_p           (o) coupling element matrix for fluid pressure
    \param dphi_dp            (i) derivative of porosity w.r.t. fluid pressure
    \param dphi_dJ            (i) derivative of porosity w.r.t. jacobian determinant
    \param refporositydot     (i) time derivative of reference porosity
    \param timefacfacpre      (i) time factor * integration factor
   */
  void LinMeshMotion_2D_Pres_OD(
      LINALG::Matrix< my::nen_, my::nsd_ * my::nen_>&                   ecoupl_p,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    refporositydot,
      const double &                                                    timefacfacpre
  );

  /*!
    \brief linearisation of continuity equation in the case of mesh motion 3-D for Poroelasticity

    \param ecoupl_p           (o) coupling element matrix for fluid pressure
    \param dphi_dp            (i) derivative of porosity w.r.t. fluid pressure
    \param dphi_dJ            (i) derivative of porosity w.r.t. jacobian determinant
    \param refporositydot     (i) time derivative of reference porosity
    \param timefacfacpre      (i) time factor * integration factor
   */
  virtual void LinMeshMotion_3D_Pres_OD(
      LINALG::Matrix< my::nen_, my::nsd_ * my::nen_>&                   ecoupl_p,
      const double &                                                    dphi_dp,
      const double &                                                    dphi_dJ,
      const double &                                                    refporositydot,
      const double &                                                    timefacfacpre
  );

  //! Compute element matrix entries: PSPG
  virtual void PSPG(
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &        estif_q_u,      ///< block (weighting function q x u)
      LINALG::Matrix<my::nen_,my::nen_> &                  ppmat,          ///< block (weighting function q x p)
      LINALG::Matrix<my::nen_,1> &                         preforce,       ///< rhs forces pressure
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &   lin_resM_Du,    ///< linearisation of the stabilization residual
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&    lin_resMRea_Du, ///< linearisation of the stabilization residual (reactive terms)
      const LINALG::Matrix<my::nsd_,my::nen_> &            lin_resM_Dp,    ///< linearisation of the stabilization residual w.r.t. pressure
      const double &                                       dphi_dp,        ///< linearisation of porosity w.r.t. pressure
      const double &                                       fac3,           ///< factor for residual in current subgrid velocities
      const double &                                       timefacfac,     ///< = timefac x fac
      const double &                                       timefacfacpre,  ///< = timefacpre x fac
      const double &                                       rhsfac          ///< right-hand-side factor for residuals
  );

  //! Compute element matrix entries: Biot
  virtual void StabBiot(
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_> &        estif_q_u,      ///< block (weighting function q x u)
      LINALG::Matrix<my::nen_,my::nen_> &                  ppmat,          ///< block (weighting function q x p)
      LINALG::Matrix<my::nen_,1> &                         preforce,       ///< rhs forces pressure
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &   lin_resM_Du,    ///< linearisation of the stabilization residual
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&    lin_resMRea_Du, ///< linearisation of the stabilization residual (reactive terms)
      const LINALG::Matrix<my::nsd_,my::nen_> &            lin_resM_Dp,    ///< linearisation of the stabilization residual w.r.t. pressure
      const double &                                       dphi_dp,        ///< linearisation of porosity w.r.t. pressure
      const double &                                       fac3,           ///< factor for residual in current subgrid velocities
      const double &                                       timefacfac,     ///< = timefac x fac
      const double &                                       timefacfacpre,  ///< = timefacpre x fac
      const double &                                       rhsfac          ///< right-hand-side factor for residuals
  );

  //! Compute element matrix entries: reactive stabilization
  virtual void ReacStab(
      LINALG::Matrix<my::nen_*my::nsd_,my::nen_*my::nsd_> &  estif_u,        ///< block (weighting function v x u)
      LINALG::Matrix<my::nen_*my::nsd_,my::nen_> &           estif_p_v,      ///< block (weighting function v x p)
      LINALG::Matrix<my::nsd_,my::nen_> &                    velforce,       ///< rhs forces velocity
      LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> &           lin_resM_Du,    ///< linearisation of the stabilization residual
      const LINALG::Matrix<my::nsd_,my::nen_> &              lin_resM_Dp,    ///< linearisation of the stabilization residual w.r.t. pressure
      const double &                                         dphi_dp,        ///< linearisation of porosity w.r.t. pressure
      const double &                                         timefacfac,     ///< = timefac x fac
      const double &                                         timefacfacpre,  ///< = timefacpre x fac
      const double &                                         rhsfac,         ///< right-hand-side factor for residuals
      const double &                                         fac3            ///< factor for residual in current subgrid velocities
  );

  //! Compute linearization of momentum residual w.r.t fluid velocities
  void ComputeLinResMDu(
      const double & timefacfac,
      LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&  lin_resM_Du,
      LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&  lin_resMRea_Du);

  //! Compute linearization of momentum residual (stabilization) w.r.t fluid velocities
  void ComputeLinResMDu_Stab(
      const double & timefacfac,
      LINALG::Matrix<my::nsd_*my::nsd_,my::nen_>&  lin_resM_Du);

  //! Compute linearization of momentum residual w.r.t fluid pressure
  void ComputeLinResMDp(
      const double & timefacfacpre,
      const double & dphi_dp,
      LINALG::Matrix<my::nsd_,my::nen_>&  lin_resM_Du);

  //! calculate div(epsilon(u))
  void CalcDivEps(const LINALG::Matrix<my::nsd_,my::nen_>& evelaf   ///< velocity at time n+alpha_f / n+1
      );

  /*!
    \brief Compute derivatives of deformation gradient

    \param edispnp         (i) nodal displacements at n+1
    \param defgrd_inv      (i) inverse deformation gradient
    \param F_x             (o) derivative of deformation gradient w.r.t. current coordinates xyz
    \param F_X             (o) derivative of deformation gradient w.r.t. material coordinates XYZ
   */
  void ComputeFDerivative(
      const LINALG::Matrix<my::nsd_, my::nen_>&         edispnp,
      const LINALG::Matrix<my::nsd_,my::nsd_>&          defgrd_inv,
      LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&       F_x,
      LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&       F_X
  );

  /*!
    \brief Compute spatial gradient of jacobian determinant and porosity

    \param J              (i) determinant of deformation gradient at gauss point
    \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
    \param dphidJ         (i) derivative of porosity w.r.t. jacobian determinant
    \param defgrd_IT_vec  (i) inverse transposed deformation gradient in vector notation
    \param F_x            (i) derivative of deformation gradient w.r.t. current coordinates xyz
    \param gradp          (i) spatial gradient of fluid pressure
    \param eporositynp    (i) nodal porosities at n+1
    \param gradJ          (o) spatial gradient of jacobian determinant
    \param grad_porosity  (o) spatial gradient of porosity
    \param refgrad_porosity  (o) reference gradient of porosity
   */
  void ComputeGradients(
      const double&                                      J,
      const double&                                      dphidp,
      const double&                                      dphidJ,
      const LINALG::Matrix<my::nsd_*my::nsd_,1>&         defgrd_IT_vec,
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>&  F_x,
      const LINALG::Matrix<my::nsd_,1>&                  gradp,
      const LINALG::Matrix<my::nen_,1>*                  eporositynp,
      LINALG::Matrix<my::nsd_,1>&                        gradJ,
      LINALG::Matrix<my::nsd_,1>&                        grad_porosity,
      LINALG::Matrix<my::nsd_,1>&                        refgrad_porosity);

  /*!
    \brief Compute spatial gradient of porosity

    \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
    \param dphidJ         (i) derivative of porosity w.r.t. jacobian determinant
    \param gradJ          (i) spatial gradient of jacobian determinant
    \param gradJ          (i) spatial gradient of fluid pressure
    \param eporositynp    (i) nodal porosities at n+1
    \param grad_porosity  (o) spatial gradient of porosity
    \param refgrad_porosity  (o) reference gradient of porosity
   */
  virtual void ComputePorosityGradient(
      const double&                                      dphidp,
      const double&                                      dphidJ,
      const LINALG::Matrix<my::nsd_,1>&                  gradJ,
      const LINALG::Matrix<my::nsd_,1>&                  gradp,
      const LINALG::Matrix<my::nen_,1>*                  eporositynp,
      LINALG::Matrix<my::nsd_,1>&                        grad_porosity,
      LINALG::Matrix<my::nsd_,1>&                        refgrad_porosity);

  /*!
    \brief Compute linearization needed for diagonal terms (lin. of porosity gradient w.r.t. fluid pressure)

    \param dphidp         (i) derivative of porosity w.r.t. fluid pressure
    \param dphi_dpp       (i) second derivative of porosity w.r.t. fluid pressure
    \param dphi_dJp       (i) mixed derivative of porosity w.r.t. fluid pressure and jacobian determinant
    \param gradJ          (i) spatial gradient of jacobian determinant
    \param dgradphi_dp    (o) derivate of spatial gradient of porosity w.r.t. fluid pressure
   */
  virtual void ComputeLinearization(
      const double&                                      dphi_dp,
      const double&                                      dphi_dpp,
      const double&                                      dphi_dJp,
      const LINALG::Matrix<my::nsd_,1>&                  gradJ,
      LINALG::Matrix<my::nsd_,my::nen_>&                 dgradphi_dp);

  /*!
    \brief Compute linearization needed for off diagonal terms
      (lin. of jacobian, porosity and porosity gradient w.r.t. structure displacement)

    \param dphi_dJ        (i) derivative of porosity w.r.t. jacobian determinant
    \param dphi_dJJ       (i) second derivative of porosity w.r.t. jacobian determinant
    \param dphi_dJp       (i) mixed derivative of porosity w.r.t. fluid pressure and jacobian determinant
    \param defgrd_inv     (i) inverse deformation gradient
    \param defgrd_IT_vec  (o) inverse transposed deformation gradient in vector notation
    \param F_x            (i) derivative of deformation gradient w.r.t. current coordinates xyz
    \param F_X            (i) derivative of deformation gradient w.r.t. material coordinates XYZ
    \param gradJ          (i) spatial gradient of jacobian determinant
    \param dJ_dus         (o) derivative of jacobian determinant w.r.t. structure displacments
    \param dphi_dus       (o) derivative of porosity determinant w.r.t. structure displacments
    \param dgradphi_dus   (o) derivate of spatial gradient of porosity w.r.t. structure displacments
   */
  virtual void ComputeLinearizationOD(
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

  /*!
    \brief Compute porosity and derivatives

    \param params         (i) element parameter list
    \param press          (i) fluid pressure at gauss point
    \param J              (i) jacobian determinant at gauss point
    \param gp             (i) number of actual gauss point
    \param shapfct        (i) shape function values at gauss point
    \param myporosity     (i) nodal porosities
    \param porosity       (o) porosity at gauss point
    \param dphi_dp        (o) derivative of porosity gradient w.r.t. fluid pressure
    \param dphi_dJ        (o) derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dJdp      (o) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
    \param dphi_dJJ       (o) second derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dpp       (o) second derivative of porosity gradient w.r.t. fluid pressure
    \param save           (i) flag for saving porosity within structure material
   */
  virtual void ComputePorosity(
      Teuchos::ParameterList&                params,
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

  /*!
   \brief evaluate pressure equation (i.e. continuity equation for standard poro elements)
          This function is overwritten by the poro_p1 element

    \param params         (i) element parameter list
    \param timefacfacpre  (i) fluid pressure at gauss point
    \param rhsfac         (i) jacobian determinant at gauss point
    \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
    \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
    \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dpp       (i) second derivative of porosity gradient w.r.t. fluid pressure
    \param eporositydot   (i) nodal values of porosity time derivative at time step n+1
    \param eporositydotn  (i) nodal values of porosity time derivative at time step n
    \param echist         (i) nodal values of history values of continuity equation
    \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
    \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
    \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
    \param preforce       (o) element rhs vector to be filled
   * */
  virtual void EvaluatePressureEquation(
      Teuchos::ParameterList&                       params,
      const double&                                 timefacfacpre,
      const double&                                 rhsfac,
      const double&                                 dphi_dp,
      const double&                                 dphi_dJ,
      const double&                                 dphi_dJdp,
      const double&                                 dphi_dpp,
      const LINALG::Matrix<my::nen_,1>*             eporositydot,
      const LINALG::Matrix<my::nen_,1>*             eporositydotn,
      const LINALG::Matrix<my::nen_,1>&             echist,
      const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
      LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
      LINALG::Matrix<my::nen_,1>&                   preforce);

  /*!
   \brief evaluate nontransient part of pressure equation (i.e. no time derivative of porosity)

    \param params         (i) element parameter list
    \param timefacfacpre  (i) fluid pressure at gauss point
    \param rhsfac         (i) jacobian determinant at gauss point
    \param dphi_dp        (i) derivative of porosity gradient w.r.t. fluid pressure
    \param dphi_dJ        (i) derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dJdp      (i) mixed derivative of porosity gradient w.r.t. fluid pressure and jacobian determinant
    \param dphi_dJJ       (i) second derivative of porosity gradient w.r.t. jacobian determinant
    \param dphi_dpp       (i) second derivative of porosity gradient w.r.t. fluid pressure
    \param echist         (i) nodal values of history values of continuity equation
    \param dgradphi_dp    (i) derivative of porosity gradient w.r.t. fluid pressure
    \param estif_q_u      (o) element matrix (pressure - fluid velocity weighting) to be filled
    \param ppmat          (o) element matrix (pressure - pressure weighting) to be filled
    \param preforce       (o) element rhs vector to be filled
   * */
  virtual void EvaluatePressureEquationNonTransient(
      Teuchos::ParameterList&                       params,
      const double&                                 timefacfacpre,
      const double&                                 rhsfac,
      const double&                                 dphi_dp,
      const double&                                 dphi_dJ,
      const double&                                 dphi_dJdp,
      const double&                                 dphi_dpp,
      const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
      LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
      LINALG::Matrix<my::nen_,1>&                   preforce);

  /*!
    \brief Gauss point loop for evaluation of diagonal terms

    \param params             (i) element parameter list
    \param evelaf             (i) nodal velocities at n+alpha_F/n+1
    \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf             (i) nodal pressure at n+alpha_F/n+1
    \param eprenp             (i) nodal pressure at n+alpha_F/n+1
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param eaccam             (i) nodal accelerations at n+alpha_M
    \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
    \param egridv             (i) grid velocity at n+1
    \param escaaf             (i) nodal scalar at n+alpha_F/n+1
    \param eporositynp        (i) nodal porosity at n+1
    \param eporositydot       (i) nodal porosity time derivative at n+1
    \param eporositydotn      (i) nodal porosity time derivative at n
    \param estif_u            (o) element matrix for fluid velocity
    \param estif_p_v          (o) element matrix for fluid velocity weighting - fluid pressure
    \param estif_q_u          (o) element matrix for pressure weighting - fluid velocity
    \param ppmat              (o) element matrix for fluid pressure
    \param preforce           (o) element rhs vector for fluid pressure
    \param velforce           (o) element rhs vector for fluid velocity
    \param material           (i) fluid material
    \param intpoints          (i) Gaussian integration points
   */
  void GaussPointLoop(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1>&                              epren,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        edispn,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>*                               eporositynp,
      const LINALG::Matrix<my::nen_,1>*                               eporositydot,
      const LINALG::Matrix<my::nen_,1>*                               eporositydotn,
      LINALG::Matrix<my::nen_*my::nsd_,my::nen_*my::nsd_>&            estif_u,
      LINALG::Matrix<my::nen_*my::nsd_,my::nen_>&                     estif_p_v,
      LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&                    estif_q_u,
      LINALG::Matrix<my::nen_,my::nen_>&                              ppmat,
      LINALG::Matrix<my::nen_,1>&                                     preforce,
      LINALG::Matrix<my::nsd_,my::nen_>&                              velforce,
      Teuchos::RCP<const MAT::Material>                               material,
      const DRT::UTILS::GaussIntegration &                            intpoints);

  /*!
    \brief Gauss point loop for evaluation of off-diagonal terms

    \param params             (i) element parameter list
    \param evelaf             (i) nodal velocities at n+alpha_F/n+1
    \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf             (i) nodal pressure at n+alpha_F/n+1
    \param eprenp             (i) nodal pressure at n+alpha_F/n+1
    \param epren              (i) nodal pressure at n
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param eaccam             (i) nodal accelerations at n+alpha_M
    \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
    \param edispn             (i) nodal displacements at n
    \param egridv             (i) grid velocity at n+1
    \param egridvn            (i) grid velocity at n
    \param escaaf             (i) nodal scalar at n+alpha_F/n+1
    \param eporositynp        (i) nodal porosity at n+alpha_F/n+1
    \param eforce             (o) coupling rhs force vector
    \param ecoupl_u           (o) coupling element matrix for fluid velocity
    \param ecoupl_p           (o) coupling element matrix for fluid pressure
    \param material           (i) fluid material
    \param intpoints          (i) Gaussian integration points
   */
  void GaussPointLoopOD(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1>&                              epren,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>*                               eporositynp,
      LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&                   eforce,
      LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&       ecoupl_u,
      LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&                  ecoupl_p,
      Teuchos::RCP<const MAT::Material>                               material,
      const DRT::UTILS::GaussIntegration &                            intpoints);

  /*!
    \brief Evaluation of gauss point values (diagonal terms)

    \param params             (i) element parameter list
    \param evelaf             (i) nodal velocities at n+alpha_F/n+1
    \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf             (i) nodal pressure at n+alpha_F/n+1
    \param eprenp             (i) nodal pressure at n+alpha_F/n+1
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
    \param edispn             (i) nodal displacements at n (on moving mesh)
    \param egridv             (i) grid velocity at n+1
    \param egridvn            (i) grid velocity at n
    \param escaaf             (i) nodal scalar at n+alpha_F/n+1
    \param eporositynp        (i) nodal porosity at n+1
    \param eporositydot       (i) nodal porosity time derivative at n+1
    \param eporositydotn      (i) nodal porosity time derivative at n
   */
  void EvaluateVariablesAtGaussPoint(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1>&                              epren,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>*                               eporositynp,
      const LINALG::Matrix<my::nen_,1>*                               eporositydot,
      const LINALG::Matrix<my::nen_,1>*                               eporositydotn);

  /*!
    \brief Evaluation of gauss point values (off-diagonal terms)

    \param params             (i) element parameter list
    \param evelaf             (i) nodal velocities at n+alpha_F/n+1
    \param evelnp             (i) nodal velocities at n+1 (np_genalpha)
    \param epreaf             (i) nodal pressure at n+alpha_F/n+1
    \param eprenp             (i) nodal pressure at n+alpha_F/n+1
    \param epren              (i) nodal pressure at n
    \param epressnp_timederiv (i) nodal pressure time derivative at n+alpha_F/n+1
    \param edispnp            (i) nodal displacements at n+1 (on moving mesh)
    \param edispn             (i) nodal displacements at n (on moving mesh)
    \param egridv             (i) grid velocity at n+1
    \param egridvn            (i) grid velocity at n
    \param escaaf             (i) nodal scalar at n+alpha_F/n+1
    \param eporositynp        (i) nodal porosity at n+alpha_F/n+1
   */
  void EvaluateVariablesAtGaussPointOD(
      Teuchos::ParameterList&                                         params,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        ebofoaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       evelnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       eveln,
      const LINALG::Matrix<my::nen_, 1>&                              epreaf,
      const LINALG::Matrix<my::nen_, 1>&                              eprenp,
      const LINALG::Matrix<my::nen_, 1>&                              epren,
      const LINALG::Matrix<my::nen_, 1> &                             epressnp_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressam_timederiv,
      const LINALG::Matrix<my::nen_, 1> &                             epressn_timederiv,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        eaccam,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       edispn,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridv,
      const LINALG::Matrix<my::nsd_, my::nen_>&                       egridvn,
      const LINALG::Matrix<my::nen_,1>&                               escaaf,
      const LINALG::Matrix<my::nsd_,my::nen_>&                        emhist,
      const LINALG::Matrix<my::nen_,1>&                               echist,
      const LINALG::Matrix<my::nen_,1>*                               eporositynp);

  /*!
    \brief Evaluate off diagonal terms in momentum equation

    \param timefacfac         (i) time factor * integration factor
    \param porosity           (i) porosity at gauss point
    \param gridvelint         (i) grid (structure) velocity at gauss point
    \param grad_porosity      (i) porosity gradient at gauss point
    \param dgradphi_dus       (i) derivative of porosity gradient w.r.t. structural displacements at gauss point
    \param dphi_dus           (i) derivative of porosity w.r.t. structural displacements at gauss point
    \param refporositydot     (i) time derivative of reference porosity
    \param lin_resM_Dus       (i) linearization of residual of momentum equation w.r.t. structureal dofs
    \param ecoupl_u           (o) coupling element matrix of momentum equation
   */
  void FillMatrixMomentumOD(
      const double&                                               timefacfac,
      const LINALG::Matrix<my::nsd_, my::nen_>&                   evelaf,
      const LINALG::Matrix<my::nsd_, my::nen_>&                   egridv,
      const LINALG::Matrix<my::nen_, 1>&                          epreaf,
      const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
      const double &                                              dphi_dp,
      const double &                                              dphi_dJ,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
      const double &                                              refporositydot,
      LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&              lin_resM_Dus,
      LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&              lin_resM_Dus_gridvel,
      LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&   ecoupl_u);

  /*!
    \brief Evaluate off diagonal terms in continuity equation

    \param timefacfacpre      (i) time factor (pressure) * integration factor
    \param dphi_dJ            (i) derivative of porosity gradient w.r.t. jacobian determinant at gauss point
    \param dphi_dJJ           (i) second derivative of porosity gradient w.r.t. jacobian determinant at gauss point
    \param dphi_dJdp          (i) mixed derivative of porosity gradient w.r.t. jacobian and pressure at gauss point
    \param dgradphi_dus       (i) derivative of porosity gradient w.r.t. structural displacements at gauss point
    \param dphi_dus           (i) derivative of porosity w.r.t. structural displacements at gauss point
    \param dJ_dus             (i) derivative of jacobian w.r.t. structural displacements at gauss point
    \param egridv             (i) nodal grid velocities
    \param lin_resM_Dus       (i) linearization of residual of momentum equation w.r.t. structureal dofs
    \param ecoupl_u           (o) coupling element matrix of continuity equation
   */
  virtual void FillMatrixContiOD(
      const double&                                               timefacfacpre,
      const double &                                              dphi_dp,
      const double &                                              dphi_dJ,
      const double&                                               dphi_dJJ,
      const double&                                               dphi_dJdp,
      const double &                                              refporositydot,
      const LINALG::Matrix<my::nsd_,my::nen_*my::nsd_>&           dgradphi_dus,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dphi_dus,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&                  dJ_dus,
      const LINALG::Matrix<my::nsd_, my::nen_>&                   egridv,
      const LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&        lin_resM_Dus,
      const LINALG::Matrix<my::nsd_, my::nen_ * my::nsd_>&        lin_resM_Dus_gridvel,
      LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&              ecoupl_p);

  //! do some evaluation before actual element matrix assembly
  void PreEvaluate
  (
      Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::ELEMENTS::Fluid*                 ele, ///< fluid element
      DRT::Discretization&      discretization  ///< pointer to discretization for de-assembly
  );

  /// computation of material derivatives
  double SetupMaterialDerivatives();

  /// access structure material of corresponding solid (poro) element
  void GetStructMaterial(DRT::ELEMENTS::Fluid*                 ele);

  /// get material parameters of poro fluid element
  void GetMaterialParamters(Teuchos::RCP<const MAT::Material>   material);

  /// compute spatial reactive term (darcy term)
  void ComputeSpatialReactionTerms(
      Teuchos::RCP<const MAT::Material>        material,   //< fluid material
      const LINALG::Matrix<my::nsd_,my::nsd_>& invdefgrd   //< inverse of deformationgradient at gausspoint
      );

  /// compute linearization of spatial reactive term (darcy term) w.r.t to structural displacements
  void ComputeLinSpatialReactionTerms(
      Teuchos::RCP<const MAT::Material>              material,   //< fluid material
      const LINALG::Matrix<my::nsd_,my::nsd_>&       defgrd_inv, //< inverse of deformationgradient at gausspoint
      const LINALG::Matrix<1,my::nsd_*my::nen_>*     dJ_dus,     //< derivative of jacobian w.r.t. structural displacements at gauss point
      const LINALG::Matrix<1,my::nsd_*my::nen_>*     dphi_dus    //< derivative of porosity w.r.t. structural displacements at gauss point
      );

  /// get compute RHS of momentum equation of time step n and subgrid-scale velocity
  void ComputeOldRHSAndSubgridScaleVelocity();

  /// get stabilization paramters
  void ComputeStabilizationParameters(const double& vol);

  /// get compute RHS of contiuity equation of time step n
  void ComputeOldRHSConti(double dphi_dp);

  void ComputeMixtureStrongResidual(
      Teuchos::ParameterList&                           params,
      const LINALG::Matrix<my::nsd_,my::nsd_>&          defgrd,
      const LINALG::Matrix<my::nsd_, my::nen_>&         edispnp,
      const LINALG::Matrix<my::nsd_, my::nen_>&         edispn,
      const LINALG::Matrix<my::nsd_*my::nsd_,my::nsd_>& F_X,
      const bool                                        computeLinOD
      );

  virtual int ComputeVolume(
      Teuchos::ParameterList&         params,           //< paramters
      DRT::ELEMENTS::Fluid*           ele,              //< current fluid element
      DRT::Discretization&            discretization,   //< fluid discretization
      std::vector<int>&               lm,               //< location vector for DOF management
      Epetra_SerialDenseVector&       elevec1           //< reference to element vector to be filled
  );

  //! Compute deformation gradient
  void ComputeDefGradient(
      LINALG::Matrix<my::nsd_,my::nsd_>&       defgrd,   ///<<    (i) deformation gradient at gausspoint
      const LINALG::Matrix<my::nsd_,my::nen_>& N_XYZ,    ///<<    (i) derivatives of shape functions w.r.t. reference coordinates
      const LINALG::Matrix<my::nsd_,my::nen_>& xcurr     ///<<    (i) current position of gausspoint
    );

  //! Compute Jacobian Determinant and the volume change
  void ComputeJacobianDeterminantVolumeChange(
     double & J,
     double & volchange,
     const LINALG::Matrix<my::nsd_,my::nsd_>& defgrd,
     const LINALG::Matrix<my::nsd_,my::nen_>& N_XYZ,
     const LINALG::Matrix<my::nsd_,my::nen_>& nodaldisp
     );

  //! Compute deformation gradient
  double ComputeEffectiveStiffness(
    );

  /// Evaluate element ERROR
  /*!
      general function to compute the error (analytical solution) for particular problem type
   */
  int ComputeError(
      DRT::ELEMENTS::Fluid*         ele,
      Teuchos::ParameterList&       params,
      Teuchos::RCP<MAT::Material>&  mat,
      DRT::Discretization&          discretization,
      std::vector<int>&             lm,
      Epetra_SerialDenseVector&     elevec);

  int ComputeError(
      DRT::ELEMENTS::Fluid*                ele,
      Teuchos::ParameterList&              params,
      Teuchos::RCP<MAT::Material>&         mat,
      DRT::Discretization&                 discretization,
      std::vector<int>&                    lm,
      Epetra_SerialDenseVector&            elevec1,
      const DRT::UTILS::GaussIntegration & intpoints2);

  ///first derivatives of shape functions w.r.t. material coordinates
  LINALG::Matrix<my::nsd_,my::nen_> N_XYZ_;
  /*!
  \brief second derivatives of shape functions w.r.t. material coordinates XYZ
    3D: (N,XX ; N,YY ; N,ZZ ; N,XY ; N,XZ ; N,YX ; N,YZ;  N,ZX ; N,ZY)
    2D: (N,XX ; N,YY ; N,XY ; N,YX )
   */
  LINALG::Matrix<my::numderiv2_,my::nen_> N_XYZ2_;
  ///second derivatives of shape functions w.r.t. material coordinates XYZ (ordered in symmetric matrix)
  LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> N_XYZ2full_;
  /// material coordinates
  LINALG::Matrix<my::nsd_,my::nen_> xyze0_;
  /// node coordinates at time n
  LINALG::Matrix<my::nsd_,my::nen_> xyzeold_;
  //! global derivatives of shape functions w.r.t x,y,z at time n
  LINALG::Matrix<my::nsd_,my::nen_> derxyold_;

  //! vector containing all values from previous timelevel n for continuity equation
  double histcon_;

  //! porosity at gauss point (time step n+1)
  double             porosity_;
  //! time derivative of porosity at gauss point (time step n+1 or n+\alpha_M)
  double             porositydot_;
  //! porosity gradient w.r.t. spatial coordinates at gauss point
  LINALG::Matrix<my::nsd_,1>             grad_porosity_;
  //! porosity gradient w.r.t. reference coordinates at gauss point
  LINALG::Matrix<my::nsd_,1>             refgrad_porosity_;
  //! porosity gradient w.r.t. spatial coordinates at gauss point
  LINALG::Matrix<my::nsd_,1>             grad_porosityold_;
  //! porosity gradient w.r.t. reference coordinates at gauss point
  LINALG::Matrix<my::nsd_,1>             refgrad_porosityold_;
  //! grid (=structure) velocity at gauss point (time step n+1)
  LINALG::Matrix<my::nsd_,1>             gridvelint_;
  //! grid (=structure) velocity at gauss point (time step n)
  LINALG::Matrix<my::nsd_,1>             gridvelnint_;
  //! convective velocity (fluid - structure) at gauss point
  LINALG::Matrix<my::nsd_,1>             convel_;
  //! grid (=structure) velocity derivatives w.r.t. reference coordinates at integration point
  LINALG::Matrix<my::nsd_,my::nsd_>             gridvelderiv_;
  //! grid (=structure) velocity divergence at gauss point
  double gridvdiv_;
  //! determinant of deformation gradient (time step n+1)
  double J_;
  //! fluid pressure at gauss point (time step n+1)
  double press_;
  //! fluid pressure time derivative at gauss point (time step n+1 or n+\alpha_M)
  double pressdot_;
//  //! fluid pressure time derivative at gauss point (time step n)
//  double pressdotn_;
  //! fluid pressure time derivative at gauss point (time step n+1)
  //double pressdotnp_;
  //! fluid pressure gradient w.r.t to paramter space coordinates at gauss point
  LINALG::Matrix<my::nsd_,1> refgradp_;

  //! material reactive tensor
  LINALG::Matrix<my::nsd_,my::nsd_>             matreatensor_;

  //! linearisation of material reactive tensor w.r.t. porosity
  LINALG::Matrix<my::nsd_,my::nsd_>             matreatensorlinporosity_;

  //! linearisation of material reactive tensor w.r.t. J
  LINALG::Matrix<my::nsd_,my::nsd_>             matreatensorlinJ_;

  //! spatial reactive tensor
  LINALG::Matrix<my::nsd_,my::nsd_>             reatensor_;

  //! linearisation of reactive tensor w.r.t. structural displacements * fluid velocity
  LINALG::Matrix<my::nsd_,my::nsd_ * my::nen_>             reatensorlinODvel_;

  //! linearisation of reactive tensor w.r.t. structural displacements * grid velocity
  LINALG::Matrix<my::nsd_,my::nsd_ * my::nen_>             reatensorlinODgridvel_;

  //! reactive tensor x fluid velocity
  LINALG::Matrix<my::nsd_,1>             reavel_;

  //! reactive tensor x structural (grid) velocity
  LINALG::Matrix<my::nsd_,1>             reagridvel_;

  //! reactive tensor x convective velocity
  LINALG::Matrix<my::nsd_,1>             reaconvel_;

  //! linearisation of reaction tensor w.r.t. porosity * fluid velocity
  LINALG::Matrix<my::nsd_,1> lin_p_vel_;

  //! linearisation of reaction tensor w.r.t. porosity * grid velocity
  LINALG::Matrix<my::nsd_,1> lin_p_vel_grid_;

  //! linearization of stabilisation parameters w.r.t. porosity -> it is a (3,1) vector for 2D and 3D
  LINALG::Matrix<3,1> dtaudphi_;

  //! the stabilisation parameter  for biot stabilization
  double taustruct_;

  //! residual of mixture (structural) equation
  LINALG::Matrix<my::nsd_,1>            mixres_;

  //! linearisation residual of mixture (structural) equation
  LINALG::Matrix<my::nsd_,my::nsd_ * my::nen_>            mixresLinOD_;

  //! corresponding poro structure material
  Teuchos::RCP< MAT::StructPoro > structmat_;

  //! state if reaction/permeability tensor is constant
  bool const_permeability_;

  //! kinematic type
  INPAR::STR::KinemType kintype_;

  //! pointer to parameter lists
  DRT::ELEMENTS::FluidEleParameterPoro* porofldpara_;
};
}
}

#endif

