/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_H_
#define SCATRA_ELE_CALC_PORO_H_

#include "scatra_ele_calc.H"

namespace MAT
{
class ScatraMat;
}

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcPoro : public virtual ScaTraEleCalc<distype>
{
protected:
  /// (private) protected constructor, since we are a Singleton.
  ScaTraEleCalcPoro(const int numdofpernode,const int numscal);

  private:
  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcPoro() {}

  /// Singleton access method
  static ScaTraEleCalcPoro<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

  /// Evaluate the element
   /*!
     Generic virtual interface function. Called via base pointer.
    */
//   virtual int Evaluate(DRT::ELEMENTS::Transport*     ele,
//                        Teuchos::ParameterList&       params,
//                        DRT::Discretization &         discretization,
//                        const std::vector<int> &      lm,
//                        Epetra_SerialDenseMatrix&     elemat1_epetra,
//                        Epetra_SerialDenseMatrix&     elemat2_epetra,
//                        Epetra_SerialDenseVector&     elevec1_epetra,
//                        Epetra_SerialDenseVector&     elevec2_epetra,
//                        Epetra_SerialDenseVector&     elevec3_epetra);

protected:

//   int EvaluateODMesh(DRT::ELEMENTS::Transport*     ele,
//                        Teuchos::ParameterList&       params,
//                        DRT::Discretization &         discretization,
//                        const std::vector<int> &      lm,
//                        Epetra_SerialDenseMatrix&     elemat1_epetra,
//                        Epetra_SerialDenseMatrix&     elemat2_epetra,
//                        Epetra_SerialDenseVector&     elevec1_epetra,
//                        Epetra_SerialDenseVector&     elevec2_epetra,
//                        Epetra_SerialDenseVector&     elevec3_epetra);
//
//   int EvaluateODFluid(DRT::ELEMENTS::Transport*     ele,
//                        Teuchos::ParameterList&       params,
//                        DRT::Discretization &         discretization,
//                        const std::vector<int> &      lm,
//                        Epetra_SerialDenseMatrix&     elemat1_epetra,
//                        Epetra_SerialDenseMatrix&     elemat2_epetra,
//                        Epetra_SerialDenseVector&     elevec1_epetra,
//                        Epetra_SerialDenseVector&     elevec2_epetra,
//                        Epetra_SerialDenseVector&     elevec3_epetra);
//
//   //! calculate matrix and rhs. Here the whole thing is hidden.
//   virtual void SysmatODMesh(
//     DRT::Element*                         ele,       //!< the element we are dealing with
//     Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
//     const int                     numdofpernode
//   );
//
//   //! calculate matrix and rhs. Here the whole thing is hidden.
//   virtual void SysmatODFluid(
//     DRT::Element*                         ele,       //!< the element we are dealing with
//     Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
//     const int                     numdofpernode
//   );

  //! material ScaTra
  virtual void MatScaTra(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  double GetPorosityAtGP(const int  iquad);

  void SetDiffusivity(const Teuchos::RCP<const MAT::ScatraMat>& material,
                     const int                                 k,
                     Teuchos::RCP<ScaTraEleDiffManager>        diffmanager,
                     const double                              scale);

  void SetReaCoefficient(const Teuchos::RCP<const MAT::ScatraMat>& material,
                        const int                                 k,
                        Teuchos::RCP<ScaTraEleReaManager>         reamanager,
                        const double                              scale);

  void SetDensities(double  porosity,
                   double& densn,    //!< density at t_(n)
                   double& densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
                   double& densam   //!< density at t_(n+alpha_M));
                   );

};
}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */
