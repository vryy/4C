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

//forward declarations
namespace MAT
{
class ScatraMat;
}

namespace DRT
{
namespace ELEMENTS
{

class ScaTraEleDiffManagerPoro;


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

  //! read element coordinates, assuming they are all 3D and then project to the respective lower dimensional space
  virtual void ReadElementCoordinatesAndProject(
      const DRT::ELEMENTS::Transport*     ele
  );

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual const std::vector<double> ExtractElementAndNodeValues(
    DRT::ELEMENTS::Transport*  ele,
    Teuchos::ParameterList&    params,
    DRT::Discretization&       discretization,
    const std::vector<int>&    lm
  );

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );

  //! compute porosity based on solid, fluid and (potentially) scatra solution
  virtual void ComputePorosity(
      const DRT::Element* ele       //!< the element we are dealing with
    );

  //! material ScaTra
  virtual void MatScaTra(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );


  //! set diffusivity for poro scatra problem (i.e. scale by porosity)
  void SetDiffusivity(const Teuchos::RCP<const MAT::ScatraMat>& material,
                     const int                                 k,
                     const double                              scale);

  //! set reaction coefficient for poro scatra problem (i.e. scale by porosity)
  void SetReaCoefficient(const Teuchos::RCP<const MAT::ScatraMat>& material,
                        const int                                 k,
                        const double                              scale);

  //! set densisties for poro scatra problem (i.e. scale by porosity)
  void SetDensities(double  porosity,
                   double& densn,    //!< density at t_(n)
                   double& densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
                   double& densam   //!< density at t_(n+alpha_M));
                   );

  //! get poro diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerPoro> DiffManager(){return Teuchos::rcp_static_cast<ScaTraEleDiffManagerPoro>(my::diffmanager_);};

  /*========================================================================*/
  //! @name Galerkin approximation and related
  /*========================================================================*/

  //! initial node coordinates
  LINALG::Matrix<my::nsd_,my::nen_> xyze0_;

};

/// ScaTraEleDiffManagerPoro implementation
/*!
  This class keeps all poro-specific transport parameter needed for the evaluation of an element. The
  ScaTraEleDiffManagerPoro is derived from the standard ScaTraEleDiffManager.
*/
class ScaTraEleDiffManagerPoro : public ScaTraEleDiffManager
{
  public:

    ScaTraEleDiffManagerPoro(int numscal):
    ScaTraEleDiffManager(numscal),
    porosity_(0.0)
    {return;}

    ~ScaTraEleDiffManagerPoro(){};

    void SetPorosity(double porosity)
    {
      porosity_=porosity;
    }

    double GetPorosity() const
    {
      return porosity_;
    }

  protected:

    //! porosity at gauss point
    double porosity_;
};

}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */
