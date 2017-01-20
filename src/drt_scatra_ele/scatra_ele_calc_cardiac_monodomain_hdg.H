/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_cardiac_monodomain_hdg.H

\brief main file containing routines for calculation of HDG cardiac monodomain element

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_H
#define SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_H


#include "scatra_ele_calc_hdg.H"
#include "scatra_ele_hdg.H"


#include "../drt_fem_general/drt_utils_polynomial.H"

namespace DRT
{

namespace ELEMENTS
{

  /// Scatra HDG element implementation
  template<DRT::Element::DiscretizationType distype, int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
  class ScaTraEleCalcHDGCardiacMonodomain : public ScaTraEleCalcHDG<distype,probdim>
  {
  protected:

  private:

    /// (private) protected constructor, since we are a Singleton.
    /// this constructor is called from a derived class
    /// -> therefore, it has to be protected instead of private
    ScaTraEleCalcHDGCardiacMonodomain(const int numdofpernode,const int numscal,const std::string& disname);

//    typedef ScaTraEleCalc<distype,probdim> my;
//    typedef ScaTraEleCalcAniso<distype,probdim> aniso;
//    typedef ScaTraEleCalcAdvReac<distype,probdim> advreac;

    std::vector<Epetra_SerialDenseVector> values_mat_gp_all_;

  public:

    virtual ~ScaTraEleCalcHDGCardiacMonodomain() {}


    /// Singleton access method
    static ScaTraEleCalcHDGCardiacMonodomain<distype,probdim> * Instance(
      const int numdofpernode,
      const int numscal,
      const std::string& disname,
      bool create=true );

    /// called upon destruction
    virtual void Done();

    //! evaluate the element
    int EvaluateAction(DRT::Element*                 ele,
                       Teuchos::ParameterList&       params,
                       DRT::Discretization &         discretization,
                       const SCATRA::Action&         action,
                       DRT::Element::LocationArray&  la,
                       Epetra_SerialDenseMatrix&     elemat1_epetra,
                       Epetra_SerialDenseMatrix&     elemat2_epetra,
                       Epetra_SerialDenseVector&     elevec1_epetra,
                       Epetra_SerialDenseVector&     elevec2_epetra,
                       Epetra_SerialDenseVector&     elevec3_epetra) {return 0;};

  protected:

    /*========================================================================*/
    //! @name material and related and related functions
    /*========================================================================*/

    //! evaluate material
    virtual void PrepareMaterials(
      DRT::Element*                             ele,      //!< the element we are dealing with
      const Teuchos::RCP<const MAT::Material>   material, //!< pointer to current material
      const int                                 k,        //!< id of current scalar
      Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> >  difftensor
      );


    //! evaluate material
    virtual void Materials(
      const Teuchos::RCP<const MAT::Material>   material, //!< pointer to current material
      const int                                 k,        //!< id of current scalar
      Epetra_SerialDenseMatrix&                 difftensor,
      Epetra_SerialDenseVector&                 ivecn,
      Epetra_SerialDenseVector&                 ivecnp,
      Epetra_SerialDenseMatrix&                 ivecnpderiv
      );

    //! material ScaTra
    void MatMyocard(
      const Teuchos::RCP<const MAT::Material>   material, //!< pointer to current material
      const int                                 k,        //!< id of current scalar
      Epetra_SerialDenseMatrix&                 difftensor,
      Epetra_SerialDenseVector&                 ivecn,
      Epetra_SerialDenseVector&                 ivecnp,
      Epetra_SerialDenseMatrix&                 ivecnpderiv
      );

    //! update time dependent material
    void TimeUpdateMaterial(
      const DRT::Element*   ele   //!< the element we are dealing with
      );

    //! get material internal state for output
    void GetMaterialInternalState(
      const DRT::Element *              ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization &             discretization
      );

    //! set material internal state after restart
    void SetMaterialInternalState(
      const DRT::Element *              ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization &             discretization
      );

    //! project material field
    int ProjectMaterialField(
      const DRT::Element *              ele
      );

    /// polynomial space for element interior for various Gauss Points for the evaluation of the material
    Teuchos::RCP<DRT::UTILS::PolynomialSpace<probdim> > polySpace_;


    //!  Get Gauss Rule 2D
    DRT::UTILS::GaussRule2D GetGaussRule2D(
        const Teuchos::RCP<const MAT::Material>   material   //!< pointer to current material
        );

    //!  Get Gauss Rule 2D
    DRT::UTILS::GaussRule3D GetGaussRule3D(
        const Teuchos::RCP<const MAT::Material>   material   //!< pointer to current material
        );
  };

}
}

#endif /* SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_H_ */

