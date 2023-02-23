/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_cardiac_monodomain.H

\level 2

 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_H
#define SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_H

#include "scatra_ele_calc.H"
#include "scatra_ele_calc_aniso.H"
#include "scatra_ele_calc_advanced_reaction.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype, int probdim>
    class ScaTraEleCalcCardiacMonodomain : public ScaTraEleCalcAniso<distype, probdim>,
                                           public ScaTraEleCalcAdvReac<distype, probdim>
    {
     protected:
     private:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcCardiacMonodomain(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype, probdim> my;
      typedef ScaTraEleCalcAniso<distype, probdim> aniso;
      typedef ScaTraEleCalcAdvReac<distype, probdim> advreac;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      virtual ~ScaTraEleCalcCardiacMonodomain() {}


      /// Singleton access method
      static ScaTraEleCalcCardiacMonodomain<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate the element
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! extract element based or nodal values
      virtual void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! evaluate material
      virtual void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! material ScaTra
      void MatMyocard(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! calculate matrix and rhs for ep
      virtual void Sysmat(DRT::Element* ele,   ///< the element whose matrix is calculated
          Epetra_SerialDenseMatrix& emat,      ///< element matrix to calculate
          Epetra_SerialDenseVector& erhs,      ///< element rhs to calculate
          Epetra_SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
      );
    };

  }  // namespace ELEMENTS

}  // namespace DRT

#endif  // SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_H
