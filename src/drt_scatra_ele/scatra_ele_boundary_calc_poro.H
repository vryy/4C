/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_poro.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_PORO_H
#define SCATRA_ELE_BOUNDARY_CALC_PORO_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcPoro : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcPoro<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcPoro* delete_me = NULL);

      /// Empty destructor
      virtual ~ScaTraEleBoundaryCalcPoro(){};

      /// called upon destruction
      void Done();

      //! evaluate action
      int EvaluateAction(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, SCATRA::BoundaryAction action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      //! compute integral of convective mass/heat flux over boundary surface
      virtual std::vector<double> CalcConvectiveFlux(const DRT::FaceElement* ele,
          const std::vector<LINALG::Matrix<my::nen_, 1>>& ephinp,
          const LINALG::Matrix<my::nsd_ + 1, my::nen_>& evelnp, Epetra_SerialDenseVector& erhs);

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual double ComputePorosity(
          const DRT::FaceElement* ele  //!< the element we are dealing with
      );

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcPoro(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! nodal porosity values at t_(n+1)
      LINALG::Matrix<my::nen_, 1> eporosity_;

      //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
      LINALG::Matrix<my::nen_, 1> eprenp_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;

    };  // class ScaTraEleBoundaryCalcPoro
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
