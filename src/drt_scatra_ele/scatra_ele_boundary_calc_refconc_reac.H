/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_refconc_reac.H

\brief main file containing routines for calculation of scatra element formulated in reference
concentrations and with advanced reaction terms

\level 3

 <pre>
   \maintainer Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_H
#define SCATRA_ELE_BOUNDARY_CALC_REFCONC_REAC_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleBoundaryCalcRefConcReac : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

     public:
      /// Singleton access method
      static ScaTraEleBoundaryCalcRefConcReac<distype>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname,
          const ScaTraEleBoundaryCalcRefConcReac* delete_me = NULL);

      /// Empty destructor
      virtual ~ScaTraEleBoundaryCalcRefConcReac(){};

      /// called upon destruction
      void Done();

     protected:
      //! Factor needed for the calculation of reference concentrations
      virtual double FacForRefConc(const int iquad,  ///< current boundary integration point
          const DRT::FaceElement* bele,              ///< current boundary element
          Teuchos::ParameterList& params,            ///< parameter list
          DRT::Discretization& discretization        ///< discretization
      );

     private:
      /// private constructor since we are singleton
      ScaTraEleBoundaryCalcRefConcReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
      double CalcJatIntPoint(const int iquad,  ///< current boundary integration point
          const DRT::FaceElement* bele,        ///< current boundary element
          const DRT::Element* pele,            ///< current parent element
          Teuchos::ParameterList& params,      ///< parameter list
          DRT::Discretization& discretization  ///< discretization
      );

    };  // class ScaTraEleBoundaryCalcRefConcReac

  }  // namespace ELEMENTS

}  // namespace DRT
#endif
