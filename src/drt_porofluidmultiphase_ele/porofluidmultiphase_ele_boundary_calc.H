/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_boundary_calc.H

 \brief implementation of evaluation routines of porofluidmultiphase boundary element

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_H_


#include "porofluidmultiphase_ele_action.H"
#include "porofluidmultiphase_ele_boundary_interface.H"
#include "porofluidmultiphase_ele_calc_utils.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

template<DRT::Element::DiscretizationType distype>
class PoroFluidMultiPhaseEleBoundaryCalc: public PoroFluidMultiPhaseEleBoundaryInterface
{
public:

  /// Singleton access method
  static PoroFluidMultiPhaseEleBoundaryCalc<distype> * Instance(
      const int numdofpernode,
      const std::string& disname,
      const PoroFluidMultiPhaseEleBoundaryCalc* delete_me = NULL );

  /// Empty destructor
  virtual ~PoroFluidMultiPhaseEleBoundaryCalc(){return;};

  /// called upon destruction
  virtual void Done();

  //! number of element nodes (nomenclature: T. Hughes, The finite element method)
  static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of boundary(!) space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! setup element evaluation
  int SetupCalc( DRT::FaceElement*                 ele,
                 Teuchos::ParameterList&           params,
                 DRT::Discretization&              discretization);

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual void ExtractElementAndNodeValues(
    DRT::FaceElement*             ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
  );

  //! Evaluate the element (using location array)
  virtual int Evaluate(
      DRT::FaceElement*             ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      Epetra_SerialDenseMatrix&     elemat1_epetra,
      Epetra_SerialDenseMatrix&     elemat2_epetra,
      Epetra_SerialDenseVector&     elevec1_epetra,
      Epetra_SerialDenseVector&     elevec2_epetra,
      Epetra_SerialDenseVector&     elevec3_epetra
      );

  //! evaluate action
  virtual int EvaluateAction(
      DRT::FaceElement*                   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      POROFLUIDMULTIPHASE::BoundaryAction action,
      DRT::Element::LocationArray&        la,
      Epetra_SerialDenseMatrix&           elemat1_epetra,
      Epetra_SerialDenseMatrix&           elemat2_epetra,
      Epetra_SerialDenseVector&           elevec1_epetra,
      Epetra_SerialDenseVector&           elevec2_epetra,
      Epetra_SerialDenseVector&           elevec3_epetra
      );

  //! evaluate Neumann boundary condition
  virtual int EvaluateNeumann(
      DRT::FaceElement*                   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      DRT::Condition&                     condition,
      DRT::Element::LocationArray&        la,
      Epetra_SerialDenseVector&           elevec1,
      const double                        scalar
      );

protected:
  //! evaluate shape functions and derivatives at int. point
  double EvalShapeFuncAndIntFac(
      const DRT::UTILS::IntPointsAndWeights<nsd_>&   intpoints,         ///< integration points
      const int                                      iquad,             ///< id of current Gauss point
      LINALG::Matrix<1 + nsd_,1>*                    normalvec = NULL   ///< normal vector at Gauss point(optional)
  );

private:
  /// private constructor since we are singleton
  PoroFluidMultiPhaseEleBoundaryCalc(const int numdofpernode, const std::string& disname);

  //! pointer to parameter list
  DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* params_;

   //! number of dof per node
  const int numdofpernode_;

 //! node coordinates
 LINALG::Matrix<nsd_+1,nen_> xyze_;
 //! nodal displacement values for ALE
 LINALG::Matrix<nsd_+1,nen_> edispnp_;
 //! coordinates of current integration point in reference coordinates
 LINALG::Matrix<nsd_,1> xsi_;
 //! array for shape functions
 LINALG::Matrix<nen_,1> funct_;
 //! array for shape function derivatives w.r.t r,s,t
 LINALG::Matrix<nsd_,nen_> deriv_;
 //! global derivatives of shape functions w.r.t x,y,z
 LINALG::Matrix<nsd_,nen_> derxy_;
 //! unit normal vector at integration point
 LINALG::Matrix<nsd_+1,1> normal_;
 //! velocity vector in gausspoint
 LINALG::Matrix<nsd_+1,1> velint_;
 //! metric tensor at integration point
 LINALG::Matrix<nsd_,nsd_>  metrictensor_;

};
}
}



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_CALC_H_ */
