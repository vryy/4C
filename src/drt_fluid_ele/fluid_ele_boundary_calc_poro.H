/*!----------------------------------------------------------------------
\file fluid_ele_boundary_calc_poro.H

\brief evaluate boundary conditions for poroelast / fpsi

<pre>
Maintainers: Anh-Tu Vuong & Andreas Rauch
             {vuong,rauch}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264/-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef FLUID_ELE_BOUNDARY_CALC_PORO_H
#define FLUID_ELE_BOUNDARY_CALC_PORO_H

#include "fluid_ele_boundary_calc.H"

namespace DRT
{
class Condition;
class Discretization;

namespace ELEMENTS
{

class FluidBoundary;

template<DRT::Element::DiscretizationType distype>
class FluidEleBoundaryCalcPoro : public FluidBoundaryImpl<distype>
{
  typedef DRT::ELEMENTS::FluidBoundaryImpl<distype> my;

public:

  /// Singleton access method
  static FluidEleBoundaryCalcPoro<distype> * Instance( bool create=true );

  /// Empty destructor
  virtual ~FluidEleBoundaryCalcPoro() {};

  /// called upon destruction
  virtual void Done();

  /// determines which boundary integral is to be evaluated
  void EvaluateAction(DRT::ELEMENTS::FluidBoundary*  ele1,
                      Teuchos::ParameterList&         params,
                      DRT::Discretization&            discretization,
                      std::vector<int>&               lm,
                      Epetra_SerialDenseMatrix&       elemat1,
                      Epetra_SerialDenseMatrix&       elemat2,
                      Epetra_SerialDenseVector&       elevec1,
                      Epetra_SerialDenseVector&       elevec2,
                      Epetra_SerialDenseVector&       elevec3);

protected:
  /// private constructor since we are singleton
  FluidEleBoundaryCalcPoro();

  /*!
  \brief apply boundary integral for porous media problems

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  void PoroBoundary(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&        params,
      DRT::Discretization&           discretization,
      std::vector<int>&              lm,
      Epetra_SerialDenseMatrix&      elemat1,
      Epetra_SerialDenseVector&      elevec1);

  /*!
  \brief apply boundary pressure for porous media problems

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  void PressureCoupling(
      DRT::ELEMENTS::FluidBoundary*   ele,
      Teuchos::ParameterList&         params,
      DRT::Discretization&            discretization,
      std::vector<int>&               lm,
      Epetra_SerialDenseMatrix&       elemat1,
      Epetra_SerialDenseVector&       elevec1);

  /*!
  \brief apply boundary coupling terms for FPSI problems

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  void FPSICoupling(
      DRT::ELEMENTS::FluidBoundary*   ele,
      Teuchos::ParameterList&         params,
      DRT::Discretization&            discretization,
      std::vector<int>&               lm,
      Epetra_SerialDenseMatrix&       elemat1,
      Epetra_SerialDenseVector&       elevec1);

  /*!
  \brief apply no penetration boundary condition

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elemat2 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  void NoPenetration(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseMatrix&  elemat2,
      Epetra_SerialDenseVector&  elevec1,
      Epetra_SerialDenseVector&  elevec2);

  void NoPenetrationIDs(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      Epetra_SerialDenseVector&        elevec1,
      std::vector<int>&          lm);

  /*!
  \brief apply no penetration boundary condition

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  void ComputeFlowRate(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseVector&  elevec1);

  /*!
  \brief apply boundary integral (mass flux over boundary in continuity equation) for porous media problems

  \param ele (in)           : Boundary element
  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  template <DRT::Element::DiscretizationType pdistype>
  void PoroBoundary(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

  template <DRT::Element::DiscretizationType pdistype>
  void FPSICoupling(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseVector&  elevec1);

  void NoPenetrationMatAndRHS(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseMatrix&  elemat2,
      Epetra_SerialDenseVector&  elevec1);

  void NoPenetrationMatOD(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseMatrix&  elemat1,
      Epetra_SerialDenseMatrix&  elemat2);

  /*!
  \brief compute flow rate over boundary for porous media problems

  \param ele (in)           : Boundary element
  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
  */
  template <DRT::Element::DiscretizationType pdistype>
  void ComputeFlowRate(
      DRT::ELEMENTS::FluidBoundary*  ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      std::vector<int>&          lm,
      Epetra_SerialDenseVector&  elevec1);

  virtual bool ComputeNodalPorosity(DRT::ELEMENTS::FluidBoundary*  ele,
                                   const std::vector<double> &     mydispnp,
                                   LINALG::Matrix<my::bdrynen_,1>& eporosity )
  {return false;};

  virtual void ComputePorosityAtGP( Teuchos::ParameterList& params,
                                    DRT::ELEMENTS::FluidBoundary*  ele,
                                    const LINALG::Matrix<my::bdrynen_,1>& funct,
                                    const LINALG::Matrix<my::bdrynen_,1>& eporosity,
                                    double press,
                                    double J,
                                    int gp,
                                    double& porosity,
                                    double& dphi_dp,
                                    double& dphi_dJ );

}; // class FluidEleBoundaryCalcPoro

template<DRT::Element::DiscretizationType distype>
class FluidEleBoundaryCalcPoroP1 : public FluidEleBoundaryCalcPoro<distype>
{
  typedef DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype> my;

public:

  /// Singleton access method
  static FluidEleBoundaryCalcPoroP1<distype> * Instance( bool create=true );

  /// Empty destructor
  virtual ~FluidEleBoundaryCalcPoroP1() {};

  /// called upon destruction
  virtual void Done();

protected:
  bool ComputeNodalPorosity( DRT::ELEMENTS::FluidBoundary*  ele,
                             const std::vector<double> &     mydispnp,
                             LINALG::Matrix<my::bdrynen_,1>& eporosity);

  void ComputePorosityAtGP( Teuchos::ParameterList& params,
                                    DRT::ELEMENTS::FluidBoundary*  ele,
                                    const LINALG::Matrix<my::bdrynen_,1>& funct,
                                    const LINALG::Matrix<my::bdrynen_,1>& eporosity,
                                    double press,
                                    double J,
                                    int gp,
                                    double& porosity,
                                    double& dphi_dp,
                                    double& dphi_dJ );
};// class FluidEleBoundaryCalcPoroP1

} // namespace ELEMENTS
} // namespace DRT

#endif
