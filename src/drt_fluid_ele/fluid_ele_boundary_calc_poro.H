/*----------------------------------------------------------------------*/
/*! \file

\brief evaluate boundary conditions for poroelast / fpsi

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_BOUNDARY_CALC_PORO_H
#define FLUID_ELE_BOUNDARY_CALC_PORO_H

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "fluid_ele_boundary_calc.H"

/*---------------------------------------------------------------------*
 | forward declarations                                                 |
 *---------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;

  namespace ELEMENTS
  {
    class FluidBoundary;

    /// Class for Evaluating boundary integrals for porous media problems
    /*!
     This class is derived from the FluidBoundaryImpl class, i.e. it is capable of
     evaluated all integrals implemented there. It will do so, if the evaluate action
     given by the control routine is not known (see method EvaluateAction).
     Otherwise, it can evaluate integrals for special poro boundary conditions (as
     no penetration constraint terms or pressure coupling) or overwrite existing methods
     that need to be reimplemented for porous flow (like flow rate calculation).

     This a calculation class implemented as a singleton, like all calc classes in fluid
     (see comments on base classes for more details). In short this means that on instance
     exists for every discretization type of the boundary element (because of the template).

     Also for those cases where the boundary element needs to assemble into dofs of its parent
     element, the corresponding methods are implemented twice. The first just contains a
     switch identifying the discretization type of the parent element and the second methods
     is templated by the discretization type of the parent element (and also by the boundary
     distype) and does the actual work. Note that all those methods need the location vector
     of the parent(!) element to fill the matrixes. Therefore the evaluate actions
     corresponding to these boundaries need to be listed in the FluidPoroBoundary::LocationVector()
     method. Otherwise one might get an dserror or even no error at all and very strange results
     instead.

     \author vuong 10/14
     */

    template <DRT::Element::DiscretizationType distype>
    class FluidEleBoundaryCalcPoro : public FluidBoundaryImpl<distype>
    {
      typedef DRT::ELEMENTS::FluidBoundaryImpl<distype> my;

     public:
      /// Singleton access method
      static FluidEleBoundaryCalcPoro<distype>* Instance(bool create = true);

      /// Empty destructor
      virtual ~FluidEleBoundaryCalcPoro(){};

      /// called upon destruction
      virtual void Done();

      /// determines which boundary integral is to be evaluated
      void EvaluateAction(DRT::ELEMENTS::FluidBoundary* ele1, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3);

     protected:
      /// protected constructor since we are singleton
      FluidEleBoundaryCalcPoro();

      /*!
      \brief applies boundary integral for porous media problems
             (contains switch and calls parent distype templated version)

             This method evaluates the boundary integral appearing when integrating
             the continuity equation by parts for a porous fluid.
             I.e. when the CONTIPARTINT flag is set to 'yes' in the POROELASTICITY DYNAMIC section,
             it will be evaluated on the PORO PARTIAL INTEGRATION condition

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
      void PoroBoundary(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply boundary integral (mass flux over boundary in continuity equation) for porous
      media problems (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void PoroBoundary(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

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
      void PressureCoupling(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply boundary coupling terms for FPSI problems
             (contains switch and calls parent distype templated version)

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
      void FPSICoupling(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply boundary coupling terms for FPSI problems
             (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void FPSICoupling(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply no penetration boundary condition
             This method applies no penetration boundary condition as a strong nodal constraint.
             It works for plain surfaces, but may lead to problems for curved boundaries!!

             Alternative (and probably better!!) ways to apply the no penetration condition
             are in a weak substitution sense (PoroBoundary method) or with Lagrange multiplier
             (implemented in poro_monolithicsplit_nopenetration.cpp and methods
      NoPenetrationMatAndRHS and NoPenetrationMatOD)

       Note: This method is called with one additional action added to the parameter list.
             It can either be POROELAST::fluidfluid or POROELAST::fluidstructure.

             In case of POROELAST::fluidfluid it fill elevec1 with the residual of the no
      penetration constraint and elemat1 with the linearization of the constraint w.r.t. fluid
      velocities. elemat2 will not be filled.

             In case of POROELAST::fluidstructure it fill elemat1 with the linearization of the
             constraint w.r.t. structure displacements and elemat2 with the linearization of the
             constraint w.r.t. structure velocities .
             elevec1 will not be filled.


      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element.
      \param elemat2 (out)      : matrix to be filled by element.
      \param elevec1 (out)      : vector to be filled by element.
      */
      void NoPenetration(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply no penetration boundary condition
             This method is to be called before the NoPenetration method to find the dof IDs
             that will be subject to the nodal no penetration constraint. I.e. find the
             normal direction of a potentially curved boundary.
             This a quick and dirty version. I would recommend to do it in another way (see
             comment to NoPenetration method)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. This vector is used as
                                  toggle vector for the constraint dofs, i.e. dofs
                                  where a constraint is to be applied will be marked with 1.0
                                  in this vector otherwise 0.0
      */
      void NoPenetrationIDs(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1,
          std::vector<int>& lm);

      /*!
      \brief compute flow rate over boundary for porous media problems
            (contains switch and calls parent distype templated version)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      */
      void ComputeFlowRate(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1);

      /*!
      \brief compute flow rate over boundary for porous media problems
             (templated by discretization type of parent element)

      */
      template <DRT::Element::DiscretizationType pdistype>
      void ComputeFlowRate(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply no penetration boundary condition
             This method applies no penetration boundary condition using Lagrange Multipliers.
             (see poro_monolithicsplit_nopenetration.cpp for control method)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element. Linearization of constraint w.r.t
      to fluid dofs \param elemat2 (out)      : matrix to be filled by element. D matrix (mortar
      nomenclature) \param elevec1 (out)      : vector to be filled by element. Residual of the
      constraint.
      */
      void NoPenetrationMatAndRHS(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply no penetration boundary condition
             (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void NoPenetrationMatAndRHS(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseVector& elevec1);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             This method applies no penetration boundary condition using Lagrange Multipliers.
             (see poro_monolithicsplit_nopenetration.cpp for control method)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element. Linearization of constraint w.r.t
      to structure dofs \param elemat2 (out)      : matrix to be filled by element. Linearization of
      constraint w.r.t to lagrange multiplier
      */
      void NoPenetrationMatOD(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void NoPenetrationMatOD(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             This method applies no penetration boundary condition using Lagrange Multipliers.
             (see poro_monolithicsplit_nopenetration.cpp for control method)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element. Linearization of constraint w.r.t
      pressure part ofporosity
      */
      void NoPenetrationMatODPoroPres(DRT::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void NoPenetrationMatODPoroPres(DRT::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             This method applies no penetration boundary condition using Lagrange Multipliers.
             (see poro_monolithicsplit_nopenetration.cpp for control method)

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element. Linearization of constraint w.r.t
      displacement part ofporosity
      */
      void NoPenetrationMatODPoroDisp(DRT::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1);

      /*!
      \brief apply no penetration boundary condition (Off Diagonal terms)
             (templated by discretization type of parent element)
      */
      template <DRT::Element::DiscretizationType pdistype>
      void NoPenetrationMatODPoroDisp(DRT::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1);


      /*!
      \brief compute the porosity at nodes

      This is a dummy class and will be reimplemented in FluidEleBoundaryCalcPoroP1

      \param ele (in)            : current element
      \param mydispnp (in)       : nodal displacements
      \param eporosity (in/out)  : matrix to be filled with nodal porosities
      */
      virtual bool ComputeNodalPorosity(DRT::ELEMENTS::FluidBoundary* ele,
          const std::vector<double>& mydispnp, LINALG::Matrix<my::bdrynen_, 1>& eporosity)
      {
        return false;
      };

      /*!
      \brief compute the porosity at guass point

      This class will access the poro structure material and evalute the porosity.
      reimplemented in FluidEleBoundaryCalcPoroP1 class.

      \param params (in)         : ParameterList for communication between control routine
      \param ele (in)            : current element
      \param funct (in)          : nodal displacements
      \param eporosity (in)      : matrix filled nodal porosities (only needed in
      FluidEleBoundaryCalcPoroP1) \param press (in)          : fluid pressure at gauss point \param
      J (in)              : Jacobian determinant of deformation gradient in gauss point \param gp
      (in)             : current gauss point number (needed for saving) \param porosity (in/out)   :
      porosity at gauss point \param dphi_dp (in/out)    : derivative of porosity w.r.t. fluid
      pressure \param dphi_dJ (in/out)    : derivative of porosity w.r.t. Jacobian determinant
      \param save (in)           : flag whether the porosity is to be saved for later access
      */
      virtual void ComputePorosityAtGP(Teuchos::ParameterList& params,
          DRT::ELEMENTS::FluidBoundary* ele, const LINALG::Matrix<my::bdrynen_, 1>& funct,
          const LINALG::Matrix<my::bdrynen_, 1>& eporosity, double press, double J, int gp,
          double& porosity, double& dphi_dp, double& dphi_dJ, bool save);

    };  // class FluidEleBoundaryCalcPoro

    /// Class for Evaluating boundary integrals for porous media problems (P1 approach)
    /*!
      This class is implements the poro boundary condition for poro P1 elements, i.e. with
     additional nodal porosity degree of freedom. Therefore the only difference should be the way
     the porosity is evaluated.

     \author vuong 10/14
     */
    template <DRT::Element::DiscretizationType distype>
    class FluidEleBoundaryCalcPoroP1 : public FluidEleBoundaryCalcPoro<distype>
    {
      typedef DRT::ELEMENTS::FluidEleBoundaryCalcPoro<distype> my;

     public:
      /// Singleton access method
      static FluidEleBoundaryCalcPoroP1<distype>* Instance(bool create = true);

      /// Empty destructor
      virtual ~FluidEleBoundaryCalcPoroP1(){};

      /// called upon destruction
      virtual void Done();

     protected:
      /*!
      \brief compute the porosity at nodes

      As we have nodal porosities in this case, they can be read from the nodal state vector

      \param ele (in)            : current element
      \param mydispnp (in)       : nodal displacements/porosities
      \param eporosity (in/out)  : matrix to be filled with nodal porosities
      */
      bool ComputeNodalPorosity(DRT::ELEMENTS::FluidBoundary* ele,
          const std::vector<double>& mydispnp, LINALG::Matrix<my::bdrynen_, 1>& eporosity);

      /*!
      \brief compute the porosity at nodes

      As we have porosities as primary variables here, they can be evaluated with the shape
      functions

      \param ele (in)            : current element
      \param mydispnp (in)       : nodal displacements/porosities
      \param eporosity (in/out)  : matrix to be filled with nodal porosities
      */
      void ComputePorosityAtGP(Teuchos::ParameterList& params, DRT::ELEMENTS::FluidBoundary* ele,
          const LINALG::Matrix<my::bdrynen_, 1>& funct,
          const LINALG::Matrix<my::bdrynen_, 1>& eporosity, double press, double J, int gp,
          double& porosity, double& dphi_dp, double& dphi_dJ, bool save);
    };  // class FluidEleBoundaryCalcPoroP1

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
