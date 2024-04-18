/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach.

\level 2


*/
/*---------------------------------------------------------------------------*/


#ifndef FOUR_C_W1_PORO_P1_HPP
#define FOUR_C_W1_PORO_P1_HPP

#include "baci_config.hpp"

#include "baci_w1_poro.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media

    */
    template <CORE::FE::CellType distype>
    class Wall1PoroP1 : public Wall1Poro<distype>
    {
      using Base = DRT::ELEMENTS::Wall1Poro<distype>;

     public:
      //!@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1PoroP1(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1PoroP1(const Wall1PoroP1& old);

      //!@}

      //! number of dofs per node
      static constexpr int noddof_ = Base::noddof_ + 1;

      //! total dofs per element
      static constexpr int numdof_ = noddof_ * Base::numnod_;

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override;

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      //! Get vector of Teuchos::RCPs to the lines of this element
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      //! Get vector of Teuchos::RCPs to the surfaces of this element
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //! @name Access methods

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override;

      //!@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate element stiffness, mass, internal forces, etc.

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  //!< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  //!< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  //!< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   //!< vector to be filled by element
          ) override;

      //! Evaluate a Neumann boundary condition
      //!
      //! this method evaluates a surfaces Neumann condition on the wall element
      //!
      //! \return 0 if successful, negative otherwise
      int EvaluateNeumann(
          Teuchos::ParameterList& params,  //!< (in/out) ParameterList for communication between
                                           //!< control routine and elements
          DRT::Discretization& discretization,  //!< A reference to the underlying discretization
          DRT::Condition& condition,            //!<  The condition to be evaluated
          std::vector<int>& lm,                 //!< location vector of this element
          CORE::LINALG::SerialDenseVector& elevec1,  //!< vector to be filled by element.
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      //!@}

      //! initialize the inverse of the jacobian and its determinant in the material configuration
      void InitElement() override;

      //! don't want = operator
      Wall1PoroP1& operator=(const Wall1PoroP1& old) = delete;

     protected:
      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int MyEvaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  //!< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  //!< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  //!< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   //!< vector to be filled by element
          ) override;

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      void NonlinearStiffnessPoroelast(std::vector<int>& lm,         //!< location matrix
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,  //!< current displacements
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,   //!< current velocities
          CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,      //!< porosity at element nodes
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>&
              evelnp,                                           //!< fluid velocity of element
          CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,       //!< fluid pressure of element
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<numdof_, numdof_>* reamatrix,    //!< element reactive matrix
          CORE::LINALG::Matrix<numdof_, 1>* force,              //!< element internal force vector
          Teuchos::ParameterList& params  //!< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems
      void CouplingPoroelast(std::vector<int>& lm,                   //!< location matrix
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,  //!< current displacements
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,   //!< current velocities
          CORE::LINALG::Matrix<Base::numnod_, 1>* porosity,          //!< porosity value
          CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>&
              evelnp,                                      //!< fluid velocity of element
          CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,  //!< fluid pressure of element
          CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
              stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
              reamatrix,                            //!< element reactive matrix
          CORE::LINALG::Matrix<numdof_, 1>* force,  //!< element internal force vector
          Teuchos::ParameterList& params);          //!< algorithmic parameters e.g. time

      //! compute porosity at gausspoint and linearization of porosity w.r.t. structural
      //! displacements
      void ComputePorosityAndLinearization(Teuchos::ParameterList& params, const double& press,
          const double& J, const int& gp, const CORE::LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const CORE::LINALG::Matrix<Base::numnod_, 1>* myporosity,
          const CORE::LINALG::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
          CORE::LINALG::Matrix<1, Base::numdof_>& dphi_dus) override;

      //! compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      void ComputePorosityAndLinearizationOD(Teuchos::ParameterList& params, const double& press,
          const double& J, const int& gp, const CORE::LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const CORE::LINALG::Matrix<Base::numnod_, 1>* myporosity, double& porosity,
          double& dphi_dp) override;

      //! evaluate gauss points for diagonal terms
      void GaussPointLoopP1(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
          const CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
          const CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
          CORE::LINALG::Matrix<Base::numdof_, Base::numdof_>& erea_v,
          CORE::LINALG::Matrix<Base::numdof_, Base::numdof_>* sub_stiff,
          CORE::LINALG::Matrix<Base::numdof_, 1>* sub_force,
          CORE::LINALG::Matrix<Base::numdof_, Base::numnod_>& ecoupl_p1,
          CORE::LINALG::Matrix<Base::numnod_, numdof_>& estiff_p1,
          CORE::LINALG::Matrix<Base::numnod_, 1>& ecoupl_force_p1);

      //! evaluate gauss points for off-diagonal terms
      void GaussPointLoopP1OD(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
          const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
          const CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
          const CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
          CORE::LINALG::Matrix<Base::numnod_, Base::numnod_>& ecoupl_p1,
          CORE::LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>& sub_stiff);
    };
  }  // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
