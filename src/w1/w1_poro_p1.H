/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach.

\level 2


*/
/*---------------------------------------------------------------------------*/


#ifndef WALL1_PORO_P1_H_
#define WALL1_PORO_P1_H_

#include "w1_poro.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media

    */
    template <DRT::Element::DiscretizationType distype>
    class Wall1_PoroP1 : public Wall1_Poro<distype>
    {
      using Base = DRT::ELEMENTS::Wall1_Poro<distype>;

     public:
      //!@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1_PoroP1(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1_PoroP1(const Wall1_PoroP1& old);

      //!@}

      //! number of dofs per node
      static const int noddof_ = Base::noddof_ + 1;

      //! total dofs per element
      static const int numdof_ = noddof_ * Base::numnod_;

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
      void Pack(DRT::PackBuffer& data) const override;

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
          Epetra_SerialDenseMatrix& elemat1,    //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     //!< vector to be filled by element
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
          Epetra_SerialDenseVector& elevec1,    //!< vector to be filled by element.
          Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

      //!@}

      //! initialize the inverse of the jacobian and its determinant in the material configuration
      virtual void InitElement();

      //! don't want = operator
      Wall1_PoroP1& operator=(const Wall1_PoroP1& old) = delete;

     protected:
      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int MyEvaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     //!< vector to be filled by element
      );

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      void NonlinearStiffnessPoroelast(std::vector<int>& lm,     //!< location matrix
          LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,    //!< current displacements
          LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,     //!< current velocities
          LINALG::Matrix<Base::numnod_, 1>* porosity_dof,        //!< porosity at element nodes
          LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,  //!< fluid velocity of element
          LINALG::Matrix<Base::numnod_, 1>& epreaf,              //!< fluid pressure of element
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix,         //!< element stiffness matrix
          LINALG::Matrix<numdof_, numdof_>* reamatrix,           //!< element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,                     //!< element internal force vector
          Teuchos::ParameterList& params  //!< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems
      void CouplingPoroelast(std::vector<int>& lm,               //!< location matrix
          LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,    //!< current displacements
          LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,     //!< current velocities
          LINALG::Matrix<Base::numnod_, 1>* porosity,            //!< porosity value
          LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,  //!< fluid velocity of element
          LINALG::Matrix<Base::numnod_, 1>& epreaf,              //!< fluid pressure of element
          LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
              stiffmatrix,  //!< element stiffness matrix
          LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
              reamatrix,                      //!< element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,  //!< element internal force vector
          Teuchos::ParameterList& params);    //!< algorithmic parameters e.g. time

      //! compute porosity at gausspoint and linearization of porosity w.r.t. structural
      //! displacements
      virtual void ComputePorosityAndLinearization(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const LINALG::Matrix<Base::numnod_, 1>* myporosity,
          const LINALG::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
          LINALG::Matrix<1, Base::numdof_>& dphi_dus);

      //! compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      virtual void ComputePorosityAndLinearizationOD(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const LINALG::Matrix<Base::numnod_, 1>* myporosity, double& porosity, double& dphi_dp);

      //! evaluate gauss points for diagonal terms
      void GaussPointLoopP1(Teuchos::ParameterList& params,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
          const LINALG::Matrix<Base::numnod_, 1>& epreaf,
          const LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
          LINALG::Matrix<Base::numdof_, Base::numdof_>& erea_v,
          LINALG::Matrix<Base::numdof_, Base::numdof_>* sub_stiff,
          LINALG::Matrix<Base::numdof_, 1>* sub_force,
          LINALG::Matrix<Base::numdof_, Base::numnod_>& ecoupl_p1,
          LINALG::Matrix<Base::numnod_, numdof_>& estiff_p1,
          LINALG::Matrix<Base::numnod_, 1>& ecoupl_force_p1);

      //! evaluate gauss points for off-diagonal terms
      void GaussPointLoopP1OD(Teuchos::ParameterList& params,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
          const LINALG::Matrix<Base::numnod_, 1>& epreaf,
          const LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
          LINALG::Matrix<Base::numnod_, Base::numnod_>& estiff_p1,
          LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>& sub_stiff);
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* WALL1_PORO_P1_H_ */
