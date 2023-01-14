/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach)

 \level 2

 *----------------------------------------------------------------------*/


#ifndef SO3_PORO_P1_H_
#define SO3_PORO_P1_H_

#include "so3_poro.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 3 dimensional solid element with modifications for porous media

    A structural 3 dimensional solid displacement element for large deformations
    and (near)-incompressibility.

    */
    template <class so3_ele, DRT::Element::DiscretizationType distype>
    class So3_Poro_P1 : public So3_Poro<so3_ele, distype>
    {
      //! @name Friends
      friend class So_hex8PoroP1Type;

      using Base = So3_Poro<so3_ele, distype>;

     public:
      //!@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3_Poro_P1(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3_Poro_P1(const So3_Poro_P1& old);

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

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the volumes of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Volumes() override;

      //! @name Access methods

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override { return 4; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override;

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.

      If NULL on input, the controlling method does not expect the element
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

      //!@}

      //! initialize the inverse of the jacobian and its determinant in the material configuration
      void InitElement() override;

      void PreEvaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       //!< location array for de-assembly
          ) override;

      //! don't want = operator
      So3_Poro_P1& operator=(const So3_Poro_P1& old) = delete;

     protected:
      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int MyEvaluate(
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

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      void NonlinearStiffnessPoroelast(std::vector<int>& lm,     //!< location matrix
          LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,    //!< current displacements
          LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,     //!< current velocities
          LINALG::Matrix<Base::numnod_, 1>* porosity,            //!< porosity value
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
      void ComputePorosityAndLinearization(Teuchos::ParameterList& params, const double& press,
          const double& J, const int& gp, const LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const LINALG::Matrix<Base::numnod_, 1>* myporosity,
          const LINALG::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
          LINALG::Matrix<1, Base::numdof_>& dphi_dus) override;

      //! compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      void ComputePorosityAndLinearizationOD(Teuchos::ParameterList& params, const double& press,
          const double& J, const int& gp, const LINALG::Matrix<Base::numnod_, 1>& shapfct,
          const LINALG::Matrix<Base::numnod_, 1>* myporosity, double& porosity,
          double& dphi_dp) override;

      //! gauss point loop for evaluation of stiffness and rhs vector
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

      //! gauss point loop for evaluation of stiffness (off diagonal)
      void GaussPointLoopP1OD(Teuchos::ParameterList& params,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
          const LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
          const LINALG::Matrix<Base::numnod_, 1>& epreaf,
          const LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
          LINALG::Matrix<Base::numnod_, Base::numnod_>& estiff_p1,
          LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>* sub_stiff);

      Teuchos::RCP<LINALG::Matrix<Base::numnod_, 1>> init_porosity_;

      bool is_init_porosity_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* SO3_PORO_P1_H_ */
