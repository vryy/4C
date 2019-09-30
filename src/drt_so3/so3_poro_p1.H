/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach)

 \level 2

\maintainer  Johannes Kremheller
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

      typedef So3_Poro<so3_ele, distype> my;

     public:
      //@}
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

      /*!
      \brief Destructor

      */
      virtual ~So3_Poro_P1(){};

      //@}

      //! number of dofs per node
      static const int noddof_ = my::noddof_ + 1;

      //! total dofs per element
      static const int numdof_ = noddof_ * my::numnod_;

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const;

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Unpack(const std::vector<char>& data);

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volumes of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      //! @name Access methods

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const { return 4; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const;

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
      );

      //@}

      //! initialize the inverse of the jacobian and its determinant in the material configuration
      virtual void InitElement();

      virtual void PreEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       ///< location array for de-assembly
      );

     protected:
      //! don't want = operator
      So3_Poro_P1& operator=(const So3_Poro_P1& old);

      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int MyEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
      );

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      void nlnstiff_poroelast(std::vector<int>& lm,        ///< location matrix
          LINALG::Matrix<my::numdim_, my::numnod_>& disp,  // current displacements
          LINALG::Matrix<my::numdim_, my::numnod_>& vel,   // current velocities
          LINALG::Matrix<my::numnod_, 1>* porosity,
          //    std::vector<double>&      residual,       ///< current residual displ
          LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,  //< fluid velocity of element
          LINALG::Matrix<my::numnod_, 1>& epreaf,            //< fluid pressure of element
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix,     ///< element stiffness matrix
          LINALG::Matrix<numdof_, numdof_>* reamatrix,       // element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,                 ///< element internal force vector
          Teuchos::ParameterList& params                     ///< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems
      void coupling_poroelast(std::vector<int>& lm,        ///< location matrix
          LINALG::Matrix<my::numdim_, my::numnod_>& disp,  // current displacements
          LINALG::Matrix<my::numdim_, my::numnod_>& vel,   // current velocities
          LINALG::Matrix<my::numnod_, 1>* porosity,
          LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,  //< fluid velocity of element
          LINALG::Matrix<my::numnod_, 1>& epreaf,            //< fluid pressure of element
          LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>*
              stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>*
              reamatrix,                      //< element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,  ///< element internal force vector
          Teuchos::ParameterList& params);    ///< algorithmic parameters e.g. time

      //! compute porosity at gausspoint and linearization of porosity w.r.t. structural
      //! displacements
      virtual void ComputePorosityAndLinearization(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<my::numnod_, 1>& shapfct,
          const LINALG::Matrix<my::numnod_, 1>* myporosity,
          const LINALG::Matrix<1, my::numdof_>& dJ_dus, double& porosity,
          LINALG::Matrix<1, my::numdof_>& dphi_dus);

      //! compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      virtual void ComputePorosityAndLinearizationOD(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<my::numnod_, 1>& shapfct,
          const LINALG::Matrix<my::numnod_, 1>* myporosity, double& porosity, double& dphi_dp);

      //! gauss point loop for evaluation of stiffness and rhs vector
      void GaussPointLoopP1(Teuchos::ParameterList& params,
          const LINALG::Matrix<my::numdim_, my::numnod_>& xrefe,
          const LINALG::Matrix<my::numdim_, my::numnod_>& xcurr,
          const LINALG::Matrix<my::numdim_, my::numnod_>& nodaldisp,
          const LINALG::Matrix<my::numdim_, my::numnod_>& nodalvel,
          const LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,
          const LINALG::Matrix<my::numnod_, 1>& epreaf,
          const LINALG::Matrix<my::numnod_, 1>* porosity_dof,
          LINALG::Matrix<my::numdof_, my::numdof_>& erea_v,
          LINALG::Matrix<my::numdof_, my::numdof_>* sub_stiff,
          LINALG::Matrix<my::numdof_, 1>* sub_force,
          LINALG::Matrix<my::numdof_, my::numnod_>& ecoupl_p1,
          LINALG::Matrix<my::numnod_, numdof_>& estiff_p1,
          LINALG::Matrix<my::numnod_, 1>& ecoupl_force_p1);

      //! gauss point loop for evaluation of stiffness (off diagonal)
      void GaussPointLoopP1OD(Teuchos::ParameterList& params,
          const LINALG::Matrix<my::numdim_, my::numnod_>& xrefe,
          const LINALG::Matrix<my::numdim_, my::numnod_>& xcurr,
          const LINALG::Matrix<my::numdim_, my::numnod_>& nodaldisp,
          const LINALG::Matrix<my::numdim_, my::numnod_>& nodalvel,
          const LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,
          const LINALG::Matrix<my::numnod_, 1>& epreaf,
          const LINALG::Matrix<my::numnod_, 1>* porosity_dof,
          LINALG::Matrix<my::numnod_, my::numnod_>& estiff_p1,
          LINALG::Matrix<my::numdof_, (my::numdim_ + 1) * my::numnod_>* sub_stiff);

      Teuchos::RCP<LINALG::Matrix<my::numnod_, 1>> init_porosity_;

      bool is_init_porosity_;
    };  // class So3_Poro_P1


  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* SO3_PORO_P1_H_ */
