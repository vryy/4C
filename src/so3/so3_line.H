/*----------------------------------------------------------------------*/
/*! \file

\brief line element

\level 1

*----------------------------------------------------------------------*/
#ifndef SO3_LINE_H
#define SO3_LINE_H


#include <Teuchos_SerialDenseMatrix.hpp>
#include "lib_node.H"
#include "lib_element.H"
#include "lib_elementtype.H"
#include "discretization_fem_general_utils_integration.H"



namespace DRT
{
  namespace ELEMENTS
  {
    class StructuralLineType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "StructuralLineType"; }

      static StructuralLineType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Teuchos::SerialDenseMatrix<int, double> nullspace;
        dserror("method ComputeNullSpace not implemented!");
        return nullspace;
      }

     private:
      static StructuralLineType instance_;
    };

    /*!
    \brief An element representing a line edge of any 3D structural element
           Note that this element should not be used in 2D cases!

    */
    class StructuralLine : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this line
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent shell element of this line
      \param lline: the local line number of this line w.r.t. the parent element
      */
      StructuralLine(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::Element* parent, const int lline);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      StructuralLine(const StructuralLine& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      inline int UniqueParObjectId() const override
      {
        return StructuralLineType::Instance().UniqueParObjectId();
      }

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
      \brief Destructor

      */
      ~StructuralLine() override = default;

      //@}

      //! @name Acess methods

      /*!
      \brief Get shape type of element
      */
      DiscretizationType Shape() const override;

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      inline int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      inline int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return StructuralLineType::Instance(); }

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the shell element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     private:
      // don't want = operator
      StructuralLine& operator=(const StructuralLine& old);

      //! gaussian integration to be used
      CORE::DRT::UTILS::GaussRule1D gaussrule_;

      //! line integration
      void LineIntegration(double& dL, const CORE::LINALG::SerialDenseMatrix& x,
          const CORE::LINALG::SerialDenseMatrix& deriv);

      /*!
      \brief Create matrix with material configuration

      \param x  (out)  : nodal coords in material frame
      */
      inline void MaterialConfiguration(CORE::LINALG::SerialDenseMatrix& x) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0];
          x(i, 1) = Nodes()[i]->X()[1];
          x(i, 2) = Nodes()[i]->X()[2];
        }
        return;
      }

    };  // class StructuralLine


  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SO3_LINE_H
