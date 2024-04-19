/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_BELE_VELE3_HPP
#define FOUR_C_BELE_VELE3_HPP


#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;


  namespace ELEMENTS
  {
    class Vele3Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Vele3Type"; }

      static Vele3Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Vele3Type instance_;
    };

    /*!
     * A 3D volume element with no physics attached
     *
     * This element is meant to have no physics. It can be used to have a volume discretization
     * of volume elements. They can be of any 3d shape
     *
     * The number of dof per node is set to 3, so we can define displacement vectors by
     * using FillComplete on the discretization.
     *
     */


    /*!
    \brief A register for bele3 element

    */
    class Vele3 : public DRT::Element
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      explicit Vele3(int id,  ///< A unique global id
          int owner           ///< proc num that owns this element
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit Vele3(const Vele3& old);


      /*!
      \brief Deep copy this instance of Fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override
      {
        if (NumNode() == 8 || NumNode() == 20 || NumNode() == 27)
          return 12;
        else if (NumNode() == 4 || NumNode() == 10)
          return 6;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override
      {
        if (NumNode() == 8 || NumNode() == 20 || NumNode() == 27)
          return 6;
        else if (NumNode() == 4 || NumNode() == 10)
          return 4;
        else if (NumNode() == 6 || NumNode() == 15 || NumNode() == 5)
          return 5;
        else
        {
          FOUR_C_THROW("Could not determine number of surfaces");
          return -1;
        }
      }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int NumVolume() const override { return 1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return Vele3Type::Instance().UniqueParObjectId(); }

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



      //@}

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
    \brief Get number of degrees of freedom per element
           (implements pure virtual DRT::Element)

    The element decides how many element degrees of freedom it has.
    It can redecide along the way of a simulation.

    \note Element degrees of freedom mentioned here are dofs that are visible
          at the level of the total system of equations. Purely internal
          element dofs that are condensed internally should NOT be considered.
    */
      int NumDofPerElement() const override { return 0; }

      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return Vele3Type::Instance(); }


      //@}

      //@}

      //! @name Evaluation
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

     private:
      //! action parameters recognized by bele3
      enum ActionType
      {
        none
      };


      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool isHigherOrderElement(const CORE::FE::CellType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case CORE::FE::CellType::hex8:
          case CORE::FE::CellType::hex20:
          case CORE::FE::CellType::hex27:
          case CORE::FE::CellType::tet10:
          case CORE::FE::CellType::wedge15:
          case CORE::FE::CellType::nurbs8:
          case CORE::FE::CellType::nurbs27:
            hoel = true;
            break;
          case CORE::FE::CellType::tet4:
          case CORE::FE::CellType::wedge6:
          case CORE::FE::CellType::pyramid5:  //!!!TODO:  wedge und pyramid have 2nd
                                              //! derivatives!!!!!!!!!!!!!!!!!!!!!!!!
            hoel = false;
            break;
          default:
            FOUR_C_THROW("distype unknown!");
        }
        return hoel;
      };


      //! don't want = operator
      Vele3& operator=(const Vele3& old);

      //! set number of gauss points to element shape default
      CORE::FE::GaussRule3D getOptimalGaussrule(const CORE::FE::CellType& distype) const;

    };  // class Bele3



    class Vele3SurfaceType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Vele3SurfaceType"; }

      static Vele3SurfaceType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented for element type vele3!");
        return nullspace;
      }

     private:
      static Vele3SurfaceType instance_;
    };

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================
    class Vele3Surface : public DRT::FaceElement
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      explicit Vele3Surface(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Vele3* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit Vele3Surface(const Vele3Surface& old);

      DRT::Element* Clone() const override;
      CORE::FE::CellType Shape() const override;
      int NumLine() const override
      {
        if (NumNode() == 9 || NumNode() == 8 || NumNode() == 4)
          return 4;
        else if (NumNode() == 3 || NumNode() == 6)
          return 3;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }
      int NumSurface() const override { return 1; }
      int NumVolume() const override { return -1; }
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;
      int UniqueParObjectId() const override
      {
        return Vele3SurfaceType::Instance().UniqueParObjectId();
      }
      void Pack(CORE::COMM::PackBuffer& data) const override;
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods

      int NumDofPerNode(const DRT::Node&) const override { return 3; }
      int NumDofPerElement() const override { return 0; }
      void Print(std::ostream& os) const override;
      DRT::ElementType& ElementType() const override { return Vele3SurfaceType::Instance(); }

      //@}

      //@}

      //! @name Evaluation
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;


     private:
      //! action parameters recognized by bele3
      enum ActionType
      {
        none
      };

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool isHigherOrderElement(const CORE::FE::CellType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case CORE::FE::CellType::quad4:
          case CORE::FE::CellType::quad8:
          case CORE::FE::CellType::quad9:
          case CORE::FE::CellType::tri6:
            hoel = true;
            break;
          case CORE::FE::CellType::tri3:
            hoel = false;
            break;
          default:
            FOUR_C_THROW("distype unknown!");
        }
        return hoel;
      };


      //! don't want = operator
      Vele3Surface& operator=(const Vele3Surface& old);

      //! set number of gauss points to element shape default
      CORE::FE::GaussRule2D getOptimalGaussrule(const CORE::FE::CellType& distype) const;

    };  // class Vele3Surface


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    class Vele3LineType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Vele3LineType"; }

      static Vele3LineType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented for element type vele3!");
        return nullspace;
      }

     private:
      static Vele3LineType instance_;
    };


    /*!
    \brief An element representing a line of a vele3 element

    */
    class Vele3Line : public DRT::FaceElement
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
      \param parent: The parent fluid element of this line
      \param lline: the local line number of this line w.r.t. the parent element
      */
      Vele3Line(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::Element* parent, const int lline);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Vele3Line(const Vele3Line& old);

      DRT::Element* Clone() const override;
      CORE::FE::CellType Shape() const override;
      int UniqueParObjectId() const override
      {
        return Vele3LineType::Instance().UniqueParObjectId();
      }
      void Pack(CORE::COMM::PackBuffer& data) const override;
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      For this 3D boundary element, we have 3 displacements, if needed
      */
      int NumDofPerNode(const DRT::Node&) const override { return 3; }

      int NumDofPerElement() const override { return 0; }

      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return Vele3LineType::Instance(); }

      //@}


      //! @name Evaluation
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


      //! @name Evaluate methods
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     private:
      //! action parameters recognized by Vele3Line
      enum ActionType
      {
        none
      };

      //! don't want = operator
      Vele3Line& operator=(const Vele3Line& old);

      //! Get Rule for Gaussintegration according to DRT::UTIL
      CORE::FE::GaussRule1D getOptimalGaussrule(const CORE::FE::CellType& distype);

    };  // class Vele3Line


  }  // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
