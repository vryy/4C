/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the thermo element

This file contains the element specific service routines like
Pack, Unpack, NumDofPerNode etc.

\level 1
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_ELEMENT_H
#define BACI_THERMO_ELEMENT_H

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_serialdensematrix.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class ThermoBoundary;

    class ThermoType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ThermoType"; }

      static ThermoType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static ThermoType instance_;

    };  // class: ThermoType

    //!
    //! \brief A C++ wrapper for the thermo element
    //!
    class Thermo : public DRT::Element
    {
     public:
      //! @name Friends
      friend class ThermoBoundary;

      //@}
      //! @name Constructors and destructors and related methods

      //! \brief Standard Constructor
      Thermo(int id,  ///< A unique global id
          int owner   ///< processor id who owns a certain instance of this class
      );

      //! \brief Copy Constructor
      //!
      //! Makes a deep copy of a Element
      Thermo(const Thermo& old);

      //! \brief Deep copy this instance of Thermo and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const override;

      //! \brief Get shape type of element
      CORE::FE::CellType Shape() const override;

      //! \brief set discretization type of element
      virtual void SetDisType(CORE::FE::CellType shape)
      {
        distype_ = shape;
        return;
      };

      //! \brief Return number of lines of this element
      int NumLine() const override { return CORE::DRT::UTILS::getNumberOfElementLines(distype_); }

      //! \brief Return number of surfaces of this element
      int NumSurface() const override
      {
        switch (distype_)
        {
          case CORE::FE::CellType::hex8:
          case CORE::FE::CellType::hex20:
          case CORE::FE::CellType::hex27:
          case CORE::FE::CellType::nurbs27:
            return 6;
            break;
          case CORE::FE::CellType::tet4:
          case CORE::FE::CellType::tet10:
            return 4;
            break;
          case CORE::FE::CellType::wedge6:
          case CORE::FE::CellType::wedge15:
          case CORE::FE::CellType::pyramid5:
            return 5;
            break;
          case CORE::FE::CellType::quad4:
          case CORE::FE::CellType::quad8:
          case CORE::FE::CellType::quad9:
          case CORE::FE::CellType::nurbs4:
          case CORE::FE::CellType::nurbs9:
          case CORE::FE::CellType::tri3:
          case CORE::FE::CellType::tri6:
            return 1;
            break;
          case CORE::FE::CellType::line2:
          case CORE::FE::CellType::line3:
            return 0;
            break;
          default:
            dserror("discretization type not yet implemented");
            break;
        }
        return 0;
      }

      //! \brief Return number of volumes of this element
      int NumVolume() const override
      {
        switch (distype_)
        {
          case CORE::FE::CellType::hex8:
          case CORE::FE::CellType::hex20:
          case CORE::FE::CellType::hex27:
          case CORE::FE::CellType::tet4:
          case CORE::FE::CellType::tet10:
          case CORE::FE::CellType::wedge6:
          case CORE::FE::CellType::wedge15:
          case CORE::FE::CellType::pyramid5:
            return 1;
            break;
          case CORE::FE::CellType::quad4:
          case CORE::FE::CellType::quad8:
          case CORE::FE::CellType::quad9:
          case CORE::FE::CellType::nurbs4:
          case CORE::FE::CellType::nurbs9:
          case CORE::FE::CellType::tri3:
          case CORE::FE::CellType::tri6:
          case CORE::FE::CellType::line2:
          case CORE::FE::CellType::line3:
            return 0;
            break;
          default:
            dserror("discretization type not yet implemented");
            break;
        }
        return 0;
      }

      //! \brief Get vector of Teuchos::RCPs to the lines of this element
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      //! \brief Get vector of Teuchos::RCPs to the surfaces of this element
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //! \brief Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the
      //! top of this file.
      int UniqueParObjectId() const override { return ThermoType::Instance().UniqueParObjectId(); }

      //! \brief Pack this class so it can be communicated
      //! \ref Pack and \ref Unpack are used to communicate this element
      void Pack(CORE::COMM::PackBuffer& data) const override;

      //! \brief Unpack data from a char vector into this class
      //!
      //! \ref Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods

      //! \brief Get number of degrees of freedom of a certain node
      //!        (implements pure virtual DRT::Element)
      //!
      //! The element decides how many degrees of freedom its nodes must have.
      //! As this may vary along a simulation, the element can redecide the
      //! number of degrees of freedom per node along the way for each of it's nodes
      //! separately.
      int NumDofPerNode(const DRT::Node& node) const override { return numdofpernode_; }

      //!
      //! \brief Get number of degrees of freedom per element
      //!        (implements pure virtual DRT::Element)
      //!
      //! The element decides how many element degrees of freedom it has.
      //! It can redecide along the way of a simulation.
      //!
      //! \note Element degrees of freedom mentioned here are dofs that are visible
      //!       at the level of the total system of equations. Purely internal
      //!       element dofs that are condensed internally should NOT be considered.
      int NumDofPerElement() const override { return 0; }

      //! \brief Print this element
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ThermoType::Instance(); }

      //! \brief Query names of element data to be visualized using BINIO
      //!
      //! The element fills the provided map with key names of
      //! visualization data the element wants to visualize AT THE CENTER
      //! of the element geometry. The value is supposed to be dimension of the
      //! data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym.
      //! tensor) or 9 (nonsym. tensor)
      //!
      //! Example:
      //! \code
      //!  // Name of data is 'Owner', dimension is 1 (scalar value)
      //!  names.insert(std::pair<std::string,int>("Owner",1));
      //!  // Name of data is 'HeatfluxXYZ', dimension is 3 (vector value)
      //!  names.insert(std::pair<std::string,int>("HeatfluxXYZ",3));
      //! \endcode
      //!
      //! \param names (out): On return, the derived class has filled names with
      //!                     key names of data it wants to visualize and with int
      //!                     dimensions of that data.
      void VisNames(std::map<std::string, int>& names) override;

      //! \brief Query data to be visualized using BINIO of a given name
      //!
      //! The method is supposed to call this base method to visualize the owner of
      //! the element.
      //! If the derived method recognizes a supported data name, it shall fill it
      //! with corresponding data.
      //! If it does NOT recognizes the name, it shall do nothing.
      //!
      //! \warning The method must not change size of data
      //!
      //! \param name (in):   Name of data that is currently processed for visualization
      //! \param data (out):  data to be filled by element if element recognizes the name
      bool VisData(const std::string& name, std::vector<double>& data) override;

      //! \brief Add values to the data_ container
      virtual void AddToData(const std::string& name, std::vector<double> data)
      {
        data_.Add(name, data);
      };

      //@}

      //! @name Input and Creation

      //! \brief Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      //! \brief Evaluate an element, i.e. call the implementation to evaluate element
      //! tangent, capacity, internal forces or evaluate errors, statistics or updates
      //! etc. directly.
      //!
      //! Following implementations of the element are allowed:
      //!       //!  o Evaluation of thermal system matrix and residual for the One-Step-Theta
      //!
      //!  o Evaluation of thermal system matrix and residual for the stationary thermal solver
      //!       //!
      //! \param params (in/out): ParameterList for communication between control routine
      //!                         and elements
      //! \param discretization (in): A reference to the underlying discretization
      //! \param la (in)        : location array of this element
      //! \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //! \brief Evaluate a Neumann boundary condition
      //!
      //! this method evaluates a surfaces Neumann condition on the shell element
      //!
      //! \param params (in/out)    : ParameterList for communication between control
      //!                             routine and elements
      //! \param discretization (in): A reference to the underlying discretization
      //! \param condition (in)     : The condition to be evaluated
      //! \param lm (in)            : location vector of this element
      //! \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      //!
      //! \return 0 if successful, negative otherwise
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1) override;

      //@}

      //! kinematic type passed from structural element
      virtual void SetKinematicType(INPAR::STR::KinemType kintype)
      {
        kintype_ = kintype;
        return;
      };
      //! kinematic type
      INPAR::STR::KinemType kintype_;

      INPAR::STR::KinemType KinType() const { return kintype_; }

     private:
      //! number of dofs per node (for systems of thermo equations)
      //! (storage neccessary because we don't know the material in the post filters anymore)
      static const int numdofpernode_ = 1;
      //! the element discretization type
      CORE::FE::CellType distype_;
      //! container containing all kind of element this-and-that
      DRT::Container data_;

      //! don't want = operator
      Thermo& operator=(const Thermo& old);

    };  // class Thermo


    ////=======================================================================
    ////=======================================================================
    ////=======================================================================
    ////=======================================================================
    class ThermoBoundaryType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ThermoBoundaryType"; }

      static ThermoBoundaryType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented");
        return nullspace;
      }

     private:
      static ThermoBoundaryType instance_;
    };

    //! \brief An element representing a boundary element of a thermo element
    //!
    //! \note This is a pure boundary condition element. It's only
    //!       purpose is to evaluate certain boundary conditions that might be
    //!       adjacent to a parent Thermo element.
    class ThermoBoundary : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      //! \brief Standard Constructor
      //!
      //! \param id : A unique global id
      //! \param owner: Processor owning this surface
      //! \param nnode: Number of nodes attached to this element
      //! \param nodeids: global ids of nodes attached to this element
      //! \param nodes: the discretization map of nodes to build ptrs to nodes from
      //! \param parent: The parent fluid element of this surface
      //! \param lsurface: the local surface number of this surface w.r.t. the parent element
      ThermoBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Thermo* parent, const int lsurface);

      //! \brief Copy Constructor
      //!
      //! Makes a deep copy of a Element
      ThermoBoundary(const ThermoBoundary& old);

      //! \brief Deep copy this instance of an element and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-constructor is needed
      DRT::Element* Clone() const override;

      //! \brief Get shape type of element
      CORE::FE::CellType Shape() const override;

      //! \brief Return number of lines of boundary element
      int NumLine() const override
      {
        // get spatial dimension of boundary
        const int nsd = CORE::DRT::UTILS::getDimension(ParentElement()->Shape()) - 1;

        if ((NumNode() == 4) or (NumNode() == 8) or (NumNode() == 9))
          return 4;
        else if (NumNode() == 6)
          return 3;
        else if ((NumNode() == 3) and (nsd == 2))
          return 3;
        else if ((NumNode() == 3) and (nsd == 1))
          return 1;
        else if (NumNode() == 2)
          return 1;
        else
        {
          dserror("Could not determine number of lines");
          return -1;
        }
      }

      //! \brief Return number of surfaces of boundary element
      int NumSurface() const override
      {
        // get spatial dimension of parent element
        const int nsd = CORE::DRT::UTILS::getDimension(ParentElement()->Shape());

        if (nsd == 3)
          return 1;
        else
          return 0;
      }

      //! \brief Get vector of Teuchos::RCPs to the lines of this element
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      //! \brief Get vector of Teuchos::RCPs to the surfaces of this element
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //! \brief Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the
      //! top of the parobject.H file.
      int UniqueParObjectId() const override
      {
        return ThermoBoundaryType::Instance().UniqueParObjectId();
      }

      //! \brief Pack this class so it can be communicated
      //!
      //! \ref Pack and \ref Unpack are used to communicate this element
      virtual void Pack(std::vector<char>& data) const;

      //! \brief Unpack data from a char vector into this class
      //!
      //! \ref Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods

      //! \brief Get number of degrees of freedom of a certain node
      //!       (implements pure virtual DRT::Element)
      //!
      //! The element decides how many degrees of freedom its nodes must have.
      //! As this may vary along a simulation, the element can redecide the
      //! number of degrees of freedom per node along the way for each of it's nodes
      //! separately.
      int NumDofPerNode(const DRT::Node& node) const override
      {
        return ParentElement()->NumDofPerNode(node);
      }

      /*
      //! Return a pointer to the parent element of this boundary element
      virtual DRT::ELEMENTS::Thermo* ParentElement()
      {
        return parent_;
      }
      */

      //! \brief Get number of degrees of freedom per element
      //!       (implements pure virtual DRT::Element)
      //!
      //! The element decides how many element degrees of freedom it has.
      //! It can redecide along the way of a simulation.
      //!
      //! \note Element degrees of freedom mentioned here are dofs that are visible
      //!      at the level of the total system of equations. Purely internal
      //!      element dofs that are condensed internally should NOT be considered.
      int NumDofPerElement() const override { return 0; }

      //! \brief Print this element
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ThermoBoundaryType::Instance(); }

      //@}

      //! @name Evaluation

      //! \brief Evaluate an element
      //!
      //! Evaluate Thermo element tangent, capacity, internal forces etc
      //!
      //! \param params (in/out): ParameterList for communication between control routine
      //!                         and elements
      //! \param discretization (in): A reference to the underlying discretization
      //! \param la (in):         location array of this element, vector of
      //!                         degrees of freedom adressed by this element
      //! \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //@}

      //! @name Evaluate methods

      //! \brief Evaluate a Neumann boundary condition
      //!
      //! this method evaluates a surface Neumann condition on the thermo element
      //!
      //! \param params (in/out)    : ParameterList for communication between control routine
      //!                             and elements
      //! \param discretization (in): A reference to the underlying discretization
      //! \param condition (in)     : The condition to be evaluated
      //! \param lm (in)            : location vector of this element
      //! \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      //!
      //! \return 0 if successful, negative otherwise
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1) override;

      //@}

     private:
      // don't want = operator
      ThermoBoundary& operator=(const ThermoBoundary& old);

    };  // class ThermoBoundary

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // THERMO_ELEMENT_H
