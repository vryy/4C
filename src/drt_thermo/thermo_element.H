/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the thermo element

This file contains the element specific service routines like
Pack, Unpack, NumDofPerNode etc.

\level 1
\maintainer Sebastian Proell
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

#ifndef THERMO_ELEMENT_H
#define THERMO_ELEMENT_H

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

//#include "../drt_tsi/kinematictype.H"
#include "../drt_inpar/inpar_structure.H"


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
      std::string Name() const { return "ThermoType"; }

      static ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static ThermoType instance_;

    };  // class: ThermoType

    //!
    //! \brief A C++ wrapper for the thermo element
    //!
    //! \author danowski (danowski@lnm.mw.tum.de)
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
      DRT::Element* Clone() const;

      //! \brief Get shape type of element
      virtual DiscretizationType Shape() const;

      //! \brief set discretization type of element
      virtual void SetDisType(DRT::Element::DiscretizationType shape)
      {
        distype_ = shape;
        return;
      };

      //! \brief Return number of lines of this element
      virtual int NumLine() const { return DRT::UTILS::getNumberOfElementLines(distype_); }

      //! \brief Return number of surfaces of this element
      virtual int NumSurface() const
      {
        switch (distype_)
        {
          case DRT::Element::hex8:
          case DRT::Element::hex20:
          case DRT::Element::hex27:
          case DRT::Element::nurbs27:
            return 6;
            break;
          case DRT::Element::tet4:
          case DRT::Element::tet10:
            return 4;
            break;
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          case DRT::Element::pyramid5:
            return 5;
            break;
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::quad9:
          case DRT::Element::nurbs4:
          case DRT::Element::nurbs9:
          case DRT::Element::tri3:
          case DRT::Element::tri6:
            return 1;
            break;
          case DRT::Element::line2:
          case DRT::Element::line3:
            return 0;
            break;
          default:
            dserror("discretization type not yet implemented");
            break;
        }
        return 0;
      }

      //! \brief Return number of volumes of this element
      virtual int NumVolume() const
      {
        switch (distype_)
        {
          case DRT::Element::hex8:
          case DRT::Element::hex20:
          case DRT::Element::hex27:
          case DRT::Element::tet4:
          case DRT::Element::tet10:
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          case DRT::Element::pyramid5:
            return 1;
            break;
          case DRT::Element::quad4:
          case DRT::Element::quad8:
          case DRT::Element::quad9:
          case DRT::Element::nurbs4:
          case DRT::Element::nurbs9:
          case DRT::Element::tri3:
          case DRT::Element::tri6:
          case DRT::Element::line2:
          case DRT::Element::line3:
            return 0;
            break;
          default:
            dserror("discretization type not yet implemented");
            break;
        }
        return 0;
      }

      //! \brief Get vector of Teuchos::RCPs to the lines of this element
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      //! \brief Get vector of Teuchos::RCPs to the surfaces of this element
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      //! \brief Get vector of Teuchos::RCPs to the volume of this element
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      //! \brief Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the
      //! top of this file.
      virtual int UniqueParObjectId() const { return ThermoType::Instance().UniqueParObjectId(); }

      //! \brief Pack this class so it can be communicated
      //! \ref Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      //! \brief Unpack data from a char vector into this class
      //!
      //! \ref Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      //! \brief Destructor
      virtual ~Thermo();

      //@}

      //! @name Acess methods

      //! \brief Get number of degrees of freedom of a certain node
      //!        (implements pure virtual DRT::Element)
      //!
      //! The element decides how many degrees of freedom its nodes must have.
      //! As this may vary along a simulation, the element can redecide the
      //! number of degrees of freedom per node along the way for each of it's nodes
      //! separately.
      virtual int NumDofPerNode(const DRT::Node& node) const { return numdofpernode_; }

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
      virtual int NumDofPerElement() const { return 0; }

      //! \brief Print this element
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return ThermoType::Instance(); }

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
      virtual void VisNames(std::map<std::string, int>& names);

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
      virtual bool VisData(const std::string& name, std::vector<double>& data);

      //! \brief Add values to the data_ container
      virtual void AddToData(const std::string& name, std::vector<double> data)
      {
        data_.Add(name, data);
      };

      //@}

      //! @name Input and Creation

      //! \brief Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

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
      //! \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elevec1 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec2 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec3 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //! \brief Evaluate a Neumann boundary condition
      //!
      //! this method evaluates a surfaces Neumann condition on the shell element
      //!
      //! \param params (in/out)    : ParameterList for communication between control
      //!                             routine and elements
      //! \param discretization (in): A reference to the underlying discretization
      //! \param condition (in)     : The condition to be evaluated
      //! \param lm (in)            : location vector of this element
      //! \param elevec1 (out)      : vector to be filled by element. If NULL on input,
      //!
      //! \return 0 if successful, negative otherwise
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1);

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
      //! (storage neccessary because we don´t know the material in the post filters anymore)
      static const int numdofpernode_ = 1;
      //! the element discretization type
      DRT::Element::DiscretizationType distype_;
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
      std::string Name() const { return "ThermoBoundaryType"; }

      static ThermoBoundaryType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
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
      DRT::Element* Clone() const;

      //! \brief Get shape type of element
      virtual DiscretizationType Shape() const;

      //! \brief Return number of lines of boundary element
      virtual int NumLine() const
      {
        // get spatial dimension of boundary
        const int nsd = DRT::UTILS::getDimension(ParentElement()->Shape()) - 1;

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
      virtual int NumSurface() const
      {
        // get spatial dimension of parent element
        const int nsd = DRT::UTILS::getDimension(ParentElement()->Shape());

        if (nsd == 3)
          return 1;
        else
          return 0;
      }

      //! \brief Get vector of Teuchos::RCPs to the lines of this element
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      //! \brief Get vector of Teuchos::RCPs to the surfaces of this element
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      //! \brief Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the
      //! top of the drt_parobject.H file.
      virtual int UniqueParObjectId() const
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
      virtual void Unpack(const std::vector<char>& data);

      //! brief Destructor
      virtual ~ThermoBoundary();

      //@}

      //! @name Acess methods

      //! \brief Get number of degrees of freedom of a certain node
      //!       (implements pure virtual DRT::Element)
      //!
      //! The element decides how many degrees of freedom its nodes must have.
      //! As this may vary along a simulation, the element can redecide the
      //! number of degrees of freedom per node along the way for each of it's nodes
      //! separately.
      virtual int NumDofPerNode(const DRT::Node& node) const
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
      virtual int NumDofPerElement() const { return 0; }

      //! \brief Print this element
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return ThermoBoundaryType::Instance(); }

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
      //! \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element to fill
      //!                         this matrix.
      //! \param elevec1 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec2 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \param elevec3 (out)  : vector to be filled by element. If NULL on input,
      //!                         the controlling method does not expect the element
      //!                         to fill this vector
      //! \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

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
      //! \param elevec1 (out)      : vector to be filled by element. If NULL on input,
      //!
      //! \return 0 if successful, negative otherwise
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1);

      //@}

     private:
      // don't want = operator
      ThermoBoundary& operator=(const ThermoBoundary& old);

    };  // class ThermoBoundary

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // #ifndef THERMO_ELEMENT_H
#endif  // #ifdef D_THERMO
