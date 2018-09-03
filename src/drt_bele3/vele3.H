/*----------------------------------------------------------------------*/
/*!
\file vele3.H

\brief volume element

\maintainer Jonas Eichinger

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef VELE3_H
#define VELE3_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
  // forward declarations
  class Discretization;


  namespace ELEMENTS
  {
    class Vele3Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Vele3Type"; }

      static Vele3Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

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
     * \author u.may (mayer@lnm.mw.tum.de)
     */


    /*!
    \brief A register for bele3 element

    \author u.may (mayer@lnm.mw.tum.de)
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
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const
      {
        if (NumNode() == 8 || NumNode() == 20 || NumNode() == 27)
          return 12;
        else if (NumNode() == 4 || NumNode() == 10)
          return 6;
        else
        {
          dserror("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const
      {
        if (NumNode() == 8 || NumNode() == 20 || NumNode() == 27)
          return 6;
        else if (NumNode() == 4 || NumNode() == 10)
          return 4;
        else if (NumNode() == 6 || NumNode() == 15 || NumNode() == 5)
          return 5;
        else
        {
          dserror("Could not determine number of surfaces");
          return -1;
        }
      }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      virtual int NumVolume() const { return 1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volume of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return Vele3Type::Instance().UniqueParObjectId(); }

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
      \brief Destructor

      */
      virtual ~Vele3();

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
      virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

      /*!
    \brief Get number of degrees of freedom per element
           (implements pure virtual DRT::Element)

    The element decides how many element degrees of freedom it has.
    It can redecide along the way of a simulation.

    \note Element degrees of freedom mentioned here are dofs that are visible
          at the level of the total system of equations. Purely internal
          element dofs that are condensed internally should NOT be considered.
    */
      virtual int NumDofPerElement() const { return 0; }

      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Vele3Type::Instance(); }


      //@}

      //@}

      //! @name Evaluation
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      /// Read input for this element
      virtual bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

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
      bool isHigherOrderElement(const DRT::Element::DiscretizationType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case hex8:
          case hex20:
          case hex27:
          case tet10:
          case wedge15:
          case nurbs8:
          case nurbs27:
            hoel = true;
            break;
          case tet4:
          case wedge6:
          case pyramid5:  //!!!TODO:  wedge und pyramid have 2nd derivatives!!!!!!!!!!!!!!!!!!!!!!!!
            hoel = false;
            break;
          default:
            dserror("distype unknown!");
        }
        return hoel;
      };


      //! don't want = operator
      Vele3& operator=(const Vele3& old);

      //! set number of gauss points to element shape default
      DRT::UTILS::GaussRule3D getOptimalGaussrule(const DiscretizationType& distype) const;

    };  // class Bele3



    class Vele3SurfaceType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Vele3SurfaceType"; }

      static Vele3SurfaceType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
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

      DRT::Element* Clone() const;
      virtual DiscretizationType Shape() const;
      virtual int NumLine() const
      {
        if (NumNode() == 9 || NumNode() == 8 || NumNode() == 4)
          return 4;
        else if (NumNode() == 3 || NumNode() == 6)
          return 3;
        else
        {
          dserror("Could not determine number of lines");
          return -1;
        }
      }
      virtual int NumSurface() const { return 1; }
      virtual int NumVolume() const { return -1; }
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();
      virtual int UniqueParObjectId() const
      {
        return Vele3SurfaceType::Instance().UniqueParObjectId();
      }
      virtual void Pack(DRT::PackBuffer& data) const;
      virtual void Unpack(const std::vector<char>& data);

      /*!
      \brief Destructor

      */
      virtual ~Vele3Surface();

      //@}

      //! @name Access methods

      virtual int NumDofPerNode(const DRT::Node&) const { return 3; }
      virtual int NumDofPerElement() const { return 0; }
      void Print(std::ostream& os) const;
      virtual DRT::ElementType& ElementType() const { return Vele3SurfaceType::Instance(); }

      //@}

      //@}

      //! @name Evaluation
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);


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
      bool isHigherOrderElement(const DRT::Element::DiscretizationType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case quad4:
          case quad8:
          case quad9:
          case tri6:
            hoel = true;
            break;
          case tri3:
            hoel = false;
            break;
          default:
            dserror("distype unknown!");
        }
        return hoel;
      };


      //! don't want = operator
      Vele3Surface& operator=(const Vele3Surface& old);

      //! set number of gauss points to element shape default
      DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype) const;

    };  // class Vele3Surface


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    class Vele3LineType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Vele3LineType"; }

      static Vele3LineType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
      }

     private:
      static Vele3LineType instance_;
    };


    /*!
    \brief An element representing a line of a vele3 element

    \author mayer (mayer@lnm.mw.tum.de)
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

      DRT::Element* Clone() const;
      virtual DiscretizationType Shape() const;
      virtual int UniqueParObjectId() const
      {
        return Vele3LineType::Instance().UniqueParObjectId();
      }
      virtual void Pack(DRT::PackBuffer& data) const;
      virtual void Unpack(const std::vector<char>& data);

      /*!
      \brief Destructor

      */
      virtual ~Vele3Line();

      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      For this 3D boundary element, we have 3 displacements, if needed
      */
      virtual int NumDofPerNode(const DRT::Node&) const { return 3; }

      virtual int NumDofPerElement() const { return 0; }

      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Vele3LineType::Instance(); }

      //@}


      //! @name Evaluation
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


      //! @name Evaluate methods
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

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
      DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

    };  // class Vele3Line


  }  // namespace ELEMENTS
}  // namespace DRT


#endif
