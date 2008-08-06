/*!
 * \file bele3.H
 * 
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 */
#ifdef CCADISCRET
#ifndef BELE3_H
#define BELE3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
// forward declarations
class Discretization;


namespace ELEMENTS
{

/*!
 * A 3D boundary element with no physics attached
 * 
 * This element is meant to have no physics. It can be used to have a boundary discretization
 * of surface/boundary elements. They can be of any 2d shape (quad4,quad9,tri3,...)
 * 
 * The number of dof per node is set to 3, so we can define displacement vectors by
 * using FillComplete on the boundary discretization.
 * 
 * \author axelchen (gerstenberger@lnm.mw.tum.de)
 */
class Bele3 : public DRT::Element
{
public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  explicit Bele3(
      int id,     ///< A unique global id
      int owner   ///< proc num that owns this element
      );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  explicit Bele3(const Bele3& old);

  DRT::Element* Clone() const;
  virtual DiscretizationType Shape() const;
  virtual int NumLine() const
  {
    if (NumNode()==9 || NumNode()==8 || NumNode()==4) return 4;
    else if (NumNode()==3 || NumNode()==6) return 3;
    else
    {
      dserror("Could not determine number of lines");
      return -1;
    }
  }
  virtual int NumSurface() const {return 1;}
  virtual int NumVolume() const {return -1;}
  virtual vector<RCP<DRT::Element> > Lines();
  virtual vector<RCP<DRT::Element> > Surfaces();
  virtual int UniqueParObjectId() const { return ParObject_Bele3;}
  virtual void Pack(vector<char>& data) const;
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3();

  //@}

  //! @name Access methods

  virtual int NumDofPerNode(const DRT::Node& node) const {return 3;}
  virtual int NumDofPerElement() const { return 0; }
  void Print(ostream& os) const;
  virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const;

  //@}

  //! @name Input and Creation

  /*!
   * \brief Read input for this element
   * \note Bele3 will be created on the fly for various tasks and no read in is planned
  */
  bool ReadElement()
  {
      dserror("No reading for this element! Will be created on the fly, not from a .dat file.");
      return false;
  }

  //@}

  //! @name Evaluation

  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other
  //! does this element have non-zero displacements or not
  bool IsMoving() const { return is_moving_; }

  //@}


private:

  //! action parameters recognized by bele3
  enum ActionType
  {
    none
  };

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype
    ) const
    {
      bool hoel = true;
      switch (distype)
      {
          case quad4: case quad8: case quad9: case tri6:
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
  
  //! container containing all kind of this and that
  DRT::Container                                  data_;
  
  //! flag for fixed or moving boundary
  bool                                            is_moving_;

  //! don't want = operator
  Bele3& operator = (const Bele3& old);
  
  //! set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype);

}; // class Bele3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for bele3 element

\author axelchen (gerstenberger@lnm.mw.tum.de)
*/
class Bele3Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  explicit Bele3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Bele3Register

  */
  Bele3Register(const DRT::ELEMENTS::Bele3Register& old);

  /*!
  \brief Deep copy this instance of Bele3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Bele3Register* Clone() const;

  virtual int UniqueParObjectId() const { return ParObject_Bele3Register;}
  virtual void Pack(vector<char>& data) const;
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3Register();

  //@}

  //! @name Acess methods

  void Print(ostream& os) const;

  //@}

  //! @name Construction

  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Bele3Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a line of a bele3 element

\author gerstenberger (gerstenberger@lnm.mw.tum.de)
*/
class Bele3Line : public DRT::Element
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
  Bele3Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Bele3* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Bele3Line(const Bele3Line& old);

  DRT::Element* Clone() const;
  virtual DiscretizationType Shape() const;
  virtual int UniqueParObjectId() const { return ParObject_Bele3Line;}
  virtual void Pack(vector<char>& data) const;
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3Line();

  //@}

  //! @name Access methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  For this 3D boundary element, we have 3 displacements, if needed
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return 3;}

  virtual int NumDofPerElement() const { return 0; }

  void Print(ostream& os) const;

  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const { return null; }

  //@}

    //! @name Evaluation

  virtual int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //! @name Evaluate methods

  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //@}

private:

  //! action parameters recognized by Bele3Line
  enum ActionType
  {
     none,
     integrate_Shapefunction
  };

  //! don't want = operator
  Bele3Line& operator = (const Bele3Line& old);

  //! The parent element of this line
  DRT::ELEMENTS::Bele3* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;



  //! compute infintesimal line element dr for integration along the line
  double  f2_substitution(
		      const Epetra_SerialDenseMatrix  xyze,
		      const Epetra_SerialDenseMatrix  deriv,
		      const int iel);

  //! Get Rule for Gaussintegration according to DRT::UTIL
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  //! integrate shape functions over a line
  void IntegrateShapeFunction(
      ParameterList&             params,
      DRT::Discretization&       discretization,
      const vector<int>&         lm,
      Epetra_SerialDenseVector&  elevec1,
      const std::vector<double>& edispnp
      );


}; // class Bele3Line




} // namespace ELEMENTS
} // namespace DRT


#endif
#endif
