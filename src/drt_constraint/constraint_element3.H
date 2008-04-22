/*!----------------------------------------------------------------------
\file constraint_element3.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONSTRAINT_ELEMENT3_H
#define CONSTRAINT_ELEMENT3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"



namespace DRT
{
// forward declarations
class Discretization;
class ConstraintElement3Register;

namespace ELEMENTS
{
// forward declarations
//class ConstraintElement3Line;

/*!
\brief A 3D constraint element with no physics attached

\author tk (kloeppel@lnm.mw.tum.de)
*/
class ConstraintElement3 : public DRT::Element
{
public:

  //! @name Friends
  friend class ConstraintElement3Register;
  friend class ConstraintElement3Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  ConstraintElement3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  ConstraintElement3(const ConstraintElement3& old);

  /*!
  \brief Deep copy this instance of ConstraintElement3 and return pointer to the copy

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
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of ptrs to the volume of this element

  */
  virtual DRT::Element** Lines();

  /*!
  \brief Get vector of ptrs to the surfaces of this element

  */
  virtual DRT::Element** Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ConstraintElement3;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~ConstraintElement3();

  //@}

  //! @name Acess methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 3;}

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

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Allocate and return a DesignElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const;

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement()
  {
      dserror("No reading for this element! Will be created on the fly, not from a .dat file.");
      return false;
  }

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate fluid2 element stiffness, mass, internal forces etc

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \return 0 if successful, negative otherwise
  */
  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the fluid2 element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


  //@}

private:

  //! action parameters recognized by bele3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual
  };

  //! container containing all kind of fluid2 thisandthat
  DRT::Container                                  data_;

  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surface_;

//  //! vector of lines adjacent to this element (fluid2: 3 or 4 lines)
//  vector<RefCountPtr<DRT::ELEMENTS::ConstraintElement3Line> > lines_;
//  //! vector of lines adjacent to this element (fluid2: 3 or 4 lines)
//  vector<DRT::ELEMENTS::ConstraintElement3Line*>              lineptrs_;

  // internal calculation methods

  // don't want = operator
  ConstraintElement3& operator = (const ConstraintElement3& old);

}; // class ConstraintElement3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for bele3 element

\author axelchen (gerstenberger@lnm.mw.tum.de)
*/
class ConstraintElement3Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  ConstraintElement3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an ConstraintElement3Register

  */
  ConstraintElement3Register(const DRT::ELEMENTS::ConstraintElement3Register& old);

  /*!
  \brief Deep copy this instance of ConstraintElement3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::ConstraintElement3Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ConstraintElement3Register;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~ConstraintElement3Register();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}

  //! @name Construction

  /*!
  \brief Call Element init routines

  The method performs all steps needed to init the fluid2 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class ConstraintElement3Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

//
///*!
//\brief An element representing a line of a bele3 element
//
//\note This is a pure Neumann boundary condition element. It's only
//      purpose is to evaluate line Neumann boundary conditions that might be
//      adjacent to a parent fluid2 element. It therefore does not implement
//      the DRT::Element::Evaluate method and does not have its own ElementRegister class.
//
//\author gamnitzer (gamnitzer@lnm.mw.tum.de)
//*/
//class ConstraintElement3Line : public DRT::Element
//{
//public:
//
//  //! @name Constructors and destructors and related methods
//
//  /*!
//  \brief Standard Constructor
//
//  \param id : A unique global id
//  \param owner: Processor owning this line
//  \param nnode: Number of nodes attached to this element
//  \param nodeids: global ids of nodes attached to this element
//  \param nodes: the discretizations map of nodes to build ptrs to nodes from
//  \param parent: The parent fluid element of this line
//  \param lline: the local line number of this line w.r.t. the parent element
//  */
//  ConstraintElement3Line(int id, int owner, int nnode, const int* nodeids,
//             DRT::Node** nodes, DRT::ELEMENTS::ConstraintElement3* parent, const int lline);
//
//  /*!
//  \brief Copy Constructor
//
//  Makes a deep copy of a Element
//
//  */
//  ConstraintElement3Line(const ConstraintElement3Line& old);
//
//  /*!
//  \brief Deep copy this instance of an element and return pointer to the copy
//
//  The Clone() method is used from the virtual base class Element in cases
//  where the type of the derived class is unknown and a copy-ctor is needed
//
//  */
//  DRT::Element* Clone() const;
//
//  /*!
//  \brief Get shape type of element
//  */
//  virtual DiscretizationType Shape() const;
//
//  /*!
//  \brief Return unique ParObject id
//
//  every class implementing ParObject needs a unique id defined at the
//  top of the drt_parobject.H file.
//  */
//  virtual int UniqueParObjectId() const { return ParObject_ConstraintElement3Line;}
//
//  /*!
//  \brief Pack this class so it can be communicated
//
//  \ref Pack and \ref Unpack are used to communicate this element
//
//  */
//  virtual void Pack(vector<char>& data) const;
//
//  /*!
//  \brief Unpack data from a char vector into this class
//
//  \ref Pack and \ref Unpack are used to communicate this element
//
//  */
//  virtual void Unpack(const vector<char>& data);
//
//  /*!
//  \brief Destructor
//
//  */
//  virtual ~ConstraintElement3Line();
//
//  //@}
//
//  //! @name Access methods
//
//
//  /*!
//  \brief Get number of degrees of freedom of a certain node
//         (implements pure virtual DRT::Element)
//
//  The element decides how many degrees of freedom its nodes must have.
//  As this may vary along a simulation, the element can redecide the
//  number of degrees of freedom per node along the way for each of it's nodes
//  separately.
//  */
//  virtual int NumDofPerNode(const DRT::Node& node) const
//  {return 3;}
//
//  /*!
//  \brief Get number of degrees of freedom per element
//         (implements pure virtual DRT::Element)
//
//  The element decides how many element degrees of freedom it has.
//  It can redecide along the way of a simulation.
//
//  \note Element degrees of freedom mentioned here are dofs that are visible
//        at the level of the total system of equations. Purely internal
//        element dofs that are condensed internally should NOT be considered.
//  */
//  virtual int NumDofPerElement() const { return 0; }
//
//  /*!
//  \brief Print this element
//  */
//  void Print(ostream& os) const;
//
//  /*!
//  \brief Allocate and return a ElementRegister class
//
//  Every implementation of DRT::Element must also implement
//  an DRT::ElementRegister. This method allocates and returns
//  an instance of a matching DRT::ElementRegister.
//  */
//  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const { return null; }
//
//  //@}
//
//    //! @name Evaluation
//
//  /*!
//  \brief Evaluate an element
//
//  Evaluate fluid2 line element stiffness, mass, internal forces etc
//
//  \param params (in/out): ParameterList for communication between control routine
//                          and elements
//  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
//                          the controling method does not epxect the element to fill
//                          this matrix.
//  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
//                          the controling method does not epxect the element to fill
//                          this matrix.
//  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
//                          the controlling method does not epxect the element
//                          to fill this vector
//  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
//                          the controlling method does not epxect the element
//                          to fill this vector
//  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
//                          the controlling method does not epxect the element
//                          to fill this vector
//  \return 0 if successful, negative otherwise
//  */
//  virtual int Evaluate(ParameterList& params,
//               DRT::Discretization&      discretization,
//               vector<int>&              lm,
//               Epetra_SerialDenseMatrix& elemat1,
//               Epetra_SerialDenseMatrix& elemat2,
//               Epetra_SerialDenseVector& elevec1,
//               Epetra_SerialDenseVector& elevec2,
//               Epetra_SerialDenseVector& elevec3);
//
//  //! @name Evaluate methods
//
//  /*!
//  \brief Evaluate a Neumann boundary condition
//
//  this method evaluates a line Neumann condition on the fluid2 element
//
//  \param params (in/out)    : ParameterList for communication between control routine
//                              and elements
//  \param discretization (in): A reference to the underlying discretization
//  \param condition (in)     : The condition to be evaluated
//  \param lm (in)            : location vector of this element
//  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
//
//  \return 0 if successful, negative otherwise
//  */
//  virtual int EvaluateNeumann(ParameterList& params,
//                              DRT::Discretization&      discretization,
//                              DRT::Condition&           condition,
//                              vector<int>&              lm,
//                              Epetra_SerialDenseVector& elevec1);
//
//  //@}
//
//private:
//
//  //! action parameters recognized by ConstraintElement3Line
//  enum ActionType
//  {
//     none,
//     integrate_Shapefunction
//  };
//
//  // don't want = operator
//  ConstraintElement3Line& operator = (const ConstraintElement3Line& old);
//
//  //! The parent element of this line
//  DRT::ELEMENTS::ConstraintElement3* parent_;
//  //! The local line number of this line w.r.t to the parent_ element
//  int                    lline_;
//
//
//
//  // compute infintesimal line element dr for integration along the line
//  double  f2_substitution(
//		      const Epetra_SerialDenseMatrix  xyze,
//		      const Epetra_SerialDenseMatrix  deriv,
//		      const int iel);
//
//  // Get Rule for Gaussintegration according to DRT::UTIL
//  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);
//
//  // integrate shape functions over a line
//  void IntegrateShapeFunction(ParameterList& params,
//                  DRT::Discretization&       discretization,
//                  vector<int>&               lm,
//                  Epetra_SerialDenseVector&  elevec1,
//                  const std::vector<double>& edispnp);
//
//
//}; // class ConstraintElement3Line




} // namespace ELEMENTS
} // namespace DRT


#endif
#endif
