/*!----------------------------------------------------------------------
\file drt_celement.H

\class CONTACT::CElement

\brief A contact element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_CELEMENT_H
#define DRT_CELEMENT_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "drt_cnode.H"
using namespace std;
using namespace Teuchos;

// forward declarations

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class CElement

\brief A contact element


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CElement : public DRT::Element
{
public:

  //! @name Enums and Friends
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  \param shape (in): shape of this element
  \param numnode (in): Number of nodes to this element
  \param nodeids (in): ids of nodes adjacent to this element
  */
  explicit CElement(int id, ElementType etype, int owner, 
                    const DRT::Element::DiscretizationType& shape, 
                    const int numnode,
                    const int* nodeids, 
                    const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this class

  */
  CElement(const CONTACT::CElement& old);

  /*!
  \brief Destructor

  */
  virtual ~CElement() { return; }

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CElement* Clone() const;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return ParObject_CElement; }

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
  \brief Allocate and return a ElementRegister class

   Assuming the contact elements do not need initialization, we return a 
   dummy base class here.

  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const
  {
    RCP<DRT::ElementRegister> reg = rcp(new DRT::ElementRegister(Type()));
    return reg;
  }

  //@}

  //! @name Query methods


  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const { return shape_;}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces to this element
  */
  virtual int NumSurface() const {return 0;}


  /*!
  \brief Get vector of ptrs to the lines of this element

  This is a base class dummy routine that always returns NULL.
  The derived element class is expected to allocate and store
  a vector of elements that represent the edges of this element.
  These edges are then used to create and evaluate boundary conditions.

  \note A 1D type of element (e.g. a beam) may return
       a vector of length one pointing to itself. It then does not need
       to be able to spin of an explicit separate line element.
       In this case, the element itself must be able to evaluate the
       Neumann boundary conditions on a line itself.

  */
  virtual Element** Lines() { return NULL; }

  /*!
  \brief Get vector of ptrs to the surfaces of this element

  This is a base class dummy routine that always returns NULL.
  The derived element class is expected to allocate and store
  a vector of elements that represent the surfaces of this element.
  These surfaces are then used to create and evaluate boundary conditions.

  \note A 2D type of element (e.g. a shell, wall, fluid2 etc) may return
       a vector of length one pointing to itself. It then does not need
       to be able to spin of an explicit separate surface element.
       In this case, the element itself must be able to evaluate the
       Neumann boundary conditions on a surface itself.

  */
  virtual Element** Surfaces() { return NULL; }


  /*!
  \brief Get number of degrees of freedom of a certain node

  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    const CONTACT::CNode* cnode = dynamic_cast<const CONTACT::CNode*>(&node);
    if (!cnode) dserror("Node is not a CNode");
    return cnode->NumDof();
  }

  /*!
  \brief Get number of degrees of freedom per element

  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element

  */
  virtual void Print(ostream& os) const;

  //@}

  //! @name Construction


  //@}

  //! @name Conditions

  //@}

  //! @name Evaluation methods



  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elemat2 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elevec1 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec2 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec3 (out)      : vector to be filled by element depending on commands
                              given in params
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);



  //@}


private:

  DRT::Element::DiscretizationType shape_;   // shape of this element
  bool                             isslave_; // indicating slave or master side

}; // class CElement
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CElement& ele);









#endif  // #ifndef DRT_CELEMENT_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
