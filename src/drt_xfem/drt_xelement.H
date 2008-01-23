/*!----------------------------------------------------------------------
\file drt_xelement.H

\class DRT::ELEMENTS::XElement

\brief A pure virtual XFEM element class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef XELEMENT_H
#define XELEMENT_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "physics.H"
#include "dof_management.H"

using namespace Teuchos;

class XFEM::Enrichment;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Discretization;
class Node;
class Element;
class ElementRegister;

namespace ELEMENTS
{

/*!
\class XElement

\brief A virtual class all elements that are used in DRT have to implement

This is the pure virtual base class for all finite elements to be used with
the DRT::Discretization. Every element (and boundary condition) to be used
with the discretization management module has to implement this class.
It implements various basic element methods and
stores basic information such as element to node connectivity.

\author gerstenberger (gerstenberger@lnm.mw.tum.de)
*/
class XElement : public Element
{
public:

  typedef int EleId;    ///< unique elemend id

  typedef int ProcId;   ///< processor id

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  XElement(const int id, const ElementType etype, const int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  XElement(const XElement& old);

  /*!
  \brief Destructor

  */
  virtual ~XElement();

  /*!
  \brief Deep copy the derived class and return pointer to it

  This method is sort of a copy constructor for a class derived from Element.
  It allows to copy construct the derived class without knowing what it
  actually is using the base class Element.

  */
  virtual XElement* Clone() const = 0;


  //@}

  //! @name Query methods

  /*!
  \brief Get number of degrees of freedom of a certain node

  Implementation for XFEM elements: the dof management for XFEM elements is generalized
  such that it can handle the dof management without knowing about the particular
  physics used. All constraints that the standard element \ref Element asks for are fulfilled
  */
  virtual int NumDofPerNode(const DRT::Node& node) const = 0;

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
  virtual int EvaluateShape(const ParameterList&            params,
                            const DRT::Discretization&      discretization,
                            const vector<int>&              lm,
                            Epetra_SerialDenseMatrix& elemat1,
                            Epetra_SerialDenseMatrix& elemat2,
                            Epetra_SerialDenseVector& elevec1,
                            Epetra_SerialDenseVector& elevec2,
                            Epetra_SerialDenseVector& elevec3);

  //@}



private:

}; // class XElement
}; // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef XELEMENT_H
#endif  // #ifdef CCADISCRET
