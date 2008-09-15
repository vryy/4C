/*!----------------------------------------------------------------------
\file so_line.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_LINE_H
#define SO_LINE_H


#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_fem_general/drt_utils_integration.H"



namespace DRT
{

namespace ELEMENTS
{

/*!
\brief An element representing a line edge of any 3D structural element
       Note that this element should not be used in 2D cases!

\author gee (gee@lnm.mw.tum.de)
*/
class StructuralLine : public DRT::Element
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
  StructuralLine(int id, int owner, int nnode, const int* nodeids,
                 DRT::Node** nodes, DRT::Element* parent, const int lline);

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
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual inline int UniqueParObjectId() const
  { return ParObject_StructuralLine;}

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
  virtual ~StructuralLine() {}

  //@}

  //! @name Acess methods

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual inline int NumDofPerNode(const DRT::Node& node) const
  {return 3;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual inline int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Allocate and return a ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual inline RCP<DRT::ElementRegister> ElementRegister() const
  { return rcp(new DRT::ElementRegister(Type())); }

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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //@}

private:

  // don't want = operator
  StructuralLine& operator = (const StructuralLine& old);

  //! The parent element of this line
  DRT::Element* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int           lline_;
  //! gaussian integration to be used
  DRT::UTILS::GaussRule1D gaussrule_;

  //! line integration
  void LineIntegration(double&                          dL,
                       const LINALG::SerialDenseMatrix& x,
                       const LINALG::SerialDenseMatrix& deriv);

  /*!
  \brief Create matrix with material configuration

  \param x  (out)  : nodal coords in material frame
  */
  inline void MaterialConfiguration(LINALG::SerialDenseMatrix& x) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = Nodes()[i]->X()[0];
      x(i,1) = Nodes()[i]->X()[1];
      x(i,2) = Nodes()[i]->X()[2];
    }
    return;
  }


}; // class StructuralLine


} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SO_LINE_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
