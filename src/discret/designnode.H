/*!----------------------------------------------------------------------
\file designnode.H

\class CCADISCRETIZATION::DesignNode

\brief A node that is part of a CAD design description

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DESIGN_NODE_H
#define DESIGN_NODE_H



#include "node.H"

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class DesignNode

\brief A node that is part of a CAD design description

\author gee (gee@lnm.mw.tum.de)
*/
class DesignNode : public CCADISCRETIZATION::Node
{
public:

  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id    : A globally unique node id
  \param coords: vector of nodal coordinates, length 3
  */
  DesignNode(int id, const double* coords, const int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Design_Node
  
  */
  DesignNode(const CCADISCRETIZATION::DesignNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it
  
  */
  virtual CCADISCRETIZATION::DesignNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~DesignNode();

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this node
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual bool Unpack(const char* data);

  //@}

  // @{ \name Acess methods

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const ;

  //@}

  // @{ \name Construction


  //@}

protected:


}; // class DesignNode
} // namespace CCADISCRETIZATION











#endif  // #ifndef DESIGN_NODE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
