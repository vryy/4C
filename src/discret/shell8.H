/*!----------------------------------------------------------------------
\file shell8.H

\class DRT::Elements::Shell8

\brief A C++ wrapper for the shell8 element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SHELL8_H
#define SHELL8_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "drt_element.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
namespace Elements
{

/*!
\class Shell8

\brief A C++ wrapper for the shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Shell8 : public DRT::Element
{
public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param id : A unique global id
  */
  Shell8(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  
  */
  Shell8(const Shell8& old);

  /*!
  \brief Deep copy this instance of Shell8 and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::Element* Clone() const;

  /*!
  \brief Pack this class so it can be communicated
  
  \ref Pack and \ref Unpack are used to communicate this element
  
  */
  virtual const char* Pack(int& size) const;

  /*!
  \brief Unpack data from a char vector into this class
  
  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual bool Unpack(const char* data);

  /*!
  \brief Destructor

  */
  virtual ~Shell8();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)
  
  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 6;}

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

  //@}

  // @{ \name Other

  /*!
  \brief Read input for this element
  */
  bool ReadElement();

  /*!
  \brief Evaluate an element
  
  Evaluate shell8 element stuffness, mass, internal forces etc
  
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

  //@}


private:

  // don't want = operator
  Shell8 operator = (const Shell8& old);

  enum ActionType
  {
    none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload
  };

  // element specific data
  enum ForceType
  {
    s8_none,
    s8_xyz,
    s8_rst,
    s8_rst_ortho
  };
  ForceType        forcetype_;  // type of force calculation
  double           thickness_;  // shell thickness
  int              ngp_[3];     // no. gaussian points in each spatial dimension
  int              ngptri_;     // no. gaussian points for triangle elements
  int              nhyb_;       // no. of internal dofs
  int              eas_[5];     // eas parameters
  int              ans_;        // type of ANS: 0=none 1=Q 2=T 3=QT  
  double           sdc_;        // scaled director conditioning
  int              material_;   // number of the material law

}; // class Shell8
} // namespace Elements
} // namespace DRT











#endif  // #ifndef SHELL8_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
