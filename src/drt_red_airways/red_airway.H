/*!----------------------------------------------------------------------
\file red_airway.H

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_RED_AIRWAYS
#ifdef CCADISCRET
#ifndef RED_AIRWAY_H_
#define RED_AIRWAY_H_

#include <blitz/array.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


// forward declarations
struct _MATERIAL;


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class RedAirwayImplInterface;
template <DRT::Element::DiscretizationType distype> class AirwayImpl;
class RedAirwayRegister;

/*!
\brief A C++ wrapper for the airway element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirway : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayRegister;
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedAirway(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedAirway(const RedAirway& old);

  /*!
  \brief Deep copy this instance of RedAirway and return pointer to the copy

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
  {if (NumNode()==2) return 1;
      else {
	  dserror("Could not determine number of lines");
	  return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
  */
  virtual int NumSurface() const {return -1;}

  /*!
  \brief Return number of volumes of this element (always 1) 
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_RedAirway;}

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
  virtual ~RedAirway();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

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
  //bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate airway element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the airway element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(ParameterList&            params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                vector<int>&              lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other
  std::string type(){return elemType_;}


  //@}

  //! get cross-sectional area
  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  double         getEa(){ return Ea_;}

  //! get wall thickness
  double         gettw(){ return tw_;}

private:

  //! action parameters recognized by airway
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    get_initial_state,
    set_bc
  };


  //! Element Type
  std::string                                     elemType_;

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;

  //! Wall Compliance
  double                                          Ew_;

  //! Air Compliance
  double                                          Ea_;

  //! Wall Thickness
  double                                          tw_;

  //! cross-sectional area
  double                                          A_;

  // internal calculation methods

  // don't want = operator
  RedAirway& operator = (const RedAirway& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;


}; // class RedAirway


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for airway element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirwayRegister : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  RedAirwayRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an RedAirwayRegister

  */
  RedAirwayRegister(const DRT::ELEMENTS::RedAirwayRegister& old);

  /*!
  \brief Deep copy this instance of RedAirwayRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::RedAirwayRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_RedAirwayRegister;}

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
  virtual ~RedAirwayRegister();

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

  The method performs all steps needed to init the airway elements.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class RedAirwayRegister


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef RED_AIRWAY_H_
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_RED_AIRWAYS
