/*!----------------------------------------------------------------------
\file artery.H

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ARTNET
#ifdef CCADISCRET
#ifndef ARTERY_H
#define ARTERY_H

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
class ArteryExpInterface;
template <DRT::Element::DiscretizationType distype> class ArteryLinExp;
class ArteryRegister;

/*!
\brief A C++ wrapper for the artery element

\author gammi (ismail@lnm.mw.tum.de)
*/
class Artery : public DRT::Element
{
public:

  //! @name Friends
  friend class ArteryRegister;
  friend class ArteryExpInterface;
  friend class ArteryLinExp<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Artery(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Artery(const Artery& old);

  /*!
  \brief Deep copy this instance of Artery and return pointer to the copy

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
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 0;}

  /*!
  \brief Return number of volumes of this element (always 1) even though it's a 1D element
  */
  virtual int NumVolume() const {return 0;}

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Artery;}

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
  virtual ~Artery();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 2;}

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
  bool ReadElement();

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate artery element stiffness, mass, internal forces etc

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

  this method evaluates a line Neumann condition on the artery element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(ParameterList&            params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                vector<int>&              lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other

  bool IsAle() const { return is_ale_; }

  //@}


private:

  //! action parameters recognized by Artery
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_stiffness_matrix,
    calc_mass_matrix,
    calc_charcteristic_variables,
    get_density,
    get_viscosity,
    get_stiffness_matrix,
    get_mass_matrix,
    get_beta,
    get_sound_speed,
    get_initail_artery_state,
    set_term_bc
  };


  //! Gaussrule
  DRT::UTILS::GaussRule1D                         gaussrule_;

  //! flag for euler/ale net algorithm
  bool                                            is_ale_;

  //! container containing all kind of Artery thisandthat
  DRT::Container                                  data_;


  // internal calculation methods

  // don't want = operator
  Artery& operator = (const Artery& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;


}; // class Artery


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for artery element

\author ismail (ismail@lnm.mw.tum.de)
*/
class ArteryRegister : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  ArteryRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an ArteryRegister

  */
  ArteryRegister(const DRT::ELEMENTS::ArteryRegister& old);

  /*!
  \brief Deep copy this instance of ArteryRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::ArteryRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ArteryRegister;}

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
  virtual ~ArteryRegister();

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

  The method performs all steps needed to init the artery elements. 
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class ArteryRegister


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef ARTERY_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ARTNET
