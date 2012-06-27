/*!----------------------------------------------------------------------
\file So3_volcoupl.H

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/


#ifndef SO3_VOLCOUPL_H_
#define SO3_VOLCOUPL_H_

#include "so3_poro.H"
#include "so3_scatra.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief A C++ version of a general solid element with volumetric coupling

*/

template <DRT::Element::DiscretizationType distype> class coupltype;
template<class so3_ele, class coupltype>
class So3_volcoupl : public so3_ele,
                     public coupltype
{
public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */

  So3_volcoupl(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So3_volcoupl(const So3_volcoupl<so3_ele,coupltype> & old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const;

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~So3_volcoupl();

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const ;

  //@}

  //! @name Access methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;


  virtual DRT::ElementType & ElementType() const;

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  virtual bool ReadElement(const std::string& eletype,
                   const std::string& eledistype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.

  If NULL on input, the controlling method does not expect the element
  to fill these matrices or vectors.

  \return 0 if successful, negative otherwise
  */
  int Evaluate
  (
    ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );


protected:

  //! don't want = operator
  So3_volcoupl& operator = (const So3_volcoupl& old);

  //! special action parameters for volume coupled solid elements
  enum ActionType
  {
    none,
    calc_struct_multidofsetcoupling,//!< structure-fluid coupling: internal force, stiffness for other dofset (structural part)
  };

}; // class So3_volcoupl

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

} // namespace ELEMENTS
} // namespace DRT


#endif /* SO3_VOLCOUPL_H_ */
