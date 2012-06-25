/*
 * so_test4_volcoupl.H
 *
 *  Created on: Jun 22, 2012
 *      Author: bertoglio
 */

#ifndef SO_TEST4_VOLCOUPL_H_
#define SO_TEST4_VOLCOUPL_H_

#include "so_tet4.H"
#include "so3_poro.H"
#include "so3_poro_eletypes.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief A C++ version of a 4-node tet solid element with modifications for porous media

A structural 4-node hexahedral solid displacement element for large deformations
and (near)-incompressibility.

*/

template <DRT::Element::DiscretizationType distype> class coupltype;
template< class coupltype>
class So_tet4_volcoupl : public So_tet4,
                     public coupltype
{
public:

  //! @name Friends
  friend class So_tet4PoroType;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */

  So_tet4_volcoupl(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_tet4_volcoupl(const So_tet4_volcoupl<coupltype> & old);

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
  virtual int UniqueParObjectId() const { return So_tet4PoroType::Instance().UniqueParObjectId(); }

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
  virtual ~So_tet4_volcoupl();

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


  virtual DRT::ElementType & ElementType() const
  { return So_tet4PoroType::Instance(); }

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
  So_tet4_volcoupl& operator = (const So_tet4_volcoupl& old);

  //! special action parameters for volume coupled solid elements
  enum ActionType
  {
    none,
    calc_struct_multidofsetcoupling,//!< structure-fluid coupling: internal force, stiffness for other dofset (structural part)
  };

}; // class So_tet4_volcoupl

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

} // namespace ELEMENTS
} // namespace DRT


#endif /* SO_TEST4_VOLCOUPL_H_ */
