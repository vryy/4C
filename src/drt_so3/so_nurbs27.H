/*!----------------------------------------------------------------------
\file so_nurbs27.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SONURBS27_H
#define SONURBS27_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

//@=
//#include "so_surface.H"
//#include "so_line.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_mat/material.H"

#include "../drt_lib/linalg_fixedsizematrix.H"

/// Several parameters which are fixed for Solid nurbs27
const int NUMNOD_SONURBS27 = 27;      ///< number of nodes
const int NODDOF_SONURBS27 = 3;      ///< number of dofs per node
const int NUMDOF_SONURBS27 = 81;     ///< total dofs per element
const int NUMSTR_SONURBS27 = 6;      ///< number of strains/stresses
const int NUMGPT_SONURBS27 = 27;      ///< total gauss points per element
const int NUMDIM_SONURBS27 = 3;      ///< number of dimensions


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief NURBS: namespace of all isogeometric analysis tools

*/
namespace NURBS
{  


// forward declarations
class Sonurbs27Register;

/*!
\brief A C++ version of a 27-control point nurbs solid element

A structural 27-control points nurbs solid displacement element for large deformations. 
As its discretization is fixed many data structures are evaluated just once and kept
for performance.

\author gammi
*/
class So_nurbs27 : public DRT::Element
{
public:

  //! @name Friends
  friend class Sonurbs27Register;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */
  So_nurbs27(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_nurbs27(const So_nurbs27& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of volumes of this element
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 6;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 12;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volumes of this element

  */
  virtual vector<RCP<DRT::Element> > Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_So_Nurbs27;}

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
  virtual ~So_nurbs27();

  //@}

  //! @name Access methods


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

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */

  virtual bool ReadElement();

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate So_nurbs27 element stiffness, mass, internal forces, etc.

  If NULL on input, the controling method does not expect the element 
  to fill these matrices or vectors.
  
  \return 0 if successful, negative otherwise
  */
  int Evaluate
  (
    ParameterList&            params        , ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization, ///< pointer to discretization for de-assembly
    vector<int>&              lm            , ///< location matrix for de-assembly
    Epetra_SerialDenseMatrix& elemat1       , ///< (stiffness-)matrix to be filled by element. 
    Epetra_SerialDenseMatrix& elemat2       , ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1       , ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2       , ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3         ///< vector to be filled by element
  );


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surface Neumann condition on the solid3 element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann
  (
    ParameterList&            params        ,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition     ,
    vector<int>&              lm            ,
    Epetra_SerialDenseVector& elevec1
  );

  //@}

protected:

  //! action parameters recognized by So_nurbs27
  enum ActionType
  {
    none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass
  };

  //! kinematic types of So_nurbs27; only total-lagrangian 
  //! is implemented
  enum KinematicType
  {
    sonurbs27_totlag
  };

   //! kinematic type
  KinematicType                 kintype_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                data_;
  //! vector of inverses of the jacobian in material frame
  vector<LINALG::Matrix<3,3> >  invJ_;
  //! determinant of Jacobian in material frame
  vector<double>                detJ_;


  // internal calculation methods

  //! don't want = operator
  So_nurbs27& operator = (const So_nurbs27& old);

  //! init the inverse of the jacobian and its determinant in 
  //! the material configuration
  virtual void InitJacobianMapping(DRT::Discretization& dis);

  //! Calculate nonlinear stiffness and mass matrix
  virtual void sonurbs27_nlnstiffmass(
      vector<int>&           lm            , ///< location matrix
      DRT::Discretization&   discretization, ///< discretisation to extract knot vector
      vector<double>&        disp          , ///< current displacements
      vector<double>&        residual      , ///< current residual displ
      LINALG::Matrix<81,81>* stiffmatrix   , ///< element stiffness matrix
      LINALG::Matrix<81,81>* massmatrix    , ///< element mass matrix
      LINALG::Matrix<81,1>*  force         , ///< element internal force vector
      ParameterList&         params        );

  /*!
   \brief Interface to material laws 
   
   Here the interface to any material takes place. Stress and C-mat must be retrieved
   whereas all necessary data is passed. Add whatever your material needs, but
   make sure that exchange is not overdone performance-wise.

   \param *stress (out): Voigt-Vector of stresses at current gp
   \param *cmat (out): Elasticity matrix at current gp
   \param *density (out): density of material
   \param *glstrain (in): Voigt-Vector of GL-strains at current gp
   \param *defgrd (in): F at current gp, CAUTION! only for disp-based soh8
   \param gp (in): current gp
   \param params (in): parameterlist to access time, etc. in materials
   */

  void sonurbs27_mat_sel(
    LINALG::Matrix<6,1>* stress  ,
    LINALG::Matrix<6,6>* cmat    ,
    double*              density ,
    LINALG::Matrix<6,1>* glstrain,
    LINALG::Matrix<3,3>* defgrd  ,
    const int            gp      ,
    ParameterList&       params  );


  //! Evaluate Nurbs27 Shapefcts to keep them static
  const vector<LINALG::Matrix<27, 1> > sonurbs27_shapefcts(
    const std::vector<Epetra_SerialDenseVector> & myknots,
    const LINALG::Matrix<27,1>                  & weights
    );
  //! Evaluate Nurbs27 Derivs to keep them static
  const vector<LINALG::Matrix< 3,27> > sonurbs27_derivs(
    const std::vector<Epetra_SerialDenseVector> & myknots,
    const LINALG::Matrix<27,1>                  & weights
    );
  //! Evaluate Nurbs27 Weights to keep them static
  const vector<double> sonurbs27_gpweights();


}; // class So_nurbs27


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief Just another register 

\author gammi (gammi@lnm.mw.tum.de)
*/
class Sonurbs27Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Sonurbs27Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Sonurbs27Register

  */
  Sonurbs27Register(const DRT::ELEMENTS::NURBS::Sonurbs27Register& old);

  /*!
  \brief Deep copy this instance of Sonurbs27Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::NURBS::Sonurbs27Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return -1;}

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
  virtual ~Sonurbs27Register();

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

  The method performs all steps needed to init the elements.

  */
  int Initialize(DRT::Discretization& dis);

  //@}


private:


}; // class Sonurbs27Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace NURBS
} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SONURBS27_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
