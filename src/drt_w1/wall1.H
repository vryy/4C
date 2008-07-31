/*----------------------------------------------------------------------*/
/*!
\file wall1.H

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifndef WALL1_H
#define WALL1_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"

// element constants
#ifndef MAXNOD_WALL1
#define MAXNOD_WALL1 (9)  /* maximum of element nodes : max quad9 */
#endif

// parameters which are fixed for Wall1
const int NUMDIM_W1 = 2;   ///< number of dimensions
const int NODDOF_W1 = 2;   ///< number od DOFs at each node
const int NUMSTR_W1 = 3;   ///< number of strains/stresses
const int neas = 4;        ///< number of eas parameters

// forward declarations
struct _MATERIAL;
struct _W1_DATA;


namespace DRT
{
// forward declarations
class Discretization;
class Wall1Register;

namespace ELEMENTS
{
// forward declarations
class Wall1Line;

/*!
\brief A C++ wrapper for the wall1 element

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class Wall1 : public DRT::Element
{
public:

  //! @name Friends
  friend class Wall1Register;
  friend class Wall1Line;

  //@}
  /// @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Wall1(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Wall1(const Wall1& old);

  /*!
  \brief Deep copy this instance of Wall1 and return pointer to the copy

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
  {if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4; else return 3;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Wall1;}

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
  virtual ~Wall1();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return NODDOF_W1;}

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
  \brief Allocate and return an ElementRegister class

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

  Evaluate Wall1 element stiffness, mass, internal forces etc

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

  this method evaluates a surfaces Neumann condition on the wall element

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

  //! action parameters recognized by wall1
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep  //!< reset elementwise internal variables
                             //!< during iteration to last converged state
  };

  /// 
  enum DimensionalReduction
  {
    plane_none,
    plane_stress,
    plane_strain
  };

  //! type of stresses to be calculatecd in postprocessing
  enum StressType
  {
    w1_none,
    w1_xy,
    w1_rs
  };

//  //! couptyp for thermo-structure-interaction
//  enum ThermoStructureType
//  {
//    tsi_coup_none,
//    tsi_coup_thermconf,
//    tsi_coup_thermcreate
//  };

  /// container containing all kind of wall1 thisandthat
  DRT::Container data_;

  /// number of the material law
  int material_; 
  /// element thickness 
  double thickness_; 
  /// gaussian points
  DRT::UTILS::GaussRule2D gaussrule_;
  /// problem type
  DimensionalReduction wtype_;
  /// type of stress calculation
  StressType stresstype_;
  /// eas or not
  bool iseas_;

//  //! thermo-structure-interaction couptype
//  ThermoStructureType                         tsi_couptyp_;


  /// @name Internal calculation methods
  //@{

  /// evaluate the element forces and stiffness and mass
  /// \author mgit \date 03/07
  void w1_nlnstiffmass
  (
    const vector<int>& lm,  ///< location vector
    const vector<double>& disp,  ///< element displacements
    const vector<double>& residual,  ///< residual displacements 
    Epetra_SerialDenseMatrix* stiffmatrix,  ///< element stiffness matrix
    Epetra_SerialDenseMatrix* massmatrix,  ///< element mass matrix
    Epetra_SerialDenseVector* force,  ///< element internal force vector
    Epetra_SerialDenseMatrix* elestress,  ///< element stresses
    Epetra_SerialDenseMatrix* elestrain,  ///< element strains
    struct _MATERIAL* material,  ///< element material
    const bool cauchy = false  ///< if true calculate Cauchy instead of PK2 stresses
  );

 // void w1_integration_points(struct _W1_DATA& data);

 // void w1_gpdom(int& totngp,
 //               double** gpcr, double** gpcs,
 //               double** gpwr);

//  void w1_shapefunctions(vector<double>& funct, Epetra_SerialDenseMatrix& deriv,
//                               const double r, const double s, const int numnode,
//                               const int doderiv) const;

  /// Jacobian matrix for mapping from parameter space in physical material space
  /// at point parameter space
  /// \author mgit \date 04/07
  void w1_jacobianmatrix
  (
    double xrefe[2][MAXNOD_WALL1],  ///< reference/material co-ordinates of element nodes
    const Epetra_SerialDenseMatrix& deriv,  ///< derivatives of shape functions at parameter point
    Epetra_SerialDenseMatrix& xjm,  ///< Jacobi matrix
    double* det,  ///< determinant of Jacobi matrix
    const int iel  ///< actual number of element nodes
  );

  /// Linear B-operator in reference configuration at point parameter space
  /// \author mgit \date 04/07
  void w1_boplin
  (
    Epetra_SerialDenseMatrix& boplin,  ///< the B-operator
    Epetra_SerialDenseMatrix& deriv,  ///< derivatives of shape functions at parameter point
    Epetra_SerialDenseMatrix& xjm,  ///< Jacobian at parameter point
    double& det,  ///< Jacobian determinant at parameter point
    const int iel  ///< number of element nodes
  );

  /// (Material) Deformation gradient \f$F\f$ and Green-Lagrange strains \f$E\f$
  /// at parameter point
  /// \author mgit \date 04/07
  void w1_defgrad
  (
    Epetra_SerialDenseVector& F,  ///< deformation gradient
    Epetra_SerialDenseVector& strain,  ///< GL strain \f$E^T=[E_{11} \; E_{22} \; E_{12}]\f$
    const double xrefe[][MAXNOD_WALL1],  ///< reference/material co-ordinates of element nodes
    const double xcure[][MAXNOD_WALL1], ///< current/spatial co-ordinates of element nodes
    Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const int iel  ///< number of element nodes
  );

  /// Non-linear B-operator in reference configuration
  /// \author mgit \date 04/07
  void w1_boplin_cure
  (
    Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-operator
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const Epetra_SerialDenseVector& F,  ///< deformation gradient as Voigt-vector
    const int numeps,  ///< number of strains
    const int nd  ///< number of element nodes
  );

  /// Constitutive matrix \f$C\f$ and stresses
  /// \author mgit \date 05/07
  void w1_call_matgeononl
  (
    const Epetra_SerialDenseVector& strain,  ///< Green-Lagrange strain vector
    Epetra_SerialDenseMatrix& stress,  ///< stress vector
    Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
    const int numeps,  ///< number of strains
    const struct _MATERIAL* material  ///< the material data
  );

  /// Geometric stiffness constribution (total Lagrange) 
  /// \author mgit \date 05/07
  void w1_kg
  (
    Epetra_SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix 
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    const double fac,  ///< integration factor
    const int nd,  ///< number of element DOFs
    const int numeps  ///< number of strains
  );
  
  /// elastic and initial displacement stiffness contribution (total Lagrange)
  /// \author mgit \date 05/07
  void w1_keu
  (
    Epetra_SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix 
    const Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-operator
    const Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
    const double fac,  ///< integration factor
    const int nd,  ///< number of element DOFs
    const int numeps  ///< number of strains
  );

  /// Evaluate internal element forces for large def (total Lagr)
  /// \author mgit \date 05/07
  void w1_fint
  (
    const Epetra_SerialDenseMatrix& stress, ///< PK2 stress vector
    const Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-op
    Epetra_SerialDenseVector& intforce,  ///< internal force vector
    const double fac,  ///< integration factor
    const int nd  ///< number of element DOFs
  );

  /// setup of constant EAS data
  void w1_eassetup
  (
    Epetra_SerialDenseMatrix& boplin0,  ///< linear B-op at origin
    Epetra_SerialDenseVector& F0,  ///< deformation gradient at origin
    Epetra_SerialDenseMatrix& xjm0,  ///< jacobian matrix at origin
    double& detJ0,  ///< det of Jacobian at origin
    double xrefe[2][MAXNOD_WALL1],  ///< material element coords
    double xcure[2][MAXNOD_WALL1],  ///< current element coords
    const DRT::Element::DiscretizationType& distype  ///< discretisation type
  );

  /// Get the enhanced deformation gradient and
  /// also the operators G,W0 and Z
  /// at point in parameter space
  /// \author mgit \date 01/08
  void w1_call_defgrad_enh
  (
    Epetra_SerialDenseMatrix& F_enh,  ///< enhanced deformation gradient
    const Epetra_SerialDenseMatrix xjm0,  ///< Jacobian at origin
    const Epetra_SerialDenseMatrix xjm,  ///< Jacobian at parameter point
    const double detJ0,  ///< det of Jacobian at origin
    const double det,  ///< det of Jacobian at parameter point
    const Epetra_SerialDenseVector F0,  ///< deformation gradient at origin
    const Epetra_SerialDenseMatrix alpha,  ///< alpha parameters (EAS params)
    const double e1,  ///< \f$\xi\f$ co-ordinate of parameter point
    const double e2,  ///< \f$\eta\f$ co-ordinate of parameter point
    Epetra_SerialDenseMatrix& G,  ///< G-operator
    Epetra_SerialDenseMatrix& W0,  ///< W_0-operator
    const Epetra_SerialDenseMatrix boplin0,  ///< linear B-op at origin
    Epetra_SerialDenseMatrix& Z  ///< Z-operator
  ); 

  /// Total deformation gradient and (enhanced) Green-Lagrange strain
  /// \author mgit \date 01/08
  void w1_call_defgrad_tot
  (
    const Epetra_SerialDenseMatrix F_enh,  ///< enhanced def.grad.
    Epetra_SerialDenseMatrix& F_tot,  ///< total def.grad.
    const Epetra_SerialDenseVector F,  ///< displ.-based def.grad.
    Epetra_SerialDenseVector& strain  ///< GL strains
  ); 
  
  /// first Piola-Kirchhoff stress vector  
  /// \author mgit \author 02/08
  void w1_stress_eas
  (
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    const Epetra_SerialDenseMatrix& F_tot,  ///< total def.grad.
    Epetra_SerialDenseMatrix& p_stress  ///< PK1 stress vector
  ); 
  
  /// calculate stiffness matrix kdd \f$\partial f_{int}/\partial d\f$
  /// \author mgit \date 03/08
  void w1_kdd(
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator at origin
    const Epetra_SerialDenseMatrix& F,  ///< total def.grad.
    const Epetra_SerialDenseMatrix& C,  ///< consititutive matrix
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    Epetra_SerialDenseMatrix& FCF,  ///< \f$F^T \dot C \dot F\f$
    Epetra_SerialDenseMatrix& estif,  ///< element stiff matrix
    const double fac  ///< integration factor
  );
  
  /// calculate tangential matrix kda \f$\partial f_{int}/\partial \alpha\f$
  /// \author mgit \date 03/08
  void w1_kda
  (
    const Epetra_SerialDenseMatrix& FCF,  ///< a product
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator at origin
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress
    const Epetra_SerialDenseMatrix& G,  ///< G-operator
    const Epetra_SerialDenseMatrix& Z,  ///< Z-operator
    Epetra_SerialDenseMatrix& Kda,  ///< target: kda
    const Epetra_SerialDenseMatrix& p_stress,  ///< PK1 stress
    const double fac  ///< integration factor
  );
  
  /// calculate tangential matrix kaa \f$\partial s/\partial \alpha\f$
  /// \author mgit \date 03/08
  void w1_kaa
  (
    const Epetra_SerialDenseMatrix& FCF,  ///< a product \f$F^T \dot C \dot F\f$
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress
    const Epetra_SerialDenseMatrix& G,  ///< G-op
    Epetra_SerialDenseMatrix& Kaa,  ///< target: kaa
    const double fac  ///< integration factor
  );
  
  /// calculate internal forces fint(displacements u) and feas
  void w1_fint_eas
  (
    const Epetra_SerialDenseMatrix& W0,  ///< W-op at origin
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& G,  ///< G-op
    const Epetra_SerialDenseMatrix& p_stress,  ///< PK1 stress
    Epetra_SerialDenseVector& intforce,  ///< internal force vector
    Epetra_SerialDenseVector& feas,  ///< internal EAS constraint condition
    const double fac  ///< integration factor
  );

  /// lump mass matrix
  void w1_lumpmass
  (
    Epetra_SerialDenseMatrix* emass  ///< (in/out) element mass matrix
  );

  /// don't want = operator
  Wall1& operator = (const Wall1& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getGaussrule
  (
    int* ngp  ///< number of Gauss points
  );

  //@}


}; // class Wall1


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for Wall1 element

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class Wall1Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Wall1Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Wall1Register

  */
  Wall1Register(const DRT::ELEMENTS::Wall1Register& old);

  /*!
  \brief Deep copy this instance of Wall1Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Wall1Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Wall1Register;}

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
  virtual ~Wall1Register();

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

  The method performs all steps needed to init the Wall1 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Wall1Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a line edge of a Wall1 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent wall1 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class Wall1Line : public DRT::Element
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
  \param parent: The parent wall element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Wall1Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Wall1* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Wall1Line(const Wall1Line& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;


  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Wall1Line;}

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
  virtual ~Wall1Line();

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
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return 2;}

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
  \brief Allocate and return a ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const { return null; }

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the wall element

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
  
  /*!
    \brief Evaluate a boundary condition

    this method evaluates a line Neumann condition on the wall element

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param condition (in)     : The condition to be evaluated
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,

    \return 0 if successful, negative otherwise
    */
    virtual int Evaluate(ParameterList& params,
                                DRT::Discretization&      discretization,
                                vector<int>&              lm,
                                Epetra_SerialDenseMatrix& elematrix1,
                                Epetra_SerialDenseMatrix& elematrix2,
                                Epetra_SerialDenseVector& elevector1,
                                Epetra_SerialDenseVector& elevector2,
                                Epetra_SerialDenseVector& elevector3);

  //@}
    
    
private:

  enum ActionType
  {
    none,
    calc_struct_areaconstrstiff,
    calc_struct_constrarea
  };
  
  // don't want = operator
  Wall1Line& operator = (const Wall1Line& old);

  //! The parent element of this line
  DRT::ELEMENTS::Wall1* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

  //! Submethod to compute necessary change to stiffness matrix due to the constraint
  void ComputeAreaConstrStiff(Epetra_SerialDenseMatrix xscurr,
      Epetra_SerialDenseMatrix& elematrix);
  
  //! Submethod to compute first derivatives of constraint area w.r.t. the displacements
  void ComputeAreaConstrDeriv(Epetra_SerialDenseMatrix xscurr,
      Epetra_SerialDenseVector& elevector);
  
  // compute infintesimal line element dr for integration along the line
  double  w1_substitution(const Epetra_SerialDenseMatrix& xyze,
		                  const Epetra_SerialDenseMatrix& deriv,
		                  vector<double>* unrm,
		                  const int iel);

  // set number of gauss points to element shape default
   DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

}; // class Wall1Line





} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef WALL1_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
