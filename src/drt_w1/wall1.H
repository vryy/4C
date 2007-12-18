/*!----------------------------------------------------------------------
\file wall1.H

\class DRT::Elements::Wall1

\brief A C++ wrapper for the wall1 element

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifndef WALL1_H
#define WALL1_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"

// element constants
#ifndef MAXNOD_WALL1
#define MAXNOD_WALL1 (9)  /* maximum of element nodes : max quad9 */
#endif

// forward declarations
struct _MATERIAL;
struct _W1_DATA;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class Wall1Register;

namespace Elements
{
// forward declarations
class Wall1Line;

/*!
\class wall1

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
  //! @name Constructors and destructors and related methods

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
  \brief Get vector of ptrs to the lines of this element
  
  */
  virtual DRT::Element** Lines();

  /*!
  \brief Get vector of ptrs to the surfaces of this element
  
  */
  virtual DRT::Element** Surfaces();

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
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep
  };

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

  //! couptyp for thermo-structure-interaction
  enum ThermoStructureType
  {
    tsi_coup_none,
    tsi_coup_thermconf,
    tsi_coup_thermcreate
  }; 

  //! container containing all kind of wall1 thisandthat
  DRT::Container                                  data_;      

  //! number of the material law
  int                                             material_; 
  //! element thickness 
  double                                          thickness_; 
  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surfaces_;  
  //! vector of lines adjacent to this element (wall1: 3 or 4 lines)
  vector<RefCountPtr<DRT::Elements::Wall1Line> >  lines_;     
  //! vector of lines adjacent to this element (wall1: 3 or 4 lines)
  vector<DRT::Elements::Wall1Line*>               lineptrs_;
  //! gaussian points
  DRT::Utils::GaussRule2D                         gaussrule_;
  //! problem type
  DimensionalReduction                            wtype_;
  //! type of stress calculation
  StressType                                   stresstype_;
  //! eas or not
  bool                                              iseas_;
  //! thermo-structure-interaction couptype
  ThermoStructureType                         tsi_couptyp_;
  
  
  // internal calculation methods

  void w1_nlnstiffmass(vector<int>&              lm, 
                       vector<double>&           disp, 
                       vector<double>&           residual,
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force,
                       struct _MATERIAL*         material);

 // void w1_integration_points(struct _W1_DATA& data); 

 // void w1_gpdom(int& totngp, 
 //               double** gpcr, double** gpcs,
 //               double** gpwr);
 
//  void w1_shapefunctions(vector<double>& funct, Epetra_SerialDenseMatrix& deriv,
//                               const double r, const double s, const int numnode, 
//                               const int doderiv) const;

  void w1_jacobianmatrix(double xrefe[2][MAXNOD_WALL1], const Epetra_SerialDenseMatrix& deriv,
                         Epetra_SerialDenseMatrix& xjm,
			 double* det,
                         const int iel);

  void w1_boplin(Epetra_SerialDenseMatrix& boplin,
                 Epetra_SerialDenseMatrix& deriv,
                 Epetra_SerialDenseMatrix& xjm,
                 double& det,  
                 const int iel);

  void w1_defgrad(Epetra_SerialDenseVector& F,
                  Epetra_SerialDenseVector& strain, 
                  const double xrefe[][MAXNOD_WALL1],
                  const double xcure[][MAXNOD_WALL1], 
                  Epetra_SerialDenseMatrix& boplin,
                  const int iel);

  void w1_boplin_cure(Epetra_SerialDenseMatrix& b_cure,
                      Epetra_SerialDenseMatrix& boplin,
                      Epetra_SerialDenseVector& F,
                      int numeps,
                      int nd);  

  void w1_call_matgeononl(Epetra_SerialDenseVector& strain,
                          Epetra_SerialDenseMatrix& stress,
                          Epetra_SerialDenseMatrix& C,                    
                          const int numeps,
                          struct _MATERIAL* material);

  void w1_kg(Epetra_SerialDenseMatrix& estif,
             Epetra_SerialDenseMatrix& boplin,
             Epetra_SerialDenseMatrix& stress,
             double fac,
             int nd,
             int numeps);

  void w1_keu(Epetra_SerialDenseMatrix& estif,
              Epetra_SerialDenseMatrix& b_cure,
              Epetra_SerialDenseMatrix& C,
              double fac,
              int nd,
              int numeps);  
  
  void w1_fint(Epetra_SerialDenseMatrix& stress,
               Epetra_SerialDenseMatrix& b_cure,
               Epetra_SerialDenseVector& intforce,
               double fac,
               int nd);				    

  // don't want = operator
  Wall1& operator = (const Wall1& old);
  
  
  void CreateLinesTri(const int& nline,
                      const int& nnode);
  
  void CreateLinesQuad(const int& nline,
                       const int& nnode);
  
  // set number of gauss points to element shape default
   DRT::Utils::GaussRule2D getGaussrule(int* ngp);
  

}; // class Wall1


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Wall1Register

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
  Wall1Register(const DRT::Elements::Wall1Register& old);

  /*!
  \brief Deep copy this instance of Wall1Register and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::Elements::Wall1Register* Clone() const;

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
\class Wall1Line

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
             DRT::Node** nodes, DRT::Elements::Wall1* parent, const int lline);

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

  //@}

private:
	
  // don't want = operator
  Wall1Line& operator = (const Wall1Line& old);
  
  //! The parent element of this line
  DRT::Elements::Wall1* parent_;  
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;   
  
  // compute infintesimal line element dr for integration along the line
  double  w1_substitution(const Epetra_SerialDenseMatrix& xyze,
		                  const Epetra_SerialDenseMatrix& deriv,
		                  const int iel);
  
  // set number of gauss points to element shape default
   DRT::Utils::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

}; // class Wall1Line





} // namespace Elements
} // namespace DRT

#endif  // #ifndef WALL1_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
