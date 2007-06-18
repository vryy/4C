/*!----------------------------------------------------------------------
\file fluid2.H

\class DRT::Elements::Fluid2

\brief A C++ wrapper for the fluid2 element

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef FLUID2_H
#define FLUID2_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"


// forward declarations
struct _MATERIAL;
struct _FLUID_DATA;
/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class Fluid2Register;

namespace Elements
{
// forward declarations
class Fluid2Line;

/*!
\class Fluid2

\brief A C++ wrapper for the fluid2 element

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2 : public DRT::Element
{
public:

  //! @name Friends
  friend class Fluid2Register;
  friend class Fluid2Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Fluid2(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid2(const Fluid2& old);

  /*!
  \brief Deep copy this instance of Fluid2 and return pointer to the copy

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
  virtual int NumLine()
  {if (NumNode()==9 || NumNode()==8 || NumNode()==4) return 4;
      else if (NumNode()==3 || NumNode()==6) return 3;
      else {
	  dserror("Could not determine number of lines");
	  return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of ptrs to the volume of this element

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid2;}

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
  virtual ~Fluid2();

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

  Evaluate fluid2 element stiffness, mass, internal forces etc

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

  this method evaluates a line Neumann condition on the fluid2 element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other

  bool IsAle() const { return is_ale_; }

  //@}


private:

  //! action parameters recognized by fluid2
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual
  };

  //! no. gaussian points in each spatial dimension
  int                                             ngp_[2];
  //! number of the material law
  int                                             material_;
  //! flag for euler/ale net algorithm
  bool                                            is_ale_;

  //! container containing all kind of fluid2 thisandthat
  DRT::Container                                  data_;

  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surface_;

  //! vector of lines adjacent to this element (fluid2: 3 or 4 lines)
  vector<RefCountPtr<DRT::Elements::Fluid2Line> > lines_;
  //! vector of lines adjacent to this element (fluid2: 3 or 4 lines)
  vector<DRT::Elements::Fluid2Line*>              lineptrs_;


  // internal calculation methods

  // don't want = operator
  Fluid2& operator = (const Fluid2& old);

  // calculating system matrix
  void f2_sys_mat(vector<int>&              lm,
                  vector<double>&           evelnp,
		  vector<double>&           eprenp,
                  vector<double>&           ehist,
                  vector<double>&           edispnp,
                  vector<double>&           egridv,
                  Epetra_SerialDenseMatrix* sys_mat,
                  Epetra_SerialDenseVector* residual,
                  struct _MATERIAL*         material,
		  ParameterList&            params
		  );

  // evaluate element integration points (and storing them in data)
  void f2_integration_points(struct _FLUID_DATA& data);

  // calculate Jacobian matrix and it's determinant
  void f2_jaco(const Epetra_SerialDenseMatrix& xyze,
	       const Epetra_SerialDenseMatrix& deriv,
               Epetra_SerialDenseMatrix& xjm,
	       double* 		det,
               const int 	iel
	       );

  // shape functions and natural derivatives for hexaeder
  void f2_shape_function(vector<double>& funct,
               		Epetra_SerialDenseMatrix& deriv,
               		Epetra_SerialDenseMatrix& deriv2,
               		const double&      r,
               		const double&      s,
               		const int&         iel,
               		const int          icode
          		);

  // calculate global derivatives w.r.t. x,y,z at point r,s,t are
    void f2_gder(Epetra_SerialDenseMatrix& derxy,
		 const Epetra_SerialDenseMatrix& deriv,
                 Epetra_SerialDenseMatrix& xjm,
		 double& det,
                 const int iel
		 );

  // calculate second global derivatives w.r.t x/y/z at point r,s,t
    void f2_gder2(const Epetra_SerialDenseMatrix& xyze,
				    const Epetra_SerialDenseMatrix& xjm,
                                    const Epetra_SerialDenseMatrix& derxy,
				    Epetra_SerialDenseMatrix& derxy2,
				    const Epetra_SerialDenseMatrix& deriv2,
                                    const int iel
				    );

  // evaluate fluid coefficient matrix
    void f2_calmat( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  	eforce,
                vector<double>&  	velint,
                vector<double>&   	histvec,
                vector<double>&   	gridvint,
                double&   		press,
                Epetra_SerialDenseMatrix& vderxy,
                Epetra_SerialDenseMatrix& vderxy2,
                vector<double>&   	gradp,
                vector<double>&   	funct,
		vector<double>&   	tau,
                Epetra_SerialDenseMatrix& derxy,
                Epetra_SerialDenseMatrix& derxy2,
                double&   		fac,
                const double&   	visc,
                const int&      	iel,
		ParameterList& 		params
      		);

  // evaluate fluid coefficient matrix for stationary problems
    void f2_calmat_stationary( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  	eforce,
                vector<double>&  	velint,
                vector<double>&   	histvec,
                vector<double>&   	gridvint,
                double&   		press,
                Epetra_SerialDenseMatrix& vderxy,
                Epetra_SerialDenseMatrix& vderxy2,
                vector<double>&   	gradp,
                vector<double>&   	funct,
		vector<double>&   	tau,
                Epetra_SerialDenseMatrix& derxy,
                Epetra_SerialDenseMatrix& derxy2,
                double&   		fac,
                const double&   	visc,
                const int&      	iel,
		ParameterList& 		params
      		);

}; // class Fluid2


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Fluid2Register

\brief A register for fluid2 element

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Fluid2Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Fluid2Register

  */
  Fluid2Register(const DRT::Elements::Fluid2Register& old);

  /*!
  \brief Deep copy this instance of Fluid2Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Fluid2Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid2Register;}

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
  virtual ~Fluid2Register();

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

  The method performs all steps needed to init the fluid2 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Fluid2Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Fluid2Line

\brief An element representing a line of a fluid2 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent fluid2 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gamnitzer (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2Line : public DRT::Element
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
  \param parent: The parent fluid element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Fluid2Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::Fluid2* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid2Line(const Fluid2Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid2Line;}

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
  virtual ~Fluid2Line();

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
  virtual int NumDofPerNode(const DRT::Node& node) const
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

  this method evaluates a line Neumann condition on the fluid2 element

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
  Fluid2Line& operator = (const Fluid2Line& old);

  //! The parent element of this line
  DRT::Elements::Fluid2* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

}; // class Fluid2Line




} // namespace Elements
} // namespace DRT




#endif  // #ifndef FLUID2_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
