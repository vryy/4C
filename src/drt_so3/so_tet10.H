/*!----------------------------------------------------------------------
\file so_tet10.H

\class DRT::ELEMENTS::So_tet10

\brief A C++ version of the 10-node tet solid element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by: Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_TET10_H_
#define SO_TET10_H_


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

#include "../drt_mat/material.H"

#define NUMNOD_SOTET10      10      // number of nodes
#define NODDOF_SOTET10       3      // number of dofs per node
#define NUMDOF_SOTET10      30      // total dofs per element
#define NUMSTR_SOTET10       6      // number of strains/stresses 	/****/
#define NUMGPT_SOTET10       4      // total gauss points per element	/****/
#define NUMDIM_SOTET10       3      // number of dimensions 			/****/
#define NUMCOORD_SOTET10	 4 		// number of shape function cooordinates (ksi1-ksi4) 
#define NUMNOD_SOTET10_FACE  6		// number of nodes on a TET10 face (which is a TRI6) 
#define NUMGPT_SOTET10_FACE	 3		// number of GP    on a TET10 face (which is a TRI6) 

// forward declarations
struct _MATERIAL;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Sotet10Volume;
class Sotet10Surface;
class Sotet10Line;
class Sotet10Register;

/*!
\class So_tet10

\brief A C++ version of the 10-node tet solid element

\author Alexander Volf (alexander.volf@mytum.de)
*/

class So_tet10 : public DRT::Element
{
public:

  //! @name Friends
  friend class Sotet10Register;
  friend class Sotet10Surface;
  friend class Sotet10Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_tet10(int id, int owner);
  
  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_tet10(const So_tet10& old);
  
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
  virtual int NumSurface() const {return 4;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 6;}

  /*!
  \brief Get vector of ptrs to the lines of this element

  */
  virtual DRT::Element** Lines();

  /*!
  \brief Get vector of ptrs to the surfaces of this element

  */
  virtual DRT::Element** Surfaces();

  /*!
  \brief Get vector of ptrs to the volumes of this element

  */
  virtual DRT::Element** Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_So_tet10;}

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
  virtual ~So_tet10();
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

  Evaluate so_tet10 element stiffness, mass, internal forces, etc.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                          the controling method does not expect the element to fill
                          this matrix.
  \param elemat2 (out)  : (mass-)matrix to be filled by element. If NULL on input,
                          the controling method does not expect the element to fill
                          this matrix.
  \param elevec1 (out)  : (internal force-)vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not expect the element
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

  this method evaluates a surface Neumann condition on the solid3 element

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

protected:

  //! action parameters recognized by so_tet10
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
    calc_struct_fsiload,
    calc_struct_update_istep
  };

  //! kinematic types of so_tet10
  enum KinematicType
  {
    so_tet10_geolin,
    so_tet10_totlag,
    so_tet10_updlag
  };
  enum StressType//###################
  {
    so_tet10_stress_none,
    so_tet10_stress_gpxyz,
    so_tet10_stress_gprst,
    so_tet10_stress_gp123,
    so_tet10_stress_ndxyz,
    so_tet10_stress_ndrst,
    so_tet10_stress_nd123
  };

  //! kinematic type
  KinematicType                                   kintype_;
  //! type of stress evaluation
  StressType                                      stresstype_;
  //! no. gaussian points in each spatial dimension
  int                                             ngp_[3];
  //! number of the material law
  int                                             material_;

  DRT::Container                                  data_;
  //! vector of volume (length 1) of this element
  vector<RefCountPtr<DRT::ELEMENTS::Sotet10Volume> > volume_;
  //! vector of surfaces of this element
    vector<RefCountPtr<DRT::ELEMENTS::Sotet10Surface> > surfaces_;
  //! vector of surfaces of this element
  vector<DRT::ELEMENTS::Sotet10Surface*>           surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::ELEMENTS::Sotet10Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::ELEMENTS::Sotet10Line*>              lineptrs_;
  //! vector of history variables for each Gauss point
  vector<RefCountPtr<MAT::Material> >                mat_;
  
  
  //! some distypes need rewinding of NodeIds due to bad input from mesher
  //! flag whether process of rewinding is done
  bool                                            donerewinding_;

  // internal calculation methods

  // don't want = operator
  So_tet10& operator = (const So_tet10& old);
  virtual void so_tet10_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 struct _MATERIAL*         material);

  void so_tet10_stress(struct _MATERIAL* material,
                   vector<double>& mydisp,
                   Epetra_SerialDenseMatrix* stresses);
  
/*!
\class DRT::ELEMENTS::So_tet10::Discrete_integrator

\brief A general class agregating all necessery data for discrete integration

<h3>Purpose</h3>
The class offers a standard interface for discrete integration, mainly the agregated
shape functions and shape function derivatives, it is a purely general class and 
does NOT own cnstructor. For use constructors of derived classes. It has however  
a destructor (not virtual!(to save on vtable* pointers)), that should be the same 
for all derived classes. Derived classes are not supposed to add extra attriutes,
 they may however add extra methods.

pointer to (static) array of shape function stored in Epetra_Vectors			
evaluated at each gp , for each node it may look like this:

shapefct[gp] = [N_1(gp)  N_2(gp)  -- N_X(gp)]

pointer to (static) array of shape function derivatives stored in
Epetra_Matrices evaluated at each gp, for xsi1-4(or 1-3) it may look like this:	
													    |
               [  N1_,xsi1   N1_,xsi2    N1_,xsi3   N1_,xsi4 ]
deriv_gp[gp] = [  N2_,xsi1   N2_,xsi2    N2_,xsi3   N2_,xsi4 ] (at a given gp)
               [    |           |           |          |     ]
               [  NX_,xsi1   NX_,xsi2    NX_,xsi3   NX_,xsi4 ]
               
NOTE: X stands for the number of NODES (and shape functions) in a element


\author a.vlf
*/

  class Discrete_integrator
  {
	public: 
		Epetra_SerialDenseVector* shapefct_gp;  
      	Epetra_SerialDenseMatrix* deriv_gp;  
      	Epetra_SerialDenseVector weights;
      	int num_gp;
      	int num_nodes;
      	int num_coords;
      	~Discrete_integrator();
  };
 /*---------------------------------------------------------------------------------*
  | class Tet_integrator_4point											vlf 09/07   |
  |																					|
  | a special class for quadratic 10-node tetrahedra using a 4 point quadrature   	|
  |	it is exaxt for linear and quadratic polynomials								|
  |																					|
  *---------------------------------------------------------------------------------*/
  class Tet_integrator_4point : public Discrete_integrator
  {
  	public:
  		Tet_integrator_4point(void);
  };
 
  /*---------------------------------------------------------------------------------*
  | class Tet_integrator_14point										vlf 09/07   |
  |																					|
  | a special class for quadratic 10-node tetrahedra using a 14 point quadrature   	|
  |	it is exaxt for 0-4th order polynomials											|
  |																					|
  *---------------------------------------------------------------------------------*/
  
  class Tet_integrator_14point : public Discrete_integrator
  {
  	public:
  		Tet_integrator_14point(void);
  };
  
  // became obsolete
  void so_tet10_shapederiv(
  	  Epetra_SerialDenseVector** shapefct_ptr,  
      Epetra_SerialDenseMatrix** deriv_gp_ptr,  
      Epetra_SerialDenseVector** weights);       

  // #### hier deklariert, in implementierung soll einfach soh8_mat_sel gerufen werden
  void so_tet10_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    int gp);       

}; // class So_tet10

/*!
\class Sotet10Register

\brief A register for tet10 element

\author gee (gee@lnm.mw.tum.de)
*/

class Sotet10Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Sotet10Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Sotet10Register

  */
  Sotet10Register(const DRT::ELEMENTS::Sotet10Register& old);

  /*!
  \brief Deep copy this instance of Sotet10Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Sotet10Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Sotet10Register;}

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
  virtual ~Sotet10Register();

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

  The method performs all steps needed to init the shell8 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Sotet10Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Sotet10Surface

\brief An element representing a surface of a So_tet10 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author maf (frenzel@lnm.mw.tum.de)
*/

class Sotet10Surface : public DRT::Element
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this surface
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent fluid element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element
  */
  Sotet10Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::So_tet10* parent, const int lsurface);


  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Sotet10Surface(const Sotet10Surface& old);
  
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
  virtual int UniqueParObjectId() const { return ParObject_Sotet10Surface;}

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
  virtual ~Sotet10Surface();

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

  this method evaluates a surface Neumann condition on the fluid3 element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList&            params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //@}

private:

  // don't want = operator
  Sotet10Surface& operator = (const Sotet10Surface& old);

  //! The parent element of this surface
  DRT::ELEMENTS::So_tet10* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;

  void sotet10_surface_shapefunc(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseVector** weights);
  
  /*void sotet10_surface_integ(	//became obsolete !!!!!!!!!!!!!!!!
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      const Epetra_SerialDenseMatrix* xsrefe,// (i) material element coords
      const double r,                        // (i) coord in r-direction
      const double s);                       // (i) coord in s-direction*/
}; // class Sotet10Surface

/*!
\class Sotet10Line

\brief An element representing a line edge of a So_tet10 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent So_tet10 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/

class Sotet10Line : public DRT::Element
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
  \param parent: The parent shell element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Sotet10Line(int id, int owner, int nnode, const int* nodeids,
               DRT::Node** nodes, DRT::ELEMENTS::So_tet10* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Sotet10Line(const Sotet10Line& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Sotet10Line;}
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
  virtual ~Sotet10Line();

  //@}

  //! @name Acess methods

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

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

  this method evaluates a line Neumann condition on the shell element

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
  Sotet10Line& operator = (const Sotet10Line& old);
  //! The parent element of this line
  DRT::ELEMENTS::So_tet10* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

}; // class Sotet10Line




} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef SO_TET10_H_
#endif  // #ifdef  CCADISCRET
#endif  // #ifdef  D_SOLID3
