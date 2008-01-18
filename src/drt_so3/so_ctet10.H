/*!----------------------------------------------------------------------
\file so_ctet10.H

\class DRT::Elements::So_ctet10

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
#ifdef D_SOTET
#ifdef CCADISCRET
#ifndef SO_CTET10_H_
#define SO_CTET10_H_

#include "so_integrator.H"
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

#define NUMNOD_SOCTET10      10     // number of nodes
#define NODDOF_SOCTET10       3     // number of dofs per node
#define NUMDOF_SOCTET10      30     // total dofs per element
#define NUMSTR_SOCTET10       6     // number of strains/stresses 	/****/
#define NUMGPT_SOCTET10       5     // total gauss points per element	/****/
#define NUMDIM_SOCTET10       3     // number of dimensions 			/****/
#define NUMCOORD_SOCTET10     4     // number of shape function cooordinates (ksi1-ksi4) 
#define NUMNOD_SOCTET10_FACE  6	    // number of nodes on a CTET10 face (which is a TRI6) 
#define NUMGPT_SOCTET10_FACE  3	    // number of GP    on a CTET10 face (which is a TRI6) 

#define NUMSUBELEMS_SOCTET10 12		// number of sub elements in a CTET10

#define NUMGPT_SUBTET4        1
#define NUMNOD_SUBTET4        4
#define NUMCOORD_SUBTET4      4
#define NUMDIM_SUBTET4        3

// forward declarations
struct _MATERIAL;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;

namespace Elements
{
// forward declarations
class Soctet10Volume;
class Soctet10Surface;
class Soctet10Line;
class Soctet10Register;

/*!
\class So_ctet10

\brief A C++ version of the 10-node tet solid element

\author Alexander Volf (alexander.volf@mytum.de)
*/
class So_ctet10 : public DRT::Element
{
public:

  //! @name Friends
  friend class Soctet10Register;
  friend class Soctet10Surface;
  friend class Soctet10Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_ctet10(int id, int owner);
  
  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_ctet10(const So_ctet10& old);
  
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
  virtual ~So_ctet10();
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

  Evaluate so_ctet10 element stiffness, mass, internal forces, etc.

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

  //! action parameters recognized by so_ctet10
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

  //! kinematic types of so_ctet10
  enum KinematicType
  {
    so_ctet10_geolin,
    so_ctet10_totlag,
    so_ctet10_updlag,
  };
  enum StressType//###################
  {
    so_ctet10_stress_none,
    so_ctet10_stress_gpxyz,
    so_ctet10_stress_gprst,
    so_ctet10_stress_gp123,
    so_ctet10_stress_ndxyz,
    so_ctet10_stress_ndrst,
    so_ctet10_stress_nd123,
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
  vector<RefCountPtr<DRT::Elements::Soctet10Volume> > volume_;
  //! vector of surfaces of this element
    vector<RefCountPtr<DRT::Elements::Soctet10Surface> > surfaces_;
  //! vector of surfaces of this element
  vector<DRT::Elements::Soctet10Surface*>           surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::Elements::Soctet10Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::Elements::Soctet10Line*>              lineptrs_;
  //! vector of history variables for each Gauss point
  vector<RefCountPtr<MAT::Material> >                mat_;


  // internal calculation methods

  // don't want = operator
  So_ctet10& operator = (const So_ctet10& old);
  virtual void so_ctet10_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 struct _MATERIAL*         material);

  void so_ctet10_stress(struct _MATERIAL* material,
                   vector<double>& mydisp,
                   Epetra_SerialDenseMatrix* stresses);
                   
  // SUB_NODE class provides information needed by a node of a sub-element 
  class SUB_NODE
  {	
  	public:
      int local_id;
  	  int global_id;
  	  double my_x[3];
  	  
  	  SUB_NODE();
  	  void init(const int local_id,const int global_id,const Epetra_SerialDenseMatrix& xrefe);
  	  ~SUB_NODE();   				
  };

  // TET4_SUB class provides the information needed to integrate over a sub-element
  class TET4_SUB
  {
    public: 
  	  // RefCountPtr<DRT::Elements::So_ctet10::SUB_NODE> my_nodes;
  	  SUB_NODE my_nodes[4];
  	  double global_baricentre_xi[4];
  	  
  	  TET4_SUB();
  	  void init(const int&, const int&, const int&, const int&,const Epetra_SerialDenseMatrix& xrefe,
  	    const double& xi_babicentre1, const double& xi_babicentre2,
  	    const double& xi_babicentre3, const double& xi_babicentre4);
  	  			
	  //const Epetra_SerialDenseMatrix& integrate(Epetra_SerialDenseVector& xi_c);
	  void integrate(Epetra_SerialDenseVector& xi_c,  Epetra_SerialDenseMatrix& in_Naj);
	  double my_volume();
  	  ~TET4_SUB(); 				
  };

  /* SUB_STRUCTURE holds the sub elements of a CTET10 element 
   * and provides a constructor, to build these from the element reference coordinates
   * NOTE: these already should have the extra 11th node
   */
  class SUB_STRUCTURE
  {
  	public:
  	TET4_SUB my_elements[12];

  	SUB_STRUCTURE(const Epetra_SerialDenseMatrix& xrefe);
  	~SUB_STRUCTURE();
  };
  

/*!
\class DRT::Elements::So_ctet10::Tet_integrator_14point	

\brief A spezialisation of the DRT::Elements::Discrete_integrator class for the quadratic
10-node tetrahedra using a  14 point quadrature, it is exaxt for 0-4th order polynomials.
\author a.vlf
*/

  class Tet10c_integrator_5point : public DRT::Elements::So_integrator
  {
  	public:
  		Tet10c_integrator_5point(SUB_STRUCTURE& );
  		//~Tet10c_integrator_5point();
  };
  

  class L_AJ_integrator : public DRT::Elements::So_integrator
  {
  	public:
  		L_AJ_integrator(SUB_STRUCTURE&);
  		//~L_AJ_integrator();
  };
  
  
  double ctet10_midpoint(int coord);

  // #### hier deklariert, in implementierung soll einfach soh8_mat_sel gerufen werden
  void so_ctet10_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    int gp);       

}; // class So_ctet10

/*!
\class Soctet10Register

\brief A register for tet10 element

\author gee (gee@lnm.mw.tum.de)
*/

class Soctet10Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Soctet10Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soctet10Register

  */
  Soctet10Register(const DRT::Elements::Soctet10Register& old);

  /*!
  \brief Deep copy this instance of Soctet10Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Soctet10Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Soctet10Register;}

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
  virtual ~Soctet10Register();

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


}; // class Soctet10Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Soctet10Surface

\brief An element representing a surface of a So_ctet10 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author maf (frenzel@lnm.mw.tum.de)
*/

class Soctet10Surface : public DRT::Element
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
  Soctet10Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::So_ctet10* parent, const int lsurface);


  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Soctet10Surface(const Soctet10Surface& old);
  
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
  virtual int UniqueParObjectId() const { return ParObject_Soctet10Surface;}

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
  virtual ~Soctet10Surface();

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
  Soctet10Surface& operator = (const Soctet10Surface& old);

  //! The parent element of this surface
  DRT::Elements::So_ctet10* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;

  class SUB_NODE
  {	
  	public:
      int local_id;
  	  int global_id;
  	  double my_x[3];
  	  
  	  SUB_NODE();
  	  void init(const int local_id,const int global_id,const Epetra_SerialDenseMatrix& xrefe);
  	  ~SUB_NODE();   				
  };

  class TRI3_SUB
  {
    public: 
  	  // RefCountPtr<DRT::Elements::So_ctet10::SUB_NODE> my_nodes;
  	  SUB_NODE my_nodes[3];
  	  
  	  TRI3_SUB();
  	  void init(const int&, const int&, const int&,const Epetra_SerialDenseMatrix& xrefe);
  	  			
	  double my_area();
  	  ~TRI3_SUB(); 				
  };

  class SUB_STRUCTURE_SURF
  {
  	public:
  	TRI3_SUB my_surface[4];

  	SUB_STRUCTURE_SURF(const Epetra_SerialDenseMatrix& xrefe);
  	~SUB_STRUCTURE_SURF();
  };

  void sotet4_surface_shapefunc(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseVector** weights);
  

}; // class Soctet10Surface

/*!
\class Soctet10Line

\brief An element representing a line edge of a So_ctet10 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent So_ctet10 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/

class Soctet10Line : public DRT::Element
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
  Soctet10Line(int id, int owner, int nnode, const int* nodeids,
               DRT::Node** nodes, DRT::Elements::So_ctet10* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Soctet10Line(const Soctet10Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Soctet10Line;}
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
  virtual ~Soctet10Line();

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
  Soctet10Line& operator = (const Soctet10Line& old);
  //! The parent element of this line
  DRT::Elements::So_ctet10* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

}; // class Soctet10Line




} // namespace Elements
} // namespace DRT

#endif  // #ifndef SO_CTET10_H_
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET
