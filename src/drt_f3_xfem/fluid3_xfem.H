/*!----------------------------------------------------------------------
\file fluid3_xfem.H

\class DRT::Elements::XFluid3

\brief A XFEM enabled copy of the fluid3 element

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_FLUID3_XFEM_H
#define DRT_FLUID3_XFEM_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "fluid3_xfem_integration.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"


// forward declarations
struct _MATERIAL;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class XFluid3Register;

namespace Elements
{
// forward declarations
class XFluid3Surface;
class XFluid3Line;

/*!
\class XFluid3

\brief A C++ wrapper for the xfluid3 element

\author g.bau (bauer@lnm.mw.tum.de)
*/
class XFluid3 : public DRT::Element
{
public:

  //! @name Friends
  friend class XFluid3Register;
  friend class XFluid3Surface;
  friend class XFluid3Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  XFluid3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  XFluid3(const XFluid3& old);

  /*!
  \brief Deep copy this instance of Fluid3 and return pointer to the copy

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
  {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 12;
   else if (NumNode()==4 || NumNode()==10) return 6;
   else {
    dserror("Could not determine number of lines");
    return -1;
    }
  }   

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 6;
   else if (NumNode()==4 || NumNode()==10) return 4;
   else {
    dserror("Could not determine number of surfaces");
    return -1;
    }
  }

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Get vector of ptrs to the lines of this element

  */
  virtual DRT::Element** Lines();
  
  /*!
  \brief Get vector of ptrs to the surfaces of this element

  */
  virtual DRT::Element** Surfaces();

  /*!
  \brief Get vector of ptrs to the volume of this element

  */
  virtual DRT::Element** Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_XFluid3;}

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
  virtual ~XFluid3();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 4;}

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

  Evaluate fluid3 element stiffness, mass, internal forces etc

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

  this method evaluates a surfaces Neumann condition on the shell element

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

  //! @name Other


  //@}


private:

  //! action parameters recognized by XFluid3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_fluid_beltrami_error,
    calc_ShapefunctDeriv1Deriv2
  };
  
  // constants
  static const int                                NSD_ = 3;
  
  // constant mappings between volume and its surfaces and lines
  // they are used as the one and only reference for such numbering
  static const int                                hex27_surfaces_[6][9];
  static const int                                hex27_lines_[12][3];
  
  static const int                                tet10_surfaces_[4][6];
  static const int                                tet10_lines_[6][3];

  //! gaussrule used for weak form integration
  GaussRule3D                                     gaussrule_;
  //! number of the material law
  int                                             material_;
  //! flag for euler/ale net algorithm
  bool                                            is_ale_;

  //! container containing all kind of fluid3 thisandthat
  DRT::Container                                  data_;

  //! vector of volumes of this element (length 1)
  vector<DRT::Element*>                   volume_;

  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<RefCountPtr<DRT::Elements::XFluid3Surface> > surfaces_;
  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<DRT::Elements::XFluid3Surface*>              surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::Elements::XFluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::Elements::XFluid3Line*>              lineptrs_;




  // internal calculation methods

  // don't want = operator
  XFluid3& operator = (const XFluid3& old);

  // calculating system matrix
  void f3_sys_mat(const vector<int>&              lm,
                  const vector<double>&           evelnp,
                  const vector<double>&           eprenp,
                  const vector<double>&           ehist,
                  Epetra_SerialDenseMatrix* sys_mat,
                  Epetra_SerialDenseVector* residual,
                  const struct _MATERIAL*         material,
                  ParameterList&            params
                 );

  // calculate Jacobian matrix and it's determinant
  void f3_jaco(const Epetra_SerialDenseMatrix& xyze,
               const Epetra_SerialDenseMatrix& deriv,
               Epetra_SerialDenseMatrix& xjm,
               double*      det,
               const int    iel
              );

  // get the body force from the nodes condition
  void f3_getbodyforce (Epetra_SerialDenseMatrix& edeadng,
                        const int                 iel     ,
                        ParameterList&        params);

  // calculate global derivatives w.r.t. x,y,z at point r,s,t are
    void f3_gder(Epetra_SerialDenseMatrix& derxy,
                 const Epetra_SerialDenseMatrix& deriv,
                 const Epetra_SerialDenseMatrix& xjm,
                 const double& det,
                 const int iel
         );

  // calculate second global derivatives w.r.t x/y/z at point r,s,t
    void f3_gder2(const Epetra_SerialDenseMatrix& xyze,
                  const Epetra_SerialDenseMatrix& xjm,
                  const Epetra_SerialDenseMatrix& derxy,
                  Epetra_SerialDenseMatrix& derxy2,
                  const Epetra_SerialDenseMatrix& deriv2,
                  const int iel
                  );
  // evaluate error for beltrami test problem
  void f3_int_beltrami_err(const vector<double>&       myvelnp,
                           const vector<double>&       myprenp,
                           const struct _MATERIAL*     material,
                           ParameterList&        params
                           );

  // evaluate fluid coefficient matrix (instationary case)
    void f3_calmat(Epetra_SerialDenseMatrix&           estif,
                   Epetra_SerialDenseVector&           eforce,
                   const vector<double>&               velint,
                   const vector<double>&               histvec,
                   const vector<double>&               gridvint,
                   const double&                       press,
                   const Epetra_SerialDenseMatrix&     vderxy,
                   const Epetra_SerialDenseMatrix&     vderxy2,
                   const vector<double>&               gradp,
                   const Epetra_SerialDenseVector&     funct,
                   const vector<double>&               tau,
                   const Epetra_SerialDenseMatrix&     derxy,
                   const Epetra_SerialDenseMatrix&     derxy2,
                   const vector<double>&               edeadng,
                   const double&                       fac,
                   const double&                       visc,
                   const int&                          iel,
                   ParameterList&                      params
            );

  // evaluate fluid coefficient matrix (stationary case)
    void f3_calmat_stationary( 
                Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&           eforce,
                const vector<double>&               velint,
                const vector<double>&               histvec,
                const vector<double>&               gridvint,
                const double&                       press,
                const Epetra_SerialDenseMatrix&     vderxy,
                const Epetra_SerialDenseMatrix&     vderxy2,
                const vector<double>&               gradp,
                const Epetra_SerialDenseVector&     funct,
                const vector<double>&               tau,
                const Epetra_SerialDenseMatrix&     derxy,
                const Epetra_SerialDenseMatrix&     derxy2,
                const vector<double>&               edeadng,
                const double&                       fac,
                const double&                       visc,
                const int&                          iel,
                ParameterList&                      params
            );
  
  // get optimal gaussrule for given discretization type
  GaussRule3D get_optimal_gaussrule(const DiscretizationType& distype);
  
  inline vector<double> f3_caltau(
    Epetra_SerialDenseVector&    funct,
    Epetra_SerialDenseMatrix&    deriv,
    const Epetra_SerialDenseMatrix&    xyze,
    Epetra_SerialDenseMatrix&    xjm,
    Epetra_SerialDenseMatrix&    derxy,
    const vector<double>&        evelnp,
    const DRT::Element::DiscretizationType           distype,
    const double                 visc,
    const int                    numnode,
    const double                 timefac,
    const bool                   is_stationary
    );

  inline bool is_higher_order_element(
    const DRT::Element::DiscretizationType  distype);
    
  void CreateLinesTet(const int& nline,
                      const int& nnode);
  
  void CreateLinesHex(const int& nline,
                      const int& nnode);
                       
  void CreateSurfacesTet(const int& nline,
                         const int& nnode);
  
  void CreateSurfacesHex(const int& nline,
                         const int& nnode);

}; // class Fluid3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Fluid3Register

\brief A register for fluid3 element

\author gee (gee@lnm.mw.tum.de)
*/
class XFluid3Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  XFluid3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Fluid3Register

  */
  XFluid3Register(const DRT::Elements::XFluid3Register& old);

  /*!
  \brief Deep copy this instance of Fluid3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::XFluid3Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_XFluid3Register;}

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
  virtual ~XFluid3Register();

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

  The method performs all steps needed to init the fluid3 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Fluid3Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class XFluid3Surface

\brief An element representing a surface of a fluid3 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
class XFluid3Surface : public DRT::Element
{
public:

friend class XFluid3Line;

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
  XFluid3Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::XFluid3* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  XFluid3Surface(const XFluid3Surface& old);

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
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4;
   	else if (NumNode()==3 || NumNode()==6) return 3;
   	else 
   	{
    		dserror("Could not determine number of surfaces");
    	return -1;
    	}
  }    
  
  /*!
  \brief Get vector of ptrs to the lines of this element

  */
  virtual DRT::Element** Lines();
  
  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return ParObject_XFluid3Surface;}

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
  virtual ~XFluid3Surface();

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
  {return 4;}

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

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate fluid3 element stiffness, mass, internal forces etc

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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //@}

private:

  //! action parameters recognized by XFluid3Surface
  enum ActionType
  {
  	 none,
    calc_ShapefunctDeriv1Deriv2
  };

  // don't want = operator
  XFluid3Surface& operator = (const XFluid3Surface& old);

  // constants
  static const int                                NSD_ = 3;
  
  static const int                                quad9_lines_[4][3];
  
  static const int                                tri6_lines_[3][3];
  
  //! The parent element of this surface
  DRT::Elements::XFluid3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;
   //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::Elements::XFluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::Elements::XFluid3Line*>              lineptrs_;

  //  compute kovariant metric tensor G for fluid surface element
  //                                                  gammi 04/07

  void  f3_metric_tensor_for_surface(
    const Epetra_SerialDenseMatrix  xyze,
    const Epetra_SerialDenseMatrix  deriv,
    Epetra_SerialDenseMatrix&       metrictensor,
    double                         *drs);

  GaussRule2D get_optimal_gaussrule(const DiscretizationType& distype);

  void CreateLinesTri(const int& nline,
                      const int& nnode);
  
  void CreateLinesQuad(const int& nline,
                       const int& nnode);

}; // class XFluid3Surface




class XFluid3Line : public DRT::Element
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
  
  XFluid3Line(  int id, 
                int owner, 
                int nnode, 
                const int* nodeids,
                DRT::Node** nodes, 
                DRT::Elements::XFluid3Surface* surfaceParent,
                DRT::Elements::XFluid3* parent,  
                const int lline);   

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  XFluid3Line(const XFluid3Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_XFluid3Line;}

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
  virtual ~XFluid3Line();

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
  {return 4;}

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

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate fluid3 element stiffness, mass, internal forces etc

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
  virtual int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);
               
               
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

  //! action parameters recognized by XFluid3Line
  enum ActionType
  {
	 none, 		
    calc_ShapefunctDeriv1Deriv2
  };
  
  // don't want = operator
  XFluid3Line& operator = (const XFluid3Line& old);

  //! The surface parent element of this line
  DRT::Elements::XFluid3Surface* surfaceParent_;

  //! The volume parent element of this line
  DRT::Elements::XFluid3* parent_;					// problem !!! wth initialization
  

  
  
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;
  
  void f3_shapefunction_for_line(	Epetra_SerialDenseVector&	funct ,
  												Epetra_SerialDenseMatrix& 	deriv1,
  												Epetra_SerialDenseMatrix& 	deriv2,
  												int&                 		iel   ,
  												double&              		r     );

}; // class XFluid3Line



} // namespace Elements
} // namespace DRT




#endif  // #ifndef FLUID3_XFEM_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
