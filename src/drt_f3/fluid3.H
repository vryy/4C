/*!
\file fluid3.H

\class DRT::Elements::Fluid3

\brief A C++ wrapper for the fluid3 element

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_FLUID3_H
#define DRT_FLUID3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"


// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;
class Fluid3Register;

namespace Elements
{
// forward declarations
class Fluid3Surface;
class Fluid3Line;

/*!
\class Fluid3

\brief A C++ wrapper for the fluid3 element

\author g.bau (bauer@lnm.mw.tum.de)
*/
class Fluid3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Fluid3Register;
  friend class Fluid3Surface;
  friend class Fluid3Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  Fluid3(
    int id,    ///< A unique global id
    int owner  ///< ???
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3(const Fluid3& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid3;}

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
  virtual ~Fluid3();

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

  bool IsAle() const { return is_ale_; }

  //@}


private:

  //! action parameters recognized by fluid3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_beltrami_error,
    calc_turbulence_statistics,
    calc_Shapefunction,
    calc_ShapeDeriv1,
    calc_ShapeDeriv2
  };

  //! Gaussrule
  DRT::Utils::GaussRule3D                         gaussrule_;
  //! flag for euler/ale net algorithm
  bool                                            is_ale_;
  //! container containing all kind of fluid3 thisandthat
  DRT::Container                                  data_;

  //! vector of volumes of this element (length 1)
  vector<DRT::Element*>                           volume_;

  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<RefCountPtr<DRT::Elements::Fluid3Surface> > surfaces_;
  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<DRT::Elements::Fluid3Surface*>              surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::Elements::Fluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::Elements::Fluid3Line*>              lineptrs_;


  // internal calculation methods

  // don't want = operator
  Fluid3& operator = (const Fluid3& old);

  // calculating system matrix
  void f3_sys_mat(const vector<int>&        lm,
                  const vector<double>&     evelnp,
                  const vector<double>&     eprenp,
                  const vector<double>&     ehist,
                  const vector<double>&     edispnp,
                  const vector<double>&     egridv,
                  Epetra_SerialDenseMatrix* sys_mat,
                  Epetra_SerialDenseVector* residual,
                  struct _MATERIAL*         material,
                  ParameterList&            params
                  );

  // calculating system matrix for genalpha
  void f3_genalpha_sys_mat(vector<int>&              lm,
                           vector<double>&           evelnp,
                           vector<double>&           eprenp,
                           vector<double>&           evelaf,
                           vector<double>&           eaccam,
                           Epetra_SerialDenseMatrix* sys_mat,
                           Epetra_SerialDenseVector* elevec,
                           struct _MATERIAL*         material,
                           ParameterList&            params
                           );

  // calculating system rhs for genalpha
  void f3_genalpha_rhs(vector<int>&              lm,
                       vector<double>&           myvelnp,
                       vector<double>&           myprenp,
                       vector<double>&           myvelaf,
                       vector<double>&           myaccam,
                       Epetra_SerialDenseVector& elevec,
                       struct _MATERIAL*         material,
                       ParameterList& 	         params);

  // calculate the stabilisation parameter
  void f3_calc_stabpar(vector<double>           &tau,
                       int                       iel,
                       Epetra_SerialDenseMatrix &xyze,
                       vector<double>           &myvelnp,
                       double                    visc,
                       ParameterList            &params,
                       int                       version);

  /*!
  \brief calculate stabilization parameter
  \return vector with 3 stabilization parameters
   */
  vector<double> f3_caltau(
    const Epetra_SerialDenseMatrix&         xyze,
    const vector<double>&                   evelnp,
    const DRT::Element::DiscretizationType  distype,
    const double                            visc,
    const int                               numnode,
    const double                            timefac,
    const bool                              is_stationary
    );

  /// create array with element corner positions
  Epetra_SerialDenseMatrix f3_getPositionArray(
          const vector<double>&   edispnp);

  /*!
    \brief calculate Jacobian matrix
    \return 3x3 Jacobi matrix
   */
  Epetra_SerialDenseMatrix getJacobiMatrix(
    const Epetra_SerialDenseMatrix& xyze,   ///< positions of element nodes as array
    const Epetra_SerialDenseMatrix& deriv,  ///< shape function derivatives with respect to element coordinates
    const int                       iel     ///< number of nodes in element
    ) const;

  /*!
   \brief calculate determinant of a given Jacobimatrix

   \return determinante of a 3x3 matrix
   */
  double getDeterminante(
    const Epetra_SerialDenseMatrix&  xjm  ///< a 3x3 matrix
    ) const;

  /*!
   * \brief get the body force from the nodes condition
   * \return ???
   */
  Epetra_SerialDenseMatrix f3_getbodyforce (
    const double          time,   ///< ???
    const ParameterList&  params  ///< ???
    );

  /*!
   * \brief calculate global derivatives w.r.t. x,y,z at point r,s,t are
   */
  void f3_gder(
    Epetra_SerialDenseMatrix&        derxy,  ///< first derivative with respect to global coordinates x,y,z
    const Epetra_SerialDenseMatrix&  deriv,  ///< first derivative with respect to local coordinates r,s,t
    const Epetra_SerialDenseMatrix&  xjm,    ///< 3x3 Jacobi matrix dx/dxi
    const double&                    det,    ///< determinante of Jacobi matrix
    const int                        iel     ///< number of nodes for this element
    );

  /*!
   * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
   */
  void f3_gder2(
    const Epetra_SerialDenseMatrix&  xyze,    ///< nodal position array
    const Epetra_SerialDenseMatrix&  xjm,     ///< 3x3 Jacobi matrix dx/dxi
    const Epetra_SerialDenseMatrix&  derxy,   ///< first derivative with respect to global coordinates x,y,z
    Epetra_SerialDenseMatrix&        derxy2,  ///< second derivative with respect to global coordinates x,y,z
    const Epetra_SerialDenseMatrix&  deriv2,  ///< secand derivative with respect to element coordinates r,s,t
    const int                        iel      ///< number of nodes for this element
    );

  /*!
   * \brief evaluate error for beltrami test problem
   */
  void f3_int_beltrami_err(
    vector<double>&       myvelnp,
    vector<double>&       myprenp,
    struct _MATERIAL*     material,
    ParameterList&    	  params
    );

  /*!
    \brief calculate spatial mean values for channel flow
    (requires cartesian mesh)
   
                                                             gammi 07/07

    <pre>
                                                             
    this method assumes that each 2 dimensional integration element
    in the homogeneous plane has the same size!!!

    The necessary element integration is performed in here. The element
    is cut into two (HEX8) or three (quadratic elements) planes, the
    spatial functions (velocity, pressure etc.) are integrated over this
    plane and this element contribution is added to a processor local
    vector (see formulas below for a exact description of the output).
    The method assumes, that all elements are of the same rectangular
    shape in the "inplanedirection".
 
 
                       ^ normdirect       integration plane       
                       |                /         
                       |               /          
                       |                          
                 +-----|-------------+            
                /|     |            /|            
               / |     |           / |            
              /  |     |          /  |            
             /   |     |         /   |            
            /    +-----|--------/----+ ---- additional integration
           /    /|     |       /    /|      plane for quadratic elements
          /    / |     |      /    / |            
         +-------------------+    /  |            
         |   /   |     *-----|---+------------>   
         |  /    +----/------|--/----+         inplanedirect[1]
         | /    /    /       | /    /            
         |/    /    /        |/    /   \         
         +---------+---------+    /     \        
         |   /    /          |   /       integration plane
         |  /    /           |  /                 
         | /    /            | /                  
         |/    /             |/                   
         +----/--------------+                    
             /                                    
            /   inplanedirect[0]                                
               
               
    Example for a mean value evaluation:
               
                            1.0       /                   
                     _               |                    
                     u = -------- *  | u(x,y,z) dx dy dz =
                         +---        |                    
                          \         / A                   
                          / area                          
                         +---                             
                                                                  
                                                               
                           1.0      /                             
                                   |            area           
                     =  -------- * | u(r,s,t) * ---- dr ds dt  
                        +---       |              4            
                         \        /  [-1:1]^2                  
                         / area                                                                          
                        +---                                                                             
                                
                   
                                                               
                            1.0      /                         
                                    |            1             
                     =   -------- * | u(r,s,t) * - dr ds dt    
                                    |            4             
                          numele   /  [-1:1]^2                 
                   
                                   |                        |
                                   +------------------------+
                                this is the integral we compute!
 
    The factor 1/4 is necessary since we use a reference element of
    size 2x2
 
    The method computes:
                       _             _             _             _
              numele * u  , numele * v  , numele * w  , numele * p
                       ___           ___           ___           ___
                        ^2            ^2            ^2            ^2
    and       numele * u  , numele * v  , numele * w  , numele * p
 
    as well as numele.
    All results are communicated vi the parameter list!
    
    </pre>
 
   */

  void f3_calc_means(
    vector<double>&           vel  ,
    vector<double>&           velsq,
    ParameterList& 	      params
    );
    
  /*!
   * \brief evaluate fluid coefficient matrix (instationary case)
   */
  void f3_calmat(
    Epetra_SerialDenseMatrix&        estif,
    Epetra_SerialDenseVector&        eforce,
    const vector<double>&            velint,
    const vector<double>&            histvec,
    const vector<double>&            gridvint,
    const double&                    press,
    const Epetra_SerialDenseMatrix&  vderxy,
    const Epetra_SerialDenseMatrix&  vderxy2,
    const vector<double>&            gradp,
    const Epetra_SerialDenseVector&  funct,
    const vector<double>&            tau,
    const Epetra_SerialDenseMatrix&  derxy,
    const Epetra_SerialDenseMatrix&  derxy2,
    const vector<double>&            edeadng,
    const double&                    fac,
    const double&                    visc,
    const int&                       iel,
    ParameterList&                   params
    );

  /*!
   * \brief evaluate fluid coefficient matrix (stationary case)
   */
  void f3_calmat_stationary(
    Epetra_SerialDenseMatrix&        estif,
    Epetra_SerialDenseVector&        eforce,
    const vector<double>&            velint,
    const vector<double>&            histvec,
    const vector<double>&            gridvint,
    const double&                    press,
    const Epetra_SerialDenseMatrix&  vderxy,
    const Epetra_SerialDenseMatrix&  vderxy2,
    const vector<double>&            gradp,
    const Epetra_SerialDenseVector&  funct,
    const vector<double>&            tau,
    const Epetra_SerialDenseMatrix&  derxy,
    const Epetra_SerialDenseMatrix&  derxy2,
    const vector<double>&            edeadng,
    const double&                    fac,
    const double&                    visc,
    const int&                       iel,
    ParameterList&                   params
    );

  void CreateLinesTet(const int& nline,
                      const int& nnode);

  void CreateLinesHex(const int& nline,
                      const int& nnode);

  void CreateSurfacesTet(const int& nline,
                         const int& nnode);

  void CreateSurfacesHex(const int& nline,
                         const int& nnode);

  /*!
   \brief get optimal gaussrule for given discretization type
   \return enum with optimal gausrule for this element
   */
  DRT::Utils::GaussRule3D getOptimalGaussrule(
    const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
    );

  // evaluate fluid coefficient matrix and rhs for a generalised
  // alpha problem
  void f3_genalpha_calmat(
    Epetra_SerialDenseMatrix&  elemat,
    Epetra_SerialDenseVector&  elevec,
    vector<double>&            accintam,
    vector<double>&            velintaf,
    Epetra_SerialDenseMatrix&  vderxyaf,
    Epetra_SerialDenseMatrix&  vderxyaf2,
    vector<double>&            velintnp,
    Epetra_SerialDenseMatrix&  vderxynp,
    double&                    prenp,
    vector<double>&            pderxynp,
    vector<double>&            edeadaf,
    Epetra_SerialDenseVector&  funct,
    Epetra_SerialDenseMatrix&  derxy,
    Epetra_SerialDenseMatrix&  derxy2,
    vector<double>&            tau,
    const double&              fac,
    const double&              visc,
    const int&                 iel,
    ParameterList&             params);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;

  ActionType convertStringToActionType(const string& action) const;

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
class Fluid3Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Fluid3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Fluid3Register

  */
  Fluid3Register(const DRT::Elements::Fluid3Register& old);

  /*!
  \brief Deep copy this instance of Fluid3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Fluid3Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid3Register;}

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
  virtual ~Fluid3Register();

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
\class Fluid3Surface

\brief An element representing a surface of a fluid3 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
class Fluid3Surface : public DRT::Element
{
public:

  friend class Fluid3Line;

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
  Fluid3Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::Fluid3* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3Surface(const Fluid3Surface& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid3Surface;}

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
  virtual ~Fluid3Surface();

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

#ifdef XFEM
  // to be removed soon and without warning!!! a.gerstenberger
  DRT::Elements::Fluid3* GetParent(){return parent_;}
#endif

private:

  //! action parameters recognized by Fluid3Surface
  enum ActionType
  {
    none,
    calc_Shapefunction,
    calc_ShapeDeriv1,
    calc_ShapeDeriv2
  };

  // don't want = operator
  Fluid3Surface& operator = (const Fluid3Surface& old);

  //! The parent element of this surface
  DRT::Elements::Fluid3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::Elements::Fluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::Elements::Fluid3Line*>              lineptrs_;

  //  compute kovariant metric tensor G for fluid surface element
  //                                                  gammi 04/07

  void  f3_metric_tensor_for_surface(
    const Epetra_SerialDenseMatrix  xyze,
    const Epetra_SerialDenseMatrix  deriv,
    Epetra_SerialDenseMatrix&       metrictensor,
    double                         *drs);

  void CreateLinesTri(const int& nline,
                      const int& nnode);

  void CreateLinesQuad(const int& nline,
                       const int& nnode);

}; // class Fluid3Surface

class Fluid3Line : public DRT::Element
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

  Fluid3Line(   int id,
                int owner,
                int nnode,
                const int* nodeids,
                DRT::Node** nodes,
                DRT::Elements::Fluid3Surface* surfaceParent,
                DRT::Elements::Fluid3* parent,
                const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3Line(const Fluid3Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid3Line;}

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
  virtual ~Fluid3Line();

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

  //! action parameters recognized by Fluid3Line
  enum ActionType
  {
     none,
     calc_Shapefunction,
     calc_ShapeDeriv1,
     calc_ShapeDeriv2
  };

  // don't want = operator
  Fluid3Line& operator = (const Fluid3Line& old);

  //! The surface parent element of this line
  DRT::Elements::Fluid3Surface* surfaceParent_;

  //! The volume parent element of this line
  DRT::Elements::Fluid3* parent_;                  // problem !!! wth initialization




  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;


}; // class Fluid3Line


} // namespace Elements
} // namespace DRT




#endif  // #ifndef FLUID3_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
