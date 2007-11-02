/*!----------------------------------------------------------------------
\file condif2.H

\class DRT::Elements::Condif2

\brief A C++ wrapper for the condif2 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifndef CONDIF2_H
#define CONDIF2_H


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
struct _FLUID_DATA;
/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
// forward declarations
class Discretization;
class Condif2Register;

namespace Elements
{
// forward declarations
class Condif2Line;

/*!
\class Condif2

\brief A C++ wrapper for the condif2 element

\author vg (vgravem@lnm.mw.tum.de)
*/
class Condif2 : public DRT::Element
{
public:

  //! @name Friends
  friend class Condif2Register;
  friend class Condif2Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Condif2(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Element

  */
  Condif2(const Condif2& old);

  /*!
  \brief Deep copy this instance of Condif2 and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return ParObject_Condif2;}

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
  virtual ~Condif2();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

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

  Evaluate condif2 element stiffness, mass, internal forces etc

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

  this method evaluates a line Neumann condition on the condif2 element

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

  //! action parameters recognized by condif2
  enum ActionType
  {
    none,
    calc_condif_systemmat_and_residual
  };

  bool is_higher_order_element(const DRT::Element::DiscretizationType  distype) const;

  //! container containing all kind of condif2 thisandthat
  DRT::Container                                  data_;

  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surface_;

  //! vector of lines adjacent to this element (condif2: 3 or 4 lines)
  vector<RefCountPtr<DRT::Elements::Condif2Line> > lines_;
  //! vector of lines adjacent to this element (condif2: 3 or 4 lines)
  vector<DRT::Elements::Condif2Line*>              lineptrs_;


  // internal calculation methods

  // don't want = operator
  Condif2& operator = (const Condif2& old);

  // calculating system matrix
  void condif2_sys_mat(vector<int>&              lm,
                   vector<double>&           ehist,
                   Epetra_SerialDenseMatrix* sys_mat,
                   Epetra_SerialDenseMatrix* sys_mat_dc,
                   Epetra_SerialDenseVector* residual,
                   struct _MATERIAL*         material,
                   double                    time,
                   double                    timefac,
                   int                       vel_field,
                   int                       discap,
                   bool                      is_stationary
		   );

  // calculate Jacobian matrix and it's determinant
  void condif2_jaco(const Epetra_SerialDenseMatrix& xyze,
	       const Epetra_SerialDenseMatrix& deriv,
               Epetra_SerialDenseMatrix& xjm,
	       double* 		det,
               const int 	iel
	       );

  /*!
   * \brief get the body force from the nodes condition
   * \return ???
   */
  Epetra_SerialDenseVector condif2_getbodyforce();

  // calculate global derivatives w.r.t. x,y,z at point r,s,t are
  void condif2_gder(Epetra_SerialDenseMatrix& derxy,
		 const Epetra_SerialDenseMatrix& deriv,
                 Epetra_SerialDenseMatrix& xjm,
		 double& det,
                 const int iel
		 );

  // calculate second global derivatives w.r.t x/y/z at point r,s,t
  void condif2_gder2(const Epetra_SerialDenseMatrix& xyze,
				    const Epetra_SerialDenseMatrix& xjm,
                                    const Epetra_SerialDenseMatrix& derxy,
				    Epetra_SerialDenseMatrix& derxy2,
				    const Epetra_SerialDenseMatrix& deriv2,
                                    const int iel
				    );

  // evaluate convection-diffusion coefficient matrix
  void condif2_calmat(Epetra_SerialDenseMatrix& estif,
                    Epetra_SerialDenseMatrix& edc,
                    Epetra_SerialDenseVector& eforce,
                    vector<double>&           velint,
                    const double&             hist,
                    Epetra_SerialDenseVector& funct,
                    Epetra_SerialDenseMatrix& derxy,
                    Epetra_SerialDenseMatrix& derxy2,
                    const double&             edeadng,
                    const double&             tau,
                    const double&             kart,
                    const double&             fac,
                    const double&             diffus,
                    const int&                iel,
                    const int&                discap,
                    double                    timefac);
		
  // evaluate convection-diffusion coefficient matrix for stationary problems
  void condif2_calmat_stat(Epetra_SerialDenseMatrix& estif,
                         Epetra_SerialDenseMatrix& edc,
                         Epetra_SerialDenseVector& eforce,
                         vector<double>&           velint,
                         const double&             hist,
                         Epetra_SerialDenseVector& funct,
                         Epetra_SerialDenseMatrix& derxy,
                         Epetra_SerialDenseMatrix& derxy2,
                         const double&             edeadng,
                         const double&             tau,
                         const double&             kart,
                         const double&             fac,
                         const double&             diffus,
                         const int&                iel,
                         const int&                discap);

  // set number of gauss points to element shape default
  DRT::Utils::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype);

  void CreateLinesTri(const int& nline,
                      const int& nnode);

  void CreateLinesQuad(const int& nline,
                       const int& nnode);

}; // class Condif2


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Condif2Register

\brief A register for condif2 element

\author vg (vgravem@lnm.mw.tum.de)
*/
class Condif2Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Condif2Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Condif2Register

  */
  Condif2Register(const DRT::Elements::Condif2Register& old);

  /*!
  \brief Deep copy this instance of Condif2Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Elements::Condif2Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Condif2Register;}

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
  virtual ~Condif2Register();

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

  The method performs all steps needed to init the condif2 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Condif2Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Condif2Line

\brief An element representing a line of a condif2 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent condif2 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author vg (gravem@lnm.mw.tum.de)
*/
class Condif2Line : public DRT::Element
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
  \param parent: The parent condif element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Condif2Line(int id, int owner, int nnode, const int* nodeids,
              DRT::Node** nodes, DRT::Elements::Condif2* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Element

  */
  Condif2Line(const Condif2Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Condif2Line;}

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
  virtual ~Condif2Line();

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
  {return 1;}

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

  /*!
  \brief Evaluate an element

  Evaluate condif2 line element stiffness, mass, internal forces etc

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

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the condif2 element

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

  //! action parameters recognized by Condif2Line
  enum ActionType
  {
     none,
     integrate_Shapefunction
  };

  // don't want = operator
  Condif2Line& operator = (const Condif2Line& old);

  //! The parent element of this line
  DRT::Elements::Condif2* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

  // compute infintesimal line element dr for integration along the line
  double  condif2_substitution(const Epetra_SerialDenseMatrix  xyze,
		               const Epetra_SerialDenseMatrix  deriv,
		               const int iel);

  // Get Rule for Gaussintegration according to DRT::UTIL
  DRT::Utils::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  // integrate shape functions over a line
  void IntegrateShapeFunction(ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm,
                  Epetra_SerialDenseVector&  elevec1);


}; // class Condif2Line




} // namespace Elements
} // namespace DRT




#endif  // #ifndef FLUID2_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
