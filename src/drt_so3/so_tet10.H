/*!----------------------------------------------------------------------
\file so_tet10.H

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
#ifndef SO_TET10_H
#define SO_TET10_H


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
#include "so_surface.H"
#include "so_line.H"

#include "../drt_lib/linalg_fixedsizematrix.H"

// gee: note that these are also defined in so_integrator.H
#define NUMNOD_SOTET10      10      ///< number of nodes
#define NODDOF_SOTET10       3      ///< number of dofs per node
#define NUMDOF_SOTET10      30      ///< total dofs per element
#define NUMSTR_SOTET10       6      ///< number of strains/stresses 	/****/
#define NUMGPT_SOTET10       4      ///< total gauss points per element	/****/
#define NUMDIM_SOTET10       3      ///< number of dimensions 			/****/
#define NUMCOORD_SOTET10     4      ///< number of shape function cooordinates (ksi1-ksi4)
#define NUMNOD_SOTET10_FACE  6	    ///< number of nodes on a TET10 face (which is a TRI6)
#define NUMGPT_SOTET10_FACE  3      ///< number of GP    on a TET10 face (which is a TRI6)

// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Sotet10Register;

/*!
\brief A C++ version of the 10-node tet solid element

\author Alexander Volf (alexander.volf@mytum.de)
*/

class So_tet10 : public DRT::Element
{
public:

  //! @name Friends
  friend class Sotet10Register;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owning processor
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
  \param discretization : pointer to discretization for de-assembly
  \param lm (in)        : location matrix for de-assembly
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
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    postprocess_stress
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
  KinematicType                                            kintype_;
  //! type of stress evaluation
  StressType                                               stresstype_;
  //! no. gaussian points in each spatial dimension
  int                                                      ngp_[3];
  //! number of the material law
  int                                                      material_;

  DRT::Container                                           data_;
  //! vector of history variables for each Gauss point
  vector<RefCountPtr<MAT::Material> >                      mat_;




  // internal calculation methods

  // don't want = operator
  So_tet10& operator = (const So_tet10& old);
  virtual void so_tet10_nlnstiffmass(vector<int>&          lm,
                                     vector<double>&           disp,
                                     vector<double>&           residual,
                                     LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* stiffmatrix,
                                     LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* massmatrix,
                                     LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET10,1>* force,
                                     LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOTET10,NUMSTR_SOTET10>* stresses,
                                     LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOTET10,NUMSTR_SOTET10>* elestrain,
                                     struct _MATERIAL*         material,
                                     bool cauchy = false);

  //! lump mass matrix (bborn 07/08)
  void so_tet10_lumpmass(LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* emass);

  void so_tet10_expol(LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOTET10,NUMSTR_SOTET10>& stresses,
                      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOTET10,NUMSTR_SOTET10>& nodalstresses);

  // These functions are basically copied from the So_integrator,
  // along with the rather inconsitent matrix sizes.
  //! Evaluate Tet10 Shapefcts at 4 gausspoints to keep them static
  const vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOTET10,1> >& so_tet10_4gp_shapefcts();
  //! Evaluate Tet10 Derivs at 4 gausspoints to keep them static
  const vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOTET10,NUMCOORD_SOTET10> >& so_tet10_4gp_derivs();
  //! Evaluate Tet10 Weights at 4 gausspoints to keep them static
  const vector<double>&so_tet10_4gp_weights();

  //! Evaluate Tet10 Shapefcts at 10 gausspoints to keep them static
  const vector<LINALG::FixedSizeSerialDenseMatrix<NUMNOD_SOTET10,1> >& so_tet10_10gp_shapefcts();
  //! Evaluate Tet10 Derivs at 10 gausspoints to keep them static
  const vector<LINALG::FixedSizeSerialDenseMatrix<NUMDIM_SOTET10,NUMNOD_SOTET10> >& so_tet10_10gp_derivs();
  //! Evaluate Tet10 Weights at 10 gausspoints to keep them static
  const vector<double>& so_tet10_10gp_weights();

//   void so_tet10_stress(struct _MATERIAL* material,
//                    vector<double>& mydisp,
//                    Epetra_SerialDenseMatrix* stresses);

/*!
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
  | class Tet_integrator_4point 	                        	vlf 09/07   |
  |																	|
  | a special class for quadratic 10-node tetrahedra using a 4 point quadrature 	|
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

  /// became obsolete
  void so_tet10_shapederiv(
  	  Epetra_SerialDenseVector** shapefct_ptr,
      Epetra_SerialDenseMatrix** deriv_gp_ptr,
      Epetra_SerialDenseVector** weights);

  void so_tet10_mat_sel(
                      LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
                      LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                      double* density,
                      LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                      LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
                      const int gp);
  
  /// material selection
  void so_tet10_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    int gp);

}; // class So_tet10

/*!
\brief A register for tet10 element

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




} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef SO_TET10_H_
#endif  // #ifdef  CCADISCRET
#endif  // #ifdef  D_SOLID3
