/*!----------------------------------------------------------------------
\file so_tet4.H

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_TET4_H
#define SO_TET4_H


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

#define NUMNOD_SOTET4       4      ///< number of nodes
#define NODDOF_SOTET4       3      ///< number of dofs per node
#define NUMDOF_SOTET4      12      ///< total dofs per element
#define NUMSTR_SOTET4       6      ///< number of strains/stresses 	/****/
#define NUMGPT_SOTET4       1      ///< total gauss points per element	/****/
#define NUMDIM_SOTET4       3      ///< number of dimensions 			/****/
#define NUMCOORD_SOTET4	    4 		///< number of shape function cooordinates (ksi1-ksi4)
#define NUMNOD_SOTET4_FACE  3		///< number of nodes on a TET4 face (which is a TRI3)
#define NUMGPT_SOTET4_FACE  1		///< number of GP    on a TET4 face (which is a TRI3)

// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Sotet4Volume;
class Sotet4Register;

/*!
\brief A C++ version of the 4-node tet solid element

*/

class So_tet4 : public DRT::Element
{
public:

  //! @name Friends
  friend class Sotet4Register;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_tet4(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_tet4(const So_tet4& old);

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
  virtual int UniqueParObjectId() const { return ParObject_So_tet4;}

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
  virtual ~So_tet4();
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

  Evaluate so_tet4 element stiffness, mass, internal forces, etc.

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
    calc_struct_nlnstifflmass,
    calc_struct_stress,
    postprocess_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep, //!< reset elementwise internal variables
                             //!< during iteration to last converged state
    update_gl
  };

  //! kinematic types of so_tet10
  enum KinematicType
  {
    so_tet4_geolin,
    so_tet4_totlag,
    so_tet4_updlag
  };
  enum StressType
  {
    so_tet4_stress_none,
    so_tet4_stress_gpxyz,
    so_tet4_stress_gprst,
    so_tet4_stress_gp123,
    so_tet4_stress_ndxyz,
    so_tet4_stress_ndrst,
    so_tet4_stress_nd123
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
  //! volume of the element
  double                                                   V_;
  //! vector of partial derivatives in material frame
  vector<Epetra_SerialDenseMatrix>                         nxyz_;

  //! green lagrange strains storage needed for prestressing structures
#if defined(PRESTRESS) || defined(POSTSTRESS)
  virtual inline Teuchos::RCP<Epetra_SerialDenseMatrix>& PreStrains()
  { return glprestrain_; }
  RCP<Epetra_SerialDenseMatrix>                    glprestrain_;
#endif

  // internal calculation methods

  // don't want = operator
  So_tet4& operator = (const So_tet4& old);
  
  // init the inverse of the jacobian and its determinant
  // in the material configuration
  virtual void InitJacobianMapping();

  virtual void so_tet4_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 Epetra_SerialDenseMatrix* stresses,
                                 Epetra_SerialDenseMatrix* elestrain,
                                 struct _MATERIAL*         material,
                                 bool                      cauchy = false);

  //! lump mass matrix (bborn 07/08)
  void so_tet4_lumpmass(Epetra_SerialDenseMatrix* emass);  //!< element mass matrix


  void so_tet4_expol(Epetra_SerialDenseMatrix& stresses,
                     Epetra_SerialDenseMatrix& nodalstresses);


  // #### hier deklariert, in implementierung soll einfach soh8_mat_sel gerufen werden
  void so_tet4_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    int gp);

}; // class So_tet4

/*!
\brief A register for tet4 element

\author gee (gee@lnm.mw.tum.de)
*/

class Sotet4Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Sotet4Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  Sotet4Register(const DRT::ELEMENTS::Sotet4Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Sotet4Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Sotet4Register;}

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
  virtual ~Sotet4Register();

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


}; // class Sotet4Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================





} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SO_TET4_H_
#endif  // #ifdef  CCADISCRET
#endif  // #ifdef  D_SOLID3


