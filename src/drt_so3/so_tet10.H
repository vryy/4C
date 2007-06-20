/*-----------------------------------------**
*/
#ifdef D_SOTET10
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SO_TET10_H_
#define SO_TET10_H_

#endif /*SO_TET10_H_*/


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

// Several parameters which are fixed for Solid Hex8
////#define NUMNOD_SOH8 8      // number of nodes
////#define NODDOF_SOH8 3      // number of dofs per node
////#define NUMDOF_SOH8 24     // total dofs per element
////#define NUMSTR_SOH8 6      // number of strains/stresses
////#define NUMGPT_SOH8 8      // total gauss points per element
////#define NUMDIM_SOH8 3      // number of dimensions

#define NUMNOD_SOTET10 10      // number of nodes
#define NODDOF_SOTET10 3      // number of dofs per node
#define NUMDOF_SOTET10 30     // total dofs per element
/**/ #define NUMSTR_SOTET10 6      // number of strains/stresses
/**/ #define NUMGPT_SOTET10 8      // total gauss points per element
/**/ #define NUMDIM_SOTET10 3      // number of dimensions


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
class Sotet10Volume;
class Sotet10Surface;
class Sotet10Line;
class Sotet10Register;

/*!
\class So_tet10

\brief A C++ version of the 10-node tet solid element

\author volf (alexander.volf@mytum.de)
*/

//class So_hex8 : public DRT::Element
class So_tet10 : public DRT::Element
{
public:

  //! @name Friends
  //friend class Soh8Register;
  friend class Tet10Register;
  //friend class Soh8Surface;
  friend class Tet10Surface;
  //friend class Soh83Line;
  friend class Tet103Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  //So_hex8(int id, int owner);
  So_tet10(int id, int owner);
  
  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  //So_hex8(const So_hex8& old);
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
  //virtual ~So_hex8();
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

  Evaluate so_hex8 element stiffness, mass, internal forces, etc.

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

  // //! action parameters recognized by so_hex8
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

  // //! kinematic types of so_hex8
  //! kinematic types of so_tet10
  enum KinematicType
  {
    //soh8_geolin,
    sot10_geolin,
    //soh8_totlag,
    sot10_totlag,
    //soh8_updlag,
    sot10_updlag,
  };
  enum StressType
  {
    //soh8_stress_none,
    sot8_stress_none,
    //soh8_stress_gpxyz,
    sot8_stress_gpxyz,
    //soh8_stress_gprst,
    sot8_stress_gprst,
    //soh8_stress_gp123,
    sot8_stress_gp123,
    //soh8_stress_ndxyz,
    sot8_stress_ndxyz,
    //soh8_stress_ndrst,
    sot8_stress_ndrst,
    //soh8_stress_nd123,
    sot8_stress_nd123,
  };

  //! kinematic type
  KinematicType                                   kintype_;
  //! type of stress evaluation
  StressType                                      stresstype_;
  //! no. gaussian points in each spatial dimension
  int                                             ngp_[3];
  //! number of the material law
  int                                             material_;
  //! container containing so_hex3 eas stuff//ommited
  //DRT::Container                                  easdata_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                  data_;
  //! vector of volume (length 1) of this element
  //vector<RefCountPtr<DRT::Elements::Soh8Volume> > volume_;
  vector<RefCountPtr<DRT::Elements::Sotet10Volume> > volume_;
  //! vector of surfaces of this element
  //vector<RefCountPtr<DRT::Elements::Soh8Surface> > surfaces_;
  vector<RefCountPtr<DRT::Elements::Sotet10Surface> > surfaces_;
  //! vector of surfaces of this element
  //vector<DRT::Elements::Soh8Surface*>           surfaceptrs_;
  vector<DRT::Elements::Sotet10Surface*>           surfaceptrs_;
  //! vector of lines adjacent to this element
  //vector<RefCountPtr<DRT::Elements::Soh8Line> > lines_;
  vector<RefCountPtr<DRT::Elements::Sotet10Line> > lines_;
  //! vector of lines adjacent to this element
  //vector<DRT::Elements::Soh8Line*>              lineptrs_;
  vector<DRT::Elements::Sotet10Line*>              lineptrs_;
  //! vector of history variables for each Gauss point
  vector<RefCountPtr<Material> >                mat_;


  // internal calculation methods

  // don't want = operator
  //So_hex8& operator = (const So_hex8& old);
  So_tet10& operator = (const So_tet10& old);
  //virtual void soh8_nlnstiffmass(vector<int>&              lm,
  virtual void sot10_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 struct _MATERIAL*         material);

  //void soh8_stress(struct _MATERIAL* material,
  void sot10_stress(struct _MATERIAL* material,
                   vector<double>& mydisp,
                   Epetra_SerialDenseMatrix* stresses);

  //void soh8_shapederiv(Epetra_SerialDenseMatrix** shapefct,
  void sot10_shapederiv(Epetra_SerialDenseMatrix** shapefct,
                       Epetra_SerialDenseMatrix** deriv,
                       Epetra_SerialDenseVector** weights);

  //void soh8_mat_sel(Epetra_SerialDenseVector* stress,
  // #### hier deklariert, in implementierung soll einfach soh8_mat_sel gerufen werden
  void sot10_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    struct _MATERIAL* material,
                    int gp);

}; // class So_tet10

/*!
\class Sotet10Register

\brief A register for shell8 element

\author gee (gee@lnm.mw.tum.de)
*/

//class Soh8Register : public DRT::ElementRegister
class Sotet10Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  //Soh8Register(DRT::Element::ElementType etype);
  Sotet10Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  //Soh8Register(const DRT::Elements::Soh8Register& old);
  Sotet10Register(const DRT::Elements::Sotet10Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  //DRT::Elements::Soh8Register* Clone() const;
  DRT::Elements::Sotet10Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
 // virtual int UniqueParObjectId() const { return ParObject_Soh8Register;}
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
  //virtual ~Soh8Register();
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


}; // class Soh8Register


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

//class Soh8Surface : public DRT::Element
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
//  Soh8Surface(int id, int owner, int nnode, const int* nodeids,
  //           DRT::Node** nodes, DRT::Elements::So_hex8* parent, const int lsurface);
  Sotet10Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::Elements::So_tet10* parent, const int lsurface);


  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
 // Soh8Surface(const Soh8Surface& old);
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
  //virtual int UniqueParObjectId() const { return ParObject_Soh8Surface;}
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
  //virtual ~Soh8Surface();
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
  //Soh8Surface& operator = (const Soh8Surface& old);
  Sotet10Surface& operator = (const Sotet10Surface& old);

  //! The parent element of this surface
  //DRT::Elements::So_hex8* parent_;
  DRT::Elements::So_tet10* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;

  //void soh8_surface_integ(
  void sotet10_surface_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      const Epetra_SerialDenseMatrix* xsrefe,// (i) material element coords
      const double r,                        // (i) coord in r-direction
      const double s);                       // (i) coord in s-direction
}; // class Sotet10Surface

/*!
\class Sotet10Line

\brief An element representing a line edge of a So_tet10 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent So_hex8 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
//class Soh8Line : public DRT::Element
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
  //Soh8Line(int id, int owner, int nnode, const int* nodeids,
    //         DRT::Node** nodes, DRT::Elements::So_hex8* parent, const int lline);
  Sotet10Line(int id, int owner, int nnode, const int* nodeids,
               DRT::Node** nodes, DRT::Elements::So_tet10* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  //Soh8Line(const So8Line& old);
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
  // virtual int UniqueParObjectId() const { return ParObject_Soh8Line;}
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
  //virtual ~Soh8Line();
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
  //Soh8Line& operator = (const Soh8Line& old);

  Sotet10Line& operator = (const Sotet10Line& old);
  //! The parent element of this line
  //DRT::Elements::So_hex8* parent_;
  DRT::Elements::So_tet10* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

}; // class Soh8Line




} // namespace Elements
} // namespace DRT

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET10
