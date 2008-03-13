/*!----------------------------------------------------------------------
\file so_hex8.H

\class DRT::ELEMENTS::So_hex8

\brief A C++ version of the 8-node hex solid element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SOH8_H
#define SOH8_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
using namespace std;
using namespace Teuchos;

#include "so_integrator.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

#include "../drt_mat/material.H"

// Several parameters which are fixed for Solid Hex8
const int NUMNOD_SOH8 = 8;      // number of nodes
const int NODDOF_SOH8 = 3;      // number of dofs per node
const int NUMDOF_SOH8 = 24;     // total dofs per element
const int NUMSTR_SOH8 = 6;      // number of strains/stresses
const int NUMGPT_SOH8 = 8;      // total gauss points per element
const int NUMDIM_SOH8 = 3;      // number of dimensions

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
class Soh8Volume;
class Soh8Surface;
class Soh8Line;
class Soh8Register;

/*!
\class So_hex8

\brief A C++ version of the 8-node hex solid element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_hex8 : public DRT::Element
{
public:

  //! @name Friends
  friend class Soh8Register;
  friend class Soh8Surface;
  friend class Soh83Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  So_hex8(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_hex8(const So_hex8& old);

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
  virtual int NumSurface() const {return 6;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 12;}

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
  virtual int UniqueParObjectId() const { return ParObject_So_hex8;}

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
  virtual ~So_hex8();

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
  /*!
  \brief Query names of element data to be visualized using BINIO

  The element fills the provided map with key names of
  visualization data the element wants to visualize AT THE CENTER
  of the element geometry. The values is supposed to be dimension of the
  data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
  or 9 (nonsym. tensor)

  Example:
  \code
    // Name of data is 'Owner', dimension is 1 (scalar value)
    names.insert(pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(pair<string,int>("StressesXYZ",6));
  \endcode

  \param names (out): On return, the derived class has filled names with
                      key names of data it wants to visualize and with int dimensions
                      of that data.
  */
  virtual void VisNames(map<string,int>& names);

  /*!
  \brief Query data to be visualized using BINIO of a given name

  The method is supposed to call this base method to visualize the owner of
  the element.
  If the derived method recognizes a supported data name, it shall fill it
  with corresponding data.
  If it does NOT recognizes the name, it shall do nothing.

  \warning The method must not change size of data

  \param name (in):   Name of data that is currently processed for visualization
  \param data (out):  data to be filled by element if element recognizes the name
  */
  virtual void VisData(const string& name, vector<double>& data);

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */

  virtual bool ReadElement();

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


  const vector<double> GetThickvec(){return thickvec_;};
  const vector<double> GetFibervec(){return fiberdirection_;};
  const vector<double> soh8_ElementCenterRefeCoords();

  //@}

protected:

  //! action parameters recognized by so_hex8
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
    calc_struct_update_istep,
    calc_homog_stressdens,
    postprocess_stress
  };

  //! kinematic types of so_hex8
  enum KinematicType
  {
    soh8_geolin,
    soh8_totlag,
    soh8_updlag
  };
  //! stress evaluation types of so_hex8
  enum StressType
  {
    soh8_stress_none,
    soh8_stress_gpxyz,
    soh8_stress_gprst,
    soh8_stress_gp123,
    soh8_stress_ndxyz,
    soh8_stress_ndrst,
    soh8_stress_nd123,
    soh8_stress_cxyz     // interpolated element center stresses
  };
  //! EAS technology types of so_hex8
  enum EASType
  {
    soh8_easfull,
    soh8_easmild,
    soh8_eassosh8,
    soh8_easnone
  };
  //! kinematic type
  KinematicType                                   kintype_;
  //! type of EAS technology
  EASType                                         eastype_;
  //! number of EAS parameters (alphas)
  int                                             neas_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                  data_;
  //! vector of volumes of this element (length 1)
  vector<DRT::Element*>                           volume_;
  //! vector of surfaces of this element
  vector<RefCountPtr<DRT::ELEMENTS::Soh8Surface> > surfaces_;
  //! vector of surfaces of this element
  vector<DRT::ELEMENTS::Soh8Surface*>           surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::ELEMENTS::Soh8Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::ELEMENTS::Soh8Line*>              lineptrs_;
  //! vector of element stress tensor in Voigt notation

  //! some distypes need rewinding of NodeIds due to bad input from mesher
  //! flag whether process of rewinding is done
  bool                                            donerewinding_;

  //! vector in thickness direction for compatibility with sosh8
  vector<double>                                   thickvec_;
  vector<double>                                   fiberdirection_;

  // internal calculation methods

  // don't want = operator
  So_hex8& operator = (const So_hex8& old);

  virtual void soh8_nlnstiffmass(vector<int>&              lm,
                                 vector<double>&           disp,
                                 vector<double>&           residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 Epetra_SerialDenseMatrix* stresses,
                                 Epetra_SerialDenseMatrix* elestrain,
                                 const double              time);

  void soh8_expol(Epetra_SerialDenseMatrix& stresses,
                  Epetra_SerialDenseMatrix& nodalstresses);

  /*!
  \brief Create a integrator for the SoHex8 element

  this constructor provides Object Oriented So_integrator functionality by adapting the strictly procedural
  DRT::UTILS::IntegrationPoints3D struct.
  */
  class Integrator_So_hex8 : public DRT::ELEMENTS::So_integrator
  {
  	public:
  		Integrator_So_hex8();
  };
  /***************OBSOLETE**********************************/
  void soh8_shapederiv(Epetra_SerialDenseMatrix** shapefct,
                       Epetra_SerialDenseMatrix** deriv,
                       Epetra_SerialDenseVector** weights);
  /***************OBSOLETE**********************************/
  // the deformation gradient, Gauss point number, element ID and the
  // action parameter (default: empty string) is needed only for multi-scale
  // materials
  void soh8_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    const Epetra_SerialDenseMatrix* defgrd,
                    const int gp,
                    const int ele_ID,
                    const double time,
                    const string action = string());

  // EAS technology section
  void soh8_easinit();

  void soh8_eassetup(
          Epetra_SerialDenseMatrix** M_GP,    // M-matrix evaluated at GPs
          double& detJ0,                      // det of Jacobian at origin
          Epetra_SerialDenseMatrix& T0invT,   // maps M(origin) local to global
          const Epetra_SerialDenseMatrix& xrefe);     // material element coords

  void soh8_initvol(ParameterList& params);

  void soh8_homog(ParameterList&  params,
                  vector<double>& disp,
                  const double    time,
                  vector<double>& residual);

  vector<double> getthicknessvector() const {dserror("not implemented"); return vector<double>(3);};


}; // class So_hex8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class Soh8Register

\brief A register for shell8 element

\author gee (gee@lnm.mw.tum.de)
*/
class Soh8Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Soh8Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Soh8Register

  */
  Soh8Register(const DRT::ELEMENTS::Soh8Register& old);

  /*!
  \brief Deep copy this instance of Soh8Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Soh8Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Soh8Register;}

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
  virtual ~Soh8Register();

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

  The method performs all steps needed to init the wedge elements,
  e.g. defines shell-thickness direction
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
\class Soh8Surface

\brief An element representing a surface of a So_hex8 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author maf (frenzel@lnm.mw.tum.de)
*/

class Soh8Surface : public DRT::Element
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
  Soh8Surface(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::So_hex8* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Soh8Surface(const Soh8Surface& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Soh8Surface;}

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
  virtual ~Soh8Surface();

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

  /*!
    \brief Evaluate a boundary condition

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): A reference to the underlying discretization
    \param condition (in)     : The condition to be evaluated
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,

    \return 0 if successful, negative otherwise
    */
  virtual int Evaluate(ParameterList& params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elematrix1,
                       Epetra_SerialDenseMatrix& elematrix2,
                       Epetra_SerialDenseVector& elevector1,
                       Epetra_SerialDenseVector& elevector2,
                       Epetra_SerialDenseVector& elevector3);
  //@}


protected:
  //! action parameters recognized by so_hex8_surface
  enum ActionType
  {
    none,
    calc_struct_constrvol,
    calc_struct_volconstrstiff,
    calc_struct_constrarea,
    calc_init_vol,
    calc_surfstress_stiff,
    calc_struct_areaconstrstiff,
    calc_struct_monitarea
  };


private:

  // don't want = operator
  Soh8Surface& operator = (const Soh8Surface& old);

  //! The parent element of this surface
  DRT::ELEMENTS::So_hex8* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;


  // Submethod to compute the enclosed volume for vlume constraint boundary condition
  double ComputeConstrVols(Epetra_SerialDenseMatrix xscurr);

  void ComputeVolconstrStiff(Epetra_SerialDenseMatrix xscurr,
		  Epetra_SerialDenseMatrix& elematrix);

  void ComputeVolconstrVolDeriv(Epetra_SerialDenseMatrix xscurr,
		  Epetra_SerialDenseVector& elevector1);


  void soh8_surface_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      vector<double>* unrm,                  // (o) unit normal
      const Epetra_SerialDenseMatrix* xs,    // (i) material element coords
      const double r,                        // (i) coord in r-direction
      const double s);                       // (i) coord in s-direction
}; // class Soh8Surface

/*!
\class Soh8Line

\brief An element representing a line edge of a So_hex8 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent So_hex8 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gee (gee@lnm.mw.tum.de)
*/
class Soh8Line : public DRT::Element
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
  Soh8Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::So_hex8* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Soh8Line(const Soh8Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Soh8Line;}

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
  virtual ~Soh8Line();

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
  Soh8Line& operator = (const Soh8Line& old);

  //! The parent element of this line
  DRT::ELEMENTS::So_hex8* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;

  void soh8_line_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      vector<double>* unrm,                  // (o) unit normal
      const Epetra_SerialDenseMatrix* xs,    // (i) material element coords
      const double r);                       // (i) coord in r-direction

}; // class Soh8Line


} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SOH8_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
