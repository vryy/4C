/*!----------------------------------------------------------------------
\file so_hex8.H

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
#include "so_surface.H"
#include "so_line.H"

// prestressing object
#if defined(PRESTRESS) || defined(POSTSTRESS)
#include "prestress.H"
#endif

#include "../drt_mat/material.H"

// Several parameters which are fixed for Solid Hex8
#define NUMNOD_SOH8  8      ///< number of nodes
#define NODDOF_SOH8  3      ///< number of dofs per node
#define NUMDOF_SOH8  24     ///< total dofs per element
#define NUMSTR_SOH8  6      ///< number of strains/stresses
#define NUMGPT_SOH8  8      ///< total gauss points per element
#define NUMDIM_SOH8  3      ///< number of dimensions

// forward declarations
struct _MATERIAL;


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Soh8Register;

/*!
\brief A C++ version of the 8-node hex solid element

A structural 8-node hexahedral solid element for large deformations. As its
discretization is fixed many data structures are evaluated just once and kept
for performance (e.g. shape functions, derivatives, etc.,
see DRT::ELEMENTS::So_hex8::Integrator_So_hex8). It heavily uses
Epetra objects and methods and therefore relies on their performance.

There are 2 sets of EAS enhancements for GL-strains to alleviate locking
(see DRT::ELEMENTS::So_hex8::EASType).

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_hex8 : public DRT::Element
{
public:

  //! @name Friends
  friend class Soh8Register;
  friend class Sosh8Register;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
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


  //const vector<double> GetFibervec(){return fiberdirection_;};
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
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    calc_homog_dens,
    eas_init_multi,
    eas_set_multi,
    postprocess_stress,
    prestress_update,
    multi_readrestart
  };

  //! kinematic types of so_hex8, but yet only total-lagrangian is implemented
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

  /*!
   * \brief EAS technology enhancement types of so_hex8
   *
   * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
   * \param soh8_easmild: 9 parameters consisting of modes to alleviate
   * shear locking (bending) and main incompressibility modes
   * \param soh8_easfull: 21 parameters to prevent almost all locking modes.
   * Equivalent to all 30 parameters to fully complete element
   * with quadratic modes (see Andelfinger 1993 for details) and therefore
   * also suitable for distorted elements.
   * \param soh8_eassoh8: related to Solid-Shell, 7 parameters to alleviate
   * inplane (membrane) locking and main modes for Poisson-locking
   * \param soh8_easnone: no EAS i.e. displacement based with tremendous locking
   */
  enum EASType
  {
    soh8_easfull,
    soh8_easmild,
    soh8_eassosh8,
    soh8_easnone
  };
  //! kinematic type
  KinematicType                                    kintype_;
  //! type of EAS technology
  EASType                                          eastype_;
  //! number of EAS parameters (alphas), defined by 'EASType'
  int                                              neas_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                   data_;
  //! vector of inverses of the jacobian in material frame
  vector<Epetra_SerialDenseMatrix>                 invJ_;
  //! determinant of Jacobian in material frame
  vector<double>                                   detJ_;

#if defined(PRESTRESS) || defined(POSTSTRESS)
  /// Prestressing object
  RCP<DRT::ELEMENTS::PreStress>                    prestress_;

  // compute Jacobian mapping wrt to deformed configuration
  void UpdateJacobianMapping(const vector<double>& disp,
                             DRT::ELEMENTS::PreStress& prestress);

  // compute defgrd in all gp for given disp
  void DefGradient(const vector<double>& disp,
                   Epetra_SerialDenseMatrix& gpdefgrd,
                   DRT::ELEMENTS::PreStress& prestress);
#endif



  //! vector for fiber direction, needs better place in future
  //vector<double>                                   fiberdirection_;


  // internal calculation methods

  //! don't want = operator
  So_hex8& operator = (const So_hex8& old);


  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping();

  //! Calculate nonlinear stiffness and mass matrix
  virtual void soh8_nlnstiffmass(
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy=false,   ///< stress output option
      const bool                euler_almansi=false);  ///< strain output option

  //! Lump mass matrix (bborn 07/08)
  void soh8_lumpmass(Epetra_SerialDenseMatrix* emass);

  //! Extrapolate Gauss-point values (e.g. stresses) to nodes
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

  /*!
   * \brief Interface to material laws working with soh8.
   *
   * Here the interface to any material takes place. Stress and C-mat must be retrieved
   * whereas all necessary data is passed. Add whatever your material needs, but
   * make sure that exchange is not overdone performance-wise.
   * Right now EAS is based on GL-Strains which is ok for any hyperelastic material.
   * This means that e.g. the deformation gradient is not alleviated from locking
   * and therefore only to be used with disp-based soh8.
   *
   * \param *stress (out): Voigt-Vector of stresses at current gp
   * \param *cmat (out): Elasticity matrix at current gp
   * \param *density (out): density of material
   * \param *glstrain (in): Voigt-Vector of GL-strains at current gp
   * \param *defgrd (in): F at current gp, CAUTION! only for disp-based soh8
   * \param gp (in): current gp
   * \param params (in): parameterlist to access time, etc. in materials
   * \param action (in): action, multiscale-related
   * \param ele_ID (in): current element Id
   * \param time (in): current algorithmic time
   * \param dt (in): current algorithmic timestepsize
   * \param action (in): action for material law (default empty)
   * */
  void soh8_mat_sel(Epetra_SerialDenseVector* stress,
                    Epetra_SerialDenseMatrix* cmat,
                    double* density,
                    const Epetra_SerialDenseVector* glstrain,
                    Epetra_SerialDenseMatrix* defgrd,
                    const int gp,
                    ParameterList&  params,
                    const string action = string());

  //! EAS technology, init
  void soh8_easinit();

  //! Re-initialize EAS data, needed for sosh8 morphing
  void soh8_reiniteas(const DRT::ELEMENTS::So_hex8::EASType EASType);

  //! EAS technology, setup necessary data
  void soh8_eassetup(
          vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
          double& detJ0,                      // det of Jacobian at origin
          Epetra_SerialDenseMatrix& T0invT,   // maps M(origin) local to global
          const Epetra_SerialDenseMatrix& xrefe);     // material element coords


  // multi-scale section
  void soh8_initvol(ParameterList& params);

  void soh8_homog(ParameterList& params);

  void soh8_set_eas_multi(ParameterList&  params);

  void soh8_eas_init_multi(ParameterList&  params);

  void soh8_read_restart_multi(ParameterList& params);


  /// temporary method for compatibility with solidshell, needs clarification
  vector<double> getthicknessvector() const {dserror("not implemented"); return vector<double>(3);};


}; // class So_hex8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
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

  The method performs all steps needed to init the elements,
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



} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SOH8_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
