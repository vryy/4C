/*!
\file xfluid3.H

\class DRT::ELEMENTS::XFluid3

\brief A C++ wrapper for the fluid3 element

This file contains the element specific service routines like
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and gausspoint loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifndef DRT_XFLUID3_H
#define DRT_XFLUID3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "fluid3.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/dof_management.H"
#include "xfluid3_interpolation.H"


// forward declarations
struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class XFluid3Surface;
class XFluid3Line;
class XFluid3Impl;
class XFluid3Register;

/*!
\class XFluid3

\brief A C++ wrapper for the fluid3 element

\author g.bau (bauer@lnm.mw.tum.de)
*/
class XFluid3 : public DRT::Element
{
public:

  //! @name Friends
  friend class XFluid3Register;
  friend class XFluid3Surface;
  friend class XFluid3Line;
  friend class XFluid3Impl;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  XFluid3(
    int id,    ///< A unique global id
    int owner  ///< ???
    );

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
   else if (NumNode()==6 || NumNode()==15 || NumNode()==5) return 5;
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
  virtual int NumDofPerNode(const DRT::Node& node) const {return this->eleDofManager_.NumDofPerNode(node.Id(), 0);}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const {return this->eleDofManager_.NumDofPerElement();}

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
  \brief Evaluate an element, i.e. call the implementation to evaluate fluid3
  element stiffness, mass, internal forces or evaluate errors, statistics or updates
  etc. directly.

  Following implementations of the element are allowed:
  <pre>
    o Evaluation of fluid system matrix and residual for the One-Step-Theta or BDF2
      method

    o Evaluation of fluid system matrix and residual for the stationary fluid solver

    o Evaluation of fluid system matrix and residual for the Generalised alpha method
      ---
      the element matrix computation could be skipped, this is useful for a
      speedup if itemax is reached in the nonlinear iteration and the matrix is not
      used in further steps
  </pre>

  Up to now, forther evaluation of element integrals is done to
  <pre>
    o evaluate the L2 error for beltrami flows

    o compute integrals over homogeneous planes for turbulent channel flow simulations
  </pre>

  

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
  
  DRT::UTILS::GaussRule3D GetGaussrule() const { return gaussrule_; }

  template <DRT::Element::DiscretizationType  DISTYPE>
  bool isHigherOrderElement() const
  {
    bool hoel = true;
    switch (DISTYPE)
    {
    case hex8: case hex20: case hex27: case tet10: case wedge15:
      hoel = true;
      break;
    case tet4: case wedge6: case pyramid5:
      hoel = false;
      break;
    default:
      dserror("DISTYPE unknown!");
    }
    return hoel;
  };
  
  //@}


private:

  //! action parameters recognized by fluid3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_genalpha_average_for_subscales_and_residual,
    calc_fluid_stationary_systemmat_and_residual,    
    calc_fluid_beltrami_error,
    calc_turbulence_statistics,
    calc_fluid_box_filter,
    calc_smagorinsky_const,
    calc_Shapefunction,
    calc_ShapeDeriv1,
    calc_ShapeDeriv2,
    store_xfem_info
  };

  //! stabilisation actions recognized by fluid3 (genalpha implementation)
  enum StabilisationAction
  {
    stabaction_unspecified,
    subscales_quasistatic,
    subscales_time_dependent,
    inertia_stab_drop,
    inertia_stab_keep,
    pstab_assume_inf_sup_stable,
    pstab_use_pspg,
    convective_stab_none,
    convective_stab_supg,
    viscous_stab_none,
    viscous_stab_gls,
    viscous_stab_gls_only_rhs,
    viscous_stab_agls,
    viscous_stab_agls_only_rhs,
    continuity_stab_none,
    continuity_stab_yes,
    cross_stress_stab_none,
    cross_stress_stab,
    cross_stress_stab_only_rhs,
    reynolds_stress_stab_none,
    reynolds_stress_stab_only_rhs
  };

  //! physical turbulence models recognized by fluid3_genalpha implementation
  enum TurbModelAction
  {
    no_model,
    smagorinsky,
    smagorinsky_with_wall_damping,
    dynamic_smagorinsky
  };
  
  //! map that converts strings to actions (stabilisation)
  static map<string,StabilisationAction>          stabstrtoact_;

  //! Gaussrule
  DRT::UTILS::GaussRule3D                         gaussrule_;
  //! flag for euler/ale net algorithm
  bool                                            is_ale_;
  //! container containing all kind of fluid3 thisandthat
  DRT::Container                                  data_;

  //! vector of volumes of this element (length 1)
  vector<DRT::Element*>                           volume_;

  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<RefCountPtr<DRT::ELEMENTS::XFluid3Surface> > surfaces_;
  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  vector<DRT::ELEMENTS::XFluid3Surface*>              surfaceptrs_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::ELEMENTS::XFluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::ELEMENTS::XFluid3Line*>              lineptrs_;
  
 
  //! some distypes need rewinding of NodeIds due to bad input from mesher
  //! flag to rewind this element, must be availabale for friend class fluid3register
protected:
  bool                                            rewind_;
  //! flag whether process of rewinding is done
  bool                                            donerewinding_;

private:
  //! matrices of subscale acceleration values in gausspoints of this element
  blitz::Array<double,2>                        sub_acc_old_;

  //! matrices of subscale velocity values, current iteration value, in gausspoints of this element
  blitz::Array<double,2>                        sub_vel_;
  //! matrices of subscale velocity values, last timestep, in gausspoints of this element
  blitz::Array<double,2>                        sub_vel_old_;

  //! vector of subscale pressure values, current iteration value, in gausspoints
  blitz::Array<double,1>                        sub_pre_;
  //! vector of subscale pressure values (history values) in gausspoints
  blitz::Array<double,1>                        sub_pre_old_;

  //! store information about the nodal-unknowns of this element
  XFEM::ElementDofManager                         eleDofManager_;

  // internal calculation methods

  // don't want = operator
  XFluid3& operator = (const XFluid3& old);

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
    is cut into two (HEX8) or three (quadratic elements) planes (plus
    additional planes for visualisation purposes, defined by planes
    vector), the spatial functions (velocity, pressure etc.) are
    integrated over this plane and this element contribution is added
    to a processor local vector (see formulas below for a exact
    description of the output).
    The method assumes, that all elements are of the same rectangular
    shape in the "inplanedirection". In addition, it is assumed that
    the sampling planes are distributed equidistant in the element.


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
           /    /|     |       /    /|      plane (for quadratic elements)
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
    and        numele * u  , numele * v  , numele * w  , numele * p

                        _ _           _ _           _ _  
   as well as  numele * u*v, numele * u*w, numele * v*w
  

    as well as numele.
    All results are communicated via the parameter list!

    </pre>

   */

  void f3_calc_means(
    vector<double>&           vel  ,
    ParameterList& 	      params
    );

  /*!
    \brief one point integration of convolutions with box filter
    function. The filter function is not normalized yet, this has
    to be done as a 'postprocessing' step after the patchvolume
    is calculated, i.e. after the volume contribution of all
    elements adjacent to a node has been added to the
    normalisation constant.

    See reference

    A.E. Tejada-Martinez, K.E. Jansen
    Spatial test filters for dynamic model large-eddy simulation with
    finite elements.
    Communications in numerical methods in engineering, 2000

    for details

    \param evelaf                        (in )  : nodal velocities
    \param vel_hat                       (out)  : integrated convolution
                                                  with velocity
    \param reystr_hat                    (out)  : integrated convolution
                                                  with reynolds stresses 
    \param modeled_stress_grid_scale_hat (out)  : integrated convolution
                                                  with subgrid scale
                                                  stress tensor model
    \param volume                        (out)  : element volume
    \return void
    
   */
  void f3_apply_box_filter(
    blitz::Array<double, 2>&  evelaf,
    blitz::Array<double, 1>&  vel_hat,
    blitz::Array<double, 2>&  reystr_hat,
    blitz::Array<double, 2>&  modeled_stress_grid_scale_hat,
    double&                   volume
    );

  /*!
    \brief compute the quantities necessary to determine Cs for the dynamic
    model.

    References are

    
    M. Germano, U. Piomelli, P. Moin, W.H. Cabot:
    A dynamic subgrid-scale eddy viscosity model
    (Phys. Fluids 1991)
      
    or

    D.K. Lilly:
    A proposed modification of the Germano subgrid-scale closure method
    (Phys. Fluids 1992)

    or
    A.E. Tejada-Martinez
    Dynamic subgrid-scale modeling for large eddy simulation of turbulent
    flows with a stabilized finite element method
    (Phd thesis, Rensselaer Polytechnic Institute, Troy, New York)

    <pre>
    
    L_ij : the tensor of resolved components of the stress tensor associated
           with the scales of motion between the resolved and the test scale

    M_ij : the tensor of modeled components of the stress tensor associated
           with the scales of motion between the resolved and the test scale
           (modulo 2Cs, a constant)
    
    </pre>

    
    \param evelaf                        (in )  : nodal velocities
    \param vel_hat                       (in)   : integrated convolution
                                                  with velocity
    \param reystr_hat                    (in)   : integrated convolution
                                                  with reynolds stresses 
    \param modeled_stress_grid_scale_hat (in)   : integrated convolution
                                                  with subgrid scale
                                                  stress tensor model
    \param LijMij                        (out)  : trace of product of resolved
                                                  stress tensor with modeled
                                                  stress tensor
    \param MijMij                        (out)  : trace of product of modeled
                                                  resolved stress tensor with
                                                  itself
    \param center                        (out)  : element center
    \return void
    
   */

void f3_calc_smag_const_LijMij_and_MijMij(
  blitz::Array<double, 2> evelaf_hat,
  blitz::Array<double, 3> ereynoldsstress_hat,
  blitz::Array<double, 3> efiltered_modeled_subgrid_stress_hat,
  double&                 LijMij,
  double&                 MijMij,
  double&                 center);

  
  void CreateLinesTet(const int& nline,
                      const int& nnode);

  void CreateLinesHex(const int& nline,
                      const int& nnode);

  void CreateSurfacesTet(const int& nline,
                         const int& nnode);

  void CreateSurfacesHex(const int& nline,
                         const int& nnode);

  void CreateSurfaceWedge(const int& nsurf,
                          const int& wedgetype);

  void CreateSurfacesPyramid();

  /*!
   \brief get optimal gaussrule for given discretization type
   \return enum with optimal gausrule for this element
   */
  DRT::UTILS::GaussRule3D getOptimalGaussrule(
    const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
    );

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;

  ActionType convertStringToActionType(const string& action) const;

  //! convert string to stabilisation action
  StabilisationAction ConvertStringToStabAction (const string& action) const;

  /// Internal implementation class for specific fluid element
  XFluid3Impl* Impl(
          const int numnode,
          const int numparamvelx,
          const int numparamvely,
          const int numparamvelz,
          const int numparampres);

protected:
  /*!
   * \brief check, whether NodeIds need rewinding due to wrong input from mesher,
   * based on negative Jacobian determinant.
   * \return boolean rewind or not rewind
   */
  bool checkRewinding();

}; // class Fluid3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\class XFluid3Register

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
  XFluid3Register(const DRT::ELEMENTS::XFluid3Register& old);

  /*!
  \brief Deep copy this instance of Fluid3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::XFluid3Register* Clone() const;

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


}; // class XFluid3Register


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
  XFluid3Surface(int id, int owner, const int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::XFluid3* parent, const int lsurface);

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
  {return parent_->NumDofPerNode(node);}

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

  // to be removed soon or moved to the DRT::Element class
  // please ask me before using it! a.gerstenberger
  const DRT::ELEMENTS::XFluid3* GetParent(){return parent_;}

private:

  //! action parameters recognized by Fluid3Surface
  enum ActionType
  {
    none,
    calc_Shapefunction,
    calc_ShapeDeriv1,
    calc_ShapeDeriv2,
    integrate_Shapefunction
  };

  // don't want = operator
  XFluid3Surface& operator = (const XFluid3Surface& old);

  //! The parent element of this surface
  DRT::ELEMENTS::XFluid3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  const int                    lsurface_;
  //! vector of lines adjacent to this element
  vector<RefCountPtr<DRT::ELEMENTS::XFluid3Line> > lines_;
  //! vector of lines adjacent to this element
  vector<DRT::ELEMENTS::XFluid3Line*>              lineptrs_;

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

  /*!
  \brief integrate shapefunction over surface element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
  */

  virtual void IntegrateShapeFunction(ParameterList&         params,
                                  DRT::Discretization&       discretization,
                                  vector<int>&               lm,
                                  Epetra_SerialDenseVector&  elevec1,
                                      const std::vector<double>& edispnp);


}; // class Fluid3Surface

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

  XFluid3Line(   int id,
                int owner,
                int nnode,
                const int* nodeids,
                DRT::Node** nodes,
                DRT::ELEMENTS::XFluid3Surface* surfaceParent,
                DRT::ELEMENTS::XFluid3* parent,
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
  {return parent_->NumDofPerNode(node);}

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
  XFluid3Line& operator = (const XFluid3Line& old);

  //! The surface parent element of this line
  DRT::ELEMENTS::XFluid3Surface* surfaceParent_;

  //! The volume parent element of this line
  DRT::ELEMENTS::XFluid3* parent_;                  // problem !!! wth initialization




  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;


}; // class XFluid3Line


} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef XFLUID3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
