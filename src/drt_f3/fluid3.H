/*!
\file fluid3.H

\brief A C++ wrapper for the fluid3 element

This file contains the element specific service routines like
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and gausspoint loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifndef DRT_FLUID3_H
#define DRT_FLUID3_H

#include "../drt_lib/linalg_serialdensematrix.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/linalg_fixedsizematrix.H"


// forward declarations
//struct _MATERIAL;

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Fluid3Surface;

template<DRT::Element::DiscretizationType distype,DRT::Element::DiscretizationType pdistype> class Fluid3SurfaceWeakDBC;

class Fluid3Line;
template<DRT::Element::DiscretizationType distype> class Fluid3Impl;
template<DRT::Element::DiscretizationType distype> class Fluid3StationaryImpl;
template<DRT::Element::DiscretizationType distype> class Fluid3lin_Impl;
template<DRT::Element::DiscretizationType distype> class Fluid3GenalphaResVMM;
class Fluid3Register;

/*!
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

  // All the template specializations that are used, new ones must be
  // added here.
  friend class Fluid3Impl<DRT::Element::tet4>;
  friend class Fluid3Impl<DRT::Element::pyramid5>;
  friend class Fluid3Impl<DRT::Element::wedge6>;
  friend class Fluid3Impl<DRT::Element::hex8>;
  friend class Fluid3Impl<DRT::Element::tet10>;
  friend class Fluid3Impl<DRT::Element::wedge15>;
  friend class Fluid3Impl<DRT::Element::hex20>;
  friend class Fluid3Impl<DRT::Element::hex27>;
  friend class Fluid3StationaryImplInterface;

  // All the template specializations that are used, new ones must be
  // added here.
  friend class Fluid3StationaryImpl<DRT::Element::tet4>;
  friend class Fluid3StationaryImpl<DRT::Element::pyramid5>;
  friend class Fluid3StationaryImpl<DRT::Element::wedge6>;
  friend class Fluid3StationaryImpl<DRT::Element::hex8>;
  friend class Fluid3StationaryImpl<DRT::Element::tet10>;
  friend class Fluid3StationaryImpl<DRT::Element::wedge15>;
  friend class Fluid3StationaryImpl<DRT::Element::hex20>;
  friend class Fluid3StationaryImpl<DRT::Element::hex27>;
  friend class Fluid3lin_ImplInterface;
  // All the template specializations that are used, new ones must be
  // added here.
  friend class Fluid3lin_Impl<DRT::Element::tet4>;
  friend class Fluid3lin_Impl<DRT::Element::pyramid5>;
  friend class Fluid3lin_Impl<DRT::Element::wedge6>;
  friend class Fluid3lin_Impl<DRT::Element::hex8>;
  friend class Fluid3lin_Impl<DRT::Element::tet10>;
  friend class Fluid3lin_Impl<DRT::Element::wedge15>;
  friend class Fluid3lin_Impl<DRT::Element::hex20>;
  friend class Fluid3lin_Impl<DRT::Element::hex27>;
  friend class Fluid3ImplInterface;

  friend class Fluid3GenalphaResVMM<DRT::Element::tet4>;
  friend class Fluid3GenalphaResVMM<DRT::Element::pyramid5>;
  friend class Fluid3GenalphaResVMM<DRT::Element::wedge6>;
  friend class Fluid3GenalphaResVMM<DRT::Element::hex8>;
  friend class Fluid3GenalphaResVMM<DRT::Element::tet10>;
  friend class Fluid3GenalphaResVMM<DRT::Element::wedge15>;
  friend class Fluid3GenalphaResVMM<DRT::Element::hex20>;
  friend class Fluid3GenalphaResVMM<DRT::Element::hex27>;
  friend class Fluid3GenalphaResVMM<DRT::Element::nurbs8>;
  friend class Fluid3GenalphaResVMM<DRT::Element::nurbs27>;
  friend class Fluid3GenalphaResVMMInterface;

  friend class Fluid3SystemEvaluator;

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
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element
  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volume of this element
  */
  virtual vector<RCP<DRT::Element> > Volumes();

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

  //! @name Geometry related methods

  /*!
   * \brief calculate volume of current element
   * @return volume of current element (out double)
   */
  double CalculateVolume(ParameterList& params, DRT::Discretization& discretization,vector<int>& lm) const;

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
  //bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

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
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


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
    calc_fluid_stationary_systemmat_and_residual,
    calc_fluid_afgenalpha_systemmat_and_residual,
    calc_linear_fluid,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_genalpha_average_for_subscales_and_residual,
    calc_fluid_beltrami_error,
    calc_turbulence_statistics,
    calc_loma_statistics,
    calc_fluid_box_filter,
    calc_smagorinsky_const,
    get_density,
    get_gas_constant,
    integrate_shape,
    calc_node_normal,
    calc_fluid_elementvolume
  };

  //! stabilisation actions recognized by fluid3 (genalpha implementation)
  enum StabilisationAction
  {
    stabaction_unspecified,
    subscales_quasistatic,
    subscales_time_dependent,
    inertia_stab_drop,
    inertia_stab_keep,
    inertia_stab_keep_complete,
    pstab_assume_inf_sup_stable,
    pstab_use_pspg,
    convective_stab_none,
    convective_stab_supg,
    viscous_stab_none,
    viscous_stab_gls,
    viscous_stab_gls_only_rhs,
    viscous_stab_usfem,
    viscous_stab_usfem_only_rhs,
    continuity_stab_none,
    continuity_stab_yes,
    cross_stress_stab_none,
    cross_stress_stab,
    cross_stress_stab_only_rhs,
    reynolds_stress_stab_none,
    reynolds_stress_stab,
    reynolds_stress_stab_only_rhs
  };

  //! linearisation actions recognized by fluid3 (genalpha implementation)
  enum LinearisationAction
  {
    no_linearisation,
    Newton          ,
    fixed_point_like,
    minimal
  };

  //! possible options for the stabilisation parameter
  enum TauType
  {
    smoothed_franca_barrenechea_valentin_wall,
    franca_barrenechea_valentin_wall         ,
    franca_barrenechea_valentin_codina       ,
    fbvw_gradient_based_hk                   ,
    codina                                   ,
    bazilevs                                 ,
    fbvw_wo_dt                               ,
    tau_not_defined
  };


  //! physical turbulence models recognized by fluid3_genalpha implementation
  enum TurbModelAction
  {
    no_model,
    smagorinsky,
    smagorinsky_with_wall_damping,
    dynamic_smagorinsky
  };

  //! options for fine-scale subgrid viscosity
  enum FineSubgridVisc
  {
    no_fssgv,
    smagorinsky_all,
    smagorinsky_small
  };

  //! map that converts strings to actions (stabilisation)
  static map<string,StabilisationAction>          stabstrtoact_;

  //! Gaussrule
  DRT::UTILS::GaussRule3D                         gaussrule_;
  //! flag for euler/ale net algorithm
  bool                                            is_ale_;
  //! container containing all kind of fluid3 thisandthat
  DRT::Container                                  data_;

private:
  //! factor including Smagorinsky's constant
  double                                          Cs_delta_sq_;

  //! matrices of subscale acceleration values in gausspoints of this element
  Epetra_SerialDenseMatrix                        saccn_;
  //! matrices of subscale velocity values, current iteration value, in gausspoints of this element
  Epetra_SerialDenseMatrix                        svelnp_;
  //! matrices of subscale velocity values, last timestep, in gausspoints of this element
  Epetra_SerialDenseMatrix                        sveln_;

  // internal calculation methods

  // don't want = operator
  Fluid3& operator = (const Fluid3& old);

  /*!
   * \brief evaluate error for beltrami test problem
   */
  void f3_int_beltrami_err(
    vector<double>&       myvelnp,
    vector<double>&       myprenp,
    Teuchos::RCP<const MAT::Material> material,
    ParameterList&    	  params
    );

  /*!
    \brief calculate spatial mean values for channel flow
    (requires wall parallel layers of elements)

                                                             gammi 07/07

    <pre>

    this method assumes that each 2 dimensional integration element
    in the homogeneous plane is parallel to the wall!!!

    The necessary element integration is done in here. The element
    is cut into two (HEX8) or three (quadratic elements) planes (plus
    additional planes for visualisation purposes, defined by planes
    vector), the spatial functions (velocity, pressure etc.) are
    integrated over this plane and this element contribution is added
    to a processor local vector (see formulas below for a exact
    description of the output).
    It is assumed that the sampling planes are distributed equidistant
    in the element. The result is normalized by the area afterwards


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

           1.0       /                     1.0      /
    _               |                              |            detJ
    u = -------- *  | u(x,y,z) dx dz =  -------- * | u(r,s,t) * ---- dr ds
        +---        |                   +---       |             h
         \         / A                   \        /  [-1:1]^2     y
         / area                          / area
        +---                            +---                    +--+
                                                               Jacobi-
                                                             determinant
                                                              computed
                                                           from 3d mapping
                                                         h  is the (constant)
                                                          y
                                                        height of the element

    The method computes:
                        _             _             _             _
               numele * u  , numele * v  , numele * w  , numele * p
                        ___           ___           ___           ___
                         ^2            ^2            ^2            ^2
    and        numele * u  , numele * v  , numele * w  , numele * p

                        _ _           _ _           _ _
    as well as  numele * u*v, numele * u*w, numele * v*w


    as well as numele and the element area.
    All results are communicated via the parameter list!

    </pre>

   */
  template<int iel>
  void f3_calc_means(
    DRT::Discretization&      discretization,
    vector<double>&           vel   ,
    vector<double>&           disp  ,
    ParameterList& 	      params
    );

  /*!
    \brief same as above for low-Mach-number flow
   */
  template<int iel>
  void f3_calc_loma_means(
    DRT::Discretization&      discretization,
    vector<double>&           velocitypressure ,
    vector<double>&           temperature  ,
    ParameterList&            params,
    const double              eosfac
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
  template<int iel>
  void f3_apply_box_filter(
    vector<double>&           myvel,
    double*                   vel_hat,
    double*                   reystr_hat,
    double*                   modeled_stress_grid_scale_hat,
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
template<int iel>
void f3_calc_smag_const_LijMij_and_MijMij(
  RCP<Epetra_MultiVector>& filtered_vel                       ,
  RCP<Epetra_MultiVector>& col_filtered_reynoldsstress        ,
  RCP<Epetra_MultiVector>& col_filtered_modeled_subgrid_stress,
  double&                  LijMij,
  double&                  MijMij,
  double&                  center);

  /*!
   \brief get optimal gaussrule for given discretization type
   \return enum with optimal gausrule for this element
   */
  DRT::UTILS::GaussRule3D getOptimalGaussrule(
    const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
    );

  /*!
    \brief integrate shape functions over volume element

    \param discretization (i) the discretisation (required in the nurbs case)
    \param lm             (i) location vector (to get ale displacements from global vectors)
    \param w              (o) shape functions integrated over element (a so called 'weight' vector)

    \return
   */
  template<int iel>
  void integrateShapefunction(
    DRT::Discretization&      discretization,
    vector<int>&              lm            ,
    Epetra_SerialDenseVector& w             );

  /*!
    \brief Calculate node normals using the volume integral in Wall (7.13).
    That formula considers all surfaces of the element, not only the
    free surfaces.

    Acc. to Wall (7.13)

    \param elevec1 (out)      : grad(N) integrated over Element
    \param edispnp (in)       : Displacement-vector
  */
  template<int iel>
  void ElementNodeNormal(ParameterList& 		 params,
                         DRT::Discretization&       discretization,
                         vector<int>&               lm,
                         Epetra_SerialDenseVector&  elevec1);

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

}; // class Fluid3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
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
  Fluid3Register(const DRT::ELEMENTS::Fluid3Register& old);

  /*!
  \brief Deep copy this instance of Fluid3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Fluid3Register* Clone() const;

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

  friend class Fluid3SurfaceWeakDBC<DRT::Element::quad4 ,DRT::Element::hex8   >;
  friend class Fluid3SurfaceWeakDBC<DRT::Element::nurbs9,DRT::Element::nurbs27>;
  friend class Fluid3SurfaceWeakDBCInterface;

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
             DRT::Node** nodes, DRT::ELEMENTS::Fluid3* parent, const int lsurface);

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
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

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
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

private:

  //! action parameters recognized by Fluid3Surface
  enum ActionType
  {
    none,
    integrate_Shapefunction,
    areacalc,
    flowratecalc,
    Outletimpedance,
    calc_node_normal,
    calc_node_curvature,
    calc_surface_tension,
    enforce_weak_dbc,
    conservative_outflow_bc,
    calc_Neumann_inflow
  };

  // don't want = operator
  Fluid3Surface& operator = (const Fluid3Surface& old);

  //! The parent element of this surface
  DRT::ELEMENTS::Fluid3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;


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

  virtual void IntegrateShapeFunction(ParameterList& 		 params,
                          	      DRT::Discretization&       discretization,
                          	      vector<int>&               lm,
                          	      Epetra_SerialDenseVector&  elevec1,
                                      const std::vector<double>& edispnp);

  /*!
    \brief Calculate node normals

    Acc. to Wall (7.13), surface integral.

    \param elevec1 (out)      : Weighted element normal for assembly
    \param edispnp (in)       : Displacement-vector
  */
  void ElementNodeNormal(ParameterList& 		 params,
                       DRT::Discretization&       discretization,
                       vector<int>&               lm,
                       Epetra_SerialDenseVector&  elevec1,
                       const std::vector<double>& edispnp);
  /*!
    \brief Calculate mean curvature H. Interpolate the results to achieve better
    results in the surface tension algorithm (c0 field -> c1 field).

    \param elevec1  (out)     : Nodal values of mean curvature
    \param edispnp  (in)      : Displacement-vector
    \param enormals (in)      : Node normals
  */
  void ElementMeanCurvature(ParameterList& params,
                            DRT::Discretization& discretization,
                            vector<int>& lm,
                            Epetra_SerialDenseVector& elevec1,
                            const std::vector<double>& edispnp,
                            std::vector<double>& enormals);


  /*!
    \brief Integrate surface tension

    \param elevec1   (out)    : RHS contribution of surface tension effect
    \param edispnp   (in)     : Displacement-vector
    \param enormals  (in)     : Node normals
    \param curvature (in)     : Interpolated curvature
  */
  void ElementSurfaceTension(ParameterList& params,
                             DRT::Discretization& discretization,
                             vector<int>& lm,
                             Epetra_SerialDenseVector& elevec1,
                             const std::vector<double>& edispnp,
                             std::vector<double>& enormals,
                             std::vector<double>& ecurvature);

  /*!
  brief integrate elemental areas over a surface

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  */
  void AreaCaculation(ParameterList& params);

 /*!
  brief integrate elemental flow rates over a surface

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  */
  void FlowRateParameterCalculation(ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm);

  /*!
  \brief apply impedance boundary condition (outlet pressures)

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
  */
  void ImpedanceIntegration(ParameterList& params,
          DRT::Discretization&       discretization,
          vector<int>&               lm,
          Epetra_SerialDenseVector&  elevec1);

  /*!
  \brief apply outflow boundary condition which is necessary for the conservative
  element formulation (since the convective term was partially integrated)

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element.
  \param elevec1 (out)      : vector to be filled by element.

  */
  void SurfaceConservativeOutflowConsistency(
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1,
    Epetra_SerialDenseVector&  elevec1);

  /*!
  \brief compute potential Neumann inflow

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.                          and elements
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  */
  void NeumannInflow(
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1,
    Epetra_SerialDenseVector&  elevec1);

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
                DRT::Element* parent,
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
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

private:

  //! action parameters recognized by Fluid3Line
  enum ActionType
  {
     none
  };

  // don't want = operator
  Fluid3Line& operator = (const Fluid3Line& old);

  //! The volume parent element of this line
  DRT::Element* parent_;




  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;


}; // class Fluid3Line


} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef FLUID3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
