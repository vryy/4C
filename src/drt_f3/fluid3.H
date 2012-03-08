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

#include "Teuchos_RefCountPtr.hpp"


#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations

class Fluid3ImplParameter;
template<DRT::Element::DiscretizationType distype> class Fluid3Impl;
template<DRT::Element::DiscretizationType distype> class Fluid3GenalphaResVMM;
class Fluid3Boundary;

class Fluid3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Fluid3Type"; }

  static Fluid3Type & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  /// pre-evaluation
  virtual void PreEvaluate(DRT::Discretization& dis,
                          Teuchos::ParameterList& p,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                          Teuchos::RCP<Epetra_Vector>          systemvector1,
                          Teuchos::RCP<Epetra_Vector>          systemvector2,
                          Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  static Fluid3Type instance_;
};

/*!
\brief A C++ wrapper for the fluid3 element

\author g.bau (bauer@lnm.mw.tum.de)
*/
class Fluid3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Fluid3Boundary;

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
  virtual DiscretizationType Shape() const
  {
    return distype_;
  };

  /*!
  \brief set discretization type of element
  */
  virtual void SetDisType(DRT::Element::DiscretizationType shape)
  {
    distype_ = shape;

    if (distype_ == DRT::Element::tet10)
    {
      cout << endl;
      cout << "Warning:The integration rule of a tet10 element was changed from tet-5-point to tet-11-point integration rule, "
      " since the deg. of precision 3 is to small as a standard"<< endl << endl;
    }
    return;
  };

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {
    return DRT::UTILS::getNumberOfElementLines(distype_);
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    return DRT::UTILS::getNumberOfElementSurfaces(distype_);
  }

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const
  {
    return DRT::UTILS::getNumberOfElementVolumes(distype_);
  }

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
  virtual int UniqueParObjectId() const { return Fluid3Type::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

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
  {
    // number of Dof's is fluid specific.
    // Therefore, it is not pushed into drt_utilis_local_connectivity
    switch(distype_)
    {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet4:
    case DRT::Element::tet10:
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
    case DRT::Element::pyramid5:
    case DRT::Element::nurbs8:
    case DRT::Element::nurbs27:
      return 4;
      break;
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs4:
    case DRT::Element::nurbs9:
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      return 3;
      break;
    case DRT::Element::line2:
    case DRT::Element::line3:
      dserror("1D Fluid elements are not supported");
    default:
      dserror("discretization type %s not yet implemented", (DRT::DistypeToString(distype_)).c_str());
    }
    return 0;
  }

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.

  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const;

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

  virtual DRT::ElementType & ElementType() const
  { return Fluid3Type::Instance(); }

  //@}

  //! @name Input and Creation

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

  //! @name Time-dependent subgrid scales
  /*!
  \brief Memory allocation for subgrid-scale arrays
  */
  void ActivateTDS(int nquad,int nsd, double** saccn=NULL, double** sveln=NULL, double** svelnp=NULL);

  /*!
  \brief Nonlinear update for current subgrid velocities according to the current
         residual --- reduced version for afgenalpha and one-step-theta
  */
  void UpdateSvelnpInOneDirection(
    const double  fac1   ,
    const double  fac2   ,
    const double  fac3   ,
    const double  resM   ,
    const double  alphaF ,
    const int     dim    ,
    const int     iquad  ,
    double&       svelaf );

  /*!
  \brief Nonlinear update for current subgrid velocities according to the current
         residual --- svelnp as additional return value
  */
  void UpdateSvelnpInOneDirection(
    const double  fac1   ,
    const double  fac2   ,
    const double  fac3   ,
    const double  resM   ,
    const double  alphaF ,
    const int     dim    ,
    const int     iquad  ,
    double&       svelnp,
    double&       svelaf);

  //@}


  //! @name Other

  /*!
  \brief Flag for ALE form of equations
  */
  bool IsAle() const { return is_ale_; }

  virtual void SetIsAle(bool is_ale)
  {
	  is_ale_ = is_ale;

    return;
  };

  /*!
  \brief Returns the quantity Cs*h^2 from the dynamic Smagorinsky model
  */
  double CsDeltaSq() const {return Cs_delta_sq_; }

  /*!
  \brief Returns the subgrid velocity at time n (sveln_)
  */
  Epetra_SerialDenseMatrix Sveln() const {return sveln_; }

  /*!
  \brief Returns the subgrid velocity at time n+1 (svelnp_)
  */
  Epetra_SerialDenseMatrix Svelnp() const {return svelnp_; }

private:

  //! action parameters recognized by fluid3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_porousflow_fluid_coupling,
    calc_loma_mono_odblock,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_genalpha_average_for_subscales_and_residual,
    calc_fluid_error,
    calc_turbulence_statistics,
    calc_loma_statistics,
    calc_turbscatra_statistics,
    calc_dissipation,
    calc_model_params_mfsubgr_scales,
    calc_fluid_box_filter,
    calc_smagorinsky_const,
    get_gas_constant,
    integrate_shape,
    calc_node_normal,
    calc_fluid_elementvolume,
    set_general_fluid_parameter,
    set_time_parameter,
    set_turbulence_parameter,
    set_loma_parameter,
    calc_adjoint_systemmat_and_residual,
    set_general_adjoint_parameter,
    set_adjoint_time_parameter
  };

  //! discretization type
  DRT::Element::DiscretizationType distype_;

  //! flag for euler/ale net algorithm
  bool                                            is_ale_;

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
    \brief same as above for turbulent scalar transport
   */
  template<int iel>
  void f3_calc_scatra_means(
    DRT::Discretization&      discretization,
    vector<double>&           velocitypressure,
    vector<double>&           scalar,
    ParameterList&            params
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
    bool                      dyn_smagorinsky,
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
  double&                  xcenter,
  double&                  ycenter,
  double&                  zcenter);

/*!
   \brief compute parameters of multifractal subgrid-scale model
*/
template<int NEN, int NSD, DRT::Element::DiscretizationType DISTYPE>
void f3_get_mf_params(
  ParameterList&       params,
  RCP<MAT::Material>   mat   ,
  vector<double>&      vel   ,
  vector<double>&      fsvel );

  /*!
    \brief Calculate node normals using the volume integral in Wall (7.13).
    That formula considers all surfaces of the element, not only the
    free surfaces.

    Acc. to Wall (7.13)

    \param elevec1 (out)      : grad(N) integrated over Element
    \param edispnp (in)       : Displacement-vector
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  void ElementNodeNormal(ParameterList& 		 params,
                         DRT::Discretization&       discretization,
                         vector<int>&               lm,
                         Epetra_SerialDenseVector&  elevec1);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */

  ActionType convertStringToActionType(const string& action) const;

}; // class Fluid3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a boundary element of a fluid3 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid3 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

*/

class Fluid3BoundaryType : public DRT::ElementType
{
public:

  std::string Name() const { return "Fluid3BoundaryType"; }

  static Fluid3BoundaryType & Instance() { return instance_; };

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Fluid3BoundaryType instance_;
};


// class Fluid3Boundary

class Fluid3Boundary : public DRT::Element
{
public:

  //! @name Constructors and destructors and related methods

  //! number of space dimensions
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
  Fluid3Boundary(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Fluid3* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3Boundary(const Fluid3Boundary& old);

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
  {
    return DRT::UTILS::getNumberOfElementLines(Shape());
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    return DRT::UTILS::getNumberOfElementSurfaces(Shape());
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

  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return Fluid3BoundaryType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Fluid3Boundary();

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
  {
    return parent_->NumDofPerNode(node);
  }

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

  virtual DRT::ElementType & ElementType() const
  { return Fluid3BoundaryType::Instance(); }

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

  virtual DRT::ELEMENTS::Fluid3* ParentElement()
  {
    return parent_;
  }

  virtual int SurfaceNumber()
  {
    return lsurface_;
  }

  //@}


  enum ActionType
  {
    none,
    integrate_Shapefunction,
    areacalc,
    flowratecalc,
    calc_flowrate,
    flowratederiv,
    Outletimpedance,
    calc_node_normal,
    calc_node_curvature,
    calc_surface_tension,
    enforce_weak_dbc,
    mixed_hybrid_dbc,
    conservative_outflow_bc,
    calc_Neumann_inflow,
    integ_pressure_calc,
    center_of_mass_calc,
    traction_velocity_component,
    traction_Uv_integral_component,
    edge_based_stabilization,
    ghost_penalty
  };

private:

  // don't want = operator
  Fluid3Boundary& operator = (const Fluid3Boundary& old);

  //! The parent element of this surface
  DRT::ELEMENTS::Fluid3* parent_;

  //! The local surface number of this surface w.r.t to the parent_ element
  int                    lsurface_;


}; // class Fluid3Boundary

} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef FLUID3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
