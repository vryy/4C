/*!----------------------------------------------------------------------
\file fluid2.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifndef FLUID2_H
#define FLUID2_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


// forward declarations
//struct _MATERIAL;
struct _FLUID_DATA;


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Fluid2Line;
class Fluid2ImplInterface;
class Fluid2GenalphaResVMMInterface;
template <DRT::Element::DiscretizationType distype> class Fluid2Impl;
class Fluid2StationaryInterface;
template <DRT::Element::DiscretizationType distype> class Fluid2Stationary;
template <DRT::Element::DiscretizationType distype> class Fluid2GenalphaResVMM;
class Fluid2Register;
class Fluid2THInterface;	// for Taylor-Hood element
template <DRT::Element::DiscretizationType distype> class Fluid2TH;

template<DRT::Element::DiscretizationType distype,DRT::Element::DiscretizationType pdistype> class Fluid2LineWeakDBC;

/*!
\brief A C++ wrapper for the fluid2 element

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2 : public DRT::Element
{
public:

  //! @name Friends
  friend class Fluid2Register;
  friend class Fluid2Line;
  friend class Fluid2ImplInterface;
  friend class Fluid2GenalphaResVMMInterface;
  friend class Fluid2Impl<DRT::Element::tri3>;
  friend class Fluid2Impl<DRT::Element::tri6>;
  friend class Fluid2Impl<DRT::Element::quad4>;
  friend class Fluid2Impl<DRT::Element::quad8>;
  friend class Fluid2Impl<DRT::Element::quad9>;
  friend class Fluid2StationaryInterface;
  friend class Fluid2Stationary<DRT::Element::tri3>;
  friend class Fluid2Stationary<DRT::Element::tri6>;
  friend class Fluid2Stationary<DRT::Element::quad4>;
  friend class Fluid2Stationary<DRT::Element::quad8>;
  friend class Fluid2Stationary<DRT::Element::quad9>;
  friend class Fluid2GenalphaResVMM<DRT::Element::tri3>;
  friend class Fluid2GenalphaResVMM<DRT::Element::tri6>;
  friend class Fluid2GenalphaResVMM<DRT::Element::quad4>;
  friend class Fluid2GenalphaResVMM<DRT::Element::quad8>;
  friend class Fluid2GenalphaResVMM<DRT::Element::quad9>;
  friend class Fluid2GenalphaResVMM<DRT::Element::nurbs4>;
  friend class Fluid2GenalphaResVMM<DRT::Element::nurbs9>;
  friend class Fluid2THInterface;	// Taylor-Hood element (only quad)
  friend class Fluid2TH<DRT::Element::quad4>;
  friend class Fluid2TH<DRT::Element::quad9>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Fluid2(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid2(const Fluid2& old);

  /*!
  \brief Deep copy this instance of Fluid2 and return pointer to the copy

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
  \brief Get vector of RefCountPtrs to the volume of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid2;}

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
  virtual ~Fluid2();

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
  virtual int NumDofPerNode(const DRT::Node& node) const;

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
  \brief Evaluate an element

  Evaluate fluid2 element stiffness, mass, internal forces etc

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

  this method evaluates a line Neumann condition on the fluid2 element

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

  //! discretization layouts supported by fluid2 (equal-order (quad or tri) or Taylor-Hood Q2Q1 (quad9))
  enum DiscretizationMode
  {
	  dismod_equalorder,		// equal-order elements (Q1Q1)
	  dismod_taylorhood	// non equal-order Taylor-Hood (Q2Q1)
  };


   /*
    * \brief get discretization mode (dismod_equalorder or dismod_taylorhood for Taylor-Hood)
    */
   inline DiscretizationMode DisMode() const {return dismode_; }


private:



  //! action parameters recognized by fluid2
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_fluid_stationary_systemmat_and_residual,
    calc_fluid_afgenalpha_systemmat_and_residual,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_genalpha_average_for_subscales_and_residual,
    get_density,
    get_gas_constant,
    integrate_shape,
    calc_gradop_and_massmatrix,
    calc_impulseeqn_implicit,
    calc_impulseeqn_semiimplicit,
    calc_fluid_residual
  };

  //! linearisation actions recognized by fluid3 (genalpha implementation)
  enum LinearisationAction
  {
    no_linearisation,
    Newton          ,
    fixed_point_like,
    minimal
  };

  //! stabilisation actions recognized by fluid2 (genalpha implementation)
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
    smagorinsky
  };

  //! options for fine-scale subgrid viscosity
  enum FineSubgridVisc
  {
    no_fssgv,
    artificial_all,
    artificial_small,
    smagorinsky_all,
    smagorinsky_small
  };

  //! map that converts strings to actions (stabilisation)
  static map<string,StabilisationAction>          stabstrtoact_;

  //! convert string to stabilisation action
  StabilisationAction ConvertStringToStabAction (const string& action) const;

  //! Gaussrule
  DRT::UTILS::GaussRule2D                         gaussrule_;

  //! flag for euler/ale net algorithm
  bool                                            is_ale_;

  //! discretization mode: equal-order element or Taylor-Hood element (nonequal)
  DiscretizationMode dismode_;

  //! container containing all kind of fluid2 thisandthat
  DRT::Container                                  data_;

  //! matrices of subscale acceleration values in gausspoints of this element
  Epetra_SerialDenseMatrix                        saccn_;
  //! matrices of subscale velocity values, current iteration value, in gausspoints of this element
  Epetra_SerialDenseMatrix                        svelnp_;
  //! matrices of subscale velocity values, last timestep, in gausspoints of this element
  Epetra_SerialDenseMatrix                        sveln_;

  // internal calculation methods

  // don't want = operator
  Fluid2& operator = (const Fluid2& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;

  /*!
    \brief integrate shape functions over volume element

    \param discretization (i) the discretisation (required in the nurbs case)
    \param lm             (i) location vector (to get ale displacements from global vectors)
    \param w              (o) shape functions integrated over element

    \return
   */
  void integrateShapefunction(
    DRT::Discretization&      discretization,
    vector<int>&              lm            ,
    Epetra_SerialDenseVector& w             );

}; // class Fluid2


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for fluid2 element

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2Register : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Fluid2Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Fluid2Register

  */
  Fluid2Register(const DRT::ELEMENTS::Fluid2Register& old);

  /*!
  \brief Deep copy this instance of Fluid2Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Fluid2Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid2Register;}

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
  virtual ~Fluid2Register();

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

  The method performs all steps needed to init the fluid2 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Fluid2Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a line of a fluid2 element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent fluid2 element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author gamnitzer (gamnitzer@lnm.mw.tum.de)
*/
class Fluid2Line : public DRT::Element
{
public:

  friend class Fluid2LineWeakDBC<DRT::Element::line2 ,DRT::Element::quad4 >;
  friend class Fluid2LineWeakDBC<DRT::Element::nurbs3,DRT::Element::nurbs9>;
  friend class Fluid2LineWeakDBCInterface;

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
  Fluid2Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Fluid2* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid2Line(const Fluid2Line& old);

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
  virtual int UniqueParObjectId() const { return ParObject_Fluid2Line;}

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
  virtual ~Fluid2Line();

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

    //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate fluid2 line element stiffness, mass, internal forces etc

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

  this method evaluates a line Neumann condition on the fluid2 element

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
                              Epetra_SerialDenseMatrix* elemat1 = NULL
			      );

  //@}

private:

  //! action parameters recognized by Fluid2Line
  enum ActionType
  {
     none,
     integrate_Shapefunction,
     calc_node_normal,
     calc_surface_tension,
     enforce_weak_dbc,
     conservative_outflow_bc,
     calc_Neumann_inflow
  };

  // don't want = operator
  Fluid2Line& operator = (const Fluid2Line& old);

  //! The parent element of this line
  DRT::ELEMENTS::Fluid2* parent_;
  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;



  // compute infintesimal line element dr for integration along the line
  double  f2_substitution(
		      const Epetra_SerialDenseMatrix  xyze,
		      const Epetra_SerialDenseMatrix  deriv,
		      const int iel);

  // Get Rule for Gaussintegration according to DRT::UTIL
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  // integrate shape functions over a line
  void IntegrateShapeFunction(ParameterList& params,
                              DRT::Discretization&       discretization,
                              vector<int>&               lm,
                              Epetra_SerialDenseVector&  elevec1,
                              const std::vector<double>& edispnp);

  /*!
    \brief Calculate node normals

    This method calculates the node normals perpendicular to the surface-element.

    \param elevec1 (out)      : Node normals multiplied by the elements length.
    \param edispnp (in)       : Displacement-vector
  */
  void ElementNodeNormal(ParameterList& params,
                         DRT::Discretization&       discretization,
                         vector<int>&               lm,
                         Epetra_SerialDenseVector&  elevec1,
                         const std::vector<double>& edispnp);

  /*!
    \brief Integrate surface tension

    \param elevec1 (out)      : RHS contribution of surface tension effect
    \param edispnp (in)       : Displacement-vector
    \param enormals (in)      : Element node normals
  */
  void ElementSurfaceTension(ParameterList& params,
                             DRT::Discretization& discretization,
                             vector<int>& lm,
                             Epetra_SerialDenseVector& elevec1,
                             const std::vector<double>& edispnp,
                             std::vector<double>& enormals);

  /*!
    \brief Evaluate weakly imposed Dirichlet conditions

    This method calculates the contributions to rhs and matrix of
    weak Dirichlet boundary conditions for a generalized alpha
    system.

    Literature:
    Weak imposition of Dirichlet boundary conditions in fluid mechanics
    Y. Bazilevs, T.J.R. Hughes
    Computers & Fluids 36 (2007) 12–26

    \param lm        (in)  : location vector of boundary element
    \param plm       (in)  : location vector of parent element
    \param elemat    (out) : element matrix contribution
    \param elevec    (out) : element vector contribution
    \param pevelaf   (in)  : current intermediate velocities on the
                             parent elements node array
    \param pevelnp   (in)  : current new velocities on the
                             parent elements node array
    \param val       (in)  : values for weak DBC
    \param functions (in)  : spatial functions to multiply weak
                             Dirichlet condition value with
    \param curvefac  (in)  : factor from timecurve
    \param Cb        (in)  : shape funtion dependent parameter
                             for penalty term
    \param wd_gamma  (in)  : switch between adjoint-(in)consistent
    \param afgdt     (in)  : time integration dependent factor, alphaF*gamma*dt
    \param gdt       (in)  : time integration dependent factor, gamma*dt
    \param comp_lin  (in)  : linearisation of adjoint convective fluxes on/off

  */
  void EvaluateWeakDirichlet(
  vector<int>&               lm       ,
  vector<int>&               plm      ,
  Epetra_SerialDenseMatrix&  elemat   ,
  Epetra_SerialDenseVector&  elevec   ,
  Epetra_SerialDenseMatrix&  pevelaf  ,
  Epetra_SerialDenseMatrix&  pevelnp  ,
  Epetra_SerialDenseVector&  peprenp  ,
  vector<double>             val      ,
  const vector<int>*         functions,
  double                     curvefac ,
  double                     Cb       ,
  double                     wd_gamma ,
  const double               afgdt    ,
  const double               gdt      ,
  const bool                 comp_lin );


  /*!
  \brief apply outflow boundary condition which is necessary for the conservative
         element formulation (since the convective term was partially integrated)

  \param params (in)        : ParameterList for communication between control routine
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element
  \param elevec1 (out)      : vector to be filled by element

  */
  void LineConservativeOutflowConsistency(
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

}; // class Fluid2Line




} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef FLUID2_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
