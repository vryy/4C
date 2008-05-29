/*!
\file xfluid3.H

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


#ifdef PERF
#include <Teuchos_TimeMonitor.hpp>
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

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
class XFluid3Register;

/*!
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
  virtual int NumDofPerNode(const DRT::Node& node) const {return this->eleDofManager_.NumDofPerNode(node.Id(), 4);}

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
  virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const;

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

  //@}


private:

  //! action parameters recognized by fluid3
  enum ActionType
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_linear_fluid,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_genalpha_average_for_subscales_and_residual,
    calc_fluid_stationary_systemmat_and_residual,
    calc_fluid_beltrami_error,
    calc_turbulence_statistics,
    calc_fluid_box_filter,
    calc_smagorinsky_const,
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
    viscous_stab_usfem,
    viscous_stab_usfem_only_rhs,
    continuity_stab_none,
    continuity_stab_yes,
    continuity_stab_td,
    cross_stress_stab_none,
    cross_stress_stab,
    cross_stress_stab_only_rhs,
    reynolds_stress_stab_none,
    reynolds_stress_stab_only_rhs,
    fssgv_no,
    fssgv_artificial_all,
    fssgv_artificial_small,
    fssgv_Smagorinsky_all,
    fssgv_Smagorinsky_small,
    fssgv_mixed_Smagorinsky_all,
    fssgv_mixed_Smagorinsky_small,
    fssgv_scale_similarity
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

  //! flag for euler/ale net algorithm
  bool                                            is_ale_;
  //! container containing all kind of fluid3 thisandthat
  DRT::Container                                  data_;

  //! vector of volumes of this element (length 1)
  std::vector<DRT::Element*>                                volume_;

  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  std::vector<Teuchos::RCP<DRT::ELEMENTS::XFluid3Surface> > surfaces_;
  //! vector of surfaces adjacent to this element (fluid3: 4 or 6 surfaces)
  std::vector<DRT::ELEMENTS::XFluid3Surface*>               surfaceptrs_;
  //! vector of lines adjacent to this element
  std::vector<Teuchos::RCP<DRT::ELEMENTS::XFluid3Line> >    lines_;
  //! vector of lines adjacent to this element
  std::vector<DRT::ELEMENTS::XFluid3Line*>                  lineptrs_;

private:

  //! store information about the unknowns of this element
  XFEM::ElementDofManager                                   eleDofManager_;
  
  //! store information about the intersection handle
  Teuchos::RCP<XFEM::InterfaceHandle>                       ih_;

  // internal calculation methods

  // don't want = operator
  XFluid3& operator = (const XFluid3& old);

  /*!
   * \brief evaluate error for beltrami test problem
   */
  void f3_int_beltrami_err(
    vector<double>&       myvelnp,
    vector<double>&       myprenp,
    const struct _MATERIAL*     material,
    ParameterList&    	  params
    );

  /*!
   \brief get optimal gaussrule for given discretization type
   \return enum with optimal gausrule for this element
   */
  DRT::UTILS::GaussRule3D getOptimalGaussrule(
    const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
    );

  ActionType convertStringToActionType(const string& action) const;

  //! convert string to stabilisation action
  StabilisationAction ConvertStringToStabAction (const string& action) const;

}; // class XFluid3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
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

private:

  //! action parameters recognized by Fluid3Surface
  enum ActionType
  {
    none,
    integrate_Shapefunction
  };

  // don't want = operator
  XFluid3Surface& operator = (const XFluid3Surface& old);

  //! The parent element of this surface
  DRT::ELEMENTS::XFluid3* parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  const int                    lsurface_;
  //! vector of lines adjacent to this element
  std::vector<Teuchos::RCP<DRT::ELEMENTS::XFluid3Line> > lines_;
  //! vector of lines adjacent to this element
  std::vector<DRT::ELEMENTS::XFluid3Line*>              lineptrs_;

  //  compute kovariant metric tensor G for fluid surface element
  //                                                  gammi 04/07

  void  f3_metric_tensor_for_surface(
    const Epetra_SerialDenseMatrix  xyze,
    const Epetra_SerialDenseMatrix  deriv,
    Epetra_SerialDenseMatrix&       metrictensor,
    double                         *drs);

  /*!
  \brief integrate shapefunction over surface element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
  */

  virtual void IntegrateShapeFunction(
          ParameterList&                params,
          DRT::Discretization&          discretization,
          vector<int>&                  lm,
          Epetra_SerialDenseVector&     elevec1,
          const std::vector<double>&    edispnp);


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

  XFluid3Line(  int id,
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
     none
  };

  // don't want = operator
  XFluid3Line& operator = (const XFluid3Line& old);

  //! The volume parent element of this line
  DRT::Element* parent_;




  //! The local line number of this line w.r.t to the parent_ element
  int                    lline_;


}; // class XFluid3Line


} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef XFLUID3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
