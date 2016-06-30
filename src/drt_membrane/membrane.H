/*!----------------------------------------------------------------------
\file membrane.H

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

\brief Nonlinear Membrane Finite Element

The input line (in the header file) should read
MAT x KINEM nonlinear THICK x STRESS_STRAIN [plane_stress/plane_strain]


*----------------------------------------------------------------------*/
#ifndef MEMBRANE_H
#define MEMBRANE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_so3/so_base.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "membrane_eletypes.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
template<DRT::Element::DiscretizationType distype2>
class MembraneLine;

/*!
\brief A C++ wrapper for the membrane element
*/
template<DRT::Element::DiscretizationType distype>
class Membrane : public So_base
{

public:

  //! @name Friends
  friend class Membrane_tri3Type;
  friend class Membrane_tri9Type;
  friend class Membrane_quad4Type;
  friend class Membrane_quad9Type;
  template<DRT::Element::DiscretizationType distype2> friend class MembraneLine;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */
  Membrane(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Membrane(const Membrane<distype>& old);

  /*!
  \brief Deep copy this instance of Membrane and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-constructor is needed

  */
  DRT::Element* Clone() const;

  //! number of element nodes
  static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int numdim_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of dofs per node
  static const int noddof_ = 3;

  //! total dofs per element
  static const int numdof_ = noddof_ * numnod_;

  //!total gauss points per element
  static const int numgpt_ = numnod_;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const;

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const
  {
    switch(distype)
    {
    case DRT::Element::tri3:
    {
      return Membrane_tri3Type::Instance().UniqueParObjectId();
    }
    case DRT::Element::tri6:
    {
      return Membrane_tri6Type::Instance().UniqueParObjectId();
    }
    case DRT::Element::quad4:
    {
      return Membrane_quad4Type::Instance().UniqueParObjectId();
    }
    case DRT::Element::quad9:
    {
      return Membrane_quad9Type::Instance().UniqueParObjectId();
    }
    default: dserror("unknown element type!");
      break;
    }
    // Intel compiler needs a return so
    return -1;
  };

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Membrane();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return noddof_;}

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
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  {
    switch(distype)
    {
    case DRT::Element::tri3:
    {
       return Membrane_tri3Type::Instance();
    }
    break;
    case DRT::Element::tri6:
    {
       return Membrane_tri6Type::Instance();
    }
    break;
    case DRT::Element::quad4:
    {
       return Membrane_quad4Type::Instance();
    }
    break;
    case DRT::Element::quad9:
    {
       return Membrane_quad9Type::Instance();
    }
    break;
    default:
      dserror("unknown element type!");
    break;
    }
    // Intel compiler needs a return so
    return Membrane_quad4Type::Instance();
  };

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
    names.insert(std::pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(std::pair<string,int>("StressesXYZ",6));
  \endcode

  \param names (out): On return, the derived class has filled names with
                      key names of data it wants to visualize and with int dimensions
                      of that data.
  */
  virtual void VisNames(std::map<std::string,int>& names);

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
  virtual bool VisData(const std::string& name, std::vector<double>& data);

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& eledistype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate Membrane element stiffness, mass, internal forces etc

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param discretization : pointer to discretization for de-assembly
  \param lm (in)        : location matrix for de-assembly
  \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                          the controlling method does not expect the element to fill
                          this matrix.
  \param elemat2 (out)  : (mass-)matrix to be filled by element. If NULL on input,
                          the controlling method does not expect the element to fill
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
  int Evaluate(Teuchos::ParameterList& params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1_epetra,
               Epetra_SerialDenseMatrix& elemat2_epetra,
               Epetra_SerialDenseVector& elevec1_epetra,
               Epetra_SerialDenseVector& elevec2_epetra,
               Epetra_SerialDenseVector& elevec3_epetra);


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surfaces Neumann condition on the membrane element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  //@}

  //! @name Other

  /*!
  \brief Get ptr to nodal thicknesses
  */
  const std::vector<double>* thickness() { return data_.Get<std::vector<double> >("thick");}

  //@}


private:

  //! action parameters recognized by membrane
  enum ActionType
  {
    none,
    calc_struct_nlnstiff,
    calc_struct_nlnstiffmass,
    calc_struct_update_istep,
    calc_struct_reset_istep,
    calc_struct_stress,
    postprocess_stress,
    calc_cur_normal_at_point
  };

  /// type of 2D dimension reduction
  enum DimensionalReduction
  {
    plane_stress,  ///< plane stress, i.e. lateral stress is zero \f$S_{33}=S_{13}=S_{23}=0\f$

    //Membrane not intended for plane strain evaluation (mentioned here for completeness)
    plane_strain   ///< plane strain, i.e. lateral strain is zero \f$E_{33}=E_{13}=E_{23}=0\f$
  };

  //! membrane thickness
  double                                          thickness_;

  //! current membrane thickness at gauss point
  LINALG::Matrix<numgpt_,1>                       curr_thickness_;

  //! container containing all kind of membrane thisandthat
  DRT::Container                                  data_;

  //! gaussian integration to be used
  DRT::UTILS::GaussRule2D                         gaussrule_;

  //! membrane stress/strain state
  DimensionalReduction                            planetype_;

  //! current time
  double                                          time_;


  // internal calculation methods

  // don't want = operator
  Membrane<distype>& operator = (const Membrane<distype>& old);

  //! calculate nonlinear stiffness and mass matrix
  void mem_nlnstiffmass(std::vector<int>&                      lm,                    // location matrix
                        std::vector<double>&                   disp,                  // current displacements
                        LINALG::Matrix<numdof_,numdof_>*       stiffmatrix,           // element stiffness matrix
                        LINALG::Matrix<numdof_,numdof_>*       massmatrix,            // element mass matrix
                        LINALG::Matrix<numdof_,1>*             force,                 // element internal force vector
                        LINALG::Matrix<numgpt_,6>*             elestress,             // stresses at GP
                        LINALG::Matrix<numgpt_,6>*             elestrain,             // strains at GP
                        Teuchos::ParameterList&                params,                // algorithmic parameters e.g. time
                        const INPAR::STR::StressType           iostress,              // stress output option
                        const INPAR::STR::StrainType           iostrain);             // strain output option

  //! get reference and current configuration
  void mem_configuration(const std::vector<double>&             disp,
                         LINALG::Matrix<numnod_,noddof_>&       xrefe,
                         LINALG::Matrix<numnod_,noddof_>&       xcurr);

  //! introduce orthonormal base in the undeformed configuration at current Gauss point
  void mem_orthonormalbase(const LINALG::Matrix<numnod_,noddof_>&                      xrefe,
                           const LINALG::Matrix<numnod_,noddof_>&                      xcurr,
                           const LINALG::Matrix<numdim_,numnod_>&                      derivs,
                           LINALG::Matrix<numdim_,numnod_>&                            derivs_ortho,
                           double&                                                     G1G2_cn,
                           LINALG::Matrix<noddof_,1>&                                  dXds1,
                           LINALG::Matrix<noddof_,1>&                                  dXds2,
                           LINALG::Matrix<noddof_,1>&                                  dxds1,
                           LINALG::Matrix<noddof_,1>&                                  dxds2,
                           LINALG::Matrix<noddof_,noddof_>&                            Q_trafo
                           );

  //! pushforward of 2nd Piola-Kirchhoff stresses to Cauchy stresses at Gauss point
  void mem_PK2toCauchy(const LINALG::Matrix<3,3>&                      pkstress_global,
                       const LINALG::Matrix<noddof_,noddof_>&          defgrd,
                       LINALG::Matrix<3,3>&                            cauchy);

  // pushforward of Green-Lagrange to Euler-Almansi strains at Gauss point
  void mem_GLtoEA(const LINALG::Matrix<3,3>&                      glstrain_global,
                  const LINALG::Matrix<noddof_,noddof_>&          defgrd,
                  LINALG::Matrix<3,3>&                            euler_almansi);

  // transforms local membrane surface tensor to tensor in global coordinates
  void mem_localtoglobal(const LINALG::Matrix<noddof_,noddof_>&    Q_trafo,
                         const LINALG::Matrix<noddof_,noddof_>&    local,
                         LINALG::Matrix<noddof_,noddof_>&          global);

  // transforms tensor in global coords to local membrane surface tensor
  void mem_globaltolocal(const LINALG::Matrix<noddof_,noddof_>&    Q_trafo,
                         const LINALG::Matrix<noddof_,noddof_>&    global,
                         LINALG::Matrix<noddof_,noddof_>&          local);

  // evaluate 3D Material law and reduce to plane stress state
  void mem_Material3dPlane(const LINALG::Matrix<noddof_,1>&          dXds1,
                           const LINALG::Matrix<noddof_,1>&          dXds2,
                           const LINALG::Matrix<noddof_,1>&          dxds1,
                           const LINALG::Matrix<noddof_,1>&          dxds2,
                           LINALG::Matrix<noddof_,noddof_>&          defgrd,
                           LINALG::Matrix<noddof_,noddof_>&          cauchygreen_local,
                           LINALG::Matrix<3,1>&                      pkstress,
                           LINALG::Matrix<3,3>&                      cmat,
                           const LINALG::Matrix<noddof_,noddof_>&    Q_trafo,
                           Teuchos::ParameterList&                   params);

  // transformation matrix for material tangent tensor
  void mem_cmat_trafo(const LINALG::Matrix<noddof_,noddof_>&    Q_trafo,
                      LINALG::Matrix<6,6>&                      cmat_trafo);

  // transforms material tangent tensor from global coordinates to local membrane surface
  void mem_cmat_globaltolocal(const LINALG::Matrix<6,6>&    cmat_trafo,
                              const LINALG::Matrix<6,6>&    cmat_global,
                              LINALG::Matrix<6,6>&          cmat_local);

  // transforms matrix to vector in stress-like voigt notation
  void mem_matrixtovoigt(const LINALG::Matrix<3,3>&    matrix,
                         LINALG::Matrix<6,1>&          voigt);

  // transforms vector in stress-like voigt notation to matrix
  void mem_voigttomatrix(const LINALG::Matrix<6,1>&    voigt,
                         LINALG::Matrix<3,3>&          matrix);

  // determine deformation gradient in global frame on membrane surface
  void mem_defgrd_global(const LINALG::Matrix<noddof_,1>&          dXds1,
                         const LINALG::Matrix<noddof_,1>&          dXds2,
                         const LINALG::Matrix<noddof_,1>&          dxds1,
                         const LINALG::Matrix<noddof_,1>&          dxds2,
                         const double&                             lambda3,
                         LINALG::Matrix<noddof_,noddof_>&          defgrd_global);

}; // class Membrane


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  LINE 2 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_line2Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_line2Type"; }

  static Membrane_line2Type & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Membrane_line2Type instance_;
};

/*----------------------------------------------------------------------*
 |  LINE 3 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_line3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_line3Type"; }

  static Membrane_line3Type & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Membrane_line3Type instance_;
};

/*!
\brief An element representing a line edge of a membrane element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate line Neumann boundary conditions that might be
      adjacent to a parent membrane element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

\author braeu (braeu@lnm.mw.tum.de)
*/
template<DRT::Element::DiscretizationType distype2>
class MembraneLine : public DRT::FaceElement
{
public:

  //! @name Friends
  friend class Membrane_line2Type;
  friend class Membrane_line3Type;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes
  \param parent: The parent shell element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  MembraneLine(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes, DRT::ELEMENTS::Membrane<distype2>* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  MembraneLine(const MembraneLine<distype2>& old);

  //! number of nodes per line
  static const int numnod_line_ = DRT::UTILS::DisTypeToNumNodePerEle<DRT::UTILS::DisTypeToFaceShapeType<distype2>::shape>::numNodePerElement;
  //! total gauss points per line
  static const int numgpt_line_ = numnod_line_;

  static const int noddof_ = 3;

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
  virtual int UniqueParObjectId() const
  {
    switch(DRT::UTILS::DisTypeToFaceShapeType<distype2>::shape)
    {
    case DRT::Element::line2:
    {
      return Membrane_line2Type::Instance().UniqueParObjectId();
    }
    case DRT::Element::line3:
    {
      return Membrane_line3Type::Instance().UniqueParObjectId();
    }
    default: dserror("unknown line type!");
      break;
    }
    // Intel compiler needs a return so
    return -1;
  };


  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~MembraneLine();

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
  virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

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
   * \brief Return pointer to the parent element
  */
  virtual DRT::ELEMENTS::Membrane<distype2>* ParentElement() const
  {
    DRT::Element* parent = this->DRT::FaceElement::ParentElement();
    // make sure the static cast below is really valid
    dsassert(dynamic_cast<DRT::ELEMENTS::Membrane<distype2>* >(parent) != NULL, "Parent element is no membrane element");
    return static_cast<DRT::ELEMENTS::Membrane<distype2>* >(parent);
  }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  {
    switch(DRT::UTILS::DisTypeToFaceShapeType<distype2>::shape)
    {
    case DRT::Element::line2:
    {
       return Membrane_line2Type::Instance();
    }
    case DRT::Element::line3:
    {
       return Membrane_line3Type::Instance();
    }
    default:
      dserror("unknown line type!");
    break;
    }
    // Intel compiler needs a return so
    return Membrane_line2Type::Instance();
  };

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the membrane element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

private:

  // don't want = operator
  MembraneLine<distype2>& operator = (const MembraneLine<distype2>& old);

  //! gaussian integration to be used
  DRT::UTILS::GaussRule1D                         gaussrule_;

}; // class MembraneLine

} // namespace ELEMENTS
} // namespace DRT

#endif
