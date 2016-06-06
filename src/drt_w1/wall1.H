/*----------------------------------------------------------------------*/
/*!
\file wall1.H

\brief ToDo Add meaningful comment.

\level 1

\maintainer Markus Gitterle

*/
/*----------------------------------------------------------------------*/
/* macros */
#ifndef WALL1_H
#define WALL1_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_so3/so_base.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class Wall1Line;

class Wall1Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Wall1Type"; }

  static Wall1Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Wall1Type instance_;
};

/*======================================================================*/
/*!
 * \brief A C++ wrapper for the wall1 element
 *
 * A planar solid element for geometrically non-linear analyses.
 *
 * \author gitterle (gitterle@lnm.mw.tum.de)
 *
 * <b>Typical quantities and their storage pattern</b>
 * - Green-Lagrange strain vector
 *   \f$E = \left[\begin{array}{c} E_{11} \\ E_{22} \\ 2E_{12} \end{array}\right]\f$
 * - Deformation gradient vector
 *   \f$\tilde{F} = \left[\begin{array}{c} F_{11} \\ F_{22} \\ F_{12} \\ F_{21} \end{array}\right]\f$
 * - Constitutive matrix
 *   \f$C = \left[\begin{array}{ccc} C_{11} & C_{12} & C_{13} \\ C_{21} & C_{22} & C_{23} \\ C_{31} & C_{32} & C_{33} \end{array}\right]\f$
 * - Second Piola--Kirchhoff stress vector
 *   \f$S = \left[\begin{array}{c} S_{11} \\ S_{22} \\ S_{12} \end{array}\right]\f$
 * - Second Piola--Kirchhoff stress matrix
 *   \f$\bar{S} = \left[\begin{array}{cccc} S_{11} & 0 & S_{12} & 0 \\ 0 & S_{22} & 0 & S_{12} \\ S_{12} & 0 & S_{22} & 0 \\ 0 & S_{12} & 0 & S_{11} \end{array}\right]\f$
 * - First Piola-Kirchhoff stress vector
 *   \f$\tilde{P} = \left[\begin{array}{c} P_{11} \\ P_{22} \\ P_{12} \\ P_{21} \end{array}\right]\f$
 *
 * <b>References</b>
 * - [1] JC Simo and F Armero, Geometrically non-linear enhanced strain
 *   mixed methods and the method of incompatible modes, International
 *   Journal for Numerical Methods in Engineering, 33:1413-1449, 1992.
 * - [2] S Glaser and F Armero, On the formulation of enhanced strain finite
 *   elements in finite deformations, Engineering Computations, 14:759-791,
 *   1997.
 * - [3] P Wriggers, Nichtlineare Finite-Element-Methoden, Springer, 2001.
 * - [4] WA Wall and B Bornemann, Nichtlineare Finite-Element-Methoden,
 *   Vorlesungsskript, Lehrstuhl fuer Numerische Mechanik, SS 2007.
 */
class Wall1 : public So_base
{
public:

  /// @name Friends
  //@{
  friend class Wall1Line;
  //@}

  /// @name object-wide constants
  //@{
  static const int maxnod_ = 9;  ///< maximally permitted number of nodes
  static const int numdim_ = 2;  ///< 2-dimensional/planar element
  static const int noddof_ = 2;  ///< number of nodal DOFs
  static const int numstr_ = 4;  ///< number of (symmetric) strain/stress components --- THIS SHOULD BE 3
  static const int numnstr_ = 4; ///< number of non-symmetric strain/stress components
  static const int neas_ = 4;    ///< number of EAS parameters
  //@}

  /// @name Constructors and destructors and related methods

  /// Standard Constructor
  Wall1(
    int id,  ///< A unique global id
    int owner
  );

  /// Copy Constructor
  ///
  /// Makes a deep copy of a Element
  Wall1(const Wall1& old);

  /// Deep copy this instance of Wall1 and return pointer to the copy
  ///
  /// The Clone() method is used from the virtual base class Element in cases
  /// where the type of the derived class is unknown and a copy-ctor is needed
  DRT::Element* Clone() const;

  /// Get shape type of element
  virtual DiscretizationType Shape() const;

  /// Set discretization type of element
  virtual void SetDisType(DRT::Element::DiscretizationType shape)
  {
    distype_ = shape;
    return;
  };


  /// Return number of lines of this element
  virtual int NumLine() const
  {
    if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4;
    else return 3;
  }

  /// Return number of surfaces of this element
  virtual int NumSurface() const { return 1; }

  /// Get vector of Teuchos::RCPs to the lines of this element
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /// Get vector of Teuchos::RCPs to the surfaces of this element
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /// Return unique ParObject id
  ///
  /// every class implementing ParObject needs a unique id defined at the
  /// top of this file.
  virtual int UniqueParObjectId() const { return Wall1Type::Instance().UniqueParObjectId(); }

  /// Pack this class so it can be communicated
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Pack(DRT::PackBuffer& data) const;

  /// Unpack data from a char vector into this class
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Unpack(const std::vector<char>& data);

  /// Destructor
  virtual ~Wall1();

  /// get Cauchy stress
  virtual void GetCauchyAtXi(const LINALG::Matrix<2,1>& xi,
        const std::vector<double>& disp,
        LINALG::Matrix<2,2>& cauchy,
        Epetra_SerialDenseMatrix& dsdd);
  //@}

  /// @name Acess methods
  //@{
    /*!
    \brief Does this element use EAS?
    */
    bool HaveEAS() const { return (eastype_!=eas_vague); };

  /// Get number of degrees of freedom of a certain node
  /// (implements pure virtual DRT::Element)
  ///
  /// The element decides how many degrees of freedom its nodes must have.
  /// As this may vary along a simulation, the element can redecide the
  /// number of degrees of freedom per node along the way for each of it's nodes
  /// separately.
  virtual int NumDofPerNode(const DRT::Node& node) const { return Wall1::noddof_; }

  /// Get number of degrees of freedom per element
  /// (implements pure virtual DRT::Element)
  ///
  /// The element decides how many element degrees of freedom it has.
  /// It can redecide along the way of a simulation.
  ///
  /// \note Element degrees of freedom mentioned here are dofs that are visible
  ///      at the level of the total system of equations. Purely internal
  ///      element dofs that are condensed internally should NOT be considered.
  virtual int NumDofPerElement() const { return 0; }

  /// Print this element
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return Wall1Type::Instance(); }

  //@}

  /// @name Input and Creation
  //@{

  /// Read input for this element
  virtual bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  /// @name Evaluation
  //@{

  /// Evaluate an element
  ///
  /// Evaluate Wall1 element stiffness, mass, internal forces etc
  ///
  /// \return 0 if successful, negative otherwise
  virtual int Evaluate(
    Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication between control routine
                            ///< and elements
    DRT::Discretization& discretization,  ///< A reference to the underlying discretization
    std::vector<int>& lm,  ///< location vector of this element
    Epetra_SerialDenseMatrix& elemat1,  ///< matrix to be filled by element. If NULL on input,
                                        ///< the controling method does not epxect the element to fill
                                        ///< this matrix.
    Epetra_SerialDenseMatrix& elemat2,  ///< matrix to be filled by element. If NULL on input,
                                        ///< the controling method does not epxect the element to fill
                                        ///< this matrix.
    Epetra_SerialDenseVector& elevec1,  ///< vector to be filled by element. If NULL on input,
                                        ///< the controlling method does not epxect the element
                                        ///< to fill this vector
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element. If NULL on input,
                                        ///< the controlling method does not epxect the element
                                        ///< to fill this vector
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element. If NULL on input,
                                       ///< the controlling method does not epxect the element
                                       ///< to fill this vector
  );

  /// Evaluate a Neumann boundary condition
  ///
  /// this method evaluates a surfaces Neumann condition on the wall element
  ///
  /// \return 0 if successful, negative otherwise
  virtual int EvaluateNeumann(
    Teuchos::ParameterList& params,   ///< (in/out) ParameterList for communication between control routine
                             ///< and elements
    DRT::Discretization& discretization,  ///< A reference to the underlying discretization
    DRT::Condition& condition,  ///<  The condition to be evaluated
    std::vector<int>& lm,  ///< location vector of this element
    Epetra_SerialDenseVector& elevec1,  ///< vector to be filled by element. If NULL on input
    Epetra_SerialDenseMatrix* elemat1 = NULL);

  /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store them in elevecs
  void w1_expol(
      Epetra_SerialDenseMatrix& stresses, ///< gp stresses
      Epetra_MultiVector&       expolstress ///< expol stresses
  );

  //@}

protected:

  /// action parameters recognized by wall1
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,
    calc_struct_nlnstiff_gemm,  ///< internal force, stiffness and mass
                                ///< for generalised energy-momentum method (GEMM)
    calc_struct_stress,
    postprocess_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_reset_istep, ///< reset elementwise internal variables
                             ///< during iteration to last converged state
    calc_struct_energy,      ///< calculate internal energy of element
    calc_struct_errornorms,  ///< compute error norms (L2,H1,energy)
    calc_potential_stiff,    ///< calculate volume potential

    calc_struct_mass_volume  ///< calc volume and mass
  };

  /// type of 2D dimension reduction
  enum DimensionalReduction
  {
    plane_none,  ///< undetermined
    plane_stress,  ///< plane stress, i.e. lateral stress is zero \f$S_{33}=S_{13}=S_{23}=0\f$
    plane_strain  ///< plane strain, i.e. lateral strain is zero \f$E_{33}=E_{13}=E_{23}=0\f$
  };

  /// type of stresses to be calculatecd in postprocessing
  ///
  /// describes in co-ordinate system the stress components are expressed
  enum StressType
  {
    w1_none,  ///< undetermined
    w1_xy,  ///< Cartesian
    w1_rs  ///< local/parametric/natural
  };

  /// EAS type
  enum EasType
  {
    eas_vague,  ///< unknown type
    eas_q1e4,  ///< Q1E4 formulation due to Simo and Armero [1]
    eas_q1et4  ///< Q1ET4 formulation due to Glaser and Armero [2]
  };

  /// container containing all kind of wall1 thisandthat
  DRT::Container data_;

  /// number of the material law
  int material_;
  /// element thickness
  double thickness_;
  // line search parameter (old step length)
  double old_step_length_;
  /// gaussian points
  DRT::UTILS::GaussRule2D gaussrule_;
  /// problem type
  DimensionalReduction wtype_;
  /// type of stress calculation
  StressType stresstype_;
  /// eas or not
  bool iseas_;
  /// EAS type
  enum EasType eastype_;
  //! struct_ale
  bool structale_;

  //! the element discretization type (shape)
  DRT::Element::DiscretizationType distype_;

  /// @name Internal calculation methods
  //@{
  /** recover elementwise stored stuff
   *
   * \author hiermeier
   * \date 05/16 */
  void w1_recover(
    const std::vector<int>&    lm,
    const std::vector<double>& disp,
    const std::vector<double>& residual);

  /// evaluate the element forces and stiffness and mass
  /// \author mgit \date 03/07
  void w1_nlnstiffmass(
    const std::vector<int>               & lm         , ///< location vector
    const std::vector<double>            & disp       , ///< element displacements
    const std::vector<double>            & residual   , ///< residual displacements
    const std::vector<double>            & dispmat    , ///< residual displacements
    std::vector<Epetra_SerialDenseVector>& myknots    , ///< knot vector for nurbs elements
    Epetra_SerialDenseMatrix             * stiffmatrix, ///< element stiffness matrix
    Epetra_SerialDenseMatrix             * massmatrix , ///< element mass matrix
    Epetra_SerialDenseVector             * force      , ///< element internal force vector
    Epetra_SerialDenseMatrix             * elestress  , ///< element stresses
    Epetra_SerialDenseMatrix             * elestrain  , ///< element strains
    Teuchos::RCP<const MAT::Material>      material   , ///< element material
    Teuchos::ParameterList&                params     , ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType           iostress   , ///< stress output option
    const INPAR::STR::StrainType           iostrain     ///< strain output option
  );

  /// evaluate the geometrically linear element forces and stiffness and mass
  void w1_linstiffmass(
    const std::vector<int>               & lm         , ///< location vector
    const std::vector<double>            & disp       , ///< element displacements
    const std::vector<double>            & residual   , ///< residual displacements
    const std::vector<double>            & dispmat    , ///< residual displacements
    std::vector<Epetra_SerialDenseVector>& myknots    , ///< knot vector for nurbs elements
    Epetra_SerialDenseMatrix             * stiffmatrix, ///< element stiffness matrix
    Epetra_SerialDenseMatrix             * massmatrix , ///< element mass matrix
    Epetra_SerialDenseVector             * force      , ///< element internal force vector
    Epetra_SerialDenseMatrix             * elestress  , ///< element stresses
    Epetra_SerialDenseMatrix             * elestrain  , ///< element strains
    Teuchos::RCP<const MAT::Material>      material   , ///< element material
    Teuchos::ParameterList&                params     , ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType           iostress   , ///< stress output option
    const INPAR::STR::StrainType           iostrain     ///< strain output option
  );

  /// Jacobian matrix for mapping from parameter space in physical material space
  /// at point parameter space
  /// \author mgit \date 04/07
  void w1_jacobianmatrix(
    const Epetra_SerialDenseMatrix& xrefe,  ///< reference/material co-ordinates of element nodes
    const Epetra_SerialDenseMatrix& deriv,  ///< derivatives of shape functions at parameter point
    Epetra_SerialDenseMatrix& xjm,  ///< Jacobi matrix
    double* det,  ///< determinant of Jacobi matrix
    const int iel  ///< actual number of element nodes
  );

  /// Linear B-operator in reference configuration at point parameter space
  /// \author mgit \date 04/07
  void w1_boplin(
    Epetra_SerialDenseMatrix& boplin,  ///< the B-operator
    Epetra_SerialDenseMatrix& deriv,  ///< derivatives of shape functions at parameter point
    Epetra_SerialDenseMatrix& xjm,  ///< Jacobian at parameter point
    double& det,  ///< Jacobian determinant at parameter point
    const int iel  ///< number of element nodes
  );

  /// (Material) Deformation gradient \f$F\f$ and Green-Lagrange strains \f$E\f$
  /// at parameter point
  /// \author mgit \date 04/07
  void w1_defgrad(
    Epetra_SerialDenseVector& F,  ///< deformation gradient
    Epetra_SerialDenseVector& strain,  ///< GL strain \f$E^T=[E_{11} \; E_{22} \; E_{12}]\f$
    const Epetra_SerialDenseMatrix& xrefe,  ///< reference/material co-ordinates of element nodes
    const Epetra_SerialDenseMatrix& xcure, ///< current/spatial co-ordinates of element nodes
    Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const int iel  ///< number of element nodes
  );

  /// Deformation gradient, measures distortion of the mesh of the
  /// material configuration with respect to the referential configuration
  /// in structure with ale approaches (fractional step method)
  /// \author mgit \date 04/11
  void w1_defgradmat(
    Epetra_SerialDenseVector& F,  ///< deformation gradient (mesh distortion of material configuration)
    Epetra_SerialDenseVector& Fmat,  ///< deformation gradient (mesh distortion of spatial configuration)
    Epetra_SerialDenseVector& FFmatinv,  ///< product of F and Fmat(inv)
    Epetra_SerialDenseVector& strain,  ///< GL strain \f$E^T=[E_{11} \; E_{22} \; E_{12}]\f$
    const Epetra_SerialDenseMatrix& xrefe,  ///< reference coordinates of element nodes
    const Epetra_SerialDenseMatrix& xmat, ///< material coordinates of element nodes
    Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const int iel  ///< number of element nodes
  );

  /// Non-linear B-operator in reference configuration
  /// \author mgit \date 04/07
  void w1_boplin_cure(
    Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-operator
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const Epetra_SerialDenseVector& F,  ///< deformation gradient as Voigt-vector
    const int numeps,  ///< number of strains
    const int nd  ///< number of element nodes
  );

  /// Geometric stiffness constribution (total Lagrange)
  /// \author mgit \date 05/07
  void w1_kg(
    Epetra_SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    const double fac,  ///< integration factor
    const int nd,  ///< number of element DOFs
    const int numeps  ///< number of strains
  );

  /// elastic and initial displacement stiffness contribution (total Lagrange)
  /// \author mgit \date 05/07
  void w1_keu(
    Epetra_SerialDenseMatrix& estif,  ///< (in/out) element stiffness matrix
    const Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-operator
    const Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
    const double fac,  ///< integration factor
    const int nd,  ///< number of element DOFs
    const int numeps  ///< number of strains
  );

  /// Evaluate internal element forces for large def (total Lagr)
  /// \author mgit \date 05/07
  void w1_fint(
    const Epetra_SerialDenseMatrix& stress, ///< PK2 stress vector
    const Epetra_SerialDenseMatrix& b_cure,  ///< non-linear B-op
    Epetra_SerialDenseVector& intforce,  ///< internal force vector
    const double fac,  ///< integration factor
    const int nd  ///< number of element DOFs
  );

  /// lump mass matrix
  void w1_lumpmass(
    Epetra_SerialDenseMatrix* emass  ///< (in/out) element mass matrix
  );

  /// determine cauchy stress and store it
  void StressCauchy(
    const int ip,  ///< Gauss point index
    const double& F11,  ///< F_{11} component of def.grad.
    const double& F22,  ///< F_{22} component of def.grad.
    const double& F12,  ///< F_{12} component of def.grad.
    const double& F21,  ///< F_{21} component of def.grad.
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress matrix
    Epetra_SerialDenseMatrix* elestress  ///< element stress array
  );

  /// determine internal energy
  void Energy(
    Teuchos::ParameterList& params,  ///<  a list containing GEMM coefficients
    const std::vector<int>& lm,  ///< location vector
    const std::vector<double>& dis,  ///< element displacements \f$d_{n}^{(e)}\f$ at \f$t_{n}\f$
    Epetra_SerialDenseVector* energies,  ///< (in/out) energies
    Teuchos::RCP<const MAT::Material> material  ///< element material
  );

  //@}

  /// @name EAS tools
  //@{
  /// return the number of eas dofs
  int w1_neas() const
  { return neas_; };

  /// setup of constant EAS data
  void w1_eassetup(
    Epetra_SerialDenseMatrix& boplin0,  ///< linear B-op at origin
    Epetra_SerialDenseVector& F0,  ///< deformation gradient at origin
    Epetra_SerialDenseMatrix& xjm0,  ///< jacobian matrix at origin
    double& detJ0,  ///< det of Jacobian at origin
    const Epetra_SerialDenseMatrix& xrefe,  ///< material element coords
    const Epetra_SerialDenseMatrix& xcure,  ///< current element coords
    const DRT::Element::DiscretizationType& distype  ///< discretisation type
  );

  /// Get the enhanced deformation gradient and
  /// also the operators G,W0 and Z
  /// at point in parameter space
  /// \author mgit \date 01/08
  void w1_call_defgrad_enh(
    Epetra_SerialDenseMatrix& F_enh,  ///< enhanced deformation gradient vector
    const Epetra_SerialDenseMatrix xjm0,  ///< Jacobian at origin
    const Epetra_SerialDenseMatrix xjm,  ///< Jacobian at parameter point
    const double detJ0,  ///< det of Jacobian at origin
    const double det,  ///< det of Jacobian at parameter point
    const Epetra_SerialDenseVector F0,  ///< deformation gradient at origin
    const Epetra_SerialDenseMatrix alpha,  ///< alpha parameters (EAS params)
    const double e1,  ///< \f$\xi\f$ co-ordinate of parameter point
    const double e2,  ///< \f$\eta\f$ co-ordinate of parameter point
    Epetra_SerialDenseMatrix& G,  ///< G-operator
    Epetra_SerialDenseMatrix& W0,  ///< W_0-operator
    const Epetra_SerialDenseMatrix boplin0,  ///< linear B-op at origin
    Epetra_SerialDenseMatrix& Z  ///< Z-operator
  );

  /// Total deformation gradient and (enhanced) Green-Lagrange strain
  /// \author mgit \date 01/08
  void w1_call_defgrad_tot(
    const Epetra_SerialDenseMatrix& F_enh,  ///< enhanced def.grad.
    Epetra_SerialDenseMatrix& F_tot,  ///< total def.grad.
    const Epetra_SerialDenseVector& F,  ///< displ.-based def.grad.
    Epetra_SerialDenseVector& strain  ///< GL strains
  );

  /// first Piola-Kirchhoff stress vector
  /// \author mgit \author 02/08
  void w1_stress_eas(
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    const Epetra_SerialDenseMatrix& F_tot,  ///< total def.grad.
    Epetra_SerialDenseMatrix& p_stress  ///< PK1 stress vector
  );

  /// calculate stiffness matrix kdd \f$\partial f_{int}/\partial d\f$
  /// \author mgit \date 03/08
  void w1_kdd(
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator at origin
    const Epetra_SerialDenseMatrix& F,  ///< total def.grad.
    const Epetra_SerialDenseMatrix& C,  ///< consititutive matrix
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress vector
    Epetra_SerialDenseMatrix& FCF,  ///< \f$F^T \dot C \dot F\f$
    Epetra_SerialDenseMatrix& estif,  ///< element stiff matrix
    const double fac  ///< integration factor
  );

  /// calculate tangential matrix kda \f$\partial f_{int}/\partial \alpha\f$
  /// \author mgit \date 03/08
  void w1_kda(
    const Epetra_SerialDenseMatrix& FCF,  ///< a product
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator at origin
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress
    const Epetra_SerialDenseMatrix& G,  ///< G-operator
    const Epetra_SerialDenseMatrix& Z,  ///< Z-operator
    Epetra_SerialDenseMatrix& Kda,  ///< target: kda
    const Epetra_SerialDenseMatrix& p_stress,  ///< PK1 stress
    const double fac  ///< integration factor
  );

  /// calculate tangential matrix kaa \f$\partial s/\partial \alpha\f$
  /// \author mgit \date 03/08
  void w1_kaa(
    const Epetra_SerialDenseMatrix& FCF,  ///< a product \f$F^T \dot C \dot F\f$
    const Epetra_SerialDenseMatrix& stress,  ///< PK2 stress
    const Epetra_SerialDenseMatrix& G,  ///< G-op
    Epetra_SerialDenseMatrix& Kaa,  ///< target: kaa
    const double fac  ///< integration factor
  );

  /// calculate internal forces fint(displacements u) and feas
  void w1_fint_eas(
    const Epetra_SerialDenseMatrix& W0,  ///< W-op at origin
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-op
    const Epetra_SerialDenseMatrix& G,  ///< G-op
    const Epetra_SerialDenseMatrix& p_stress,  ///< PK1 stress
    Epetra_SerialDenseVector& intforce,  ///< internal force vector
    Epetra_SerialDenseVector& feas,  ///< internal EAS constraint condition
    const double fac  ///< integration factor
  );

  //@}

  /// @name Material matters
  //@{

  /// Constitutive matrix \f$C\f$ and stresses
  /// \author mgit \date 05/07
  void w1_call_matgeononl(
    const Epetra_SerialDenseVector& strain,  ///< Green-Lagrange strain vector
    Epetra_SerialDenseMatrix& stress,  ///< stress matrix
    Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
    const int numeps,  ///< number of strains
    Teuchos::RCP<const MAT::Material> material,  ///< the material data
    Teuchos::ParameterList& params ///< element parameter list
  );

  /// Stress and constitutive matrix mapper from 3d to 2d
  /// due to dimensional reduction #wtype_
  ///
  /// \author bborn \date 06/09
  void MaterialResponse3dPlane(
    Epetra_SerialDenseMatrix& stress,  ///< stress matrix
    Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
    const Epetra_SerialDenseVector& strain,  ///< Green-Lagrange strain vector
    Teuchos::ParameterList& params ///< element parameter list
  );

  /// Generic 3D stress response
  ///
  /// \author bborn \date 06/09
  void MaterialResponse3d(
    LINALG::Matrix<6,1>* stress,  ///< 3D 2nd Piola-Kirchhoff stress vector
    LINALG::Matrix<6,6>* cmat,  ///< 3D elasticity matrix
    const LINALG::Matrix<6,1>* glstrain,  ///< 3D Green-Lagrange strain vector
    Teuchos::ParameterList& params ///< element parameter list
  );

  /// Map plane Green-Lagrange strains to 3d
  void GreenLagrangePlane3d(
    const Epetra_SerialDenseVector& glplane,  ///< 2d version of GL strain (Voigt notation)
    LINALG::Matrix<6,1>& gl3d ///< 3d version of GL strain (Voigt notation)
  );

  /// Internal/strain energy
  double EnergyInternal(
    Teuchos::RCP<const MAT::Material> material,  ///< element material
    Teuchos::ParameterList& params,///< element parameter list
    const Epetra_SerialDenseVector& Ev  ///< Green-Lagrange strain vector
  );

  /// Kinetic Energy
  /// \return kinetic energy of element
  double EnergyKinetic(
    const Epetra_SerialDenseMatrix& mass,  ///< element mass matrix
    const std::vector<double>& vel  ///< element velocity vector
  );

  //@}



  /// @name Tools for generalised energy-momentum method
  //@{

  /// evaluate the element forces and stiffness and mass for GEMM
  /// \author bborn \date 08/08
  void FintStiffMassGEMM(
    Teuchos::ParameterList& params,  ///<  a list containing GEMM coefficients
    const std::vector<int>& lm,  ///< location vector
    const std::vector<double>& dispo,  ///< element displacements \f$d_{n+1}^{(e)}\f$ at \f$t_{n}\f$
    const std::vector<double>& disp,  ///< element displacements \f$d_{n}^{(e)}\f$ at \f$t_{n+1}\f$
    const std::vector<double>& residual,  ///< residual displacements \f$\Delta d_{n+1}^{<k>(e)}\f$
    Epetra_SerialDenseMatrix* stiffmatrix,  ///< element stiffness matrix \f$k^{(e)}\f$
    Epetra_SerialDenseMatrix* massmatrix,  ///< element mass matrix \f$m^{(e)}\f$
    Epetra_SerialDenseVector* force,  ///< element internal mid-force vector \f$f_{int;m}^{(e)}\f$
    Epetra_SerialDenseMatrix* elestress,  ///< element stresses
    Epetra_SerialDenseMatrix* elestrain,  ///< element strains
    Teuchos::RCP<const MAT::Material> material,  ///< element material
    const INPAR::STR::StressType   iostress,   ///< stress output option
    const INPAR::STR::StrainType   iostrain    ///< strain output option
  );

  /// Tangent of internal force vector \f$f_{int;m}\f$ with respect to displacements \f$d_{n+1}\f$
  /// in pure displacement-based approach
  void TangFintByDispGEMM(
    const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
    const double& xigemm,  ///< GEMM's \f$\xi\f$ coefficient
    const double& fac,  ///< integration factor
    const Epetra_SerialDenseMatrix& bopm,  ///< non-lin. mid-B-op. \f$B_{m}\f$
    const Epetra_SerialDenseMatrix& bop,  ///< non-linear B-operator \f$B_{n+1}\f$ at \f$t_{n+1}\f$
    const Epetra_SerialDenseMatrix& C,  ///< constitutive matrix
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator
    const Epetra_SerialDenseMatrix& Smm,  ///< PK2 mid-stress matrix \f$\bar{S}_{m}\f$
    Epetra_SerialDenseMatrix& estif  ///< (in/out) element stiffness matrix
  );

  /// Tangent of internal force vector \f$f_{int;m}\f$ with respect to displacements \f$d_{n+1}\f$
  /// in EAS approach
  void TangFintByDispGEMM(
    const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
    const double& xigemm,  ///< GEMM's \f$\xi\f$ coefficient
    const double& fac,  ///< integration factor
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
    const Epetra_SerialDenseMatrix& W0m,  ///< mid-W-operator \f$W_{0;m}\f$
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator \f$W_{0;n+1}\f$
    const Epetra_SerialDenseMatrix& Fmm,  ///< mid def.grad. matrix \f$\bar{F}_m\f$
    const Epetra_SerialDenseMatrix& Fm,  ///< def.grad. matrix \f$\bar{F}_{n+1}\f$
    const Epetra_SerialDenseMatrix& C,  ///< constititive matrix \f$C\f$
    const Epetra_SerialDenseMatrix& Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
    Epetra_SerialDenseMatrix& FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
    Epetra_SerialDenseMatrix& estif  ///< (in/out) element stiffness matrix
  );

  /// Tangent of internal force vector \f$f_{int,m}\f$ with respect to enhanced strain parameters \f$\alpha_{n+1}\f$
  void TangFintByEnhGEMM(
    const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
    const double& xigemm,  ///< GEMM's \f$\xi\f$ coefficient
    const double& fac,  ///< integration factor
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
    const Epetra_SerialDenseMatrix& W0m,  ///< mid-W-operator \f$W_{0;m}\f$
    const Epetra_SerialDenseMatrix& FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
    const Epetra_SerialDenseMatrix& Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
    const Epetra_SerialDenseMatrix& G,  ///< G-operator \f$G_{n+1}\f$
    const Epetra_SerialDenseMatrix& Z,  ///< Z-operator \f$Z_{n+1}\f$
    const Epetra_SerialDenseMatrix& Pvmm,  ///< mid-mid 1st Piola-Kirchhoff stress vector \f$\tilde{P}_{mm} = \bar{F}_m \cdot S_m\f$
    Epetra_SerialDenseMatrix& kda  ///< (in/out) tangent \f$\partial f_{int;m}/\partial \alpha_{n+1}\f$
  );

  /// Tangent of enhancemed constraint vector \f$s_m\f$ with respect to displacements \f$d_{n+1}\f$
  void TangEconByDispGEMM(
    const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
    const double& xigemm,  ///< GEMM's \f$\xi\f$ coefficient
    const double& fac,  ///< integration factor
    const Epetra_SerialDenseMatrix& boplin,  ///< linear B-operator \f$B_L\f$
    const Epetra_SerialDenseMatrix& W0,  ///< W-operator \f$W_{0;n+1}\f$
    const Epetra_SerialDenseMatrix& FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
    const Epetra_SerialDenseMatrix& Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
    const Epetra_SerialDenseMatrix& Gm,  ///< mid-G-operator \f$G_{m}\f$
    const Epetra_SerialDenseMatrix& Z,  ///< Z-operator \f$Z_{n+1}\f$
    const Epetra_SerialDenseMatrix& Pvmm,  ///< mid-mid 1st Piola-Kirchhoff stress vector \f$\tilde{P}_{mm} = \bar{F}_m \cdot S_m\f$
    Epetra_SerialDenseMatrix& kad  ///< (in/out) tangent \f$\partial s_{m}/\partial d_{n+1}\f$
  );

  /// Tangent of enhancemed constraint vector \f$s_m\f$ with respect to enhanced strain parameters \f$\alpha_{n+1}\f$
  void TangEconByEnhGEMM(
    const double& alphafgemm,  ///< GEMM's \f$\alpha_f\f$ coefficient
    const double& xigemm,  ///< GEMM's \f$\xi\f$ coefficient
    const double& fac,  ///< integration factor
    const Epetra_SerialDenseMatrix& FmCF,  ///< the product \f$\bar{F}_m \cdot C \cdot \bar{F}_{n+1}^T\f$
    const Epetra_SerialDenseMatrix& Smm,  ///< 2nd Piola-Kirchhof mid-stress matrix \f$\bar{S}_m\f$
    const Epetra_SerialDenseMatrix& G,  ///< G-operator \f$G_{n+1}\f$
    const Epetra_SerialDenseMatrix& Gm,  ///< mid-G-operator \f$G_{m}\f$
    Epetra_SerialDenseMatrix& kaa  ///< (in/out) tangent \f$\partial s_{m}/\partial \alpha_{n+1}\f$
  );

  //@}

  /// don't want = operator
  Wall1& operator = (const Wall1& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getGaussrule
  (
    int* ngp  ///< number of Gauss points
  );
}; // class Wall1



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

class Wall1LineType : public DRT::ElementType
{
public:

  std::string Name() const { return "Wall1LineType"; }

  static Wall1LineType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Wall1LineType instance_;
};


/*======================================================================*/
/*!
 * \brief An element representing a line edge of a Wall1 element
 *
 * \note This is a pure Neumann boundary condition element. It's only
 *      purpose is to evaluate line Neumann boundary conditions that might be
 *      adjacent to a parent wall1 element. It therefore does not implement
 *      the DRT::Element::Evaluate method and does not have its own ElementRegister class.
 *
 * \author gitterle (gitterle@lnm.mw.tum.de)
 */
class Wall1Line : public DRT::FaceElement
{
public:

  /// @name Constructors and destructors and related methods
  //@{

  /// Standard Constructor
  Wall1Line(
    int id,  ///< A unique global id
    int owner,  ///< Processor owning this line
    int nnode,  ///< Number of nodes attached to this element
    const int* nodeids,  ///< global ids of nodes attached to this element
    DRT::Node** nodes,  ///<  the discretizations map of nodes to build ptrs to nodes from
    DRT::ELEMENTS::Wall1* parent, ///< The parent wall element of this line
    const int lline  ///< the local line number of this line w.r.t. the parent element
  );

  /// Copy Constructor
  ///
  /// Makes a deep copy of a Element
  Wall1Line(const Wall1Line& old);


  /// Deep copy this instance of an element and return pointer to the copy
  ///
  /// The Clone() method is used from the virtual base class Element in cases
  /// where the type of the derived class is unknown and a copy-ctor is needed
  DRT::Element* Clone() const;


  /// Get shape type of element
  virtual DiscretizationType Shape() const;


  /// Return unique ParObject id
  ///
  /// every class implementing ParObject needs a unique id defined at the
  /// top of the drt_parobject.H file.
  virtual int UniqueParObjectId() const { return Wall1LineType::Instance().UniqueParObjectId(); }

  /// Pack this class so it can be communicated
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Pack(DRT::PackBuffer& data) const;

  /// Unpack data from a char vector into this class
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Unpack(const std::vector<char>& data);

  /// Destructor
  virtual ~Wall1Line();

  //@}

  /// @name Access methods
  //@{

  /// Get number of degrees of freedom of a certain node
  /// (implements pure virtual DRT::Element)
  ///
  /// The element decides how many degrees of freedom its nodes must have.
  /// As this may vary along a simulation, the element can redecide the
  /// number of degrees of freedom per node along the way for each of it's nodes
  /// separately.
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return ParentMasterElement()->NumDofPerNode(node);}

  /// Get number of degrees of freedom per element
  /// (implements pure virtual DRT::Element)
  ///
  /// The element decides how many element degrees of freedom it has.
  /// It can redecide along the way of a simulation.
  ///
  /// \note Element degrees of freedom mentioned here are dofs that are visible
  ///       at the level of the total system of equations. Purely internal
  ///       element dofs that are condensed internally should NOT be considered.
  virtual int NumDofPerElement() const { return 0; }

  /// Print this element
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return Wall1LineType::Instance(); }

  //@}

  /// @name Evaluate methods
  //@{

  /// Evaluate a Neumann boundary condition
  ///
  /// this method evaluates a line Neumann condition on the wall element
  ///
  /// \return 0 if successful, negative otherwise
  virtual int EvaluateNeumann(
    Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication between control routine
                            ///<   and elements
    DRT::Discretization& discretization,  ///< (in) A reference to the underlying discretization
    DRT::Condition& condition,  ///< (in) The condition to be evaluated
    std::vector<int>& lm,  ///< (in) location vector of this element
    Epetra_SerialDenseVector& elevec1,  ///< (out) vector to be filled by element. If NULL on input
    Epetra_SerialDenseMatrix* elemat1 = NULL
  );

  /// Evaluate a boundary condition
  ///
  /// this method evaluates a line Neumann condition on the wall element
  ///
  /// \return 0 if successful, negative otherwise
  virtual int Evaluate(
    Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication between control routine
                            ///< and elements
    DRT::Discretization& discretization,  ///< (in) A reference to the underlying discretization
    std::vector<int>& lm,  ///< (in) location vector of this element
    Epetra_SerialDenseMatrix& elematrix1,
    Epetra_SerialDenseMatrix& elematrix2,
    Epetra_SerialDenseVector& elevector1,
    Epetra_SerialDenseVector& elevector2,
    Epetra_SerialDenseVector& elevector3
  );

  //! Evaluate method on mutliple dofsets for wall element
  virtual int Evaluate(Teuchos::ParameterList& params,
                       DRT::Discretization&      discretization,
                       LocationArray&            la,
                       Epetra_SerialDenseMatrix& elematrix1,
                       Epetra_SerialDenseMatrix& elematrix2,
                       Epetra_SerialDenseVector& elevector1,
                       Epetra_SerialDenseVector& elevector2,
                       Epetra_SerialDenseVector& elevector3);

  //@}

private:

  /// digestible actions
  enum ActionType
  {
    none,  ///< undetermined
    calc_struct_areaconstrstiff,
    calc_struct_centerdisp,
    calc_struct_constrarea,
    calc_potential_stiff,
    calc_struct_area_poro
  };

  /// don't want = operator
  Wall1Line& operator = (const Wall1Line& old);

  /// Submethod to compute necessary change to stiffness matrix due to the constraint
  void ComputeAreaConstrStiff(
    Epetra_SerialDenseMatrix xscurr,
    Epetra_SerialDenseMatrix& elematrix
  );

  /// Submethod to compute first derivatives of constraint area w.r.t. the displacements
  void ComputeAreaConstrDeriv(
    Epetra_SerialDenseMatrix xscurr,
    Epetra_SerialDenseVector& elevector
  );

  /// compute infintesimal line element dr for integration along the line
  double w1_substitution(
    const Epetra_SerialDenseMatrix& xyze,
    const Epetra_SerialDenseMatrix& deriv,
    std::vector<double>* unrm,
    const int iel
  );

  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(
    const DiscretizationType& distype
  );

}; // class Wall1Line





} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef WALL1_H
