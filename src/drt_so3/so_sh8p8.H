/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef D_SOLID3
#ifdef CCADISCRET

#ifndef SO_SH8P8_H
#define SO_SH8P8_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "so_sh8.H"
#include "../drt_inpar/inpar_structure.H"


// forward declarations
struct _SOH8_DATA;

namespace DRT
{

// forward declarations
class Discretization;

namespace ELEMENTS
{

// forward declarations
class Sosh8p8Register;
class So_sh8;

/// An incompressible 8-node solid shell element inherited from #DRT::ELEMENTS::so_sh8
/// utilising a Bochev-stabilised equal-order approach for 
/// tri-linearly Lagrangean interpolated displacement and pressure fields
///
/// <h3>References</h3>
/// The equal-order Bochev-stabilised approach for incompressible solid
/// is based on
/// - [1] C.R. Dohrmann, P.B. Bochev, "A stabilized finite element for the Stokes problem based on plynomial projections", Int. J. Numer. Fluids, 2000.
/// - [2] P.B. Bochev, C.R. Dohrmann, "Stabilization of low-order mixed finite elements for the Stokes equations".
/// - [3] C. Foerster, "Zur Bochev-Stabilisierung", internal note, 16.7.2008.
/// - [4] B. Bornemann, "Zum Bochev-stabilisierten Strukturelement fuer
///     geometrisch nicht-lineare Verformungen", internal note, 5.11.2008.
///
/// The solid shell/ANS  element technology is based on
/// - [5] M.A. Frenzel, "An advanced ...", PhD thesis, LNM, TU Muenchen, 2009.
/// - [6] T. Vu-Quoc, "Optimal solid shells for non-linear analyses
///     of multilayer composites", CMAME 2003
/// - [7] S. Klinkel, Gruttmann, W. Wagner, "A robust non-linear solid shell element
///     based on a mixed variational fromulation"
///
/// \author bborn
/// \date 03/09
class So_sh8p8 : public So_sh8
{
public:

  /// @name Friends
  //@{
  friend class Sosh8p8Register;
  friend class Soh8Surface;
  friend class Soh8Line;
  //@}

  /// @name Constructors and destructors and related methods
  //@{

  /// Standard Constructor
  So_sh8p8(
    int id,  ///<  A unique global ID
    int owner  ///< elements owning processor
    );

  /// Copy Constructor
  ///
  /// Makes a deep copy of a Element
  So_sh8p8(const So_sh8p8& old);

  /// Deep copy this instance of Solid3 and return pointer to the copy
  ///
  /// The Clone() method is used from the virtual base class Element in cases
  /// where the type of the derived class is unknown and a copy-ctor is needed
  DRT::Element* Clone() const;

  /// Return unique ParObject id
  ///
  /// every class implementing ParObject needs a unique id defined at the
  /// top of this file.
  virtual int UniqueParObjectId() const { return ParObject_So_sh8p8; }

  /// Pack this class so it can be communicated
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Pack(std::vector<char>& data) const;

  /// Unpack data from a char vector into this class
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Unpack(const std::vector<char>& data);

  /// Print this element
  void Print(std::ostream& os) const;

  /// Allocate and return a DesignElementRegister class
  ///
  /// Every implementation of DRT::Element must also implement
  /// an DRT::ElementRegister. This method allocates and returns
  /// an instance of a matching DRT::ElementRegister.
  virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const;

  //@}

  /// @name Input and Creation
  //@{

  /// Read input for this element
  bool ReadElement();

  /// Destructor
  virtual ~So_sh8p8();

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
  virtual int NumDofPerNode(const DRT::Node& node) const { return 4; }

  /// Get number of degrees of freedom per element
  /// (implements pure virtual DRT::Element)
  ///
  /// The element decides how many element degrees of freedom it has.
  /// It can redecide along the way of a simulation.
  ///
  /// \note Element degrees of freedom mentioned here are dofs that are visible
  /// at the level of the total system of equations. Purely internal
  /// element dofs that are condensed internally should NOT be considered.
  virtual int NumDofPerElement() const { return 0; }

  //@}

  /// @name Evaluation
  //@{

  /// Evaluate an element
  ///
  /// Evaluate so_sh8p8 element stiffness, mass, internal forces, etc.
  ///
  /// \return 0 if successful, negative otherwise
  int Evaluate(
    Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication between control routine
                                     ///< and elements
    DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
    std::vector<int>& lm,  ///< (in) location matrix for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (out) (stiffness-)matrix to be filled by element. If NULL on input,
                                        ///< the controling method does not expect the element to fill
                                        ///< this matrix.
    Epetra_SerialDenseMatrix& elemat2,  ///< (out) (mass-)matrix to be filled by element. If NULL on input,
                                        ///< the controling method does not expect the element to fill
                                        ///< this matrix.
    Epetra_SerialDenseVector& elevec1,  ///< (out) (internal force-)vector to be filled by element. If NULL on input,
                                        ///< the controlling method does not expect the element
                                        ///< to fill this vector
    Epetra_SerialDenseVector& elevec2,  ///< (out)  vector to be filled by element. If NULL on input,
                                        ///< the controlling method does not expect the element
                                        ///< to fill this vector
    Epetra_SerialDenseVector& elevec3  ///< (out) vector to be filled by element. If NULL on input,
                                       ///< the controlling method does not expect the element
                                       ///< to fill this vector
    );

  /// Kind of stabilisation for mixed, equal-order displacement-pressure appraoch
  enum StabilisationType
  {
    stab_affine,  ///< Bochev stabilisation for affine elements
    stab_nonaffine,  ///< Bochev stabilisation for non-affine elements
    stab_puredisp  ///< DEBUG ONLY: recover pure displacement-based approach
  };

  /// Kind of ANS anti-locking technique
  enum AnsType
  {
    ans_lateral,  ///< ANS active in t-/out-of-plane/thickness direction
    ans_none   ///< DEBUG ONLY: ANS switched off
  };

  /// Type of Voigt vector
  enum VoigtType
  {
    voigt6_strain,  ///< strain-like vector [ E11, E22, E33, 2*E12, 2*E23, 2*E31 ]
    voigt6_stress,  ///< stress-like vector [ S11, S22, S33, S12, S23, S31 ]
    voigt9  ///< non-symmetric vector [ F11, F22, F33, F12, F23, F31, F13, F32, F21 ]
  };

  // LinearizationType
  enum LinearizationType
  {
    lin_one,  ///< full linearisation
    lin_half,  ///< half-way linearisation
    lin_third ///< etc.
  };

private:

  // don't want = operator
  So_sh8p8& operator = (const So_sh8p8& old);

  /// @name Element constants
  //@{

  static const int NUMNOD_SOSH8P8 = NUMNOD_SOH8; //8  ///< number of nodes
  static const int NODDOF_SOSH8P8 = 4;  ///< number of dofs per node
  static const int NODDISP_SOSH8P8 = 3;  ///< number of displacements per node
  static const int NODPRES_SOSH8P8 = 1;  ///< number of pressures per node
  static const int NUMDOF_SOSH8P8 = 32;  ///< total dofs per element
  static const int NUMDISP_SOSH8P8 = 24;  ///< total discrete displacements per element
  static const int NUMPRES_SOSH8P8 = 8;  ///< total discrete pressures per element
  static const int NUMSTR_SOSH8P8 = NUMSTR_SOH8; //6  ///< number of strains/stresses
  static const int NUMDFGR_SOSH8P8 = 9;  ///< number of deformation gradient components
                                         ///< (deformation gradient is non-symmetric)
  static const int NUMGPT_SOSH8P8 = NUMGPT_SOH8; //8  ///< total gauss points per element
  static const int NUMDIM_SOSH8P8 = NUMDIM_SOH8; //3  ///< number of dimensions

  /// number of ANS sampling points, here 8
  static const int NUMSP_SOSH8P8 = 8;
  /// number of modified ANS strains (E_rt,E_st,E_tt), here 3
  static const int NUMANS_SOSH8P8 = 3;

  //@}

  /// kind of stabilisation, cf. #StabilisationType
  StabilisationType stab_;

  /// kind of ANS, cf. #AnsType
  AnsType ans_;

  /// kind of linearization
  LinearizationType lin_;

  /// @name Evaluate force and stiffness
  //@{

  /// Compute stiffness and mass matrix
  void ForceStiffMass(
    const std::vector<int>& lm,  ///< location matrix
    const LINALG::Matrix<NUMDISP_SOSH8P8,1>& disp, ///< current displacements
    const LINALG::Matrix<NUMPRES_SOSH8P8,1>& pres, ///< current pressures
    LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* massmatrix,  ///< element mass matrix
    LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* gradmatrix,  ///< element gradient matrix
    LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8>* stabmatrix,  ///< element stabilisation matrix
    LINALG::Matrix<NUMDISP_SOSH8P8,1>* force,  ///< element internal force vector
    LINALG::Matrix<NUMPRES_SOSH8P8,1>* incomp,  ///< incompressibility residual
    LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestrain,  ///< strains at GP
    double* volume,  ///< current element volume
    Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,   ///< stress output option
    const INPAR::STR::StrainType iostrain    ///< strain output option
    );

  /// Return stress at Gauss point
  void Stress(
    LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestress,  ///< store the stress herein
    const INPAR::STR::StressType iostress,  ///< stress type
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrd,  ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain,  ///< Green-Lagrange strain vector
    const LINALG::Matrix<NUMSTR_SOSH8P8,1>& stress,  ///< (deviatoric) 2nd Piola-Kirchhoff stress vector
    const double& pressure  ///< true pressure
    );

  /// Return strain at Gauss point
  void Strain(
    LINALG::Matrix<NUMGPT_SOSH8P8,NUMSTR_SOSH8P8>* elestrain,  ///< store the strain herein
    const INPAR::STR::StrainType iostrain,  ///< strain type
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrd,  ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invdefgrd,  ///< (assumed) inverted deformation gradient
    const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain  ///< Green-Lagrange strain vector
    );

  /// Recover deformation gradient incoperating assumed natural GL strain
  static void AssDefGrad(
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgrad,  ///< deformation gradient \f$[\boldsymbol{F}\f$
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& rgtstr,  ///< right stretch tensor \f$[\boldsymbol{U}\f$
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& defgradD,  ///< pure disp-based deformation gradient \f$[\boldsymbol{F}^d\f$
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& rgtstrD,  ///< pure disp-based right stretch tensor \f$[\boldsymbol{U}^d\f$
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& Jinv,  ///< inverse of transposed Jacobi matrix X_{,xi}
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& Jac,  ///< transposed Jacobi matrix X_{,xi}
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& jac,  ///< transposed Jacobi matrix x_{,xi}
    const LINALG::Matrix<NUMSTR_SOSH8P8,1>& glstrain  ///< material Green-Lagrange strain vector
                                                      ///< in global Cartesian components
    );

  /// Retrieve shear modulus
  ///
  /// Shear modulus is needed for stabilisation
  double ShearMod() const;

  //@}

  /// @name Voigt vector/matrix converters
  //@{

  /// Provide tensor indices of Voigt 6-vector notation
  ///
  /// These arrays are used to get from Voigt index to tensor index couple.
  static void Indices6VoigtTo2Tensor(
    const int*& voigt6row,  ///< tensor row indices
    const int*& voigt6col,  ///< tensor column indices
    const bool transpose = true  ///< deliver tranposed content if true
    );

  /// Provide tensor indices of Voigt 9-vector notation
  ///
  /// These arrays are used to get from Voigt index to tensor index couple.
  static void Indices9VoigtTo2Tensor(
    const int*& voigt9row,  ///< tensor row indices
    const int*& voigt9col,  ///< tensor column indices
    const bool transpose = false  ///< deliver tranposed content if true
    );

  /// Provide Voigt 9-indices 
  ///
  /// This array is used get from tensor couple to Voigt index
  /// access is via
  /// (i,j) -> 3*i+j
  static void Indices2TensorTo9Voigt(
    const int*& voigt3x3  ///< Voigt indices
    );

  /// Provide Voigt 6-indices 
  ///
  /// This array is used get from symmetric tensor couple to Voigt index
  /// access is via
  /// (i,j) -> 3*i+j
  static void Indices2TensorTo6Voigt(
    const int*& voigt3x3  ///< Voigt indices
    );

  /// Convert 3x3 matrix of non-sym 2-tensor in 9x1 Voigt vector
  ///
  /// From
  /// \f[
  ///   [\boldsymbol{F}] = \left[\begin{array}{ccc}
  ///   F_{11} & F_{12} & F_{13} \\ F_{21} & F_{22} & F_{23} \\ F_{31} & F_{32} & F_{33}
  ///   \end{array}\right]
  /// \f]
  /// to
  /// \f[
  ///   \tilde{\mathbf{F}} = \left[\begin{array}{c}
  ///   F_{11} \\ F_{22} \\ F_{33} \\ F_{12} \\ F_{23} \\ F_{31} \\ F_{13} \\ F_{32} \\ F_{21}
  ///   \end{array}\right]
  /// \f]
  static void Matrix2TensorToVector9Voigt(
    LINALG::Matrix<NUMDFGR_SOSH8P8,1>& fvct,  ///< (out) 9x1 Voigt vector
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat,  ///< (in) 3x3 matrix
    const bool transpose = false  ///< use transposed input 3x3 matrix
    );

  /// Convert 3x3 matrix of sym 2-tensor in 9-Voigt vector
  ///
  /// From
  /// \f[
  ///   [\boldsymbol{B}] = \left[\begin{array}{ccc}
  ///   B_{11} & B_{12} & B_{31} \\ B_{12} & B_{22} & B_{23} \\ B_{31} & B_{23} & B_{33}
  ///   \end{array}\right]
  /// \f]
  /// to if \i strain-like vector
  /// \f[
  ///   \tilde{\mathbf{B}} = \left[\begin{array}{c}
  ///   B_{11} \\ B_{22} \\ B_{33} \\ 2 B_{12} \\ 2 B_{23} \\ 2 B_{31}
  ///   \end{array}\right]
  /// \f]
  /// or to if \i stress-like vector
  /// \f[
  ///   \tilde{\mathbf{B}} = \left[\begin{array}{c}
  ///   B_{11} \\ B_{22} \\ B_{33} \\ B_{12} \\ B_{23} \\ B_{31}
  ///   \end{array}\right]
  /// \f]
  static void Matrix2TensorToVector6Voigt(
    LINALG::Matrix<NUMSTR_SOSH8P8,1>& bvct,  ///< 6-Voigt vector
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bmat,  ///< 3x3 sym 2-tensor
    const VoigtType outvoigt6 = voigt6_strain  ///< type of 6-Voigt vector
    );

  /// convert 6-Voigt vector into matrix of 2-tensor
  ///
  /// opposite of #Matrix2TensorToVector6Voigt
  static void Vector6VoigtToMatrix2Tensor(
    LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bmat,  ///< 3x3 sym 2-tensor
    const LINALG::Matrix<NUMSTR_SOSH8P8,1>& bvct,  ///< 6-Voigt vector
    const VoigtType invoigt6 = voigt6_strain   ///< type of 6-Voigt vector
    );

  /// Derivative of inverse of non-symmetric 2-tensor with respect to itself
  ///
  /// In index notation:
  /// \f[
  ///   \frac{\partial (F^{-1})_{ij}}{\partial F_{kl}}
  ///   = -(F^{-1})_{ik} \cdot (F^{-1})_{lj}
  /// \f]
  static void InvVector9VoigtDiffByItself(
    LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8>& invfderf,
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat,
    const bool transpose = false  ///< use transposed input 3x3 matrix
    );

  /// Derivative of inverse of symmetric 2-tensor with respect to itself
  ///
  /// In index notation:
  /// \f[
  ///   \frac{\partial (F^{-1})_{ij}}{\partial F_{kl}}
  ///   = -1/2*\big( (F^{-1})_{ik} \cdot (F^{-1})_{lj} 
  ///                + (F^{-1})_{il} \cdot (F^{-1})_{kj} \big)
  /// \f]
  static void InvVector6VoigtDiffByItself(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& invfderf,  ///< 6x6 derivative
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat  ///< 3x3 inverse symmetric 2-tensor
    );

  /// 2nd derivative of inverse of symmetric 2-tensor with respect to itself
  static void InvVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8*NUMSTR_SOSH8P8>& invbvdderb,  ///< 6x36 matrix
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt  ///< 3x3 inverse symmetric 2-tensor
    );

  /// Derivative of square of symmetric 2-tensor by itself
  ///
  /// In index notation:
  /// \f[
  ///   \frac{\partial F_{im} F_{mj}}{\partial F_{kl}}
  ///   = \delta_{ik}\cdot F_{il} + \delta_{jl}\cdot F_{ik}
  /// \f]
  static void SqVector6VoigtDiffByItself(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& sqfderf,  ///< diff. squared matrix 
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat  ///< symmetric 2-tensor
    );

  /// Derivative of "square" of non-symmetric 2-tensor by itself
  ///
  /// cf. #SqVector6VoigtDiffByItself
  ///
  /// In compact tensor notation
  /// \f[
  ///   \big( \boldsymbol{F}^T \cdot \boldsymbol{F} \big)_{,\boldsymbol{F}}
  /// \f]
  static void SqVector9VoigtDiffByItself(
    LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8>& sqfderf,  ///< diff. squared matrix 
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat,  ///< non-symmetric 2-tensor
    const bool transpose = false  ///< transpose non-symmetric 2-tensor
    );

  /// 2nd derivative of square of symmetric 2-tensor with respect to itself
  static void SqVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8*NUMSTR_SOSH8P8>& sqfdderf,  ///< 2nd derivative
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat  ///< symmetric 2-tensor (not needed)
    );

  /// Build 6x9 Voigt matrix for 2-tensor-dot-2-tensor product
  ///
  /// This method builds a matrix $\mathbf{B}$ which allows to perform the 
  /// following operation utilising a Voigt matrix.
  ///
  /// The operation in
  /// (1) compact tensor notation
  /// \f[
  ///   \boldsymbol{A} = \boldsymbol{A}^T = \boldsymbol{B}^T \cdot \boldsymbol{C}
  /// \f]
  /// or (2) in index tensor notation
  /// \f[
  ///   A_{ij} = A_{ji} = B_{ki} \cdot C_{kj}
  /// \f]
  /// or (3) in Voigt 9-vector by 6x9-matrix product, i.e.
  /// \f[
  ///   \mathbf{A} = \mathbf{B} \; \tilde{\mathbf{C}}
  /// \f]
  static void Matrix2TensorToMatrix6x9Voigt(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMDFGR_SOSH8P8>& bm,  ///< (out) 6x9 Voigt matrix
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt,  ///< (in) 3x3 matrix
    const bool transpose = true  ///< transpose entries
    );

  /// Build 6x6-Voigt matrix to multiply 6-Voigt with which is equivalent to
  /// multiplication with non-sym 2-tensor \f$\boldsymbol{F}\f$ from left and right 
  /// of a symmetric 2-tensor \f$\boldsymbol{e}\f$
  ///
  /// The operation in
  /// (1) compact tensor notation
  /// \f[
  ///    \boldsymbol{E} = \boldsymbol{F}^T \cdot \boldsymbol{e} \cdot \boldsymbol{F}
  /// \f]
  /// or (2) in index tensor notation
  /// \f[
  ///    E_AB = F_A^a  \cdot e_ab  \cdot F^b_B
  /// \f]
  /// or  (3) in Voigt 6-vector by 6x6-matrix product, i.e.
  /// \f[
  ///    \mathbf{E} = \mathbf{B} \; \mathbf{e}
  /// \f]
  ///
  /// This operation is also known as 'pull back'
  static void Matrix2TensorToLeftRightProductMatrix6x6Voigt(
    LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& bm,  ///< (out) 6x6 Voigt matrix
    const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt,  ///< (in) 3x3 matrix of 2-tensor
    const bool transpose = false, ///< 3x3 input matrix is transposed
    const VoigtType outvoigt6 = voigt6_strain,  ///< 6-Voigt vector layout on rows of 6x6 matrix
    const VoigtType invoigt6 = voigt6_strain  ///< 6-Voigt vector layout on columns of 6x6 matrix
    );

  //@}

  /// @name Local assemble of displacement and pressure quantities
  //@{

  /// Extract 24x1 displacement and 8x1 pressure 
  /// of 32x1 displacement-pressure vector
  static void ExtractDispAndPres(
    std::vector<double>& mystat,  ///< 32x1 element displacement-pressure vector
    LINALG::Matrix<NUMDISP_SOSH8P8,1>& mydisp,  ///< 24x1 element displacement vector
    LINALG::Matrix<NUMPRES_SOSH8P8,1>& mypres  ///< 8x1 element pressure vector
    );

  /// Build 32x32 element matrix
  static void BuildElementMatrix(
    LINALG::Matrix<NUMDOF_SOSH8P8,NUMDOF_SOSH8P8>* mat,  ///< 32x32 matrix
    const LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* matdd,  ///< 24x24 sub-matrix, if NULL then 0s are inserted
    const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matdp,  ///< 24x8 sub-matrix, if NULL then 0s are inserted
    const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matpd,  ///< 8x24 sub-matrix, if NULL then transpose of matdp is inserted or 0s
    const LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8>* matpp  ///< 8x8 sub-matrix, if NULL then 0s are inserted
    );

  /// Build 32x1 element vector
  static void BuildElementVector(
    LINALG::Matrix<NUMDOF_SOSH8P8,1>* vct,  ///< 32x1 vector
    const LINALG::Matrix<NUMDISP_SOSH8P8,1>* vctd,  ///< 24x1 sub-vector, if NULL then 0s are inserted
    const LINALG::Matrix<NUMPRES_SOSH8P8,1>* vctp  ///< 8x1 sub-vector, if NULL then 0s are inserted
    );

  /// Assemble global volume
  static void AssembleVolume(
    Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
    const double& volume  ///< current element volume
    );

  //@}

}; // class So_sh8p8


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/// A register for incompressible Solid-Shell8 element, inherited from Sosh8
///
/// \author bborn
/// \date 03/09
class Sosh8p8Register : public Sosh8Register
{
public:

  /// @name Constructors and destructors and related methods
  //@{

  /// Standard Constructor
  Sosh8p8Register(DRT::Element::ElementType etype);

  /// Copy Constructor
  ///
  /// Makes a deep copy of an Soh8Register
  Sosh8p8Register(const DRT::ELEMENTS::Sosh8p8Register& old);

  /// Deep copy this instance of Sosh8p8Register and return pointer to the copy
  ///
  /// The Clone() method is used from the virtual base class Element in cases
  /// where the type of the derived class is unknown and a copy-ctor is needed
  DRT::ELEMENTS::Sosh8p8Register* Clone() const;

  /// Return unique ParObject id
  ///
  /// every class imploementing ParObject needs a unique id defined at the
  /// top of this file.
  virtual int UniqueParObjectId() const { return ParObject_Sosh8p8Register;}

  /// Pack this class so it can be communicated
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Pack(std::vector<char>& data) const;

  /// Unpack data from a char vector into this class
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Unpack(const std::vector<char>& data);

  /// Destructor
  virtual ~Sosh8p8Register();

  //@}

  /// @name Acess methods
  //@{

  /// Print this element
  void Print(std::ostream& os) const;

  /// Call Element init routines
  ///
  /// The method performs all steps needed to init the solid-shell8 elements.
  /// This includes defining the thickness direction of the solid-shell by looking
  /// at its geometric shape in reference configuration. Thus, the loop looks a bit
  /// complicated because a "not-thin" So_sh8p8 is morphed to a So_hex8 where the
  /// discretization is not filled-complete anymore.
  int Initialize(DRT::Discretization& dis);

  //@}

}; // class Sosh8p8Register


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef SO_SH8P8_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
