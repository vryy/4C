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
  
private:

  /// @name Properties
  //@{
  
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
    lin_sixth = 0,  ///< minimum linearisation to achieve anti-locking
    lin_half,  ///< half-way linearisation
    lin_one  ///< full linearisation
  };
  
  //@}

  /// @name Element constants
  //@{

  static const int NUMNOD_ = NUMNOD_SOH8; //8  ///< number of nodes
  static const int NODDOF_ = 4;  ///< number of dofs per node
  static const int NODDISP_ = 3;  ///< number of displacements per node
  static const int NODPRES_ = 1;  ///< number of pressures per node
  static const int NUMDOF_ = 32;  ///< total dofs per element
  static const int NUMDISP_ = 24;  ///< total discrete displacements per element
  static const int NUMDISPSQSYM_ = 300;  ///< (NUMDISP_*NUMDISP_ + NUMDISP_)/2, number of relevant entries if symmetric in element displacements
  static const int NUMPRES_ = 8;  ///< total discrete pressures per element
  static const int NUMSTR_ = NUMSTR_SOH8; //6  ///< number of strains/stresses
  static const int NUMDFGR_ = 9;  ///< number of deformation gradient components
                                         ///< (deformation gradient is non-symmetric)
  static const int NUMGPT_ = NUMGPT_SOH8; //8  ///< total gauss points per element
  static const int NUMDIM_ = NUMDIM_SOH8; //3  ///< number of dimensions

  /// number of ANS sampling points, here 8
  static const int NUMSP_ = 8;
  /// number of modified ANS strains (E_rt,E_st,E_tt), here 3
  static const int NUMANS_ = 3;

  //@}

public:

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
  
  /// Set ANS type
  void SetANS(
     const AnsType& newans  ///< ANS to set
     ) 
   { 
      ans_ = newans;
   }

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



private:

  /// don't want = operator
  So_sh8p8& operator = (const So_sh8p8& old);

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
    const LINALG::Matrix<NUMDISP_,1>& disp, ///< current displacements
    const LINALG::Matrix<NUMPRES_,1>& pres, ///< current pressures
    LINALG::Matrix<NUMDISP_,NUMDISP_>* massmatrix,  ///< element mass matrix
    LINALG::Matrix<NUMDISP_,NUMDISP_>* stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<NUMDISP_,NUMPRES_>* gradmatrix,  ///< element gradient matrix
    LINALG::Matrix<NUMPRES_,NUMPRES_>* stabmatrix,  ///< element stabilisation matrix
    LINALG::Matrix<NUMDISP_,1>* force,  ///< element internal force vector
    LINALG::Matrix<NUMPRES_,1>* incomp,  ///< incompressibility residual
    LINALG::Matrix<NUMGPT_,NUMSTR_>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_,NUMSTR_>* elestrain,  ///< strains at GP
    double* volume,  ///< current element volume
    Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,   ///< stress output option
    const INPAR::STR::StrainType iostrain    ///< strain output option
    );

  /// Return stress at Gauss point
  void Stress(
    LINALG::Matrix<NUMGPT_,NUMSTR_>* elestress,  ///< store the stress herein
    const INPAR::STR::StressType iostress,  ///< stress type
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrd,  ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMSTR_,1>& glstrain,  ///< Green-Lagrange strain vector
    const LINALG::Matrix<NUMSTR_,1>& stress,  ///< (deviatoric) 2nd Piola-Kirchhoff stress vector
    const double& pressure  ///< true pressure
    );

  /// Return strain at Gauss point
  void Strain(
    LINALG::Matrix<NUMGPT_,NUMSTR_>* elestrain,  ///< store the strain herein
    const INPAR::STR::StrainType iostrain,  ///< strain type
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrd,  ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& invdefgrd,  ///< (assumed) inverted deformation gradient
    const LINALG::Matrix<NUMSTR_,1>& glstrain  ///< Green-Lagrange strain vector
    );

  /// Recover deformation gradient incoperating assumed natural GL strain
  static void AssDefGrad(
    double& detdefgrad,  ///< determinat of deformation gradient
    LINALG::Matrix<NUMDIM_,NUMDIM_>& defgrad,  ///< deformation gradient \f$[\boldsymbol{F}\f$
    LINALG::Matrix<NUMDIM_,NUMDIM_>& invdefgrad,  ///< inverse deformation gradient \f$[\boldsymbol{F}^{-1}\f$
    LINALG::Matrix<NUMDIM_,NUMDIM_>& rgtstr,  ///< right stretch tensor \f$[\boldsymbol{U}\f$
    LINALG::Matrix<NUMDIM_,NUMDIM_>& defgradD,  ///< pure disp-based deformation gradient \f$[\boldsymbol{F}^d\f$
    LINALG::Matrix<NUMDIM_,NUMDIM_>& rgtstrD,  ///< pure disp-based right stretch tensor \f$[\boldsymbol{U}^d\f$
    LINALG::Matrix<NUMDIM_,NUMDIM_>& invrgtstrD,  ///< inverted pure disp-based right stretch tensor \f$[\boldsymbol{U}^{d-1}\f$
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& Jinv,  ///< inverse of transposed Jacobi matrix X_{,xi}
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& Jac,  ///< transposed Jacobi matrix X_{,xi}
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& jac,  ///< transposed Jacobi matrix x_{,xi}
    const LINALG::Matrix<NUMSTR_,1>& glstrain  ///< material Green-Lagrange strain vector
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
    LINALG::Matrix<NUMDFGR_,1>& fvct,  ///< (out) 9x1 Voigt vector
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,  ///< (in) 3x3 matrix
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
    LINALG::Matrix<NUMSTR_,1>& bvct,  ///< 6-Voigt vector
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& bmat,  ///< 3x3 sym 2-tensor
    const VoigtType outvoigt6 = voigt6_strain  ///< type of 6-Voigt vector
    );

  /// convert 6-Voigt vector into matrix of 2-tensor
  ///
  /// opposite of #Matrix2TensorToVector6Voigt
  static void Vector6VoigtToMatrix2Tensor(
    LINALG::Matrix<NUMDIM_,NUMDIM_>& bmat,  ///< 3x3 sym 2-tensor
    const LINALG::Matrix<NUMSTR_,1>& bvct,  ///< 6-Voigt vector
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
    LINALG::Matrix<NUMDFGR_,NUMDFGR_>& invfderf,
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& invfmat,
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
    LINALG::Matrix<NUMSTR_,NUMSTR_>& invfderf,  ///< 6x6 derivative
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& invfmat  ///< 3x3 inverse symmetric 2-tensor
    );

  /// 2nd derivative of inverse of symmetric 2-tensor with respect to itself
  static void InvVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_>& invbvdderb,  ///< 6x36 matrix
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& ibt  ///< 3x3 inverted symmetric 2-tensor
    );

  /// Derivative of square of symmetric 2-tensor by itself
  ///
  /// In index notation:
  /// \f[
  ///   \frac{\partial F_{im} F_{mj}}{\partial F_{kl}}
  ///   = \delta_{ik}\cdot F_{il} + \delta_{jl}\cdot F_{ik}
  /// \f]
  static void SqVector6VoigtDiffByItself(
    LINALG::Matrix<NUMSTR_,NUMSTR_>& sqfderf,  ///< diff. squared matrix 
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,  ///< symmetric 2-tensor
    const VoigtType outvoigt6 = voigt6_strain  ///< 6-Voigt row type of resulting matrix
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
    LINALG::Matrix<NUMDFGR_,NUMDFGR_>& sqfderf,  ///< diff. squared matrix 
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,  ///< non-symmetric 2-tensor
    const bool transpose = false  ///< transpose non-symmetric 2-tensor
    );

  /// 2nd derivative of square of symmetric 2-tensor with respect to itself
  static void SqVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_>& sqfdderf,  ///< 2nd derivative
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat  ///< symmetric 2-tensor (not needed)
    );

  /// 2nd derivative of square of symmetric 2-tensor with respect to itself
  ///
  /// cf. #SqVector6VoigtTwiceDiffByItself bit with sparse storage of 6-tensor
  static void SqVector6VoigtTwiceDiffByItself(
    int* isqfdderf,  ///< indices
    LINALG::Matrix<NUMSTR_,6>& sqfdderf  ///< 2nd derivative data (sparsely stored)
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
    LINALG::Matrix<NUMSTR_,NUMDFGR_>& bm,  ///< (out) 6x9 Voigt matrix
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& bt,  ///< (in) 3x3 matrix
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
    LINALG::Matrix<NUMSTR_,NUMSTR_>& bm,  ///< (out) 6x6 Voigt matrix
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& bt,  ///< (in) 3x3 matrix of 2-tensor
    const bool transpose = false, ///< 3x3 input matrix is transposed
    const VoigtType outvoigt6 = voigt6_strain,  ///< 6-Voigt vector layout on rows of 6x6 matrix
    const VoigtType invoigt6 = voigt6_strain  ///< 6-Voigt vector layout on columns of 6x6 matrix
    );

  //@}

  /// determine (inverse) right stretch tensor
  ///
  /// The right stretch tensor is obtained by polar decomposition 
  /// of the right Cauchy-Green 2-tensor in 3 dimensions.
  ///
  /// About:
  ///   Polar decomposition is often applied to the material deformation
  /// tensor F, i.e.
  ///      F = R . U = v . R
  /// in which R is the rotation matrix (two-point tensor), U the material
  /// stretch tensor (refers to the undeformed configuration) and v the
  /// spatial stretch tensor (refers to the deformed configuration).
  ///   This polar decomposition is also applied to the isoparametric
  /// Jacobian tensor (mapping quantities in parameter space to material
  /// configuration). However, the local variables are denoted for the
  /// case of a deformation tensor.
  ///
  /// References:
  /// [1] A. Hoger & D.E. Carlson, "Determination of the stretch and
  ///        rotation in the polar decomposition of the deformation
  ///        gradient", Quart. Appl. Math., 42(2):113-117, 1984.
  /// [2] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  ///        esp. Section 2.6
  static void StretchTensor(
    double* detut,  ///< determinant of material stretch tensor
    LINALG::Matrix<NUMDIM_,NUMDIM_>* ut,  ///< material stretch tensor
    LINALG::Matrix<NUMDIM_,NUMDIM_>* invut,  ///< inverse material stretch tensor
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& ct  ///< right Cauchy-Green tensor
    );

  /// Spectral decomposition of symmetric 2-tensor calculated
  /// iteratively using Jacobi's method
  ///
  /// References:
  /// [1] I.N. Bronstein & K.A. Semendjajew, "Taschenbuch der
  ///     Mathematik", Teubner, 25.ed, 1991.
  ///     esp. p. 741 ff.
  ///
  /// \return error: 0=success, 1=failure
  static int SymSpectralDecompJacIter(
    LINALG::Matrix<NUMDIM_,NUMDIM_>& ew,  ///< diagional matrix of eigenvalues
    LINALG::Matrix<NUMDIM_,NUMDIM_>& ev,  ///< orthornormal eigenvectors (direction cosines)
    const LINALG::Matrix<NUMDIM_,NUMDIM_>& at,  ///< <i>symmetric</i> input matrix
    const double itertol,  ///< tolerance
    const int itermax  ///< maximally allowd iteration steps
    );

  /// @name Local assemble of displacement and pressure quantities
  //@{

  /// Extract 24x1 displacement and 8x1 pressure 
  /// of 32x1 displacement-pressure vector
  static void ExtractDispAndPres(
    std::vector<double>& mystat,  ///< 32x1 element displacement-pressure vector
    LINALG::Matrix<NUMDISP_,1>& mydisp,  ///< 24x1 element displacement vector
    LINALG::Matrix<NUMPRES_,1>& mypres  ///< 8x1 element pressure vector
    );

  /// Build 32x32 element matrix
  static void BuildElementMatrix(
    LINALG::Matrix<NUMDOF_,NUMDOF_>* mat,  ///< 32x32 matrix
    const LINALG::Matrix<NUMDISP_,NUMDISP_>* matdd,  ///< 24x24 sub-matrix, if NULL then 0s are inserted
    const LINALG::Matrix<NUMDISP_,NUMPRES_>* matdp,  ///< 24x8 sub-matrix, if NULL then 0s are inserted
    const LINALG::Matrix<NUMDISP_,NUMPRES_>* matpd,  ///< 8x24 sub-matrix, if NULL then transpose of matdp is inserted or 0s
    const LINALG::Matrix<NUMPRES_,NUMPRES_>* matpp  ///< 8x8 sub-matrix, if NULL then 0s are inserted
    );

  /// Build 32x1 element vector
  static void BuildElementVector(
    LINALG::Matrix<NUMDOF_,1>* vct,  ///< 32x1 vector
    const LINALG::Matrix<NUMDISP_,1>* vctd,  ///< 24x1 sub-vector, if NULL then 0s are inserted
    const LINALG::Matrix<NUMPRES_,1>* vctp  ///< 8x1 sub-vector, if NULL then 0s are inserted
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
