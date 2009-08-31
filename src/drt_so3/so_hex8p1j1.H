/*----------------------------------------------------------------------*/
/*!
\file so_q1p0hex8.H

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_HEX8P1J1_H
#define SO_HEX8P1J1_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "so_hex8.H"
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
class SoHex8P1J1Register;

/// The so-called 'Q1P0' element in 8-node hexahedron shape
///
/// <h3>About</h3>
///   The element is a mixed method based on a three-field principle.
/// Firstly, the displacement field, which is C^0-discretised with
/// tri-linear Lagrangean polynomials. Secondly, the pressure #p_,
/// which is discretised discontinuously across element boundaries
/// in a constant manner, ie one pressure DOF per element. Thirdly,
/// the determinant of the deformation gradient (Jacobian) #t_, which is as
/// well discretised discontinuously across element boundaries
/// in a constant manner, ie one Jacobian DOF per element.
///   The approach tackles volumetric locking. However, it does
/// not prevent shear locking etc.
///
/// <h3>References</h3>
/// - [1] OC Zienkiewicz, RL Taylor, The Finite Element Method for Solid
///       and Structural Mechanics, Butterworth Heinemann, 6th edition,
///       2005. Especially Section 5.5.
///
/// \author lw
/// \date spring/09
class So_Hex8P1J1 : public So_hex8
{
public:

  //! @name Friends
  friend class SoHex8P1J1Register;
  friend class Soh8Surface;
  friend class Soh8Line;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owning processor
  */
  So_Hex8P1J1(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_Hex8P1J1(const So_Hex8P1J1& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_So_Hex8P1J1;}

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
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Calculate current from reference configuration moduli

  \param cmat (in): reference configuration moduli
  \param F (in): modified deformation gradient
  \param D_T_bar (out): current reference moduli
  \param t (in): current volume theta
  */
  void ConvertMat(const LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8>& cmat,
                  const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& F,
                  LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8>& D_T_bar,
                  const double t);

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

  void InitKpt()
  {
    // K_pt = N_p * N_t * detJ * w(gp)

    K_pt_ = 0.0;

    const static vector<double> gpweights = soh8_weights();

    for (int gp=0; gp<NUMGPT_SOH8; ++gp)
    {
      K_pt_ -= detJ_[gp] * gpweights[gp];
    }
  };

  /*!
  \brief Destructor

  */
  virtual ~So_Hex8P1J1();


  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate element stiffness, mass, internal forces, etc.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param discretization : pointer to discretization for de-assembly
  \param lm (in)        : location matrix for de-assembly
  \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                          the controlling method does not expect the element to fill
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
  int Evaluate(ParameterList&            params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //! Compute internal force, its stiffness and mass matrix
  void ForceStiffMass(
    const std::vector<int>& lm,  ///< location matrix
    const std::vector<double>& disp,  ///< current displacements
    const std::vector<double>& residual,  ///< current residual displ
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,    ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,     ///< element mass matrix
    LINALG::Matrix<NUMDOF_SOH8,1>* force,  ///< element internal force vector
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,   ///< stress output option
    const INPAR::STR::StrainType iostrain   ///< strain output option
    );

  /// Return stress at Gauss point
  void Stress(
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  ///< store the stress herein
    const INPAR::STR::StressType iostress,  ///< stress type
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& defgrd,  ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMSTR_SOH8,1>& stress  ///< Cauchy stress vector
    );

  /// Return strain at Gauss point
  void Strain(
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  ///< store the strain herein
    const INPAR::STR::StrainType iostrain,  ///< strain type to store for post-proc
    const int gp,  ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& defgrd,  ///< deformation gradient
    const LINALG::Matrix<NUMSTR_SOH8,1>& glstrain  ///< Green-Lagrange strain vector
    );

  /// Push-pull-operator
  static void PushPullOperator(
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8>& g, ///< (G_ab^AB or G_AB^ab) or (G^ab_AB or G^AB_ab)
    const LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& f,  ///< [F^-1]=[F^B_b] or [F]=[F^b_B]
    const bool& transpose,  ///< co-variant if true
    const double& fac  //  a scaling factor, eg det(F)
    );

#if 0
  void test_stiffmat(
    const std::vector<int>& lm,  ///< location matrix
    const std::vector<double>& disp,  ///< current displacements
    const std::vector<double>& residual,   ///< current residual displ
    ParameterList& params  ///< algorithmic parameters e.g. time
    );
#endif


private:

  // don't want = operator
  So_Hex8P1J1& operator = (const So_Hex8P1J1& old);

  LINALG::Matrix<1,NUMDOF_SOH8> K_pu_;
  LINALG::Matrix<1,NUMDOF_SOH8> K_tu_;
  LINALG::Matrix<1,1> R_t_;
  LINALG::Matrix<1,1> R_p_;

  LINALG::Matrix<6,1> m_;
  LINALG::Matrix<6,6> Identity6_;
  LINALG::Matrix<6,6> I_d_;
  LINALG::Matrix<6,6> I_0_;


  /// @name Discontinuous primary field variables
  //@{
  LINALG::Matrix<1,1> p_;  ///< pressure at current time/load step
  LINALG::Matrix<1,1> p_o_;  ///< (old) pressure at last converged time/load step
  LINALG::Matrix<1,1> t_;  ///< determinant of deformation gradient at current time/load step
  LINALG::Matrix<1,1> t_o_;  ///< (old) Jacobian of deformation gradient at last converged time/load step
  //@}

  double K_pt_;
  double K_tt_;
  LINALG::Matrix<1,1> p_temp_;
  LINALG::Matrix<1,1> t_temp_;

  LINALG::Matrix<24,24> K_uu_;
  LINALG::Matrix<24,1> F_u_;

};


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


class SoHex8P1J1Register : public Soh8Register
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  SoHex8P1J1Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an SoHex8P1J1Register

  */
  SoHex8P1J1Register(const DRT::ELEMENTS::SoHex8P1J1Register& old);

  /*!
  \brief Deep copy this instance of SoHex8P1J1Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::SoHex8P1J1Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_SoHex8P1J1Register;}

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
  virtual ~SoHex8P1J1Register();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Call Element init routines
  */
  int Initialize(DRT::Discretization& dis);


}; // class SoHex8P1J1Register

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SO_HEX8P1J1_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

