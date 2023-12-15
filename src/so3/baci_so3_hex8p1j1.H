/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SO3_HEX8P1J1_H
#define BACI_SO3_HEX8P1J1_H

#include "baci_linalg_serialdensematrix.H"
#include "baci_so3_hex8.H"

BACI_NAMESPACE_OPEN

// forward declarations
struct _SOH8_DATA;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class So_Hex8P1J1Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "So_Hex8P1J1Type"; }

      static So_Hex8P1J1Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_Hex8P1J1Type instance_;

      std::string GetElementTypeString() const { return "SOLIDH8P1J1"; }
    };

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
      friend class So_Hex8P1J1Type;
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
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return So_Hex8P1J1Type::Instance().UniqueParObjectId();
      }

      /*!
     \brief Pack this class so it can be communicated

     \ref Pack and \ref Unpack are used to communicate this element

     */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      /*!
      \brief Calculate current from reference configuration moduli

      \param cmat (in): reference configuration moduli
      \param F (in): modified deformation gradient
      \param D_T_bar (out): current reference moduli
      \param t (in): current volume theta
      */
      void ConvertMat(const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat,
          const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& F,
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& D_T_bar, const double t);

      So_Hex8P1J1Type& ElementType() const override { return So_Hex8P1J1Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      void InitKpt()
      {
        // K_pt = N_p * N_t * detJ * w(gp)

        K_pt_ = 0.0;

        const static std::vector<double> gpweights = soh8_weights();

        for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
        {
          K_pt_ -= detJ_[gp] * gpweights[gp];
        }
      };



      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //! Compute internal force, its stiffness and mass matrix
      void ForceStiffMass(const std::vector<int>& lm,  ///< location matrix
          const std::vector<double>& disp,             ///< current displacements
          const std::vector<double>& residual,         ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>*
              stiffmatrix,                                             ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,      ///< element internal force vector
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,  // structure force
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
      );

      /// Return stress at Gauss point
      void Stress(CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>*
                      elestress,                  ///< store the stress herein
          const INPAR::STR::StressType iostress,  ///< stress type
          const int gp,                           ///< Gauss point index
          const double& detdefgrd,                ///< determinant of (assumed) deformation gradient
          const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>&
              defgrd,  ///< (assumed) deformation gradient
          const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress  ///< Cauchy stress vector
      );

      /// Return strain at Gauss point
      void Strain(CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>*
                      elestrain,                  ///< store the strain herein
          const INPAR::STR::StrainType iostrain,  ///< strain type to store for post-proc
          const int gp,                           ///< Gauss point index
          const double& detdefgrd,                ///< determinant of (assumed) deformation gradient
          const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,  ///< deformation gradient
          const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>&
              glstrain  ///< Green-Lagrange strain vector
      );

      /// Push-pull-operator
      static void PushPullOperator(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
                                       g,  ///< (G_ab^AB or G_AB^ab) or (G^ab_AB or G^AB_ab)
          const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>&
              f,                  ///< [F^-1]=[F^B_b] or [F]=[F^b_B]
          const bool& transpose,  ///< co-variant if true
          const double& fac       //  a scaling factor, eg det(F)
      );


     private:
      std::string GetElementTypeString() const { return "SOLIDH8P1J1"; }

      /** recover elementwise stored stuff */
      void soh8P1J1_recover(const std::vector<double>& residual);

      // don't want = operator
      So_Hex8P1J1& operator=(const So_Hex8P1J1& old);

      CORE::LINALG::Matrix<1, NUMDOF_SOH8> K_pu_;
      CORE::LINALG::Matrix<1, NUMDOF_SOH8> K_tu_;
      CORE::LINALG::Matrix<1, 1> R_t_;
      CORE::LINALG::Matrix<1, 1> R_p_;

      CORE::LINALG::Matrix<6, 1> m_;
      CORE::LINALG::Matrix<6, 6> Identity6_;
      CORE::LINALG::Matrix<6, 6> I_d_;
      CORE::LINALG::Matrix<6, 6> I_0_;


      /// @name Discontinuous primary field variables
      //@{
      CORE::LINALG::Matrix<1, 1> p_;    ///< pressure at current time/load step
      CORE::LINALG::Matrix<1, 1> p_o_;  ///< (old) pressure at last converged time/load step
      CORE::LINALG::Matrix<1, 1> dp_;   ///< pressure increment
      CORE::LINALG::Matrix<1, 1>
          t_;  ///< determinant of deformation gradient at current time/load step
      CORE::LINALG::Matrix<1, 1>
          t_o_;  ///< (old) Jacobian of deformation gradient at last converged time/load step
      CORE::LINALG::Matrix<1, 1> dt_;  ///< jacobian increment
      //@}

      double K_pt_{};
      double K_tt_{};
      CORE::LINALG::Matrix<1, 1> p_temp_;
      CORE::LINALG::Matrix<1, 1> t_temp_;

      CORE::LINALG::Matrix<24, 24> K_uu_;
      CORE::LINALG::Matrix<24, 1> F_u_;
    };


  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // SO3_HEX8P1J1_H
