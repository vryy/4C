/*----------------------------------------------------------------------*/
/*! \file
\brief a in-plane bi-quadratic solid-like 7-parameter shell
\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_SH18_HPP
#define FOUR_C_SO3_SH18_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_mat_material.hpp"
#include "baci_so3_hex18.hpp"

FOUR_C_NAMESPACE_OPEN

const int num_eas = 9;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class SoSh18Type : public SoHex18Type
    {
     public:
      //! @name Friends
      friend class SoSh18PlastType;

      std::string Name() const override { return "SoSh18Type"; }

      static SoSh18Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoSh18Type instance_;

      std::string GetElementTypeString() const { return "SOLIDSH18"; }
    };

    /*!
    \brief A C++ version of the 18-node hex solid-shell element

    A structural 18-node hexahedral solid element for large deformations.
    Quadratic interpolation within a plane and linear interpolation in the third direction
    */
    class SoSh18 : public virtual SoHex18
    {
     public:
      //! @name Friends
      friend class SoSh18Type;
      friend class SoSh18PlastType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoSh18(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoSh18(const SoSh18& old);

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
      int UniqueParObjectId() const override { return SoSh18Type::Instance().UniqueParObjectId(); }

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


      //@}

      //! @name Acess methods

      /*!
       \brief Does this element use EAS?
       */
      bool HaveEAS() const override { return eas_; };

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;



      DRT::ElementType& ElementType() const override { return SoSh18Type::Instance(); }

      //@}

      //! @name Input and Creation

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      //@}


     protected:
      //! use DSG for transverse shear locking
      bool dsg_shear_{};
      //! use DSG for membrane locking
      bool dsg_membrane_{};
      //! use DSG for curvature thickness locking
      bool dsg_ctl_{};
      //! use EAS
      bool eas_{};

      // EAS stuff *****************************************************
      CORE::LINALG::Matrix<num_eas, num_eas> KaaInv_;
      CORE::LINALG::Matrix<num_eas, 1> feas_;
      CORE::LINALG::Matrix<num_eas, NUMDOF_SOH18> Kad_;
      CORE::LINALG::Matrix<num_eas, 1> alpha_eas_;
      CORE::LINALG::Matrix<num_eas, 1> alpha_eas_last_timestep_;
      CORE::LINALG::Matrix<num_eas, 1> alpha_eas_delta_over_last_timestep_;
      CORE::LINALG::Matrix<num_eas, 1> alpha_eas_inc_;
      double old_step_length_{};
      // EAS stuff *****************************************************

      // DSG factors ***************************************************
      // every DSG-modification is equivalent to a certain
      // linear combination of nodal coordinates / displacements
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_shear_r_;
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_shear_s_;
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_membrane_r_;
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_membrane_s_;
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_membrane_rs_;
      std::vector<CORE::LINALG::Matrix<9, 9>> dsg_transverse_t_;
      // DSG factors ***************************************************


      // internal calculation methods

      //! don't want = operator
      SoSh18& operator=(const SoSh18& old);


      //! init the inverse of the jacobian and its determinant in the material configuration
      //! return the number of negative eigenvalues
      int InitJacobianMapping() override;

      //! Calculate nonlinear stiffness and mass matrix
      void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,           ///< current displacements
          std::vector<double>& residual,       ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>*
              stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          CORE::LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
          ) override;

      // parameter space coords of one node
      CORE::LINALG::Matrix<3, 1> NodeParamCoord(const int node);
      // parameter space coords of all nodes
      CORE::LINALG::Matrix<18, 3> NodeParamCoord();

      void FlipT();

      //! Lump mass matrix
      void soh18_lumpmass(CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* emass);

      void EvaluateT(const CORE::LINALG::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>& jac,
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT);

      void EasSetup(
          std::vector<CORE::LINALG::Matrix<6, num_eas>>& M_gp,  // M-matrix evaluated at GPs
          CORE::LINALG::Matrix<3, 1>& G3_contra,  // contravariant basis vector G3 at element center
          const CORE::LINALG::Matrix<NUMNOD_SOH18, 3>& xrefe);  // material element coords

      void SetupDSG();
      void Integrate_dsg_shear_r(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_shear_r);
      void Integrate_dsg_shear_s(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_shear_s);
      void Integrate_dsg_membrane_r(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_membrane_r);
      void Integrate_dsg_membrane_s(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_membrane_s);
      void Integrate_dsg_membrane_rs(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_membrane_rs);
      void Integrate_dsg_transverse_t(const int gp, CORE::LINALG::Matrix<9, 9>& dsg_transverse_t);

      void CalcConsistentDefgrd(const CORE::LINALG::Matrix<3, 3>& defgrd_disp,
          CORE::LINALG::Matrix<6, 1> glstrain_mod, CORE::LINALG::Matrix<3, 3>& defgrd_mod);

      void CalculateBopLoc(const CORE::LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
          const CORE::LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
          const CORE::LINALG::Matrix<9, 1>& shape_q9, const CORE::LINALG::Matrix<2, 9>& deriv_q9,
          const int gp, CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH18>& bop_loc);

      void CalculateLocStrain(const CORE::LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
          const CORE::LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
          const CORE::LINALG::Matrix<9, 1>& shape_q9, const CORE::LINALG::Matrix<2, 9>& deriv_q9,
          const int gp, CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& lstrain);

      void CalculateGeoStiff(const CORE::LINALG::Matrix<9, 1>& shape_q9,
          const CORE::LINALG::Matrix<2, 9>& deriv_q9,
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT, const int gp,
          const double detJ_w, const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress,
          CORE::LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix);

      void Update();


      /** recover elementwise stored stuff */
      void Recover(const std::vector<double>& residual) override;

      //@}

     private:
      std::string GetElementTypeString() const { return "SOLIDSH18"; }
    };  // class So_sh18

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
