/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for aneurysmatic artery wall
damaging following Simo approach explained in Holzapfel (p298).
SEF by Raghavan and Vorp[2000](Isochoric) and Simo&Miehe version of
Ogden (1972) (Volumetric) with beta=-2

The material is a special case of a generalised pover law neo-Hookean material
with a damage model inside relative only to the isochoric part.

the input line should read
///  MAT 1 MAT_Raghavan_Damage BULK 1.044E7 ALPHA 6.86E-2 BETA 188.1E5 INDAMSTR 1,25654 A 0.3 B 121
DENS 1.0 MAT 3 MAT_Raghavan_Damage BULK 0.120755 ALPHA 0.068632  BETA 5.799445 EQSTRMIN 0.206141
EQSTRMAX 0.6444974  A 2.4816557 B  0.478626783 DENS 0.001

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_AAARAGHAVANVORP_DAMAGE_HPP
#define FOUR_C_MAT_AAARAGHAVANVORP_DAMAGE_HPP

#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for aneurysm wall material
    class AAAraghavanvorp_damage : public Parameter
    {
     public:
      /// standard constructor
      AAAraghavanvorp_damage(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      const double bulk_;   /// Bulk's modulus (Volumetric)
      const double alpha_;  /// 1st parameter, alpha (Isochoric)
      const double beta_;   /// 2nd parameter, beta (Isochoric)

      /// damage parameters

      const double eqstrmin_;  /// equivalent strain initial damage
      const double a_;         /// 1st parameter, a
      const double b_;         /// 2nd parameter, b
      const double density_;   /// mass density

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class AAAraghavanvorp_damage

  }  // namespace PAR

  class AAAraghavanvorp_damageType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "AAAraghavanvorp_damageType"; }

    static AAAraghavanvorp_damageType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static AAAraghavanvorp_damageType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// aneurysm wall material with damage according Simo approach explained in Holzapfel (p298)
  class AAAraghavanvorp_damage : public So3Material
  {
   public:
    // empty constructor
    AAAraghavanvorp_damage();

    // constructor with given material parameters
    AAAraghavanvorp_damage(MAT::PAR::AAAraghavanvorp_damage* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return AAAraghavanvorp_damageType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
      History data is unpacked in restart.

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material mass density
    double Density() const override { return params_->density_; }

    /// shear modulus
    double ShearMod() const
    {
      dserror("");
      return (6.0 * params_->bulk_ * params_->alpha_) /
             (3.0 * params_->bulk_ - 2 * params_->alpha_);
    }

    /// Check if history variables are already initialized

    bool Initialized() const { return isinit_ && (histgcurr_ != Teuchos::null); }

    // material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_aaaraghavanvorp_damage;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new AAAraghavanvorp_damage(*this));
    }

    /// Initialize internal stress variables
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /// Update internal stress variables
    void Update() override;

    void StressTensTransfSPKtoCauchy(
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& f,  ///< deformation gradient tensor
        const double detf,                          ///< determinant of deformation gradient tensor
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& pktwo,   ///< 2nd PK-stress
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& cstress  ///< Cauchy-stress
    );

    void StressTensTransfCauchytoSPK(
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& invf,  ///< deformation gradient tensor
        const double detf,  ///< determinant of deformation gradient tensor
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& cstress,  ///< Cauchy-stress
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& pktwo     ///< 2nd PK-stress
    );

    // THE material routine
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,  ///< green lagrange strain
        Teuchos::ParameterList& params,                  ///< parameter list for communication
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  ///< 2nd PK-stress
        CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
        int gp,                                                    ///< Gauss point
        int eleGID) override;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::AAAraghavanvorp_damage* params_;
    // damage history parameter g
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<1, 1>>> histgcurr_;  ///< current damage parameter
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<1, 1>>>
        histglast_;  ///< damage of last converged state

    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<1, 1>>>
        histeqstrmaxcurr_;  ///< current damage parameter
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<1, 1>>>
        histeqstrmaxlast_;             ///< damage of last converged state
    Teuchos::RCP<double> elstrength_;  ///< strength of the element

    bool isinit_;  ///< indicates if material is initialized by calling the #Initialized routine
  };
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
