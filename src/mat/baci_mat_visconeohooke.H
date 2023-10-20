/*----------------------------------------------------------------------*/
/*! \file
\brief
Viscohyperelastic material model containing the following parts:
IsoNeohooke + VolSussmannBathe with Generalized Maxwell (just on isochoric part).
The input line should read
MAT 1 MAT_VISCONEOHOOKE YOUNGS_SLOW 1.0 POISSON 0.499 DENS 0.1 YOUNGS_FAST 100.0 RELAX 10.0 THETA
0.5

\level 2


*----------------------------------------------------------------------*/
#ifndef BACI_MAT_VISCONEOHOOKE_H
#define BACI_MAT_VISCONEOHOOKE_H

#include "baci_lib_parobjectfactory.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_so3_material.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class ViscoNeoHooke : public Parameter
    {
     public:
      /// standard constructor
      ViscoNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      const double youngs_slow_;
      const double poisson_;
      const double density_;
      const double youngs_fast_;
      const double relax_;
      const double theta_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class ViscoNeoHooke

  }  // namespace PAR

  class ViscoNeoHookeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ViscoNeoHookeType"; }

    static ViscoNeoHookeType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ViscoNeoHookeType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Visco-NeoHooke material
  class ViscoNeoHooke : public So3Material
  {
   public:
    /// construct empty material object
    ViscoNeoHooke();

    /// construct the material object given material parameters
    explicit ViscoNeoHooke(MAT::PAR::ViscoNeoHooke* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ViscoNeoHookeType::Instance().UniqueParObjectId();
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
    void Pack(DRT::PackBuffer& data) const override;

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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_visconeohooke; }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new ViscoNeoHooke(*this)); }

    /// Initialize internal stress variables
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /// Update internal stress variables
    void Update() override;

    /// Evaluate material
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,      ///< deformation gradient
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,  ///< green lagrange strain
        Teuchos::ParameterList& params,                  ///< parameter list for communication
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  ///< 2nd PK-stress
        CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
        int gp,                                                    ///< Gauss point
        int eleGID) override;

    /// Return density
    double Density() const override { return params_->density_; }

    /// Return shear modulus
    double ShearMod() const { return 0.5 * params_->youngs_slow_ / (1.0 + params_->poisson_); };

    /// Check if history variables are already initialized
    bool Initialized() const { return isinit_ && (histstresscurr_ != Teuchos::null); }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ViscoNeoHooke* params_;

    /// visco history stresses
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histstresscurr_;  ///< current stress
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histstresslast_;  ///< stress of last converged state
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        artstresscurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        artstresslast_;  ///< artificial stress in last converged state

    bool isinit_;  ///< indicates if #Initialized routine has been called
  };
}  // namespace MAT

#endif  // MAT_VISCONEOHOOKE_H
