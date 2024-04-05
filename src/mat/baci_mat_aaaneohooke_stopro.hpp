/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for aneurysmatic artery wall following
Raghavan and Vorp [2000]

The material is a special case of a generalised pover law neo-Hookean material

This material law accounts for spatial variation of the material parameters young and beta
by using an random field, such that every element gets its own set of mat parameters.

The sample of the random field is not calculated in the material routine but globally in MLMC
because it cannot be stored in each element. At the moment beta is the only stochastic parameter.

the input line should read
  MAT 1 MAT_Struct_AAANeoHooke_Stopro YOUNG 1.044E7 BETA 188.1E5 NUE 0.3 DENS 1.0 SIGMA 0.25
CORRLENGTH 5.0

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_AAANEOHOOKE_STOPRO_HPP
#define FOUR_C_MAT_AAANEOHOOKE_STOPRO_HPP

#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for aneurysm wall material
    class AAAneohooke_stopro : public Parameter
    {
     public:
      /// standard constructor
      AAAneohooke_stopro(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_mean_;
      /// Possion's ratio
      const double nue_;
      /// 2nd parameter
      const double beta_mean_;
      /// mass density
      const double density_;

      //
      int init_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class AAAneohooke
  }     // namespace PAR


  class AAAneohooke_stoproType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "AAAneohooke_stoproType"; }

    static AAAneohooke_stoproType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static AAAneohooke_stoproType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// aneurysm wall material according to Raghavan and Vorp [2000]
  class AAAneohooke_stopro : public So3Material
  {
   public:
    // empty constructor
    AAAneohooke_stopro();

    // constructor with given material parameters
    AAAneohooke_stopro(MAT::PAR::AAAneohooke_stopro* params);

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return AAAneohooke_stoproType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

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

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    // stochastic beta
    double Beta() const { return beta_; }

    // stochastic youngs
    double Youngs() const { return youngs_; }

    /// material mass density
    double Density() const override { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5 * params_->youngs_mean_ / (1.0 + params_->nue_); }

    // material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_aaaneohooke_stopro;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new AAAneohooke_stopro(*this));
    }

    // Init stochastic beta
    void Init(double value_stopro, std::string stochpar);

    // THE material routine
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
        const int eleGID) override;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// evaluate strain energy function
    void StrainEnergy(const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, const int gp,
        const int eleGID) override;

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// my material parameters
    MAT::PAR::AAAneohooke_stopro* params_;

    // stochastic beta
    // put beta and youngs not in params because it will be elementspecific
    double beta_;
    /// init flag for stochastic beta
    bool isinit_beta_;

    double youngs_;
    bool isinit_youngs_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
