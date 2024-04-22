/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MICROMATERIAL_HPP
#define FOUR_C_MAT_MICROMATERIAL_HPP


#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN



namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for micro material
    class MicroMaterial : public Parameter
    {
     public:
      /// standard constructor
      MicroMaterial(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// inputfile for microstructure
      const std::string microfile_;
      /// Number of microscale discretization
      const int microdisnum_;
      ///
      const double initvol_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class MicroMaterial

  }  // namespace PAR

  class MicroMaterialType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "MicroMaterialType"; }

    static MicroMaterialType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MicroMaterialType instance_;
  };

  /*----------------------------------------------------------------------*/
  // forward
  class MicroMaterialGP;

  /*----------------------------------------------------------------------*/
  /// class for handling of micro-macro transitions
  class MicroMaterial : public So3Material
  {
   public:
    /// construct empty material object
    MicroMaterial();

    /// construct the material object given material parameters
    explicit MicroMaterial(MAT::PAR::MicroMaterial* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return MicroMaterialType::Instance().UniqueParObjectId();
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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_struct_multiscale;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new MicroMaterial(*this)); }

    /// evaluate micro material on a processor with macro scale
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, int gp,
        int eleGID) override;

    /// evaluate micro material on a processor which only knows about the micro scale (supporting
    /// proc)
    void Evaluate(CORE::LINALG::Matrix<3, 3>* defgrd, CORE::LINALG::Matrix<6, 6>* cmat,
        CORE::LINALG::Matrix<6, 1>* stress, const int gp, const int ele_ID, const int microdisnum,
        double V0, bool eleowner);

    double Density() const override;

    /// Calculate stresses and strains on the micro-scale
    void PrepareOutput();

    /// Write output on micro-scale
    void Output();

    // Update state vectors
    void Update() override;

    /// Read restart of micro scale on a processor with macro scale
    void ReadRestart(const int gp, const int eleID, const bool eleowner);

    /// restart micro material on a processor which only knows about the micro scale (supporting
    /// proc)
    void ReadRestart(
        const int gp, const int eleID, const bool eleowner, int microdisnum, double V0);

    /// @name Access parameters
    //@{
    std::string MicroInputFileName() const { return params_->microfile_; }
    int MicroDisNum() const { return params_->microdisnum_; }
    double InitVol() const { return params_->initvol_; }
    //@}

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    std::map<int, Teuchos::RCP<MicroMaterialGP>> matgp_;

    double density_;

    /// my material parameters
    MAT::PAR::MicroMaterial* params_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
