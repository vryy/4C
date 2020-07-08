/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIAL_H
#define MICROMATERIAL_H


#include <Teuchos_RCP.hpp>

#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"



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

      /// destructor
      virtual ~MicroMaterial() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class MicroMaterial

  }  // namespace PAR

  class MicroMaterialType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MicroMaterialType"; }

    static MicroMaterialType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
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
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_struct_multiscale;
    }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new MicroMaterial(*this)); }

    /// evaluate micro material on a processor with macro scale
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int gp, int eleGID);

    /// evaluate micro material on a processor which only knows about the micro scale (supporting
    /// proc)
    void Evaluate(LINALG::Matrix<3, 3>* defgrd, LINALG::Matrix<6, 6>* cmat,
        LINALG::Matrix<6, 1>* stress, const int gp, const int ele_ID, const int microdisnum,
        double V0, bool eleowner);

    virtual double Density() const;

    /// Calculate stresses and strains on the micro-scale
    void PrepareOutput();

    /// Write output on micro-scale
    void Output();

    // Update state vectors
    virtual void Update();

    /// Read restart of micro scale on a processor with macro scale
    void ReadRestart(const int gp, const int eleID, const bool eleowner);

    /// restart micro material on a processor which only knows about the micro scale (supporting
    /// proc)
    void ReadRestart(
        const int gp, const int eleID, const bool eleowner, int microdisnum, double V0);

    /// Prepare new run of inverse analysis
    void InvAnaInit(bool eleowner, int eleID);

    /// @name Access parameters
    //@{
    std::string MicroInputFileName() const { return params_->microfile_; }
    int MicroDisNum() const { return params_->microdisnum_; }
    double InitVol() const { return params_->initvol_; }
    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    std::map<int, Teuchos::RCP<MicroMaterialGP>> matgp_;

    double density_;

    /// my material parameters
    MAT::PAR::MicroMaterial* params_;
  };

}  // namespace MAT

#endif
