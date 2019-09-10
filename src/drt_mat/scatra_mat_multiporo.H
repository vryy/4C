/*----------------------------------------------------------------------*/
/*! \file
 \brief scatra material for transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_
#define SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_


#include "scatra_mat.H"

namespace MAT
{
  namespace ScatraMatMultiPoro
  {
    /*--------------------------------------------------------------------------*/
    /*!
     * \brief enum that provides all possible MatScatraMultiPoro species types
     *///                                                        kremheller 03/18
    /*--------------------------------------------------------------------------*/
    enum SpeciesType
    {
      // species type
      species_in_fluid,    /*!< species in fluid */
      species_in_volfrac,  /*!< species in volume fraction */
      species_in_solid,    /*!< species in solid */
      species_temperature, /*!< temperature */
      species_undefined
    };  // enum SpeciesType

    /**
     *  \brief a structure for the mapping of scalars to the transporting phase
     */
    struct ScalarToPhaseMap
    {
      int phaseID;              /*!< ID of the fluid phase containing the scalar*/
      SpeciesType species_type; /*!< species type */
    };

  }  // namespace ScatraMatMultiPoro


  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material (species in fluid)
    class ScatraMatMultiPoroFluid : public ScatraMat
    {
     public:
      /// standard constructor
      ScatraMatMultiPoroFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMatMultiPoroFluid() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// ID of fluid phase the scalar is associated with
      const int phaseID_;
      /// delta used for modelling dependency of diffusivity on (saturation*porosity)^delta
      //  as in G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and B. A. Schrefler.
      //  Three phase flow dynamics in tumor growth. Computational Mechanics, 53:465–484,
      //  2014.
      const double delta_;
      // minimum saturation under which also corresponding mass fraction is equal to zero
      const double min_sat_;
    };
    // class Scatra

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material (species in volume fraction)
    class ScatraMatMultiPoroVolFrac : public ScatraMat
    {
     public:
      /// standard constructor
      ScatraMatMultiPoroVolFrac(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMatMultiPoroVolFrac() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// ID of fluid phase the scalar is associated with
      const int phaseID_;
      /// delta used for modelling dependency of diffusivity on volfrac^delta
      const double delta_;
    };
    // class Scatra
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material (species in solid)
    class ScatraMatMultiPoroSolid : public ScatraMat
    {
     public:
      /// standard constructor
      ScatraMatMultiPoroSolid(Teuchos::RCP<MAT::PAR::Material> matdata);
      /// destructor
      virtual ~ScatraMatMultiPoroSolid() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// delta used for modelling dependency of diffusivity on
      /// (saturation*porosity)^delta
      //  as in G. Sciume, William G. Gray, F. Hussain, M. Ferrari, P. Decuzzi, and
      //  B. A. Schrefler. Three phase flow dynamics in tumor growth. Computational
      //  Mechanics, 53:465–484, 2014.
      const double delta_;
    };
    // class Scatra
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material (temperature)
    class ScatraMatMultiPoroTemperature : public ScatraMat
    {
     public:
      /// standard constructor
      ScatraMatMultiPoroTemperature(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMatMultiPoroTemperature() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      const int numfluidphases_;
      const int numvolfrac_;
      // heat capacities for fluids, volfracs and solid
      const std::vector<double> cp_fluid_;
      const std::vector<double> cp_volfrac_;
      const double cp_solid_;

      // thermal diffusivity for fluids, volfracs and solid
      const std::vector<double> kappa_fluid_;
      const std::vector<double> kappa_volfrac_;
      const double kappa_solid_;
    };
    // class Scatra

  }  // namespace PAR

  class ScatraMatMultiPoroFluidType : public ScatraMatType
  {
   public:
    std::string Name() const { return "ScatraMatMultiPoroFluidType"; }

    static ScatraMatMultiPoroFluidType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatMultiPoroFluidType instance_;
  };

  class ScatraMatMultiPoroVolFracType : public ScatraMatType
  {
   public:
    std::string Name() const { return "ScatraMatMultiPoroVolFracType"; }

    static ScatraMatMultiPoroVolFracType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatMultiPoroVolFracType instance_;
  };

  class ScatraMatMultiPoroSolidType : public ScatraMatType
  {
   public:
    std::string Name() const { return "ScatraMatMultiPoroSolidType"; }

    static ScatraMatMultiPoroSolidType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatMultiPoroSolidType instance_;
  };
  class ScatraMatMultiPoroTemperatureType : public ScatraMatType
  {
   public:
    std::string Name() const { return "ScatraMatMultiPoroTemperatureType"; }

    static ScatraMatMultiPoroTemperatureType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatMultiPoroTemperatureType instance_;
  };
  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material (species in fluid)
  class ScatraMatMultiPoroFluid : public ScatraMat
  {
   public:
    /// construct empty material object
    ScatraMatMultiPoroFluid();

    /// construct the material object given material parameters
    explicit ScatraMatMultiPoroFluid(MAT::PAR::ScatraMatMultiPoroFluid* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatMultiPoroFluidType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_scatra_multiporo_fluid;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraMatMultiPoroFluid(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// return phase ID
    virtual int PhaseID() const { return params_->phaseID_; }

    /// return delta
    virtual double Delta() const { return params_->delta_; }

    /// return minimum saturation
    virtual double MinSat() const { return params_->min_sat_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatMultiPoroFluid* params_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material (species in volfrac)
  class ScatraMatMultiPoroVolFrac : public ScatraMat
  {
   public:
    /// construct empty material object
    ScatraMatMultiPoroVolFrac();

    /// construct the material object given material parameters
    explicit ScatraMatMultiPoroVolFrac(MAT::PAR::ScatraMatMultiPoroVolFrac* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatMultiPoroVolFracType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_scatra_multiporo_volfrac;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraMatMultiPoroVolFrac(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// return phase ID
    virtual int PhaseID() const { return params_->phaseID_; }

    /// return delta
    virtual double Delta() const { return params_->delta_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatMultiPoroVolFrac* params_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material (species in solid)
  class ScatraMatMultiPoroSolid : public ScatraMat
  {
   public:
    /// construct empty material object
    ScatraMatMultiPoroSolid();

    /// construct the material object given material parameters
    explicit ScatraMatMultiPoroSolid(MAT::PAR::ScatraMatMultiPoroSolid* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatMultiPoroSolidType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_scatra_multiporo_solid;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraMatMultiPoroSolid(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// return delta
    virtual double Delta() const { return params_->delta_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatMultiPoroSolid* params_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material (species in temperature)
  class ScatraMatMultiPoroTemperature : public ScatraMat
  {
   public:
    /// construct empty material object
    ScatraMatMultiPoroTemperature();

    /// construct the material object given material parameters
    explicit ScatraMatMultiPoroTemperature(MAT::PAR::ScatraMatMultiPoroTemperature* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatMultiPoroTemperatureType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_scatra_multiporo_temperature;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraMatMultiPoroTemperature(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    std::vector<double> CP_Fluid() const { return params_->cp_fluid_; }
    double CP_Fluid(int phase) const { return params_->cp_fluid_[phase]; }
    std::vector<double> CP_Volfrac() const { return params_->cp_volfrac_; }
    double CP_Volfrac(int phase) const { return params_->cp_volfrac_[phase]; }
    double CP_Solid() const { return params_->cp_solid_; }

    std::vector<double> KAPPA_Fluid() const { return params_->kappa_fluid_; }
    double KAPPA_Fluid(int phase) const { return params_->kappa_fluid_[phase]; };
    std::vector<double> KAPPA_Volfrac() const { return params_->kappa_volfrac_; }
    double KAPPA_Volfrac(int phase) const { return params_->kappa_volfrac_[phase]; }
    double KAPPA_Solid() const { return params_->kappa_solid_; }

   private:
    /// my material parameters
    MAT::PAR::ScatraMatMultiPoroTemperature* params_;
  };

}  // namespace MAT



#endif /* SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_ */
