/*----------------------------------------------------------------------*/
/*!
 \brief scatra material for chemical diffusion under a variational framework

   \level 3

 \maintainer  Martin Kronbichler
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_SCATRA_MAT_VAR_CHEMDIFFUSION_H_
#define SRC_DRT_MAT_SCATRA_MAT_VAR_CHEMDIFFUSION_H_


#include "scatra_mat.H"

namespace MAT
{
  namespace PAR
  {
    //! models for half cell open circuit potential of electrode
    enum Models
    {
      chemdiff_undefined,
      chemdiff_linear,
      chemdiff_fickean
    };  // enum Models
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraMatVarChemDiffusion : public ScatraMat
    {
     public:
      //! model for chemical diffusion
      const Models chemdiff_model_;

      /// standard constructor
      ScatraMatVarChemDiffusion(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ScatraMatVarChemDiffusion() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// Reference chemical potential
      const double refMu_;
      /// Reference concentration
      const double refC_;
      /// Reference Temperature
      const double refTemp_;
      /// Gas constant
      const double gasConstant_;

     private:
      //! convert string to model for constitutive description
      Models StringToModel(const std::string& modelstring) const;
    };
    // class Scatra

  }  // namespace PAR

  class ScatraMatVarChemDiffusionType : public ScatraMatType
  {
   public:
    std::string Name() const { return "ScatraMatVarChemDiffusionType"; }

    static ScatraMatVarChemDiffusionType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ScatraMatVarChemDiffusionType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMatVarChemDiffusion : public ScatraMat
  {
   private:
    /// my material parameters
    MAT::PAR::ScatraMatVarChemDiffusion* params_;
    //  int nsd_; //DRT::Problem::Instance()->ScalarTransportDynamicParams().get<int>("NSD")
    //  int numscal_; //DRT::Problem::Instance()->ScalarTransportDynamicParams().get<int>("NSD")

   public:
    /// construct empty material object
    ScatraMatVarChemDiffusion();  //:
                                  //    ScatraMatVarChemDiffusion()
    //  ,nsd_(DRT::Problem::Instance()->ScalarTransportDynamicParams().get<int>("NSD"))
    //  ,numscal_(DRT::Problem::Instance()->ScalarTransportDynamicParams().get<int>("NUMSCAL"))
    //  {return;}

    /// construct the material object given material parameters
    explicit ScatraMatVarChemDiffusion(MAT::PAR::ScatraMatVarChemDiffusion* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return ScatraMatVarChemDiffusionType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_var_chemdiffusion;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ScatraMatVarChemDiffusion(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// returns constitutive model type
    virtual MAT::PAR::Models ModelType() const { return params_->chemdiff_model_; }

    /// return reference Chemical potential
    virtual double RefMu() const { return params_->refMu_; }

    /// return reference concentration
    virtual double RefC() const { return params_->refC_; }

    /// return reference temperature
    virtual double RefTemp() const { return params_->refTemp_; }

    /// return gas Constant
    virtual double GasConstant() const { return params_->gasConstant_; }

    //! Compute the internal energy of the system and its derivatives
    double ComputeInternalEnergy(const double concentration,  //!< Concentration
        const double refMu,                                   //!< Reference chemical potential
        const double refC,                                    //!< Reference concentration
        const double rt,                                      //!< factor RT
        const int orderderivative                             //!< Order of the derivative to use
        ) const;

    double ComputeDualInternalEnergy(const double ChemPot,  //!< Chemical Potential
        const double refMu,                                 //!< Reference chemical potential
        const double refC,                                  //!< Reference concentration
        const double rt,                                    //!< factor RT
        const int orderderivative                           //!< Order of the derivative to use
        ) const;

    template <int NSD>
    void ComputeDissipationPot(double concentration  //!< Concentration at t_(n+1)
        ,
        const double& refC  //!< Reference concentration
        ,
        const LINALG::Matrix<NSD, 1>& chemicalfield  //!< Chemical field at t_(n+1)
        ,
        const double& mobility  //!< mobility
        ,
        double& DissipationPot  //!< Output
        ) const;

    template <int NSD>
    void ComputeDissipationPot_D1(double concentration  //!< Concentration at t_(n+1)
        ,
        const double& refC  //!< Reference concentration
        ,
        const LINALG::Matrix<NSD, 1>& chemicalfield  //!< Chemical field at t_(n+1)
        ,
        const double& mobility  //!< mobility
        ,
        LINALG::Matrix<NSD, 1>& DissipationPot1deriv  //!< Output
        ) const;

    template <int NSD>
    void ComputeDissipationPot_D2(double concentration  //!< Concentration at t_(n+1)
        ,
        const double& refC  //!< Reference concentration
        ,
        const LINALG::Matrix<NSD, 1>& chemicalfield  //!< Chemical field at t_(n+1)
        ,
        const double& mobility  //!< mobility
        ,
        LINALG::Matrix<NSD, NSD>& DissipationPot2deriv  //!< Output
        ) const;
  };

}  // namespace MAT

#endif /* SRC_DRT_MAT_SCATRA_MAT_MULTIPORO_H_ */
