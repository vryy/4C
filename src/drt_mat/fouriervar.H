/*----------------------------------------------------------------------*/
/*! \file
\brief material for heat conduction according to fourier's law with variable conductivity

\level 1
\maintainer Sebastian Proell
*/

/*----------------------------------------------------------------------*
 |  definitions                                                         |
 *----------------------------------------------------------------------*/
#ifndef FOURIERVAR_H_
#define FOURIERVAR_H_

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "material_thermo.H"
#include "consolidation.H"

#include "../drt_lib/drt_parobjectfactory.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for FourierVar material
    ///
    /// <h3>Input line</h3>
    /// MAT 2 THERM_FourierVar CAPAFUNCT 1 2 3 CONDUCTFUNCT 4 5 6 CONSOLMAT 3
    class FourierVar : public Parameter
    {
     public:
      /// standard constructor
      explicit FourierVar(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      ~FourierVar() override = default;

      /// @name material parameters
      //@{

      /// functions for volumetric heat capacity
      const std::vector<int> capafunct_;
      /// functions for heat conductivity evaluation
      const std::vector<int> conductfunct_;
      /// material specifying the consolidation implementation
      const int consolmat_;

      //@}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class FourierVar

  }  // namespace PAR

  class FourierVarType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "FourierVarType"; }

    static FourierVarType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FourierVarType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// FourierVar material
  ///
  /// This is a an extension of Fourier's isotropic law of heat conduction to
  /// multiple phases and varying heat capacity and conductivity.
  ///
  /// the input line should read something like:
  /// MAT 2 THERM_FourierVar CAPANUM § CAPAFUNCT 5 6 6 CONDUCTNUM 3 CONDUCTFUNCT 4 3 4 CONSOLMAT 3
  ///
  /// \author proell
  /// \date 04/18
  class FourierVar : public ThermoMaterial
  {
   public:
    /// empty constructor
    FourierVar();

    /// constructor with given material parameters
    explicit FourierVar(MAT::PAR::FourierVar* params);

    /// @name Packing and Unpacking
    //@{

    int UniqueParObjectId() const override
    {
      return FourierVarType::Instance().UniqueParObjectId();
    }

    void Pack(DRT::PackBuffer& data) const override;

    void Unpack(const std::vector<char>& data) override;

    //@}

    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_th_fourier_var; }

    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new FourierVar(*this)); }


    /// setup rotine to create history variables based on number of Gauss points
    void Setup(const int numgp);

    /*!
     * @brief set Consolidation object to be used by this material.
     *
     * This method is intended for classes that want to reuse this implementation with a shared
     * consolidation history.
     *
     * @param consolidation the Consolidation object to be used
     */
    void InjectConsolidation(Teuchos::RCP<MAT::Consolidation> consolidation);

    /// @name Material evaluation and parameters
    //! @{

    void Evaluate(const LINALG::Matrix<3, 1>& gradtemp, LINALG::Matrix<3, 3>& cmat,
        LINALG::Matrix<3, 1>& heatflux) const override
    {
      EvaluateInternal(gradtemp, cmat, heatflux);
    }

    void Evaluate(const LINALG::Matrix<2, 1>& gradtemp, LINALG::Matrix<2, 2>& cmat,
        LINALG::Matrix<2, 1>& heatflux) const override
    {
      EvaluateInternal(gradtemp, cmat, heatflux);
    }

    void Evaluate(const LINALG::Matrix<1, 1>& gradtemp, LINALG::Matrix<1, 1>& cmat,
        LINALG::Matrix<1, 1>& heatflux) const override
    {
      EvaluateInternal(gradtemp, cmat, heatflux);
    }

    void ConductivityDerivT(LINALG::Matrix<3, 3>& dCondDT) const override
    {
      GetConductivityDerivTInternal(dCondDT);
    }

    void ConductivityDerivT(LINALG::Matrix<2, 2>& dCondDT) const override
    {
      GetConductivityDerivTInternal(dCondDT);
    }

    void ConductivityDerivT(LINALG::Matrix<1, 1>& dCondDT) const override
    {
      GetConductivityDerivTInternal(dCondDT);
    }

    double Capacity() const override;

    double CapacityDerivT() const override;

    //! @}

    void Reinit(double temperature, unsigned gp) override
    {
      currentTemperature = temperature;
      currentGp = gp;
    };

    void ResetCurrentState() override { consolidation_->ResetLatentHeatLastStep(); }

    void CommitCurrentState() override { consolidation_->Update(currentTemperature, currentGp); }

    /// Get a capacity value useable for heat integration method
    double HeatIntegrationCapacity() const;

    /// get a pointer to the (converged!) consolidation fraction, for postproc
    template <unsigned int nen>
    void NodalPhase(LINALG::Matrix<nen, 1>& nodaltemp,  ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalphase              ///< (return) fraction of solid material
    )
    {
      consolidation_->NodalPhase(nodaltemp, nodalphase);
    }


    /// evaluate conductivity and capacity at nodes
    ///
    /// only for post-processing purposes
    template <unsigned int nen>
    void NodalConductivityCapacity(
        LINALG::Matrix<nen, 1>& nodaltemp,  ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalcond,  ///< (return) conductivity at nodes
        LINALG::Matrix<nen, 1>& nodalcapa   ///< (return) capacity at nodes
    )
    {
      consolidation_->EvaluateTempDependentFunctionAtNodes(
          nodaltemp, nodalcond, params_->conductfunct_);
      consolidation_->EvaluateTempDependentFunctionAtNodes(
          nodaltemp, nodalcapa, params_->capafunct_);
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Get direct access to the Consolidation object
    virtual Teuchos::RCP<MAT::Consolidation> Consolidation() { return consolidation_; }

   private:
    //! templated internal evaluate implementation
    template <unsigned nsd>
    void EvaluateInternal(const LINALG::Matrix<nsd, 1>& gradtemp, LINALG::Matrix<nsd, nsd>& cmat,
        LINALG::Matrix<nsd, 1>& heatflux) const;

    template <unsigned nsd>
    void GetConductivityDerivTInternal(LINALG::Matrix<nsd, nsd>& dCondDT) const;

    /// my material parameters
    MAT::PAR::FourierVar* params_;

    /// pointer to the consolidation material managing the actual evaluation of function
    Teuchos::RCP<MAT::Consolidation> consolidation_;

    /// temperature used for evaluation
    double currentTemperature{};

    /// current evaluation Gauss point
    unsigned currentGp{};

  };  // FourierVar

}  // namespace MAT
/*----------------------------------------------------------------------*/
#endif  // FOURIERVAR_H_
