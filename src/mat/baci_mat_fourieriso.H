/*----------------------------------------------------------------------*/
/*! \file
\brief material for heat conduction according to fourier's law
\level 1
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/
#ifndef BACI_MAT_FOURIERISO_H
#define BACI_MAT_FOURIERISO_H

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include "baci_lib_parobjectfactory.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_thermo.H"

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for FourierIso material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 THERM_FourierIsoIso CAPA 1.0 COND 1.0
    class FourierIso : public Parameter
    {
     public:
      /// standard constructor
      FourierIso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// volumetric heat capacity
      const double capa_;
      /// heat conductivity
      const double conduct_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class FourierIso

  }  // namespace PAR

  class FourierIsoType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "FourierIsoType"; }

    static FourierIsoType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FourierIsoType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// FourierIso material according to [1]
  ///
  /// This is a FourierIso's law of isotropic, instationary heat conduction
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author dano
  /// \date 09/09
  class FourierIso : public ThermoMaterial
  {
   public:
    /// empty constructor
    FourierIso();

    /// constructor with given material parameters
    FourierIso(MAT::PAR::FourierIso* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return FourierIsoType::Instance().UniqueParObjectId();
    }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    void Pack(DRT::PackBuffer& data  ///< (in/out): char vector to store class information
    ) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    void Unpack(const std::vector<char>& data  ///< vector storing all data to be unpacked into this
        ) override;

    //@}

    /// @name Access material constants
    //@{

    /// conductivity
    double Conductivity() const { return params_->conduct_; }

    /// volumetric heat capacity
    double Capacity() const override { return params_->capa_; }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_th_fourier_iso; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new FourierIso(*this)); }

    //@}

    void Evaluate(const CORE::LINALG::Matrix<1, 1>& gradtemp, CORE::LINALG::Matrix<1, 1>& cmat,
        CORE::LINALG::Matrix<1, 1>& heatflux) const override;

    void Evaluate(const CORE::LINALG::Matrix<2, 1>& gradtemp, CORE::LINALG::Matrix<2, 2>& cmat,
        CORE::LINALG::Matrix<2, 1>& heatflux) const override;

    void Evaluate(const CORE::LINALG::Matrix<3, 1>& gradtemp, CORE::LINALG::Matrix<3, 3>& cmat,
        CORE::LINALG::Matrix<3, 1>& heatflux) const override;

    void ConductivityDerivT(CORE::LINALG::Matrix<3, 3>& dCondDT) const override { dCondDT.Clear(); }

    void ConductivityDerivT(CORE::LINALG::Matrix<2, 2>& dCondDT) const override { dCondDT.Clear(); }

    void ConductivityDerivT(CORE::LINALG::Matrix<1, 1>& dCondDT) const override { dCondDT.Clear(); }

    double CapacityDerivT() const override { return 0; }

    void Reinit(double temperature, unsigned gp) override
    {
      // do nothing
    }

    void ResetCurrentState() override
    {
      // do nothing
    }

    void CommitCurrentState() override
    {
      // do nothing
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::FourierIso* params_;

  };  // FourierIso

}  // namespace MAT

/*----------------------------------------------------------------------*/
#endif  // MAT_FOURIERISO_H
