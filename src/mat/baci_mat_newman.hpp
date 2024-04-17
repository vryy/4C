/*----------------------------------------------------------------------------*/
/*! \file
\brief material stores parameters for ion species in electrolyte solution. The newman material is
derived for a binary electrolyte using the electroneutrality condition to condense the non-reacting
species

\level 2


*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_NEWMAN_HPP
#define FOUR_C_MAT_NEWMAN_HPP

#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_elchsinglemat.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class Newman : public ElchSingleMat
    {
     public:
      /// standard constructor
      Newman(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// valence (= charge number)
      const double valence_;

      /// definition of transference number
      /// (by curve number or implemented concentration dependence)
      const int transnrcurve_;

      /// definition of thermodynamic factor
      /// (by curve number or implemented concentration dependence)
      const int thermfaccurve_;

      // Important:
      // pointer to vectors containing parameter for predefined curves
      // -> faster than saving the actual vector

      /// number of parameter needed for implemented concentration dependence
      const int transnrparanum_;
      /// parameter needed for implemented concentration dependence
      const std::vector<double> transnrpara_;

      /// number of parameter needed for implemented concentration dependence
      const int thermfacparanum_;
      /// parameter needed for implemented concentration dependence
      const std::vector<double> thermfacpara_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };  // class Newman

  }  // namespace PAR

  class NewmanType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "NewmanType"; }

    static NewmanType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static NewmanType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Newman : public ElchSingleMat
  {
   public:
    /// construct empty material object
    Newman();

    /// construct the material object given material parameters
    explicit Newman(MAT::PAR::Newman* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return NewmanType::Instance().UniqueParObjectId(); }

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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_newman; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Newman(*this)); }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }

    /// computation of the transference number based on the defined curve
    double ComputeTransferenceNumber(const double cint) const;
    /// computation of the first derivative of the transference number based on the defined curve
    double ComputeFirstDerivTrans(const double cint) const;

    /// computation of the thermodynamic factor based on the defined curve
    double ComputeThermFac(const double cint) const;
    /// computation of the first derivative of the transference number based on the defined curve
    double ComputeFirstDerivThermFac(const double cint) const;

   private:
    /// return curve defining the transference number
    int TransNrCurve() const { return params_->transnrcurve_; }
    /// return curve defining the thermodynamic factor
    int ThermFacCurve() const { return params_->thermfaccurve_; }

    /// parameter needed for implemented concentration dependence
    const std::vector<double>& TransNrParams() const { return params_->transnrpara_; }

    /// parameter needed for implemented concentration dependence
    const std::vector<double>& ThermFacParams() const { return params_->thermfacpara_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// my material parameters
    MAT::PAR::Newman* params_;
  };
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
