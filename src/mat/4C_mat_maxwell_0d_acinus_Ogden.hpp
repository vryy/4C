/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional
acinus elements with non-linear Ogden-like spring, inherits from Maxwell_0d_acinus

The originally linear spring (Stiffness1) of the 4-element Maxwell model is substituted by a
non-linear pressure-volume relation derived from the Ogden strain energy function considering pure
volumetric expansion (derivation: see Christian Roth's dissertation, Appendix B)


\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_MAXWELL_0D_ACINUS_OGDEN_HPP
#define FOUR_C_MAT_MAXWELL_0D_ACINUS_OGDEN_HPP


#include "4C_config.hpp"

#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_red_airways_elem_params.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell0dAcinusOgden : public Maxwell0dAcinus
    {
     public:
      /// standard constructor
      Maxwell0dAcinusOgden(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// enum for mapping between material parameter and entry in the matparams_ vector
      enum Matparamnames
      {
        kappa,
        beta,
        first = kappa,
        last = beta
      };

    };  // class Maxwell_0d_acinus_Ogden
  }     // namespace PAR


  class Maxwell0dAcinusOgdenType : public Maxwell0dAcinusType
  {
   public:
    std::string Name() const override { return "maxwell_0d_acinusOgdenType"; }

    static Maxwell0dAcinusOgdenType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Maxwell0dAcinusOgdenType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Maxwell 0D acinar material
  ///
  /// This object exists (several times) at every element
  class Maxwell0dAcinusOgden : public Maxwell0dAcinus
  {
   public:
    /// construct empty material object
    Maxwell0dAcinusOgden();

    /// construct the material object given material parameters
    Maxwell0dAcinusOgden(Mat::PAR::Maxwell0dAcinus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return Maxwell0dAcinusOgdenType::Instance().UniqueParObjectId();
    }


    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(Core::Communication::PackBuffer& data) const override;

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
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_0d_maxwell_acinus_ogden;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new Maxwell0dAcinus(*this));
    }

    /*!
      \brief
    */
    void Setup(Input::LineDefinition* linedef) override;

    /*!
       \brief
     */
    void Evaluate(Core::LinAlg::SerialDenseVector& epnp, Core::LinAlg::SerialDenseVector& epn,
        Core::LinAlg::SerialDenseVector& epnm, Core::LinAlg::SerialDenseMatrix& sysmat,
        Core::LinAlg::SerialDenseVector& rhs, const Discret::ReducedLung::ElemParams& params,
        const double NumOfAcini, const double Vo, double time, double dt) override;

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int eleGID) override;

    /// Return value of class parameter
    double GetParams(std::string parametername) override;

    /// Set value of class parameter
    void SetParams(std::string parametername, double new_value) override;

   private:
    double kappa_;
    double beta_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
