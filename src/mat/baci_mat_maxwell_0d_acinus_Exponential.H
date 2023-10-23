/*----------------------------------------------------------------------*/
/*! \file

\brief Four-element Maxwell material model for reduced dimensional acinus elements with non-linear
spring with exponential behaviour, inherits from Maxwell_0d_acinus

The originally linear spring (Stiffness1) of the 4-element Maxwell model is substituted by an
exponential pressure-volume relation (derivation: see Ismail Mahmoud's dissertation, chapter 3.4)


\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_MAXWELL_0D_ACINUS_EXPONENTIAL_H
#define BACI_MAT_MAXWELL_0D_ACINUS_EXPONENTIAL_H


#include "baci_mat_maxwell_0d_acinus.H"
#include "baci_red_airways_elem_params.h"
#include "baci_red_airways_elementbase.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell_0d_acinus_Exponential : public Maxwell_0d_acinus
    {
     public:
      /// standard constructor
      Maxwell_0d_acinus_Exponential(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class Maxwell_0d_acinus_Exponential
  }     // namespace PAR


  class Maxwell_0d_acinusExponentialType : public Maxwell_0d_acinusType
  {
   public:
    std::string Name() const override { return "maxwell_0d_acinusExponentialType"; }

    static Maxwell_0d_acinusExponentialType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Maxwell_0d_acinusExponentialType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Maxwell 0D acinar material
  ///
  /// This object exists (several times) at every element
  class Maxwell_0d_acinus_Exponential : public Maxwell_0d_acinus
  {
   public:
    /// construct empty material object
    Maxwell_0d_acinus_Exponential();

    /// construct the material object given material parameters
    Maxwell_0d_acinus_Exponential(MAT::PAR::Maxwell_0d_acinus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return Maxwell_0d_acinusExponentialType::Instance().UniqueParObjectId();
    }


    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

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

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;
    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_0d_maxwell_acinus_exponential;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Maxwell_0d_acinus(*this));
    }

    /*!
      \brief
    */
    void Setup(DRT::INPUT::LineDefinition* linedef) override;

    /*!
       \brief
     */
    void Evaluate(CORE::LINALG::SerialDenseVector& epnp, CORE::LINALG::SerialDenseVector& epn,
        CORE::LINALG::SerialDenseVector& epnm, CORE::LINALG::SerialDenseMatrix& sysmat,
        CORE::LINALG::SerialDenseVector& rhs, const DRT::REDAIRWAYS::ElemParams& params,
        const double NumOfAcini, const double Vo, double time, double dt) override;

   private:
    double e1_0_;
    double e1_lin_;
    double e1_exp_;
    double tau_;
  };

}  // namespace MAT

#endif
