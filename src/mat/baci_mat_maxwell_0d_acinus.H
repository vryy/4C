/*----------------------------------------------------------------------*/
/*! \file

\brief Base of four-element Maxwell material model for reduced dimensional
acinus elements

Four-element Maxwell model consists of a parallel configuration of a spring (Stiffness1),
spring-dashpot (Stiffness2 and Viscosity1) and dashpot (Viscosity2) element
(derivation: see Ismail Mahmoud's dissertation, chapter 3.4)


\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_MAXWELL_0D_ACINUS_H
#define BACI_MAT_MAXWELL_0D_ACINUS_H


#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"
#include "baci_red_airways_elem_params.h"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell_0d_acinus : public Parameter
    {
     public:
      /// standard constructor
      Maxwell_0d_acinus(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// first stiffness of the Maxwell model
      const double stiffness1_;
      /// first stiffness of the Maxwell model
      const double stiffness2_;
      /// first viscosity of the Maxwell model
      const double viscosity1_;
      /// first viscosity of the Maxwell model
      const double viscosity2_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class Maxwell_0d_acinus

  }  // namespace PAR

  class Maxwell_0d_acinusType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "maxwell_0d_acinusType"; }

    static Maxwell_0d_acinusType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Maxwell_0d_acinusType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Maxwell 0D acinar material
  ///
  /// This object exists (several times) at every element
  class Maxwell_0d_acinus : public Material
  {
   public:
    /// construct empty material object
    Maxwell_0d_acinus();

    /// construct the material object given material parameters
    explicit Maxwell_0d_acinus(MAT::PAR::Maxwell_0d_acinus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return Maxwell_0d_acinusType::Instance().UniqueParObjectId();
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

    /*!
      \brief
    */
    virtual void Setup(INPUT::LineDefinition* linedef)
    {
      dserror(
          "Setup not implemented yet! Check your material type, "
          "maybe you are still using the base class MAT_0D_MAXWELL_ACINUS.");
    }

    /*!
       \brief
     */
    virtual void Evaluate(CORE::LINALG::SerialDenseVector& epnp,
        CORE::LINALG::SerialDenseVector& epn, CORE::LINALG::SerialDenseVector& epnm,
        CORE::LINALG::SerialDenseMatrix& sysmat, CORE::LINALG::SerialDenseVector& rhs,
        const DRT::REDAIRWAYS::ElemParams& params, const double NumOfAcini, const double Vo,
        double time, double dt)
    {
      dserror("Evaluate not implemented yet !");
    }

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_0d_maxwell_acinus;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Maxwell_0d_acinus(*this));
    }

    /// return density
    double Density() const override { return -1; }

    /// return first stiffness of the Maxwell model
    double Stiffness1() const { return params_->stiffness1_; }

    /// return first stiffness of the Maxwell model
    double Stiffness2() const { return params_->stiffness2_; }

    /// return first viscosity of the Maxwell model
    double Viscosity1() const { return params_->viscosity1_; }

    /// return first viscosity of the Maxwell model
    double Viscosity2() const { return params_->viscosity2_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Return value of class parameter
    virtual double GetParams(std::string parametername);

    /// Set value of class parameter
    virtual void SetParams(std::string parametername, double new_value);

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names){
        /* do nothing for simple material models */};

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int eleID)
    { /* do nothing for simple material models */
      return false;
    };

   protected:
    /// my material parameters
    MAT::PAR::Maxwell_0d_acinus* params_;
  };

}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
