/*----------------------------------------------------------------------*/
/*! \file
\brief
Former file of Ursula Mayer

\level 3

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef CARREAUYASUDA_H
#define CARREAUYASUDA_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class CarreauYasuda : public Parameter
    {
     public:
      /// standard constructor
      CarreauYasuda(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~CarreauYasuda() { ; }

      /// @name material parameters
      //@{

      const double nu_0_;    /* zero-shear viscosity */
      const double nu_inf_;  /* infinite-shear viscosity */
      const double lambda_;  /* characteristic time */
      const double a_param_; /* constant parameter */
      const double b_param_; /* constant parameter */
      const double density_; /* density */

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class CarreauYasuda

  }  // namespace PAR

  class CarreauYasudaType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "CarreauYasudaType"; }

    static CarreauYasudaType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static CarreauYasudaType instance_;
  };

  /// Nonlinear viscosity according to Carreau-Yasuda
  class CarreauYasuda : public Material
  {
   public:
    /// construct empty material object
    CarreauYasuda();

    /// construct the material object given material parameters
    explicit CarreauYasuda(MAT::PAR::CarreauYasuda* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const
    {
      return CarreauYasudaType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_carreauyasuda; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new CarreauYasuda(*this)); }

    /// return material parameters for element calculation
    //@{

    /// return parameter for zero-shear viscosity
    double Nu0() const { return params_->nu_0_; }
    /// return parameter for infinite-shear viscosity
    double NuInf() const { return params_->nu_inf_; }
    /// parameter for characteristic time
    double Lambda() const { return params_->lambda_; }
    /// constant parameter
    double AParam() const { return params_->a_param_; }
    /// constant parameter
    double BParam() const { return params_->b_param_; }
    /// density
    virtual double Density() const { return params_->density_; }

    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::CarreauYasuda* params_;
  };

}  // namespace MAT

#endif
