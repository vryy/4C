/*----------------------------------------------------------------------*/
/*! \file
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       example input line:
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732

\level 2


 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ELECTROMAGNETIC_HPP
#define FOUR_C_MAT_ELECTROMAGNETIC_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for electromagnetic wave propagation
    class ElectromagneticMat : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      ElectromagneticMat(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      enum Matparamnames
      {
        sigma_,
        epsilon_,
        mu_,
        first = sigma_,
        last = mu_
      };

      double get_parameter(int parametername, const int EleId)
      {
        // check if we have an element based value via size
        if (matparams_[parametername]->GlobalLength() == 1)
        {
          // we have a global value hence we directly return the first entry
          return (*matparams_[parametername])[0];
        }
        // If someone calls this functions without a valid EleID and we have element based values
        // throw error
        else if (EleId < 0 && matparams_[parametername]->GlobalLength() > 1)
        {
          FOUR_C_THROW("Global mat parameter requested but we have elementwise mat params");
          return 0.0;
        }
        // otherwise just return the element specific value
        else
        {
          // calculate LID here, instead of before each call
          return (*matparams_[parametername])[matparams_[parametername]->Map().LID(EleId)];
        }
      }

     private:
      std::vector<Teuchos::RCP<Epetra_Vector>> matparams_;

    };  // class ElectromagneticMat

  }  // namespace PAR

  class ElectromagneticMatType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "ElectromagneticMatType"; }

    static ElectromagneticMatType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    static ElectromagneticMatType instance_;

  };  // class ElectromagneticMatType

  /*----------------------------------------------------------------------*/
  /// Wrapper for Sound propagation material
  class ElectromagneticMat : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    ElectromagneticMat();

    /// construct the material object given material parameters
    explicit ElectromagneticMat(Mat::PAR::ElectromagneticMat* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of parobject.H (this file) and should return it in this method.
     */
    int unique_par_object_id() const override
    {
      return ElectromagneticMatType::instance().unique_par_object_id();
    }

    /*!
        \brief Pack this class so it can be communicated

        Resizes the vector data and stores all information of a class in it.
        The first information to be stored in data has to be the
        unique parobject id delivered by unique_par_object_id() which will then
        identify the exact class on the receiving processor.

        \param data (in/out): char vector to store class information
     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
        \brief Unpack data from a char vector into this class

        The vector data contains all information to rebuild the
        exact copy of an instance of a class on a different processor.
        The first entry in data has to be an integer which is the unique
        parobject id defined at the top of this file and delivered by
        unique_par_object_id().

        \param data (in) : vector storing all data to be unpacked into this
        instance.
     */
    void unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_electromagneticmat;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> clone() const override
    {
      return Teuchos::rcp(new ElectromagneticMat(*this));
    }

    /// conductivity
    double sigma(int eleid = -1) const { return params_->get_parameter(params_->sigma_, eleid); }

    /// permittivity coefficient
    double epsilon(int eleid = -1) const
    {
      return params_->get_parameter(params_->epsilon_, eleid);
    }

    /// permeability coefficient
    double mu(int eleid = -1) const { return params_->get_parameter(params_->mu_, eleid); }



    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    //@}


   private:
    /// my material parameters
    Mat::PAR::ElectromagneticMat* params_;
  };

}  // namespace Mat



FOUR_C_NAMESPACE_CLOSE

#endif
