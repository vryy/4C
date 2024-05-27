/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of ion species in electrolyte solutions

\level 2


*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ELCHPHASE_HPP
#define FOUR_C_MAT_ELCHPHASE_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class ElchPhase : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      ElchPhase(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> create_material() override;

      /// provide ids of the individual mat
      const std::vector<int>& MatIds() const { return matids_; }

      /// provide access to phases by its ID
      Teuchos::RCP<CORE::MAT::Material> MatById(const int id) const
      {
        if (not local_)
        {
          std::map<int, Teuchos::RCP<CORE::MAT::Material>>::const_iterator m = mat_.find(id);

          if (m == mat_.end())
          {
            FOUR_C_THROW("Material %d could not be found", id);
            return Teuchos::null;
          }
          else
            return m->second;
        }
        else
          FOUR_C_THROW("This is not allowed");

        return Teuchos::null;
      }

      /// @name material parameters
      //@{

      /// porosity
      const double epsilon_;

      /// tortuosity
      const double tortuosity_;

      /// number of materials
      const int nummat_;

      /// the list of material IDs
      const std::vector<int> matids_;

      /// flag for individual materials or only one at global scope
      bool local_;

      //@}

     private:
      /// map to materials (only used for local_==true)
      std::map<int, Teuchos::RCP<CORE::MAT::Material>> mat_;

    };  // class ElchPhase

  }  // namespace PAR

  class ElchPhaseType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ElchPhaseType"; }

    static ElchPhaseType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ElchPhaseType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class ElchPhase : public CORE::MAT::Material
  {
   public:
    /// construct empty material object
    ElchPhase();

    /// construct the material object given material parameters
    explicit ElchPhase(MAT::PAR::ElchPhase* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return ElchPhaseType::Instance().UniqueParObjectId(); }

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
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_elchphase;
    }

    /// return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new ElchPhase(*this));
    }

    /// return constant porosity
    double Epsilon() const { return params_->epsilon_; }
    /// return constant tortuosity
    double Tortuosity() const { return params_->tortuosity_; }

    int NumMat() const { return params_->nummat_; }

    /// material ID by Index
    int MatID(const unsigned index) const
    {
      if ((int)index < params_->nummat_)
        return params_->matids_.at(index);
      else
      {
        FOUR_C_THROW("Index too large");
        return -1;
      }
    }

    /// provide access to material by its ID
    Teuchos::RCP<CORE::MAT::Material> MatById(const int id) const
    {
      if (params_->local_)
      {
        std::map<int, Teuchos::RCP<CORE::MAT::Material>>::const_iterator m = mat_.find(id);
        if (m == mat_.end())
        {
          FOUR_C_THROW("Material %d could not be found", id);
          return Teuchos::null;
        }
        else
          return m->second;
      }
      else  // material is global (stored in material parameters)
        return params_->MatById(id);
    }

    /// Return quick accessible material parameter data
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// setup of material map
    void setup_mat_map();

    /// clear everything
    void clear();

    /// my material parameters
    MAT::PAR::ElchPhase* params_;

    /// map to materials
    std::map<int, Teuchos::RCP<CORE::MAT::Material>> mat_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
