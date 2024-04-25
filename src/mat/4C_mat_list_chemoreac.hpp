/*----------------------------------------------------------------------*/
/*! \file
 \brief

This file contains the material for reactive AND chemotactic scalars. It is
in diamond inheritance with MatListReactions and MatListChemotaxis,
which govern the actual doings

\level 3
*----------------------------------------------------------------------*/


#ifndef FOUR_C_MAT_LIST_CHEMOREAC_HPP
#define FOUR_C_MAT_LIST_CHEMOREAC_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_chemotaxis.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_material.hpp"
#include "4C_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatListChemoReac : public MatListReactions, public MatListChemotaxis
    {
     public:
      /// standard constructor
      MatListChemoReac(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// @name material parameters

    };  // class MatListReactions

  }  // namespace PAR

  class MatListChemoReacType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "MatListChemoReacType"; }

    static MatListChemoReacType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MatListChemoReacType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatListChemoReac : public MatListChemotaxis, public MatListReactions
  {
   public:
    /// construct empty material object
    MatListChemoReac();

    /// construct the material object given material parameters
    explicit MatListChemoReac(MAT::PAR::MatListChemoReac* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return MatListChemoReacType::Instance().UniqueParObjectId();
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

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_matlist_chemoreac;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new MatListChemoReac(*this));
    }

    /// Return quick accessible material parameter data
    MAT::PAR::MatListChemoReac* Parameter() const override { return paramsreachemo_; }

   private:
    /// setup of material map
    void SetupMatMap() override;

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatListChemoReac* paramsreachemo_;
  };

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif
