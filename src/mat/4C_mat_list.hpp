/*----------------------------------------------------------------------------*/
/*! \file
\brief generic material that stores a list of materials, where each material itself defines the
properties of e.g. one species in a scalar transport problem, or one phase in a fluid problem

\level 1


*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_LIST_HPP
#define FOUR_C_MAT_LIST_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatList : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      MatList(Teuchos::RCP<Core::Mat::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// @name material parameters
      //@{

      /// provide ids of the individual materials
      const std::vector<int>* MatIds() const { return &matids_; }

      /// provide access to material by its ID
      Teuchos::RCP<Core::Mat::Material> MaterialById(const int id) const;

      std::map<int, Teuchos::RCP<Core::Mat::Material>>* material_map_write() { return &mat_; }

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int> matids_;

      /// flag for individual materials or only one at global scope
      bool local_;

     private:
      /// map to materials (only used for local_==true)
      std::map<int, Teuchos::RCP<Core::Mat::Material>> mat_;

      //@}

    };  // class MatList

  }  // namespace PAR

  class MatListType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "MatListType"; }

    static MatListType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MatListType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatList : public Core::Mat::Material
  {
   public:
    /// construct empty material object
    MatList();

    /// construct the material object given material parameters
    explicit MatList(Mat::PAR::MatList* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return MatListType::Instance().UniqueParObjectId(); }

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
      return Core::Materials::m_matlist;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new MatList(*this));
    }

    /// number of materials
    int NumMat() const { return params_->nummat_; }

    /// material ID by Index
    int MatID(const unsigned index) const;

    /// provide access to material by its ID
    virtual Teuchos::RCP<Core::Mat::Material> MaterialById(const int id) const;

    /// Return quick accessible material parameter data
    Mat::PAR::MatList* Parameter() const override { return params_; }

   protected:
    /// return pointer to the materials map, which has read-only access.
    const std::map<int, Teuchos::RCP<Core::Mat::Material>>* material_map_read() const
    {
      return &mat_;
    }

    /// return pointer to the materials map, which has read and write access.
    std::map<int, Teuchos::RCP<Core::Mat::Material>>* material_map_write() { return &mat_; }

   private:
    /// setup of material map
    void setup_mat_map();

    /// clear everything
    void clear();

    /// my material parameters
    Mat::PAR::MatList* params_;

    /// map to materials
    std::map<int, Teuchos::RCP<Core::Mat::Material>> mat_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
