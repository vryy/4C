/*----------------------------------------------------------------------*/
/*! \file
 \brief a fluid material for porous multiphase flow with reactions (mass sources and sinks)

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H
#define BACI_MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H



#include "baci_lib_parobjectfactory.H"
#include "baci_mat_fluidporo_multiphase.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class FluidPoroMultiPhaseReactions : public FluidPoroMultiPhase
    {
     public:
      /// standard constructor
      FluidPoroMultiPhaseReactions(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// @name material parameters
      //@{

      /// provide ids of the individual reaction materials
      const std::vector<int>* ReacIds() const { return reacids_; }

      /// length of reaction list
      const int numreac_;

      /// the list of reaction IDs
      const std::vector<int>* reacids_;

      //@}

    };  // class FluidPoroMultiPhaseReactions

  }  // namespace PAR

  class FluidPoroMultiPhaseReactionsType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "FluidPoroMultiPhaseReactions"; }

    static FluidPoroMultiPhaseReactionsType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FluidPoroMultiPhaseReactionsType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class FluidPoroMultiPhaseReactions : public FluidPoroMultiPhase
  {
   public:
    /// construct empty material object
    FluidPoroMultiPhaseReactions();

    /// construct the material object given material parameters
    explicit FluidPoroMultiPhaseReactions(MAT::PAR::FluidPoroMultiPhaseReactions* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return FluidPoroMultiPhaseReactionsType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_fluidporo_multiphase_reactions;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new FluidPoroMultiPhaseReactions(*this));
    }

    /// number of reactions
    int NumReac() const { return paramsreac_->numreac_; }

    /// reaction ID by Index
    int ReacID(const unsigned index) const;

    /// Return quick accessible material parameter data
    MAT::PAR::FluidPoroMultiPhaseReactions* Parameter() const override { return paramsreac_; }

    /// return whether reaction terms need to be evaluated
    bool IsReactive() const override { return true; };

   protected:
    /// setup of material map
    virtual void SetupMatMap();

   private:
    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::FluidPoroMultiPhaseReactions* paramsreac_;
  };

}  // namespace MAT


#endif  // MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H
