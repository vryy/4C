/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_multiphase_reactions.H

 \brief a fluid material for porous multiphase flow with reactions (mass sources and sinks)

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H_
#define SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H_



#include "matpar_parameter.H"
#include "material.H"
#include "fluidporo_multiphase.H"
#include "../drt_lib/drt_parobjectfactory.H"

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

      /// destructor
      virtual ~FluidPoroMultiPhaseReactions() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

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
    std::string Name() const { return "FluidPoroMultiPhaseReactions"; }

    static FluidPoroMultiPhaseReactionsType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
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
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_fluidporo_multiphase_reactions;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new FluidPoroMultiPhaseReactions(*this));
    }

    /// number of reactions
    int NumReac() const { return paramsreac_->numreac_; }

    /// reaction ID by Index
    int ReacID(const unsigned index) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::FluidPoroMultiPhaseReactions* Parameter() const { return paramsreac_; }

    /// return whether reaction terms need to be evaluated
    virtual bool IsReactive() const { return true; };

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


#endif /* SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_REACTIONS_H_ */
