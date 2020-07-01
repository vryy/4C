/*----------------------------------------------------------------------*/
/*! \file
 \brief a fluid material for porous multiphase flow defining one reaction (mass sources and sinks)

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_SINGLEREACTION_H_
#define SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_SINGLEREACTION_H_


#include "fluidporo_singlephase.H"

/*---------------------------------------------------------------------*
 | forward declarations                                                |
 *---------------------------------------------------------------------*/
namespace DRT
{
  namespace UTILS
  {
    class VariableExprFunction;
  }
}  // namespace DRT

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for a single phase of porous multiphase fluid
    ///
    /// This object exists only once for each read fluid.
    class FluidPoroSingleReaction : public Parameter
    {
      enum porofluid_reaction_coupling
      {
        porofluid_reac_coup_none,              ///< no coupling, initialization value
        porofluid_reac_coup_scalarsbyfunction  ///< reaction depending on scalars defined by
                                               ///< function
      };

     public:
      /// standard constructor
      FluidPoroSingleReaction(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FluidPoroSingleReaction() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// initialize
      void Initialize();

      /// evaluate reaction for by-function definition
      void EvaluateFunction(std::vector<double>& reacval,
          std::vector<std::vector<double>>& reacderivspressure,
          std::vector<std::vector<double>>& reacderivssaturation,
          std::vector<double>& reacderivsporosity,
          std::vector<std::vector<double>>& reacderivsvolfrac,
          std::vector<std::vector<double>>& reacderivsvolfracpressure,
          std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
          const std::vector<double>& saturation, const double& porosity,
          const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
          const std::vector<double>& scalar);

      /// Check sizes of vectors
      void CheckSizes(std::vector<double>& reacval,
          std::vector<std::vector<double>>& reacderivspressure,
          std::vector<std::vector<double>>& reacderivssaturation,
          std::vector<double>& reacderivsporosity,
          std::vector<std::vector<double>>& reacderivsvolfrac,
          std::vector<std::vector<double>>& reacderivsvolfracpressure,
          std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
          const std::vector<double>& saturation, const double& porosity,
          const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
          const std::vector<double>& scalar);

      /// @name material parameters
      //@{
      /// number of scalars in this problem (not only in this reaction but total problem!)
      const int numscal_;

      /// number of additional volume fractions of this problem
      const int numvolfrac_;

      /// total number of multiphase-dofs = number of fluid phases + number of additional volume
      /// fractions
      const int totalnummultiphasedof_;

      /// number of fluid phases in this problem
      const int numfluidphases_;

      /// the list of material IDs
      const std::vector<int>* scale_;

      /// type of coupling
      const MAT::PAR::FluidPoroSingleReaction::porofluid_reaction_coupling coupling_;

      /// ID of the function defining the reaction
      const int functID_;
      //@}

     private:
      /// returns the enum of the current coupling type
      MAT::PAR::FluidPoroSingleReaction::porofluid_reaction_coupling SetCouplingType(
          Teuchos::RCP<MAT::PAR::Material> matdata);

      /// return variable expression function holding the reaction equation
      inline DRT::UTILS::VariableExprFunction& Function(int functnum) const;

      /// flag if Initialize() has been called
      bool isinit_;

      /// string name used for scalars in function parser
      std::vector<std::string> scalarnames_;
      /// string name used for pressure in function parser
      std::vector<std::string> pressurenames_;
      /// string name used for saturation in function parser
      std::vector<std::string> saturationnames_;
      /// string name used for porosity in function parser
      const std::string porosityname_;
      /// string name used for volume fractions in function parser
      std::vector<std::string> volfracnames_;
      /// string name used for volume fraction pressures in function parser
      std::vector<std::string> volfracpressurenames_;

    };  // class FluidPoroSinglePhaseReaction

  }  // namespace PAR

  /*----------------------------------------------------------------------*
   | instance access method                                   vuong 08/16 |
   *----------------------------------------------------------------------*/
  class FluidPoroSingleReactionType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "FluidPoroSingleReactionType"; }

    static FluidPoroSingleReactionType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static FluidPoroSingleReactionType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a single porous fluid phase within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroSingleReaction : public FluidPoroSinglePhaseBase
  {
   public:
    /// construct empty material object
    FluidPoroSingleReaction();

    /// construct the material object given material parameters
    explicit FluidPoroSingleReaction(MAT::PAR::FluidPoroSingleReaction* params);

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_fluidporo_singlereaction;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new FluidPoroSingleReaction(*this));
    }

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return FluidPoroSinglePhaseType::Instance().UniqueParObjectId();
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

    /// initialize
    void Initialize();

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// evaluate reaction
    void EvaluateReaction(std::vector<double>& reacval,
        std::vector<std::vector<double>>& reacderivspressure,
        std::vector<std::vector<double>>& reacderivssaturation,
        std::vector<double>& reacderivsporosity,
        std::vector<std::vector<double>>& reacderivsvolfrac,
        std::vector<std::vector<double>>& reacderivsvolfracpressure,
        std::vector<std::vector<double>>& reacderivsscalar, const std::vector<double>& pressure,
        const std::vector<double>& saturation, const double& porosity,
        const std::vector<double>& volfracs, const std::vector<double>& volfracpressures,
        const std::vector<double>& scalar);

    /// return whether phase 'phasenum' is involved in this reaction
    bool IsReactive(int phasenum) const { return params_->scale_->at(phasenum) != 0; };

    /// Check sizes of vectors
    int TotalNumDof() { return params_->totalnummultiphasedof_; }

   private:
    /// my material parameters
    MAT::PAR::FluidPoroSingleReaction* params_;
  };

}  // namespace MAT

#endif /* SRC_DRT_MAT_FLUIDPORO_MULTIPHASE_SINGLEREACTION_H_ */
