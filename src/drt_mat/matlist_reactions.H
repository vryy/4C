/*!----------------------------------------------------------------------
 \brief This file contains the material for reactive scalars. It derives from MAT_matlist
and adds everything to supervise all the MAT_scatra_raction materials. The reactions
itself are defined inside the MAT_scatra_raction materials. So MAT_matlist_reactions
is just a "control instance".

\level 2

\maintainer Johannes Kremheller
*----------------------------------------------------------------------*/

#ifndef MatListReactions_H
#define MatListReactions_H



#include "matpar_parameter.H"
#include "material.H"
#include "matlist.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatListReactions : public virtual MatList
    {
     public:
      /// standard constructor
      MatListReactions(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~MatListReactions() { ; }

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

    };  // class MatListReactions

  }  // namespace PAR

  class MatListReactionsType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MatListReactionsType"; }

    static MatListReactionsType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MatListReactionsType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatListReactions : public virtual MatList
  {
   public:
    /// construct empty material object
    MatListReactions();

    /// construct the material object given material parameters
    explicit MatListReactions(MAT::PAR::MatListReactions* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return MatListReactionsType::Instance().UniqueParObjectId();
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
    virtual void Initialize();

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_matlist_reactions;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new MatListReactions(*this));
    }

    /// number of reactions
    int NumReac() const { return paramsreac_->numreac_; }

    /// reaction ID by Index
    int ReacID(const unsigned index) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::MatListReactions* Parameter() const { return paramsreac_; }

    /// advanced reaction terms
    double CalcReaBodyForceTerm(const int k, const std::vector<double>& phinp,
        const double* gpcoord,  //!< current Gauss-point coordinates
        const double scale = 1.0) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k, std::vector<double>& derivs,
        const std::vector<double>& phinp,
        const double* gpcoord,  //!< current Gauss-point coordinates
        const double scale = 1.0) const;

    /// advanced reaction terms
    double CalcReaBodyForceTerm(const int k, const std::vector<double>& phinp,
        const std::vector<std::pair<std::string, double>>& constants,
        const double* gpcoord,  //!< current Gauss-point coordinates
        const double scale = 1.0) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k, std::vector<double>& derivs,
        const std::vector<double>& phinp,
        const std::vector<std::pair<std::string, double>>& constants,
        const double* gpcoord,  //!< current Gauss-point coordinates
        const double scale = 1.0) const;

    // add additional variables to the reaction (only for by-function coupling)
    void AddAdditionalVariables(
        const int k, const std::vector<std::pair<std::string, double>>& variables) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrixAddVariables(const int k, std::vector<double>& derivs,
        const std::vector<double>& phinp,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants, const double* gpcoord,
        const double scale = 1.0) const;


   protected:
    /// setup of material map
    virtual void SetupMatMap();

   private:
    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatListReactions* paramsreac_;
  };

}  // namespace MAT

#endif
