/*!----------------------------------------------------------------------
\file scatra_reaction_mat.H

\brief This file contains the base material for reactive scalars.

\level 2
<pre>
\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/


#ifndef SCATRA_REACTION_MAT_H
#define SCATRA_REACTION_MAT_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  //forward declaration
  class MatListReactions;

  namespace PAR
  {

  enum reaction_coupling
  {
    reac_coup_none,  ///< no coupling, initialization value
    reac_coup_simple_multiplicative, ///< coupling of type A*B
    reac_coup_power_multiplicative, ///< coupling of type A*B
    reac_coup_constant, //< no coupling, constant increase/decrease
    reac_coup_michaelis_menten, //coupling of type (B/(const+B))*A
  };

  /*----------------------------------------------------------------------*/
  /// parameters for scalar transport material
  class ScatraReactionMat : public Parameter
  {
  public:

    /// standard constructor
    ScatraReactionMat(
      Teuchos::RCP<MAT::PAR::Material> matdata
      );

    /// destructor
    virtual ~ScatraReactionMat() { ; }

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    void OptParams(std::map<std::string, int>* pnames);

    /// returns the enum of the current coupling type
    MAT::PAR::reaction_coupling SetCouplingType( Teuchos::RCP<MAT::PAR::Material> matdata );

    /// number of scalars in this reaction
    const int numscal_;

    /// the list of material IDs
    const std::vector<int>* stoich_;

    /// reaction coefficient
    const double reaccoeff_;

    /// type of coupling
    const MAT::PAR::reaction_coupling coupling_;

    /// specifies scalar type in reaction
    const std::vector<double>* couprole_;

    /// parameter to define start of reaction
    const double reacstart_;

  };  // class Scatra

} // namespace PAR

class ScatraReactionMatType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ScatraReactionMatType"; }

  static ScatraReactionMatType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ScatraReactionMatType instance_;
};

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraReactionMat : public Material
  {
  friend class MAT::MatListReactions;

  public:

    /// construct empty material object
    ScatraReactionMat();

    /// construct the material object given material parameters
    explicit ScatraReactionMat(MAT::PAR::ScatraReactionMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ScatraReactionMatType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_reaction; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ScatraReactionMat( *this ) ); }

    /// return number of scalars for this reaction
    int NumScal() const { return params_->numscal_; }

  protected:

    /// return reaction coefficient
    virtual double ReacCoeff() const { return params_->reaccoeff_; }

    /// return stoichometrie
    std::vector<int>* Stoich() const { return (std::vector<int>*)params_->stoich_; }

    /// return type of coupling
    MAT::PAR::reaction_coupling Coupling() const { return params_->coupling_; }

    /// return role in coupling
    std::vector<double>* Couprole() const { return (std::vector<double>*)params_->couprole_;}

    /// delayed reaction start coefficient
    double ReacStart() const { return params_->reacstart_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Calculate K(c)
    double CalcReaCoeff(
        const int k,
        const std::vector<double>& phinp,
        const double scale
        ) const;

    /// calculate \frac{partial}{\partial c} K(c)
    double CalcReaCoeffDerivMatrix(
        const int k,
        const int j,
        const std::vector<double>& phinp,
        const double scale
        ) const;

    /// Calculate f(c)
    double CalcReaBodyForceTerm(
        const int k,
        const std::vector<double>& phinp,
        const double scale
        ) const;

    /// calculate \frac{partial}{\partial c} f(c)
    double CalcReaBodyForceDerivMatrix(
        const int k,
        const int j,
        const std::vector<double>& phinp,
        const double scale
        ) const;

  private:

    ///helper for calculating K(c)
    double CalcReaCoeffFac(
        const int k,
        const std::vector<double>& phinp,
        const double reacstart,              //!<reaction start coefficient
        const double scale
        ) const;

    ///helper for calculating \frac{partial}{\partial c} K(c)
    double CalcReaCoeffDerivFac(
        const int k,
        const int toderive,
        const std::vector<double>& phinp,
        const double reacstart,              //!<reaction start coefficient
        const double scale
        ) const;

    ///helper for calculating f(c)
    double CalcReaBodyForceTermFac(
        const int k,
        const std::vector<double>& phinp,
        const double reacstart,              //!<reaction start coefficient
        const double scale ) const;

    ///helper for calculating calculate \frac{partial}{\partial c} f(c)
    double CalcReaBodyForceDerivFac(
        const int k, const int toderive,
        const std::vector<double>& phinp,
        const double reacstart,              //!<reaction start coefficient
        const double scale ) const;

  protected:

    /// my material parameters
    MAT::PAR::ScatraReactionMat* params_;

  };

}

#endif
