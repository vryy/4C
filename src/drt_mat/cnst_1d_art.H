/*----------------------------------------------------------------------*/
/*!
\file cnst_1d_art.H

\level 3

\maintainer Lena Yoshihara
*/
/*----------------------------------------------------------------------*/
#ifndef CNST_1D_ART_H
#define CNST_1D_ART_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for constant 1D_Artery
    ///
    // This object exists only once for each read Newton fluid. ???
    class Cnst_1d_art : public Parameter
    {
    public:

      /// standard constructor
      Cnst_1d_art(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Cnst_1d_art() { }

      /// @name material parameters
      //@{
      /// Newtonian viscosity of blood
      const double viscosity_;
      /// density of blood
      const double density_;
      /// Artery Youngs modulus of elasticity */
      const double young_;
      /// Artery Poisson's ratio
      const double nue_;
      /// Artery initial diameter
      const double diam_;
      /// Artery wall thickness
      const double th_;
      /// Fixed external pressure at node 1
      const double pext1_;
      /// Fixed external pressure at node 2
      const double pext2_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class NewtonianFluid

  }  // namespace PAR

class Cnst_1d_artType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "Cnst_1d_artType"; }

  static Cnst_1d_artType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static Cnst_1d_artType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for constant 1D_Artery material
  ///
  /// This object exists (several times) at every element
  class Cnst_1d_art : public Material
  {
  public:

    /// construct empty material object
    Cnst_1d_art();

    /// construct the material object given material parameters
    explicit Cnst_1d_art(MAT::PAR::Cnst_1d_art* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return Cnst_1d_artType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_cnst_art; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Cnst_1d_art( *this ) ); }

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return density
    virtual double Density() const { return params_->density_; }

    /// return Youngs modulus
    double Young() const { return params_->young_; }

    /// return Poisson's ratio
    double Nue() const { return params_->nue_; }

    /// return artery diameter
    double Diam() const { return params_->diam_; }

    /// return artery wall thickness
    double Th() const { return params_->th_; }

    /// return artery external pressure
    double pext(int i) const
           {if(i==0)
              return params_->pext1_;
            else if(i==1)
              return params_->pext2_;
            else dserror("There is no pressure with id %d",i);
              return 0.0;
           }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::Cnst_1d_art* params_;
  };

}

#endif

