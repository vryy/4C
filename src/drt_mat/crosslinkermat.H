/*-----------------------------------------------------------*/
/*!
\file crosslinkermat.H

\brief A class for a crosslinker material

\maintainer Jonas Eichinger

\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef CROSSLINKER_H_
#define CROSSLINKER_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Teuchos_ParameterList.hpp>

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for de St. Venant--Kirchhoff
    class CrosslinkerMat : public Parameter
    {
    public:

      /// standard constructor
      CrosslinkerMat(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CrosslinkerMat() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
//      const double rlink_;
//      const double deltarlink_;
//      const double alink_;
//      const double iplink;


      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class CrosslinkerMat

  }  // namespace PAR

class CrosslinkerMatType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "CrosslinkerMatType"; }

  static CrosslinkerMatType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static CrosslinkerMatType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for Crosslinker Material
  class CrosslinkerMat : public Material
  {
  public:

    /// construct empty material object
    CrosslinkerMat();

    /// construct the material object given material parameters
    explicit CrosslinkerMat(MAT::PAR::CrosslinkerMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return CrosslinkerMatType::Instance().UniqueParObjectId(); }

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

    //! @name Access methods

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_crosslinkermat; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new CrosslinkerMat( *this ) ); }

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    /// density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->poissonratio_); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}

  private:

    /// my material parameters
    MAT::PAR::CrosslinkerMat* params_;

  };
}





#endif /* CROSSLINKER_H_ */
