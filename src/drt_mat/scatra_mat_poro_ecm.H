/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_poro_ecm.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_MAT_PORO_ECM_H_
#define SCATRA_MAT_PORO_ECM_H_


#include "scatra_mat.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraMatPoroECM : public ScatraMat
    {
    public:

      /// standard constructor
      ScatraMatPoroECM(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ScatraMatPoroECM() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Scatra

  } // namespace PAR

class ScatraMatPoroECMType : public ScatraMatType
{
public:

  std::string Name() const { return "ScatraMatPoroECMType"; }

  static ScatraMatPoroECMType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ScatraMatPoroECMType instance_;
};

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMatPoroECM : public ScatraMat
  {
  public:

    /// construct empty material object
    ScatraMatPoroECM();

    /// construct the material object given material parameters
    explicit ScatraMatPoroECM(MAT::PAR::ScatraMatPoroECM* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ScatraMatPoroECMType::Instance().UniqueParObjectId(); }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra_poroECM; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ScatraMatPoroECM( *this ) ); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ScatraMatPoroECM* params_;

  };

}



#endif /* SCATRA_MAT_PORO_ECM_H_ */
