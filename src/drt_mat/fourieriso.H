/*----------------------------------------------------------------------*/
/*!
\brief material for heat conduction according to fourier's law
\level 1
<pre>
\maintainer Christoph Meier
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/
#ifndef FOURIERISO_H_
#define FOURIERISO_H_

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"

#include "../drt_lib/drt_parobjectfactory.H"

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for FourierIso material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 THERM_FourierIsoIso CAPA 1.0 COND 1.0
    class FourierIso : public Parameter
    {
     public:
      /// standard constructor
      FourierIso(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FourierIso() { ; }

      /// @name material parameters
      //@{

      /// volumetric heat capacity
      const double capa_;
      /// heat conductivity
      const double conduct_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class FourierIso

  }  // namespace PAR

  class FourierIsoType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "FourierIsoType"; }

    static FourierIsoType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static FourierIsoType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// FourierIso material according to [1]
  ///
  /// This is a FourierIso's law of isotropic, instationary heat conduction
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author dano
  /// \date 09/09
  class FourierIso : public Material
  {
   public:
    /// empty constructor
    FourierIso();

    /// constructor with given material parameters
    FourierIso(MAT::PAR::FourierIso* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return FourierIsoType::Instance().UniqueParObjectId(); }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    virtual void Pack(DRT::PackBuffer& data  ///< (in/out): char vector to store class information
        ) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    virtual void Unpack(
        const std::vector<char>& data  ///< vector storing all data to be unpacked into this
    );

    //@}

    /// @name Access material constants
    //@{

    /// conductivity
    double Conductivity() const { return params_->conduct_; }

    /// volumetric heat capacity
    double Capacity() const { return params_->capa_; }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_th_fourier_iso; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new FourierIso(*this)); }

    //@}

    /// Main 1D material call to determine heat flux and constitutive tensor
    void Evaluate(const LINALG::Matrix<1, 1>& gradtemp,  ///< temperature gradient (strain tensor)
        LINALG::Matrix<1, 1>& cmat,                      ///< constitutive matrix
        LINALG::Matrix<1, 1>& heatflux                   ///< 3D heatflux (2nd PK stress)
        ) const;

    /// Main 2D material call to determine heat flux and constitutive tensor
    void Evaluate(const LINALG::Matrix<2, 1>& gradtemp,  ///< temperature gradient (strain tensor)
        LINALG::Matrix<2, 2>& cmat,                      ///< constitutive matrix
        LINALG::Matrix<2, 1>& heatflux                   ///< 3D heatflux (2nd PK stress)
        ) const;

    /// Main 3D material call to determine heat flux and constitutive tensor
    void Evaluate(const LINALG::Matrix<3, 1>& gradtemp,  ///< temperature gradient (strain tensor)
        LINALG::Matrix<3, 3>& cmat,                      ///< constitutive matrix
        LINALG::Matrix<3, 1>& heatflux                   ///< 3D heatflux (2nd PK stress)
        ) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::FourierIso* params_;

  };  // FourierIso

}  // namespace MAT

/*----------------------------------------------------------------------*/
#endif  // FOURIERISO_H_
