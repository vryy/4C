/*----------------------------------------------------------------------*/
/*!
\file logneohooke.H
\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

*----------------------------------------------------------------------*/
#ifndef LOGNEOHOOKE_H_
#define LOGNEOHOOKE_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for logarithmic neo-Hooke material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 MAT_Struct_LogNeoHooke YOUNG 1.044E7 NUE 0.3 DENS 1.0
    class LogNeoHooke : public Parameter
    {
    public:

      /// standard constructor
      LogNeoHooke(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~LogNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double nue_;
      /// mass density
      const double density_;
      /// model type
      const int model_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class LogNeoHooke

  }  // namespace PAR

class LogNeoHookeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "LogNeoHookeType"; }

  static LogNeoHookeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static LogNeoHookeType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Logarithmic neo-Hooke material according to [1]
  ///
  /// This is a compressible, hyperelastic, isotropic material
  /// of the most simple kind.
  ///
  /// You can easily get a purely incompressible (or isochoric) response
  /// by setting Poisson's ration to 0.5. However be warned, this will only
  /// work with an elment capable of purely incompressible response without
  /// addressing a volumetric strain-stress relation.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] J Bonet and RD Wood, "Nonlinear continuum mechanics for finite
  ///          element analysis", Cambridge, 1997.
  /// <li> [2] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author bborn
  /// \date 04/09
  class LogNeoHooke : public Material
  {

  public:

    /// empty constructor
    LogNeoHooke();

    /// constructor with given material parameters
    LogNeoHooke(MAT::PAR::LogNeoHooke* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return LogNeoHookeType::Instance().UniqueParObjectId(); }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    virtual void Pack(
      DRT::PackBuffer& data  ///< (in/out): char vector to store class information
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

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->nue_); }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_logneohooke; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new LogNeoHooke( *this ) ); }

    //@}

    /// Main <i>0-type</i> material call to determine  2nd PK stress and material constitutive tensor
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = \frac{\mu}{2} (I_{\boldsymbol{C}} - 3)
    ///       - \mu \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
    ///      + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
    ///\f]
    /// which was taken from [1].
    ///
    /// Here is \f$I_{\boldsymbol{C}}\f$ the first principal invariant
    /// of the right Cauchy--Green strain tensor \f$\boldsymbol{C}\f$
    /// and \f$I\!I\!I_{\boldsymbol{C}}\f$ its third.
    /// The isochoric part is proportional to \f$\mu\f$,
    /// whereas the \f$\lambda\f$-proportional part constitutes the volumetric contribution.
    ///
    /// The parameters are the shear modulus
    ///\f[
    /// \mu = \frac{E}{2(1+\nu)}
    ///\f]
    /// and Lame's coefficient
    ///\f[
    /// \lambda = \left\{\begin{array}{ll}
    ///      \frac{E \nu}{(1+\nu) (1-2\nu)} & \text{if $\nu \neq 1/2$}
    ///   \\ 0                              & \text{else}
    /// \end{array}\right.
    ///\f]
    ///
    /// The 2nd Piola--Kirchhoff stress is
    ///\f[
    /// \boldsymbol{S} = \mu \big( \boldsymbol{1} - \boldsymbol{C}^{-1} \big)
    ///                + \lambda \, \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \,  \boldsymbol{C}^{-1}
    ///\f]
    ///
    /// The material constitutive 4-tensor
    /// \f$\boldsymbol{C}_\text{m}=C_{IJKL}\boldsymbol{E}^I\otimes\boldsymbol{E}^J\otimes\boldsymbol{E}^K\otimes\boldsymbol{E}^L\f$
    /// is determined by
    ///\f[
    /// C_{IJKL}
    /// = \lambda (\boldsymbol{C}^{-1})_{IJ} \, (\boldsymbol{C}^{-1})_{KL}
    /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
    ///   \, \frac{1}{2} \big( \boldsymbol{C}^{-1})_{IK} \, (\boldsymbol{C}^{-1})_{JL}
    ///                        + (\boldsymbol{C}^{-1})_{IL} \, (\boldsymbol{C}^{-1})_{JK} \big)
    ///\f]
    /// or
    ///\f[
    /// \boldsymbol{C}_\text{m}
    /// = \lambda \boldsymbol{C}^{-1} \otimes  \boldsymbol{C}^{-1}
    /// + 2\big(\mu-\lambda \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})\big)
    ///   \, \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
    ///\f]
    void EvaluateModelZero(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain tensor
      LINALG::Matrix<6,6>& cmat,  ///< constitutive matrix
      LINALG::Matrix<6,1>& stress  ///< 2nd PK stress
      );

    /// Main <i>1-type</i> material call to determine 2nd PK stress and material constitutive tensor
    ///
    /// The material strain energy density function is
    ///\f[
    /// \Psi = \frac{\mu_1}{2} \Big(
    ///         \frac{I_{\boldsymbol{C}}}{\sqrt[3]{I\!I\!I_{\boldsymbol{C}}}} - 3
    ///      \Big)
    ///      - \mu_3 \log(\sqrt{I\!I\!I_{\boldsymbol{C}}})
    ///      + \frac{\lambda}{2} \big( \log(\sqrt{I\!I\!I_{\boldsymbol{C}}}) \big)^2
    ///\f]
    /// which was taken from [1] and adapted to a clean volumetric-isochoric decomposition.
    ///
    /// Here \f$I_{\boldsymbol{C}}\f$ is the first principal invariant
    /// of the right Cauchy--Green strain tensor \f$\boldsymbol{C}\f$
    /// and \f$I\!I\!I_{\boldsymbol{C}}\f$ its third.
    /// The isochoric part is proportional to \f$\mu_1\f$,
    /// whereas the \f$\mu_3,\lambda\f$-proportional parts constitute the volumetric contribution.
    ///
    /// The parameters are the '1st shear modulus'
    ///\f[
    /// \mu_1 = \frac{E}{2(1+\nu)}
    ///\f]
    /// and '3rd shear modulus'
    ///\f[
    /// \mu_3 = \left\{\begin{array}{ll}
    ///      \mu_1 & \text{if $\nu \neq 1/2$}
    ///   \\ 0     & \text{else}
    /// \end{array}\right.
    ///\f]
    /// and Lame's coefficient
    ///\f[
    /// \lambda = \left\{\begin{array}{ll}
    ///      \frac{E \nu}{(1+\nu) (1-2\nu)} & \text{if $\nu \neq 1/2$}
    ///   \\ 0                              & \text{else}
    /// \end{array}\right.
    ///\f]
    ///
    /// The 2nd Piola--Kirchhoff stress is
    ///\f[
    /// \boldsymbol{S} = \frac{\mu_1}{\sqrt[3]{I\!I\!I_{\boldsymbol{C}}}} \Big(
    ///                  \boldsymbol{1} - \frac{I_{\boldsymbol{C}}}{3} \boldsymbol{C}^{-1}
    ///                \Big)
    ///                - \Big(
    ///                  \mu_3 - \frac{\lambda}{2} \, \log(I\!I\!I_{\boldsymbol{C}})
    ///                \Big) \,  \boldsymbol{C}^{-1}
    ///\f]
    ///
    /// The material constitutive 4-tensor
    ///\f[
    /// \boldsymbol{C}_\text{m}
    /// = -\frac{2 \mu_1 I_{\boldsymbol{C}}}{3 \sqrt[3]{I\!I\!I_{\boldsymbol{C}}}} \Big(
    ///    \boldsymbol{1} \otimes \boldsymbol{C}^{-1}
    ///    + \boldsymbol{C}^{-1} \otimes \boldsymbol{I}
    /// \Big)
    /// + \Big( \frac{2 \mu_1 I_{\boldsymbol{C}}}{9 \sqrt[3]{I\!I\!I_{\boldsymbol{C}}}} + \lambda \Big)
    ///    \, \boldsymbol{C}^{-1} \otimes  \boldsymbol{C}^{-1}
    /// + \Big( \frac{2 \mu_1 I_{\boldsymbol{C}}}{3 \sqrt[3]{I\!I\!I_{\boldsymbol{C}}}}
    ///    + 2 \mu_3 - \lambda\log(I\!I\!I_{\boldsymbol{C}}) \Big)
    ///    \, \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
    ///\f]
    void EvaluateModelOne(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain tensor
      LINALG::Matrix<6,6>& cmat,  ///< constitutive matrix
      LINALG::Matrix<6,1>& stress  ///< 2nd PK stress
      );

    /// Main material call to determine 2nd PK stress and material constitutive tensor
    ///
    /// The method branches according to material sub-model #model_.
    void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain tensor
      LINALG::Matrix<6,6>& cmat,  ///< constitutive matrix
      LINALG::Matrix<6,1>& stress  ///< 2nd PK stress
      );

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::LogNeoHooke* params_;

  };
}

#endif  // LOGNEOHOOKE_H_
