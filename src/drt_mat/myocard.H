/*----------------------------------------------------------------------*/
/*!
\file myocard.H
\brief

<pre>
Maintainer:Cristobal Bertoglio
            bertoglio@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_H
#define MYOCARD_H

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 08/13 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "myocard_general.H"

/*----------------------------------------------------------------------*
 |                                                           cbert 08/13 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class Myocard : public Parameter
    {
    public:
      /// standard constructor
      Myocard( Teuchos::RCP<MAT::PAR::Material> matdata );

      /// destructor
      virtual ~Myocard() { ; }

      /// @name material parameters
      //@{

      /// Diffusivity
      const double maindirdiffusivity;
      const double offdirdiffusivity;

      /// Pertubation for calculation of derivative
      const double dt_deriv;

      /// Model type. Possibilities are: "MV", "FHN", "INADA", "TNNP", "SAN"
      const std::string* model;

      /// Tissue type. Possibilities are: "M", "ENDO", "EPI" for "TNNP"
      ///                                 "AN", "N", "NH" for "INADA"
      const std::string* tissue;

      /// Time factor to correct for different Model specific time units
      const double time_scale;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    }; // class myocard
  } // namespace PAR

class MyocardType : public DRT::ParObjectType
{
public:
  std::string Name() const { return "MyocardType"; }
  static MyocardType & Instance() { return instance_; };
  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:
  static MyocardType instance_;
};

/*----------------------------------------------------------------------*/
/// Myocard material
///
/// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac muscle tissue
///
/// \author cbert
///
/// \date 08/13

class Myocard : public Material

{

public:
  /// construct empty material object
  Myocard();

  /// constructor with given material parameters
  Myocard(MAT::PAR::Myocard* params);

  /// @name Packing and Unpacking
  //@{

  /// Return unique ParObject id
  ///
  ///  every class implementing ParObject needs a unique id defined at the
  ///  top of drt_parobject.H (this file) and should return it in this method.
  virtual int UniqueParObjectId() const { return MyocardType::Instance().UniqueParObjectId(); }

  /// Pack this class so it can be communicated
  ///
  /// Resizes the vector data and stores all information of a class in it.
  /// The first information to be stored in data has to be the
  /// unique parobject id delivered by UniqueParObjectId() which will then
  /// identify the exact class on the receiving processor.
  virtual void Pack(DRT::PackBuffer& data) const;  ///< (in/out): char vector to store class information

  /// \brief Unpack data from a char vector into this class
  ///
  /// The vector data contains all information to rebuild the
  /// exact copy of an instance of a class on a different processor.
  /// The first entry in data has to be an integer which is the unique
  /// parobject id defined at the top of this file and delivered by
  /// UniqueParObjectId().
  ///
  virtual void Unpack(const std::vector<char>& data);  ///< vector storing all data to be unpacked into this

  //@}

  /// material type
  virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_myocard; }

  /// return copy of this material object
  virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Myocard( *this ) ); }

  /// material call from DRT::ELEMENTS::Transport::ReadElement function
  /// to setup conductivity tensor for each element
  void Setup(const std::vector<double> &fiber1);
  void Setup(DRT::INPUT::LineDefinition* linedef);
  void SetupDiffusionTensor(const std::vector<double> &fiber1);

  /// diffusivity
  void ComputeDiffusivity(LINALG::Matrix<1,1>& diffus3) const;
  void ComputeDiffusivity(LINALG::Matrix<2,2>& diffus3) const;
  void ComputeDiffusivity(LINALG::Matrix<3,3>& diffus3) const;

  /// compute reaction coefficient
  double ComputeReactionCoeff(const double phi, const double dt);

  /// compute reaction coefficient derivative
  double ComputeReactionCoeffDeriv(const double csnp, const double dt);

  /// compute Heaviside step function
  double GatingFunction(const double Gate1, const double Gate2, const double p, const double var, const double thresh) const;

  /// compute gating variable 'y' from dy/dt = (y_inf-y)/y_tau
  double GatingVarCalc(const double dt, double y_0, const double y_inf, const double y_tau) const;

  ///  returns number of internal state variables of the material
  int GetNumberOfInternalStateVariables() const;

  ///  returns current internal state of the material
  double GetInternalState(const int k) const;

  ///  set internal state of the material
  void SetInternalState(const int k, const double val);

  /// initialize internal variables (called by constructors)
  void Initialize();

  /// time update for this material
  void Update(const double phi, const double dt);

  /// @name Access material constants
  //@{

  //@}

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const { return params_; }

private:

  /// my material parameters
  MAT::PAR::Myocard* params_;

  /// conductivity tensor
  LINALG::Matrix<3,3> difftensor_;

  /// number of internal state variables
  int nb_state_variables_;

  // Type of material model
  Teuchos::RCP<Myocard_General> myocard_mat_;

}; // Myocard
} // namespace MAT

#endif
