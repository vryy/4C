/*----------------------------------------------------------------------*/
/*!
\file elast_aniso_structuraltensor_strategy.H

\brief strategy for evaluation of the structural tensor in anisotropic materials

\level 2

\maintainer Andreas Rauch

*/
/*----------------------------------------------------------------------*/
#ifndef ELAST_STRUCTURALTENSOR_STRATEGY_
#define ELAST_STRUCTURALTENSOR_STRATEGY_

#define numbgp 32
#define twice 64

#include "../drt_mat/matpar_parameter.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MAT
{
namespace ELASTIC
{
namespace PAR
{

typedef enum _DISTR_TYPE_
{
  distr_type_undefined=0,
  distr_type_bingham=1,
  distr_type_vonmisesfisher=2
} DISTR_TYPE_;

typedef enum _STRATEGY_TYPE_
{
  strategy_type_undefined=0,
  strategy_type_standard=1,
  strategy_type_bydistributionfunction=2
} STRATEGY_TYPE_;

/*----------------------------------------------------------------------*/
/// material parameters for generalized structural tensor with distribution
/// around mean fiber direction M.
///
/// We assume only one symmetry for the distribution function \rho(M), i.e.:
/// \rho(M) = \rho(-M) .
///
/// Example input line in .dat file:
///
/// MAT 2 ELAST_IsoAnisoExpo  K1 1.0E6 K2 100.0 GAMMA 0.0 K1COMP 0.0 K2COMP 0.0 STR_TENS_ID 100
/// MAT 100 ELAST_StructuralTensor STRATEGY ByDistributionFunction DISTR vonMisesFisher C1 500.0
///
class StructuralTensorParameter : public MAT::PAR::Parameter
{
public:

  /// standard constructor
  StructuralTensorParameter(
      Teuchos::RCP<MAT::PAR::Material> matdata
  );

  /// destructor
  virtual ~StructuralTensorParameter() { ; }

  /// @name material parameters
  //@{
  double c1_;     //!< constant 1 for distribution function
  double c2_;     //!< constant 2 for distribution function
  double c3_;     //!< constant 3 for distribution function
  double c4_;     //!< constant 4 for distribution function
  //@}

  /// type of distribution function around mean fiber direction
  int distribution_type_;

  /// type of structural tensor strategy
  int strategy_type_;

  /// pure virtual in base class
  virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

};  // class StructuralTensorParameter
}  // namespace PAR

/*----------------------------------------------------------------------*/
/// Base class for evaluation strategy of structural tensor for anisotropic
/// materials.
///
/// \author rauch
/// \date 10/17
class StructuralTensorStrategyBase
{

public:

  /// constructor
  StructuralTensorStrategyBase(MAT::ELASTIC::PAR::StructuralTensorParameter* params);

  /// destructor
  virtual ~StructuralTensorStrategyBase(){ ; };

  /// method for computing the structural tensor for anisotropic materials.
  /// this is the core functionality of this object.
  /// each derived class has to implement this method (pure virtual).
  virtual void SetupStructuralTensor(
      LINALG::Matrix<3,1>  &fiber_vector,
      LINALG::Matrix<6,1>  &structural_tensor
  ) = 0;

protected:

  /// return residual tolerance of structural problem
  double GetResidualTol();

  /// my material parameters
  MAT::ELASTIC::PAR::StructuralTensorParameter* params_;

}; // class StructuralTensorStrategyBase


/*----------------------------------------------------------------------*/
/// Class for computing the standard structural tensor for anisotropic
/// materials via a dyadic product of the current fiber direction.
///
/// <h3>Definiton of Structural Tensor H</h3>
///
/// By H we denote the structural tensor. M is the direction of the fiber.
///
/// H = M x M
///
/// Here, x denotes the dyadic product.
///
/// \author rauch
/// \date 10/17
class StructuralTensorStrategyStandard : public StructuralTensorStrategyBase
{
public:

  /// constructor with given material parameters
  StructuralTensorStrategyStandard(MAT::ELASTIC::PAR::StructuralTensorParameter* params);

  /// \brief method for computing the structural tensor for anisotropic materials
  ///
  /// Simplest variant assuming perfect alignment of fiber with a given fiber direction M.
  ///
  /// <h3>Definiton of Structural Tensor H</h3>
  ///
  ///  H = M x M (x denotes the dydic product)
  ///
  /// \param fiber_vector (in) : direction of fiber 'M'
  /// \param structural_tensor (out) : structural tensor is 'H' filled inside
  ///
  virtual void SetupStructuralTensor(
      LINALG::Matrix<3,1>& fiber_vector,
      LINALG::Matrix<6,1>& structural_tensor
  );

}; // class StructuralTensorStrategyStandard


/*----------------------------------------------------------------------*/
/// \brief Class for computing the generalized structural tensor
///
/// Class for computing the generalized structural tensor for anisotropic
/// materials, by integrating a given distribution function for fiber
/// orientations.
///
///  See:
///  T. C. Gasser, R. W. Ogden, and G. A. Holzapfel. Hyperelastic modelling of arterial layers
///  with distributed collagen fibre orientations. J. R. Soc. Interface, 3:15–35, 2006. doi:
///  10.1098/rsif.2005.0073.
///
/// <h3>Definiton of Structural Tensor H</h3>
///
/// H = \frac{1}{4\pi}\int \rho(Theta,Phi) [m x m] sin(Theta) dTheta dPhi
/// with
/// m = \sin(\Theta)\cos(\Phi)e_1 + \sin(\Theta)\sin(\Phi)e_2 + \cos(\Theta)e_3 .
///
///
/// <h3>Distribution functions rho</h3>
///
/// By \rho we denote the distribution function. For now we can choose between:
///
/// 1) vonMisesFisher distribution:
///
///    This distribution function has only one parameter C1 determining the dispersion
///    of fibers around the given mean direction M.
///
///     \rho = \frac{C1}{4\pi\sinh(C1)} \exp(C1 M \cdot m)
///
///     with 0 \leq C1 , |M| = 1 , and |m| = 1
///
///    See also: https://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution
///
///  2) Bingham distribution:
///
///     This distribution function requires four parameters C1, C2, C3, and C4.
///
///     \rho(\Theta,\Phi) = \frac{1}{C4} \exp(C1 X_1 + C2 X_2 + C3 X_3)
///
///     with X_1 = m_1^2
///          X_2 = m_2^2
///          X_3 = m_3^2
///
///     See:
///     T. Gasser, S. Gallinetti, X. Xing, C. Forsell, J. Swedenborg, and J. Roy. Spatial orientation
///     of collagen fibers in the abdominal aortic aneurysm’s wall and its relation to wall mechanics
///     Acta Biomater., 8:3091–3103, 2012. doi: http://dx.doi.org/10.1016/j.actbio.2012.04.044.
///
/// The integration is performed by an Gauss integration rule in accordance with:
///
/// Atkinson K.; Numerical Integration on the Sphere.; J. Austral. Math. Soc.B; 23; 1982
///
/// \author rauch
/// \date 10/17
class StructuralTensorStrategyByDistributionFunction : public StructuralTensorStrategyBase
{
public:

  /// constructor with given material parameters
  StructuralTensorStrategyByDistributionFunction(MAT::ELASTIC::PAR::StructuralTensorParameter* params);

  /*! \brief evaluate generalized structural tensor with given distribution function

    Here we integrate the general integral for the structural tensor:

    H = \frac{1}{4\pi}\int \rho(Theta,Phi) [m x m] sin(Theta) dTheta dPhi

    wherein H denotes the structural tensor tensor, \rho is the distribution
    function for the dispersion of fibers around the mean fiber direction.
    Here, M denotes the direction vector in each direction (Theta,Phi).
    In case of a perfect alignment of all fibres the integral yields the standard

    \param fiber_vector (in) : mean direction of fiber
    \param structural_tensor (out) : generalized structural tensor

    \author rauch
    \date 10/17    */
  virtual void SetupStructuralTensor(
      LINALG::Matrix<3,1>& fiber_vector,
      LINALG::Matrix<6,1>& structural_tensor
  );
}; // class StructuralTensorStrategyByDistributionFunction

}  // namespace ELASTIC
}  // namespace MAT

#endif  // ELAST_STRUCTURALTENSOR_STRATEGY_
