/*! \file
\brief Interface for every material that can evaluate thermo material laws

\level 3

\maintainer Sebastian Proell
*/

#ifndef TRAIT_THERMO_H
#define TRAIT_THERMO_H

namespace MAT
{
  namespace TRAIT
  {
    class Thermo
    {
     public:
      virtual ~Thermo() = default;

      //! Main material call to determine heat flux and constitutive tensor in 3D
      virtual void Evaluate(
          const LINALG::Matrix<3, 1>& gradtemp,  ///< temperature gradient (strain tensor)
          LINALG::Matrix<3, 3>& cmat,            ///< constitutive matrix
          LINALG::Matrix<3, 1>& heatflux         ///< heatflux
          ) const = 0;

      //! Main material call to determine heat flux and constitutive tensor in 2D
      virtual void Evaluate(
          const LINALG::Matrix<2, 1>& gradtemp,  ///< temperature gradient (strain tensor)
          LINALG::Matrix<2, 2>& cmat,            ///< constitutive matrix
          LINALG::Matrix<2, 1>& heatflux         ///< heatflux
          ) const = 0;

      //! Main material call to determine heat flux and constitutive tensor in 1D
      virtual void Evaluate(
          const LINALG::Matrix<1, 1>& gradtemp,  ///< temperature gradient (strain tensor)
          LINALG::Matrix<1, 1>& cmat,            ///< constitutive matrix
          LINALG::Matrix<1, 1>& heatflux         ///< heatflux
          ) const = 0;


      //! @name Derivatives of conductivity tensor
      //! @{

      //! @brief Derivative of conductivity tensor wrt to temperature
      virtual void ConductivityDerivT(LINALG::Matrix<3, 3>& dCondDT) const = 0;

      virtual void ConductivityDerivT(LINALG::Matrix<2, 2>& dCondDT) const = 0;

      virtual void ConductivityDerivT(LINALG::Matrix<1, 1>& dCondDT) const = 0;

      //! @}

      //! @brief get volumetric heat capacity
      //!
      //! @pre the state must be set by Reinit() if necessary
      virtual double Capacity() const = 0;

      //! @brief get derivative of volumetric heat capacity wrt temperature
      //!
      //! @pre state must be set by Reinit() if necessary.
      virtual double CapacityDerivT() const = 0;

      //! Set necessary variables for Evaluation
      virtual void Reinit(double temperature, unsigned gp) = 0;

      //! reset current state e.g. due to Newton failed
      virtual void ResetCurrentState() = 0;

      //! persist currently set state to history
      virtual void CommitCurrentState() = 0;
    };
  }  // namespace TRAIT
}  // namespace MAT

#endif
