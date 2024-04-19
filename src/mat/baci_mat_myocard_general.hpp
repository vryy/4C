/*----------------------------------------------------------------------*/
/*! \file
\brief general myocard material model

\level 3

*/

#ifndef FOUR_C_MAT_MYOCARD_GENERAL_HPP
#define FOUR_C_MAT_MYOCARD_GENERAL_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

class MyocardGeneral

{
 public:
  /// construct empty material object
  MyocardGeneral(){};

  /// construct empty material object
  explicit MyocardGeneral(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~MyocardGeneral() = default;
  /// compute reaction coefficient for multiple Gauss points
  virtual double ReaCoeff(const double phi, const double dt, int gp)
  {
    if (gp > 0)
    {
      FOUR_C_THROW("Multiple Gauss points only implemented for MV and FHN model");
      return 0.;
    }

    else
      return ReaCoeff(phi, dt);
  };

  /// compute reaction coefficient
  virtual double ReaCoeff(const double phi, const double dt) = 0;

  /// compute reaction coefficient at timestep n
  virtual double ReaCoeffN(const double phi, const double dt) { return 0; };

  /// compute reaction coefficient for multiple Gauss points at timestep n
  virtual double ReaCoeffN(const double phi, const double dt, int gp)
  {
    if (gp > 0)
    {
      FOUR_C_THROW("Multiple Gauss points only implemented for MV and FHN model");
      return 0.;
    }

    else
      return ReaCoeffN(phi, dt);
  };

  ///  returns number of internal state variables of the material
  virtual int GetNumberOfInternalStateVariables() const = 0;

  ///  returns current internal state of the material
  virtual double GetInternalState(const int k) const = 0;

  ///  returns current internal state of the material
  virtual double GetInternalState(const int k, int gp) const
  {
    if (gp > 0)
    {
      FOUR_C_THROW("Multiple Gauss points only implemented for MV and FHN model");
      return 0.;
    }
    else
      return GetInternalState(k);
  };

  ///  set internal state of the material
  virtual void SetInternalState(const int k, const double val) = 0;

  ///  set internal state of the material for multiple Gauss points
  virtual void SetInternalState(const int k, const double val, int gp)
  {
    if (gp > 0)
      FOUR_C_THROW("Multiple Gauss points only implemented for MV and FHN model");
    else
      SetInternalState(k, val);

    return;
  };

  ///  return number of ionic currents
  virtual int GetNumberOfIonicCurrents() const = 0;

  ///  return ionic currents
  virtual double GetIonicCurrents(const int k) const = 0;

  ///  return ionic currents
  virtual double GetIonicCurrents(const int k, int gp) const
  {
    FOUR_C_THROW("Multiple Gauss Points only implemented for MV and FHN model");
    return 0.;
  };

  // resize internal state variables if number of Gauss point changes
  virtual void ResizeInternalStateVariables(int gp)
  {
    FOUR_C_THROW("Multiple Gauss Points only implemented for MV and FHN model");
    return;
  };

  /// get number of Gauss points
  virtual int GetNumberOfGP() const { return 1; };

  /// time update for this material
  virtual void Update(const double phi, const double dt) = 0;

};  // Myocard_general


FOUR_C_NAMESPACE_CLOSE

#endif
