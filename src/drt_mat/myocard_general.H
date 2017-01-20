/*----------------------------------------------------------------------*/
/*!
\file myocard_general.H
\brief general myocard material model

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_GENERAL_H
#define MYOCARD_GENERAL_H

/*----------------------------------------------------------------------*
 |  headers                                                  cbert 08/13 |
 *----------------------------------------------------------------------*/
/*#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
*/
// #include "myocard_minimal.H"

/// \date 08/13

class Myocard_General

{

public:

  /// construct empty material object
  Myocard_General(){};

  /// construct empty material object
  explicit Myocard_General(const double eps_deriv_myocard, const std::string tissue);

  /// destructor
  virtual ~Myocard_General() { ; }

  /// compute reaction coefficient for multiple Gauss points
  virtual double ReaCoeff(const double phi, const double dt, int gp)
  {
    if(gp>0)
    {
      dserror("Multiple Gauss points only implemented for MV and FHN model");
      return 0.;
    }

    else
      return ReaCoeff(phi,dt);
  };

  /// compute reaction coefficient
  virtual double ReaCoeff(const double phi, const double dt) = 0;

  /// compute reaction coefficient at timestep n
  virtual double ReaCoeffN(const double phi, const double dt) {return 0;};

  /// compute reaction coefficient for multiple Gauss points at timestep n
  virtual double ReaCoeffN(const double phi, const double dt, int gp)
  {
    if(gp>0)
  {
    dserror("Multiple Gauss points only implemented for MV and FHN model");
    return 0.;
  }

  else
    return ReaCoeffN(phi,dt);
  };

  ///  returns number of internal state variables of the material
  virtual int GetNumberOfInternalStateVariables() const = 0;

  ///  returns current internal state of the material
  virtual double GetInternalState(const int k) const = 0;

  ///  returns current internal state of the material
  virtual double GetInternalState(const int k, int gp) const
  {
    if(gp > 0)
    {
      dserror("Multiple Gauss points only implemented for MV and FHN model");
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
    if(gp > 0)
      dserror("Multiple Gauss points only implemented for MV and FHN model");
    else
      SetInternalState(k,val);

    return;
  };

  ///  return number of ionic currents
  virtual int GetNumberOfIonicCurrents() const = 0;

  ///  return ionic currents
  virtual double GetIonicCurrents(const int k) const = 0;

  ///  return ionic currents
  virtual double GetIonicCurrents(const int k, int gp) const
  {
    dserror("Multiple Gauss Points only implemented for MV and FHN model");
    return 0.;
  };

  // resize internal state variables if number of Gauss point changes
  virtual void ResizeInternalStateVariables(int gp)
  {
    dserror("Multiple Gauss Points only implemented for MV and FHN model");
    return ;
  };

  /// get number of Gauss points
  virtual int GetNumberOfGP() const
  {
    return 1;
  };

  /// time update for this material
  virtual void Update(const double phi, const double dt) = 0;

}; // Myocard_general


#endif
