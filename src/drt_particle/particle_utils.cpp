/*----------------------------------------------------------------------*/
/*!
\file particle_utils.cpp

\brief General functions for the Particle-MeshFree dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#include "particle_utils.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/extparticle_mat.H"

/*-----------------------------------------------------------------------------*/
/* Compute the inertia vector */
void PARTICLE::Utils::ComputeInertia(
  const Teuchos::RCP<const Epetra_Vector> radius,
  const Teuchos::RCP<const Epetra_Vector> mass,
  Teuchos::RCP<Epetra_Vector> &inertia,
  bool trg_createInertiaVector)
{
  // checks
  if (radius == Teuchos::null ||
      mass == Teuchos::null)
    dserror("radius or mass vectors are empty");

  // rebuild the inertia vector
  if (trg_createInertiaVector || inertia == Teuchos::null)
    inertia = Teuchos::rcp(new Epetra_Vector(mass->Map(), true));

  // compute inertia for every particle
  for (int lidNode = 0; lidNode < mass->MyLength(); ++lidNode)
    (*inertia)[lidNode] = ComputeInertia((*radius)[lidNode], (*mass)[lidNode]);
}


/*----------------------------------------------------------------------*/
/* compute temperature from the specEnthalpy
 * TODO: (see reference... equation...)
 * */
double PARTICLE::Utils::SpecEnthalpy2Temperature(
    const double& specEnthalpy,
    const MAT::PAR::ExtParticleMat* extParticleMat)
{
  //extract the interesting parameters
  const double specEnthalpyST = extParticleMat->SpecEnthalpyST();
  const double specEnthalpyTL = extParticleMat->SpecEnthalpyTL();
  const double transitionTemperature = extParticleMat->transitionTemperature_;
  const double CPS = extParticleMat->CPS_;
  const double inv_CPS = 1/CPS;
  const double CPL = extParticleMat->CPL_;
  const double inv_CPL = 1/CPL;

  // compute temperature of the node
  if (specEnthalpy < specEnthalpyST)
    return specEnthalpy * inv_CPS;
  else if (specEnthalpy > specEnthalpyTL)
    return transitionTemperature + (specEnthalpy - specEnthalpyTL) * inv_CPL;
  else
    return transitionTemperature;
}


/*----------------------------------------------------------------------*/
/* compute temperature from the specEnthalpy - (vector version) */
void PARTICLE::Utils::SpecEnthalpy2Temperature(
    Teuchos::RCP<Epetra_Vector> temperature,
    const Teuchos::RCP<const Epetra_Vector>& specEnthalpy,
    const MAT::PAR::ExtParticleMat* extParticleMat)
{
  // check: no specEnthalpy? no temperature :)
  if (specEnthalpy == Teuchos::null)
    dserror("specEnthalpy is a null pointer!");
  if (temperature == Teuchos::null)
    dserror("temperature is a null pointer!");
  if (!(temperature->Map().SameAs(specEnthalpy->Map())))
    dserror("temperature map and specEnthalpy map missmatch!");

  for (int lidNode = 0; lidNode < specEnthalpy->MyLength(); ++lidNode)
  {
    (*temperature)[lidNode] = PARTICLE::Utils::SpecEnthalpy2Temperature((*specEnthalpy)[lidNode],extParticleMat);
  }
}

/*-----------------------------------------------------------------------------*/
// Compute pressure vector (see Antoci2007-E4)
void PARTICLE::Utils::Density2Pressure(
  const double& restDensity,
  const double& refdensfac,
  const Teuchos::RCP<const Epetra_Vector> density,
  const Teuchos::RCP<const Epetra_Vector> specEnthalpy,
  Teuchos::RCP<Epetra_Vector> &pressure,
  const MAT::PAR::ExtParticleMat* extParticleMat,
  bool trg_createPressureVector,
  bool solve_thermal_problem)
{
  // checks
  if (density == Teuchos::null)
  {
    pressure = Teuchos::null;
    return;
  }
  if (specEnthalpy == Teuchos::null)
    dserror("specEnthalpy is a null pointer!");

  // rebuild the pressure vector
  if (trg_createPressureVector || pressure == Teuchos::null)
  {
    pressure = Teuchos::rcp(new Epetra_Vector(density->Map(), true));
  }

  // compute inertia for every particle
  for (int lidNode = 0; lidNode < density->MyLength(); ++lidNode)
  {
    double speedOfSound = 0.0;
    //If no thermal problem is solved, only liquid phases are considered
    if(solve_thermal_problem)
      speedOfSound = SpeedOfSound((*specEnthalpy)[lidNode], extParticleMat);
    else
      speedOfSound = extParticleMat->SpeedOfSoundL();

    (*pressure)[lidNode] = Density2Pressure(speedOfSound, restDensity, refdensfac, (*density)[lidNode],extParticleMat->exponent_);
  }
}


/*-----------------------------------------------------------------------------*/
// compute the intersection area of two particles that are in contact. It returns 0 if there is no contact
// TODO: (see reference... equation...)
double PARTICLE::Utils::IntersectionAreaPvsP(const double& radius1, const double& radius2, const double& dis)
{
  //checks
  if (radius1<=0 || radius2<=0 || dis<=0)
    dserror("input parameters are unacceptable");
  if (dis >= radius1 + radius2)
    return 0;

  return M_PI * (- 0.25 * dis * dis
                 - radius1 * radius1 * radius1 * radius1 / (4 * dis * dis)
                 - radius2 * radius2 * radius2 * radius2 / (4 * dis * dis)
                 + 0.5 * ( radius1 * radius1 + radius2 * radius2 + radius1 * radius1 * radius2 * radius2));

}

/*-----------------------------------------------------------------------------*/
// Provide the correct speed of sound
double PARTICLE::Utils::SpeedOfSound(
    const double &specEnthalpy,
    const MAT::PAR::ExtParticleMat* extParticleMat)
{

  if (specEnthalpy <= extParticleMat->SpecEnthalpyST())
  {
    return extParticleMat->SpeedOfSoundS();
  }
  else if (specEnthalpy >= extParticleMat->SpecEnthalpyTL())
  {
    return extParticleMat->SpeedOfSoundL();
  }
  else
  {
    return extParticleMat->SpeedOfSoundT(specEnthalpy);
  }
}

