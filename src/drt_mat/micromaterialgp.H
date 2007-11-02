/*----------------------------------------------------------------------*/
/*!
\file micromaterial.H

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIALGP_H
#define MICROMATERIALGP_H

#ifdef D_SOH8
#ifdef CCADISCRET

#include "micromaterial.H"

using namespace std;
using namespace Teuchos;
using namespace IO;

namespace MAT
{

/// one Gauss point of the great micro material
class MicroMaterialGP
{
public:

  /// construct an instance of MicroMaterial for a given Gauss point
  /// and microscale discretization
  MicroMaterialGP(const int gp, const int ele_ID);

  /// destructor
  ~MicroMaterialGP();


  /// Perform microscale simulation
  void PerformMicroSimulation(const Epetra_SerialDenseMatrix* defgrd,
                              Epetra_SerialDenseVector* stress,
                              Epetra_SerialDenseMatrix* cmat,
                              double* density,
                              const double time,
                              string action);

  /// Set up microscale generalized alpha
  void SetUpMicroGenAlpha();


private:

  /// my Gauss point
  const int gp_;

  /// my macroscopic element
  const int ele_ID_;

  /// microstructure generalized alpha class
  static RefCountPtr<MicroStruGenAlpha> microgenalpha_;

  /// microstructure discretization writer
  RefCountPtr<MicroDiscretizationWriter> micro_output_;

  /// my vector of displacements
  RefCountPtr<Epetra_Vector> disp_;

  /// my vector of velocities
  RefCountPtr<Epetra_Vector> vel_;

  /// my vector of accelerations
  RefCountPtr<Epetra_Vector> acc_;

  /// my vector of residual displacements
  RefCountPtr<Epetra_Vector> disi_;

  /// current absolute time
  double timen_;

  /// current step
  int istep_;

  /// timestep size
  double dt_;

};
}

#endif
#endif
#endif
