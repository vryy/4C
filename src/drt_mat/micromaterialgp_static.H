/*----------------------------------------------------------------------*/
/*!
\file micromaterialgp_static.H

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIALGP_STAT_H
#define MICROMATERIALGP_STAT_H

#ifdef CCADISCRET

#include "micromaterial.H"
#include "../drt_surfstress/drt_surfstress_manager.H"

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
  MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner,
                  const double time, const int microdisnum, const double V0);

  /// destructor
  ~MicroMaterialGP();

  /// Read restart
  void ReadRestart();

  /// Perform microscale simulation
  void PerformMicroSimulation(LINALG::Matrix<3,3>* defgrd,
                              LINALG::Matrix<6,1>* stress,
                              LINALG::Matrix<6,6>* cmat,
                              double* density,
                              const double time,
                              const double dt,
                              const bool eleowner);

  /// Set up microscale static "time integrator"
  void SetUpMicroStatic(const double V0);

  /// Create and initialize "empty" EAS history map
  void EasInit();


private:

  /// corresponding macroscale Gauss point
  const int gp_;

  /// corresponding macroscale element
  const int ele_ID_;

  /// corresponding microstructure discretization number
  const int microdisnum_;

  /// microstructure "time integration" classes (one for each micro-discretization)
  static std::map<int, RefCountPtr<STRUMULTI::MicroStatic> > microstaticmap_;

  static std::map<int, int> microstaticcounter_;

  /// microstructure discretization writer
  RefCountPtr<DiscretizationWriter> micro_output_;

  /// my vector of old displacements
  RefCountPtr<Epetra_Vector> dis_;

  /// my vector of displacements at the generalized midpoint
  RefCountPtr<Epetra_Vector> dism_;

  /// my vector of new displacements
  RefCountPtr<Epetra_Vector> disn_;

  /// my manager of surface stresses
  RefCountPtr<UTILS::SurfStressManager> surf_stress_man_;

  // my EAS history data -> note that microstructure is not parallel
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > lastalpha_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda_;

  /// current absolute time
  double timen_;

  /// current step
  int istep_;

  /// timestep size
  double dt_;

  /// restart name
  string restartname_;

  /// flag for modified Newton on macroscale
  bool mod_newton_;

  /// flag for build of stiffness matrix
  bool build_stiff_;

};
}

#endif
#endif
