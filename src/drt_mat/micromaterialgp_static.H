/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIALGP_STAT_H
#define MICROMATERIALGP_STAT_H



#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations

namespace UTILS
{
  class SurfStressManager;
}

namespace STRUMULTI
{
  class MicroStatic;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace MAT
{
  /// one Gauss point of the micro material
  class MicroMaterialGP
  {
   public:
    /// construct an instance of MicroMaterial for a given Gauss point
    /// and microscale discretization
    MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner, const int microdisnum,
        const double V0);

    /// destructor
    ~MicroMaterialGP();

    /// Read restart
    void ReadRestart();

    /// New result file
    void NewResultFile(bool eleowner, std::string& newfilename);

    /// Perform microscale simulation
    void PerformMicroSimulation(
        LINALG::Matrix<3, 3>* defgrd, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat);

    void Update();

    /// Calculate stresses and strains on the micro-scale
    void PrepareOutput();

    void Output();

    /// Create and initialize "empty" EAS history map
    void EasInit();

    /// Reset global time and step number (needed for multi-scale inverse
    /// analyses with multiple runs)
    void ResetTimeAndStep();

    /// get ele id
    int eleID() { return ele_ID_; }

    /// get density
    double Density() const { return density_; }


   private:
    /// corresponding macroscale Gauss point
    const int gp_;

    /// corresponding macroscale element
    const int ele_ID_;

    /// corresponding microstructure discretization number
    const int microdisnum_;

    /// microstructure "time integration" classes (one for each micro-discretization)
    static std::map<int, Teuchos::RCP<STRUMULTI::MicroStatic>> microstaticmap_;

    static std::map<int, int> microstaticcounter_;

    /// microstructure discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> micro_output_;

    /// homogenized density
    double density_;

    /// my vector of old displacements
    Teuchos::RCP<Epetra_Vector> dis_;

    /// my vector of new displacements
    Teuchos::RCP<Epetra_Vector> disn_;

    /// my manager of surface stresses
    Teuchos::RCP<UTILS::SurfStressManager> surf_stress_man_;

    // my EAS history data -> note that microstructure is not parallel
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> lastalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> oldalpha_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> oldfeas_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> oldKaainv_;
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> oldKda_;

    /// my stresses and strains
    Teuchos::RCP<std::vector<char>> stress_;
    Teuchos::RCP<std::vector<char>> strain_;
    Teuchos::RCP<std::vector<char>> plstrain_;

    /// old absolute time
    double time_;

    /// current absolute time
    double timen_;

    /// old step
    int step_;

    /// current step
    int stepn_;

    /// timestep size
    double dt_;

    /// restart name
    std::string restartname_;

    /// flag for modified Newton on macroscale
    bool mod_newton_;

    /// flag for build of stiffness matrix
    bool build_stiff_;
  };
}  // namespace MAT

#endif
