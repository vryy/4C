/*---------------------------------------------------------------------*/
/*! \file
\brief Data for airway elements
\level 3
*/
/*---------------------------------------------------------------------*/

#ifndef RED_AIRWAYS_EVALUATION_DATA_H
#define RED_AIRWAYS_EVALUATION_DATA_H

#include <Teuchos_RCPDecl.hpp>
#include <Epetra_Vector.h>

namespace DRT::REDAIRWAYS
{
  /**
   * Store all global vectors that are required to evaluate different kind of reduced airway
   * elements.
   */
  struct EvaluationData
  {
   protected:
    EvaluationData() = default;

   public:
    Teuchos::RCP<Epetra_Vector> acinar_vnp_strain;
    Teuchos::RCP<Epetra_Vector> acinar_vnp;
    Teuchos::RCP<Epetra_Vector> acinar_vn;
    Teuchos::RCP<Epetra_Vector> acinar_v;

    Teuchos::RCP<Epetra_Vector> qin_nm;
    Teuchos::RCP<Epetra_Vector> qin_n;
    Teuchos::RCP<Epetra_Vector> qin_np;

    Teuchos::RCP<Epetra_Vector> x_n;
    Teuchos::RCP<Epetra_Vector> x_np;
    Teuchos::RCP<Epetra_Vector> open;

    Teuchos::RCP<Epetra_Vector> p_extn;
    Teuchos::RCP<Epetra_Vector> p_extnp;
    Teuchos::RCP<Epetra_Vector> airway_acinus_dep;
    bool compute_awacinter{};

    Teuchos::RCP<Epetra_Vector> qout_np;
    Teuchos::RCP<Epetra_Vector> qout_n;
    Teuchos::RCP<Epetra_Vector> qout_nm;

    Teuchos::RCP<Epetra_Vector> p0np;
    Teuchos::RCP<Epetra_Vector> p0n;
    Teuchos::RCP<Epetra_Vector> p0nm;
    Teuchos::RCP<Epetra_Vector> elemArea0;

    Teuchos::RCP<Epetra_Vector> acini_e_volume;
    Teuchos::RCP<Epetra_Vector> elemVolume;
    Teuchos::RCP<Epetra_Vector> elemVolumen;
    Teuchos::RCP<Epetra_Vector> elemVolumenp;

    Teuchos::RCP<Epetra_Vector> generations;

    bool solveScatra{};

    Teuchos::RCP<Epetra_Vector> junVolMix_Corrector;
    Teuchos::RCP<Epetra_Vector> scatran;
    Teuchos::RCP<Epetra_Vector> scatranp;

    Teuchos::RCP<Epetra_Vector> e1scatran;
    Teuchos::RCP<Epetra_Vector> e2scatran;

    Teuchos::RCP<Epetra_Vector> e1scatranp;
    Teuchos::RCP<Epetra_Vector> e2scatranp;

    Teuchos::RCP<Epetra_Vector> dscatranp;

    Teuchos::RCP<Epetra_Vector> elemRadiusnp;

    Teuchos::RCP<Epetra_Vector> cfl;

    Teuchos::RCP<Epetra_Vector> po2;

    Teuchos::RCP<Epetra_Vector> bcval;
    Teuchos::RCP<Epetra_Vector> dbctog;

    Teuchos::RCP<Epetra_Vector> acini_bc;

    double lungVolume_np{};
    double lungVolume_n{};
    double lungVolume_nm{};

    double time{};
    double dt{};

    static EvaluationData& get()
    {
      static EvaluationData evaluation_data;
      return evaluation_data;
    }

    EvaluationData(const EvaluationData&) = delete;
    EvaluationData(EvaluationData&&) = delete;
    EvaluationData& operator=(const EvaluationData&) = delete;
    EvaluationData& operator=(EvaluationData&&) = delete;
  };
}  // namespace DRT::REDAIRWAYS

#endif
