/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I) specifically related to one-way-coupled problem
       configurations

\level 2



*----------------------------------------------------------------------*/


#ifndef BACI_FS3I_PARTITIONED_1WC_H
#define BACI_FS3I_PARTITIONED_1WC_H


#include "baci_fs3i_partitioned.H"


namespace FS3I
{
  class PartFS3I_1WC : public PartFS3I
  {
   public:
    PartFS3I_1WC(const Epetra_Comm& comm);

    void Init() override;

    void Setup() override;

    void Timeloop() override;

    void DoFSIStep();

    void DoScatraStep();

    void PrepareTimeStep() override;

    bool ScatraConvergenceCheck(int itnum) override;
  };

}  // namespace FS3I

#endif
