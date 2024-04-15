/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I) specifically related to one-way-coupled problem
       configurations

\level 2



*----------------------------------------------------------------------*/


#ifndef FOUR_C_FS3I_PARTITIONED_1WC_HPP
#define FOUR_C_FS3I_PARTITIONED_1WC_HPP


#include "baci_config.hpp"

#include "baci_fs3i_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FS3I
{
  class PartFS3I1Wc : public PartFS3I
  {
   public:
    PartFS3I1Wc(const Epetra_Comm& comm);

    void Init() override;

    void Setup() override;

    void Timeloop() override;

    void DoFSIStep();

    void DoScatraStep();

    void PrepareTimeStep() override;

    bool ScatraConvergenceCheck(int itnum) override;
  };

}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
