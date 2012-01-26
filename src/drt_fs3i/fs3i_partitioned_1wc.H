#ifdef CCADISCRET

#ifndef FS3I_PART_1WC_H
#define FS3I_PART_1WC_H

#include "fs3i_partitioned.H"

namespace FS3I
{
  class PartFS3I_1WC : public PartFS3I
  {
  public:

    PartFS3I_1WC(const Epetra_Comm& comm);

    virtual ~PartFS3I_1WC(){};

    virtual void Timeloop();

    void DoFSIStep();

    void DoScatraStep();

    virtual bool ScatraConvergenceCheck(int itnum);
  };

}

#endif
#endif
