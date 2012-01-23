#ifdef CCADISCRET

#ifndef FS3I_PART_2WC_H
#define FS3I_PART_2WC_H

#include "fs3i_partitioned.H"

namespace FS3I
{
  class PartFS3I_2WC : public PartFS3I
  {
  public:

    PartFS3I_2WC(const Epetra_Comm& comm);

    virtual ~PartFS3I_2WC(){};

    virtual void Timeloop();

    virtual bool ScatraConvergenceCheck(int itnum);

  private:

    //! @name  parameters controlling iteration between FSI and ScaTra
    int itermax_;
    double convtol_;
    //@}

    /// flag for constant thermodynamic pressure
    string consthermpress_;
  };
}

#endif
#endif
