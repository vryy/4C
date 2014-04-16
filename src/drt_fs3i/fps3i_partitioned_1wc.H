/*!----------------------------------------------------------------------
\file fps3i_partitioned_1wc.H
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I) specifically related to one-way-coupled problem
       configurations

<pre>
Maintainers: Moritz Thon
             thon@mhpc.mw.tum.de
             089/289-10324
</pre>

*----------------------------------------------------------------------*/


#ifndef FPS3I_PART_1WC_H
#define FPS3I_PART_1WC_H

#include "fs3i_partitioned.H"

namespace FS3I
{
  class PartFPS3I_1WC : public PartFS3I
  {
  public:

    PartFPS3I_1WC(const Epetra_Comm& comm);

    virtual ~PartFPS3I_1WC(){};

    virtual void Timeloop();

    //! This is actually a fpsi step
    void DoFSIStep();

    void DoScatraStep();

    virtual void PrepareTimeStep();

    virtual bool ScatraConvergenceCheck(int itnum);
  };

}

#endif
