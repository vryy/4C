/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I) specifically related to one-way-coupled problem
       configurations

\level 2

\maintainer Christoph Schmidt


*----------------------------------------------------------------------*/


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

    virtual void Init();

    virtual void Setup();

    virtual void Timeloop();

    void DoFSIStep();

    void DoScatraStep();

    virtual void PrepareTimeStep();

    virtual bool ScatraConvergenceCheck(int itnum);
  };

}  // namespace FS3I

#endif
