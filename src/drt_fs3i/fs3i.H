#ifdef CCADISCRET

#ifndef FS3I_H
#define FS3I_H

#include <Epetra_Comm.h>

namespace FS3I
{
  class FS3I_Base
  {

  public:

    /// constructor of base class
    FS3I_Base(){};

    /// destructor of base class
    virtual ~FS3I_Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) = 0;

    /// read restart
    virtual void ReadRestart() = 0;
  };
}

#endif
#endif
