/*----------------------------------------------------------------------------*/
/**
\file xstr_xstructure_state.H

\brief XStructure state handling ( only one single XFEM-discretization )

\maintainer Michael Hiermeier

\date Jun 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_
#define SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_

#include "../drt_xfem/xfield_state.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

namespace XSTR {
  class XStructureState : public XFEM::XFieldState,
                          public STR::TIMINT::BaseDataGlobalState
  {
  public:
    /// constructor
    XStructureState();

    /// destructor
    virtual ~XStructureState() {};

    /// setup of member variables
    virtual void Setup();

    virtual bool Destroy();

  protected:
    inline void CheckInit() const
    {
      XFEM::XFieldState::CheckInit();
      STR::TIMINT::BaseDataGlobalState::CheckInit();
    }

    inline void CheckInitSetup() const
    {
      XFEM::XFieldState::CheckInitSetup();
      STR::TIMINT::BaseDataGlobalState::CheckInitSetup();
    }
  };  // class XStructureState
} // XSTR


#endif /* SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_ */
