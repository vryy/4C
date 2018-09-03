/*----------------------------------------------------------------------------*/
/**
\file xstr_xstructure_structure_state.H

\brief XStructure/Structure state handling ( combination of XFEM and std.
       discretizations )

\maintainer Michael Hiermeier

\date Jun 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_
#define SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_

#include "xstr_xstructure_state.H"

// forward declarations
namespace XFEM
{
  class MultiFieldMapExtractor;
}  // namespace XFEM
namespace XSTR
{
  class XStructureStructureState : public XStructureState
  {
   public:
    /// constructor
    XStructureStructureState();

    /// destructor
    virtual ~XStructureStructureState(){};

    virtual void Setup();

  };  // class XStructureStructureState
}  // namespace XSTR


#endif /* SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_ */
