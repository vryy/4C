/*----------------------------------------------------------------------------*/
/** \file

\brief XStructure/Structure state handling ( combination of XFEM and std.
       discretizations )


\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef STRUCTURE_XSTRUCTURE_STATE_H
#define STRUCTURE_XSTRUCTURE_STATE_H

#include "structure_xstructure_xstr_state.H"

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


#endif  // STRUCTURE_XSTRUCTURE_STATE_H
