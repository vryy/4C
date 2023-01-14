/*----------------------------------------------------------------------------*/
/** \file

\brief XStructure/Structure state handling ( combination of XFEM and std.
       discretizations )


\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef SRC_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_
#define SRC_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_

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


#endif /* SRC_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STRUCTURE_STATE_H_ */
