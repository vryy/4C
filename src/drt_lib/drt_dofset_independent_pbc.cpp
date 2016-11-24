/*!----------------------------------------------------------------------
\file drt_dofset_independent_pbc.cpp

\brief Independent dof set for periodic boundary conditions

\level 2

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "drt_dofset_independent_pbc.H"


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::IndependentPBCDofSet::IndependentPBCDofSet(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes)
  :DRT::PBCDofSet(couplednodes)
{
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                                       |
 *----------------------------------------------------------------------*/
DRT::IndependentPBCDofSet::~IndependentPBCDofSet()
{
  return;
}


int DRT::IndependentPBCDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // assign dofs for the standard dofset, that is without periodic boundary
  // conditions and using the independent dofset's ADOF
  int count = DRT::IndependentDofSet::AssignDegreesOfFreedom(dis,dspos,start);

  myMaxGID_ = DRT::DofSet::MaxAllGID();

  // loop all master nodes and set the dofs of the slaves to the dofs of the master
  // remark: the previously assigned dofs of slave nodes are overwritten here
  for(std::map<int,std::vector<int> >::iterator master = perbndcouples_->begin();
      master != perbndcouples_->end();
      ++master )
  {
    int master_lid=dis.NodeColMap()->LID(master->first);

    if (master_lid<0)
    {
      dserror("master gid %d not on proc %d, but required by slave %d",
              master->first,
              dis.Comm().MyPID(),
              master->second[0]);
    }

    for (std::vector<int>::iterator slave=master->second.begin();
         slave!=master->second.end();
         ++slave)
    {
      int slave_lid=dis.NodeColMap()->LID(*slave);

      if (slave_lid>-1)
      {
        (*numdfcolnodes_)[slave_lid] = (*numdfcolnodes_)[master_lid];
        (*idxcolnodes_)  [slave_lid] = (*idxcolnodes_)  [master_lid];
      }
      else
      {
#ifdef DEBUG
        if(dis.NodeRowMap()->MyGID(master->first))
        {
          dserror("slave not on proc but master owned by proc\n");
        }
#endif
      }
    }
  }

  return count;
}
