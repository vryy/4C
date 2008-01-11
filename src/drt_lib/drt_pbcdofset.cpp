/*!----------------------------------------------------------------------
\file drt_pbcdofset.H

\brief A modified set of degrees of freedom for periodic boundary
       conditions

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_pbcdofset.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "../headers/define_sizes.h"





/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::PBCDofSet::PBCDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes)
{
  perbndcouples_=couplednodes;
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::PBCDofSet::~PBCDofSet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  setup everything  (public)                               gammi 05/07|
 |  this function is a specialisation of the AssignDegreesOfFreedom of  |
 |  the base class DofSet                                               |
 *----------------------------------------------------------------------*/
int DRT::PBCDofSet::AssignDegreesOfFreedom(const Discretization& dis, const int start)
{
  if (!dis.Filled()) dserror("discretization Filled()==false");
  if (!dis.NodeRowMap()->UniqueGIDs()) dserror("Nodal row map is not unique");
  if (!dis.ElementRowMap()->UniqueGIDs()) dserror("Element row map is not unique");

  // Now this is tricky. We have to care for nodes and elements, both
  // row and column maps. In general both nodes and elements can have
  // dofs. In both cases these dofs might be shared with other nodes
  // or elements. (The very general case. For elements we'd probably
  // don't need that.)
  //
  // The point is that we have to make sure the dof numbering of a
  // mesh is independent of its parallel distribution. Otherwise we
  // could not redistribute a mesh. We would not be able to use old
  // distributed vectors afterwards.
  //
  // Each object (node or element) could have a different number of
  // dofs. The parallel distribution is arbitrary. So we fall back to
  // two redundant vectors here to gather the number of dofs per node
  // or element.

  // numdf for all nodes and elements
  numdfcolnodes_ = rcp(new Epetra_IntVector(*dis.NodeColMap()));
  numdfcolelements_ = rcp(new Epetra_IntVector(*dis.ElementColMap()));

  // index of first dof for all nodes and elements
  idxcolnodes_ = rcp(new Epetra_IntVector(*dis.NodeColMap()));
  idxcolelements_ = rcp(new Epetra_IntVector(*dis.ElementColMap()));

  //////////////////////////////////////////////////////////////////

  // do the nodes first

  map<int,int> nidx;
  Utils::AllreduceEMap(nidx, *dis.NodeRowMap());

  // build a redundant that holds all the node's numdof
  vector<int> sredundantnodes(dis.NumGlobalNodes());
  fill(sredundantnodes.begin(), sredundantnodes.end(), 0);

  int lnumdof = 0;

  vector<int> localnodenumdf(dis.NumMyRowNodes());

  // loop my row nodes and remember the number of degrees of freedom
  for (int i=0; i<dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const int numele = actnode->NumElement();
    DRT::Element** myele = actnode->Elements();
    int numdf=0;
    for (int j=0; j<numele; ++j)
      numdf = max(numdf,myele[j]->NumDofPerNode(*actnode));
    const int gid = actnode->Id();
    sredundantnodes[nidx[gid]] = numdf;
    lnumdof += numdf;
  }

  if(!perbndcouples_->empty())
  {
    // loop all slave nodes and set the number of dofs to zero
    for(map<int,vector<int> >::iterator curr = perbndcouples_->begin();
        curr != perbndcouples_->end();
        ++curr )
    {

      if(!dis.NodeColMap()->MyGID(curr->first))
      {
        dserror("master not ghosted on proc but in connectivity");
      }

      for(vector<int>::iterator slave=curr->second.begin();
          slave!=curr->second.end();
          ++slave)
      {
#ifdef DEBUG
        if((!dis.NodeColMap()->MyGID(*slave)) &&
           (dis.NodeRowMap()->MyGID(curr->first)))
        {
          char *string;
          sprintf(string,"slave %d to master %d  expected to be on that proc %d\n",*slave,curr->first,dis.Comm().MyPID());
          cout << string;
          dserror(string);
        }
#endif
        // get the number of dofs associated with this slave
        int numdf = sredundantnodes[nidx[*slave]];
        // reset them to zero
        sredundantnodes[nidx[*slave]] =0;
        // reduce the number of local degrees of freedom
        lnumdof-=numdf;
      }
    }
  }

  // make all nodal numdfs globally known
  vector<int> rredundantnodes(dis.NumGlobalNodes());
  dis.Comm().SumAll(&sredundantnodes[0],&rredundantnodes[0],dis.NumGlobalNodes());

  int count=start;
  int localcolpos=0;

  // Vectors to keep the dof gid of both nodes and elements. We need
  // both the row and column fashion to create the row and column map
  // later on.
  vector<int> localrowdofs;
  localrowdofs.reserve(lnumdof); // exact

  vector<int> localcoldofs;
  localcoldofs.reserve(lnumdof); // just a guess

  // We know all the numdfs of all the nodes and we even know the
  // nodes from our local row map. Use that information.
  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::iterator i=nidx.begin(); i!=nidx.end(); ++i)
  {
    int numdf = sredundantnodes[i->second];

    // we know our row nodes from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col nodes we might have to look up
    if ((numdf!=0) || dis.NodeColMap()->MyGID(i->first))
    {
      numdf = rredundantnodes[i->second];

      // remember numdf and index for nodal based lookup
      int lid = dis.NodeColMap()->LID(i->first);
      (*numdfcolnodes_)[lid] = numdf;
      (*idxcolnodes_)[lid] = localcolpos;
      localcolpos += numdf;

      for (int j=0; j<numdf; ++j)
      {
        localcoldofs.push_back(count+j);
      }
    }
    else
    {
      // but in any case we need to increase our dof count
      numdf = rredundantnodes[i->second];
    }

    count += numdf;
  }

  if(!perbndcouples_->empty())
  {
    // loop all master nodes and set the degrees of freedom of
    // the slaves to the degrees of freedom of the master
    for(map<int,vector<int> >::iterator master = perbndcouples_->begin();
        master != perbndcouples_->end();
        ++master )
    {


      int master_lid=dis.NodeColMap()->LID(master->first);

      if(master_lid<0)
      {
        dserror("master not on proc\n");
      }


      for(vector<int>::iterator slave=master->second.begin();
          slave!=master->second.end();
          ++slave)
      {
        int slave_lid=dis.NodeColMap()->LID(*slave);

        if(slave_lid>-1)
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
  }

  sredundantnodes.clear();
  rredundantnodes.clear();
  nidx.clear();

  //////////////////////////////////////////////////////////////////

  // Now do all this fun again for the elements

  //////////////////////////////////////////////////////////////////

  // Now do all this fun again for the elements

  map<int,int> eidx;
  Utils::AllreduceEMap(eidx, *dis.ElementRowMap());

  vector<int> sredundantelements(dis.NumGlobalElements());
  fill(sredundantelements.begin(), sredundantelements.end(), 0);

  vector<int> localelenumdf(dis.NumMyRowElements());

  // loop my row elements and set number of degrees of freedom
  for (int i=0; i<dis.NumMyRowElements(); ++i)
  {
    DRT::Element* actele = dis.lRowElement(i);
    const int gid = actele->Id();
    int numdf = actele->NumDofPerElement();
    sredundantelements[eidx[gid]] = numdf;
    localelenumdf[i] = numdf;
    lnumdof += numdf;
  }

  vector<int> rredundantelements(dis.NumGlobalElements());
  dis.Comm().SumAll(&sredundantelements[0],&rredundantelements[0],dis.NumGlobalElements());

  // enlarge the big dof vectors
  localrowdofs.reserve(lnumdof); // exact
  localcoldofs.reserve(lnumdof); // guess

  // We have to loop the gids in order. This way we will get an
  // ordered set of dofs.
  for (map<int,int>::iterator i=eidx.begin(); i!=eidx.end(); ++i)
  {
    int numdf = sredundantelements[i->second];

    // we know our row elements from the redundant vector
    if (numdf!=0)
    {
      for (int j=0; j<numdf; ++j)
      {
        localrowdofs.push_back(count+j);
      }
    }

    // the col elements we might have to look up
    if ((numdf!=0) || dis.ElementColMap()->MyGID(i->first))
    {
      numdf = rredundantelements[i->second];

      // remember numdf and index for elemental based lookup
      int lid = dis.ElementColMap()->LID(i->first);
      (*numdfcolelements_)[lid] = numdf;
      (*idxcolelements_)[lid] = localcolpos;
      localcolpos += numdf;

      for (int j=0; j<numdf; ++j)
      {
        localcoldofs.push_back(count+j);
      }
    }
    else
    {
      // but in any case we need to increase our dof count
      numdf = rredundantelements[i->second];
    }

    count += numdf;
  }

  sredundantelements.clear();
  rredundantelements.clear();
  eidx.clear();

  // Now finally we have everything in place to build the maps.

  dofrowmap_ = rcp(new Epetra_Map(-1,localrowdofs.size(),&localrowdofs[0],0,dis.Comm()));
  if (!dofrowmap_->UniqueGIDs()) dserror("Dof row map is not unique");

  dofcolmap_ = rcp(new Epetra_Map(-1,localcoldofs.size(),&localcoldofs[0],0,dis.Comm()));

  return count;
}


#endif  // #ifdef CCADISCRET
