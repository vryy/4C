/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Epetra_SerialDenseMatrix.h"
#include "global_inp_control2.H"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

#ifdef DEBUG
/*!----------------------------------------------------------------------
  \brief the tracing variable

  <pre>                                                         m.gee 8/00
  defined in pss_ds.c, declared in tracing.h
  </pre>
 *----------------------------------------------------------------------*/
extern struct _CCA_TRACE         trace;
#endif

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*----------------------------------------------------------------------*
 * Read any topology and make sure node GIDs are sorted
 *----------------------------------------------------------------------*/
static void input_design_read(const std::string& name,
                              std::vector<std::vector<int> >& dobj_fenode,
                              std::vector<int>& ndobj_fenode)
{
  std::map<int,std::set<int> > topology;

  std::string sectionname = name + "-NODE TOPOLOGY";
  std::string marker = std::string("--") + sectionname;
  if (frfind(const_cast<char*>(marker.c_str())))
  {
    frread();

    // read the whole thing
    while (strncmp(allfiles.actplace,"------",6)!=0)
    {
      int dobj;
      int ierr;
      int nodeid;
      frint(const_cast<char*>(name.c_str()),&dobj,&ierr);
      if (ierr!=1)
        dserror("Cannot read %s", sectionname.c_str());
      frint("NODE",&nodeid,&ierr);
      if (ierr!=1)
        dserror("Cannot read %s", sectionname.c_str());
      topology[dobj-1].insert(nodeid-1);
      frread();
    }

    // copy all design object entries
    for (std::map<int,std::set<int> >::iterator i=topology.begin();
         i!=topology.end();
         ++i)
    {
      if (i->first >= static_cast<int>(dobj_fenode.size()))
      {
        dserror("Illegal design object number %d in section '%s'", i->first+1, sectionname.c_str());
      }

      // this is probably obsolete
      ndobj_fenode[i->first] = i->second.size();

      // we copy from a std::set, thus the gids are sorted
      dobj_fenode[i->first].reserve(i->second.size());
      dobj_fenode[i->first].assign(i->second.begin(),i->second.end());
    }
  }
}


/*----------------------------------------------------------------------*
 | input of design volumes to fe-node topology            m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_design_dvol_fenode_read(vector<vector<int> >& dvol_fenode,
                                   vector<int>& ndvol_fenode)
{
  input_design_read("DVOL",dvol_fenode,ndvol_fenode);
}


/*----------------------------------------------------------------------*
 | input of design surface to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
void input_design_dsurf_fenode_read(vector<vector<int> >& dsurf_fenode,
                                     vector<int>& ndsurf_fenode)
{
  input_design_read("DSURF",dsurf_fenode,ndsurf_fenode);
}


/*----------------------------------------------------------------------*
 | input of design line  to fe-node topology             m.gee 11/06    |
 *----------------------------------------------------------------------*/
void input_design_dline_fenode_read(vector<vector<int> >& dline_fenode,
                                     vector<int>& ndline_fenode)
{
  input_design_read("DLINE",dline_fenode,ndline_fenode);
}


/*----------------------------------------------------------------------*
 | input of design nodes  to fe-node topology             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_design_dpoint_fenode_read(vector<vector<int> >& dnode_fenode,
                                     vector<int>& ndnode_fenode)
{
  input_design_read("DNODE",dnode_fenode,ndnode_fenode);
}


#endif  // #ifdef CCADISCRET
