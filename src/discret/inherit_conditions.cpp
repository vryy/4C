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
#ifdef TRILINOS_PACKAGE

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


// static methods
static void inherit_dirichlet_high_to_low_entity_elements(
                            RefCountPtr<DRT::DesignDiscretization> highdis,
                            RefCountPtr<DRT::DesignDiscretization> lowdis);
static void inherit_dirichlet_element_to_node(
                            RefCountPtr<DRT::DesignDiscretization> dis);
static void inherit_dirichlet_design_to_discretization(
                                        DRT::DesignDiscretization& ddis,
                                        DRT::Discretization&        dis);

static void inherit_neumann_designelements_to_discretization(
                                      DRT::DesignDiscretization& ddis,
                                      DRT::Discretization&        dis,
                                      const int ndele,
                                      const vector<int>& ndele_fenode,
                                      const vector<vector<int> >& dele_fenode);

static void inherit_neumann_designnodes_to_discretization(
                                      DRT::DesignDiscretization& ddis,
                                      DRT::Discretization&        dis,
                                      const int ndnode,
                                      const vector<int>& ndnode_fenode,
                                      const vector<vector<int> >& dnode_fenode);

/*----------------------------------------------------------------------*
 | input of conditions                                    m.gee 11/06   |
 *----------------------------------------------------------------------*/
void inherit_conditions()
{
  DSTraceHelper dst("inherit_conditions");

  if (design)
  {
    RefCountPtr<DRT::Design>* tmp = (RefCountPtr<DRT::Design>*)design->ccadesign;
    DRT::Design& ccadesign = **tmp;
    RefCountPtr<DRT::DesignDiscretization> designlines = ccadesign[0];
    RefCountPtr<DRT::DesignDiscretization> designsurfs = ccadesign[1];
    RefCountPtr<DRT::DesignDiscretization> designvols  = ccadesign[2];
  
    // number of design points, lines, surfaces and volumes
    const int ndnode = designlines->NumGlobalNodes();
    const int ndline = designlines->NumGlobalElements();
    const int ndsurf = designsurfs->NumGlobalElements();
    const int ndvol  = designvols->NumGlobalElements();

    if (ccadesign.Comm().MyPID()==0)
    {
      // dirichlet conditions are inherited as follows:
      // DVOL inherits to its DSURFS if the DSURF does not have its own
      // DSURF inherits to its DLINEs if the DLINE does not have its own
      // DLINE inherits to its DNODEs if the DNODE does not have its own
      inherit_dirichlet_high_to_low_entity_elements(designvols,designsurfs);
      inherit_dirichlet_high_to_low_entity_elements(designsurfs,designlines);
      inherit_dirichlet_element_to_node(designlines);
      
      // read design nodes <-> nodes
      vector<int> ndnode_fenode(ndnode);
      vector<vector<int> > dnode_fenode(ndnode);
      for (int i=0; i<ndnode; ++i) 
        ndnode_fenode[i] = 0; 
      input_design_dpoint_fenode_read(dnode_fenode,ndnode_fenode);
    
      // read design lines <-> nodes
      vector<int> ndline_fenode(ndline);
      vector<vector<int> > dline_fenode(ndline);
      for (int i=0; i<ndline; ++i) 
        ndline_fenode[i] = 0;
      input_design_dline_fenode_read(dline_fenode,ndline_fenode); 
  
      // read design surfaces <-> nodes
      vector<int> ndsurf_fenode(ndsurf);
      vector<vector<int> > dsurf_fenode(ndsurf);
      for (int i=0; i<ndsurf; ++i) 
        ndsurf_fenode[i] = 0;
      input_design_dsurf_fenode_read(dsurf_fenode,ndsurf_fenode);

      // read design volumes <-> nodes
      vector<int> ndvol_fenode(ndvol);
      vector<vector<int> > dvol_fenode(ndvol);
      for (int i=0; i<ndvol; ++i) 
        ndvol_fenode[i] = 0;
      input_design_dvol_fenode_read(dvol_fenode,ndvol_fenode);

      // inherit all conditions from design to discretization
      for (int i=0; i<genprob.numfld; ++i)
      {
        vector<RefCountPtr<DRT::Discretization> >* discretization =
                      (vector<RefCountPtr<DRT::Discretization> >*)field[i].ccadis;
        for (int j=0; j<field[i].ndis; ++j)
        {
          // inherit the dirichlet boundary conditions from design
          // to the discretization
          RefCountPtr<DRT::Discretization> actdis = (*discretization)[j];
          inherit_dirichlet_design_to_discretization(*designlines,*actdis);
          inherit_dirichlet_design_to_discretization(*designsurfs,*actdis);
          inherit_dirichlet_design_to_discretization(*designvols,*actdis);
          // inherit neumann conditions from design to discretization 
          // based on the sets read above
          inherit_neumann_designelements_to_discretization(
                                             *designvols,*actdis,
                                             ndvol,ndvol_fenode,dvol_fenode);
          inherit_neumann_designelements_to_discretization(
                                             *designsurfs,*actdis,
                                             ndsurf,ndsurf_fenode,dsurf_fenode);
          inherit_neumann_designelements_to_discretization(
                                             *designlines,*actdis,
                                             ndline,ndline_fenode,dline_fenode);
          inherit_neumann_designnodes_to_discretization(
                                             *designlines,*actdis,
                                             ndnode,ndnode_fenode,dnode_fenode);
        }
      } // for (int i=0; i<genprob.numfld; ++i)
    } // if (ccadesign.Comm().MyPID()==0)
  } // if (design)
  else
  {
    dserror("Design free boundary conditions not yet impl.");
  }

  return;
} /* end of inherit_conditions */


/*----------------------------------------------------------------------*
 | inherit Neumann conditions from design set to nodes      m.gee 12/06 |
 *----------------------------------------------------------------------*/
void inherit_neumann_designnodes_to_discretization(
                                      DRT::DesignDiscretization& ddis,
                                      DRT::Discretization&        dis,
                                      const int ndnode,
                                      const vector<int>& ndnode_fenode,
                                      const vector<vector<int> >& dnode_fenode)
{
  DSTraceHelper dst("inherit_neumann_designnodes_to_discretization");
  if (!dis.Filled()) dserror("FillComplete() must have been called before");
  if (!ddis.Filled()) dserror("FillComplete() must have been called before");
  
  for (int i=0; i<ddis.NumMyColNodes(); ++i)
  {
    DRT::Node* dnode = ddis.lColNode(i);
    if (!dnode) dserror("Cannot get lColNode");
    const int dnodeid = dnode->Id();
    
    // check design node for PointNeumann condition
    DRT::Condition* neumpoint  = dnode->GetCondition("PointNeumann");
    if (!neumpoint) continue;
    
    const int nfenode = ndnode_fenode[dnodeid];
    for (int j=0; j<nfenode; ++j)
    {
      const int gid = dnode_fenode[dnodeid][j];
      DRT::Node* node = dis.gNode(gid);
      if (!node) dserror("Cannot find fe node");
      node->SetCondition("PointNeumann",rcp(new DRT::Condition(*neumpoint)));
    }
  }
  
  return;
}                                      



/*----------------------------------------------------------------------*
 | inherit Neumann conditions from design set to nodes      m.gee 12/06 |
 *----------------------------------------------------------------------*/
void inherit_neumann_designelements_to_discretization(
                                      DRT::DesignDiscretization& ddis,
                                      DRT::Discretization&        dis,
                                      const int ndele,
                                      const vector<int>& ndele_fenode,
                                      const vector<vector<int> >& dele_fenode)
{
  DSTraceHelper dst("inherit_neumann_designelements_to_discretization");
  if (!dis.Filled()) dserror("FillComplete() must have been called before");
  if (!ddis.Filled()) dserror("FillComplete() must have been called before");
  
  for (int i=0; i<ddis.NumMyColElements(); ++i)
  {
    DRT::Element* ele = ddis.lColElement(i);
    if (!ele) dserror("Cannot get lColElement");
    const int eleid = ele->Id();
    
    // check design Element for VolumeNeumann, SurfaceNeumann,
    //                          LineNeumann condition
    DRT::Condition* neumvol  = ele->GetCondition("VolumeNeumann");
    DRT::Condition* neumsurf = ele->GetCondition("SurfaceNeumann");
    DRT::Condition* neumline = ele->GetCondition("LineNeumann");
    
    const int nfenode = ndele_fenode[eleid];
    
    for (int j=0; j<nfenode; ++j)
    {
      const int gid = dele_fenode[eleid][j];
      DRT::Node* node = dis.gNode(gid);
      if (!node) dserror("Cannot find fe node");
      
      if (neumvol)
        node->SetCondition("VolumeNeumann",rcp(new DRT::Condition(*neumvol)));
      if (neumsurf)
        node->SetCondition("SurfaceNeumann",rcp(new DRT::Condition(*neumsurf)));
      if (neumline)
        node->SetCondition("LineNeumann",rcp(new DRT::Condition(*neumline)));
    }
  }
  
  return;
}                                      



/*----------------------------------------------------------------------*
 | inherit Dirichlet conditions from eles to nodes          m.gee 11/06 |
 *----------------------------------------------------------------------*/
void inherit_dirichlet_design_to_discretization(
                                        DRT::DesignDiscretization& ddis,
                                        DRT::Discretization&        dis)
{
  DSTraceHelper dst("inherit_dirichlet_design_to_discretization");
  
  // it is implicitly assumed that all neccessary design objects are
  // available on this processor.
  // Normally, when this method is called, everything is still on proc 0
  
  if (!ddis.Filled()) dserror("FillComplete was not called on design discretization");
  if (!dis.Filled()) dserror("FillComplete was not called on discretization");

  // loop nodes in discretization 
  for (int i=0; i<dis.NumMyColNodes(); ++i)
  {
    DRT::Node* actnode = dis.lColNode(i);
    
    // get design entity type and id this node is on
    DRT::Node::OnDesignEntity type = actnode->GetDesignEntityType();
    if (type==DRT::Node::on_none) continue;
    int dgid = actnode->GetDesignEntityId();
    
    // get the design entity and inherit its dirichlet condition
    if (type==DRT::Node::on_dnode)
    {
      if (!ddis.HaveGlobalNode(dgid)) continue;
      DRT::Node* dnode = ddis.gNode(dgid);
      DRT::Condition* dirich = dnode->GetCondition("Dirichlet");
      if (!dirich) continue;
      // deep copy and set in actnode
      RefCountPtr<DRT::Condition> newdirich = rcp(new DRT::Condition(*dirich));
      actnode->SetCondition("Dirichlet",newdirich);
    }
    else
    {
      if (!ddis.HaveGlobalElement(dgid)) continue;
      DRT::Element* dele = ddis.gElement(dgid);
      
      // make sure its the correct type
      DRT::Element::ElementType dtype = dele->Type();
      if (type  == DRT::Node::on_dline && 
          dtype != DRT::Element::element_designline) continue;
      if (type  == DRT::Node::on_dsurface && 
          dtype != DRT::Element::element_designsurface) continue;
      if (type  == DRT::Node::on_dvolume && 
          dtype != DRT::Element::element_designvolume) continue;
      
      DRT::Condition* dirich = dele->GetCondition("Dirichlet");
      if (!dirich) continue;
      // deep copy and set in actnode
      RefCountPtr<DRT::Condition> newdirich = rcp(new DRT::Condition(*dirich));
      actnode->SetCondition("Dirichlet",newdirich);
    }
  } // for (int i=0; i<dis.NumMyColNodes(); ++i)
  
  return;
} // inherit_dirichlet_design_to_discretization


/*----------------------------------------------------------------------*
 | inherit Dirichlet conditions from eles to nodes          m.gee 11/06 |
 *----------------------------------------------------------------------*/
void inherit_dirichlet_element_to_node(
                               RefCountPtr<DRT::DesignDiscretization> dis)
{
  DSTraceHelper dst("inherit_dirichlet_element_to_node");
  if (!dis->Filled()) dserror("FillComplete() must have been called before");

  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    // get high entity element
    DRT::Element* ele = dis->lColElement(i);
    if (!ele) dserror("Cannot get lElement");
    
    // get dirichlet condition from high entity element
    const DRT::Condition* highcond = ele->GetCondition("Dirichlet");
    if (!highcond) continue;
    
    // get number and ptrs to its nodes
    int nnode = ele->NumNode();
    DRT::Node** nodes = ele->Nodes();
    if (!nodes) dserror("node ptrs not set though Fillcomplete() was called");
    
    // loop lower entities and put condition if there's none there already
    for (int j=0; j<nnode; ++j)
    {
      const DRT::Condition* lowcond = nodes[j]->GetCondition("Dirichlet");
      if (lowcond) continue;
      RefCountPtr<DRT::Condition> newcond = rcp(new DRT::Condition(*highcond));
      nodes[j]->SetCondition("Dirichlet",newcond);
    }
  }
  
  return;
} // inherit_dirichlet_element_to_node



/*----------------------------------------------------------------------*
 | inherit Dirichlet conditions from higher to lower eles   m.gee 11/06 |
 *----------------------------------------------------------------------*/
void inherit_dirichlet_high_to_low_entity_elements(
                          RefCountPtr<DRT::DesignDiscretization> highdis,
                          RefCountPtr<DRT::DesignDiscretization> lowdis)
{
  DSTraceHelper dst("inherit_dirichlet_high_to_low_entity_elements");
  if (!highdis->Filled()) dserror("FillComplete() must have been called before");
  if (!lowdis->Filled())  dserror("FillComplete() must have been called before");

  for (int i=0; i<highdis->NumMyColElements(); ++i)
  {
    // get high entity element
    DRT::DesignElement* highele = 
                   dynamic_cast<DRT::DesignElement*>(highdis->lColElement(i));
    if (!highele) dserror("dynamic_cast to DesignElement failed");
    
    // get dirichlet condition from high entity element
    const DRT::Condition* highcond = highele->GetCondition("Dirichlet");
    if (!highcond) continue;
    
    // get number and ptrs to its lower entity elements
    int numlower = highele->NumLowerEntityIds();
    if (!numlower) continue;
    DRT::Element** lowele = highele->LowerEntities();
    if (!lowele) dserror("low entity ptrs not set though Fillcomplete() was called");
    
    // loop lower entities and put condition if there's none there already
    for (int j=0; j<numlower; ++j)
    {
      const DRT::Condition* lowcond = lowele[j]->GetCondition("Dirichlet");
      if (lowcond) continue;
      RefCountPtr<DRT::Condition> newcond = rcp(new DRT::Condition(*highcond));
      lowele[j]->SetCondition("Dirichlet",newcond);
    }
  }
  
  return;
} // inherit_dirichlet_high_to_low_entity_elements


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
