/*!----------------------------------------------------------------------
\file contact_augmented_interface_tools.cpp

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Created on: Apr 11, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#include "contact_augmented_interface.H"
#include "../drt_contact/contact_integrator.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_integrator.H"
#include "../drt_mortar/mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"


/*----------------------------------------------------------------------*
 | Finite difference check for KappaLin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FDCheckKappaLin()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  std::map<int,double> refKappa;
  std::map<int,std::map<int,double> > refKappaLin;

  double newKappa = 0.0;

  int dim = Dim();

  // print reference to screen (kappaLin-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if ((int) cnode->CoData().GetKappaLin().size()==0)
      continue;

    refKappa[gid]    = cnode->CoData().GetKappa();
    refKappaLin[gid] = cnode->CoData().GetKappaLin();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd=0;fd<dim*snodefullmap->NumMyElements();++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd/dim);
    CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if(!snode)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    int sdof = snode->Dofs()[fd%dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newKappa = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetKappaLin().size()==0)
        continue;

      newKappa = kcnode->CoData().GetKappa();

      // print results (derivatives) to screen
      if (abs(newKappa-refKappa[kgid]) > 1e-12)
      {
        double finit = (newKappa-refKappa[kgid])/delta;

        double analy = refKappaLin[kgid][sdof];
        double dev = finit - analy;

        // p->first: currently tested dof of slave node p->first/Dim()
        // (p->first)/Dim(): paired slave
        // sdof: currently modified slave dof
        std::cout << "(" << kgid << "," << sdof << ") :"
            "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
            "   kappaLin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
            "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
            "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // global loop to apply FD scheme to all MASTER dofs (=dim*nodes)
  for (int fd=0; fd<dim*mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd/dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* mnode = static_cast<CoNode*>(node);

    int mdof = mnode->Dofs()[fd%dim];

    std::cout << "\nDEVIATION FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newKappa = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetKappaLin().size()==0)
              continue;

      newKappa = kcnode->CoData().GetKappa();

      // print results (derivatives) to screen

      if (abs(newKappa-refKappa[kgid]) > 1e-12)
      {
        double finit = (newKappa-refKappa[kgid])/delta;
        double analy = refKappaLin[kgid][mdof];
        double dev = finit - analy;

        // p->first: currently tested dof of master node p->first/Dim()
        // (p->first)/Dim(): paired master
        // mdof: currently modified master dof
        std::cout << "(" << kgid << "," << mdof << ") :"
        "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
        "   kappaLin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
        "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
        "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
     mnode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
     mnode->xspatial()[1] -= delta;
    }
    else
    {
     mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  Evaluate();
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for AWGapLin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FDCheckAWGapLin()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  std::map<int,double> refAWGap;
  std::map<int,std::map<int,double> > refAWGapLin;

  double newAWGap = 0.0;

  int dim = Dim();

  // print reference to screen (kappaLin-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i=0;i<augActiveSlaveNodes_->NumMyElements();++i)
  {
    int gid = augActiveSlaveNodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if ((int) cnode->CoData().GetAWGapLin().size()==0)
      dserror("Linearization map shouldn't be empty for active nodes!");

    double wg = cnode->CoData().GetWGap();
    double kappa = cnode->CoData().GetKappa();
    if (wg==1.0e12 or wg==0.0)
      dserror("The weighted gap is at the active node % equal 0.0 or 1.0e12!",gid);
    if (kappa==1.0e12 or kappa==0.0)
      dserror("The scaling factor kappa is at the active node % equal 0.0 or 1.0e12!",gid);

    refAWGap[gid]    = wg/kappa;
    refAWGapLin[gid] = cnode->CoData().GetAWGapLin();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd=0;fd<dim*snodefullmap->NumMyElements();++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd/dim);
    CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if(!snode)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    int sdof = snode->Dofs()[fd%dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();
    WGap();
    AWGapLin();

    // compute finite difference derivative
    for (int k=0; k<augActiveSlaveNodes_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAWGap = 0.0;
      int kgid = augActiveSlaveNodes_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetAWGapLin().size()==0)
        dserror("Linearization map shouldn't be empty for active nodes!");

      double wg = kcnode->CoData().GetWGap();
      double kappa = kcnode->CoData().GetKappa();
      if (wg==1.0e12 or wg==0.0)
        dserror("The weighted gap is at the active node % equal 0.0 or 1.0e12!",kgid);
      if (kappa==1.0e12 or kappa==0.0)
        dserror("The scaling factor kappa is at the active node % equal 0.0 or 1.0e12!",kgid);

      newAWGap = wg/kappa;

      // print results (derivatives) to screen
      if (abs(newAWGap-refAWGap[kgid]) > 1e-12)
      {
        double finit = (newAWGap-refAWGap[kgid])/delta;

        double analy = refAWGapLin[kgid][sdof];
        double dev = finit - analy;

        // p->first: currently tested dof of slave node p->first/Dim()
        // (p->first)/Dim(): paired slave
        // sdof: currently modified slave dof
        std::cout << "(" << kgid << "," << sdof << ") :"
            "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
            "   AWGapLin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
            "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
            "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // global loop to apply FD scheme to all MASTER dofs (=dim*nodes)
  for (int fd=0; fd<dim*mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd/dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* mnode = static_cast<CoNode*>(node);

    int mdof = mnode->Dofs()[fd%dim];

    std::cout << "\nDEVIATION FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();
    WGap();
    AWGapLin();

    // compute finite difference derivative
    for (int k=0; k<augActiveSlaveNodes_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAWGap = 0.0;
      int kgid = augActiveSlaveNodes_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetAWGapLin().size()==0)
        dserror("Linearization map shouldn't be empty for active nodes!");

      double wg = kcnode->CoData().GetWGap();
      double kappa = kcnode->CoData().GetKappa();
      if (wg==1.0e12 or wg==0.0)
        dserror("The weighted gap is at the active node % equal 0.0 or 1.0e12!",kgid);
      if (kappa==1.0e12 or kappa==0.0)
        dserror("The scaling factor kappa is at the active node % equal 0.0 or 1.0e12!",kgid);

      newAWGap = wg/kappa;

      // print results (derivatives) to screen

      if (abs(newAWGap-refAWGap[kgid]) > 1e-12)
      {
        double finit = (newAWGap-refAWGap[kgid])/delta;
        double analy = refAWGapLin[kgid][mdof];
        double dev = finit - analy;

        // p->first: currently tested dof of master node p->first/Dim()
        // (p->first)/Dim(): paired master
        // mdof: currently modified master dof
        std::cout << "(" << kgid << "," << mdof << ") :"
        "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
        "   AWGapLin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
        "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
        "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
     mnode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
     mnode->xspatial()[1] -= delta;
    }
    else
    {
     mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  Evaluate();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for VarWGapLinSl              hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FDCheckVarWGapLinSl()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // first integration loop has to be activated again
  IParams().set<int>("AugLagStep",0);

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // create storage for Dn-Matrix entries
  // stores the nodal vector:
  std::map<int,std::map<int,std::pair<int,double> > > refVarWGapSl;

  std::map<int,std::map<int,std::map<int,double> > > refVarWGapLinSl; // stores old derivdn for every node
  int dim = Dim();

  // store the current values for later comparison
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if ((int) cnode->CoData().GetVarWGapSl().size()==0)
      continue;

    // get reference values
    refVarWGapSl[gid].insert(cnode->CoData().GetVarWGapSl().begin(),cnode->CoData().GetVarWGapSl().end());

    refVarWGapLinSl[gid] = cnode->CoData().GetVarWGapLinSl();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd=0;fd<dim*snodefullmap->NumMyElements();++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd/dim);
    CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if(!snode)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    int sdof = snode->Dofs()[fd%dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
//      newVarWGapSl.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int)(kcnode->CoData().GetVarWGapSl().size())==0)
        continue;

      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

      GEN::pairedvector<int,std::pair<int,double> >& newVarWGapSl = kcnode->CoData().GetVarWGapSl();

      // print results (derivatives) to screen
      for (CI p=newVarWGapSl.begin(); p!=newVarWGapSl.end(); ++p)
      {
        if (abs(newVarWGapSl[p->first].second-refVarWGapSl[kgid][p->first].second) > 1e-12)
        {
          double finit = (newVarWGapSl[p->first].second-refVarWGapSl[kgid][p->first].second)/delta;

          double analy = ((refVarWGapLinSl[kgid])[p->first])[sdof];
          double dev = finit - analy;

          // p->first: currently tested dof of slave node p->first/Dim()
          // (p->first)/Dim(): paired slave
          // sdof: currently modified slave dof
          std::cout << "(" << p->first << "," << p->first/dim <<
              "(" << kgid << ")" << "," << sdof << ") :"
              "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
              "   varWGapLinSl=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
              "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
              "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

          if( abs(dev) > 1e-4 )
          {
            std::cout << " ***** WARNING ***** ";
            w++;
          }
          else if( abs(dev) > 1e-5 )
          {
            std::cout << " ***** warning ***** ";
            w++;
          }

          std::cout << std::endl;
        }
      }
    }

    // undo finite difference modification
    if (fd%dim==0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // global loop to apply FD scheme to all MASTER dofs (=dim*nodes)
  for (int fd=0; fd<dim*mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd/dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* mnode = static_cast<CoNode*>(node);

    int mdof = mnode->Dofs()[fd%dim];

    std::cout << "\nDEVIATION FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
//      newVarWGapSl.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      int dim = kcnode->NumDof();

      if ((int)(kcnode->CoData().GetVarWGapLinSl().size())==0)
        continue;

      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

      GEN::pairedvector<int,std::pair<int,double> >& newVarWGapSl = kcnode->CoData().GetVarWGapSl();

      // print results (derivatives) to screen
      for (CI p=newVarWGapSl.begin(); p!=newVarWGapSl.end(); ++p)
      {
        if (abs(newVarWGapSl[p->first].second-refVarWGapSl[kgid][p->first].second) > 1e-12)
        {
          double finit = (newVarWGapSl[p->first].second-refVarWGapSl[kgid][p->first].second)/delta;
          double analy = ((refVarWGapLinSl[kgid])[(p->first)])[mdof];
          double dev = finit - analy;

          // p->first: currently tested dof of master node p->first/Dim()
          // (p->first)/Dim(): paired master
          // mdof: currently modified master dof
          std::cout << "(" << p->first << "," << p->first/dim <<
          "(" << kgid << ")" << "," << mdof << ") :"
          "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
          "   varWGapLinSl=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
          "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
          "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

          if( abs(dev) > 1e-4 )
          {
            std::cout << " ***** WARNING ***** ";
            w++;
          }
          else if( abs(dev) > 1e-5 )
          {
            std::cout << " ***** warning ***** ";
            w++;
          }

          std::cout << std::endl;
        }
      }
    }

    // undo finite difference modification
    if (fd%dim==0)
    {
     mnode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
     mnode->xspatial()[1] -= delta;
    }
    else
    {
     mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  Evaluate();
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for VarWGapLinMa              hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FDCheckVarWGapLinMa()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // first integration loop has to be activated again

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // create storage for Dn-Matrix entries
  // stores the nodal vector:
  std::map<int,std::map<int,std::pair<int,double> > > refVarWGapMa;
  std::map<int,std::pair<int,double> > newVarWGapMa;

  std::map<int,std::map<int,std::map<int,double> > > refVarWGapLinMa; // stores old derivdn for every node
  int dim = Dim();

  // store the current values for later comparison
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if ((int) cnode->CoData().GetVarWGapSl().size()==0)
      continue;

    // get reference values
    refVarWGapMa[gid]    = cnode->CoData().GetVarWGapMa();

    refVarWGapLinMa[gid] = cnode->CoData().GetVarWGapLinMa();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd=0;fd<dim*snodefullmap->NumMyElements();++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd/dim);
    CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if(!snode)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    int sdof = snode->Dofs()[fd%dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newVarWGapMa.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int)(kcnode->CoData().GetVarWGapMa().size())==0)
        continue;

      typedef std::map<int,std::pair<int,double> >::const_iterator CI;

      newVarWGapMa = kcnode->CoData().GetVarWGapMa();

      // print results (derivatives) to screen
      for (CI p=newVarWGapMa.begin(); p!=newVarWGapMa.end(); ++p)
      {
        if (abs(newVarWGapMa[p->first].second-refVarWGapMa[kgid][p->first].second) > 1e-12)
        {
          double finit = (newVarWGapMa[p->first].second-refVarWGapMa[kgid][p->first].second)/delta;

          double analy = ((refVarWGapLinMa[kgid])[p->first])[sdof];
          double dev = finit - analy;

          // p->first: currently tested dof of slave node p->first/Dim()
          // (p->first)/Dim(): paired slave
          // sdof: currently modified slave dof
          std::cout << "(" << p->first << "," << p->first/dim <<
              "(" << kgid << ")" << "," << sdof << ") :"
              "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
              "   varWGapLinMa=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
              "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
              "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

          if( abs(dev) > 1e-4 )
          {
            std::cout << " ***** WARNING ***** ";
            w++;
          }
          else if( abs(dev) > 1e-5 )
          {
            std::cout << " ***** warning ***** ";
            w++;
          }

          std::cout << std::endl;
        }
      }
    }

    // undo finite difference modification
    if (fd%dim==0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // global loop to apply FD scheme to all MASTER dofs (=dim*nodes)
  for (int fd=0; fd<dim*mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd/dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* mnode = static_cast<CoNode*>(node);

    int mdof = mnode->Dofs()[fd%dim];

    std::cout << "\nDEVIATION FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newVarWGapMa.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      int dim = kcnode->NumDof();

      if ((int)(kcnode->CoData().GetVarWGapMa().size())==0)
        continue;

      typedef std::map<int,std::pair<int,double> >::const_iterator CI;

      newVarWGapMa = kcnode->CoData().GetVarWGapMa();

      // print results (derivatives) to screen
      for (CI p=newVarWGapMa.begin(); p!=newVarWGapMa.end(); ++p)
      {
        if (abs(newVarWGapMa[p->first].second-refVarWGapMa[kgid][p->first].second) > 1e-12)
        {
          double finit = (newVarWGapMa[p->first].second-refVarWGapMa[kgid][p->first].second)/delta;
          double analy = ((refVarWGapLinMa[kgid])[(p->first)])[mdof];
          double dev = finit - analy;

          // p->first: currently tested dof of master node p->first/Dim()
          // (p->first)/Dim(): paired master
          // mdof: currently modified master dof
          std::cout << "(" << p->first << "," << p->first/dim <<
          "(" << kgid << ")" << "," << mdof << ") :"
          "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
          "   varWGapLinMa=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
          "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
          "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

          if( abs(dev) > 1e-4 )
          {
            std::cout << " ***** WARNING ***** ";
            w++;
          }
          else if( abs(dev) > 1e-5 )
          {
            std::cout << " ***** warning ***** ";
            w++;
          }

          std::cout << std::endl;
        }
      }
    }

    // undo finite difference modification
    if (fd%dim==0)
    {
     mnode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
     mnode->xspatial()[1] -= delta;
    }
    else
    {
     mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  Evaluate();
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for AugALin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AugmentedInterface::FDCheckAugALin()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  std::map<int,double> refAugA;
  std::map<int,std::map<int,double> > refAugALin;

  double newAugA = 0.0;

  int dim = Dim();

  // print reference to screen (kappaLin-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if ((int) cnode->CoData().GetAugALin().size()==0)
      continue;

    refAugA[gid]    = cnode->CoData().GetAugA();
    refAugALin[gid] = cnode->CoData().GetAugALin();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd=0;fd<dim*snodefullmap->NumMyElements();++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd/dim);
    CoNode* snode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if(!snode)
      dserror("ERROR: Cannot find slave node with gid %",gid);

    int sdof = snode->Dofs()[fd%dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    RedEvaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAugA = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetAugALin().size()==0)
        continue;

      newAugA = kcnode->CoData().GetAugA();

      // print results (derivatives) to screen
      if (abs(newAugA-refAugA[kgid]) > 1e-12)
      {
        double finit = (newAugA-refAugA[kgid])/delta;

        double analy = refAugALin[kgid][sdof];
        double dev = finit - analy;

        // p->first: currently tested dof of slave node p->first/Dim()
        // (p->first)/Dim(): paired slave
        // sdof: currently modified slave dof
        std::cout << "(" << kgid << "," << sdof << ") :"
            "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
            "   augALin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
            "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
            "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // global loop to apply FD scheme to all MASTER dofs (=dim*nodes)
  for (int fd=0; fd<dim*mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w=0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd/dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find slave node with gid %",gid);
    CoNode* mnode = static_cast<CoNode*>(node);

    int mdof = mnode->Dofs()[fd%dim];

    std::cout << "\nDEVIATION FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%dim==0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd%dim==1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    RedEvaluate();

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAugA = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = static_cast<CoNode*>(knode);

      if ((int) kcnode->CoData().GetAugALin().size()==0)
              continue;

      newAugA = kcnode->CoData().GetAugA();

      // print results (derivatives) to screen

      if (abs(newAugA-refAugA[kgid]) > 1e-12)
      {
        double finit = (newAugA-refAugA[kgid])/delta;
        double analy = refAugALin[kgid][mdof];
        double dev = finit - analy;

        // p->first: currently tested dof of master node p->first/Dim()
        // (p->first)/Dim(): paired master
        // mdof: currently modified master dof
        std::cout << "(" << kgid << "," << mdof << ") :"
        "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific << finit <<
        "   augALin=" << std::setw(14) << std::setprecision(5) << std::scientific << analy <<
        "   DEVIATION= " << std::setw(14) << std::setprecision(5) << std::scientific << dev <<
        "   REL-ERROR [%]= " << std::setw(14) << std::setprecision(5) << std::scientific << abs(dev/finit)*100;

        if( abs(dev) > 1e-4 )
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if( abs(dev) > 1e-5 )
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }

    }

    // undo finite difference modification
    if (fd%dim==0)
    {
     mnode->xspatial()[0] -= delta;
    }
    else if (fd%dim==1)
    {
     mnode->xspatial()[1] -= delta;
    }
    else
    {
     mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** " << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  Evaluate();
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Update of interface related quantities                hiermeier 06/14|
 | during the global FD-check                                           |
 *----------------------------------------------------------------------*/
bool CONTACT::AugmentedInterface::UpdateInterfaces(int gid,
                                                   int dof,
                                                   double delta,
                                                   bool forward)
{
  DRT::Node* node = idiscret_->gNode(gid);
  if (!node) return (false);

  CoNode* cnode = static_cast<CoNode*>(node);

  // change forward step to backward step
  if (!forward) delta *= (-1);
  if (dof==0)
  {
    cnode->xspatial()[0] += delta;
  }
  else if (dof==1)
  {
    cnode->xspatial()[1] += delta;
  }
  else
  {
    cnode->xspatial()[2] += delta;
  }
  Initialize();
  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  Evaluate();
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate();

  return (true);
}
