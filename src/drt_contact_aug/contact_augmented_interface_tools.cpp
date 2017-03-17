/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_interface_tools.cpp

\brief Tools for the augmented contact interface evaluation.

\level 2

\maintainer Michael Hiermeier

\date Apr 11, 2014

*/
/*---------------------------------------------------------------------*/
#include "contact_augmented_interface.H"
#include "../drt_contact/contact_integrator.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_integrator.H"
#include "../drt_mortar/mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"

#include "../headers/pairedvector.H"


/*----------------------------------------------------------------------*
 | Finite difference check for KappaLin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FDCheckKappaLin(
    CONTACT::ParamsInterface& cparams)
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
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
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if ((int) cnode->AugData().GetKappaLin().size()==0)
      continue;

    refKappa[gid]    = cnode->AugData().GetKappa();
    refKappaLin[gid] = cnode->AugData().GetKappaLin();
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
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newKappa = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetKappaLin().size()==0)
        continue;

      newKappa = kcnode->AugData().GetKappa();

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
    CoNode* mnode = dynamic_cast<CoNode*>(node);

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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newKappa = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetKappaLin().size()==0)
              continue;

      newKappa = kcnode->AugData().GetKappa();

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
  Evaluate(cparams_ptr);
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for AWGapLin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FDCheckAWGapLin(
    CONTACT::ParamsInterface& cparams)
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
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
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if ((int) cnode->AugData().GetAWGapLin().size()==0)
      dserror("Linearization map shouldn't be empty for active nodes!");

    double wg = cnode->AugData().GetWGap();
    double kappa = cnode->AugData().GetKappa();
    if (wg==1.0e12 or wg==0.0)
      dserror("The weighted gap is at the active node % equal 0.0 or 1.0e12!",gid);
    if (kappa==1.0e12 or kappa==0.0)
      dserror("The scaling factor kappa is at the active node % equal 0.0 or 1.0e12!",gid);

    refAWGap[gid]    = wg/kappa;
    refAWGapLin[gid] = cnode->AugData().GetAWGapLin();
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
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
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
    Evaluate(cparams_ptr);
    WGap();
    AWGapLin();

    // compute finite difference derivative
    for (int k=0; k<activenodes_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAWGap = 0.0;
      int kgid = activenodes_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetAWGapLin().size()==0)
        dserror("Linearization map shouldn't be empty for active nodes!");

      double wg = kcnode->AugData().GetWGap();
      double kappa = kcnode->AugData().GetKappa();
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
    CoNode* mnode = dynamic_cast<CoNode*>(node);

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
    Evaluate(cparams_ptr);
    WGap();
    AWGapLin();

    // compute finite difference derivative
    for (int k=0; k<activenodes_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAWGap = 0.0;
      int kgid = activenodes_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetAWGapLin().size()==0)
        dserror("Linearization map shouldn't be empty for active nodes!");

      double wg = kcnode->AugData().GetWGap();
      double kappa = kcnode->AugData().GetKappa();
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
  Evaluate(cparams_ptr);
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for VarWGapLinSl              hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FDCheckVarWGapLinSl(
    CONTACT::ParamsInterface& cparams)
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
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
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if ((int) cnode->AugData().GetVarWGapSl().size()==0)
      continue;

    // get reference values
    refVarWGapSl[gid].insert(cnode->AugData().GetVarWGapSl().begin(),cnode->AugData().GetVarWGapSl().end());

    refVarWGapLinSl[gid] = cnode->AugData().GetVarWGapLinSl();
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
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
//      newVarWGapSl.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int)(kcnode->AugData().GetVarWGapSl().size())==0)
        continue;

      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

      GEN::pairedvector<int,std::pair<int,double> >& newVarWGapSl = kcnode->AugData().GetVarWGapSl();

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
    CoNode* mnode = dynamic_cast<CoNode*>(node);

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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
//      newVarWGapSl.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      int dim = kcnode->NumDof();

      if ((int)(kcnode->AugData().GetVarWGapLinSl().size())==0)
        continue;

      typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

      GEN::pairedvector<int,std::pair<int,double> >& newVarWGapSl = kcnode->AugData().GetVarWGapSl();

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
  Evaluate(cparams_ptr);
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for VarWGapLinMa              hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FDCheckVarWGapLinMa(
    CONTACT::ParamsInterface& cparams)
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
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
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if ((int) cnode->AugData().GetVarWGapSl().size()==0)
      continue;

    // get reference values
    refVarWGapMa[gid]    = cnode->AugData().GetVarWGapMa();

    refVarWGapLinMa[gid] = cnode->AugData().GetVarWGapLinMa();
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
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newVarWGapMa.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int)(kcnode->AugData().GetVarWGapMa().size())==0)
        continue;

      typedef std::map<int,std::pair<int,double> >::const_iterator CI;

      newVarWGapMa = kcnode->AugData().GetVarWGapMa();

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
    CoNode* mnode = dynamic_cast<CoNode*>(node);

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
    Evaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newVarWGapMa.clear();
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      int dim = kcnode->NumDof();

      if ((int)(kcnode->AugData().GetVarWGapMa().size())==0)
        continue;

      typedef std::map<int,std::pair<int,double> >::const_iterator CI;

      newVarWGapMa = kcnode->AugData().GetVarWGapMa();

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
  Evaluate(cparams_ptr);
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for AugALin                  hiermeier 05/14|
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::FDCheckAugALin(
    CONTACT::ParamsInterface& cparams)
{
  // get out of here if not participating in interface
  if (!lComm()) return;
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
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
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if ((int) cnode->AugData().GetAugALin().size()==0)
      continue;

    std::map<int,double> augALinMap( cnode->AugData().GetAugALin().begin(),
        cnode->AugData().GetAugALin().end() );
    refAugA[gid]    = cnode->AugData().GetAugA();
    refAugALin[gid] = augALinMap;
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
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
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
    RedEvaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAugA = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetAugALin().size()==0)
        continue;

      newAugA = kcnode->AugData().GetAugA();

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
    CoNode* mnode = dynamic_cast<CoNode*>(node);

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
    RedEvaluate(cparams_ptr);

    // compute finite difference derivative
    for (int k=0; k<snoderowmap_->NumMyElements(); ++k)
    {
      // clear the calculated new r.h.s. map
      newAugA = 0.0;
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode)
        dserror("ERROR: Cannot find node with gid %",kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if ((int) kcnode->AugData().GetAugALin().size()==0)
              continue;

      newAugA = kcnode->AugData().GetAugA();

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
  Evaluate(cparams_ptr);
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return;
}

/*----------------------------------------------------------------------*
 | Update of interface related quantities                hiermeier 06/14|
 | during the global FD-check                                           |
 *----------------------------------------------------------------------*/
bool CONTACT::AUG::Interface::UpdateInterfaces(
    int gid,
    int dof,
    double delta,
    bool forward,
    CONTACT::ParamsInterface& cparams)
{
  DRT::Node* node = idiscret_->gNode(gid);
  if (!node) return (false);

  CoNode* cnode = dynamic_cast<CoNode*>(node);

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
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr =
      Teuchos::rcpFromRef(cparams);
  // evaluate averaged weighted gap
  Evaluate(cparams_ptr);
  WGap();
  AWGapLin();
  // evaluate remaining entities and linearization
  RedEvaluate(cparams_ptr);

  return (true);
}
