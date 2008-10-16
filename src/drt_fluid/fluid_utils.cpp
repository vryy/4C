/*!----------------------------------------------------------------------
\file fluid_utils.cpp
\brief utility functions for fluid problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <stdio.h>

#include "fluid_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidSplit(const DRT::Discretization& dis,
                                 int ndim,
                                 LINALG::MapExtractor& extractor)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    std::vector<int> dof = dis.Dof(node);
    for (unsigned j=0; j<dof.size(); ++j)
    {
      // test for dof position
      if (j<static_cast<unsigned>(ndim))
      {
        otherdofset.insert(dof[j]);
      }
      else
      {
        conddofset.insert(dof[j]);
      }
    }
  }

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(conddofset.size());
  conddofmapvec.assign(conddofset.begin(), conddofset.end());
  conddofset.clear();
  Teuchos::RCP<Epetra_Map> conddofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                conddofmapvec.size(),
                                &conddofmapvec[0],
                                0,
                                dis.Comm()));
  conddofmapvec.clear();

  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(otherdofset.size());
  otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
  otherdofset.clear();
  Teuchos::RCP<Epetra_Map> otherdofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                otherdofmapvec.size(),
                                &otherdofmapvec[0],
                                0,
                                dis.Comm()));
  otherdofmapvec.clear();

  extractor.Setup(*dis.DofRowMap(),conddofmap,otherdofmap);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidSplit(const DRT::Discretization& dis,
                                 const DRT::DofSet& dofset,
                                 int ndim,
                                 LINALG::MapExtractor& extractor)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    std::vector<int> dof = dofset.Dof(node);
    for (unsigned j=0; j<dof.size(); ++j)
    {
      // test for dof position
      if (j<static_cast<unsigned>(ndim))
      {
        otherdofset.insert(dof[j]);
      }
      else
      {
        conddofset.insert(dof[j]);
      }
    }
  }

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(conddofset.size());
  conddofmapvec.assign(conddofset.begin(), conddofset.end());
  conddofset.clear();
  Teuchos::RCP<Epetra_Map> conddofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                conddofmapvec.size(),
                                &conddofmapvec[0],
                                0,
                                dis.Comm()));
  conddofmapvec.clear();

  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(otherdofset.size());
  otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
  otherdofset.clear();
  Teuchos::RCP<Epetra_Map> otherdofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                otherdofmapvec.size(),
                                &otherdofmapvec[0],
                                0,
                                dis.Comm()));
  otherdofmapvec.clear();

  extractor.Setup(*dofset.DofRowMap(),conddofmap,otherdofmap);
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::UTILS::SetupXFluidSplit(
        const DRT::Discretization& dis,
        const RCP<XFEM::DofManager> dofman,
        LINALG::MapExtractor& extractor)
{
    // -------------------------------------------------------------------
    // get a vector layout from the discretization for a vector which only
    // contains the velocity dofs and for one vector which only contains
    // pressure degrees of freedom.
    //
    // The maps are designed assuming that every node has pressure and
    // velocity degrees of freedom --- this won't work for inf-sup stable
    // elements at the moment!
    // -------------------------------------------------------------------

    // Allocate integer vectors which will hold the dof number of the
    // velocity or pressure dofs
    vector<int> velmapdata;
    vector<int> premapdata;

    // collect global dofids for velocity and pressure in vectors
    for (int i=0; i<dis.NumMyRowNodes(); ++i) {
        const DRT::Node* node = dis.lRowNode(i);
        const std::set<XFEM::FieldEnr> enrvarset = dofman->getNodeDofSet(node->Id());
        const vector<int> dof = dis.Dof(node);
        dsassert(dof.size() == enrvarset.size(), "mismatch in length!");
        std::set<XFEM::FieldEnr>::const_iterator enrvar;
        unsigned int countdof = 0;
        for (enrvar = enrvarset.begin(); enrvar != enrvarset.end(); ++enrvar) {
            switch (enrvar->getField()) {
                case XFEM::PHYSICS::Velx:
                case XFEM::PHYSICS::Vely:
                case XFEM::PHYSICS::Velz:
                    velmapdata.push_back(dof[countdof]);
                    break;
                case XFEM::PHYSICS::Pres:
                    premapdata.push_back(dof[countdof]);
                    break;
                default:
                    break;
            }
            countdof++;
        }
    }

    // the rowmaps are generated according to the pattern provided by
    // the data vectors
    RCP<Epetra_Map> velrowmap = rcp(new Epetra_Map(-1,
            velmapdata.size(),&velmapdata[0],0,
            dis.Comm()));
    RCP<Epetra_Map> prerowmap = rcp(new Epetra_Map(-1,
            premapdata.size(),&premapdata[0],0,
            dis.Comm()));

    const Epetra_Map* map = dis.DofRowMap();
    extractor.Setup(*map, prerowmap, velrowmap);
}



// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::UTILS::LiftDrag(
  const DRT::Discretization&     dis         ,      
  const Epetra_Vector&           trueresidual,    
  ParameterList&                 params      ,
  RCP<map<int,vector<double> > > liftdragvals
  )
{
  
  const int liftdrag = params.get<int>("liftdrag");
  
  if (liftdrag == 0); // do nothing, we don't want lift & drag
  if (liftdrag == 1)
  dserror("we do not support lift&drag calculation by stresses anymore; use nodal forces!");

  if (liftdrag == 2)
  {
    int myrank=dis.Comm().MyPID();

    std::map< const int, std::set<DRT::Node* > > ldnodemap;
    std::map< const int, const std::vector<double>* > ldcoordmap;

    // allocate and initialise LiftDrag conditions
    std::vector<DRT::Condition*> ldconds;
    dis.GetCondition("LIFTDRAG",ldconds);

    // space dimension of the problem
    const int ndim = params.get<int>("number of velocity degrees of freedom");

    // there is an L&D condition if it has a size
    if( ldconds.size() )
    {

      // vector with lift&drag forces after communication
      liftdragvals = rcp(new std::map<int,std::vector<double> >);

      for( unsigned i=0; i<ldconds.size(); ++i) // loop L&D conditions (i.e. lines in .dat file)
      {
        /* get label of present LiftDrag condition  */
        const unsigned int label = ldconds[i]->Getint("label");

        ((*liftdragvals)).insert(pair<int,vector<double> >(label,vector<double> (6,0.0)));
      }

      // prepare output
      if (myrank==0)
      {
        cout << "Lift and drag calculation:" << "\n";
        if (ndim == 2)
        {
          cout << "lift'n'drag Id      F_x             F_y             M_z :" << "\n";
        }
        if (ndim == 3)
        {
          cout << "lift'n'drag Id      F_x             F_y             F_z           ";
          cout << "M_x             M_y             M_z :" << "\n";
        }
      }

      // sort data
      for( unsigned i=0; i<ldconds.size(); ++i) // loop L&D conditions (i.e. lines in .dat file)
      {
        /* get label of present LiftDrag condition  */
        const unsigned int label = ldconds[i]->Getint("label");

        /* get new nodeset for new label OR:
           return pointer to nodeset for known label ... */
        std::set<DRT::Node*>& nodes = ldnodemap[label];

        // centre coordinates to present label
        ldcoordmap[label] = ldconds[i]->Get<vector<double> >("centerCoord");

        /* get pointer to its nodal Ids*/
        const vector<int>* ids = ldconds[i]->Get<vector<int> >("Node Ids");

        /* put all nodes belonging to the L&D line or surface into 'nodes' which are
           associated with the present label */
        for (unsigned j=0; j<ids->size(); ++j)
        {
          // give me present node Id
          const int node_id = (*ids)[j];
          // put it into nodeset of actual label if node is new and mine
          if( dis.HaveGlobalNode(node_id) && dis.gNode(node_id)->Owner()==myrank )
	  nodes.insert(dis.gNode(node_id));
        }
      } // end loop over conditions


      // now step the label map
      for( std::map< const int, std::set<DRT::Node*> >::const_iterator labelit = ldnodemap.begin();
           labelit != ldnodemap.end(); ++labelit )
      {
        const std::set<DRT::Node*>& nodes = labelit->second; // pointer to nodeset of present label
        const int label = labelit->first;                    // the present label
        std::vector<double> values(6,0.0);             // vector with lift&drag forces

        // get also pointer to centre coordinates
        const std::vector<double>* centerCoord = ldcoordmap[label];

        // loop all nodes within my set
        for( std::set<DRT::Node*>::const_iterator actnode = nodes.begin(); actnode != nodes.end(); ++actnode)
        {
          const double* x = (*actnode)->X(); // pointer to nodal coordinates
          const Epetra_BlockMap& rowdofmap = trueresidual.Map();
          const std::vector<int> dof = dis.Dof(*actnode);

          std::vector<double> distances (3);
          for (unsigned j=0; j<3; ++j)
          {
            distances[j]= x[j]-(*centerCoord)[j];
          }
          // get nodal forces
          const double fx = trueresidual[rowdofmap.LID(dof[0])];
          const double fy = trueresidual[rowdofmap.LID(dof[1])];
          const double fz = trueresidual[rowdofmap.LID(dof[2])];
          values[0] += fx;
          values[1] += fy;
          values[2] += fz;

          // calculate nodal angular momenta
          values[3] += distances[1]*fz-distances[2]*fy;
          values[4] += distances[2]*fx-distances[0]*fz;
          values[5] += distances[0]*fy-distances[1]*fx;
        } // end: loop over nodes

        // care for the fact that we are (most likely) parallel
        trueresidual.Comm().SumAll (&(values[0]), &(((*liftdragvals)[label])[0]), 6);

        // do the output
        if (myrank==0)
        {
          if (ndim == 2)
          {
            cout << "     " << label << "         ";
            cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[5];
            cout << "\n";
          }
          if (ndim == 3)
          {
            cout << "     " << label << "         ";
            cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[2] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[3] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[4] << "    ";
            cout << std::scientific << ((*liftdragvals)[label])[5];
            cout << "\n";
          }
        }
      } // end: loop over L&D labels
      if (myrank== 0)
      {
        cout << "\n";
      }
    }
  }
  
  return;
}


#endif /* CCADISCRET       */
