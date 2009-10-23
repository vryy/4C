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
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"


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
    const std::set<XFEM::FieldEnr>& enrvarset(dofman->getNodeDofSet(node->Id()));
    const vector<int> dof = dis.Dof(node);
    dsassert(dof.size() == enrvarset.size(), "mismatch in length!");
    std::set<XFEM::FieldEnr>::const_iterator enrvar;
    size_t countdof = 0;
    for (enrvar = enrvarset.begin(); enrvar != enrvarset.end(); ++enrvar)
    {
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
void FLD::UTILS::SetupEnrichmentSplit(
        const DRT::Discretization& dis,
        const RCP<XFEM::DofManager> dofman,
        const set<int>& ext_enr_node_gids,
        LINALG::MapExtractor& extractor)
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization for a vector which only
  // contains the unenriched dofs and for one vector which only contains
  // enriched degrees of freedom.
  //
  // The maps are designed assuming that every node has pressure and
  // velocity degrees of freedom --- this won't work for inf-sup stable
  // elements at the moment!
  // -------------------------------------------------------------------

  // Allocate integer vectors which will hold the dof number of the
  // velocity or pressure dofs
  vector<int> normalmapdata;
  vector<int> enrichmapdata;

  const size_t numdof = 4;

  // collect global dofids for velocity and pressure in vectors
  for (int inode=0; inode<dis.NumMyRowNodes(); ++inode)
  {
    const DRT::Node* node = dis.lRowNode(inode);
    const std::set<XFEM::FieldEnr>& enrvarset(dofman->getNodeDofSet(node->Id()));
    const vector<int> dof = dis.Dof(node);
    dsassert(dof.size() == enrvarset.size(), "mismatch in length!");
    size_t countdof = 0;

    if (ext_enr_node_gids.find(node->Id()) == ext_enr_node_gids.end())
    {
      // normal node
      for (size_t i = 0;i<numdof;++i)
      {
        normalmapdata.push_back(dof[countdof]);
        countdof++;
      }
    }
    else
    {
      // enriched node
      for (size_t i = 0;i<numdof;++i)
      {
        enrichmapdata.push_back(dof[countdof]);
        countdof++;
      }
    }
  }

  // the rowmaps are generated according to the pattern provided by
  // the data vectors
  RCP<Epetra_Map> normalrowmap = rcp(new Epetra_Map(-1,
      normalmapdata.size(),&normalmapdata[0],0,
      dis.Comm()));
  RCP<Epetra_Map> enrichrowmap = rcp(new Epetra_Map(-1,
      enrichmapdata.size(),&enrichmapdata[0],0,
      dis.Comm()));

  const Epetra_Map* map = dis.DofRowMap();
  extractor.Setup(*map, enrichrowmap, normalrowmap);
}



// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::UTILS::LiftDrag(
  const DRT::Discretization      & dis         ,
  const Epetra_Vector            & trueresidual,
  const ParameterList            & params      ,
  RCP<map<int,vector<double> > > & liftdragvals
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
    std::map< const int, const std::vector<double>* > ldaxismap;

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
        const int label = ldconds[i]->GetInt("label");

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
        const int label = ldconds[i]->GetInt("label");

        /* get new nodeset for new label OR:
           return pointer to nodeset for known label ... */
        std::set<DRT::Node*>& nodes = ldnodemap[label];

        // centre coordinates to present label
        ldcoordmap[label] = ldconds[i]->Get<vector<double> >("centerCoord");

        // axis of rotation for present label (only needed for 3D)
        if(ldconds[i]->Type() == DRT::Condition::SurfLIFTDRAG)
        {
          ldaxismap[label] = ldconds[i]->Get<vector<double> >("axis");
        }

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
        const std::vector<double>* centerCoordvec = ldcoordmap[label];
        if (centerCoordvec->size() != 3) dserror("axis vector has not length 3");
          LINALG::Matrix<3,1> centerCoord(&((*centerCoordvec)[0]),false);

        // loop all nodes within my set
        for( std::set<DRT::Node*>::const_iterator actnode = nodes.begin(); actnode != nodes.end(); ++actnode)
        {
          const LINALG::Matrix<3,1> x((*actnode)->X(),false); // pointer to nodal coordinates
          const Epetra_BlockMap& rowdofmap = trueresidual.Map();
          const std::vector<int> dof = dis.Dof(*actnode);

          // get nodal forces
          const double fx = trueresidual[rowdofmap.LID(dof[0])];
          const double fy = trueresidual[rowdofmap.LID(dof[1])];
          const double fz = trueresidual[rowdofmap.LID(dof[2])];
          values[0] += fx;
          values[1] += fy;
          values[2] += fz;

          // get distance of point to center point
          LINALG::Matrix<3,1> distances;
          distances.Update(1.0, x, -1.0, centerCoord);

          // get pointer to axis vector (if available)
          const std::vector<double>* axisvecptr = ldaxismap[label];
          // we have to compute the distance of the actual node to the
          // axis of rotation:
          if (axisvecptr)
          {
            if (axisvecptr->size() != 3) dserror("axis vector has not length 3");
            const LINALG::Matrix<3,1> axisvec(&((*axisvecptr)[0]),false);
            double lambda = distances.Dot(axisvec);
            if (axisvec.Norm2() != 0.0)
              lambda = lambda / axisvec.Norm2();
            else
              dserror("zero vector is not a valid direction vector for axis of rotation.");
            distances.Update(lambda,axisvec,1.0);
          }

          // calculate nodal angular momenta
          values[3] += distances(1)*fz-distances(2)*fy;
          values[4] += distances(2)*fx-distances(0)*fz;
          values[5] += distances(0)*fy-distances(1)*fx;

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

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::UTILS::WriteLiftDragToFile(
  const double                     time,
  const int                        step,
  const map<int,vector<double> >&  liftdragvals
  )
{
  // print to file
  std::ostringstream header;
  header << right << std::setw(16) << "Time"
         << right << std::setw(10) << "Step"
         << right << std::setw(10) << "Label"
         << right << std::setw(16) << "F_x"
         << right << std::setw(16) << "F_y"
         << right << std::setw(16) << "F_z";


  for (map<int,vector<double> >::const_iterator liftdragval = liftdragvals.begin(); liftdragval != liftdragvals.end(); ++liftdragval)
  {
    std::ostringstream s;
    s << right << std::setw(16) << scientific << time
      << right << std::setw(10) << scientific << step
      << right << std::setw(10) << scientific << liftdragval->first
      << right << std::setw(16) << scientific << liftdragval->second[0]
      << right << std::setw(16) << scientific << liftdragval->second[1]
      << right << std::setw(16) << scientific << liftdragval->second[2];

    std::ostringstream slabel;
    slabel << std::setw(3) << setfill('0') << liftdragval->first;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".liftdrag_label_"+slabel.str()+".txt";

    if (step <= 1)
    {
      f.open(fname.c_str(),std::fstream::trunc);
      //f << header.str() << endl;
    }
    else
    {
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    }

    f << s.str() << "\n";
    f.close();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int,double> FLD::UTILS::ComputeSurfaceFlowRates(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp
    )
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","calc_flow_rate");

  std::map<int,double> volumeflowratepersurface;

  // get condition
  std::vector< DRT::Condition * >      conds;
  dis.GetCondition ("SurfFlowRate", conds);

  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = conds.begin(); conditer!=conds.end(); ++conditer)
  {
    const DRT::Condition* cond = *conditer;

    const int condID = cond->GetInt("ConditionID");

    // get a vector layout from the discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    const Epetra_Map* dofrowmap = dis.DofRowMap();

    // create vector (+ initialization with zeros)
    Teuchos::RCP<Epetra_Vector> flowrates = LINALG::CreateVector(*dofrowmap,true);

    // call loop over elements
    dis.ClearState();
    dis.SetState("velnp",velnp);
    dis.EvaluateCondition(eleparams,flowrates,"SurfFlowRate",condID);
    dis.ClearState();

    double locflowrate = 0.0;
    for (int i=0; i < dofrowmap->NumMyElements(); i++)
    {
      locflowrate += (*flowrates)[i];
    }

    double flowrate = 0.0;
    dofrowmap->Comm().SumAll(&locflowrate,&flowrate,1);

    volumeflowratepersurface[condID] += flowrate;
  }

  return volumeflowratepersurface;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FLD::UTILS::ComputeSurfaceImpulsRates(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp
    )
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set("action","calc_impuls_rate");

  std::map<int,LINALG::Matrix<3,1> > volumeflowratepersurface;

  // get condition
  std::vector< DRT::Condition * >      conds;
  dis.GetCondition ("SurfImpulsRate", conds);

  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = conds.begin(); conditer!=conds.end(); ++conditer)
  {
    const DRT::Condition* cond = *conditer;

    const int condID = cond->GetInt("ConditionID");

    // create vector (+ initialization with zeros)
    const Epetra_BlockMap mappy = velnp->Map();
    Teuchos::RCP<Epetra_Vector> impulsrates = Teuchos::rcp(new Epetra_Vector(mappy));
    impulsrates->PutScalar(0.0);

    // call loop over elements
    dis.ClearState();
    dis.SetState("velnp",velnp);
    dis.EvaluateCondition(eleparams,impulsrates,"SurfImpulsRate",condID);
    dis.ClearState();
    LINALG::Matrix<3,1> locflowrate(true);
    for (int inode=0; inode < dis.NumMyRowNodes(); inode++)
    {
      const DRT::Node* node = dis.lRowNode(inode);
      static std::vector<int> gdofs(4);
      dis.Dof(node,0,gdofs);
      for (size_t isd=0; isd < 3; isd++)
      {
        locflowrate(isd) += (*impulsrates)[dis.DofColMap()->LID(gdofs[isd])];
//        cout << (*impulsrates)[dis.DofColMap()->LID(gdofs[isd])] << endl;
      }
    }

//    LINALG::Matrix<3,1> flowrate(true);
//    dofrowmap->Comm().SumAll(&locflowrate(0),&flowrate(0),1);
//    dofrowmap->Comm().SumAll(&locflowrate(1),&flowrate(1),1);
//    dofrowmap->Comm().SumAll(&locflowrate(2),&flowrate(2),1);
//    cout << "locflowrate " << locflowrate << endl;
    if (volumeflowratepersurface.find(condID) == volumeflowratepersurface.end())
    {
      LINALG::Matrix<3,1> tmp(true);
      volumeflowratepersurface.insert(make_pair(condID,tmp));
    }
    LINALG::Matrix<3,1> tmp = volumeflowratepersurface[condID];
    tmp += locflowrate;
    volumeflowratepersurface[condID] = tmp;
  }

  return volumeflowratepersurface;
}

#endif /* CCADISCRET       */
