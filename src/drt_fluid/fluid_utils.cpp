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

#include <stdio.h>

#include "fluid_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_fluid_ele/fluid_ele_action.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidSplit(const DRT::Discretization& dis,
                                 int ndim,
                                 LINALG::MultiMapExtractor& extractor)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    std::vector<int> dof = dis.Dof(0,node);
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

  std::vector<Teuchos::RCP<const Epetra_Map> > maps( 2 );
  maps[0] = otherdofmap;
  maps[1] = conddofmap;
  extractor.Setup(*dis.DofRowMap(),maps);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidSplit(const DRT::Discretization& dis,
                                 int fluid_ndof,
                                 int pres_ndof,
                                 LINALG::MultiMapExtractor& extractor)
{
  unsigned fp_dim = static_cast<unsigned>(fluid_ndof + pres_ndof);

  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis.NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis.lRowNode(i);

    std::vector<int> dof = dis.Dof(node);

    if( (dof.size() % fp_dim) != 0) dserror("Fluid-Pres-Split is not unique! mismatch between number of dofs and fluid/pres dim");

    for (unsigned j=0; j<dof.size(); ++j)
    {
      // test for dof position
      if (j%fp_dim < static_cast<unsigned>(fluid_ndof))
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

  std::vector<Teuchos::RCP<const Epetra_Map> > maps( 2 );
  maps[0] = otherdofmap;
  maps[1] = conddofmap;
  extractor.Setup(*dis.DofRowMap(),maps);
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


//----------------------------------------------------------------------*/
//----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidFluidVelPresSplit(const DRT::Discretization& fluiddis,int ndim,
                                             const DRT::Discretization& alefluiddis,
                                             LINALG::MapExtractor& extractor,
                                             Teuchos::RCP<Epetra_Map> fullmap)
{
  std::set<int> veldofset;
  std::set<int> presdofset;

  // for fluid elements
  int numfluidrownodes = fluiddis.NumMyRowNodes();
  for (int i=0; i<numfluidrownodes; ++i)
  {
    DRT::Node* fluidnode = fluiddis.lRowNode(i);

    std::vector<int> fluiddof = fluiddis.Dof(fluidnode);
    for (unsigned j=0; j<fluiddof.size(); ++j)
    {
      // test for dof position
      if (j<static_cast<unsigned>(ndim))
      {
        veldofset.insert(fluiddof[j]);
      }
      else
      {
        presdofset.insert(fluiddof[j]);
      }
     }
   }

  // for ale_fluid elements
  int numalefluidrownodes = alefluiddis.NumMyRowNodes();
  for (int i=0; i<numalefluidrownodes; ++i)
  {
    DRT::Node* alefluidnode = alefluiddis.lRowNode(i);

    std::vector<int> alefluiddof = alefluiddis.Dof(alefluidnode);
    for (unsigned j=0; j<alefluiddof.size(); ++j)
    {
      // test for dof position
      if (j<static_cast<unsigned>(ndim))
      {
        veldofset.insert(alefluiddof[j]);
      }
      else
      {
        presdofset.insert(alefluiddof[j]);
      }
     }
   }

  std::vector<int> veldofmapvec;
  veldofmapvec.reserve(veldofset.size());
  veldofmapvec.assign(veldofset.begin(), veldofset.end());
  veldofset.clear();
  RCP<Epetra_Map> velrowmap = rcp(new Epetra_Map(-1,
                                  veldofmapvec.size(),&veldofmapvec[0],0,
                                  fluiddis.Comm()));
  veldofmapvec.clear();

  std::vector<int> presdofmapvec;
  presdofmapvec.reserve(presdofset.size());
  presdofmapvec.assign(presdofset.begin(), presdofset.end());
  presdofset.clear();
  RCP<Epetra_Map> presrowmap = rcp(new Epetra_Map(-1,
                                  presdofmapvec.size(),&presdofmapvec[0],0,
                                  alefluiddis.Comm()));
  extractor.Setup(*fullmap, presrowmap, velrowmap);
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

  if (liftdrag == INPAR::FLUID::liftdrag_none); // do nothing, we don't want lift & drag

  if (liftdrag == INPAR::FLUID::liftdrag_nodeforce)
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
          const unsigned numdof = dof.size();

          // extension of lift and drag function to enriched two-phase flow problems
          int dofvelx=0;
          int dofvely=0;
          int dofvelz=0;

          // numdof == 4 : standard fluid or node not enriched
          // numdof == 5 : only pressure field enriched
          if(numdof==4 or numdof==5)
          {
            dofvelx=dof[0];
            dofvely=dof[1];
            dofvelz=dof[2];
          }
          // numdof == 8 : pressure and velocity field enriched
          // numdof == 9 : only velocity field enriched (rare)
          else if(numdof==8 or numdof==7)
          {
            dofvelx=dof[0];
            dofvely=dof[2];
            dofvelz=dof[4];
          }
         else
          dserror("Number of dofs not known for lift&drag!");

          // get nodal forces
          const double fx = trueresidual[rowdofmap.LID(dofvelx)];
          const double fy = trueresidual[rowdofmap.LID(dofvely)];
          const double fz = trueresidual[rowdofmap.LID(dofvelz)];
          values[0] += fx;
          values[1] += fy;
          values[2] += fz;

          // get moment
          double mx = 0.0;
          double my = 0.0;
          double mz = 0.0;

          // get vector of point to center point
          LINALG::Matrix<3,1> distances;
          distances.Update(1.0, x, -1.0, centerCoord);

          // calculate nodal angular moment with respect to global coordinate system
          const double mx_gc = distances(1)*fz-distances(2)*fy;
          const double my_gc = distances(2)*fx-distances(0)*fz;
          const double mz_gc = distances(0)*fy-distances(1)*fx;

          // get pointer to axis vector (if available)
          const std::vector<double>* axisvecptr = ldaxismap[label];
          if (axisvecptr)
          {
            if (axisvecptr->size() != 3) dserror("axis vector has not length 3");
            LINALG::Matrix<3,1> axisvec(&((*axisvecptr)[0]),false);
            double norm = 0.0;
            if (axisvec.Norm2() != 0.0)
            {
              norm = axisvec.Norm2();
              // normed axis vector
              axisvec.Scale(1.0/norm);
            }
            else
              dserror("norm==0.0!");
            // projection of moment on given axis
            double mdir = mx_gc*axisvec(0,0) + my_gc*axisvec(1,0)+ mz_gc*axisvec(2,0);

            mx = mdir*axisvec(0,0);
            my = mdir*axisvec(1,0);
            mz = mdir*axisvec(2,0);
          }
          else
          {
            mx = mx_gc;
            my = my_gc;
            mz = mz_gc;
          }

          values[3] += mx;
          values[4] += my;
          values[5] += mz;

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
std::map<int,double> FLD::UTILS::ComputeFlowRates(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp,
    const string                   condstring)
{
  ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action",FLD::calc_flowrate);

  // note that the flowrate is not yet divided by the area
  std::map<int,double> volumeflowrateperline;

  // get condition
  std::vector< DRT::Condition* > conds;
  dis.GetCondition (condstring, conds);

  // each condition is on every proc , but might not have condition elements there
  for(vector<DRT::Condition*>::const_iterator conditer = conds.begin(); conditer!=conds.end(); ++conditer)
  {
    const DRT::Condition* cond = *conditer;
    const int condID = cond->GetInt("ConditionID");

    // get a vector layout from the discretization to construct matching
    // vectors and matrices local <-> global dof numbering
    const Epetra_Map* dofrowmap = dis.DofRowMap();

    // create vector (+ initialization with zeros)
    Teuchos::RCP<Epetra_Vector> flowrates = LINALG::CreateVector(*dofrowmap,true);

    // call loop over elements
    dis.ClearState();
    dis.SetState("velnp",velnp);

    dis.EvaluateCondition(eleparams,flowrates,condstring,condID);
    dis.ClearState();

    double local_flowrate = 0.0;
    for (int i=0; i < dofrowmap->NumMyElements(); i++)
    {
      local_flowrate +=((*flowrates)[i]);
    }

    double flowrate = 0.0;
    dofrowmap->Comm().SumAll(&local_flowrate,&flowrate,1);

    //if(dofrowmap->Comm().MyPID()==0)
    //	cout << "gobal flow rate = " << flowrate << "\t condition ID = " << condID << endl;

    //ATTENTION: new definition: outflow is positive and inflow is negative
    volumeflowrateperline[condID] = flowrate;
  }
  return volumeflowrateperline;
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::WriteFlowRatesToFile(
  const double                     time,
  const int                        step,
  const std::map<int,double>&      flowrates
  )
{
  if (flowrates.empty())
    dserror("flowratevector is empty");

  // print to file
  std::ostringstream header;
  header << right << std::setw(16) << "Time"
         << right << std::setw(10) << "Step"
         << right << std::setw(10) << "ID"
         << right << std::setw(16) << "Flowrate";

  for(map<int,double >::const_iterator flowrate = flowrates.begin(); flowrate != flowrates.end(); ++flowrate)
  {
    std::ostringstream s;
    s << right << std::setw(16) << scientific << time
      << right << std::setw(10) << scientific << step
      << right << std::setw(10) << scientific << flowrate->first
      << right << std::setw(16) << scientific << flowrate->second;

    std::ostringstream slabel;
    slabel << std::setw(3) << setfill('0') << flowrate->first;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".flowrate_ID_"+slabel.str()+".txt";

    if (step <= 1)
      f.open(fname.c_str(),std::fstream::trunc); //f << header.str() << endl;
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << s.str() << "\n";
    f.close();
  }
}

