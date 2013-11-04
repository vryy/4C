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
  RCP<Epetra_Map> velrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                  veldofmapvec.size(),&veldofmapvec[0],0,
                                  fluiddis.Comm()));
  veldofmapvec.clear();

  std::vector<int> presdofmapvec;
  presdofmapvec.reserve(presdofset.size());
  presdofmapvec.assign(presdofset.begin(), presdofset.end());
  presdofset.clear();
  RCP<Epetra_Map> presrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                  presdofmapvec.size(),&presdofmapvec[0],0,
                                  alefluiddis.Comm()));
  extractor.Setup(*fullmap, presrowmap, velrowmap);
}



// -------------------------------------------------------------------
// compute forces and moments                          rasthofer 08/13
// -------------------------------------------------------------------
void FLD::UTILS::LiftDrag(
  const DRT::Discretization      & dis         ,
  const Epetra_Vector            & trueresidual,
  const Teuchos::ParameterList   & params      ,
  Teuchos::RCP<std::map<int,std::vector<double> > > & liftdragvals
  )
{
  const int liftdrag = params.get<int>("liftdrag");

  if (liftdrag == INPAR::FLUID::liftdrag_nodeforce)
  {
    int myrank=dis.Comm().MyPID();

    std::map< const int, std::set<DRT::Node* > > ldnodemap;
    std::map< const int, const std::vector<double>* > ldcoordmap;
    std::map< const int, const std::vector<double>* > ldaxismap;
    bool axis_for_moment = false;

    // allocate and initialise LiftDrag conditions
    std::vector<DRT::Condition*> ldconds;
    dis.GetCondition("LIFTDRAG",ldconds);

    // space dimension of the problem
    const int ndim = params.get<int>("number of velocity degrees of freedom");

    // there is an L&D condition if it has a size
    if( ldconds.size() )
    {

      // vector with lift&drag forces after communication
      liftdragvals = Teuchos::rcp(new std::map<int,std::vector<double> >);

      for( unsigned i=0; i<ldconds.size(); ++i) // loop L&D conditions (i.e. lines in .dat file)
      {
        /* get label of present LiftDrag condition  */
        const int label = ldconds[i]->GetInt("label");

        ((*liftdragvals)).insert(std::pair<int,std::vector<double> >(label,std::vector<double> (6,0.0)));
      }

      // prepare output
      if (myrank==0)
      {
        std::cout << "Lift and drag calculation:" << "\n";
        if (ndim == 2)
        {
          std::cout << "lift'n'drag Id      F_x             F_y             M_z :" << "\n";
        }
        if (ndim == 3)
        {
          std::cout << "lift'n'drag Id      F_x             F_y             F_z           ";
          std::cout << "M_x             M_y             M_z :" << "\n";
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

        // center coordinates to present label
        ldcoordmap[label] = ldconds[i]->Get<std::vector<double> >("centerCoord");

        // axis of rotation for present label (only needed for 3D)
        if(ldconds[i]->Type() == DRT::Condition::SurfLIFTDRAG)
        {
          ldaxismap[label] = ldconds[i]->Get<std::vector<double> >("axis");
          // get pointer to axis vector (if available)
          const std::vector<double>* axisvecptr = ldaxismap[label];
          if (axisvecptr->size() != 3) dserror("axis vector has not length 3");
          LINALG::Matrix<3,1> axisvec(&((*axisvecptr)[0]),false);
          if (axisvec.Norm2() > 1.0e-9) axis_for_moment=true; // axis has been set
        }

        /* get pointer to its nodal Ids*/
        const std::vector<int>* ids = ldconds[i]->Get<std::vector<int> >("Node Ids");

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
        std::vector<double> myforces(3,0.0);                 // vector with lift&drag forces
        std::vector<double> mymoments(3,0.0);                // vector with lift&drag moments

        // get also pointer to center coordinates
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
          LINALG::Matrix<3,1> actforces (true);
          for (int idim=0; idim<ndim; idim++)
          {
            actforces(idim,0) = trueresidual[rowdofmap.LID(dof[idim])];
            myforces[idim] += trueresidual[rowdofmap.LID(dof[idim])];
          }
          // z-component remains zero for ndim=2

          // get moment
          LINALG::Matrix<3,1> actmoments (true);
          // get vector of point to center point
          LINALG::Matrix<3,1> distances;
          distances.Update(1.0, x, -1.0, centerCoord);

          // calculate nodal angular moment with respect to global coordinate system
          LINALG::Matrix<3,1> actmoment_gc (true);
          actmoment_gc(0,0) = distances(1)*actforces(2,0)-distances(2)*actforces(1,0); // zero for 2D
          actmoment_gc(1,0) = distances(2)*actforces(0,0)-distances(0)*actforces(2,0); // zero for 2D
          actmoment_gc(2,0) = distances(0)*actforces(1,0)-distances(1)*actforces(0,0);

          if (axis_for_moment)
          {
            const std::vector<double>* axisvecptr = ldaxismap[label];
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
            double mdir = actmoment_gc.Dot(axisvec);

           actmoments(2,0) = mdir;

          }
          else
          {
            for (int idim=0; idim<3; idim++)
              actmoments(idim,0) = actmoment_gc(idim,0);
          }

          for (int idim=0; idim<3; idim++)
            mymoments[idim] += actmoments(idim,0);

        } // end: loop over nodes

        // care for the fact that we are (most likely) parallel
        trueresidual.Comm().SumAll (&(myforces[0]), &(((*liftdragvals)[label])[0]), 3);
        trueresidual.Comm().SumAll (&(mymoments[0]), &(((*liftdragvals)[label])[3]), 3);

        // do the output
        if (myrank==0)
        {
          if (ndim == 2)
          {
            std::cout << "     " << label << "         ";
            std::cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[5];
            std::cout << "\n";
          }
          if (ndim == 3)
          {
            std::cout << "     " << label << "         ";
            std::cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[2] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[3] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[4] << "    ";
            std::cout << std::scientific << ((*liftdragvals)[label])[5];
            std::cout << "\n";
          }
        }
      } // end: loop over L&D labels
      if (myrank== 0)
      {
        std::cout << "\n";
      }
    }
  }

  return;
}

// -------------------------------------------------------------------
// write forces and moments to file                    rasthofer 08/13
// -------------------------------------------------------------------
void FLD::UTILS::WriteLiftDragToFile(
  const double                     time,
  const int                        step,
  const std::map<int,std::vector<double> >&  liftdragvals
  )
{
  // print to file
  std::ostringstream header;
  header << std::right << std::setw(16) << "Time"
         << std::right << std::setw(10) << "Step"
         << std::right << std::setw(10) << "Label"
         << std::right << std::setw(16) << "F_x"
         << std::right << std::setw(16) << "F_y"
         << std::right << std::setw(16) << "F_z"
         << std::right << std::setw(16) << "M_x"
         << std::right << std::setw(16) << "M_y"
         << std::right << std::setw(16) << "M_z";


  for (std::map<int,std::vector<double> >::const_iterator liftdragval = liftdragvals.begin(); liftdragval != liftdragvals.end(); ++liftdragval)
  {
    std::ostringstream s;
    s << std::right << std::setw(16) << std::scientific << time
      << std::right << std::setw(10) << std::scientific << step
      << std::right << std::setw(10) << std::scientific << liftdragval->first
      << std::right << std::setw(16) << std::scientific << liftdragval->second[0]
      << std::right << std::setw(16) << std::scientific << liftdragval->second[1]
      << std::right << std::setw(16) << std::scientific << liftdragval->second[2]
      << std::right << std::setw(16) << std::scientific << liftdragval->second[3]
      << std::right << std::setw(16) << std::scientific << liftdragval->second[4]
      << std::right << std::setw(16) << std::scientific << liftdragval->second[5];

    std::ostringstream slabel;
    slabel << std::setw(3) << std::setfill('0') << liftdragval->first;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".liftdrag_label_"+slabel.str()+".txt";

    if (step <= 1)
    {
      f.open(fname.c_str(),std::fstream::trunc);
      f << header.str() << std::endl;
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
    const std::string              condstring)
{
  return ComputeFlowRates(dis,velnp,Teuchos::null,Teuchos::null,condstring);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int,double> FLD::UTILS::ComputeFlowRates(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp,
    const RCP<Epetra_Vector>       gridv,
    const RCP<Epetra_Vector>       dispnp,
    const std::string              condstring)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action",FLD::calc_flowrate);

  // note that the flowrate is not yet divided by the area
  std::map<int,double> volumeflowrateperline;

  // get condition
  std::vector< DRT::Condition* > conds;
  dis.GetCondition (condstring, conds);

  // each condition is on every proc , but might not have condition elements there
  for(std::vector<DRT::Condition*>::const_iterator conditer = conds.begin(); conditer!=conds.end(); ++conditer)
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
    if(dispnp != Teuchos::null)
      dis.SetState("dispnp",dispnp);
    if(gridv != Teuchos::null)
      dis.SetState("gridv",gridv);

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
    //	std::cout << "gobal flow rate = " << flowrate << "\t condition ID = " << condID << std::endl;

    //ATTENTION: new definition: outflow is positive and inflow is negative
    volumeflowrateperline[condID] = flowrate;
  }
  return volumeflowrateperline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int,double> FLD::UTILS::ComputeVolume(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp,
    const RCP<Epetra_Vector>       gridv,
    const RCP<Epetra_Vector>       dispnp)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action",FLD::calc_volume);

  std::map<int,double> volumeperline;

  // get a vector layout from the discretization to construct matching
  // vectors and matrices local <-> global dof numbering
  const Epetra_Map* dofrowmap = dis.DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> volumes = LINALG::CreateVector(*dofrowmap,true);

  // call loop over elements
  dis.ClearState();
  dis.SetState("velnp",velnp);
  if(dispnp != Teuchos::null)
    dis.SetState("dispnp",dispnp);
  if(gridv != Teuchos::null)
    dis.SetState("gridv",gridv);

  dis.Evaluate(eleparams,Teuchos::null,volumes);
  dis.ClearState();

  double local_volume = 0.0;
  for (int i=0; i < dofrowmap->NumMyElements(); i++)
  {
    local_volume +=((*volumes)[i]);
  }

  double volume = 0.0;
  dofrowmap->Comm().SumAll(&local_volume,&volume,1);

  volumeperline[0] = volume;

  return volumeperline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FLD::UTILS::ComputeSurfaceImpulsRates(
    DRT::Discretization&           dis  ,
    const RCP<Epetra_Vector>       velnp
    )
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set("action","calc_impuls_rate");

  std::map<int,LINALG::Matrix<3,1> > volumeflowratepersurface;

  // get condition
  std::vector< DRT::Condition * >      conds;
  dis.GetCondition ("SurfImpulsRate", conds);

  // collect elements by xfem coupling label
  for(std::vector<DRT::Condition*>::const_iterator conditer = conds.begin(); conditer!=conds.end(); ++conditer)
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
//        std::cout << (*impulsrates)[dis.DofColMap()->LID(gdofs[isd])] << std::endl;
      }
    }

//    LINALG::Matrix<3,1> flowrate(true);
//    dofrowmap->Comm().SumAll(&locflowrate(0),&flowrate(0),1);
//    dofrowmap->Comm().SumAll(&locflowrate(1),&flowrate(1),1);
//    dofrowmap->Comm().SumAll(&locflowrate(2),&flowrate(2),1);
//    std::cout << "locflowrate " << locflowrate << std::endl;
    if (volumeflowratepersurface.find(condID) == volumeflowratepersurface.end())
    {
      LINALG::Matrix<3,1> tmp(true);
      volumeflowratepersurface.insert(std::make_pair(condID,tmp));
    }
    LINALG::Matrix<3,1> tmp = volumeflowratepersurface[condID];
    tmp += locflowrate;
    volumeflowratepersurface[condID] = tmp;
  }

  return volumeflowratepersurface;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::WriteDoublesToFile(
  const double                     time,
  const int                        step,
  const std::map<int,double>&      data,
  const std::string                name
  )
{
  if (data.empty())
    dserror("data vector is empty");

  // print to file
  std::ostringstream header;
  header << std::right << std::setw(16) << "Time"
         << std::right << std::setw(10) << "Step"
         << std::right << std::setw(10) << "ID"
         << std::right << std::setw(16) << name;

  for(std::map<int,double >::const_iterator iter = data.begin(); iter != data.end(); ++iter)
  {
    std::ostringstream s;
    s << std::right << std::setw(16) << std::scientific << time
      << std::right << std::setw(10) << std::scientific << step
      << std::right << std::setw(10) << std::scientific << iter->first
      << std::right << std::setw(16) << std::scientific << iter->second;

    std::ostringstream slabel;
    slabel << std::setw(3) << std::setfill('0') << iter->first;
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + "." + name + "_ID_"+slabel.str()+".txt";

    if (step <= 1)
      f.open(fname.c_str(),std::fstream::trunc); //f << header.str() << std::endl;
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << s.str() << "\n";
    f.close();
  }
}
