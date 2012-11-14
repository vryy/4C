/*!------------------------------------------------------------------------------------------------*
\file cut_parallel.cpp

\brief provides the basic parallel cut classes "Parallel"

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_pack_buffer.H"

#include "cut_volumecell.H"

#include "cut_parallel.H"



/*------------------------------------------------------------------------------------------------*
 * basic CUT parallel constructor                                                    schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
GEO::CUT::Parallel::Parallel(
    DRT::Discretization & discret,
    GEO::CUT::Mesh & mesh,
    GEO::CUT::MeshIntersection & meshintersection
) :
discret_(discret),
myrank_(discret_.Comm().MyPID()),
numproc_(discret_.Comm().NumProc()),
mesh_(mesh),
meshintersection_(meshintersection)
{
  return;
} // end constructor


/*------------------------------------------------------------------------------------------------*
 * communicate node positions                                                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::CommunicateNodePositions()
{

  // wait for all processors before the Communication starts
  discret_.Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::CommunicateNodePositions" );

//  if(myrank_==0) std::cout << "\n\t ... CommunicateNodePositions" << std::flush;
//
//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------


  // check if there are undecided node positions on this proc
  mesh_.CheckForUndecidedNodePositions(curr_undecidedNodePos_);



  int counter = 0; // loop counter to avoid infinite loops

  while(true) // loop as long as all procs finished with decided node positions for all their nodes
  {
    counter += 1;

    // avoid infinite loops, after numproc rounds all procs should know their node positions
    if(counter > numproc_+10) dserror("number of rounds for exchanging data in CommunicateNodePositions > numproc_+10");

    // check if there are undecided node positions on this proc
    bool undecided_nodes = mesh_.CheckForUndecidedNodePositions(curr_undecidedNodePos_);

    // no undecided node positions -> procDone=true
    // undecided node positions -> procDone=false
    bool procDone = !undecided_nodes;

    // export if current proc is done or not
    // explanation: export is not finished when another proc (especially previous proc) is not done
    // * after one round all procs have status procDone=false if at least one proc has undecided node positions
    // * after one round all procs have status procDone=true  only if all procs have not undecided node positions
    exportCommunicationFinished(procDone);

    // check: procDone == 1 just if all procs have finished
    if (procDone)
      break;
    else
    {
      // perform a new Robin round to gather data from other procs
      // (send from current proc to next proc and receive info from proc before)
      // fill the current maps with information (node positions) from myproc
      exportNodePositionData();

      //-----------------------------------------------------------------------
      // ... now (back to the original proc) some of the ordered data could have been obtained from other procs
      // try do distribute the information (node positions)
      distributeMyReceivedNodePositionData();
    }
  } // end while loop

#ifdef DEBUG
  if(myrank_== 0)
  {
    std::cout << "number of round Robin loops to check finished procs:\t"    << counter   << endl;
    std::cout << "number of round Robin loops to exchange node positions:\t" << counter-1 << endl;
  }
#endif

  discret_.Comm().Barrier();

  //----------------------------------------------------------

//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0 )
//  {
//    std::cout << " ... Success (" << t_diff  <<  " secs)";
//  }


  return;
}



/*------------------------------------------------------------------------------------------------*
 * export data whether proc has finished to set node positions                       schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::exportCommunicationFinished(bool & procDone)
{
  // Initialization
  int dest = myrank_+1; // destination proc (the "next" one)
  if(myrank_ == (numproc_-1))
    dest = 0;

  int source = myrank_-1; // source proc (the "last" one)
  if(myrank_ == 0)
    source = numproc_-1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc=0;iproc<numproc_-1;iproc++)
  {
    DRT::PackBuffer dataSend;

    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));

    vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;
    int allProcsDone;

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,allProcsDone);

    // if the received information is allProcsDone==false, then set the current proc also to procDone=false
    // within the next round-iteration the next proc is also set to procDone=false
    if (allProcsDone==0)
      procDone = 0;

    // processors wait for each other
    discret_.Comm().Barrier();
  }


  return;
} // end exportFinishedData



/*------------------------------------------------------------------------------------------------*
 * export position data to neighbor proc and receive data from previous proc         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::exportNodePositionData()
{

  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;


  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    //---------------------------------------------------------------------------------------------------------------
    //--------------------- send data to next proc and receive from previous proc -----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // send current node position map to next proc and receive a new map from previous proc
    {
      DRT::PackBuffer dataSend; // data to be sent

      DRT::ParObject::AddtoPack(dataSend,curr_undecidedNodePos_);

      dataSend.StartPacking();

      DRT::ParObject::AddtoPack(dataSend,curr_undecidedNodePos_);

      vector<char> dataRecv;
      sendData(dataSend,dest,source,dataRecv);

      // pointer to current position of group of cells in global string (counts bytes)
      vector<char>::size_type posinData = 0;

      // clear vector that should be filled
      curr_undecidedNodePos_.clear();

      // unpack received data
      while (posinData < dataRecv.size())
      {
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, curr_undecidedNodePos_);
      }

      discret_.Comm().Barrier(); // processors wait for each other
    }


    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- fill maps with data from current proc  -------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // find node positions for received nodes and set position if possible
    for(std::map<int,int>::iterator it=curr_undecidedNodePos_.begin(); it!=curr_undecidedNodePos_.end(); it++)
    {
      int nid = it->first;
      int pos = it->second;

      // find the node on current proc
      Node* n = mesh_.GetNode(nid);

      if(n!=NULL) // node on this proc found
      {
        Point* p = n->point();

        Point::PointPosition my_pos = p->Position();


        if( (pos == Point::undecided) and (my_pos != Point::undecided))
        {
          // set the new position
          it->second = my_pos;
        }
        // both position are undecided
        else if( (pos == Point::undecided) and (my_pos == Point::undecided) )
        {
          // no gain in information
        }
        else if( (pos != Point::undecided) and (my_pos != Point::undecided)
                  and (pos != my_pos) )
        { // only an additional check for not unique node positions
          dserror("current position on myproc and received position for node %d are set, but not equal", nid);
        }
      }

    } // end iteration map

    //---------------------------------------------------------------------------------------------------------------

  } // end loop over procs


} // end exportNodePositionData




/*------------------------------------------------------------------------------------------------*
 * distribute received node positions on my processor                                schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::distributeMyReceivedNodePositionData()
{

  // distribute the received data for myproc
  for(std::map<int,int>::iterator it=curr_undecidedNodePos_.begin(); it != curr_undecidedNodePos_.end(); it++)
  {
    int nid = it->first;
    Point::PointPosition received_pos = (Point::PointPosition)it->second;

    // find the node on current proc (remark: this node has to be on this proc, it's the original source proc)
    Node* n = mesh_.GetNode(nid);

    if(n!=NULL)
    {
      Point* p = n->point();

      // set the new position for this point and distribute the information via facets and volumecells
      if(received_pos!=Point::undecided)
      {
        //cout << "reset the position for node " << nid << endl;
        p->Position(received_pos);
      }

      if( received_pos==Point::outside or received_pos==Point::inside )
      {
        // The nodal position is already known. Set it to my facets. If the
        // facets are already set, this will not have much effect anyway. But on
        // multiple cuts we avoid unset facets this way.
        const plain_facet_set & facets = p->Facets();
        for ( plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); ++i )
        {
          Facet * f = *i;
          f->Position( received_pos );
        }
      }

    }
    else dserror(" this node should be available on its original proc that asked other procs for nodepositions");

  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 * communicate the node dofset number for single volumecells                         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::CommunicateNodeDofSetNumbers()
{

  // wait for all processors before the Communication starts
  discret_.Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::CommunicateNodeDofSetNumbers" );

//  if(myrank_==0) std::cout << "\n\t ... CommunicateNodeDofSetNumbers" << std::flush;
//
//  const double t_start = Teuchos::Time::wallTime();



  dofSetData_ = Teuchos::rcp(new std::vector<MeshIntersection::DofSetData>);

  // check if there are missing dofset for volumecells on this proc
  meshintersection_.FillParallelDofSetData(dofSetData_, discret_);


  // perform just one Robin round to gather data from other procs
  // (send from current proc to next proc and receive info from proc before)
  // fill the current maps with information (dofset number for vc and the current row node) from myproc
  exportDofSetData();

  //-----------------------------------------------------------------------
  // ... now (back to the original proc) all the ordered data should have been obtained from other procs
  // set the information
  distributeDofSetData();

  discret_.Comm().Barrier();

//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0 )
//  {
//    std::cout << " ... Success (" << t_diff  <<  " secs)\n";
//  }


  return;
}


/*------------------------------------------------------------------------------------------------*
 * export dofset data to neighbor proc and receive data from previous proc           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::exportDofSetData()
{

  // destination proc (the "next" one)
  int dest = myrank_+1;
  if(myrank_ == (numproc_-1))
    dest = 0;

  // source proc (the "last" one)
  int source = myrank_-1;
  if(myrank_ == 0)
    source = numproc_-1;


  // loop over processors
  for (int procid=0; procid<numproc_; procid++)
  {
    //---------------------------------------------------------------------------------------------------------------
    //--------------------- send data to next proc and receive from previous proc -----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // send current DofSetData to next proc and receive a new map from previous proc
    {
      DRT::PackBuffer dataSend; // data to be sent

      // packing the data
      for (vector<MeshIntersection::DofSetData>::iterator data=dofSetData_->begin(); data!=dofSetData_->end(); data++)
      {
        DRT::ParObject::AddtoPack(dataSend,data->set_index_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->inside_cell_);
        packPoints(dataSend, data->cut_points_coords_);
        DRT::ParObject::AddtoPack(dataSend,data->peid_);
        DRT::ParObject::AddtoPack(dataSend,data->node_dofsetnumber_map_);
      }

      dataSend.StartPacking();

      // packing the data
      for (vector<MeshIntersection::DofSetData>::iterator data=dofSetData_->begin(); data!=dofSetData_->end(); data++)
      {
        DRT::ParObject::AddtoPack(dataSend,data->set_index_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->inside_cell_);
        packPoints(dataSend, data->cut_points_coords_);
        DRT::ParObject::AddtoPack(dataSend,data->peid_);
        DRT::ParObject::AddtoPack(dataSend,data->node_dofsetnumber_map_);
      }

      vector<char> dataRecv;
      sendData(dataSend,dest,source,dataRecv);

      // pointer to current position of group of cells in global string (counts bytes)
      vector<char>::size_type posinData = 0;

      // clear vector that should be filled
      dofSetData_->clear();

      // unpack received data
      while (posinData < dataRecv.size())
      {
        // unpack volumecell
        int set_index = -1;                                  // set index for Volumecell
        int inside_cell = false;                             // inside or outside cell
        std::vector<LINALG::Matrix<3,1> > cut_points_coords; // coordinates of cut points
        int peid = -1;                                       // parent element id for volume cell
        std::map<int,int> node_dofsetnumber_map;             // map <nid, current dofset number>

        // unpack volumecell data
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, set_index);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, inside_cell);
        unpackPoints(posinData,dataRecv,cut_points_coords);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, peid);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, node_dofsetnumber_map);

        // create a new dofSetData object with unpacked data
        dofSetData_->push_back( GEO::CUT::MeshIntersection::DofSetData( set_index, (bool)inside_cell, cut_points_coords, peid, node_dofsetnumber_map) );
      }



      discret_.Comm().Barrier(); // processors wait for each other
    }


    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- fill maps with data from current proc  -------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(std::vector<MeshIntersection::DofSetData>::iterator vc_data=dofSetData_->begin();
        vc_data != dofSetData_->end();
        vc_data++)
    {
      bool find_volumecell = false; // do we have to identify the received volumecell on myrank?

      std::map<int,int>& node_dofsetnumber_map = vc_data->node_dofsetnumber_map_;


      // check if it necessary to find the received volumecell on myrank
      // required, if at least one node in the received node_dofsetnumber_map is a row node on myrank
      for(std::map<int,int>::iterator node_dofsetnumber_it = node_dofsetnumber_map.begin();
          node_dofsetnumber_it != node_dofsetnumber_map.end();
          node_dofsetnumber_it++)
      {
        int nid = node_dofsetnumber_it->first;

        bool haveGlobalNode = discret_.HaveGlobalNode(nid);

        if(haveGlobalNode)
        {
          DRT::Node* node = discret_.gNode(nid);

          if(node->Owner() == myrank_)
          {
            find_volumecell = true;
            //cout << "find the vc!!!" << endl;
            break; // at least one node found as row node, we have to identify the received volumecell on myrank
          }
        }
      }

      //----------------------------------------------------------------------
      if(!find_volumecell) continue; // check the next received dofSetData_
      //----------------------------------------------------------------------

      //parent element Id for current vc data
      int peid = vc_data->peid_;

      // find the volumecell on myrank
      VolumeCell * my_vc = NULL;
      my_vc = findVolumeCell( *vc_data);

      if(my_vc == NULL) dserror("no corresponding volumecell for vc in element %d found", peid);


      for(std::map<int,int>::iterator node_dofsetnumber_it = node_dofsetnumber_map.begin();
          node_dofsetnumber_it != node_dofsetnumber_map.end();
          node_dofsetnumber_it++)
      {
        //cout << "loop the node_dofsetnumber_map" << endl;
        int nid = node_dofsetnumber_it->first;
        int curr_dofset_number = node_dofsetnumber_it->second;

        bool haveGlobalNode = discret_.HaveGlobalNode(nid);

        // decide if the current proc carries the required information
        if(haveGlobalNode) // node on this proc available as row or col node
        {
          //cout << "in haveGlobalNode for node " << nid << endl;
          DRT::Node* node = discret_.gNode(nid);
          if(node->Owner() == myrank_)
          {

            int new_dofset_number = -1;

            // find the local index of the current node w.r.t the element
            int index = -1;

            index = getDofSetVecIndex(nid, peid);

            //cout << "index for node " << nid << " is: " << index << endl;

            // get the right dofsetnumber for the node
            if(my_vc != NULL)
            {
              //cout << "in my_vc != NULL" << endl;
              const std::vector<int> nds = my_vc->NodalDofSet();

              if(index >= (int)nds.size()) dserror(" index can not be read in nds vector of size %d ", index, nds.size());

              new_dofset_number = nds[index];

              if(new_dofset_number==-1) dserror("the new dofset number for node %d is not valid", nid);
            }
            else dserror("Volumecell Pointer is NULL");



            // set the new dofset-number
            if(curr_dofset_number == -1)
            {
              node_dofsetnumber_it->second = new_dofset_number;
            }
            else dserror("dofset for this node for this volumecell already set by another proc, this should be done just by the row-node proc");


          } // end if myrank_
        } // have global node on this proc
      } // end loop node_dofsetnumber_it ( some elemental nodes for the current received volumecell)
    }// end loop dofSetData_ (data for volumecells)


//    cout << "replaced data: " << endl;
//    printDofSetData();


    //---------------------------------------------------------------------------------------------------------------

  } // end loop over procs

  return;
} // end exportNodePositionData




/*------------------------------------------------------------------------------------------------*
 * distribute received dofset number for nodes of volumecell on my processor         schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::distributeDofSetData()
{

  // back to the original proc

  for(std::vector<MeshIntersection::DofSetData>::iterator data = dofSetData_->begin(); data!=dofSetData_->end(); data++)
  {
    // set data in first volumecell of set with setindex

    // get the Volumecell
    // find volumecell sets and non-row nodes for that dofset numbers has to be communicated parallel
    // the communication is done element wise for all its sets of volumecells when there is a non-row node in this element

    int peid = data->peid_ ;

    GEO::CUT::ElementHandle * e = meshintersection_.GetElement( peid );

    std::vector<int> nds;

    if(e!=NULL)
    {

      if(data->inside_cell_)
      {
        const std::vector<plain_volumecell_set> & ele_vc_sets_inside = e->GetVcSetsInside();
        std::vector<std::vector<int> > & nodaldofset_vc_sets_inside  = e->GetNodalDofSet_VcSets_Inside();

        ReplaceNdsVectors (e, ele_vc_sets_inside, nodaldofset_vc_sets_inside, data->set_index_, data->node_dofsetnumber_map_);
      }
      else
      {
        const std::vector<plain_volumecell_set> & ele_vc_sets_outside = e->GetVcSetsOutside();
        std::vector<std::vector<int> > & nodaldofset_vc_sets_outside  = e->GetNodalDofSet_VcSets_Outside();

        ReplaceNdsVectors (e, ele_vc_sets_outside, nodaldofset_vc_sets_outside,  data->set_index_, data->node_dofsetnumber_map_);
      }
    }
    else
    {
      dserror("there must be an elementhandle for element %d on my original proc %d", peid, myrank_);
    }


    //... and set data in elementhandle


    // check if all nodes are filled with a valid dofset number
    for(std::map<int,int>::iterator dofnumber_map=data->node_dofsetnumber_map_.begin();
        dofnumber_map!=data->node_dofsetnumber_map_.end();
        dofnumber_map++)
    {
      // safety check if setting dofset for data was successful
      if(dofnumber_map->second == -1)
      {
        dserror( "dofset number for node %d for vc in element %d could not be determined", dofnumber_map->first, data->peid_ );
      }


    }
  }


  return;
}

/*------------------------------------------------------------------------------------------------*
 * find the volumecell on myrank for which we received data stored in vc_data        schott 10/12 *
 *------------------------------------------------------------------------------------------------*/
GEO::CUT::VolumeCell* GEO::CUT::Parallel::findVolumeCell(
    MeshIntersection::DofSetData& vc_data           ///< volumecell data which have to be identified on myrank
    )
{
  bool vc_found = true;

  VolumeCell* my_vc = NULL;

  ElementHandle* pele = meshintersection_.GetElement(vc_data.peid_);

  if(pele == NULL) dserror("element with Id %i not found on proc %i", vc_data.peid_, myrank_);

  plain_volumecell_set my_vcs;

  pele->VolumeCells(my_vcs);

  // find all the received volumecell's points
  for(plain_volumecell_set::iterator c = my_vcs.begin(); c!=my_vcs.end(); c++)
  {
    VolumeCell* cell = *c;

    // assume that the corresponding vc on myrank_ has been found
    vc_found = true;

    // get points for my_vc
    std::vector<GEO::CUT::Point* > my_cut_points;
    const plain_facet_set facets = cell->Facets();

    for(plain_facet_set::const_iterator i=facets.begin(); i!=facets.end(); i++)
    {
      Facet* f = *i;
      std::vector<Point*> facetpoints = f->Points();
      std::copy(facetpoints.begin(), facetpoints.end(), std::inserter(my_cut_points, my_cut_points.begin()));
    }


    //-----------------------------------------------------------------------------
    // find received cut_points in my_cut_points

    // compare the size of vecs
    if(my_cut_points.size() != (vc_data.cut_points_coords_).size())
    {

      //cout << "my_cut_points.size() != (vc_data->cut_points_coords_).size()" << endl;
      // no identical number of points found!
      // this volumecell is not a candidate
      vc_found = false;
    }
    else
    {
      //cout << "my_cut_points.size() == (vc_data->cut_points_coords_).size()" << endl;
      // identical number of points found! this vc is a candidate: check the points"

      // brute force search for identical point coords
      for(std::vector<GEO::CUT::Point*>::iterator my_it=my_cut_points.begin(); my_it!=my_cut_points.end(); my_it++)
      {
        bool point_found = true;

        for(std::vector<LINALG::Matrix<3,1> >::iterator rec_it=(vc_data.cut_points_coords_).begin();
            rec_it != (vc_data.cut_points_coords_).end();
            rec_it++)
        {

          point_found = true;

          //compare the coordinates
          for(int i=0; i<3; i++)
          {
            if(fabs((*my_it)->X()[i] - (*rec_it)(i)) > PARALLEL_COORD_TOL)
            {
              point_found = false;
              break; // stop the loop over coordinates
            }
          }

          if(point_found==true) break; // stop the inner! loop over received vc's points if the point is already found

        } // inner loop over received vc's points

        if(point_found == false)
        {
          vc_found=false;
          break; // do not search for the next points, it is the wrong volumecell
        }
      } // end loop over my_points
    } // end else


    if(vc_found == true)
    {
      my_vc = cell;
      //cout << "volumecell has been found!!! Happy :-)" << endl;
      break; // stop the loop over volumecells
    }
    else
    {
      //my_vc = NULL;
      //cout << "volumecell has not been found!!! :-(" << endl;
    }


    //-----------------------------------------------------------------------------
  }// it over volumecells


  if(vc_found == false) dserror("no corresponding volumecell for vc in element %d found for node %d", vc_data.peid_);

  return my_vc;
}



void GEO::CUT::Parallel::ReplaceNdsVectors (ElementHandle*                            e,
                                            const std::vector<plain_volumecell_set> & ele_vc_sets,
                                            std::vector<std::vector<int> > &          nodaldofset_vc_sets,
                                            int                                       set_index,
                                            std::map<int,int>&                        node_dofsetnumber_map)
{

  // get the original set
  plain_volumecell_set cells = ele_vc_sets[set_index];

  // get the old nds vector for all cells in current set, represented by the first Volumecell
  const std::vector<int> nds_old = cells[0]->NodalDofSet();

  // get vector of nids in current element
  const std::vector<GEO::CUT::Node* > nodes = e->Nodes();

  // create the new dofset Vector
  std::vector<int> nds_new;

  // fill the new nds-vector
  for(int i=0; i< (int)nds_old.size(); i++)
  {

    if(nds_old[i] == -1) // entry that has to be replaced with received data
    {
      int nid = nodes[i]->Id(); // key for received node_dofsetnumber_map
      nds_new.push_back( node_dofsetnumber_map.find(nid)->second );
    }
    else // set original dofset number
    {
      nds_new.push_back(nds_old[i]);
    }
  }

  // set the new nds vector for all volumecells in the current set
  for(plain_volumecell_set::iterator c=cells.begin(); c!=cells.end(); c++)
  {
    (*c)->SetNodalDofSet(nds_new);
  }

  nodaldofset_vc_sets[set_index] = nds_new;

  return;
}


/*------------------------------------------------------------------------------------------------*
 * packing a point for parallel communication only with the basic point data                      *
 * without an underlying discretization fitting to the node's new prozessor          schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::packPoints(
    DRT::PackBuffer&              dataSend,
    vector<LINALG::Matrix<3,1> >& points_coords
) const
{
  const int nsd = 3;

  // pack number of points for current volumecell
  DRT::ParObject::AddtoPack(dataSend, (int)points_coords.size());

  for(std::vector<LINALG::Matrix<3,1> >::iterator p=points_coords.begin(); p!=points_coords.end(); p++)
  {
    // pack xyz-coordinates
    DRT::ParObject::AddtoPack(dataSend,LINALG::Matrix<nsd,1>(*p));
  }

} // end packNodes



/*------------------------------------------------------------------------------------------------*
 * unpacking a point for parallel communication only with the basic point data                    *
 * without an underlying discretization fitting to the node's new prozessor          schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::unpackPoints(
    vector<char>::size_type&       posinData,
    vector<char>&                  dataRecv,
    vector<LINALG::Matrix<3,1> >&  points_coords
) const
{
  const int nsd = 3; // dimension
//  int id; // global id
//  LINALG::Matrix<nsd,1> coords; // coordinates
//  int owner; // processor

  int num_points = 0;

  // unpack number of points for current volumecell
  DRT::ParObject::ExtractfromPack(posinData,dataRecv,num_points);


  LINALG::Matrix<nsd,1> coords(true);

  for(int i=0; i< num_points; i++)
  {
    coords.Clear();

    // pack xyz-coordinates for point
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,coords);

    points_coords.push_back(coords);
  }

} // end function unpackNodes



/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source                schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::sendData(
    DRT::PackBuffer&      dataSend,
    int&                  dest,
    int&                  source,
    std::vector<char>&    dataRecv
) const
{

  std::vector<int> lengthSend(1,0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef DEBUG
  cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << endl;
#endif

  // exporter for sending
  DRT::Exporter exporter(discret_.Comm());

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
  // ... and receive length
  std::vector<int> lengthRecv(1,0);
  exporter.Receive(source, length_tag, lengthRecv, size_one);
  exporter.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);

  // ... and receive data
  dataRecv.clear(); dataRecv.resize(lengthRecv[0]);
  exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
  exporter.Wait(req_data);

#ifdef DEBUG
  cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << endl;
#endif
} // end sendData



/*------------------------------------------------------------------------------------------------*
 * print current stored dofSetData_                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::Parallel::printDofSetData()
{
  for(std::vector<MeshIntersection::DofSetData>::iterator i=dofSetData_->begin(); i!=dofSetData_->end(); i++)
  {
    i->print();
  }
}

/*------------------------------------------------------------------------------------------------*
 * get the index of nid in the vector of elements (eid) node Ids                     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
int GEO::CUT::Parallel::getDofSetVecIndex(int nid, int eid)
{
  DRT::Element* ele = discret_.gElement(eid);

  if(ele == NULL) dserror("element %d not available on proc %d", eid, myrank_);

  int numnode = ele->NumNode();

  const int* ele_n_ids = ele->NodeIds();

  int index = 0;

  // find the local node Index for received nid w.r.t the element
  for(; index< numnode; ++index)
  {
    if(ele_n_ids[index] == nid) return index;
  }

  if(index > numnode)
  {
    dserror("nid %d not found in nodes for element %d!", nid, eid);
    return -1;
  }

  return -1;
}




