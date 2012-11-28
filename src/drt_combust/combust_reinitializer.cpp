/*!-----------------------------------------------------------------------------------------------*
 \file combust_reinitializer.cpp

  \brief reinitialization of G-function (level set function)

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "combust_reinitializer.H"
#include "combust_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_scatra/scatra_timint_implicit.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Reinitializer::Reinitializer(
    const Teuchos::ParameterList&              combustdyn,
    SCATRA::ScaTraTimIntImpl&                  scatra,
    const std::map<int,GEO::BoundaryIntCells>& boundaryintcells,
    Teuchos::RCP<Epetra_Vector>                phivector,
    const bool                                 compdist // flag for initial computation of signed distance function
  ) :
    combustdyn_(combustdyn),
    reinitaction_(DRT::INPUT::IntegralValue<INPAR::COMBUST::ReInitialActionGfunc>(combustdyn_.sublist("COMBUSTION GFUNCTION"),"REINITIALIZATION")),
    reinitband_(DRT::INPUT::IntegralValue<int>(combustdyn.sublist("COMBUSTION GFUNCTION"),"REINITBAND")),
    reinitbandwidth_(combustdyn.sublist("COMBUSTION GFUNCTION").get<double>("REINITBANDWIDTH")),
    scatra_(scatra),
    flamefront_(boundaryintcells)
{
  if (compdist)
  {
    switch (reinitaction_)
    {
    case INPAR::COMBUST::reinitaction_none:
      // do nothing
      break;
    case INPAR::COMBUST::reinitaction_byfunction:
      // read a FUNCTION from the input file and reinitialize the G-function by evaluating it
      scatra_.SetInitialField(INPAR::SCATRA::initfield_field_by_function,combustdyn_.sublist("COMBUSTION GFUNCTION").get<int>("REINITFUNCNO"));
      break;
    case INPAR::COMBUST::reinitaction_signeddistancefunction:
      SignedDistanceFunction(phivector);
      break;
    case INPAR::COMBUST::reinitaction_fastsigneddistancefunction:
      FastSignedDistanceFunction(phivector);
      break;
    default:
      dserror ("unknown option to reinitialize the G-function");
  }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Reinitializer::~Reinitializer()
{
  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: build signed distance function                                            henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Reinitializer::SignedDistanceFunction(Teuchos::RCP<Epetra_Vector> phivector)
{
  // get communicator (for output)
  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
  if (comm.MyPID()==0)
  {
    std::cout << "---  reinitializing G-function by computing distance to interface ..." << std::flush;
#ifdef COMBUST_2D
    std::cout << "\n /!\\ third component or normal vector set to 0 to keep 2D-character!" << std::endl;
#endif
  }
  // get a pointer to the G-function discretization
  Teuchos::RCP<DRT::Discretization> gfuncdis = scatra_.Discretization();

  // get map holding pairs of nodes (master-slave) on periodic boundaries (pbc nodes)
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap = scatra_.PBCmap();

  // get the pbc itself
  Teuchos::RCP<PeriodicBoundaryConditions> pbc = scatra_.PBC();

  // map holding pbc nodes (masters or slaves) <pbc node id, distance to flame front>
  std::map<int, double> pbcnodes;

  const Epetra_Map* dofrowmap = gfuncdis->DofRowMap();

  // get the following information about the pbc
  // - planenormaldirection e.g. (1,0,0)
  // - minimum in planenormaldirection
  // - maximum in planenormaldirection
  std::vector<DRT::Condition*>* surfacepbcs = pbc->ReturnSurfacePBCs();
  std::vector<int>    planenormal(0);
  std::vector<double> globalmins (0);
  std::vector<double> globalmaxs (0);
  for (size_t i = 0; i < surfacepbcs->size(); ++i)
  {
    const string* ismaster = (*surfacepbcs)[i]->Get<string>("Is slave periodic boundary condition");
    if (*ismaster == "Master")
    {
      const int masterid = (*surfacepbcs)[i]->GetInt("Id of periodic boundary condition");
      std::vector<int> nodeids(*((*surfacepbcs)[i]->Nodes()));
      for (size_t j = 0; j < surfacepbcs->size(); ++j)
      {
        const int slaveid = (*surfacepbcs)[j]->GetInt("Id of periodic boundary condition");
        if (masterid == slaveid)
        {
          const string* isslave = (*surfacepbcs)[j]->Get<string>("Is slave periodic boundary condition");
          if (*isslave == "Slave")
          {
            const std::vector<int>* slavenodeids = (*surfacepbcs)[j]->Nodes();
            // append slave node Ids to node Ids for the complete condition
            for (size_t k = 0; k < slavenodeids->size(); ++k)
              nodeids.push_back(slavenodeids->at(k));
          }
        }
      }

      // Get normal direction of pbc plane
      const string* pbcplane = (*surfacepbcs)[i]->Get<string>("degrees of freedom for the pbc plane");
      if (*pbcplane == "yz")
        planenormal.push_back(0);
      else if (*pbcplane == "xz")
        planenormal.push_back(1);
      else if (*pbcplane == "xy")
        planenormal.push_back(2);
      else
        dserror("A PBC condition could not provide a plane normal.");

      double min = +10e19;
      double max = -10e19;
      for (size_t j = 0; j < nodeids.size(); ++j)
      {

        const int gid = nodeids[j];
        const int lid = gfuncdis->NodeRowMap()->LID(gid);
        if (lid < 0)
          continue;
        const DRT::Node* lnode = gfuncdis->lRowNode(lid);
        const double* coord = lnode->X();
        if (coord[planenormal.back()] < min)
          min = coord[planenormal.back()];
        if (coord[planenormal.back()] > max)
          max = coord[planenormal.back()];
      }
      globalmins.resize(planenormal.size());
      globalmaxs.resize(planenormal.size());
      gfuncdis->Comm().MinAll(&min, &(globalmins.back()), 1);
      gfuncdis->Comm().MaxAll(&max, &(globalmaxs.back()), 1);
    }
  }//end loop over all surfacepbcs

  //  for (size_t i = 0; i < planenormal.size(); ++i)
  //  {
  //    cout << "Normal: " << planenormal[i] << " Min: " << globalmins[i] << " Max: " << globalmaxs[i] << endl;
  //  }


  // vector of coordinates of node, the G-function value of which is to be reinitialized
  static LINALG::Matrix<3,1> nodecoord(true);
  // normal vector of a flame front patch
  static LINALG::Matrix<3,1> normal(true);

  //------------------------------------
  // loop all row nodes on the processor
  //------------------------------------
  for(int lnodeid=0; lnodeid < gfuncdis->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    const DRT::Node* lnode = gfuncdis->lRowNode(lnodeid);
    //cout << "proc " << comm.MyPID() << " reinitialization for node: " << lnode->Id() << endl;

    // get the dof associated with this node
    const std::vector<int> nodedofset = gfuncdis->Dof(lnode); // this should not be a vector, but a scalar!
#ifdef DEBUG
    int numdofs = nodedofset.size(); // this should be 1 (scalar field!)
    if (numdofs != 1)
      dserror("There are more than 1 dof for a node in a scalar field!");
#endif
    const int dofgid = nodedofset[0];
    int doflid = dofrowmap->LID(dofgid);

    //---------------------------------------------------------------------------------
    // check if this node is a node with periodic boundary conditions (master or slave)
    //---------------------------------------------------------------------------------
    const int nodeid = lnode->Id();
    // boolean indicating whether this node is a pbc node
    bool pbcnode = false;
    for (std::map<int, std::vector<int>  >::const_iterator iter= (*pbcmap).begin(); iter != (*pbcmap).end(); ++iter)
    {
      if (iter->first == nodeid) // node is a pbc master node
        pbcnode = true;
      else
      {
        for (size_t islave=0;islave<iter->second.size();islave++)
        {
          if (iter->second[islave] == nodeid) // node is a pbc slave node
            pbcnode = true;
        }
      }
    }

    // get physical coordinates of this node
    nodecoord(0) = lnode->X()[0];
    nodecoord(1) = lnode->X()[1];
    nodecoord(2) = lnode->X()[2];

    //--------------------------------
    // conditions for reinitialization
    //--------------------------------
    if (!reinitband_ or // reinitialize entire level set field
        (reinitband_ and fabs((*phivector)[doflid]) <= reinitbandwidth_)) // reinitialize only within a band around the zero level set
    {
      // smallest distance to the vertex of a flame front patch
      double vertexdist = 7777.7; // default value
      // smallest distance to the edge of a flame front patch
      double edgedist = 6666.6; // default value
      // smallest distance to flame front
      double mindist = 5555.5; // default value

      //--------------------------------
      // due to the PBCs the node might actually be closer to the
      // interface then would be calculated if one only considered
      // the actual position of the node. In order to find the
      // smallest distance the node is copied along all PBC directions
      //
      //   +------------------+ - - - - - - - - - -+
      //   +             II   +
      //   +   x        I  I  +    y               +
      //   +             II   +
      //   +------------------+ - - - - - - - - - -+
      //         original           copy
      //
      //   x: current node
      //   y: copy of current node
      //   I: interface
      //   +: pbc

      if (planenormal.size() > 3)
        dserror("Sorry, but currently a maximum of three periodic boundary conditions are supported by the combustion reinitializer.");

      // since there is stl pow(INT, INT) function, we calculate it manually
      size_t looplimit = 1;
      for (size_t i = 0; i < planenormal.size(); ++i)
        looplimit *= 2;

      for (size_t ipbc = 0; ipbc < looplimit; ++ipbc)
      {
        LINALG::Matrix<3,1> tmpcoord(nodecoord);

        // determine which pbcs have to be applied
        //
        // loopcounter | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
        // ------------+---+---+---+---+---+---+---+---+
        //  first PBC  |     x       x       x       x
        // second PBC  |         x   x           x   x
        //  third PBC  |                 x   x   x   x
        //
        // this is equivalent to the binary representation
        // of the size_t
        if (ipbc & 0x01)
        {
          const double pbclength = globalmaxs[0] - globalmins[0];
          if (nodecoord(planenormal[0]) > globalmins[0] + pbclength/2.0)
            tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) - pbclength;
          else
            tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) + pbclength;
        }
        if (ipbc & 0x02)
        {
          const double pbclength = globalmaxs[1] - globalmins[1];
          if (nodecoord(planenormal[1]) > globalmins[1] + pbclength/2.0)
            tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) - pbclength;
          else
            tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) + pbclength;
        }
        if (ipbc & 0x04)
        {
          const double pbclength = globalmaxs[2] - globalmins[2];
          if (nodecoord(planenormal[2]) > globalmins[2] + pbclength/2.0)
            tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) - pbclength;
          else
            tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) + pbclength;
        }
        //-----------------------------------------------------------
        // compute smallest distance to the flame front for this node
        //-----------------------------------------------------------
        if (flamefront_.empty())
          dserror("no flamefront patches available");
        // loop groups (vectors) of flamefront patches of all elements
        for(std::map<int,GEO::BoundaryIntCells>::const_iterator elepatches = flamefront_.begin(); elepatches != flamefront_.end(); ++elepatches)
      {
        // number of flamefront patches for this element
        const std::vector<GEO::BoundaryIntCell> patches = elepatches->second;
        const int numpatch = patches.size();

        //-----------------------------------------
        // loop flame front patches of this element
        //-----------------------------------------
        for(int ipatch=0; ipatch<numpatch; ++ipatch)
        {
          // get a single patch from group of flamefront patches
          const GEO::BoundaryIntCell patch = patches[ipatch];

          // only triangles and quadrangles are allowed as flame front patches (boundary cells)
          if (!(patch.Shape() == DRT::Element::tri3 or
                patch.Shape() == DRT::Element::quad4))
          {
            dserror("invalid type of boundary integration cell for reinitialization");
          }

          // get coordinates of vertices defining flame front patch
          const LINALG::SerialDenseMatrix& patchcoord = patch.CellNodalPosXYZ();

          // compute normal vector to flame front patch
          ComputeNormalVectorToFlameFront(patch,patchcoord,normal);

          //-----------------------------------------
          // find flame front patches facing the node
          //-----------------------------------------
          // boolean indicating if facing patch was found
          bool facenode = false;
          // distance to the facing patch
          double patchdist = 7777.7; // default value
          // check if this patch faces the node
          FindFacingPatchProjCellSpace(tmpcoord,patch,patchcoord,normal,facenode,patchdist);

          // a facing patch was found
          if (facenode == true)
          {
            // overwrite smallest distance if computed patch distance is smaller
            if (fabs(patchdist) < fabs(mindist))
            {
              // if G-value at the node is negative, the minimal distance has to be negative
              if ((*phivector)[doflid] < 0.0 )
                mindist = -patchdist;
              else
                mindist = patchdist;

              if (pbcnode)
              {
                // add node to map of pbc nodes
                // remark: this node is a pbc node and it has been projected on a facing patch
                //pbcpatchnodes[nodeid] = mindist;
                pbcnodes[nodeid] = mindist;
              }

            }
          }

          //-------------------------------------------------------------
          // compute smallest distance to edges of this flame front patch
          //-------------------------------------------------------------
          ComputeDistanceToEdge(tmpcoord,patch,patchcoord,edgedist);

          //----------------------------------------------------------------
          // compute smallest distance to vertices of this flame front patch
          //----------------------------------------------------------------
          ComputeDistanceToPatch(tmpcoord,patch,patchcoord,vertexdist);
        }
      }
      }// end loop over all pbc copied nodes


      //-------------------------------------------
      // determine smallest distance to flame front
      //-------------------------------------------
      // remark: variables have the following meaning here:
      //         - "mindist" is either "patchdist" or still the default value (5555.5)
      //         - "vertexdist" is the distance to the closest vertex of any flame front patch
      //
      // case 1: a local flame front patch was found for this node -> mindist = patchdist
      // case 2: a flame front patch was found for this node, but it is not local (curved interface);
      //         a local patch was not found, because this node is located in the blind angle of all
      //         local patches -> mindist = vertexdist or edgedist
      // case 3: something went wrong:  "mindist" still has default value (5555.5)

      //cout << "minimal distance to facing patch is " << mindist << endl;
      //cout << "mindist " << std::setw(18) << std::setprecision(12) << std::scientific << mindist << endl;
      //cout << "vertexdist " << std::setw(18) << std::setprecision(12) << std::scientific << vertexdist << endl;
      if (fabs(edgedist) < fabs(mindist)) // case 2a
      {
        // if G-value at the node is negative, the minimal distance has to be negative
        if ((*phivector)[doflid] < 0.0 )
          mindist = -edgedist;
        else
          mindist = edgedist;

        if (pbcnode)
        {
          // add node to map of pbc nodes
          // remark: this node is a pbc node and the closest distance to an edge has been computed
          pbcnodes[nodeid] = mindist;
        }
      }

      if (fabs(vertexdist) < fabs(mindist)) // case 2b
      {
        // if G-value at the node is negative, the minimal distance has to be negative
        if ((*phivector)[doflid] < 0.0 )
          mindist = -vertexdist;
        else
          mindist = vertexdist;
        //cout << "node " << lnode->Id() << " does not face a patch (blind angle) or this patch is not local; distance: " << mindist << endl;

        if (pbcnode)
        {
          // add node to map of pbc nodes
          // remark: this node is a pbc node and the closest distance to a vertex has been computed
          pbcnodes[nodeid] = mindist;
        }

      }
      if (mindist == 5555.5) // case 3
        dserror ("G-fuction value of node %d was not reinitialized!", lnode->Id());

      //------------------------------
      // reinitialize G-function field
      //------------------------------
      // assign new values of signed distance function to the G-function field
      int err = phivector->ReplaceMyValues(1,&mindist,&doflid);
      if (err) dserror("this did not work");
    } // end condition for band around zero level set
  } // end loop nodes

  //-----------------------------------------------------
  // reinitialize G-function field on periodic boundaries
  //-----------------------------------------------------
  // remark: - If a node could not be projected on a facing patch, but its corresponding pbc node could,
  //           then take the projected distance to reinitialize both nodes. This happens if the
  //           flame front hits the boundary not in a 90 degree angle
  //         - previously reinitialized G-function values on periodic boundaries are overwritten

  // loop all master nodes with periodic boundary conditions
  std::map<int, std::vector<int>  >::const_iterator pbciter;
  for (pbciter = (*pbcmap).begin(); pbciter != (*pbcmap).end(); ++pbciter)
  {
    // get master node gid
    const int mastergid = pbciter->first;
    // get slave node gids
    std::vector<int> slavegid;
    for (size_t islave = 0; islave < pbciter->second.size(); islave++)
      slavegid.push_back(pbciter->second[islave]);

    // get the dof associated with these nodes (master and slave share this dof)
    DRT::Node* masternode = gfuncdis->gNode(mastergid);
    const std::vector<int> nodedofset = gfuncdis->Dof(masternode); // this should not be a vector, but a scalar!
#ifdef DEBUG
    int numdofs = nodedofset.size(); // this should be 1 (scalar field!)
    if (numdofs != 1)
      dserror("There are more than 1 dof for a node in a scalar field!");
#endif
    const int dofgid = nodedofset[0];
    int doflid = dofrowmap->LID(dofgid);

    // look for this pbc master node in map of pbc patch nodes
    std::map<int, double>::const_iterator pbcmasteriter = pbcnodes.find(mastergid);

    double mindist = pbcmasteriter->second;

    for (size_t islave = 0; islave < slavegid.size(); islave++)
    {
      // look for this pbc slave node in map of pbc patch nodes
      std::map<int, double>::const_iterator pbcslaveiter = pbcnodes.find(slavegid[islave]);

      double slavemindist = pbcslaveiter->second;
      mindist = std::min(mindist,slavemindist);
    }

    int err = phivector->ReplaceMyValues(1,&mindist,&doflid);
    if (err) dserror("this did not work");
  }


  if (comm.MyPID()==0)
    std::cout << " done" << std::endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: build signed distance function fast                                   rasthofer 08/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Reinitializer::FastSignedDistanceFunction(Teuchos::RCP<Epetra_Vector> phivector)
{
  // get communicator (for output)
  const Epetra_Comm& comm = scatra_.Discretization()->Comm();

  if (comm.MyPID()==0)
  {
    std::cout << "---  reinitializing G-function by computing distance to interface really fast..." << std::flush;
#ifdef COMBUST_2D
    std::cout << "\n /!\\ third component or normal vector set to 0 to keep 2D-character!" << std::endl;
#endif
  }

  // get a pointer to the G-function discretization
  Teuchos::RCP<DRT::Discretization> gfuncdis = scatra_.Discretization();

  // get map holding pairs of nodes (master-slave) on periodic boundaries (pbc nodes)
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmap = scatra_.PBCmap();

  // get the pbc itself
  Teuchos::RCP<PeriodicBoundaryConditions> pbc = scatra_.PBC();

  // map holding pbc nodes (masters or slaves) <pbc node id, distance to flame front>
  std::map<int, double> pbcnodes;

  const Epetra_Map* dofrowmap = gfuncdis->DofRowMap();

  // determine the number of nodes per element
  int numnodesperele = 0;
  if (gfuncdis->NumMyRowElements() <= 0)
    dserror("This discretization does not have any row elements.");
  switch (gfuncdis->lRowElement(0)->Shape())
  {
  case DRT::Element::hex8:
    numnodesperele = 8;
    break;
  case DRT::Element::hex20:
    numnodesperele = 20;
    cout << "Warning, the fast signed distance reinitialization has not been tested with hex20 elements!" << endl;
    break;
  case DRT::Element::hex27:
    numnodesperele = 27;
    cout << "Warning, the fast signed distance reinitialization has not been tested with hex27 elements!" << endl;
    break;
  default:
    dserror("The fast signed distance reinitialization only supports hex8, hex20 and hex27 elements.");
  }


  //========================================================================
  // get the following information about the pbc
  // - planenormaldirection e.g. (1,0,0)
  // - minimum in planenormaldirection
  // - maximum in planenormaldirection
  //========================================================================
  std::vector<DRT::Condition*>* surfacepbcs = pbc->ReturnSurfacePBCs();
  std::vector<int>    planenormal(0);
  std::vector<double> globalmins (0);
  std::vector<double> globalmaxs (0);

  for (size_t i = 0; i < surfacepbcs->size(); ++i)
  {
    const string* ismaster = (*surfacepbcs)[i]->Get<string>("Is slave periodic boundary condition");
    if (*ismaster == "Master")
    {
      const int masterid = (*surfacepbcs)[i]->GetInt("Id of periodic boundary condition");
      std::vector<int> nodeids(*((*surfacepbcs)[i]->Nodes()));
      for (size_t j = 0; j < surfacepbcs->size(); ++j)
      {
        const int slaveid = (*surfacepbcs)[j]->GetInt("Id of periodic boundary condition");
        if (masterid == slaveid)
        {
          const string* isslave = (*surfacepbcs)[j]->Get<string>("Is slave periodic boundary condition");
          if (*isslave == "Slave")
          {
            const std::vector<int>* slavenodeids = (*surfacepbcs)[j]->Nodes();
            // append slave node Ids to node Ids for the complete condition
            for (size_t k = 0; k < slavenodeids->size(); ++k)
              nodeids.push_back(slavenodeids->at(k));
          }
        }
      }

      // Get normal direction of pbc plane
      const string* pbcplane = (*surfacepbcs)[i]->Get<string>("degrees of freedom for the pbc plane");
      if (*pbcplane == "yz")
        planenormal.push_back(0);
      else if (*pbcplane == "xz")
        planenormal.push_back(1);
      else if (*pbcplane == "xy")
        planenormal.push_back(2);
      else
        dserror("A PBC condition could not provide a plane normal.");

      double min = +10e19;
      double max = -10e19;
      for (size_t j = 0; j < nodeids.size(); ++j)
      {

        const int gid = nodeids[j];
        const int lid = gfuncdis->NodeRowMap()->LID(gid);
        if (lid < 0)
          continue;
        const DRT::Node* lnode = gfuncdis->lRowNode(lid);
        const double* coord = lnode->X();
        if (coord[planenormal.back()] < min)
          min = coord[planenormal.back()];
        if (coord[planenormal.back()] > max)
          max = coord[planenormal.back()];
      }
      globalmins.resize(planenormal.size());
      globalmaxs.resize(planenormal.size());
      gfuncdis->Comm().MinAll(&min, &(globalmins.back()), 1);
      gfuncdis->Comm().MaxAll(&max, &(globalmaxs.back()), 1);
    }
  }//end loop over all surfacepbcs


  //=======================================================================
  // Create a vector of eleGIDs and a vector of those eles' node coords and
  // redundantly store it on each proc
  //=======================================================================
  std::vector<int> allcuteleids;
  std::vector<double> allnodecoords;
  {
    // Here we simply take the eleids from the boundaryIntCells map, which leads to our list of cut elements
    // also there is no distribution necessary, as this map is already stored on every proc
    for (std::map<int,GEO::BoundaryIntCells>::const_iterator elepatches = flamefront_.begin(); elepatches != flamefront_.end(); ++elepatches)
      allcuteleids.push_back(elepatches->first);

    // our local nodecoords
    std::vector<double> nodecoords( 3*numnodesperele*(allcuteleids.size()), 0.0);
    allnodecoords.resize(nodecoords.size(), 0.0);

    // write the node coordinates of every cut rownode of this proc into nodecoords
    for(size_t ivec = 0; ivec < allcuteleids.size(); ++ivec)
    {
      int elegid = allcuteleids[ivec];
      int elelid = gfuncdis->ElementRowMap()->LID(elegid);
      if (elelid >= 0)
      {
        const int coordbase = 3*numnodesperele*ivec;
        const DRT::Element* ele = gfuncdis->lRowElement(elelid);
        const DRT::Node* const* nodes = ele->Nodes();
        for(int inode = 0; inode < ele->NumNode(); ++inode)
        {
          const int nodecoordbase = coordbase + 3*inode;
          nodecoords[nodecoordbase + 0] = nodes[inode]->X()[0];
          nodecoords[nodecoordbase + 1] = nodes[inode]->X()[1];
          nodecoords[nodecoordbase + 2] = nodes[inode]->X()[2];
        }
      }
    }

    comm.SumAll(&(nodecoords[0]), &(allnodecoords[0]), (int)nodecoords.size());
  }

  //================================================================
  // loop all row nodes on the processor
  // those nodes will receive new phi values
  //================================================================
  for(int lnodeid=0; lnodeid < gfuncdis->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    const DRT::Node* lnode = gfuncdis->lRowNode(lnodeid);

    // get the dof associated with this node
    const int dofgid = gfuncdis->Dof(lnode, 0); // since this is a scalar field the dof is always 0
    int doflid = dofrowmap->LID(dofgid);
    if (doflid < 0)
      dserror("Proc %d: Cannot find dof gid=%d in Epetra_Vector",comm.MyPID(),dofgid);

    // get physical coordinates of this node
    LINALG::Matrix<3,1> nodecoord(false);
    nodecoord(0) = lnode->X()[0];
    nodecoord(1) = lnode->X()[1];
    nodecoord(2) = lnode->X()[2];
#ifdef ORACLES
    // skip the part of the domian left of the expansion
    if (nodecoord(0) <= 0.0)
      continue;
#endif
    //=======================================================================================
    // Build a list< pair< int eleGID, double distance > >
    // the distance is based on the distance between the current node and the closest node of
    // the cut element. This guarantees an estimated distance <= the real distance
    //=======================================================================================
    std::list< std::pair< int, double > > eledistance;

    {
      // loop all cut elements
      for (size_t ieleid = 0; ieleid < allcuteleids.size(); ++ieleid)
      {
        const size_t coordbase = 3*numnodesperele*ieleid;
        double distance = 1.0e19;

        // loop all cut element's nodes
        for (int inode = 0; inode < numnodesperele; ++inode)
        {
          const int nodecoordbase = coordbase + 3*inode;
          LINALG::Matrix<3,1> delta(false);
          delta(0) = allnodecoords[nodecoordbase + 0];
          delta(1) = allnodecoords[nodecoordbase + 1];
          delta(2) = allnodecoords[nodecoordbase + 2];

          delta.Update(1.0, nodecoord, -1.0);

          // take care of PBCs
          for (size_t ipbc = 0; ipbc < planenormal.size(); ++ipbc)
          {
            const double fulllength = (globalmaxs[ipbc] - globalmins[ipbc]);
            if (delta(planenormal[ipbc]) >= fulllength/2.0)
              delta(planenormal[ipbc]) = fulllength - delta(planenormal[ipbc]);
            else if (delta(planenormal[ipbc]) <= -fulllength/2.0)
              delta(planenormal[ipbc]) = delta(planenormal[ipbc]) + fulllength;
          }
          const double thisdistance = sqrt(delta(0)*delta(0) + delta(1)*delta(1) + delta(2)*delta(2));

          if (thisdistance < distance)
            distance = thisdistance;
        }

        std::pair< int, double> thispair;
        thispair.first  = allcuteleids[ieleid];
        thispair.second = distance;
        eledistance.push_back( thispair );
      }
    }
    if (eledistance.empty())
      dserror("No intersected elements available! G-function correct?");


    //==================================================================
    // sort the the vector in ascending order by the estimated distance
    //==================================================================
    // this is the STL sorting, which is pretty fast
    eledistance.sort(MyComparePairs);

    //--------------------------------------------------------------------------------
    // if a reinitbandwith is used the nodes not within the band will be set to the
    // estimated distance for all others the actual distance will be determined
    //--------------------------------------------------------------------------------
    if (!reinitband_ or (reinitband_ and fabs(eledistance.front().second) <= reinitbandwidth_))
    {

      //========================================================================
      // + update the eledistance vector with the real distance to the interface
      //   starting with the closest estimated element.
      // + Sort the vector by distance after every iteration.
      // + if the distance of the first element in the vector does not change
      //   any more, we have found the shortest distance
      //========================================================================
      int oldeleid = -1;
      while (oldeleid != eledistance.front().first) // this is just a safety check. usually loop should abort earlier.
      {
        oldeleid = eledistance.front().first;

        // the minimal distance, if all element patches and the PBCs are considered
        double pbcmindist = 1.0e19;

        // get patches belonging to first entry
        std::map<int,GEO::BoundaryIntCells>::const_iterator elepatches = flamefront_.find( eledistance.front().first );
        if (elepatches == flamefront_.end())
          dserror("Could not find the boundary integration cells belonging to Element %d.", eledistance.front().first);

        // number of flamefront patches for this element
        const std::vector<GEO::BoundaryIntCell> patches = elepatches->second;
        const int numpatch = patches.size();

        //--------------------------------------------------------------------
        // due to the PBCs the node might actually be closer to the
        // interface then would be calculated if one only considered
        // the actual position of the node. In order to find the
        // smallest distance the node is copied along all PBC directions
        //
        //   +------------------+ - - - - - - - - - -+
        //   +             II   +
        //   +   x        I  I  +    y               +
        //   +             II   +
        //   +------------------+ - - - - - - - - - -+
        //         original           copy
        //
        //   x: current node
        //   y: copy of current node
        //   I: interface
        //   +: pbc
        //--------------------------------------------------------------------
        if (planenormal.size() > 3)
          dserror("Sorry, but currently a maximum of three periodic boundary conditions are supported by the combustion reinitializer.");

        // since there is no stl pow(INT, INT) function, we calculate it manually
        size_t looplimit = 1;
        for (size_t i = 0; i < planenormal.size(); ++i)
          looplimit *= 2;

        for (size_t ipbc = 0; ipbc < looplimit; ++ipbc)
        {
          double mindist = 1.0e19;
          LINALG::Matrix<3,1> tmpcoord(nodecoord);

          // determine which pbcs have to be applied
          //
          // loopcounter | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
          // ------------+---+---+---+---+---+---+---+---+
          //  first PBC  |     x       x       x       x
          // second PBC  |         x   x           x   x
          //  third PBC  |                 x   x   x   x
          //
          // this is equivalent to the binary representation
          // of the size_t
          if (ipbc & 0x01)
          {
            const double pbclength = globalmaxs[0] - globalmins[0];
            if (nodecoord(0) < globalmins[0] or nodecoord(0) > globalmaxs[0])
              continue;
            if (nodecoord(planenormal[0]) > globalmins[0] + pbclength/2.0)
              tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) - pbclength;
            else
              tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) + pbclength;
          }
          if (ipbc & 0x02)
          {
            const double pbclength = globalmaxs[1] - globalmins[1];
            if (nodecoord(1) < globalmins[1] or nodecoord(1) > globalmaxs[1])
              continue;
            if (nodecoord(planenormal[1]) > globalmins[1] + pbclength/2.0)
              tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) - pbclength;
            else
              tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) + pbclength;
          }
          if (ipbc & 0x04)
          {
            const double pbclength = globalmaxs[2] - globalmins[2];
            if (nodecoord(2) < globalmins[2] or nodecoord(2) > globalmaxs[2])
              continue;
            if (nodecoord(planenormal[2]) > globalmins[2] + pbclength/2.0)
              tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) - pbclength;
            else
              tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) + pbclength;
          }

          //-----------------------------------------
          // loop flame front patches of this element
          //-----------------------------------------
          for(int ipatch=0; ipatch<numpatch; ++ipatch)
          {
            // get a single patch from group of flamefront patches
            const GEO::BoundaryIntCell patch = patches[ipatch];

            // only triangles and quadrangles are allowed as flame front patches (boundary cells)
            if (!(patch.Shape() == DRT::Element::tri3 or
                  patch.Shape() == DRT::Element::quad4))
            {
              dserror("invalid type of boundary integration cell for reinitialization");
            }

            // get coordinates of vertices defining flame front patch
            const LINALG::SerialDenseMatrix& patchcoord = patch.CellNodalPosXYZ();

            // compute normal vector to flame front patch
            LINALG::Matrix<3,1> normal(true);
            ComputeNormalVectorToFlameFront(patch,patchcoord,normal);

            //-----------------------------------------
            // find flame front patches facing the node
            //-----------------------------------------
            // boolean indicating if facing patch was found
            bool facenode = false;
            // distance to the facing patch
            double patchdist = 1.0e19; // default value
            // check if this patch faces the node
            FindFacingPatchProjCellSpace(tmpcoord,patch,patchcoord,normal,facenode,patchdist);

            // a facing patch was found
            if (facenode == true)
            {
              // overwrite smallest distance if computed patch distance is smaller
              if (fabs(patchdist) < fabs(mindist))
              {
                // if G-value at the node is negative, the minimal distance has to be negative
                if ((*phivector)[doflid] < 0.0 )
                  mindist = -patchdist;
                else
                  mindist = patchdist;
              }
            }

            //-------------------------------------------------------------
            // compute smallest distance to edges of this flame front patch
            //-------------------------------------------------------------
            // distance to the patch edge
            double edgedist = 1.0e19;
            ComputeDistanceToEdge(tmpcoord,patch,patchcoord,edgedist);

            if (fabs(edgedist) < fabs(mindist))
            {
              // if G-value at the node is negative, the minimal distance has to be negative
              if ((*phivector)[doflid] < 0.0 )
                mindist = -edgedist;
              else
                mindist = edgedist;
            }

            //----------------------------------------------------------------
            // compute smallest distance to vertices of this flame front patch
            //----------------------------------------------------------------
            // distance to the patch vertex
            double vertexdist = 1.0e19;
            ComputeDistanceToPatch(tmpcoord,patch,patchcoord,vertexdist);

            if (fabs(vertexdist) < fabs(mindist))
            {
              // if G-value at the node is negative, the minimal distance has to be negative
              if ((*phivector)[doflid] < 0.0 )
                mindist = -vertexdist;
              else
                mindist = vertexdist;
            }
          } // loop over flamefront patches

          if (fabs(mindist) < fabs(pbcmindist))
          {
            pbcmindist = mindist;
          }
        } // loop over PBCs

        // store the new distance, which is >= the estimated distance
        eledistance.front().second = pbcmindist;


        //==============================================================
        // sort the the vector in ascending order by the distance
        //==============================================================
        // here we use the fact, that everything is already sorted but the first list item
        std::pair<int,double> tmppair = eledistance.front();
        std::list<std::pair<int,double> >::iterator insertiter = eledistance.begin();
        ++insertiter;

        int loopcount = 0;
        // find the place where the item must be inserted
        while (fabs(tmppair.second) > fabs(insertiter->second) and insertiter != eledistance.end())
        {
          insertiter++;
          loopcount++;
        }

        // this removes the item from the front and inserts it where insertiter points to
        eledistance.splice(insertiter, eledistance, eledistance.begin());

        // if item was inserted at the beginnig of the list, it must be the shortest distance
        // possible and we can stop checking the other elements' distances
        if (loopcount == 0)
          break;

      } // loop over eledistance
    }
    // if outside the reinit band
    else
    {
      // correct the sign of estimated distance
      if ((*phivector)[doflid] < 0.0)
        eledistance.front().second = -eledistance.front().second;
    }

    //==========================================================
    // write the resulting minimal distance into the phivector
    //==========================================================
//#ifdef COMBUST_2D
//    if (fabs(nodecoord(0)) < 1.0e-8 and fabs(nodecoord(1)) < 1.0e-8 and nodecoord(2) > 0.0)
//#else
//    if (fabs(nodecoord(0)) < 1.0e-8 and fabs(nodecoord(1)) < 1.0e-8 and fabs(nodecoord(2)) < 1.0e-8)
//#endif
//    {
//      //cout << "alter phi wert in der Mitte " << (*phivector)[doflid] << endl;
//      //cout << "reinitilisierter phi wert in der Mitte " << eledistance.front().second << endl;
//      //---------------------
//      // write radius to file
//      //---------------------
//      {
//        std::string s = "/home/henke/simulations/results/study_timeint/2d_64_neumann_visc-3_theta05/radius";
//        //s.append(".oracles.horiz.chamber.flow_statistics");
//
//        // output to log-file
//        Teuchos::RCP<std::ofstream> log;
//        log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::app));
//        //(*log) << "\n";
//        (*log) <<  " "  << std::setw(12) << std::setprecision(8) << eledistance.front().second;
//        (*log) <<  " "  << std::setw(12) << std::setprecision(8) << (*phivector)[doflid];
//        (*log) << &endl;
//        log->flush();
//        cout << "wrote reinitialized radius at center to file" << endl;
//      }
//    }

    int err = phivector->ReplaceMyValues(1,&(eledistance.front().second),&doflid);
    if (err)
      dserror("this did not work");
  }

  if (comm.MyPID()==0)
    std::cout << " done" << std::endl;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: find a facing flame front patch by projecton of node into boundary cell space         |
 |                                                                                    henke 12/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::FindFacingPatchProjCellSpace(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3,1>&       normal,
    bool&                            facenode,
    double&                          patchdist)
{
  // indicator
  facenode = false;

  static LINALG::Matrix<2,1> eta(true);
  double alpha = 0.0;

  //-------------------------------------------------------
  // perform Newton-Raphson method to project node on patch
  //-------------------------------------------------------
  bool converged = false;
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    converged = ProjectNodeOnPatch<DRT::Element::tri3>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  case DRT::Element::quad4:
  {
    converged = ProjectNodeOnPatch<DRT::Element::quad4>(node, patch, patchcoord, normal, eta, alpha);
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
  }

  // Newton iteration converged
  //  cout << "Newton iteration converged in " << iter << " steps!" << endl;

  //----------------------------------------------------
  // check if projection lies within boundary cell space
  //----------------------------------------------------
  // remark: - tolerance has to be of same order as the tolerance that coordinates of projected nodes
  //           differ from an exact position on edges of patches (e.g. 1.0E-7 ~ 1.0E-8 -> 1.0E-6)
  //         - if this is not the case, the level set function can become tilted, since valid
  //           patches are ignored
  double TOL= 1e-6;
  
  switch(patch.Shape())
  {
  case DRT::Element::tri3:
  {
    // criteria for tri3 patch
    if ((eta(0) > -TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -TOL) and (eta(1) < 1.0+TOL) and
        (1.0-eta(0)-eta(1) > -TOL) and (1.0-eta(0)-eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = fabs(alpha);
//      cout << "facing patch found (tri3 patch)! coordinates eta(0): " << eta(0) << " eta(1) " << eta(1) << endl;
    }
    break;
  }
  case DRT::Element::quad4:
  {
    // criteria for quad4 patch
    if ((eta(0) > -1.0-TOL) and (eta(0) < 1.0+TOL) and
        (eta(1) > -1.0-TOL) and (eta(1) < 1.0+TOL) and
        converged)
    {
      facenode = true;
      patchdist = fabs(alpha);
//      cout << "facing patch found (quad4 patch)!" << endl;
    }
    break;
  }
  default:
    dserror("unknown type of boundary integration cell");
  }
//  if (!converged)
//  {
//    cout << "node x component " << node(0,0) << endl;
//    cout << "node y component " << node(1,0) << endl;
//    cout << "node z component " << node(2,0) << endl;
//    cout << "eta1 " << eta(0) << endl;
//    cout << "eta2 " << eta(1) << endl;
//    cout << "alpha " << alpha << endl;
//    cout << "patch vertices x component " << patchcoord(0,0) << " " << patchcoord(0,1) << " " << patchcoord(0,2) << endl;
//    cout << "patch vertices y component " << patchcoord(1,0) << " " << patchcoord(1,1) << " " << patchcoord(1,2) << endl;
//    cout << "patch vertices z component " << patchcoord(2,0) << " " << patchcoord(2,1) << " " << patchcoord(2,2) << endl;
//  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute distance to edge of patch                                         henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::ComputeDistanceToEdge(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    double&                          edgedist
)
{
  // set temporary edgedist to large value
  double edgedisttmp = edgedist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch (first vertex)
  static LINALG::Matrix<3,1> vertex1(true);
  // current next vertex of the patch (second vertex)
  static LINALG::Matrix<3,1> vertex2(true);
  // distance vector from first vertex to node
  static LINALG::Matrix<3,1> vertex1tonode(true);
  // distance vector from first vertex to second vertex
  static LINALG::Matrix<3,1> vertex1tovertex2(true);

  // compute distance to all vertices of patch
  for(size_t ivert = 0; ivert<numvertices; ++ivert)
  {
    // vertex1 of flame front patch
    vertex1(0) = patchcoord(0,ivert);
    vertex1(1) = patchcoord(1,ivert);
    vertex1(2) = patchcoord(2,ivert);

    if (ivert < (numvertices-1))
    {
      vertex2(0) = patchcoord(0,ivert+1);
      vertex2(1) = patchcoord(1,ivert+1);
      vertex2(2) = patchcoord(2,ivert+1);
    }
    else if (ivert == (numvertices-1))
    {
      vertex2(0) = patchcoord(0,0);
      vertex2(1) = patchcoord(1,0);
      vertex2(2) = patchcoord(2,0);
    }

    // compute distance vector from node to current first
    vertex1tonode.Update(1.0, node, -1.0, vertex1);
    // compute distance vector from current second first vertex to current frist vertex (edge)
    vertex1tovertex2.Update(1.0, vertex2, -1.0, vertex1);
    double normvertex1tovertex2 = vertex1tovertex2.Norm2();
    // normalize vector
    vertex1tovertex2.Scale(1.0/normvertex1tovertex2);

    // scalar product of vertex1tonode and the normed vertex1tovertex2
    double lotfusspointdist = vertex1tovertex2.Dot(vertex1tonode);

    if( (lotfusspointdist >= 0.0) and (lotfusspointdist <= normvertex1tovertex2) ) // lotfusspoint on edge
    {
      LINALG::Matrix<3,1> lotfusspoint(true);
      lotfusspoint.Update(1.0,vertex1,lotfusspointdist,vertex1tovertex2);
      LINALG::Matrix<3,1> nodetolotfusspoint(true);
      nodetolotfusspoint.Update(1.0,lotfusspoint,-1.0,node);

      // determine length of vector from node to lot fuss point
      edgedisttmp = nodetolotfusspoint.Norm2();
      if (edgedisttmp < edgedist)
        edgedist = edgedisttmp;
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute distance to vertex of patch                                       henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::ComputeDistanceToPatch(
    const LINALG::Matrix<3,1>&       node,
    const GEO::BoundaryIntCell&      patch,
    const LINALG::SerialDenseMatrix& patchcoord,
    double&                          vertexdist
)
{
  // set temporary vertexdist to large value
  double vertexdisttmp = vertexdist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch
  static LINALG::Matrix<3,1> vertex(true);
  // distance vector from patch to node
  static LINALG::Matrix<3,1> dist(true);

  // compute distance to all vertices of patch
  for(size_t ivert = 0; ivert<numvertices; ++ivert)
  {
    // vertex of flame front patch
    vertex(0) = patchcoord(0,ivert);
    vertex(1) = patchcoord(1,ivert);
    vertex(2) = patchcoord(2,ivert);

    // compute distance vector from flame front to node
    dist.Update(1.0, node, -1.0, vertex);

    // compute L2-norm of distance vector
    vertexdisttmp = dist.Norm2();
    if (vertexdisttmp < vertexdist)
      vertexdist = vertexdisttmp;
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: compute normal vector to flame front patch                                henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void COMBUST::Reinitializer::ComputeNormalVectorToFlameFront(
      const GEO::BoundaryIntCell&      patch,
      const LINALG::SerialDenseMatrix& patchcoord,
      LINALG::Matrix<3,1>&             normal)
{
  // first point of flame front patch
  LINALG::Matrix<3,1> point1;
  point1(0) = patchcoord(0,0);
  point1(1) = patchcoord(1,0);
  point1(2) = patchcoord(2,0);

  // second point of flame front patch
  LINALG::Matrix<3,1> point2;
  point2(0) = patchcoord(0,1);
  point2(1) = patchcoord(1,1);
  point2(2) = patchcoord(2,1);

  // first edge of flame front patch
  LINALG::Matrix<3,1> edge1;
  edge1.Update(1.0, point2, -1.0, point1);

  // third point of flame front patch
  point2(0) = patchcoord(0,2);
  point2(1) = patchcoord(1,2);
  point2(2) = patchcoord(2,2);

  // second edge of flame front patch (if patch is triangle; if not: edge 2 is secant of polygon)
  LINALG::Matrix<3,1> edge2;
  edge2.Update(1.0, point2, -1.0, point1);

  // compute normal vector of patch (cross product: edge1 x edge2)
  // remark: normal vector points into unburnt domain (G<0)
  normal(0) = (edge1(1)*edge2(2) - edge1(2)*edge2(1));
  normal(1) = (edge1(2)*edge2(0) - edge1(0)*edge2(2));
  normal(2) = (edge1(0)*edge2(1) - edge1(1)*edge2(0));
#ifdef COMBUST_2D
  normal(2) = 0.0;
#endif

//  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
//  cout << "proc " << comm.MyPID() << " normal " <<  normal << endl;
//  cout << "proc " << comm.MyPID() << " patch " <<  patchcoord << endl;

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0/norm);

#ifdef DEBUG
#ifdef COMBUST_2D
  if (!((normal(2) > 0.0-1.0E-8) and (normal(2) < 0.0+1.0E-8)))
  {
    cout << "z-component of normal: " << normal(2) << endl;
    dserror ("pseudo-3D problem not symmetric anymore!");
  }
#endif
#endif

  return;
}


/*------------------------------------------------------------------------------------------------*
 | private: find facing flame front patch by projection on the patch egdes            henke 08/09 |
 *----------------------------------------------------------------------------------------------- */
void::COMBUST::Reinitializer::FindFacingPatchProjCellEdges(
       const LINALG::Matrix<3,1>&       node,
       const LINALG::SerialDenseMatrix& patch,
       bool&                            facenode,
       double&                          vertexdist)
{
  dserror("Thou shalt not call this function! It is out of order.");
#if 0
  facenode = false;
  // number of column vectors (edges) forming flame front patch
  const int numcol = patch.N();

  // beginning point of flame front patch
  LINALG::Matrix<3,1> vecbeg(true);
  // end point of flame front patch
  LINALG::Matrix<3,1> vecend(true);
  // edge vector of flame front patch
  LINALG::Matrix<3,1> edge(true);
  // distance vector from beginning point to node
  LINALG::Matrix<3,1> dist(true);

  int counter = 0;

  for(int iedge = 0; iedge<numcol; ++iedge)
  {
    // beginning point of flame front patch
    vecbeg(0) = patch(0,iedge);
    vecbeg(1) = patch(1,iedge);
    vecbeg(2) = patch(2,iedge);

    if(iedge == (numcol-1)) // last edge
    {
      // end point of edge is first point of flame front patch
      vecend(0) = patch(0,0);
      vecend(1) = patch(1,0);
      vecend(2) = patch(2,0);
    }
    else // regular edge (not the last one)
    {
      // end point of flame front patch
      vecend(0) = patch(0,iedge+1);
      vecend(1) = patch(1,iedge+1);
      vecend(2) = patch(2,iedge+1);
    }

    // compute edge vector of flame front patch
    edge.Update(1.0, vecend, -1.0, vecbeg);
    // compute distance vector from node to flame front
    dist.Update(1.0, vecend, -1.0, node);

    // compute norms of edge and distance vectors
    double normedge = sqrt(edge(0)*edge(0) + edge(1)*edge(1) + edge(2)*edge(2));
    double normdist = sqrt(dist(0)*dist(0) + dist(1)*dist(1) + dist(2)*dist(2));
    if(normdist < fabs(vertexdist))
    {
      // compute normal vector of flame front patch
      LINALG::Matrix<3,1> normal(true);
      ComputeNormalVectorToFlameFront(patch, normal);

      // determine sign of distance vector from node to flame front
      double tmp = normal(0)*dist(0) + normal(1)*dist(1) + normal(2)*dist(2);
      // tmp < 0 if node in unburnt domain (G<0) and vice versa
      if (tmp <= 0.0)
        vertexdist = -normdist;
      else
        vertexdist = normdist;

      //cout << "distance to vertex is overwritten by: " << vertexdist << endl;
    }

    // divide vectors by their norms to get unit vectors
    if (normedge == 0.0) dserror("length of edge vector is zero");
    for (int icomp=0; icomp<3; ++icomp) edge(icomp) /= normedge;
    if (normdist == 0.0) dserror("length of distance vector is zero");
    for (int icomp=0; icomp<3; ++icomp) dist(icomp) /= normdist;

    // compute projection of distance vector in direction of edge
    double proj = edge(0)*dist(0) + edge(1)*dist(1) + edge(2)*dist(2);

    // if (0 < projection < 1), node faces this edge of the patch
    if((proj >= 0.0-1.0E-13) and (proj <= 1.0+1.0E-13))
      counter += 1;
    // if all projections face their edges, the node faces the patch
    if(counter == numcol)
      facenode = true;
  }
#endif
  return;
}

