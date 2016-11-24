/*!
\file dof_management.cpp

\brief provides the dofmanager class

<pre>
\level 2
\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */


#include "dof_management.H"
#include "dofkey.H"
#include "xdofmapcreation_combust.H"
#include "../drt_combust/combust_interface.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset.H"


/*------------------------------------------------------------------------------------------------*
 | constructor: used for combustion problems only                                     henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
XFEM::DofManager::DofManager(
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle,
    const Teuchos::RCP<Epetra_Vector>&                   phinp,
    const std::set<XFEM::PHYSICS::Field>&                fieldset,
    const Teuchos::ParameterList&                        params,
    const Teuchos::RCP<std::map<int,std::vector<int> > >           pbcmap
    ) :
  ih_(interfacehandle),
  pbcmap_(pbcmap)
{
  //if (ih_->FluidDis()->Comm().MyPID() == 0)
  //  std::cout << "Constructing DofManager for combustion problem" << std::endl;

  std::map<int, std::set<XFEM::FieldEnr> >    nodeDofMap;

  // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
  XFEM::createDofMapCombust(
      *interfacehandle,
      phinp.get(),
      nodeDofMap,
      fieldset,
      params);

  //------------------------------------
  // connect dofs on periodic boundaries
  //------------------------------------
  if (pbcmap_ != Teuchos::null)
  // remark: - for any regular call of DofManager the pbc map exists (pbcmap_ != Teuchos::null),
  //         - for the initialization call of DofManager in CombustFluidImplicitTimeInt constructor the
  //           pbc map does not make sence yet, since no XFEM dofs have been assigned
  {
    // we think that this works reliably for periodic boundary conditions in all cases (henke 18.7.2011)
    // remark: - enrichments are distributed in a loop over all column elements
    //         - periodic boundary conditions live on row (or colunm) nodes
    //         -> for a general parallel distribution 'pbcmap' and 'nodeDofMap' do not contain the same information
    //         -> simulation will crash as soon this discrepancy shows

    for (std::map<int, std::vector<int>  >::const_iterator pbciter= (*pbcmap_).begin(); pbciter != (*pbcmap_).end(); ++pbciter)
    {
      const int mastergid = pbciter->first;

      for (size_t islave = 0; islave < pbciter->second.size(); islave++)
      {
        const int slavegid  = pbciter->second[islave];
        //std::cout << "proc " << ih_->FluidDis()->Comm().MyPID() << " mastergid " << mastergid << " slavegid " << slavegid << std::endl;

        std::set<FieldEnr> masterfieldenr = emptyset_;
        std::set<FieldEnr> slavefieldenr = emptyset_;

        std::map<int, std::set<XFEM::FieldEnr> >::iterator masterentry = nodeDofMap.find(mastergid);
        if (masterentry != nodeDofMap.end())
          masterfieldenr = masterentry->second;
        else dserror("pbc master node not found in node dof map");

        std::map<int, std::set<XFEM::FieldEnr> >::iterator slaveentry = nodeDofMap.find(slavegid);
        if (slaveentry != nodeDofMap.end())
          slavefieldenr = slaveentry->second;
        else
          dserror("pbc slave node not found in node dof map");

        //if (slaveentry->second != masterentry->second)
        //{
        //  std::cout << mastergid << " das passt nicht zusammen" << std::endl;
        //  std::cout << "proc " << ih_->FluidDis()->Comm().MyPID() << " master dofs " << masterfieldenr.size() << std::endl;
        //  std::cout << "proc " << ih_->FluidDis()->Comm().MyPID() << " slave dofs " << slavefieldenr.size() << std::endl;
        //}

        // if slave is enriched, but master is not (slave has more dofs)
        if (slavefieldenr.size()>masterfieldenr.size())
          masterentry->second = slaveentry->second;
      }

      //by now the master has the most dofs and must be copied to all slaves
      std::set<FieldEnr> masterfieldenr = emptyset_;
      std::map<int, std::set<XFEM::FieldEnr> >::iterator masterentry = nodeDofMap.find(mastergid);

      for (size_t islave = 0; islave < pbciter->second.size(); islave++)
      {
        const int slavegid = pbciter->second[islave];
        std::set<FieldEnr> slavefieldenr = emptyset_;
        std::map<int, std::set<XFEM::FieldEnr> >::iterator slaveentry = nodeDofMap.find(slavegid);

        slaveentry->second = masterentry->second;
      }
    }
  }

  //------------------------------------------------------------------------------------------------
  // copy non-const dof maps to const dof maps
  // remark 1: the whole map should actually be const, not just the second element (a set)
  // remark 2: here, we need the general XFEM dof maps to be filled with the combustion dof maps,
  //           because member functions like "toGmsh" and "GatherUniqueEnrichments" rely on them
  //
  // what the following code does, is basically:   nodalDofSet_ = nodeDofMap
  //------------------------------------------------------------------------------------------------
  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodeDofMap.begin(); oneset != nodeDofMap.end(); ++oneset )
  {
    nodalDofSet_.insert( make_pair(oneset->first, oneset->second));
  };

  // collect all unique enrichments and print them on the screen (only needed for verification)
  GatherUniqueEnrichments();

  //if (ih_->FluidDis()->Comm().MyPID() == 0)
  //  std::cout << "Constructing DofManager for combustion problem done" << std::endl;
}

/*------------------------------------------------------------------------------------------------*
 | gather all unique enrichments and print them on screen (Debug)                        ag 11/07 |
 *----------------------------------------------------------- ------------------------------------*/
void XFEM::DofManager::GatherUniqueEnrichments() const
{
  // set of unique enrichments
  std::set<XFEM::Enrichment> unique_enrichments;

  // collect enrichments from nodal dofs
  for (std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator fieldenriter=nodalDofSet_.begin();
  fieldenriter!=nodalDofSet_.end(); ++fieldenriter)
  {
    const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
      enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      unique_enrichments.insert(enrfield->getEnrichment());
    }
  }


#ifdef DEBUG
  // screen output
  IO::cout << " Enrichments available:";
  for (std::set<XFEM::Enrichment>::const_iterator enr =
      unique_enrichments.begin(); enr != unique_enrichments.end(); ++enr)
  {
    IO::cout << " " << enr->toString();
  }
  IO::cout << IO::endl;
#endif
}



/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const XFEM::DofManager&)
{
  dserror("A DofManager shall not be copied!");
}



/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::~DofManager()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int XFEM::DofManager::NumNodalDof() const
{
  int locnumnodaldof = 0;

  std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator iter;
  for (iter = nodalDofSet_.begin(); iter != nodalDofSet_.end(); ++iter)
  {
    locnumnodaldof += iter->second.size();
  }

  // collect number of nodal dofs from all procs
  int numnodaldof = 0;
  ih_->FluidDis()->Comm().SumAll(&locnumnodaldof,&numnodaldof,1);

  return numnodaldof;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::DofManager::fillDofRowDistributionMaps(
    std::map<XFEM::DofKey, XFEM::DofGID>&  NodalDofDistributionMap
) const
{
  NodalDofDistributionMap.clear();
  // loop all (non-overlapping = Row)-Nodes and store the DOF information w.t.h. of DofKeys
  for (int i=0; i<ih_->FluidDis()->NumMyRowNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->FluidDis()->lRowNode(i);
    const int gid = actnode->Id();
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
    if (entry == nodalDofSet_.end())
    {
      // no dofs for this node... must be a hole or somethin'
      continue;
    }
    const std::vector<int> gdofs(ih_->FluidDis()->Dof(actnode));
    const std::set<FieldEnr> dofset = entry->second;
    if (gdofs.size() != dofset.size())
    {
      std::cout << "numdof node (Discretization): " <<  gdofs.size() << std::endl;
      std::cout << "numdof node (DofManager):     " <<  dofset.size() << std::endl;
      dserror("Bug!!! Information about nodal dofs in DofManager and Discretization does not fit together!");
    }

    int dofcount = 0;
    std::set<FieldEnr>::const_iterator fieldenr;
    for(fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
    {
      NodalDofDistributionMap.insert(std::make_pair(DofKey(gid, *fieldenr), gdofs[dofcount]));
      dofcount++;
    }
  };
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::DofManager::fillNodalDofColDistributionMap(
    std::map<XFEM::DofKey, XFEM::DofGID>&  NodalDofColDistributionMap
) const
{
  NodalDofColDistributionMap.clear();
  // loop all (overlapping = Col)-Nodes and store the DOF information w.t.h. of DofKeys
  for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->FluidDis()->lColNode(i);
    const int gid = actnode->Id();
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
    if (entry == nodalDofSet_.end())
    {
      // no dofs for this node... must be a hole or somethin'
      continue;
    }
    const std::vector<int> gdofs(ih_->FluidDis()->Dof(actnode));
    const std::set<FieldEnr> dofset = entry->second;
    if (gdofs.size() != dofset.size())
    {
      std::cout << "proc " << ih_->FluidDis()->Comm().MyPID() << " node " << actnode->Id() << std::endl;
      std::cout << "numdof node (Discretization): " <<  gdofs.size() << std::endl;
      std::cout << "numdof node (DofManager):     " <<  dofset.size() << std::endl;
      dserror("Bug!!! Information about nodal dofs in DofManager and Discretization does not fit together!");
    }

    int dofcount = 0;
    std::set<FieldEnr>::const_iterator fieldenr;
    for(fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
    {
      NodalDofColDistributionMap.insert(std::make_pair(DofKey(gid, *fieldenr), gdofs[dofcount]));
      dofcount++;
    }
  }
}


/*----------------------------------------------------------------------*
 | transform XFEM vector to (standard FEM) output vector     henke 07/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::DofManager::transformXFEMtoStandardVector(
    const Epetra_Vector&                   xfemvector,
    const DRT::DofSet&                     outdofset,
    const std::map<DofKey, DofGID>&        nodalDofDistributionMap,
    const std::set<XFEM::PHYSICS::Field>&  outputfields
) const
{
  // get DofRowMaps for output and XFEM vectors
  const Epetra_Map* outdofrowmap = outdofset.DofRowMap();
  const Epetra_BlockMap* xfemdofrowmap = &xfemvector.Map();

//  const Epetra_Map* xfemdofrowmap = ih_->FluidDis()->DofRowMap();
  // create output vector (standard FEM layout)
  Teuchos::RCP<Epetra_Vector> outvector = Teuchos::rcp(new Epetra_Vector(*outdofrowmap,true));

  // loop nodes on this processor
  for (int inode=0; inode<ih_->FluidDis()->NumMyRowNodes(); ++inode)
  {
    // get XFEM GID of this node
    const DRT::Node* xfemnode = ih_->FluidDis()->lRowNode(inode);
    const int nodegid = xfemnode->Id();
    // get vector of dof GIDs for this node according to output (standard FEM) layout
    const std::vector<int> outgid(outdofset.Dof(xfemnode));
    // find the set of field enrichments (~ XFEM dofs) for this node
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator nodeentry = nodalDofSet_.find(nodegid);

    // node was not found in nodalDofSet_
    if (nodeentry == nodalDofSet_.end())
    {
      dserror("Every node should have (at least a standard) field enrichments!");
      // one dof for every physical field (standard FEM)
      //const std::size_t numdof = outputfields.size();
      // write zero values in output vector at position of original (standard FEM) degree of freedom
      //for (std::size_t idof = 0; idof < numdof; ++idof)
      //{
      //  (*outvector)[outdofrowmap->LID(outgid[idof])] = 0.0;
      //}
    }
    else // node was found in nodalDofSet_
    {
      // get set of field enrichments for this node
      const std::set<FieldEnr> fieldenrset = nodeentry->second;
      // build a standard enrichment (label = 0)
      const XFEM::Enrichment stdenr(XFEM::Enrichment::typeStandard,0);

      size_t idof = 0;
      // loop over desired physical output fields
      for(std::set<XFEM::PHYSICS::Field>::const_iterator outputfield = outputfields.begin(); outputfield != outputfields.end(); ++outputfield)
      {
        // build a standard field enrichment with this physical field
        XFEM::FieldEnr fieldstdenr(*outputfield,stdenr);
        // check if there is a standard enrichment for this physical field available for this node
        std::set<FieldEnr>::const_iterator fieldenrentry = fieldenrset.find(fieldstdenr);

        // there is no standard enrichment for this desired output field
        if (fieldenrentry == fieldenrset.end())
        {
          dserror("There should be a standard enrichment for every physical output field!");
        }
        else // there is a standard enrichment for this desired output field
        {
          // build a dofkey (= XFEM dof)
          const XFEM::DofKey dofkey(nodegid,*fieldenrentry);
          // get dof GID (XFEM layout) corresponding to this dofkey
          const int xfemgid = nodalDofDistributionMap.find(dofkey)->second;
          if (xfemgid < 0)
            dserror("bug!");
          if (outgid[idof] < 0)
            dserror("bug!");
          // fill output vector by writing value of standard enrichment (XFEM)
          (*outvector)[outdofrowmap->LID(outgid[idof])] = xfemvector[xfemdofrowmap->LID(xfemgid)];
        }
//        // loop over field enrichments of this node
//        for(std::set<FieldEnr>::const_iterator fieldenr = fieldenrset.begin(); fieldenr != fieldenrset.end(); ++fieldenr)
//        {
//          // get the physical field and enrichment type of this field enrichment
//          const XFEM::PHYSICS::Field xfemfield = fieldenr->getField();
//          const XFEM::Enrichment::EnrType xfemenrtype = fieldenr->getEnrichment().Type();
//          // if node has a field enrichment which is also an output field and a standard dof
//          if ((xfemfield == *outputfield) and (xfemenrtype == XFEM::Enrichment::typeStandard))
//          {
//            // build a dofkey (= XFEM dof)
//            const XFEM::DofKey dofkey(nodegid,*fieldenr);
//            // get dof GID (XFEM layout) corresponding to this dofkey
//            const int xfemgid = nodalDofDistributionMap.find(dofkey)->second;
//            if (xfemgid < 0)
//              dserror("bug!");
//            if (outgid[idof] < 0)
//              dserror("bug!");
//            // fill output vector by writing value of standard enrichment (XFEM)
//            (*outvector)[outdofrowmap->LID(outgid[idof])] = xfemvector[xfemdofrowmap->LID(xfemgid)];
//          }
//        }
        idof++;
      }
#ifdef DEBUG
      if (idof != outgid.size())
        // this is not really an error, but it would be unusual to intentionally skip fields for output
        dserror("Not all available fields have been transformed to output vector!");
#endif
    }
  };
  return outvector;
}


/*----------------------------------------------------------------------*
 | transform standard vector to XFEM vector             rasthofer 03/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::DofManager::transformStandardToXFEMVector(
    const Epetra_Vector&                   stdvector,
    const DRT::DofSet&                     stddofset,
    const std::map<DofKey, DofGID>&        nodalDofDistributionMap,
    const std::set<XFEM::PHYSICS::Field>&  outputfields,
    const Epetra_BlockMap*                 xfemdofrowmap
) const
{
  // get DofRowMaps for output and XFEM vectors
  const Epetra_Map* stddofrowmap = stddofset.DofRowMap();

  // create XFEM vector (initialized with zeros)
  Teuchos::RCP<Epetra_Vector> xfemvector = Teuchos::rcp(new Epetra_Vector(*xfemdofrowmap,true));
  xfemvector->PutScalar(0.0);

  // loop nodes on this processor
  for (int inode=0; inode<ih_->FluidDis()->NumMyRowNodes(); ++inode)
  {
    // get XFEM GID of this node
    const DRT::Node* xfemnode = ih_->FluidDis()->lRowNode(inode);
    const int nodegid = xfemnode->Id();

    // get vector of dof GIDs for this node according standard FEM layout
    const std::vector<int> stdgid(stddofset.Dof(xfemnode));
    // find the set of field enrichments (~ XFEM dofs) for this node
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator nodeentry = nodalDofSet_.find(nodegid);

    // node was not found in nodalDofSet_
    if (nodeentry == nodalDofSet_.end())
    {
      dserror("Every node should have (at least a standard) field enrichments!");
      // one dof for every physical field (standard FEM)
      //const std::size_t numdof = outputfields.size();
      // write zero values in output vector at position of original (standard FEM) degree of freedom
      //for (std::size_t idof = 0; idof < numdof; ++idof)
      //{
      //  (*outvector)[outdofrowmap->LID(outgid[idof])] = 0.0;
      //}
    }
    else // node was found in nodalDofSet_
    {
      // get set of field enrichments for this node
      const std::set<FieldEnr> fieldenrset = nodeentry->second;
      // build a standard enrichment (label = 0)
      const XFEM::Enrichment stdenr(XFEM::Enrichment::typeStandard,0);

      size_t idof = 0;
      // loop over desired physical output fields
      for(std::set<XFEM::PHYSICS::Field>::const_iterator outputfield = outputfields.begin(); outputfield != outputfields.end(); ++outputfield)
      {
        // build a standard field enrichment with this physical field
        XFEM::FieldEnr fieldstdenr(*outputfield,stdenr);
        // check if there is a standard enrichment for this physical field available for this node
        std::set<FieldEnr>::const_iterator fieldenrentry = fieldenrset.find(fieldstdenr);

        // there is no standard enrichment for this desired output field
        if (fieldenrentry == fieldenrset.end())
        {
          dserror("There should be a standard enrichment for every physical output field!");
        }
        else // there is a standard enrichment for this desired output field
        {
          // build a dofkey (= XFEM dof)
          const XFEM::DofKey dofkey(nodegid,*fieldenrentry);
          // get dof GID (XFEM layout) corresponding to this dofkey
          const int xfemgid = nodalDofDistributionMap.find(dofkey)->second;
          if (xfemgid < 0)
            dserror("bug!");
          if (stdgid[idof] < 0)
            dserror("bug!");
          // fill output vector by writing value of standard enrichment (XFEM)

          (*xfemvector)[xfemdofrowmap->LID(xfemgid)]=stdvector[stddofrowmap->LID(stdgid[idof])];
        }
        idof++;
      }
#ifdef DEBUG
      // TAW: error: ‘outgid’ was not declared in this scope
      //if (idof != outgid.size())
      //  // this is not really an error, but it would be unusual to intentionally skip fields for output
      //  dserror("Not all available fields have been transformed to output vector!");
#endif
    }
  };
  return xfemvector;
}


/*----------------------------------------------------------------------*
 | write values of a certain physical field              rasthofer 09/11|
 |                                                          DA wichmann |
 *----------------------------------------------------------------------*/
void XFEM::DofManager::overwritePhysicalField(
    Teuchos::RCP<Epetra_Vector>&           vector,
    const std::map<DofKey, DofGID>&        nodalDofDistributionMap,
    const XFEM::PHYSICS::Field&            physfield,
    const XFEM::Enrichment::EnrType&       enrichment,
    const double                           value,
    const bool                             strict
) const
{
  // loop nodes on this processor
  for (int inode=0; inode < ih_->FluidDis()->NumMyRowNodes(); ++inode)
  {
    // get GID of this node
    const DRT::Node* node = ih_->FluidDis()->lRowNode(inode);
    const int nodegid = node->Id();

    // find the set of field enrichments (~ XFEM dofs) for this node
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator nodeentry = nodalDofSet_.find(nodegid);

    // node was not found in nodalDofSet_
    if (nodeentry == nodalDofSet_.end())
    {
      dserror("Every node should have (at least a standard) field enrichments!");
    }
    else // node was found in nodalDofSet_
    {
      // get set of field enrichments for this node
      const std::set<FieldEnr> fieldenrset = nodeentry->second;

      // if strict is enabled we only want dofs of nodes which use only the specified enrichment
      if (strict)
      {
        bool haswrongenrichment = false;
        // loop over physical output fields
        for(std::set<FieldEnr>::const_iterator ifield = fieldenrset.begin(); ifield != fieldenrset.end(); ++ifield)
        {
          if (ifield->getEnrichment().Type() != enrichment)
          {
            haswrongenrichment = true;
            break;
          }
        }
        // continue with next node
        if (haswrongenrichment)
          continue;
      }

      // loop over physical output fields
      for(std::set<FieldEnr>::const_iterator ifield = fieldenrset.begin(); ifield != fieldenrset.end(); ++ifield)
      {
        if (ifield->getEnrichment().Type() == enrichment)
        {
          if (ifield->getField() == physfield)
          {
            // build a dofkey (= XFEM dof)
            const XFEM::DofKey dofkey(nodegid, *ifield);
            // get dof GID (XFEM layout) corresponding to this dofkey
            std::map<DofKey, DofGID>::const_iterator idof = nodalDofDistributionMap.find(dofkey);
            if (idof == nodalDofDistributionMap.end())
              dserror("bug!");
            const int dofgid = idof->second;
            if (dofgid < 0)
              dserror("bug!");
            // fill output vector by writing value of standard enrichment (XFEM)
            vector->ReplaceGlobalValue(dofgid, 0, value);
          }
        }
      } // loop over physical output fields
    }
  } // loop nodes on this processor
  return;
}


/*----------------------------------------------------------------------*
 |  to Gmsh (Debug)                                             ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::toGmsh(
    const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
#if 1
  const bool screen_out = DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;

  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("numdof_coupled_system", step, 5, screen_out, ih_->FluidDis()->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // draw elements with associated gid
      gmshfilecontent << "View \" " << "Element->Id() \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = ih_->FluidDis()->lColElement(i);
        IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
    {
      gmshfilecontent << "View \" " << "Node->Id() \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = ih_->FluidDis()->lColNode(i);
        const LINALG::Matrix<3,1> pos(actnode->X());
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    {
      gmshfilecontent << "View \" " << "NumDof per node \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->FluidDis()->lColNode(i);

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub =
            nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fieldenrset = blub->second;
          const LINALG::Matrix<3,1> pos(xfemnode->X());
          IO::GMSH::cellWithScalarToStream(DRT::Element::point1, fieldenrset.size(), pos, gmshfilecontent);
        }
      };
      gmshfilecontent << "};\n";
    }

    {
      gmshfilecontent << "View \" " << "NumDof Jump enriched nodes \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->FluidDis()->lColNode(i);

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          double val = 0.0;
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if (f->getEnrichment().Type() == XFEM::Enrichment::typeJump)
            {
              val += 1.0;
            }
          }
          if (val > 0.5)
          {
            const LINALG::Matrix<3,1> pos(xfemnode->X());
            IO::GMSH::cellWithScalarToStream(DRT::Element::point1, val, pos, gmshfilecontent);
          }
        }
      };
      gmshfilecontent << "};\n";
    }

    {
      gmshfilecontent << "View \" " << "NumDof" << " standard enriched nodes \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->FluidDis()->lColNode(i);

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          double val = 0.0;
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if (f->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
            {
              val += 1.0;
            }
          }
          if (val > 0.5)
          {
            const LINALG::Matrix<3,1> pos(xfemnode->X());
            IO::GMSH::cellWithScalarToStream(DRT::Element::point1, val, pos, gmshfilecontent);
          }
        }
      };
      gmshfilecontent << "};\n";
    }

    {
      gmshfilecontent << "View \" " << "NumDof" << " Void enriched nodes \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->FluidDis()->lColNode(i);

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          double val = 0.0;
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if (f->getEnrichment().Type() == XFEM::Enrichment::typeVoid)
            {
              val += 1.0;
            }
          }
          if (val > 0.5)
          {
            const LINALG::Matrix<3,1> pos(xfemnode->X());
            IO::GMSH::cellWithScalarToStream(DRT::Element::point1, val, pos, gmshfilecontent);
          }
        }
      };
      gmshfilecontent << "};\n";
    }

    {
      gmshfilecontent << "View \" " << "NumDof" << " Kink enriched nodes \" {\n";
      for (int i=0; i<ih_->FluidDis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->FluidDis()->lColNode(i);

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          double val = 0.0;
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if (f->getEnrichment().Type() == XFEM::Enrichment::typeKink)
            {
              val += 1.0;
            }
          }
          if (val > 0.5)
          {
            const LINALG::Matrix<3,1> pos(xfemnode->X());
            IO::GMSH::cellWithScalarToStream(DRT::Element::point1, val, pos, gmshfilecontent);
          }
        }
      };
      gmshfilecontent << "};\n";
    }

    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << std::endl;
  }
#if 1
  if (gmshdebugout)
  {
      const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("eledofman_check", step, 5, screen_out, ih_->FluidDis()->Comm().MyPID());
      std::ofstream gmshfilecontent(filename.c_str());
      {
        gmshfilecontent << "View \" " << " NumDofPerElement() in element \" {\n";
        for (int i=0; i<ih_->FluidDis()->NumMyColElements(); ++i)
        {
          DRT::Element* actele = ih_->FluidDis()->lColElement(i);
          const double val = actele->NumDofPerElement();
          IO::GMSH::elementAtInitialPositionToStream(val, actele, gmshfilecontent);
        };
        gmshfilecontent << "};\n";
      }
      gmshfilecontent.close();
      if (screen_out) std::cout << " done" << std::endl;
  }
#endif
#endif
}
