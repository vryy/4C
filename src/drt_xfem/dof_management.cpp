/*!
\file dof_management.cpp

\brief provides the dofmanager class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "dof_management.H"
#include "xdofmapcreation.H"
#include "enrichment_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>



/*------------------------------------------------------------------------------------------------*
 | constructor: used for xfsi problems                                                   ag 11/07 |
 *------------------------------------------------------------------------------------------------*/
XFEM::DofManager::DofManager(
    const RCP<XFEM::InterfaceHandle>&      ih,
    const std::set<XFEM::PHYSICS::Field>&  fieldset,
    const XFEM::ElementAnsatz&             element_ansatz,
    const Teuchos::ParameterList&          params
    ) :
  ih_(ih)
{
  XFEM::createDofMap(*ih, nodalDofSet_, elementalDofs_, fieldset, element_ansatz, params);

  GatherUniqueEnrichments();
}

/*------------------------------------------------------------------------------------------------*
 | constructor: used for combustion problems only                                     henke 03/09 |
 *------------------------------------------------------------------------------------------------*/
XFEM::DofManager::DofManager(
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle,
    const std::set<XFEM::PHYSICS::Field>&      fieldset,
    const Teuchos::ParameterList&              params
    ) :
  ih_(interfacehandle)
{
  //if (ih_->xfemdis()->Comm().MyPID() == 0)
  //  std::cout << "Constructing DofManager for combustion problem" << std::endl;

  // build a DofMap holding dofs for all nodes including additional dofs of enriched nodes
  XFEM::createDofMapCombust(*interfacehandle, nodeDofMap_, elementDofMap_, fieldset, params);

  //------------------------------------------------------------------------------------------------
  // copy non-const dof maps to const dof maps
  // remark 1: the whole map should actually be const, not just the second element (a set)
  // remark 2: here, we need the general XFEM dof maps to be filled with the combustion dof maps,
  //           because member functions like "toGmsh" and "GatherUniqueEnrichments" rely on them
  //
  // what the following code does, is basically:   nodalDofSet_ = nodeDofMap_
  //                                             elementalDofs_ = elementDofMap_
  //------------------------------------------------------------------------------------------------
  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = nodeDofMap_.begin(); oneset != nodeDofMap_.end(); ++oneset )
  {
    nodalDofSet_.insert( make_pair(oneset->first, oneset->second));
  };
  for ( std::map<int, std::set<XFEM::FieldEnr> >::const_iterator oneset = elementDofMap_.begin(); oneset != elementDofMap_.end(); ++oneset )
  {
    elementalDofs_.insert( make_pair(oneset->first, oneset->second));
  };

  // collect all unique enrichments and print them on the screen (only needed for verification)
  GatherUniqueEnrichments();

  //if (ih_->xfemdis()->Comm().MyPID() == 0)
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
  for (map<int, const std::set<XFEM::FieldEnr> >::const_iterator fieldenriter=nodalDofSet_.begin();
  fieldenriter!=nodalDofSet_.end(); ++fieldenriter)
  {
    const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
      enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      unique_enrichments.insert(enrfield->getEnrichment());
    }
  }

  // collect enrichments from elemental dofs
  for (map<int, const std::set<XFEM::FieldEnr> >::const_iterator fieldenriter=elementalDofs_.begin();
  fieldenriter!=elementalDofs_.end(); ++fieldenriter)
  {
    const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
      enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      unique_enrichments.insert(enrfield->getEnrichment());
    }
  }

  // screen output
  if (ih_->xfemdis()->Comm().MyPID() == 0)
  {
    std::cout << " Enrichments available:";// << endl;
    for (std::set<XFEM::Enrichment>::const_iterator enr =
      unique_enrichments.begin(); enr != unique_enrichments.end(); ++enr)
    {
      std::cout << " " << enr->toString();// << endl;
    }
    std::cout << endl;
  }
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
  ih_->xfemdis()->Comm().SumAll(&locnumnodaldof,&numnodaldof,1);

  return numnodaldof;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::DofManager::fillDofDistributionMaps(
    std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofGID>&  NodalDofDistributionMap,
    std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofGID>&  ElementalDofDistributionMap
) const
{
  NodalDofDistributionMap.clear();
  // loop all (non-overlapping = Row)-Nodes and store the DOF information w.t.h. of DofKeys
  for (int i=0; i<ih_->xfemdis()->NumMyRowNodes(); ++i)
  {
    const DRT::Node* actnode = ih_->xfemdis()->lRowNode(i);
    const int gid = actnode->Id();
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
    if (entry == nodalDofSet_.end())
    {
      // no dofs for this node... must be a hole or somethin'
      continue;
    }
    const std::vector<int> gdofs(ih_->xfemdis()->Dof(actnode));
    const std::set<FieldEnr> dofset = entry->second;
    if (gdofs.size() != dofset.size())
    {
      cout << "numdof node (Discretization): " <<  gdofs.size() << endl;
      cout << "numdof node (DofManager):     " <<  dofset.size() << endl;
      dserror("Bug!!! Information about nodal dofs in DofManager and Discretization does not fit together!");
    }

    int dofcount = 0;
    std::set<FieldEnr>::const_iterator fieldenr;
    for(fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
    {
      NodalDofDistributionMap.insert(make_pair(DofKey<onNode>(gid, *fieldenr), gdofs[dofcount]));
      dofcount++;
    }
  };

  ElementalDofDistributionMap.clear();
  // loop all (non-overlapping = Row)-Elements and store the DOF information w.t.h. of DofKeys
  for (int i=0; i<ih_->xfemdis()->NumMyRowElements(); ++i)
  {
    const DRT::Element* actele = ih_->xfemdis()->lRowElement(i);
    const int gid = actele->Id();
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = elementalDofs_.find(gid);
    if (entry == elementalDofs_.end())
    {
      // no dofs for this node... must be a hole or somethin'
      continue;
    }
    const std::vector<int> gdofs(ih_->xfemdis()->Dof(actele));
    const std::set<FieldEnr> dofset = entry->second;
    if (gdofs.size() != dofset.size())
    {
//      cout << "numdof node (Discretization): " <<  gdofs.size() << endl;
//      cout << "numdof node (DofManager):     " <<  dofset.size() << endl;
//      dserror("Bug!!! Information about element dofs in DofManager and Discretization does not fit together!");
      // TODO: this mismatch is known and a better structure for element dofs should be found,
      //       e.g. DofKey(gid,fieldenr,dofperfieldenr), the latter being the number of dofs per fieldenr in the element
    }

    int dofcount = 0;
    std::set<FieldEnr>::const_iterator fieldenr;
    for(fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
    {
      ElementalDofDistributionMap.insert(make_pair(DofKey<onElem>(gid, *fieldenr), gdofs[dofcount]));
      dofcount++;
    }
  };
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::DofManager::fillPhysicalOutputVector(
    const Epetra_Vector&                    original_vector,
    const DRT::DofSet&                      dofset_out,
    const map<DofKey<onNode>, DofGID>&      nodalDofDistributionMap,
    const std::set<XFEM::PHYSICS::Field>&   fields_out
) const
{
  Teuchos::RCP<Epetra_Vector> outvec = LINALG::CreateVector(*dofset_out.DofRowMap(),true);

  const std::size_t numdof = fields_out.size();

  const Epetra_Map* dofrowmap = dofset_out.DofRowMap();
  const Epetra_Map* xdofrowmap = ih_->xfemdis()->DofRowMap();

  for (int i=0; i<ih_->xfemdis()->NumMyRowNodes(); ++i)
  {
    const DRT::Node* xfemnode = ih_->xfemdis()->lRowNode(i);
    const int gid = xfemnode->Id();
    const std::vector<int> gdofs(dofset_out.Dof(xfemnode));


    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
    if (entry == nodalDofSet_.end())
    {
      // no dofs for this node... must be a hole or somethin'
      //cout << "hole" << endl;
      for (std::size_t idof = 0; idof < numdof; ++idof)
      {
        //cout << dofrowmap->LID(gdofs[idof]) << endl;
        (*outvec)[dofrowmap->LID(gdofs[idof])] = 0.0;
      }
    }
    else
    {
      //cout << "some values available" << endl;

      //const std::vector<int> gdofs(ih_->xfemdis()->Dof(actnode));
      const std::set<FieldEnr> dofset = entry->second;

      const LINALG::Matrix<3,1> actpos(xfemnode->X());
      int idof = 0;
      for(std::set<XFEM::PHYSICS::Field>::const_iterator field_out = fields_out.begin(); field_out != fields_out.end(); ++field_out)
      {
        for(std::set<FieldEnr>::const_iterator fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
        {
          const XFEM::PHYSICS::Field fielditer = fieldenr->getField();
          if (fielditer == *field_out)
          {
            const double enrval = fieldenr->getEnrichment().EnrValue(actpos, *ih_, XFEM::Enrichment::approachUnknown);
            const XFEM::DofKey<XFEM::onNode> dofkey(gid,*fieldenr);
            const int origpos = nodalDofDistributionMap.find(dofkey)->second;
            //cout << origpos << endl;
            if (origpos < 0)
              dserror("bug!");
            if (gdofs[idof] < 0)
              dserror("bug!");
            (*outvec)[dofrowmap->LID(gdofs[idof])] += enrval * original_vector[xdofrowmap->LID(origpos)];
          }
        }
        //cout << "LID " << dofrowmap->LID(gdofs[idof]) << " -> GID " << gdofs[idof] << endl;
        idof++;
      }
    }
  };
  return outvec;
}


/*----------------------------------------------------------------------*
 | transform XFEM vector to (standard FEM) output vector     henke 07/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFEM::DofManager::transformXFEMtoOutputVector(
    const Epetra_Vector&                   xfemvector,
    const DRT::DofSet&                     outdofset,
    const map<DofKey<onNode>, DofGID>&     nodalDofDistributionMap,
    const std::set<XFEM::PHYSICS::Field>&  outputfields
) const
{
  // get DofRowMaps for output and XFEM vectors
  const Epetra_Map* outdofrowmap = outdofset.DofRowMap();
  const Epetra_Map* xfemdofrowmap = ih_->xfemdis()->DofRowMap();

  // create output vector (standard FEM layout)
  Teuchos::RCP<Epetra_Vector> outvector = LINALG::CreateVector(*outdofrowmap,true);

  // loop nodes on this processor
  for (int inode=0; inode<ih_->xfemdis()->NumMyRowNodes(); ++inode)
  {
    // get XFEM GID of this node
    const DRT::Node* xfemnode = ih_->xfemdis()->lRowNode(inode);
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
      const XFEM::Enrichment stdenr(0,XFEM::Enrichment::typeStandard);

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
          dserror("There should be a standard enrichment for every physical field!");
        }
        else // there is a standard enrichment for this desired output field
        {
          // build a dofkey (= XFEM dof)
          const XFEM::DofKey<XFEM::onNode> dofkey(nodegid,*fieldenrentry);
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
//            const XFEM::DofKey<XFEM::onNode> dofkey(nodegid,*fieldenr);
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
 |  transform  to a string (Debug)                              ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::DofManager::toString() const
{
  std::stringstream s;
  for (int i=0; i<ih_->xfemdis()->NumMyRowNodes(); ++i)
  {
    const int gid = ih_->xfemdis()->lRowNode(i)->Id();
    const set <XFEM::FieldEnr> actset = nodalDofSet_.find(gid)->second;
    for ( std::set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
    {
      s << "Node: " << gid << ", " << var->toString() << endl;
    };
  };
  return s.str();
}


/*----------------------------------------------------------------------*
 |  to Gmsh (Debug)                                             ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::toGmsh(
    const int step) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
#if 1

  const bool screen_out = false;

  const int myrank = ih_->xfemdis()->Comm().MyPID();

  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".numdof_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << DRT::Problem::Instance()->OutputControlFile()->FileName() << ".numdof_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"..."<<flush;
    std::ofstream f_system(filename.str().c_str());
    {
      // draw elements with associated gid
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Element->Id() \" {" << endl;
      for (int i=0; i<ih_->xfemdis()->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = ih_->xfemdis()->lColElement(i);
        gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(double(actele->Id()), actele);
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    {
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << " Stress unknowns in element \" {" << endl;
      for (int i=0; i<ih_->xfemdis()->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = ih_->xfemdis()->lColElement(i);
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator iter = elementalDofs_.find(actele->Id());

        if (iter != elementalDofs_.end())
        {
          const std::set<XFEM::FieldEnr> fieldenrset = iter->second;
          const double val = (double)fieldenrset.size();
          gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(val, actele);
        }
      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    {
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof per node \" {\n";
      for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
      {
        //DRT::Element* actele = ih_->xfemdis()->lColElement(i);
        const DRT::Node* xfemnode = ih_->xfemdis()->lColNode(i);
        const LINALG::Matrix<3,1> pos(xfemnode->X());

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());

        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fieldenrset = blub->second;

          gmshfilecontent << "SP(";
          gmshfilecontent << scientific << pos(0) << ",";
          gmshfilecontent << scientific << pos(1) << ",";
          gmshfilecontent << scientific << pos(2);
          gmshfilecontent << "){";
          gmshfilecontent << fieldenrset.size() << "};\n";
        }
      };
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof Jump enriched nodes \" {\n";
      for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->xfemdis()->lColNode(i);
        const LINALG::Matrix<3,1> pos(xfemnode->X());

        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeJump)
            {
              val = val+1.0;
            }
          }
          if (val > 0.5)
          {
            gmshfilecontent << "SP(";
            gmshfilecontent << scientific << pos(0) << ",";
            gmshfilecontent << scientific << pos(1) << ",";
            gmshfilecontent << scientific << pos(2);
            gmshfilecontent << "){";
            gmshfilecontent << val << "};\n";
          }
        }
      };
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof" << " standard enriched nodes \" {\n";
      for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->xfemdis()->lColNode(i);
        const LINALG::Matrix<3,1> pos(xfemnode->X());

        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeStandard)
            {
              val = val+1.0;
            }
          }
          if (val > 0.5)
          {
            gmshfilecontent << "SP(";
            gmshfilecontent << scientific << pos(0) << ",";
            gmshfilecontent << scientific << pos(1) << ",";
            gmshfilecontent << scientific << pos(2);
            gmshfilecontent << "){";
            gmshfilecontent << val << "};\n";
          }
        }
      };
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }

    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof" << " Void enriched nodes \" {\n";
      for (int i=0; i<ih_->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih_->xfemdis()->lColNode(i);
        const LINALG::Matrix<3,1> pos(xfemnode->X());

        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(xfemnode->Id());
        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeVoid)
            {
              val += 1.0;
            }
          }
          if (val > 0.5)
          {
            gmshfilecontent << "SP(";
            gmshfilecontent << scientific << pos(0) << ",";
            gmshfilecontent << scientific << pos(1) << ",";
            gmshfilecontent << scientific << pos(2);
            gmshfilecontent << "){";
            gmshfilecontent << val << "};\n";
          }
        }
      };
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }

    f_system.close();
    if (screen_out) std::cout << " done" << endl;
  }
#if 0
  if (gmshdebugout)
  {
      // debug info: print ele dofmanager information
      std::stringstream filename;
      std::stringstream filenamedel;
      const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
      filename    << filebase << "_eledofman_check_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
      filenamedel << filebase << "_eledofman_check_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
      std::remove(filenamedel.str().c_str());
      if (screen_out) std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"..."<<flush;
      std::ofstream f_system(filename.str().c_str());
      {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << " NumDofPerElement() in element \" {\n";
        for (int i=0; i<ih_->xfemdis()->NumMyColElements(); ++i)
        {
          DRT::Element* actele = ih_->xfemdis()->lColElement(i);
          //const int ele_gid = actele->Id();
          //double val = 0.0;
          //std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = elementalDofs_.find(ele_gid);
          const double val = actele->NumDofPerElement();
          gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(val, actele) << "\n";

        };
        gmshfilecontent << "};\n";
        f_system << gmshfilecontent.str();
      }
      f_system.close();
      if (screen_out) std::cout << " done" << endl;
  }
#endif
#endif
}



#endif  // #ifdef CCADISCRET
