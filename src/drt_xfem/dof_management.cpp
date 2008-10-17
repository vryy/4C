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

#include "../drt_geometry/vector_definitions.H"
#include "dof_management.H"
#include "xdofmapcreation.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const RCP<XFEM::InterfaceHandle> ih, const bool DLM_condensation) :
  ih_(ih)
{
  XFEM::createDofMap(*ih, nodalDofSet_, elementalDofs_, DLM_condensation);
  
  std::set<XFEM::Enrichment> unique_enrichments = GatherUniqueEnrichments();

  if (ih_->xfemdis()->Comm().MyPID() == 0)
  {
    std::cout << " Enrichments available:" << endl;
    for (std::set<XFEM::Enrichment>::const_iterator enr =
      unique_enrichments.begin(); enr != unique_enrichments.end(); ++enr)
    {
      std::cout << "  - " << enr->toString() << endl;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<XFEM::Enrichment> XFEM::DofManager::GatherUniqueEnrichments() const
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
  return unique_enrichments;
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
 |  transform  to a string                                      ag 11/07|
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

void XFEM::DofManager::toGmsh(
    const Teuchos::RCP<XFEM::InterfaceHandle> ih,
    const int step
) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
#if 1
  
  const int myrank = ih->xfemdis()->Comm().MyPID();
  
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename    << allfiles.outputfile_kenner << "_numdof_coupled_system_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
    filenamedel << allfiles.outputfile_kenner << "_numdof_coupled_system_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"..."<<flush;
    std::ofstream f_system(filename.str().c_str());
    //f_system << IO::GMSH::disToString("Fluid", 0.0, ih->xfemdis(), ih->elementalDomainIntCells());
    //f_system << IO::GMSH::disToString("Solid", 1.0, ih->cutterdis(), *ih->cutterposnp());
    {
      // draw elements with associated gid
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Element->Id() \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
      {
        DRT::Element* actele = ih->xfemdis()->lColElement(i);
        gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(double(actele->Id()), actele);
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    {
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << " Stress unknowns in element \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
      {
        DRT::Element* actele = ih->xfemdis()->lColElement(i);
        const int ele_gid = actele->Id();
        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = elementalDofs_.find(ele_gid);
        
        if (blub != elementalDofs_.end())
        {
          const std::set<XFEM::FieldEnr> schnapp = blub->second;
          val = schnapp.size();
          gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(val, actele);
        }
        
      };
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    {
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof per node \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
      {
        //DRT::Element* actele = ih->xfemdis()->lColElement(i);
        const DRT::Node* xfemnode = ih->xfemdis()->lColNode(i);
        const BlitzVec3 pos(toBlitzArray(xfemnode->X()));
        const int node_gid = xfemnode->Id();

        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
        
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
      gmshfilecontent << "View \" " << "NumDof Jump enriched nodes \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih->xfemdis()->lColNode(i);
        const BlitzVec3 pos(toBlitzArray(xfemnode->X()));
        const int node_gid = xfemnode->Id();
        
        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
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
      gmshfilecontent << "View \" " << "NumDof" << " standard enriched nodes \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih->xfemdis()->lColNode(i);
        const BlitzVec3 pos(toBlitzArray(xfemnode->X()));
        const int node_gid = xfemnode->Id();
        
        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
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
            gmshfilecontent << val << "};" << endl;
          }
        }
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    
    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "NumDof" << " Void enriched nodes \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColNodes(); ++i)
      {
        const DRT::Node* xfemnode = ih->xfemdis()->lColNode(i);
        const BlitzVec3 pos(toBlitzArray(xfemnode->X()));
        const int node_gid = xfemnode->Id();
        
        double val = 0.0;
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
        if (blub != nodalDofSet_.end())
        {
          const std::set<XFEM::FieldEnr> fields = blub->second;
          for (std::set<XFEM::FieldEnr>::const_iterator f = fields.begin(); f != fields.end(); ++f)
          {
            if ((f->getEnrichment().Type()) == XFEM::Enrichment::typeVoid)
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
            gmshfilecontent << val << "};" << endl;
          }
        }
      };
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
      f_system.close();
      std::cout << " done" << endl;
    }
  }
  if (gmshdebugout)
  {
    {
      // debug info: print ele dofmanager information
      std::stringstream filename;
      std::stringstream filenamedel;
      filename    << allfiles.outputfile_kenner << "_eledofman_check_" << std::setw(5) << setfill('0') << step   << ".p" << myrank << ".pos";
      filenamedel << allfiles.outputfile_kenner << "_eledofman_check_" << std::setw(5) << setfill('0') << step-5 << ".p" << myrank << ".pos";
      std::remove(filenamedel.str().c_str());
      std::cout << "writing " << std::left << std::setw(50) <<filename.str()<<"..."<<flush;
      std::ofstream f_system(filename.str().c_str());
      {
        stringstream gmshfilecontent;
        gmshfilecontent << "View \" " << " NumDofPerElement() in element \" {" << endl;
        for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
        {
          DRT::Element* actele = ih->xfemdis()->lColElement(i);
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
      std::cout << " done" << endl;
    }
  }
#endif
}

///*----------------------------------------------------------------------*
// |  construct element dof manager                               ag 11/07|
// *----------------------------------------------------------------------*/
//XFEM::ElementDofManager XFEM::DofManager::constructElementDofManager(
//    const DRT::Element&  ele,
//    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz
//) const
//{
//  // nodal dofs for ele
//  std::map<int, const set <XFEM::FieldEnr> > nodaldofset;
//  for (int inode = 0; inode < ele.NumNode(); ++inode)
//  {
//    const int gid = ele.NodeIds()[inode];
//    nodaldofset.insert(make_pair(gid,this->getNodeDofSet(gid)));
//  }
//  
//  // element dofs for ele
//  std::set<XFEM::FieldEnr> enrfieldset;
//  
//  std::map<int,const std::set<XFEM::FieldEnr> >::const_iterator enrfieldsetiter = elementalDofs_.find(ele.Id());
//  if (enrfieldsetiter != elementalDofs_.end())
//  {
//    enrfieldset = enrfieldsetiter->second;
//  }
//  else
//  {
//    // use empty set
//  }
//  
//  // return a local dofmanager
//  return XFEM::ElementDofManager(ele, nodaldofset, enrfieldset, element_ansatz);
//}


///*----------------------------------------------------------------------*
// *----------------------------------------------------------------------*/
//void XFEM::DofManager::checkForConsistency(
//    const DRT::Element& ele,
//    const XFEM::ElementDofManager& stored_eledofman
//) const
//{
//  // create local copy of current information about dofs
//  const XFEM::ElementDofManager current_eledofman = this->constructElementDofManager(ele, stored_eledofman.getDisTypePerFieldMap());
//  
//  // compare with given and report error
//  if (current_eledofman != stored_eledofman)
//  {
//    dserror("given elementdofmanager is not consistent with global dofmanger");
//  }
//  return;
//}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::DofManager::fillDofDistributionMaps(
    NodalDofPosMap&      NodalDofDistributionMap,
    ElementalDofPosMap&  ElementalDofDistributionMap
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
    const Epetra_Vector&             original_vector,
    const DRT::DofSet&               dofset_out,
    const XFEM::NodalDofPosMap&      nodalDofDistributionMap,
    const std::set<XFEM::PHYSICS::Field>&   fields_out
) const
{
  Teuchos::RCP<Epetra_Vector> outvec = LINALG::CreateVector(*dofset_out.DofRowMap(),true);
  
  const int numdof = fields_out.size();
  
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
      for (int idof = 0; idof < numdof; ++idof)
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
      
      const BlitzVec3 actpos(toBlitzArray(xfemnode->X()));
      int idof = 0;
      for(std::set<XFEM::PHYSICS::Field>::const_iterator field_out = fields_out.begin(); field_out != fields_out.end(); ++field_out)
      {
        for(std::set<FieldEnr>::const_iterator fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
        {
          const XFEM::PHYSICS::Field fielditer = fieldenr->getField(); 
          if (fielditer == *field_out)
          {
            const XFEM::Enrichment enr = fieldenr->getEnrichment();
            const double enrval = enr.EnrValue(actpos, *ih_, XFEM::Enrichment::approachUnknown);
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


#endif  // #ifdef CCADISCRET
