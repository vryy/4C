/*!
\file dof_management.cpp

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "xfem.H"
#include "dof_management.H"
#include "xdofmapcreation.H"
#include "enrichment_utils.H"
#include "dofkey.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_systemmatrix.H"



/*----------------------------------------------------------------------*
 |  default ctor                                                ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager() :
    nodalDofSet_(),
    elementnodalDofSet_()
{
    return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
        const DRT::Element& ele,
        const map<int, const set<XFEM::FieldEnr> >& nodalDofSet,
        const map<int, const set<XFEM::FieldEnr> >& elementnodalDofSet,
        const int numvirtualnodes
        ) :
            nodalDofSet_(nodalDofSet),
            elementnodalDofSet_(elementnodalDofSet)
{
    // count number of dofs for each node
    map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
    for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
        const int gid = tmp->first;
        const set<XFEM::FieldEnr> enrfieldset = tmp->second;
        const int numdof = enrfieldset.size();
        //dsassert(numdof > 0, "sollte jetzt noch nicht sein");
        nodalNumDof_[gid] = numdof;
    }

    // set number of parameters per field to zero
    for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp) {
        const set<XFEM::FieldEnr> enrfieldset = tmp->second;
        for (set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] = 0;
            paramsLocalEntries_[field] = vector<int>();
        }
    }
    for (tmp = elementnodalDofSet.begin(); tmp != elementnodalDofSet.end(); ++tmp) {
        const set<XFEM::FieldEnr> enrfieldset = tmp->second;
        for (set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield) {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] = 0;
            paramsLocalEntries_[field] = vector<int>();
        }
    }


    unique_enrichments_.clear();
    // count number of parameters per field
    // define local position of unknown by looping first over nodes and then over its unknowns!
    int counter = 0;
    const DRT::Node*const* nodes = ele.Nodes();
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
        map<int, const set <XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
        if (entry == nodalDofSet_.end())
            dserror("impossible ;-)");
        const set<XFEM::FieldEnr> enrfieldset = entry->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] += 1;
            paramsLocalEntries_[field].push_back(counter);
            unique_enrichments_.insert(enrfield->getEnrichment());
            counter++;
        }
    }
    // loop now over element dofs
    // for that we have to loop over the "virtual element nodes" ansatz function
    for (int inode=0; inode<numvirtualnodes; ++inode)
    {
        map<int, const set <XFEM::FieldEnr> >::const_iterator entry = elementnodalDofSet_.find(inode);
        if (entry == elementnodalDofSet_.end())
            dserror("impossible ;-)");
        const set<XFEM::FieldEnr> enrfieldset = entry->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            const XFEM::PHYSICS::Field field = enrfield->getField();
            numParamsPerField_[field] += 1;
            paramsLocalEntries_[field].push_back(counter);
            unique_enrichments_.insert(enrfield->getEnrichment());
            counter++;
        }
    }

    return;
}

/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(const ElementDofManager& old)
{
    dserror("no copying");
    return;
}

/*----------------------------------------------------------------------*
 |  dtor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::~ElementDofManager()
{
    return;
}

/*----------------------------------------------------------------------*
 |  transform  to a string                                      ag 11/07|
 *----------------------------------------------------------------------*/
std::string XFEM::ElementDofManager::toString() const
{
    std::stringstream s;
    map<int, const set<XFEM::FieldEnr> >::const_iterator tmp;
    for (tmp = nodalDofSet_.begin(); tmp != nodalDofSet_.end(); ++tmp)
    {
        const int gid = tmp->first;
        const set <XFEM::FieldEnr> actset = tmp->second;
        for ( set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
        {
            s << "Node: " << gid << ", " << var->toString() << endl;
        };
    };
    return s.str();
}






/*----------------------------------------------------------------------*
 |  ctor                                                        ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const RCP<XFEM::InterfaceHandle> ih) :
            xfemdis_(ih->xfemdis())
{
    XFEM::createDofMap(*ih, nodalDofSet_, elementalDofs_);

    unique_enrichments_.clear();
    for (map<int, const set<XFEM::FieldEnr> >::const_iterator fieldenriter=nodalDofSet_.begin();
            fieldenriter!=nodalDofSet_.end(); ++fieldenriter)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            unique_enrichments_.insert(enrfield->getEnrichment());
        }
    }
    for (map<int, const set<XFEM::FieldEnr> >::const_iterator fieldenriter=elementalDofs_.begin();
            fieldenriter!=elementalDofs_.end(); ++fieldenriter)
    {
        const std::set<XFEM::FieldEnr> enrfieldset = fieldenriter->second;
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            unique_enrichments_.insert(enrfield->getEnrichment());
        }
    }
    std::cout << " Enrichments available:" << endl;
    for (std::set<XFEM::Enrichment>::const_iterator enr =
        unique_enrichments_.begin(); enr != unique_enrichments_.end(); ++enr)
    {
      std::cout << "  - " << enr->toString() << endl;
    }
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                   ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::DofManager::DofManager(const XFEM::DofManager& dofman)
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
    for (int i=0; i<xfemdis_->NumMyRowNodes(); ++i)
    {
        const int gid = xfemdis_->lRowNode(i)->Id();
        const set <XFEM::FieldEnr> actset = nodalDofSet_.find(gid)->second;
        for ( set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
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
  std::stringstream filename;
  filename << "numdof_coupled_system_" << std::setw(5) << setfill('0') << step << ".pos";
  std::ofstream f_system(filename.str().c_str());
  //f_system << IO::GMSH::disToString("Fluid", 0.0, ih->xfemdis(), ih->elementalDomainIntCells());
  f_system << IO::GMSH::disToString("Solid", 1.0, ih->cutterdis(), *ih->currentcutterpositions());
  {
      // draw elements with associated gid
      std::stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << "Element->Id() \" {" << endl;
      for (int i=0; i<ih->xfemdis()->NumMyColElements(); ++i)
      {
          DRT::Element* actele = ih->xfemdis()->lColElement(i);
          gmshfilecontent << IO::GMSH::elementToString(double(actele->Id()), actele);
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
              const set<XFEM::FieldEnr> schnapp = blub->second;
              val = schnapp.size();
              gmshfilecontent << IO::GMSH::elementToString(val, actele);
          }

      };
      gmshfilecontent << "};" << endl;
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

          double val = 0.0;
          std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);

          if (blub != nodalDofSet_.end())
          {
              const set<XFEM::FieldEnr> schnapp = blub->second;
              val = schnapp.size();


          gmshfilecontent << "SP(";
          gmshfilecontent << scientific << pos(0) << ",";
          gmshfilecontent << scientific << pos(1) << ",";
          gmshfilecontent << scientific << pos(2);
          gmshfilecontent << "){";
          gmshfilecontent << val << "};" << endl;
          }
      };
      gmshfilecontent << "};" << endl;
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
          std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
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
                  gmshfilecontent << val << "};" << endl;
              }
          }
      };
      gmshfilecontent << "};" << endl;
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
          std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
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
          std::map<int, const set<XFEM::FieldEnr> >::const_iterator blub = nodalDofSet_.find(node_gid);
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
  }


  //f_system << IO::GMSH::getConfigString(2);
  f_system.close();
}

/*----------------------------------------------------------------------*
 |  construct element dof manager                               ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager XFEM::DofManager::constructElementDofManager(
        const DRT::Element&  ele,
        const int      numeleparam
        ) const
{
    // create a list with number of dofs per local node
    const int numnode = ele.NumNode();
    const int* nodegids = ele.NodeIds();

    // nodal dofs for ele
    std::map<int, const set <XFEM::FieldEnr> > nodaldofset;
    for (int inode = 0; inode < numnode; ++inode) {
        const int gid = nodegids[inode];
        nodaldofset.insert(make_pair(gid,this->getNodeDofSet(gid)));
    }

    // element dofs for ele
    std::map<int, const set <XFEM::FieldEnr> > elementnodaldofset;
    for (int inode = 0; inode < numeleparam; ++inode) {
        elementnodaldofset.insert(make_pair(inode,this->getElementDofSet(ele.Id())));
    }

    // return a local dofmanager
    return XFEM::ElementDofManager(ele, nodaldofset, elementnodaldofset, numeleparam);
}

/*----------------------------------------------------------------------*
 |  sanity check                                                ag 11/07|
 *----------------------------------------------------------------------*/
void XFEM::DofManager::checkForConsistency(
        const DRT::Element& ele,
        const XFEM::ElementDofManager& stored_eledofman
        ) const
{
    // create local copy of current information about dofs
    XFEM::ElementDofManager current_eledofman = this->constructElementDofManager(ele, stored_eledofman.NumVirtualNodes());

    // compare with given and report error
    if (current_eledofman != stored_eledofman)
    {
        dserror("given elementdofmanager is not consistent with global dofmanger");
    }
    return;
}


void XFEM::DofManager::fillDofDistributionMaps(
    NodalDofPosMap&      NodalDofDistributionMap,
    ElementalDofPosMap&  ElementalDofDistributionMap
        ) const
{
    NodalDofDistributionMap.clear();
    // loop all (non-overlapping = Row)-Nodes and store the DOF information w.t.h. of DofKeys
    for (int i=0; i<xfemdis_->NumMyRowNodes(); ++i)
    {
        const DRT::Node* actnode = xfemdis_->lRowNode(i);
        const int gid = actnode->Id();
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet_.find(gid);
        if (entry == nodalDofSet_.end())
        {
            // no dofs for this node... must be a hole or somethin'
            continue;
        }
        const std::vector<int> gdofs(xfemdis_->Dof(actnode));
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
    for (int i=0; i<xfemdis_->NumMyRowElements(); ++i)
    {
        const DRT::Element* actele = xfemdis_->lRowElement(i);
        const int gid = actele->Id();
        std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = elementalDofs_.find(gid);
        if (entry == elementalDofs_.end())
        {
            // no dofs for this node... must be a hole or somethin'
            continue;
        }
        const std::vector<int> gdofs(xfemdis_->Dof(actele));
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


XFEM::AssemblyType XFEM::CheckForStandardEnrichmentsOnly(
        const ElementDofManager&   eleDofManager_,
        const int                  numnode,
        const int*                 nodeids
        )
{
    // find out whether we can use standard assembly or need xfem assembly
    XFEM::AssemblyType assembly_type = XFEM::standard_assembly;
    for (int inode = 0; inode < numnode; ++inode)
    {
        if (assembly_type == XFEM::xfem_assembly)
        {
            break;
        }
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& fields = eleDofManager_.FieldEnrSetPerNode(gid);
        if (fields.size() != 4)
        {
            assembly_type = XFEM::xfem_assembly;
            break;
        };
        for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fields.begin(); fieldenr != fields.end(); ++fieldenr)
        {
            if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
            {
                assembly_type = XFEM::xfem_assembly;
                break;
            };
        };
    };
    const int eledof = eleDofManager_.NumDofPerElement();
    if (eledof != 0)
    {
        assembly_type = XFEM::xfem_assembly;
    }

    return assembly_type;
}


#endif  // #ifdef CCADISCRET
