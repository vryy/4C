/*!
\file dof_management_element.cpp

\brief provides the element dofmanager class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET

#include "dof_management_element.H"

#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager() :
  numElemDof_(0),
  nodalDofSet_()
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
    const DRT::Element& ele,
    const map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
) :
  nodalDofSet_(nodalDofSet),
  DisTypePerElementField_(element_ansatz)
{
  ComputeDependendInfo(ele, nodalDofSet, enrfieldset, element_ansatz);
 
  return;
}
    
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::ElementDofManager::ComputeDependendInfo(
    const DRT::Element& ele,
    const map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
)
{
  // count number of dofs for each node
  map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp;
  for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp)
  {
    const int gid = tmp->first;
    nodalNumDof_[gid] = tmp->second.size();
  }
  
  // set number of parameters per field to zero
  for (tmp = nodalDofSet.begin(); tmp != nodalDofSet.end(); ++tmp)
  {
    const std::set<XFEM::FieldEnr> enrfieldset = tmp->second;
    for (set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      numParamsPerField_[field] = 0;
      paramsLocalEntries_[field] = vector<int>();
    }
  }
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    numParamsPerField_[field] = 0;
    paramsLocalEntries_[field] = vector<int>();
  }
      
      
  unique_enrichments_.clear();
  // count number of parameters per field
  // define local position of unknown by looping first over nodes and then over its unknowns!
  int dofcounter = 0;
  const int* nodeids = ele.NodeIds();
  for (int inode=0; inode<ele.NumNode(); ++inode)
  {
    const int gid = nodeids[inode];
    map<int, const set <XFEM::FieldEnr> >::const_iterator entry = nodalDofSet.find(gid);
    if (entry == nodalDofSet_.end())
      dserror("impossible ;-)");
    const std::set<XFEM::FieldEnr> enrfieldset = entry->second;
    
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      numParamsPerField_[field] += 1;
      paramsLocalEntries_[field].push_back(dofcounter);
      unique_enrichments_.insert(enrfield->getEnrichment());
      dofcounter++;
    }
  }
  numNodeDof_ = dofcounter;
  
  // loop now over element dofs
  // we first loop over the fields and then over the params
  numElemDof_ = 0;
  enrichedFieldperPhysField_.clear();
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    //cout << physVarToString(field) << endl;
    std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator schnack = element_ansatz.find(field);
    if (schnack == element_ansatz.end())
    {
      cout << XFEM::PHYSICS::physVarToString(field) << endl;
      dserror("field not found -> bug");
    }
    
    enrichedFieldperPhysField_[field].insert(*enrfield);
    
    const DRT::Element::DiscretizationType eledofdistype = schnack->second;
    const int numparam = DRT::UTILS::getNumberOfElementNodes(eledofdistype);
    for (int inode=0; inode<numparam; ++inode)
    {
      numElemDof_ +=1;
      numParamsPerField_[field] += 1;
      paramsLocalEntries_[field].push_back(dofcounter);
      unique_enrichments_.insert(enrfield->getEnrichment());
      dofcounter++;
    }
  }
  
  for (std::map<XFEM::PHYSICS::Field, std::vector<int> >::const_iterator it = paramsLocalEntries_.begin(); 
       it != paramsLocalEntries_.end();
       ++it)
  {
    const XFEM::PHYSICS::Field field = it->first;
    std::vector<int> dofpositions = it->second;
    if (field == XFEM::PHYSICS::Velx)
    {
      paramsLocalEntriesVelx_.resize(dofpositions.size());
      std::copy (dofpositions.begin(), dofpositions.begin() + dofpositions.size(), paramsLocalEntriesVelx_.begin());
    }
    else if (field == XFEM::PHYSICS::Vely)
    {
      paramsLocalEntriesVely_.resize(dofpositions.size());
      std::copy (dofpositions.begin(), dofpositions.begin() + dofpositions.size(), paramsLocalEntriesVely_.begin());
    }
    else if (field == XFEM::PHYSICS::Velz)
    {
      paramsLocalEntriesVelz_.resize(dofpositions.size());
      std::copy (dofpositions.begin(), dofpositions.begin() + dofpositions.size(), paramsLocalEntriesVelz_.begin());
    }
    else if (field == XFEM::PHYSICS::Pres)
    {
      paramsLocalEntriesPres_.resize(dofpositions.size());
      std::copy (dofpositions.begin(), dofpositions.begin() + dofpositions.size(), paramsLocalEntriesPres_.begin());
    }
  }
  
  if (dofcounter != (numNodeDof_ + numElemDof_))
    dserror("dof number mismatch! -> bug!");
  
  return;
}
    
/*----------------------------------------------------------------------*
 |  construct element dof manager                               ag 11/07|
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
    const DRT::Element&  ele,
    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>& element_ansatz,
    const XFEM::DofManager& dofman
) :
  DisTypePerElementField_(element_ansatz)
{
  // nodal dofs for ele
  nodalDofSet_.clear();
  for (int inode = 0; inode < ele.NumNode(); ++inode)
  {
    const int gid = ele.NodeIds()[inode];
    nodalDofSet_.insert(make_pair(gid,dofman.getNodeDofSet(gid)));
  }
  
  // element dofs for ele
  const std::set<XFEM::FieldEnr>& enrfieldset(dofman.getElementDofSet(ele.Id()));
  
  ComputeDependendInfo(ele, nodalDofSet_, enrfieldset, element_ansatz);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string XFEM::ElementDofManager::toString() const
{
  std::stringstream s;
  map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp;
  for (tmp = nodalDofSet_.begin(); tmp != nodalDofSet_.end(); ++tmp)
  {
    const int gid = tmp->first;
    const set <XFEM::FieldEnr> actset = tmp->second;
    for ( std::set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
    {
      s << "Node: " << gid << ", " << var->toString() << endl;
    };
  };
  return s.str();
}



//! return reference to list of local positions in a array of dofs
template <>
const std::vector<int>& XFEM::ElementDofManager::LocalDofPosPerField<XFEM::PHYSICS::Velx>() const
{
  return paramsLocalEntriesVelx_;
}

//! return reference to list of local positions in a array of dofs
template <>
const std::vector<int>& XFEM::ElementDofManager::LocalDofPosPerField<XFEM::PHYSICS::Vely>() const
{
  return paramsLocalEntriesVely_;
}

//! return reference to list of local positions in a array of dofs
template <>
const std::vector<int>& XFEM::ElementDofManager::LocalDofPosPerField<XFEM::PHYSICS::Velz>() const
{
  return paramsLocalEntriesVelz_;
}

//! return reference to list of local positions in a array of dofs
template <>
const std::vector<int>& XFEM::ElementDofManager::LocalDofPosPerField<XFEM::PHYSICS::Pres>() const
{
  return paramsLocalEntriesPres_;
}

    
#endif  // #ifdef CCADISCRET
