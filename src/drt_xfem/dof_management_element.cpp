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
  nodalDofSet_(),
  numElemDof_(0)
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
  ComputeDependentInfo(ele, nodalDofSet, enrfieldset, element_ansatz);

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
  for (std::size_t inode = 0; inode < (size_t)ele.NumNode(); ++inode)
  {
    const int gid = ele.NodeIds()[inode];
    nodalDofSet_.insert(make_pair(gid,dofman.getNodeDofSet(gid)));
  }

  // element dofs for ele
  const std::set<XFEM::FieldEnr>& enrfieldset(dofman.getElementDofSet(ele.Id()));

  ComputeDependentInfo(ele, nodalDofSet_, enrfieldset, element_ansatz);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::ElementDofManager::ComputeDependentInfo(
    const DRT::Element& ele,
    const map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz)
{
  // count number of dofs for each node
  for (map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet.begin();
       tmp != nodalDofSet.end();
       ++tmp)
  {
    const int gid = tmp->first;
    nodalNumDof_[gid] = tmp->second.size();
  }

  // set number of parameters per field to zero
  for (map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet.begin();
       tmp != nodalDofSet.end();
       ++tmp)
  {
    const std::set<XFEM::FieldEnr> lenrfieldset = tmp->second;
    for (set<XFEM::FieldEnr>::const_iterator enrfield = lenrfieldset.begin();
         enrfield != lenrfieldset.end();
         ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      numParamsPerField_[field] = 0;
      paramsLocalEntries_[field] = vector<int>();
    }
  }
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin();
       enrfield != enrfieldset.end();
       ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    numParamsPerField_[field] = 0;
    paramsLocalEntries_[field] = vector<int>();
  }


  unique_enrichments_.clear();
  // count number of parameters per field
  // define local position of unknown by looping first over nodes and then over its unknowns!
  std::size_t dofcounter = 0;
  const int* nodeids = ele.NodeIds();
  for (std::size_t inode=0; inode<(size_t)ele.NumNode(); ++inode)
  {
    const int gid = nodeids[inode];
    map<int, const set <XFEM::FieldEnr> >::const_iterator entry = nodalDofSet.find(gid);
    if (entry == nodalDofSet.end())
      dserror("impossible ;-)");
    const std::set<XFEM::FieldEnr> lenrfieldset = entry->second;

    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = lenrfieldset.begin();
         enrfield != lenrfieldset.end();
         ++enrfield)
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

  if (dofcounter != (numNodeDof_ + numElemDof_))
    dserror("dof number mismatch! -> bug!");

  return;
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::set<XFEM::FieldEnr>& XFEM::ElementDofManager::FieldEnrSetPerNode(
    const int  gid
    ) const
{
  map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
  if (tmp == nodalDofSet_.end())
  {
    std::cout << gid << endl;
    // let this always be here so we can recognize errors early
    dserror("FieldEnrSetPerNode: global node id not found!");
    //return std::set<FieldEnr>();
  }
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::size_t XFEM::ElementDofManager::NumDofPerNode(
    const int gid             ///< unique global node id
) const
{
  map<int,std::size_t>::const_iterator tmp = nodalNumDof_.find(gid);
  dsassert(tmp != nodalNumDof_.end(), "node not found");
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::size_t XFEM::ElementDofManager::NumDofPerField(
    const XFEM::PHYSICS::Field  field  ///< field for which we seek the number of DOFs
) const
{
  map<XFEM::PHYSICS::Field, std::size_t>::const_iterator tmp = numParamsPerField_.find(field);
  if (tmp == numParamsPerField_.end()){
    //cout << XFEM::PHYSICS::physVarToString(field) << endl;
    return 0;
  }
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::ElementDofManager::getUniqueEnrichmentLabels() const
{
  std::set<int> xlabelset;
  for(std::set<XFEM::Enrichment>::const_iterator enr = unique_enrichments_.begin();
  enr != unique_enrichments_.end();
  ++enr)
  {
    xlabelset.insert(enr->XFEMConditionLabel());
  }
  return xlabelset;
}

#endif  // #ifdef CCADISCRET
