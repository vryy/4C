/*!
\file dof_management_element.cpp

\brief provides the element dofmanager class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 */


#include "dof_management.H"
#include "dof_management_element.H"
#include "field_enriched.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_combust/combust_defines.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager() :
  numNodeDof_(0),
  numElemDof_(0)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::ElementDofManager::ElementDofManager(
    const DRT::Element& ele,
    const std::map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,
    const std::set<XFEM::FieldEnr>& enrfieldset,
    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz
) :
  nodalDofSet_(nodalDofSet),
  DisTypePerElementField_(element_ansatz)
{
  ComputeDependentInfo(ele, nodalDofSet_, enrfieldset, element_ansatz);

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
  std::size_t numnode = ele.NumNode();
  for (std::size_t inode = 0; inode < numnode; ++inode)
  {
    const int gid = ele.NodeIds()[inode];
    nodalDofSet_.insert(std::make_pair(gid,dofman.getNodeDofSet(gid)));
  }

  // element dofs for ele
  const std::set<XFEM::FieldEnr>& enrfieldset(dofman.getElementDofSet(ele.Id()));

  ComputeDependentInfo(ele, nodalDofSet_, enrfieldset, element_ansatz);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::ElementDofManager::ComputeDependentInfo(
    const DRT::Element& ele,
    const std::map<int, const std::set<XFEM::FieldEnr> >& nodalDofSet,  ///< node dofs
    const std::set<XFEM::FieldEnr>& enrfieldset,                   ///< element dofs
    const std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz)
{
  // count number of dofs for each node
//   for (std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet.begin();
//        tmp != nodalDofSet.end();
//        ++tmp)
//   {
//     const int nodegid = tmp->first;
//     //nodalNumDof_[nodegid] = tmp->second.size();
//   }

  // set number of parameters per field to zero
  // for nodal dofs
  for (std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet.begin();
       tmp != nodalDofSet.end();
       ++tmp)
  {
    const std::set<XFEM::FieldEnr> lenrfieldset = tmp->second;
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = lenrfieldset.begin();
         enrfield != lenrfieldset.end();
         ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      //numParamsPerField_[field] = 0;
      paramsLocalEntries_[field] = std::vector<int>();
    }
  }
  // for element dofs
  for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin();
       enrfield != enrfieldset.end();
       ++enrfield)
  {
    const XFEM::PHYSICS::Field field = enrfield->getField();
    //numParamsPerField_[field] = 0;
    paramsLocalEntries_[field] = std::vector<int>();
  }

  unique_enrichments_.clear();

  // define local position of unknown by looping first over nodes and then over its unknowns!
  std::size_t dofcounter = 0;
  std::size_t numnode = ele.NumNode();
  const int* nodeids = ele.NodeIds();
  for (std::size_t inode=0; inode<numnode; ++inode)
  {
    const int nodegid = nodeids[inode];
    std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet.find(nodegid);
    if (entry == nodalDofSet.end())
      dserror("impossible");
    const std::set<XFEM::FieldEnr> & lenrfieldset = entry->second;

    for (std::set<XFEM::FieldEnr>::const_iterator enrfield = lenrfieldset.begin();
         enrfield != lenrfieldset.end();
         ++enrfield)
    {
      const XFEM::PHYSICS::Field field = enrfield->getField();
      paramsLocalEntries_[field].push_back(dofcounter);
      lidtofieldenr_.push_back(*enrfield);
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
    //std::cout << physVarToString(field) << std::endl;
    std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator schnack = element_ansatz.find(field);
    if (schnack == element_ansatz.end())
    {
      std::cout << XFEM::PHYSICS::physVarToString(field) << std::endl;
      dserror("field not found -> bug");
    }

    enrichedFieldperPhysField_[field].insert(*enrfield);

    std::vector<int> & local_entries = paramsLocalEntries_[field];

    const DRT::Element::DiscretizationType eledofdistype = schnack->second;
    const int numparam = DRT::UTILS::getNumberOfElementNodes(eledofdistype);
    for (int inode=0; inode<numparam; ++inode)
    {
      numElemDof_ +=1;
      //numParamsPerField_[field] += 1;
      local_entries.push_back(dofcounter);
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
  std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp;
  for (tmp = nodalDofSet_.begin(); tmp != nodalDofSet_.end(); ++tmp)
  {
    const int gid = tmp->first;
    const std::set <XFEM::FieldEnr> actset = tmp->second;
    for ( std::set<XFEM::FieldEnr>::const_iterator var = actset.begin(); var != actset.end(); ++var )
    {
      s << "Node: " << gid << ", " << var->toString() << std::endl;
    }
  }
  return s.str();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::set<XFEM::FieldEnr>& XFEM::ElementDofManager::FieldEnrSetPerNode(
    const int  gid
    ) const
{
  std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
  if (tmp == nodalDofSet_.end())
  {
    std::cout << gid << std::endl;
    // let this always be here so we can recognize errors early
    dserror("FieldEnrSetPerNode: global node id not found!");
    //return std::set<FieldEnr>();
  }
  return tmp->second;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const XFEM::FieldEnr& XFEM::ElementDofManager::FieldEnrSetPerDof(
    const int  lid
    ) const
{
  return lidtofieldenr_[lid];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::size_t XFEM::ElementDofManager::NumDofPerNode(
    const int gid             ///< unique global node id
) const
{
  std::map<int, const std::set<XFEM::FieldEnr> >::const_iterator tmp = nodalDofSet_.find(gid);
  //map<int,std::size_t>::const_iterator tmp = nodalNumDof_.find(gid);
  dsassert(tmp != nodalDofSet_.end(), "node not found");
  return tmp->second.size();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::size_t XFEM::ElementDofManager::NumDofPerField(
    const XFEM::PHYSICS::Field  field  ///< field for which we seek the number of DOFs
) const
{
  std::map<XFEM::PHYSICS::Field, std::vector<int> >::const_iterator tmp = paramsLocalEntries_.find(field);
  //map<XFEM::PHYSICS::Field, std::size_t>::const_iterator tmp = numParamsPerField_.find(field);
  if (tmp == paramsLocalEntries_.end()){
    //std::cout << XFEM::PHYSICS::physVarToString(field) << std::endl;
    return 0;
  }
  return tmp->second.size();
}


//! function that returns dofs for each node of this particular element
const std::map<int, const std::set<XFEM::FieldEnr> >& XFEM::ElementDofManager::getNodalDofSet(
) const
{
  return nodalDofSet_;
}


//! return reference to list of local positions in a array of dofs
const std::vector<int>& XFEM::ElementDofManager::LocalDofPosPerField(
    const XFEM::PHYSICS::Field field
) const
{
  std::map<XFEM::PHYSICS::Field, std::vector<int> >::const_iterator tmp = paramsLocalEntries_.find(field);
  if (tmp == paramsLocalEntries_.end()){
    dserror( "no local dofs for field %d", field );
  }
  return tmp->second;
}


//! return unique enrichments for this proc
const std::set<XFEM::Enrichment>& XFEM::ElementDofManager::getUniqueEnrichments() const
{
  return unique_enrichments_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> XFEM::ElementDofManager::getUniqueEnrichmentLabels(
) const
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


//! access to map
DRT::Element::DiscretizationType XFEM::ElementDofManager::getDisTypePerField(
    XFEM::PHYSICS::Field field
) const
{
  std::map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType>::const_iterator i = DisTypePerElementField_.find(field);
  if ( i==DisTypePerElementField_.end() )
  {
    dserror( "no field %d", field );
  }
  return i->second;
}


//! access to map
const std::set<XFEM::FieldEnr>& XFEM::ElementDofManager::getEnrichedFieldsPerEleField(
    XFEM::PHYSICS::Field field
) const
{
  std::map<XFEM::PHYSICS::Field, std::set<XFEM::FieldEnr> >::const_iterator i = enrichedFieldperPhysField_.find(field);
  if ( i==enrichedFieldperPhysField_.end() )
  {
    dserror( "no field %d", field );
  }
  return i->second;
}


//! not equal operator (for use in STL iterators)
bool XFEM::ElementDofManager::operator !=(const ElementDofManager& rhs) const
{
  // in principle, they can be only not equal, if the nodalDofSet is not equal
  // everything else is derived from that in the constructor
  return nodalDofSet_ != rhs.nodalDofSet_;
}


