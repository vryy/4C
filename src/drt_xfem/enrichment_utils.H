/*!
\file enrichment_utils.H

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef ENRICHMENT_UTILS_H
#define ENRICHMENT_UTILS_H

#include "dof_management_element.H"
#include "xfem_enums.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"


class DRT::Element;
class DRT::Discretization;

namespace XFEM
{

  class InterfaceHandle;

  /*!
   * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
   * 
   * if one already knows, an which side of the discontinuity one computes something,
   * the information can be given via approachdirection
   */
std::map<XFEM::Enrichment, double> computeEnrvalMap(
        const RCP<XFEM::InterfaceHandle>  ih,               ///< interface onformation
        const std::set<XFEM::Enrichment>&     enrset,           ///< enrichments for which one seeks enrichment values
        const LINALG::Matrix<3,1>&            actpos,           ///< a 3d position
        const XFEM::Enrichment::ApproachFrom  approachdirection ///< is location with respect to interface is known or not
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
template<class VI, class VO>
void ComputeEnrichedNodalShapefunction(
    const DRT::Element&                        ele,
    const RCP<XFEM::InterfaceHandle>           ih,
    const XFEM::ElementDofManager&             dofman,
    const XFEM::PHYSICS::Field                 field,
    const std::map<XFEM::Enrichment, double>&  enrvals,
    const VI&                                  funct,
    VO&                                        enr_funct
    )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
template<class VI, class MI, class VO, class MO>
void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const VI& funct,
        const MI& derxy,
        VO& enr_funct,
        MO& enr_derxy
        )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                for (int isd = 0; isd < 3; ++isd)
                {
                  enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
                }
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
template<class VI, class MI1, class MI2, class VO, class MO1, class MO2 >
void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const VI& funct,
        const MI1& derxy,
        const MI2& derxy2,
        VO& enr_funct,
        MO1& enr_derxy,
        MO2& enr_derxy2
        )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                for (int isd = 0; isd < 3; ++isd)
                {
                  enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
                }
                for (int isd = 0; isd < 6; ++isd)
                {
                  enr_derxy2(isd,dofcounter) = derxy2(isd,inode) * enrval;
                }
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedElementShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const LINALG::SerialDenseVector& funct,
        LINALG::SerialDenseVector& enr_funct
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedElementShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const LINALG::SerialDenseVector& funct,
        const LINALG::SerialDenseMatrix& derxy,
        LINALG::SerialDenseVector& enr_funct,
        LINALG::SerialDenseMatrix& enr_derxy
        );

//! for output, compute values on integration cells
void computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseVector& elementvalues,
  LINALG::SerialDenseVector&     cellvalues
  );

//! for output, compute values on the integration cells
void computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseVector& elementvalues,
  LINALG::SerialDenseVector&     cellvalues
  );

//! for output, compute values on the integration cells
void computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

/*! 
 * \brief calculates the fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 * 
 * \return ration between fictitious volume and element volume (can be between 0.0 and 1.0)
 */
double DomainCoverageRatio(
    const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

/*! 
 * \brief calculates the area fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 */
double BoundaryCoverageRatio(
    const DRT::Element&               ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

//! return assembly type depending on the DOF distribution
XFEM::AssemblyType CheckForStandardEnrichmentsOnly(
    const ElementDofManager&   eleDofManager,
    const int                  numnode,
    const int*                 nodeids
);

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
