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

#include "interface.H"
#include "dof_management_element.H"
#include "xfem_enums.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_discret.H"

namespace XFEM
{
  /*!
   * \brief pre-calculate enrichment values for an element and store it in this class
   */
  class ElementEnrichmentValues
  {
    
  public:
    
    /*!
    * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
    * 
    * if one already knows, an which side of the discontinuity one computes something,
    * the information can be given via approachdirection
    */
    explicit ElementEnrichmentValues(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>&     ih,               ///< interface information
        const XFEM::ElementDofManager&        dofman,
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
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field                 field,
        const VI&                                  funct,
        VO&                                        enr_funct
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const int numnode = ele_.NumNode();
      
      int dofcounter = 0;
      for (int inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }
    
    /*!
    * \brief For a given situation compute the enriched shape functions
    * 
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class MI, class VO, class MO>
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI& funct,
        const MI& derxy,
        VO& enr_funct,
        MO& enr_derxy
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const int numnode = ele_.NumNode();
      
      int dofcounter = 0;
      for (int inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            for (int isd = 0; isd < 3; ++isd)
            {
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
            }
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }
    
    /*!
    * \brief For a given situation compute the enriched shape functions
    * 
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class MI1, class MI2, class VO, class MO1, class MO2 >
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI& funct,
        const MI1& derxy,
        const MI2& derxy2,
        VO& enr_funct,
        MO1& enr_derxy,
        MO2& enr_derxy2
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const int numnode = ele_.NumNode();
      
      int dofcounter = 0;
      for (int inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
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
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }
    
    /*!
    * \brief For a given situation compute the enriched shape functions
    * 
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class VO>
    void ComputeEnrichedElementShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI&                  funct,
        VO&                        enr_funct
    ) const
    {
        int dofcounter = 0;
        
        const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
        const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
        const int numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
        dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          for (int inode = 0; inode < numvirtualnode; ++inode)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            dofcounter += 1;
          }
        }
        dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }
    
    /*!
    * \brief For a given situation compute the enriched shape functions
    * 
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class MI1, class VO, class MO1>
    void ComputeEnrichedElementShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI&                  funct,
        const MI1&                 derxy,
        VO&                        enr_funct,
        MO1&                       enr_derxy
    ) const
    {
        int dofcounter = 0;

        const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
        const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
        const int numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
        dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          for (int inode = 0; inode < numvirtualnode; ++inode)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            for (int isd = 0; isd < 3; ++isd)
            {
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
            }
            dofcounter += 1;
          }
        }
        dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }
    
  private:
    
    //! don't want default constructor
    explicit ElementEnrichmentValues();
    
    const DRT::Element&                        ele_;
    //  const RCP<XFEM::InterfaceHandle>           ih_;
    const XFEM::ElementDofManager&             dofman_;
    
    std::map<XFEM::Enrichment, double>        enrvals_;
    
  };


//! for output, compute values on integration cells
void computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&                 ele,
  const RCP<XFEM::InterfaceHandle>&   ih,
  const XFEM::ElementDofManager&      dofman,
  const GEO::DomainIntCell&           cell,
  const XFEM::PHYSICS::Field          field,
  const LINALG::SerialDenseVector&    elementvalues,
  LINALG::SerialDenseVector&          cellvalues
  );

//! for output, compute values on the integration cells
void computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&                   ele,
  const RCP<XFEM::InterfaceHandle>&     ih,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseVector&      elementvalues,
  LINALG::SerialDenseVector&            cellvalues
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
 * \brief calculates the fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 * 
 * \return ratio between fictitious volume and element volume (can be between 0.0 and 1.0)
 */
std::vector<double> DomainCoverageRatioPerNode(
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

/*! 
 * \brief calculates the volume fraction for each cell of an element .
 * 
 * \return ratio between cell volume and element volume (can be between 0.0 and 1.0)
 */
std::vector<double> DomainIntCellCoverageRatio(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
        );

//! return assembly type depending on the DOF distribution
XFEM::AssemblyType CheckForStandardEnrichmentsOnly(
    const ElementDofManager&   eleDofManager,
    const int                  numnode,
    const int*                 nodeids
);

/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
inline int getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        const int)
{
    return dofman.NumDofPerField(field);
}
template<> 
inline int getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager&,
        const XFEM::PHYSICS::Field,
        const int numnode)
{
    return numnode;
}

//! determine number of parameter
template<int NUMNODE, XFEM::AssemblyType ASSTYPE>
struct NumParam
{
  static int get(
     const XFEM::ElementDofManager& dofman,
     const XFEM::PHYSICS::Field     field
     )
  {
    return dofman.NumDofPerField(field);
  }
};

/// determine number of parameter for standard_assembly
template<int NUMNODE> struct NumParam<NUMNODE,XFEM::standard_assembly>
{
  static int get(
      const XFEM::ElementDofManager& ,
      const XFEM::PHYSICS::Field     
      )
  {
    return NUMNODE;
  }
};

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
