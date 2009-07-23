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
#include "../drt_combust/combust_interface.H"
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
     * \brief compute enrichment value for integration cell
     *
     * constructor called for combustion problems (jump enrichment)
     *
     * \author henke 06/09
     */
    template<class V1>
    explicit ElementEnrichmentValues(
            const DRT::Element&             ele,    ///< pointer to element
            const XFEM::ElementDofManager&  dofman, ///< pointer to element dof manager
            const GEO::DomainIntCell&       cell,   ///< pointer to integration cell
            const V1&                       phi     ///< vector of level set values
     ) :
       ele_(ele),
       dofman_(dofman)
    {
      enrvals_.clear();
      enrvalnodes_.clear();
      enrvalderxy_.clear();  // not used for jump enrichment
      enrvalderxy2_.clear(); // not used for jump enrichment

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();
      std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes;

      //------------------------------------------------------
      // evaluate enrichment values for this integration cell
      //------------------------------------------------------
      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const double enrvalcell = enriter->EnrValueIntCell(cell);
        enrvals_[*enriter] = enrvalcell;
      }

      //------------------------------------------------------------
      // loop over nodes to evaluate enrichment values at the nodes
      //------------------------------------------------------------
      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          const double enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
          enrvalstmp[*enriter] = enrvalnode;
        }
        enrvalnodes[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      //------------------------------------------------------------------------
      // evaluate final enrichment values for each node H(G(cell)) - H(G(node))
      //------------------------------------------------------------------------
      // literature (p. 1006, penultimate line):
      // Belytschko, T., Moës, N., Usui, S. and Parimi, C.
      // Arbitrary discontinuities in finite elements:
      // "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.

      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          enrvalstmp[*enriter] = enrvals_.find(*enriter)->second - (enrvalnodes.find(inode)->second).find(*enriter)->second;
        }
        enrvalnodes_[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      return;
    }


    /*!
    * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
    *
    * constructor called for two-phase flow problems (kink enrichment)
    *
    * \author rasthofer 06/09
    */
    template<class V1, class M1, class M2>
    explicit ElementEnrichmentValues(
            const DRT::Element&             ele,
            const XFEM::ElementDofManager&  dofman,
            const V1&                       phi,
            const V1&                       funct,   ///< vector of shape functions evaluated at Gauss point
            const M1&                       derxy,
            const M2&                       derxy2
            ) :
              ele_(ele),
              dofman_(dofman)
    {
      enrvals_.clear(); // not used for kink enrichment
      enrvalnodes_.clear();
      enrvalderxy_.clear();
      enrvalderxy2_.clear();

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());

      const int numnode = ele_.NumNode();

      for(int inode=0; inode<numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalues;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
//          const double enrval = enriter->ModifiedKinkEnrValue(actpos, inode, *ih, approachdirection, ele_);
          const double enrval = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct);
          double enrvalnode = 0.0;
          // modified enrichment in the case of jump enrichment to maintain the interpolant properties
          // Psi(xGauss)-Psi(xNode)
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
          }
          enrvalues[*enriter] = enrval - enrvalnode;
        }
        enrvalnodes_[inode] = enrvalues;
        enrvalues.clear();
      }

      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const std::vector<double> enrvalderxy = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
        enrvalderxy_[*enriter] = enrvalderxy;

        const std::vector<double> enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      return;
    }


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
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
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
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
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
            for (std::size_t isd = 0; isd < 3; ++isd)
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
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
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
            for (std::size_t isd = 0; isd < 3; ++isd)
            {
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
            }
            for (std::size_t isd = 0; isd < 6; ++isd)
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
     * \brief For a given situation compute the enriched shape functions for kink enrichment
     *
     * due to the kink enrichment, the chain rule is applied here
     *
     * \author henke (rasthofer) 06/09
     */
    template<class VI, class MI1, class MI2, class VO, class MO1, class MO2 >
    inline void ComputeKinkEnrichedNodalShapefunction(
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
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvalues.find(enrfield->getEnrichment())->second;
            const std::vector<double> enrvalderxy = enrvalderxy_.find(enrfield->getEnrichment())->second;
            const std::vector<double> enrvalderxy2 = enrvalderxy2_.find(enrfield->getEnrichment())->second;

            enr_funct(dofcounter) = funct(inode) * enrval;

            for (std::size_t isd = 0; isd < 3; ++isd)
            {
              //product rule
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval + funct(inode) * enrvalderxy[isd];
            }

//            for (std::size_t isd = 0; isd < 6; ++isd)
//            {
//              //product rule
//              // da fehlt noch was
//              enr_derxy2(isd,dofcounter) = derxy2(isd,inode) * enrval + funct(inode) * enrvalderxy2[isd];
//            }
            enr_derxy2(0,dofcounter) = derxy2(0,inode) * enrval + funct(inode) * enrvalderxy2[0] + derxy(0,inode) * enrvalderxy[0] + derxy(0,inode) * enrvalderxy[0];
            enr_derxy2(1,dofcounter) = derxy2(1,inode) * enrval + funct(inode) * enrvalderxy2[1] + derxy(1,inode) * enrvalderxy[1] + derxy(1,inode) * enrvalderxy[1];
            enr_derxy2(2,dofcounter) = derxy2(2,inode) * enrval + funct(inode) * enrvalderxy2[2] + derxy(2,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[2];
            enr_derxy2(3,dofcounter) = derxy2(3,inode) * enrval + funct(inode) * enrvalderxy2[3] + derxy(0,inode) * enrvalderxy[1] + derxy(1,inode) * enrvalderxy[0];
            enr_derxy2(4,dofcounter) = derxy2(4,inode) * enrval + funct(inode) * enrvalderxy2[4] + derxy(0,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[0];
            enr_derxy2(5,dofcounter) = derxy2(5,inode) * enrval + funct(inode) * enrvalderxy2[5] + derxy(1,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[1];
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
      std::size_t dofcounter = 0;

      const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
      const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
      const std::size_t numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
      dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
        enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        for (std::size_t inode = 0; inode < numvirtualnode; ++inode)
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
      std::size_t dofcounter = 0;

      const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
      const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
      const std::size_t numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
      dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
        enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        for (std::size_t inode = 0; inode < numvirtualnode; ++inode)
        {
          const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
          enr_funct(dofcounter) = funct(inode) * enrval;
          for (std::size_t isd = 0; isd < 3; ++isd)
          {
            enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
          }
          dofcounter += 1;
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

  private:

    //! don't want copy constructor
    explicit ElementEnrichmentValues(const ElementEnrichmentValues&);

    //! don't want assignement operator
    ElementEnrichmentValues& operator = (const ElementEnrichmentValues& old);

    //! don't want default constructor
    explicit ElementEnrichmentValues();

    const DRT::Element&                        ele_;
    //  const RCP<XFEM::InterfaceHandle>           ih_;
    const XFEM::ElementDofManager&             dofman_;

    //! map holding values for enrichments
    std::map<XFEM::Enrichment, double>        enrvals_;
    //! map holding values for enrichments for every node
    std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes_;

//URSULA
    // member variables required for kink enrichment
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy_;
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy2_;
//URSULA

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
    const std::size_t          numnode,
    const int*                 nodeids
);

/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
inline std::size_t getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        const std::size_t)
{
    return dofman.NumDofPerField(field);
}
template<>
inline std::size_t getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager&,
        const XFEM::PHYSICS::Field,
        const std::size_t numnode)
{
    return numnode;
}

//! determine number of parameter
template<std::size_t NUMNODE, XFEM::AssemblyType ASSTYPE>
struct NumParam
{
  static std::size_t get(
     const XFEM::ElementDofManager& dofman,
     const XFEM::PHYSICS::Field     field
     )
  {
    return dofman.NumDofPerField(field);
  }
};

/// determine number of parameter for standard_assembly
template<std::size_t NUMNODE> struct NumParam<NUMNODE,XFEM::standard_assembly>
{
  static std::size_t get(
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
