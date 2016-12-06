/*!-----------------------------------------------------------------------------------------------*
\file timeInt_enr.H

\brief provides the enrichment computation class(es)

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#ifndef TIMEINT_ENR_H_
#define TIMEINT_ENR_H_


#include "timeInt.H"


namespace XFEM
{

/*!
  \brief this class is used to compute enrichment values for the new timestep

 */
class EnrichmentProjection : public ENR
{

public:

  //! constructor
  explicit EnrichmentProjection(
      XFEM::TIMEINT& timeInt,
      const double& veljump,
      INPAR::COMBUST::XFEMTimeIntegrationEnr timeIntEnr,
      INPAR::COMBUST::XFEMTimeIntegrationEnrComp timeIntEnrType
  );

  //! destructor
  ~EnrichmentProjection()
  {
    return;
  }

  //! main function computing the new enrichment values
  void compute(
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
      std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
  );

  std::map<int,std::vector<LINALG::Matrix<2,4> > >* oldJumpAndKinkValues(){return &eleJumpsAndKinks_;}

protected:

  //! disabled assignment operator
  EnrichmentProjection operator = (const EnrichmentProjection& old);

  //! call computation of different values for old interface position
  void callOldValues(
  );

  //! computation of different values for old interface position
  template <const int numnode>
  void oldValues(
      const DRT::Element*& ele
  );

  //! computation of jump and kink values for jump enrichment at old interface position
  template<const int numnode>
  void oldJumpAndKinkValues(
      const DRT::Element* ele,
      LINALG::Matrix<1,numnode>& dmin,
      LINALG::Matrix<3,1> normal,
      std::vector<LINALG::Matrix<4,numnode> >& kinkEnrValues,
      int& numnodeused
  );

  //! computation of kink values for kink enrichment at old interface position
  template<const int numnode>
  void oldKinkValues(
      const DRT::Element* ele,
      LINALG::Matrix<1,numnode>& dmin,
      std::vector<LINALG::Matrix<4,numnode> >& kinkEnrValues,
      int& numnodeused
  );

  //! standard computation of new enrichment values at new interface position
  void computeNewEnrichments(
  );

  //! alternative computation of new enrichment values at new interface position
  void handleFailedNodes(
  );

  //! computation of new enrichment values at new interface position for jump enriched dofs
  void computeJumpEnrichmentValues(
      const DRT::Node* node,
      std::vector<LINALG::Matrix<2,4> > jumpsAndKinks
  );

  //! computation of new enrichment values at new interface position for kink enriched dofs
  void computeKinkEnrichmentValues(
      const DRT::Node* node,
      std::vector<LINALG::Matrix<1,4> > kinks
  );

  void setJumpEnrichmentValues(
      const DRT::Node* node
  );

  void setKinkEnrichmentValues(
      const DRT::Node* node
  );

#ifdef PARALLEL
  //! export enrichment data to the neighbour processor
  void exportEnrichmentData(
  );
#endif

  // velocity jump
  const double veljump_;

}; // class EnrichmentProjection
} // namespace XFEM

#endif /*TIMEINT_ENR_H_*/
