/*---------------------------------------------------------------------*/
/*!
\file contact_aug_steepest_ascent_interface.cpp

\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 7, 2017

*/
/*---------------------------------------------------------------------*/


#include "contact_aug_steepest_ascent_interface.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr )
    : ::CONTACT::AUG::Interface( idata_ptr )
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::STEEPESTASCENT::Interface::Interface(
    const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr,
    int id,
    const Epetra_Comm& comm,
    int dim,
    const Teuchos::ParameterList& icontact,
    bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant )
    : ::CONTACT::AUG::Interface(idata_ptr,id,comm,dim,icontact,selfcontact,redundant)
{
  /* left blank, nothing to do here */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Interface::AssembleDGGLinMatrixOnSlaveSide(
    const CoNode& cnode,
    const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap,
    const std::map<int,double>& aWGapLinMap,
    double cn,
    double aWGap,
    LINALG::SparseMatrix& dGGSlLinMatrix ) const
{
  typedef GEN::pairedvector<int,std::pair<int,double> >::const_iterator CI;

  // iteration over ALL slave Dof Ids
  for ( CI p=varWGapSlMap.begin(); p!=varWGapSlMap.end(); ++p )
  {
    const int sRow = p->first;

    // *** linearization of varWGap w.r.t. displacements ***
    // nothing to do for the steepest ascent method

    // *** linearization of the averaged weighted gap w.r.t. displacements ***
    AssembleMapIntoMatrix( sRow, cn*(p->second).second, aWGapLinMap, dGGSlLinMatrix );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Interface::AssembleDGGLinMatrixOnMasterSide(
    const CoNode& cnode,
    const std::map<int,std::pair<int,double> >& varWGapMaMap,
    const std::map<int,double>& aWGapLinMap,
    double cn,
    double aWGap,
    LINALG::SparseMatrix& dGGMaLinMatrix ) const
{
  typedef std::map<int,std::pair<int,double> >::const_iterator CI;

  // iteration over ALL master Dof Ids
  for ( CI p=varWGapMaMap.begin(); p!=varWGapMaMap.end(); ++p )
  {
    const int mRow = p->first;

    // *** linearization of varWGap w.r.t. displacements ***
    // nothing to do for the steepest ascent method

    // *** linearization of the averaged weighted gap w.r.t. displacements ***
    AssembleMapIntoMatrix( mRow, -cn*(p->second).second, aWGapLinMap, dGGMaLinMatrix );
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::STEEPESTASCENT::Interface::AddKappaLinToGapLinearization(
    const std::map<int,double>& kappaLinMap,
    double x,
    double kappainv,
    double varWGap,
    double scale,
    std::map<int,double>& aWGapLinMap ) const
{
  // do nothing for the steepest ascent method
}
