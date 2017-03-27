/*---------------------------------------------------------------------*/
/*!
\file contact_aug_steepest_ascent_interface.H

\brief Steepest ascent interface based on the augmented contact
       formulation.

\level 3

\maintainer Michael Hiermeier

\date Mar 7, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_INTERFACE_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_INTERFACE_H_

#include "contact_augmented_interface.H"

namespace CONTACT {
namespace AUG {
namespace STEEPESTASCENT {

  /** \brief Augmented steepest ascent contact interface class
   *
   *  \author hiermeier \date 03/17 */
  class Interface : public ::CONTACT::AUG::Interface
  {
  public:
    /** \brief Alternative constructor
     *
     *  A prerequisite for this constructor is, that the passed
     *  shared interface data object has been filled/initialized already.
     *
     *  \param idata_ptr (in) : filled shared augmented contact interface
     *                          data container object
     *
     *  \author hiermeier \date 03/17 */
    Interface( const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr );

    /// constructor
    Interface(
        const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr,
        int id,
        const Epetra_Comm& comm,
        int dim,
        const Teuchos::ParameterList& icontact,
        bool selfcontact,
        INPAR::MORTAR::RedundantStorage redundant );

    /// destructor
    virtual ~Interface() {};

  protected:
    /// derived
    virtual void AssembleDGGLinMatrixOnSlaveSide(
        const CoNode& cnode,
        const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGSlLinMatrix ) const;

    /// derived
    virtual void AssembleDGGLinMatrixOnMasterSide(
        const CoNode& cnode,
        const std::map<int,std::pair<int,double> >& varWGapMaMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGMaLinMatrix ) const;

    /// derived
    virtual void AddKappaLinToGapLinearization(
        const std::map<int,double>& kappaLinMap,
        double x,
        double kappainv,
        double varWGap,
        double scale,
        std::map<int,double>& aWGapLinMap ) const;

  };  // class Interface

} // namespace STEEPESTASCENT
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_STEEPEST_ASCENT_INTERFACE_H_ */
