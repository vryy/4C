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
  class Interface : public ::CONTACT::AUG::Interface
  {
  public:
    /// constructor
    Interface(
        int id,
        const Epetra_Comm& comm,
        int dim,
        const Teuchos::ParameterList& icontact,
        bool selfcontact,
        INPAR::MORTAR::RedundantStorage redundant );

    /// destructor
    virtual ~Interface() {};

  protected:
    virtual void AssembleDGGLinMatrixOnSlaveSide(
        const CoNode& cnode,
        const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGSlLinMatrix ) const;

    virtual void AssembleDGGLinMatrixOnMasterSide(
        const CoNode& cnode,
        const std::map<int,std::pair<int,double> >& varWGapMaMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGMaLinMatrix ) const;

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
