/*!----------------------------------------------------------------------
\file teko_baciepetraoperatorwrapper.H

\brief wrapper class around the epetra operator in baci for use with teko
\level 2
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef TEKO_BACIEPETRAOPERATORWRAPPER_H_
#define TEKO_BACIEPETRAOPERATORWRAPPER_H_

#ifdef HAVE_TEKO

#include <Thyra_EpetraLinearOp.hpp>

#include <Teko_EpetraOperatorWrapper.hpp>


#include "../../linalg/linalg_mapextractor.H"
#include "../../linalg/linalg_blocksparsematrix.H"

namespace LINALG
{
  namespace SOLVER
  {
    namespace TEKO
    {
      class Teko_BACIMappingStrategy : public Teko::Epetra::MappingStrategy
      {
       public:
        Teko_BACIMappingStrategy(const LINALG::MultiMapExtractor& rangemaps,
            const LINALG::MultiMapExtractor& domainmaps);

        /** Virtual function defined in MappingStrategy.  This copies
         * an Epetra_MultiVector into a Thyra::MultiVectorBase with
         * blocking handled by the strides defined in the constructor.
         *
         * \param[in]     epetra_X  source Epetra_MultiVector
         * \param[in,out]     thyra_X   destination Thyra::MultiVectorBase
         */
        virtual void copyEpetraIntoThyra(const Epetra_MultiVector& epetra_X,
            const Teuchos::Ptr<Thyra::MultiVectorBase<double>>& thyra_X) const;

        /** Virtual function defined in MappingStrategy.  This copies
         * an Epetra_MultiVector into a Thyra::MultiVectorBase with
         * blocking handled by the strides defined in the constructor.
         *
         * \param[in]     thyra_Y  source Thyra::MultiVectorBase
         * \param[in,out]     epetra_Y destination Epetra_MultiVector
         */
        virtual void copyThyraIntoEpetra(
            const Teuchos::RCP<const Thyra::MultiVectorBase<double>>& thyra_Y,
            Epetra_MultiVector& epetra_Y) const;

        /** \brief Domain map for this strategy */
        const Teuchos::RCP<const Epetra_Map> domainMap() const { return domainmaps_.FullMap(); };

        /** \brief Range map for this strategy */
        const Teuchos::RCP<const Epetra_Map> rangeMap() const { return rangemaps_.FullMap(); };

        /** \brief Identifier string */
        std::string toString() const { return std::string("Teko_BACIMappingStrategy"); }

       private:
        const LINALG::MultiMapExtractor rangemaps_;
        const LINALG::MultiMapExtractor domainmaps_;
      };

      class Teko_BACIEpetraOperatorWrapper : public Teko::Epetra::EpetraOperatorWrapper
      {
       public:
        Teko_BACIEpetraOperatorWrapper(const Teuchos::RCP<LINALG::BlockSparseMatrixBase>& mat);
        virtual ~Teko_BACIEpetraOperatorWrapper() { ; };

        void BuildBlockedOperator();

        Teuchos::RCP<const Thyra::EpetraLinearOp> GetThyraBlock(int r, int c);

       private:
        Teuchos::RCP<Teko_BACIMappingStrategy> mappingStrategy_;
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat_;
      };
    } /* namespace TEKO */
  }   /* namespace SOLVER */
} /* namespace LINALG */



#endif /* HAVE_TEKO */

#endif /* TEKO_BACIEPETRAOPERATORWRAPPER_H_ */
