/*----------------------------------------------------------------------*/
/*!
\file contact_wear_interface.cpp

\brief  ...
<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef CONTACT_WEAR_INTERFACE_H_
#define CONTACT_WEAR_INTERFACE_H_


#include "contact_interface.H"

namespace CONTACT
{

class WearInterface : public CoInterface
{

public:

  /*!
  \brief Constructor

  */
  WearInterface(const int id, const Epetra_Comm& comm, const int dim,
              const Teuchos::ParameterList& icontact,
              bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

  /*!
  \brief Destructor

  */
  virtual ~WearInterface() {}

  /*!
  \brief Assemble second mortar D matrix for both-sided wear

  */
  virtual void AssembleD2(LINALG::SparseMatrix& dglobal);

  /*!
  \brief Assemble matrix S containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the normal contact condition g~ = 0. Concretely, this
  includes assembling the linearizations of the slave side
  nodal normals and of the Mortar matrices D  and M.

  */
  virtual void AssembleS(LINALG::SparseMatrix& sglobal);

  /*!
  \brief Assemble matrix LinStick containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the tangential stick condition delta tg = 0. Concretely, this
  includes assembling the linearizations of the slave side
  nodal tangents and of the Mortar matrices D  and M.

  */
  virtual void AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                LINALG::SparseMatrix& linstickDISglobal,
                                Epetra_Vector& linstickRHSglobal);
  /*!
  \brief Assemble matrix LinSlip containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the tangential slip condition. Concretely, this
  includes assembling the linearizations of the slave side
  nodal tangents and of the Mortar matrices D  and M.

  */
  virtual void AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
                               LINALG::SparseMatrix& linslipDISglobal,
                               Epetra_Vector& linslipRHSglobal);

  /*!
  \brief Assemble matrices W containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the normal contact and slip contact condition for ~w.
  --> w.r.t. lagr. mult.

  */
  virtual void AssembleLinWLm(LINALG::SparseMatrix& sglobal);
  virtual void AssembleLinWLmSl(LINALG::SparseMatrix& sglobal);
  virtual void AssembleLinWLmSt(LINALG::SparseMatrix& sglobal);

  /*!
  \brief Assemble wear w

   This method assembles the weighted wear vector.
   */
   virtual void AssembleWear(Epetra_Vector& gglobal);

   /*!
   \brief Build active set (nodes / dofs) of this interface

   If the flag init==true, the active set is initialized (for t=0)
   according to the contact initialization defined in the input file.

   */
   virtual bool BuildActiveSet(bool init = false);

   /*!
   \brief Initialize / reset interface for contact

   Derived version!

   */
   virtual void Initialize();


   /*!
   \brief Returning dofs for both-sided wear mapping

   */
   virtual const Teuchos::RCP<Epetra_Map> InvolvedDofs() const
   {
     return involveddofs_;
   }

   virtual const Teuchos::RCP<Epetra_Map> InvolvedNodes() const
   {
     return involvednodes_;
   }

   /*!
   \brief Set element areas

   Derived version!

   */
   void SetElementAreas();

protected:

   // both-sided wear specific stuff
   Teuchos::RCP<Epetra_Map>      involvednodes_;  // row map of all involved master nodes
   Teuchos::RCP<Epetra_Map>      involveddofs_;   // row map of all involved master dofs

   bool wearimpl_;  //bool for implicit wear
}; //class


}

#endif /* CONTACT_WEAR_INTERFACE_H_ */
