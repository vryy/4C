/*----------------------------------------------------------------------*/
/*! \file
\brief Communication management for particle data

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef PARTICLE_COMM_H_
#define PARTICLE_COMM_H_

#include "../drt_comm/comm_utils.H"


namespace INVANA
{
  /*! \class ParticleComm
   * \brief Communicator to communicate between groups
   *
   * Communicator class providing the neccessary MPIComms to communicate
   * globally, locally (within a group) and accros groups to be used by
   * the INVANA::ParticleGroup
   *
   * It just provides the global communicator gcomm_ and the local communicator
   * lcomm_ from the DRT::NestedParGroup given upon Init(). Furthermore, since
   * INVANA::ParticleData might be locally distributed, it furthermore provides
   * an intergroup communicator icomm_ to communicate between the respective
   * i-th local procs of each group. Equal group size is required therefore!
   *
   * \author kehl \date 06/2016
   *
   */
  class ParticleComm
  {
   public:
    //! @name Construction, destruction, setup
    //@{

    //! constructor
    ParticleComm();

    //! destructor
    virtual ~ParticleComm() {}

    //! init
    void Init(Teuchos::RCP<Epetra_Comm> gcomm,  /// [in] the global comm
        Teuchos::RCP<Epetra_Comm> lcomm,        /// [in] the local comm
        int lnumparticles,                      /// [in] local number of particles
        int group                               /// [in] group affiliation
    );

    //! setup
    void Setup();

    //@}

    //! @name Access
    //@{

    //! local comm
    const Epetra_Comm& LComm() const { return *lcomm_; }

    //! inter comm
    const Epetra_Comm& IComm() const { return *icomm_; }

    //! global comm
    const Epetra_Comm& GComm() const { return *gcomm_; }
    //@}

   private:
    //! @name Communicators
    //@{

    //! global communicator
    Teuchos::RCP<Epetra_Comm> gcomm_;

    //! local communicator
    Teuchos::RCP<Epetra_Comm> lcomm_;

    //! inter group communicator
    Teuchos::RCP<Epetra_Comm> icomm_;

    //@}

    //! this groups number of particles
    int lnumparticles_;

    //! this group
    int group_;
  };


}  // namespace INVANA
#endif /* PARTICLE_COMM_H_ */
