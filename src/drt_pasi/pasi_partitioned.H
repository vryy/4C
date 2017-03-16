/*!----------------------------------------------------------------------
\file pasi_partitioned.H

\brief partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_H
#define PASI_PARTITIONED_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
class Particle;
class Structure;
class StructureBaseAlgorithmNew;
class PASIStructureWrapper;
} // namespace ADAPTER

namespace PARTICLE
{
class Algorithm;
} // namespace PARTICLE

/*----------------------------------------------------------------------*
 | partitioned algorithm for pasi                        sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
namespace PASI
{

//! partitioned algorithm for particle structure interaction
class PartitionedAlgo: public ADAPTER::AlgorithmBase
{

public:

  //! constructor
  explicit PartitionedAlgo(
      const Epetra_Comm&              comm,        //! communicator
      const Teuchos::ParameterList&   pasi_params  //! input parameters for particle structure interaction
      );

  //! virtual destructor
  virtual ~PartitionedAlgo(){};

  /*! \brief Initialize this object

  Hand in all objects/parameters/etc. from outside.
  Construct and manipulate internal objects.

  \note Try to only perform actions in Init(), which are still valid
        after parallel redistribution of discretizations.
        If you have to perform an action depending on the parallel
        distribution, make sure you adapt the affected objects after
        parallel redistribution.
        Example: cloning a discretization from another discretization is
        OK in Init(...). However, after redistribution of the source
        discretization do not forget to also redistribute the cloned
        discretization.
        All objects relying on the parallel distribution are supposed to
        the constructed in \ref Setup().

  \warning none
  \return int
  \date 01/17
  \author rauch  */
  virtual void Init(
      const Epetra_Comm& comm  //! communicator
      );

  /*! \brief Setup all class internal objects and members

   Setup() is not supposed to have any input arguments !

   Must only be called after Init().

   Construct all objects depending on the parallel distribution and
   relying on valid maps like, e.g. the state vectors, system matrices, etc.

   Call all Setup() routines on previously initialized internal objects and members.

  \note Must only be called after parallel (re-)distribution of discretizations is finished !
        Otherwise, e.g. vectors may have wrong maps.

  \warning none
  \return void
  \date 01/17
  \author rauch  */
  virtual void Setup();

  //! read restart data
  virtual void ReadRestart(
      int step  //! time step for restart
      );

  //! timeloop of coupled problem
  virtual void Timeloop() = 0;

  //! test results (if necessary)
  virtual void TestResults(
      const Epetra_Comm& comm  //! local epetra communicator
      );

  //! access to structural field
  const Teuchos::RCP< ::ADAPTER::PASIStructureWrapper>& StructureField() const {return structure_;};

  //! access to particle field
  const Teuchos::RCP<PARTICLE::Algorithm>&  ParticleField() const {return particles_;};

  //! returns true if Setup() was called and is still valid
  bool IsSetup(){return issetup_;};

  //! returns true if Init(..) was called and is still valid
  bool IsInit(){return isinit_;};

protected:

  //! prepare time step
  virtual void PrepareTimeStep(bool printheader);

  //! structural time step
  virtual void StructStep();

  //! particle time step
  virtual void ParticleStep();

  //! set structural displacements
  virtual void SetStructDispVel(
      Teuchos::RCP<const Epetra_Vector> dispnp,  //! structural displacements \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> velnp    //! structural velocities at \f$t_{n+1}\f$
      );

  //! structural output
  virtual void StructOutput();

  //! particle output
  virtual void ParticleOutput();

  //! check if \ref Setup() was called
  void CheckIsSetup()
  {
    if (not IsSetup())
      dserror("Setup() was not called.");
  };

  //! check if \ref Init() was called
  void CheckIsInit()
  {
    if (not IsInit())
      dserror("Init(...) was not called.");
  };

  //! ptr to underlying structure
  Teuchos::RCP< ::ADAPTER::PASIStructureWrapper> structure_;

  //! ptr to underlying particle algorithm
  Teuchos::RCP<PARTICLE::Algorithm> particles_;

private:

  //! ptr to the underlying structure problem base algorithm
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr_;

  //! flag indicating if class is setup
  bool issetup_;

  //! flag indicating if class is initialized
  bool isinit_;

  //! set flag true after setup or false if setup became invalid
  void SetIsSetup(bool trueorfalse){issetup_ = trueorfalse;};

  //! set flag true after init or false if init became invalid
  void SetIsInit(bool trueorfalse){isinit_ = trueorfalse;};

}; // class PartitionedAlgo: public ADAPTER::AlgorithmBase

} // namespace PASI

/*----------------------------------------------------------------------*/
#endif // PASI_PARTITIONED_H
