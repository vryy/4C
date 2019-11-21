/*---------------------------------------------------------------------------*/
/*! \file
\brief partitioned algorithm for particle structure interaction

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_H
#define PASI_PARTITIONED_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace ADAPTER
{
  class StructureBaseAlgorithmNew;
  class PASIStructureWrapper;
}  // namespace ADAPTER

namespace PARTICLEALGORITHM
{
  class ParticleAlgorithm;
}

namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}  // namespace STR

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PASI
{
  /*!
   * \brief partitioned algorithm for particle structure interaction
   *
   * An abstract base class for partitioned particle structure interaction problems defining
   * methods and holding members to be used in derived algorithms.
   *
   * \author Sebastian Fuchs \date 01/2017
   */
  class PartitionedAlgo : public ADAPTER::AlgorithmBase
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 01/2017
     *
     * \param[in] comm   communicator
     * \param[in] params particle structure interaction parameter list
     */
    explicit PartitionedAlgo(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PartitionedAlgo() override{};

    /*!
     * \brief init pasi algorithm
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    virtual void Init();

    /*!
     * \brief setup pasi algorithm
     *
     * \author Sebastian Fuchs \date 01/2017
     */
    virtual void Setup();

    /*!
     * \brief read restart information for given time step
     *
     * \author Sebastian Fuchs \date 01/2017
     *
     * \param[in] restartstep restart step
     */
    void ReadRestart(int restartstep) override;

    /*!
     * \brief timeloop of coupled problem
     *
     * \author Sebastian Fuchs \date 01/2017
     */
    virtual void Timeloop() = 0;

    /*!
     * \brief perform result tests
     *
     * \author Sebastian Fuchs \date 01/2017
     *
     * \param[in] comm communicator
     */
    void TestResults(const Epetra_Comm& comm);

    //! get initialization status
    bool IsInit() { return isinit_; };

    //! get setup status
    bool IsSetup() { return issetup_; };

   protected:
    /*!
     * \brief prepare time step
     *
     * \author Sebastian Fuchs \date 01/2017
     *
     * \param[in] printheader flag to control output of time step header
     */
    void PrepareTimeStep(bool printheader = true);

    /*!
     * \brief structural time step
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void StructStep();

    /*!
     * \brief particle time step
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void ParticleStep();

    /*!
     * \brief post evaluate time step
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void PostEvaluateTimeStep();

    /*!
     * \brief extract interface states
     *
     * Extract the interface states displacement, velocity, and acceleration from the structural
     * states.
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void ExtractInterfaceStates();

    /*!
     * \brief set interface states
     *
     * Set the interface states displacement, velocity, and acceleration as handed in to the
     * particle wall handler. This includes communication, since the structural discretization and
     * the particle wall discretization are in general distributed independently of each other to
     * all processors.
     *
     * \author Sebastian Fuchs \date 02/2017
     *
     * \param[in] intfdispnp interface displacement
     * \param[in] intfvelnp  interface velocity
     * \param[in] intfaccnp  interface acceleration
     */
    void SetInterfaceStates(Teuchos::RCP<const Epetra_Vector> intfdispnp,
        Teuchos::RCP<const Epetra_Vector> intfvelnp, Teuchos::RCP<const Epetra_Vector> intfaccnp);

    /*!
     * \brief output of structure field
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void StructOutput();

    /*!
     * \brief output of particle field
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void ParticleOutput();

    //! check correct setup
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("pasi algorithm not setup correctly!");
    };

    //! check correct initialization
    void CheckIsInit()
    {
      if (not IsInit()) dserror("pasi algorithm not initialized correctly!");
    };

    //! structural field
    Teuchos::RCP<::ADAPTER::PASIStructureWrapper> structurefield_;

    //! particle algorithm
    Teuchos::RCP<PARTICLEALGORITHM::ParticleAlgorithm> particlealgorithm_;

    //! communication object at the interface
    Teuchos::RCP<const STR::AUX::MapExtractor> interface_;

    //! interface displacement
    Teuchos::RCP<Epetra_Vector> intfdispnp_;

    //! interface velocity
    Teuchos::RCP<Epetra_Vector> intfvelnp_;

    //! interface acceleration
    Teuchos::RCP<Epetra_Vector> intfaccnp_;

   private:
    /*!
     * \brief init structure field
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void InitStructureField();

    /*!
     * \brief init particle algorithm
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void InitParticleAlgorithm();

    /*!
     * \brief build and register structure model evaluator
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void BuildStructureModelEvaluator();

    //! ptr to the underlying structure problem base algorithm
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr_;

    //! flag indicating correct initialization
    bool isinit_;

    //! flag indicating correct setup
    bool issetup_;

    //! set flag indicating correct initialization
    void SetIsInit(bool isinit) { isinit_ = isinit; };

    //! set flag indicating correct setup
    void SetIsSetup(bool issetup) { issetup_ = issetup; };
  };

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
