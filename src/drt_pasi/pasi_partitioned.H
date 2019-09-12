/*---------------------------------------------------------------------------*/
/*! \file
\brief partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_H
#define PASI_PARTITIONED_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 01/2017 |
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

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
namespace PASI
{
  //! partitioned algorithm for particle structure interaction
  class PartitionedAlgo : public ADAPTER::AlgorithmBase
  {
   public:
    //! constructor
    explicit PartitionedAlgo(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PartitionedAlgo() override{};

    //! init pasi algorithm
    virtual void Init();

    //! setup pasi algorithm
    virtual void Setup();

    //! read restart information for given time step
    void ReadRestart(int restartstep) override;

    //! timeloop of coupled problem
    virtual void Timeloop() = 0;

    //! perform result tests
    void TestResults(const Epetra_Comm& comm);

    //! get initialization status
    bool IsInit() { return isinit_; };

    //! get setup status
    bool IsSetup() { return issetup_; };

   protected:
    //! prepare time step
    void PrepareTimeStep(bool printheader = true);

    //! structural time step
    void StructStep();

    //! particle time step
    void ParticleStep();

    //! set structural states
    void SetStructStates();

    //! output of structure field
    void StructOutput();

    //! output of particle field
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

   private:
    //! init structure field
    void InitStructureField();

    //! init particle algorithm
    void InitParticleAlgorithm();

    //! build and register structure model evaluator
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
