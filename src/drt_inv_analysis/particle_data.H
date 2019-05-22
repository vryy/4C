/*----------------------------------------------------------------------*/
/*!
\brief Data for particles in the Sequential Monte Carlo algortihm

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#ifndef PARTICLE_DATA_H_
#define PARTICLE_DATA_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace INVANA
{
  /*! \class ParticleDataTye
   * \brief ParObjectType for the the creation of ParticleData
   *
   * \author kehl \date 06/2016
   */
  class ParticleDataType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleDataType"; }

    static ParticleDataType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleDataType instance_;
  };


  /*! \class ParticleData
   * \brief Data for each particle of the smc algorithm
   *
   * This is the data container for one particle of the INVANA::Optimzier - algorithm.
   * Therefore it stores the state and the correponding prior- and posterior-values.
   *
   * The ParticelData can be queried for the prior and posterior values only in case
   * a state has been evluated beforehand. Therefore, whenever a new state is set prior
   * and posterior have to be set accordingly afterwards.
   *
   * It derives from DRT::ParObject since the set of particles in the smc setting
   * must be resampled eventually and thus they must be sent around between the
   * INVANA::ParticleGroups.
   *
   * The particle state state_ is a locally (within a nested par group) distributed
   * quantity and is only packed locally! The INVANA::ParticleComm must be created
   * such that it can be recreated entirely in another group. The Unpack routine
   * implicitly assumes that a particle is only sent around with a properly evaluated
   * state. I cannot imagine a situation at the moment where this shoudl not apply.
   *
   * \author kehl \date 06/2016
   */
  class ParticleData : public DRT::ParObject
  {
   public:
    //! @name Construction, destruction, setup
    //@{

    //! constructor
    ParticleData();

    //! copy constructor
    ParticleData(ParticleData& data);

    //! assignment
    ParticleData& operator=(ParticleData& data);

    //! destructor
    virtual ~ParticleData() {}

    void Init(const Epetra_BlockMap& map);

    //@}

    //! @name DRT::ParObject methods
    //@{

    //! unique type id
    virtual int UniqueParObjectId() const
    {
      return ParticleDataType::Instance().UniqueParObjectId();
    }

    //! Pack
    virtual void Pack(DRT::PackBuffer& data) const;

    //! Unpack
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Data manipulation
    //@{

    //! set state, values and weight
    void SetState(const Epetra_Vector& state  ///[in] state in probability space
    );

    //! set values according to the current state
    void SetData(double& posterior,  ///[in] corresponding posterior value
        double& prior                ///[in] corresponding prior value
    );

    //@}

    //! @name Get data
    //@{

    //! Get posterior value
    double GetPosterior()
    {
      if (statechanged_) dserror("state changed!");
      return posterior_;
    }

    //! Get prior value
    double GetPrior()
    {
      if (statechanged_) dserror("state changed!");
      return prior_;
    }

    //! Get state (not to be manipulated!)
    const Epetra_Vector& GetState() const { return *state_; }

    //@}

   private:
    //! @name Data managed by this class
    //@{

    //! state of the particle in probability space
    Teuchos::RCP<Epetra_Vector> state_;

    //! posterior value
    double posterior_;

    //! prior value
    double prior_;

    //@}

    //! checking changes in state
    bool statechanged_;

    //! communicator within the group this data lives
    Teuchos::RCP<Epetra_Comm> lcomm_;
  };
}  // namespace INVANA

#endif /* PARTICLE_DATA_H_ */
