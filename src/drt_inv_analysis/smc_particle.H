/*----------------------------------------------------------------------*/
/*!
\file smc_particle.H
\brief Class to handle particles for smc algorithm

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef SMC_PARTICLE_H
#define SMC_PARTICLE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_inpar/inpar_material.H"
#include "Epetra_SerialDenseVector.h"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_parobject.H"

// forward declarations
namespace DRT
{
  class Discretization;
}
/*----------------------------------------------------------------------*/
//! Material parameters manager

namespace INVANA

{

  class SMCParticleType : public DRT::ParObjectType
  {
    public:

    std::string Name() const { return "SMCParticleType"; }

    static SMCParticleType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

    private:

    static SMCParticleType instance_;
  };

  class SMCParticle : public DRT::ParObject
  {

  public:

  //! brief standard constructor
  SMCParticle(int numparams);

  //! brief Copy Constructor Makes a deep copy of a ParticleNode
  SMCParticle(const INVANA::SMCParticle& old);

  //! brief clone particle and return rcp to this particle
  Teuchos::RCP<INVANA::SMCParticle> Clone() const;

  //! brief Return unique ParObject id. every class implementing ParObject needs a unique id defined at the
  //! top of drt_lib/drt_parobject.H.
  virtual int UniqueParObjectId() const { return SMCParticleType::Instance().UniqueParObjectId(); }

  //! brief  Pack this class so it can be communicated
  virtual void Pack(DRT::PackBuffer& data) const;

 //! brief   Unpack data from a char vector into this class
  virtual void Unpack(const std::vector<char>& data);

  //! brief  Print this Particle
  virtual void Print(std::ostream& os) const ;

  //! brief destructor
  virtual ~SMCParticle() { ; };

  //! brief return current unnormalized weight
  double GetWeight() {return weight_;};

  //! brief set new weight
  void SetWeight(double new_weight){weight_=new_weight;};

  //! brief get dimension of space the particle moves in
  int GetSizeOfPosition() {return position_.size();};

  //! brief return current position of particle
  std::vector<double> GetPosition() {return position_;};

  //! brief return proposed new position of particle
  std::vector<double> GetPropPosition() {return positionprop_;};

  //! brief set  position of particle
  virtual void SetPosition(std::vector<double> new_position);

  //! brief set new proposal position of particle
  virtual void SetPositionProp(std::vector<double>  new_position);

  //! brief comput new weight
  void ComputeWeightProp(double gamma);

  //! brief update particle data after acceptance of reweighting step
  void AcceptReweightProp(double gamma);

  //! brief update particle data after acceptance of proposed move
  void AcceptMoveProp();

  //! brief  unnormalized weight of proposal step (reweight)
  double weightprop_;

  //! brief  current unnormalized weight of particle
  double weight_;

  //! brief  old unnormalized weight of particle
  double weightold_;

  //! brief normalized weight of particle
  double normweight_;

  //! brief  log likelihood of proposal both reweight and move
  double LogLikeGammaProp_;

  //! brief  log likelihood of current position
  double LogLikeGamma_;

  //! brief  log likelihood of last position
  double LogLikeGammaOld_;

  //! brief  log likelihood of proposal both reweight and move
  double LogLikeProp_;

  //! brief  log likelihood of current position
  double LogLike_;

  //! brief  log likelihood of last position
  double LogLikeOld_;

   //! brief  log likelihood of current position
  double LogPrior_;

  //! brief  log prior of proposal  move
  double LogPriorProp_;

  protected:

  private:

  //! brief proposal position of particle in random space
  std::vector<double> positionprop_;

  //! brief current position of particle in random space
  std::vector<double> position_;



    //seed for random variable generation
    unsigned int prior_seed_;
  }; //EOF class

} // EOF INVANA

std::ostream& operator << (std::ostream& os, const INVANA::SMCParticle& particle);

#endif


#else
 // no code here
#endif
