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


/*----------------------------------------------------------------------*/
/* headers */
#include "smc_particle.H"
#include <boost/random.hpp>
// to get DRT problem
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io_pstream.H"

// only compile this on the workstation as kaisers boost version is outdated an cant run this code
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)


STR::INVANA::SMCParticleType STR::INVANA::SMCParticleType::instance_;



DRT::ParObject* STR::INVANA::SMCParticleType::Create( const std::vector<char> & data )
{
  int dummysize =1;
  STR::INVANA::SMCParticle * my_particle =new STR::INVANA::SMCParticle(dummysize);
  my_particle->Unpack(data);
  return my_particle;
}


/*----------------------------------------------------------------------*/
/* standard constructor */
/*----------------------------------------------------------------------*/
STR::INVANA::SMCParticle::SMCParticle(int numparams)
{
  weightprop_=0.0;
  weight_=0.0;
  weightold_=0.0;
  normweight_=0.0;
  LogLikeGammaProp_=0.0;
  LogLikeGamma_=0.0;
  LogLikeGammaOld_=0.0;
  LogLikeProp_=0.0;
  LogLike_=0.0;
  LogLikeOld_=0.0;
  LogPrior_= 0.0;
  LogPriorProp_= 0.0;
  position_.resize(numparams,0.0);
  positionprop_.resize(numparams,0.0);

}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       biehler 10/13|
 *----------------------------------------------------------------------*/
STR::INVANA::SMCParticle::SMCParticle(const STR::INVANA::SMCParticle& old):
ParObject(old),
weightprop_(old.weightprop_),
weight_(old.weight_),
weightold_(old.weightold_),
normweight_(old.normweight_),
LogLikeGammaProp_(old.LogLikeGammaProp_),
LogLikeGamma_(old.LogLikeGamma_),
LogLikeGammaOld_(old.LogLikeGammaOld_),
LogLikeProp_(old.LogLikeProp_),
LogLike_(old.LogLike_),
LogLikeOld_(old.LogLikeOld_),
LogPrior_(old.LogPrior_),
LogPriorProp_(old.LogPriorProp_),
positionprop_(old.positionprop_),
position_(old.position_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                         biehler 10/13|
 *----------------------------------------------------------------------*/
Teuchos::RCP<STR::INVANA::SMCParticle> STR::INVANA::SMCParticle::Clone() const
{
  Teuchos::RCP<STR::INVANA::SMCParticle> newparticle = Teuchos::rcp(new STR::INVANA::SMCParticle(*this));
  return newparticle;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         biehler 10/13|
 *----------------------------------------------------------------------*/
void STR::INVANA::SMCParticle::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);



  // get the rest
  AddtoPack(data,weightprop_);
  AddtoPack(data,weight_);
  AddtoPack(data,weightold_);
  AddtoPack(data,normweight_);
  AddtoPack(data,LogLikeGammaProp_);
  AddtoPack(data,LogLikeGamma_);
  AddtoPack(data,LogLikeGammaOld_);
  AddtoPack(data,LogLikeProp_);
  AddtoPack(data,LogLike_);
  AddtoPack(data,LogLikeOld_);
  AddtoPack(data,LogPrior_);
  AddtoPack(data,LogPriorProp_);

  // add  position of particle
  AddtoPack(data,position_);
  AddtoPack(data,positionprop_);
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         biehler 10/13|
 *----------------------------------------------------------------------*/
void STR::INVANA::SMCParticle::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");


  // get the rest
  weightprop_=ExtractDouble(position,data);
  weight_=ExtractDouble(position,data);
  weightold_=ExtractDouble(position,data);
  normweight_=ExtractDouble(position,data);
  LogLikeGammaProp_=ExtractDouble(position,data);
  LogLikeGamma_=ExtractDouble(position,data);
  LogLikeGammaOld_=ExtractDouble(position,data);
  LogLikeProp_=ExtractDouble(position,data);
  LogLike_=ExtractDouble(position,data);
  LogLikeOld_=ExtractDouble(position,data);
  LogPrior_=ExtractDouble(position,data);
  LogPriorProp_=ExtractDouble(position,data);

  // extract position of particle
  ExtractfromPack(position,data,position_);
  ExtractfromPack(position,data,positionprop_);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                             biehler 10/13|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const STR::INVANA::SMCParticle& particle)
{
  particle.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print current position an weight of this particle (public)biehler 10/13|
 *----------------------------------------------------------------------*/
void STR::INVANA::SMCParticle::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Particle   weight: " << weight_ << " normweight: "  << normweight_ << "position: ";

  for(unsigned int i=0;i<position_.size();i++)
  {
         os  << position_.at(i) << " " ;
  }

  return;
}

void STR::INVANA::SMCParticle::SetPosition(std::vector<double> new_position)
{
  position_= new_position;
}

void STR::INVANA::SMCParticle::SetPositionProp(std::vector<double> new_position)
{
  positionprop_= new_position;
}

//void STR::INVANA::SMCParticle::SetLogPriorProp(double log_prior)/
//{
//  LogPriorProp_= log_prior;
//}

//void STR::INVANA::SMCParticle::SetLogPrior(double log_prior)
//{
////  LogPrior_= log_prior;
//}

void STR::INVANA::SMCParticle::ComputeWeightProp(double gamma)
{
  weightprop_=weight_ + gamma*LogLikeProp_ - LogLikeGamma_;
}

void STR::INVANA::SMCParticle::AcceptReweightProp(double gamma)
{

  weightold_=weight_;
  weight_=weightprop_;

  LogLikeGammaOld_=LogLikeGamma_;
  LogLikeGamma_=gamma*LogLikeProp_;

  LogLikeOld_=LogLike_;
  LogLike_=LogLikeProp_;
}

void STR::INVANA::SMCParticle::AcceptMoveProp()
{
  LogLikeGammaOld_=LogLikeGamma_;
  LogLikeGamma_=LogLikeGammaProp_;

  LogLikeOld_=LogLike_;
  LogLike_=LogLikeProp_;

  LogPrior_=LogPriorProp_;
  position_=positionprop_;
}


#else
 // no code here
#endif
