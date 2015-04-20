/*----------------------------------------------------------------------*/
/*!
 \file structporo_reaction_ecm.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/



#include <vector>
#include "structporo_reaction_ecm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoroReactionECM::StructPoroReactionECM(Teuchos::RCP<MAT::PAR::Material> matdata) :
  StructPoroReaction(matdata)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoroReactionECM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoroReactionECM(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECMType MAT::StructPoroReactionECMType::instance_;

DRT::ParObject* MAT::StructPoroReactionECMType::Create(const std::vector<char> & data)
{
  MAT::StructPoroReactionECM* struct_poro = new MAT::StructPoroReactionECM();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECM::StructPoroReactionECM() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoroReactionECM::StructPoroReactionECM(MAT::PAR::StructPoroReactionECM* params) :
  StructPoroReaction(params),
  refporosity_old_(-1.0),
  refporositydot_old_(0.0),
  conc_m2c1_(0.0),
  conc_m2c1_old_(0.0),
  conc_m2c1dot_(0.0),
  conc_m2c1dot_old_(0.0),
  params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Setup(int numgp,DRT::INPUT::LineDefinition* linedef)
{
  StructPoroReaction::Setup(numgp,linedef);
  refporosity_old_ = params_->initporosity_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);

  // refporosity_
  AddtoPack(data, refporosity_old_);
  // refporositydot_old_
  AddtoPack(data, refporositydot_old_);
  // conc_m2c1_
  AddtoPack(data, conc_m2c1_);
  // conc_m2c1_old_
  AddtoPack(data, conc_m2c1_old_);
  // conc_m2c1dot_
  AddtoPack(data, conc_m2c1dot_);
  // conc_m2c1dot_old_
  AddtoPack(data, conc_m2c1dot_old_);

  // add base class material
  StructPoroReaction::Pack(data);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::StructPoroReactionECM*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  ExtractfromPack(position,data,refporosity_old_);
  ExtractfromPack(position,data,refporositydot_old_);
  ExtractfromPack(position,data,conc_m2c1_);
  ExtractfromPack(position,data,conc_m2c1_old_);
  ExtractfromPack(position,data,conc_m2c1dot_);
  ExtractfromPack(position,data,conc_m2c1dot_old_);


  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  StructPoroReaction::Unpack(basedata);

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Reaction(double cnp,
                                          double porosity,
                                          Teuchos::ParameterList& params)
{
  double dt = params.get<double>("delta time",-1.0);
  //double time = params.get<double>("total time",-1.0);

//  if(time==-1.0)
//    dserror("time step or total time not available");

  if(dt < 0.0)
    dserror("time step not available");

  double theta = 0.66;

  double k_on_m2c1 = 1.0;
  double k_off_m2c1 = 1.0;

  double k_on_c1 = 1.0;
  double k_off_c1 = 1.0;

  conc_m2c1_ = (  conc_m2c1_old_
               + dt*theta * k_on_m2c1*porosity*(1.0-porosity)*cnp
               - dt*(1.0-theta)*conc_m2c1dot_old_
               )/(1.0+dt*theta*k_off_m2c1);

  conc_m2c1dot_ = k_on_m2c1*porosity*cnp - k_off_m2c1 *conc_m2c1_;

  refporosity_ =  refporosity_old_
               - dt*theta * ( - k_on_c1*porosity*(1.0-porosity)*cnp + k_off_c1*conc_m2c1_ )
               + dt*(1.0-theta)*refporositydot_old_
               ;
  std::cout<<"refporosity_ECM Reaction: "<<refporosity_<<std::endl;
  //Todo: multiply by J!
  refporositydot_ = - k_on_c1*porosity*(1.0-porosity)*cnp + k_off_c1 *conc_m2c1_;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::StructPoroReactionECM::BodyForceTerm(double porosity) const
{
  double k_on_m2c1 = 1.0;
  double bodyforce = k_on_m2c1 * (1.0-porosity) * conc_m2c1_;
  return bodyforce;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::Update()
{
  refporosity_old_ = refporosity_;
  refporositydot_old_ = refporositydot_;
  conc_m2c1_old_ = conc_m2c1_;
  conc_m2c1dot_old_ = conc_m2c1dot_;

  StructPoroReaction::Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoroReactionECM::VisNames(std::map<std::string,int>& names)
{
  //call base class
  StructPoroReaction::VisNames(names);
  std::string name = "concentration_M2C1";
  names[name] = 1; // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool MAT::StructPoroReactionECM::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  //call base class
  if (StructPoroReaction::VisData(name,data,numgp,eleID))
    return true;
  else if (name=="concentration_M2C1")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    data[0] = conc_m2c1_;
    return true;
  }
  return false;
}


