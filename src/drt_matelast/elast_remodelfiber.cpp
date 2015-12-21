/*----------------------------------------------------------------------*/
/*!
\file elast_remodelfiber.cpp
\brief


the input line should read
  MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0

<pre>
Maintainer: Fabian Bräu
            braeu@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_remodelfiber.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::RemodelFiber::RemodelFiber(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  tdecay_(matdata->GetDouble("TDECAY")),
  sigmapre_(matdata->GetDouble("SIGMAPRE")),
  k_growth_(matdata->GetDouble("GROWTHFAC")),
  cur_w_collagen_(matdata->GetMutable<std::vector<double> >("COLMASSFRAC"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // check decay time validity
  if (tdecay_<=0.)
    dserror("decay time must be positive");
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         09/15 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::RemodelFiber::RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params)
  : params_(params),
    potsumfiber_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumfiber_.push_back(sum);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::PackSummand(DRT::PackBuffer& data) const
{
  int num_fiber = 0;
  num_fiber = last_ilambda_r_.size();

  AddtoPack(data,num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    AddtoPack(data,last_ilambda_r_[i]);
    AddtoPack(data,current_ilambda_r_[i]);
    AddtoPack(data,current_w_collagen_[i]);
  }

  if (params_ != NULL) // summands are not accessible in postprocessing mode
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumfiber_.size(); ++p)
     potsumfiber_[p]->PackSummand(data);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::UnpackSummand(const std::vector<char>& data,
                                               std::vector<char>::size_type& position)
{
  int num_fiber=0;
  ExtractfromPack(position,data,num_fiber);

  last_ilambda_r_.resize(num_fiber);
  current_ilambda_r_.resize(num_fiber);
  current_w_collagen_.resize(num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    ExtractfromPack(position,data,last_ilambda_r_[i]);
    ExtractfromPack(position,data,current_ilambda_r_[i]);
    ExtractfromPack(position,data,current_w_collagen_[i]);
  }

  // loop map of associated potential summands
  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
    potsumfiber_[p]->UnpackSummand(data,position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // setup fiber and inelastic history variable
  last_ilambda_r_.resize(potsumfiber_.size());
  current_ilambda_r_.resize(potsumfiber_.size());
  current_w_collagen_.resize(potsumfiber_.size());

  for(unsigned p=0;p<potsumfiber_.size();p++)
  {
    potsumfiber_[p]->Setup(linedef);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Update()
{
  dserror("Not implemented so far");

  return;
}


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::VisNames(std::map<std::string,int>& names)
{
  std::string inelastic_defgrd = "inelastic_defgrd_fiber";
  std::string result_inelastic_defgrad;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << inelastic_defgrd <<"_" << p;
    result_inelastic_defgrad = sstm.str();

    names[result_inelastic_defgrad] = 1;
  }


  std::string fiberdirection = "fiberdirection";
  std::string result_fiber;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << fiberdirection <<"_" << p;
    result_fiber = sstm.str();

    names[result_fiber] = 3;
  }


  std::string mass_fraction_col = "mass_fraction_col";
  std::string result_mass_fraction_col;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << mass_fraction_col <<"_" << p;
    result_mass_fraction_col = sstm.str();

    names[result_mass_fraction_col] = 1;
  }
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::ELASTIC::RemodelFiber::VisData(
  const std::string& name,
  std::vector<double>& data,
  int numgp,
  int eleID
)
{
  if (name == "inelastic_defgrd_fiber_0")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[0].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[0][gp];
    }
    data[0] = data[0]/last_ilambda_r_[0].size();

    return true;
  }
  if (name == "inelastic_defgrd_fiber_1")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[1].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[1][gp];
    }
    data[0] = data[0]/last_ilambda_r_[1].size();

    return true;
  }
  if (name == "inelastic_defgrd_fiber_2")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[2].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[2][gp];
    }
    data[0] = data[0]/last_ilambda_r_[2].size();

    return true;
  }



std::vector<LINALG::Matrix<3,1> > fiberdirection;
if(name == "fiberdirection_0")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[0]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}
if(name == "fiberdirection_1")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[1]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}
if(name == "fiberdirection_2")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[2]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}


if(name == "mass_fraction_col_0")
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<current_w_collagen_[0].size();++i)
  {
    data[0] += current_w_collagen_[0][i];
  }
  data[0] = data[0]/current_w_collagen_[0].size();


  return true;
}
if(name == "mass_fraction_col_1")
{
  if (data.size()!= 1) dserror("size mismatch");
  data[0] = current_w_collagen_[1][0];

  return true;
}
if(name == "mass_fraction_col_2")
{
  if (data.size()!= 1) dserror("size mismatch");
  data[0] = current_w_collagen_[2][0];

  return true;
}



dserror("The output is only implemented for three different fiber directions!!!");
return false;
}  // VisData()



