/*----------------------------------------------------------------------*/
/*!
\file modpowerlaw.cpp

<pre>
Maintainer: Ursula Mayer
			mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <blitz/array.h>

#include "modpowerlaw.H"

extern struct _MATERIAL *mat;


MAT::ModPowerLaw::ModPowerLaw()
  : matdata_(NULL)
{
}


MAT::ModPowerLaw::ModPowerLaw(MATERIAL* matdata)
  : matdata_(matdata)
{
}


void MAT::ModPowerLaw::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;
  AddtoPack(data,matdata);
}


void MAT::ModPowerLaw::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


void MAT::ModPowerLaw::Evaluate(const blitz::Array<double,2>& 	velderiv,
                                double& 						visc_power)
{
	
  // get material parameters
  double m  	= matdata_->m.modpowerlaw->m_cons;      // consistency constant 
  double delta 	= matdata_->m.modpowerlaw->delta;       // safety factor
  double a      = matdata_->m.modpowerlaw->a_exp;      	// exponent
 

  // compute shear rate 
  double rateofshear = 0.0;
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());   // strain rate tensor
  epsilon = 0.5 * ( velderiv(i,j) + velderiv(j,i) );
  
  for(int rr=0;rr<3;rr++)
    for(int mm=0;mm<3;mm++)
    	rateofshear += epsilon(rr,mm)*epsilon(rr,mm);                 
 
  rateofshear = sqrt(2.0*rateofshear);
  
  // compute viscosity according to a modified power law model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  visc_power = m * pow((delta + rateofshear), (-1)*a);
}

#endif
