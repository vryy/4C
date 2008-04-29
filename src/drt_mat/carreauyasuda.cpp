/*----------------------------------------------------------------------*/
/*!
\file carreauyasuda.cpp

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

#include "carreauyasuda.H"

extern struct _MATERIAL *mat;


MAT::CarreauYasuda::CarreauYasuda()
  : matdata_(NULL)
{
}


MAT::CarreauYasuda::CarreauYasuda(MATERIAL* matdata)
  : matdata_(matdata)
{
}


void MAT::CarreauYasuda::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;
  AddtoPack(data,matdata);
}


void MAT::CarreauYasuda::Unpack(const vector<char>& data)
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


void MAT::CarreauYasuda::Evaluate(const blitz::Array<double,2>& 	velderiv,
                                  double& 							visc_caryas)
{
	
  // get material parameters
  double mu_0 	= matdata_->m.carreauyasuda->mu_0;          // parameter for zero-shear viscosity
  double mu_inf = matdata_->m.carreauyasuda->mu_inf;      	// parameter for infinite-shear viscosity
  double lambda = matdata_->m.carreauyasuda->lambda;      	// parameter for characteristic time
  double a 		= matdata_->m.carreauyasuda->a_param;  		// constant parameter
  double b 		= matdata_->m.carreauyasuda->b_param;  		// constant parameter

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
  
  // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
  // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
  const double tmp = pow(lambda*rateofshear,b);
  visc_caryas = mu_inf + ((mu_0 - mu_inf)/pow((1 + tmp),a));
}

#endif
