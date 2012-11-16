/*!----------------------------------------------------------------------
\file Itskov.cpp
\brief
This file contains the routines required for incompressible Itskov material law
with a penalty-function acocrding to Balzani et al.

Itskov, M.; Ehret, A. & Mavrilas, D.
A polyconvex anisotropic strain energy function for soft collagenous tissues
Biomechanics and Modeling in Mechanobiology, 2006, 5, 17-26

Ehret, A. & Itskov, M.
A polyconvex hyperelastic model for fiber-reinforced materials in application to soft tissues
Journal of Materials Science, 2007, 42, 8853-8863


The input line should read:
MAT 1 MAT_ITSKOV  MU_GS 0.0 MU_FIBERS 2.026255352 ALPHA 20.0 BETA 20.0 EPSILON 100.0 GAMMA 10.0 C 0.0 DENS 1.0E-6

C gives the desired deviation from full incompressibility. If the Value for C
is 0.0, the Penaltyparameters EPSILON and GAMMA are used as given in the input
file. Otherwise the Penaltyparameter GAMMA is used as given and EPSILON is
computed to fullfill the desired rate of incompressibility.
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "itskov.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Itskov::Itskov( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  alpha_(matdata->GetDouble("ALPHA")),
  beta_(matdata->GetDouble("BETA")),
  mu_fibers_(matdata->GetDouble("MU_FIBERS")),
  mu_GS_(matdata->GetDouble("MU_GS")),
  epsilon_(matdata->GetDouble("EPSILON")),
  gamma_(matdata->GetDouble("GAMMA")),
  comp_(matdata->GetDouble("C")),
  density_(matdata->GetDouble("DENS"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Itskov::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Itskov(this));
}

MAT::ItskovType MAT::ItskovType::instance_;


DRT::ParObject* MAT::ItskovType::Create( const std::vector<char> & data )
{
  MAT::Itskov* its = new MAT::Itskov();
  its->Unpack(data);
  return its;
}

/*---------------------------------------------------------------------*/
MAT::Itskov::Itskov()
  : params_(NULL)
{
  dserror("This material law - ITSKOV - is not maintained anymore.");
}


/*---------------------------------------------------------------------*/
MAT::Itskov::Itskov(MAT::PAR::Itskov* params)
  : params_(params)
{
  dserror("This material law - ITSKOV - is not maintained anymore.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Itskov::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Itskov::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

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
        params_ = static_cast<MAT::PAR::Itskov*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor (Itskov material law) |
 *----------------------------------------------------------------------*/
 void MAT::Itskov::Evaluate(
            const LINALG::Matrix<6,1>& glstrain,
            const int gp, const int ele_ID,
            DRT::Container& data_,
            const double time,
                  LINALG::Matrix<6,6>& cmat,
                  LINALG::Matrix<6,1>& stress)


/*void MAT::itskov::Evaluate(
            const LINALG::Matrix<6,1> *glstrain,
            const int gp, const int ele_ID, const double time,
                  LINALG::Matrix<6,6> *cmat,
                  LINALG::Matrix<6,1> *stress)*/

{
	int i,k;
	//double energy(0.0);
	double energy_constr(0.0), delta1(0.0), delta2(0.0);
	double epsilonPen_soll(0.0), W_soll(0.0);						//for adapted penalty parameter

	// get material parameters
 	double alpha = params_->alpha_;			//parameter
	double beta = params_->beta_;				//parameter
	double mu_fibers = params_->mu_fibers_;	//parameter scaling factor fibers
	double mu_GS = params_->mu_GS_;			//parameter scaling factor fibers
	double epsilonPen = params_->epsilon_;	//parameter for penalty function
	double gammaPen = params_->gamma_;		//parameter for penalty function
	double comp = params_->comp_;				//desired incompressibility for adapted penalty parameters


	//init container
	if ((time==0.0) && (gp==0)){
		std::vector<double>time_alt(8);
		data_.Add("time", time_alt);
		Epetra_SerialDenseMatrix I3_alt(8,2);
		data_.Add("I3", I3_alt);
		Epetra_SerialDenseMatrix WPen_alt(8,2);
		data_.Add("WPen", WPen_alt);
	}

	//Update container
	std::vector<double>* his_time=data_.GetMutable<std::vector<double> >("time");
	Epetra_SerialDenseMatrix* his_I3 = data_.GetMutable<Epetra_SerialDenseMatrix>("I3");
	Epetra_SerialDenseMatrix* his_WPen = data_.GetMutable<Epetra_SerialDenseMatrix>("WPen");
	if((*his_time)[gp]!=time)
	{
		(*his_I3)(gp,1)=(*his_I3)(gp,0);
		(*his_WPen)(gp,1)=(*his_WPen)(gp,0);
	}

//-----------------------------------------------------------------------------------------------
	if(comp!=0)		//calculation of adapted parameters for the penalty function
	{
		if((*his_WPen)(gp,1)==0)
			epsilonPen_soll=1;
		else
		{
			//Value for the energy kreated ba the penalty function so that |I3-1| equals the desired comp
			W_soll=(*his_WPen)(gp,1)*comp/(fabs((*his_I3)(gp,1)-1));

			//parameter epsilonPen, for that W_soll is achieved
			if((*his_I3)(gp,1)>1)
			{   epsilonPen_soll=W_soll/(pow((1+comp),gammaPen)+pow((1+comp),-gammaPen)-2);		}
			else
			{	epsilonPen_soll=W_soll/(pow((1-comp),gammaPen)+pow((1-comp),-gammaPen)-2);}
		}

		if(epsilonPen_soll<1) epsilonPen=1;
		else
		{	epsilonPen=epsilonPen_soll;}
	}
//---------------------------------------------------------------------------------------------------

//Definition structural tensors Li (i=1,2,3) in x,y,z-direction
	LINALG::Matrix<3,3> L1(true);
	LINALG::Matrix<3,3> L2(true);
	LINALG::Matrix<3,3> L3(true);

	L1(0,0)=1.0;
	L2(1,1)=1.0;
	L3(2,2)=1.0;

	double w1=1.0/3.0;
	double w2=1.0/3.0;
	double w3=1-w1-w2;

//Definition structural tensor L (see Itskov-Paper Equation(37))
	LINALG::Matrix<3,3> L(true);	//false
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
			L(i,k)=w1*L1(i,k)+w2*L2(i,k)+w3*L3(i,k); } }

// Green-Lagrange Strain Tensor
	LINALG::Matrix<3,3> E(true);		//false
	E(0,0) = (glstrain)(0);
	E(1,1) = (glstrain)(1);
	E(2,2) = (glstrain)(2);
	E(0,1) = 0.5 * (glstrain)(3);  E(1,0) = 0.5 * (glstrain)(3);
	E(1,2) = 0.5 * (glstrain)(4);  E(2,1) = 0.5 * (glstrain)(4);
	E(0,2) = 0.5 * (glstrain)(5);  E(2,0) = 0.5 * (glstrain)(5);

// Right Cauchy-Green Tensor  C = 2 * E + I
	LINALG::Matrix<3,3> C(E);
	C.Scale(2.0);
	C(0,0) += 1.0;
	C(1,1) += 1.0;
	C(2,2) += 1.0;

	//zum testen
	/*C(0,0) = 1.410156;
	C(1,1) = 1.0;
	C(2,2) = 1.0;*/

// Principal Invariants I3 = det(C)
	const double I3 = C(0,0)*C(1,1)*C(2,2) + C(0,1)*C(1,2)*C(2,0)
                  + C(0,2)*C(1,0)*C(2,1) - (C(0,2)*C(1,1)*C(2,0)
                  + C(0,1)*C(1,0)*C(2,2) + C(0,0)*C(1,2)*C(2,1));

// Calculation of Inverse Cinv of C
	LINALG::Matrix<3,3> Cinv(C);
	Cinv.Invert();

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//------------------------------------------------------GroundSubstance
	double factors1[6];
	LINALG::Matrix<3,3> PK2_1(true);
	LINALG::Matrix<6,6> cmat_1(true);
	calc_factors (C, Cinv, L, 0, 0, mu_GS, factors1);
	calc_Itskov (C, Cinv, L, PK2_1, cmat_1, factors1);
//---------------------------------------------------------------Fibers
	double factors2[6];
	LINALG::Matrix<3,3> PK2_2(true);
	LINALG::Matrix<6,6> cmat_2(true);
	calc_factors (C, Cinv, L, alpha, beta, mu_fibers, factors2);
	calc_Itskov (C, Cinv, L, PK2_2, cmat_2, factors2);

//---------------------------------------------------------------------
// Energy
//--------------------------------constraint function incompressibility
	energy_constr= (epsilonPen*(pow(I3,gammaPen)+pow(I3,(-gammaPen))-2));
//-----------------------------------------------strain-energy-function
	//energy=energy_constr+factors1[0]+factors1[1]+factors2[0]+factors2[1];

//---------------------------------------------------------------------
// PK2 Stresses
//--------------------------------constraint function incompressibility
	LINALG::Matrix<3,3> PK2_constr(Cinv);
	PK2_constr.Scale(2.0*epsilonPen*gammaPen*(pow(I3,gammaPen)-pow(I3,(-gammaPen))));
//------------------------------------2nd Piola-Kirchhoff Stress tensor
	LINALG::Matrix<3,3> PK2(PK2_constr);
	PK2+=PK2_1;
	PK2+=PK2_2;
// Transfer PK2 tensor to stress vector
	(stress)(0) = PK2(0,0);
  	(stress)(1) = PK2(1,1);
  	(stress)(2) = PK2(2,2);
  	(stress)(3) = PK2(0,1);
  	(stress)(4) = PK2(1,2);
  	(stress)(5) = PK2(0,2);

//---------------------------------------------------------------------
// Elasticity Tensor
//--------------------------------------------------constraint function
	delta1=4.0*epsilonPen*pow(gammaPen,2)*(pow(I3,gammaPen)+pow(I3,-gammaPen));
	delta2=-4.0*epsilonPen*gammaPen*(pow(I3,gammaPen)-pow(I3,-gammaPen));

	LINALG::Matrix<6,6> cmat_constr(true); //false
	ElastSymTensorMultiply(cmat_constr,delta1,Cinv,Cinv,0.0);
	ElastSymTensor_o_Multiply(cmat_constr,delta2,Cinv,Cinv,1.0);
//----------------------------------------------------elasticity tensor
	cmat=cmat_constr;
	cmat+=cmat_1;
	cmat+=cmat_2;


//zum testen
if(gp==-1)
{
	printf("\ntime: %f, gp: %d \n", time, gp);
	/*printf("\ntime: %f, element: %f \n", time, ele_ID);  */
	printf("I3: %.10f   \n", I3);
	printf("ePenSoll: %f   ePen: %f\n", epsilonPen_soll, epsilonPen);

	/*//printf("history %f   %f   %f   %f   %f   %f   %f   %f   \n",(*his_time)[0],(*his_time)[1],(*his_time)[2],(*his_time)[3],(*his_time)[4],(*his_time)[5],(*his_time)[6],(*his_time)[7]);
	printf("history_I3[x][0]: %f   %f   %f   %f   %f   %f   %f   %f   \n",(*his_I3)(0,0),(*his_I3)(1,0),(*his_I3)(2,0),(*his_I3)(3,0),(*his_I3)(4,0),(*his_I3)(5,0),(*his_I3)(6,0),(*his_I3)(7,0));
	printf("history_I3[x][1]: %f   %f   %f   %f   %f   %f   %f   %f   \n",(*his_I3)(0,1),(*his_I3)(1,1),(*his_I3)(2,1),(*his_I3)(3,1),(*his_I3)(4,1),(*his_I3)(5,1),(*his_I3)(6,1),(*his_I3)(7,1));*/

	//printf("\n\n");

printf("\nRight cauchy green\n");
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
				printf("%.10f   ", C(i,k));}
				printf("\n");}

printf("\nSecond Piola Kirchhhoff\n");
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
				printf("%.10f   ", PK2(i,k));}
				printf("\n");}


printf("\nSecond Piola Kirchhoff_Constr\n");
	for (k=0; k<3; k++) {
		for (i=0; i<3; i++) {
				printf("%.10f   ", PK2_constr(i,k));}
				printf("\n");}

/*printf("\nCelast\n");
	for (k=0; k<6; k++) {
	  for (i=0; i<6; i++) {
			printf("%f  ", (cmat)(k,i));}
			printf("\n");}*/

//FILE *pDatei = fopen("test.txt", "at");
//	fprintf(pDatei, "C11:%f    C22:%f   C33:%f   SPK11:%f    SPK22:%f     Cmat:%f\n", C(0,0), C(1,1), C(2,2), PK2(0,0), PK2(1,1), cmat(0,0));
//	fclose(pDatei);
}
//------------------------------------------------------------------------------------
//Write I3 and energy of the penalty function to history-variables
	(*his_time)[gp]=time;
	(*his_I3)(gp,0)=I3;
	(*his_WPen)(gp,0)=energy_constr;
//-------------------------------------------------------------------------------------
	//exit(0);
  return;
}


//Calculation of the factors in front of the tensor-Products
void MAT::Itskov::calc_factors (	LINALG::Matrix<3,3>&C,
						LINALG::Matrix<3,3>&Cinv,
						LINALG::Matrix<3,3>&L,
						double alpha,
						double beta,
						double mu,
						double *factors)
{
	double Ir, Kr;
//-----modified invariants
//------------------------------------------------1. Invariant: tr(CL)
	LINALG::Matrix<3,3> CL(true); //false
	CL.Multiply(1.0,C,L,0.0);
   	Ir= CL(0,0)+CL(1,1)+CL(2,2);
//-----------------------------2. Invariant incompressible: tr[Cinv*L]
	LINALG::Matrix<3,3> CinvL(true);	//false
	CinvL.Multiply(1.0,Cinv,L,0.0);
	Kr=(CinvL(0,0)+CinvL(1,1)+CinvL(2,2));
//---------------------------------------------------------------------
if(alpha==0) factors[0]=mu/4.0*(Ir-1);
else factors[0]=mu/4.0*(1.0/alpha*(exp(alpha*(Ir-1.0))-1.0)); 	//f

if(beta==0) factors[1]=mu/4.0*(Kr-1);
else factors[1]=mu/4.0*(1.0/beta*(exp(beta*(Kr-1.0))-1.0));		//g

factors[2]=mu/4.0*(exp(alpha*(Ir-1.0)));						//fI
factors[3]=mu/4.0*(exp(beta*(Kr-1.0)));							//gI
factors[4]=mu/4.0*(alpha*exp(alpha*(Ir-1.0)));					//fII
factors[5]=mu/4.0*(beta*exp(beta*(Kr-1.0)));					//gII

} //end of calc_factors


//Berechnung von Spannung und elasticity tensor von einer schicht mit entsprecheneden Parametern
void MAT::Itskov::calc_Itskov (	LINALG::Matrix<3,3>&C,
						LINALG::Matrix<3,3>&Cinv,
						LINALG::Matrix<3,3>&L,
						LINALG::Matrix<3,3>&PK2_,
						LINALG::Matrix<6,6>&cmat_,
						double *factors)
{
//--------------------------------------------Tensorproduct Cinv*L*Cinv
	LINALG::Matrix<3,3> CinvL(true);  //false
	CinvL.Multiply(1.0,Cinv,L,0.0);

	LINALG::Matrix<3,3> CinvLCinv(true);	//false
	CinvLCinv.Multiply(1.0,CinvL,Cinv,0.0);
//---------------------------------------------------------------------
//-------------------------------------------2nd Piola-Kirchhoff Stress
//PK2=2*(fI*L-gI*CinvLCinv)
//---------------------------------------------------------------------
	PK2_=L;
	PK2_.Scale(2.0*factors[2]);
	PK2_.Multiply(-2.0*factors[3],CinvL,Cinv,1.0);
//---------------------------------------------------------------------
//--------------------------------------------Tangent elasticity tensor
//cmat=4*(m/4*fII*LxL+m/4*gII*CinvLCinvxCinvLCinv+m/4*gI*(CinvoCinvLCinv+CinvLCinvoCinv))
//---------------------------------------------------------------------
	ElastSymTensorMultiply(cmat_,4*factors[4],L,L,0.0);
	ElastSymTensorMultiply(cmat_,4*factors[5],CinvLCinv,CinvLCinv,1.0);
	ElastSymTensor_o_Multiply(cmat_,4*factors[3],Cinv,CinvLCinv,1.0);
	ElastSymTensor_o_Multiply(cmat_,4*factors[3],CinvLCinv,Cinv,1.0);
return;
} //end of calc_Itskov

