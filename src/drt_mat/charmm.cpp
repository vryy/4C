/*!----------------------------------------------------------------------
\file charmm.cpp
\brief CHARMm Interface to compute the mechanical properties of integrins

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
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "sys/types.h"
#include "sys/stat.h"
#include <math.h>
#include <string>
#include "charmm.H"
#include <time.h>
#include "../drt_so3/so_hex8.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::CHARMM::CHARMM(Teuchos::RCP<MAT::PAR::Material> matdata)
: Parameter(matdata),
origin_(matdata->GetInt("ORIGIN")),
fcl_(matdata->GetDouble("FCL")),
fcd_type_(matdata->Get<string>("FCD_TYPE")),
fcd_(matdata->Get<vector<double> >("FCD")),
fcds_(matdata->Get<vector<double> >("FCD_Space")),
scl_(matdata->GetDouble("SCL")),
scd_type_(matdata->Get<string>("SCD_TYPE")),
scd_(matdata->Get<vector<double> >("SCD")),
scds_(matdata->Get<vector<double> >("SCD_Space")),
fcdacc_(matdata->GetInt("FCD_Acceleration")),
atomicmass_(matdata->GetDouble("AtomicMass")),
facc_scale_(matdata->GetDouble("Facc_Scale")),
timeakma_(matdata->GetDouble("Time_AKMA")),
timescale_(matdata->GetDouble("Time_Scale")),
hard_(matdata->GetInt("HARD")),
c_scale_(matdata->GetDouble("c_Scale")),
path_(matdata->Get<string>("PATH")),
use_old_results_(matdata->GetInt("USE_OLD_RESULTS")),
serpar_(matdata->Get<string>("SERPAR")),
charmm_(matdata->Get<string>("CHARMM")),
input_(matdata->Get<string>("INPUT")),
nue_(matdata->GetDouble("NUE")),
density_(matdata->GetDouble("DENS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::CHARMM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CHARMM(this));
}


MAT::CHARMMType MAT::CHARMMType::instance_;


DRT::ParObject* MAT::CHARMMType::Create( const std::vector<char> & data )
{
  MAT::CHARMM* charmm = new MAT::CHARMM();
  charmm->Unpack(data);
  return charmm;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::CHARMM::CHARMM()
: params_(NULL) {
}

MAT::CHARMM::CHARMM(MAT::PAR::CHARMM* params)
: params_(params) {
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

    // pack type of this instance of ParObject
    int type = UniqueParObjectId();
    AddtoPack(data, type);

    // matid
    int matid = -1;
    if (params_ != NULL) matid = params_->Id(); // in case we are in post-process mode
    AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Unpack(const vector<char>& data) {
    vector<char>::size_type position = 0;
    // extract type
    int type = 0;
    ExtractfromPack(position, data, type);
    if (type != UniqueParObjectId()) dserror("wrong instance type data");

    // matid
    int matid;
    ExtractfromPack(position, data, matid);
    params_ = NULL;
    if (DRT::Problem::Instance()->Materials() != Teuchos::null)
      if (DRT::Problem::Instance()->Materials()->Num() != 0)
      {
        const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
        MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
        if (mat->Type() == MaterialType())
          params_ = static_cast<MAT::PAR::CHARMM*>(mat);
        else
          dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
      }

    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
//! Setup CHARMm history variables
// Actual and history variables, which need to be stored
// The updated version will be written in every iteration step.
// his_charmm = (    updated time, lasttime,
//                updated lambda(1), lambda(1)(t-dt),
//                updated lambda(2), lambda(2)(t-dt),
//                updated lambda(3), lambda(3)(t-dt),
//                I1,  I1(t-dt))
// his_mat[0] = c1 Neohookean from CHARMM for complete element
/*---------------------------------------------------------------------*/
void MAT::CHARMM::Setup(DRT::Container& data_) {

    // The following needs to come from the parameter in the final version
    vector<string> strain_type;
    strain_type.push_back("principal");
    strain_type.push_back("vector");

    vector<double> his_charmm;
    his_charmm.push_back(0.0); // actual time
    his_charmm.push_back(0.0); // time at last timestep
    for (int i = 0; i < (int) strain_type.size(); i++) {
	his_charmm.push_back(1.0); // updated lambda(1)(t)
	his_charmm.push_back(1.0); // lambda(1)(t-dt)
	his_charmm.push_back(1.0); // updated lambda(2)(t)
	his_charmm.push_back(1.0); // lambda(2)(t-dt)
	his_charmm.push_back(1.0); // updated lambda(3)(t)
	his_charmm.push_back(1.0); // lambda(3)(t-dt)
	his_charmm.push_back(3.0); // updated I1(t)
	his_charmm.push_back(3.0); // I1(t-dt)
	his_charmm.push_back(0.0); // updated v(t)
	his_charmm.push_back(0.0); // v(t-dt)
    }
    data_.Add("his_charmm", his_charmm);

    // position 1: c1 for material computation
    // position 2: energy at lambda = 0
    // position 3: # of atoms at lambda = 0
    // position 4: volume at lambda = 0
    vector<double> his_mat(4);
    his_mat.push_back(0.0);
    his_mat.push_back(0.0);
    his_mat.push_back(0.0);
    his_mat.push_back(0.0);
    data_.Add("his_mat", his_mat); // material property from CHARMm

    return;
}

/*----------------------------------------------------------------------*/
//! Compute second PK and constitutive tensor
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Evaluate(const LINALG::Matrix<NUM_STRESS_3D, 1 > * glstrain,
	LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,
	LINALG::Matrix<NUM_STRESS_3D, 1 > * stress,
	const int ele_ID,
	const int gp,
	DRT::Container& data_,
	const double time,
	const LINALG::SerialDenseMatrix& xrefe,
	const LINALG::SerialDenseMatrix& xcurr) {
#ifdef DEBUG
    if (!glstrain || !cmat || !stress)
	dserror("Data missing upon input in material CHARMm");
#endif

    //////////////////////////////////////////////////////////////////////////
    // Parameter collection
    // evaluate lamda at origin or at gp
    const bool origin = Origin(); // change only of xref and xcurr really working!!!!
    // How to treat the energy difference from MD
    // energy = store : store energy in element
    // energy = diff : use the difference between current strain and last strain
    const string energy = "store";
    // length of the protein in the main pulling direction [A]
    vector<double> characteristic_length(2);
    // Integrin length !!!!
    characteristic_length[0] = FCL(); //40.625; //50; //originally 44
    // Collagen length !!!!
    characteristic_length[1] = SCL();
    // characteristic direction of the protein
    // Possible selcetions:
    // principal = main strain direction (biggest eigenvalue)
    // vector = using the given vector
    // none = don't use the direction
    vector<string> strain_type;
    strain_type.push_back(FCDType());
    strain_type.push_back(SCDType());
    vector<LINALG::SerialDenseVector> d;
    LINALG::SerialDenseVector d_1(3);
    LINALG::SerialDenseVector d_2(3);
    d_1(0) = FCD()[0];
    d_1(1) = FCD()[1];
    d_1(2) = FCD()[2];
    d_2(0) = SCD()[0];
    d_2(1) = SCD()[1];
    d_2(2) = SCD()[2];
    d.push_back(d_1);
    d.push_back(d_2);
    // Add the directional space in case of principal direction
    vector<LINALG::SerialDenseVector> ds;
    LINALG::SerialDenseVector ds_1(3);
    LINALG::SerialDenseVector ds_2(3);
    ds_1(0) = FCDS()[0];
    ds_1(1) = FCDS()[1];
    ds_1(2) = FCDS()[2];
    ds_2(0) = SCDS()[0];
    ds_2(1) = SCDS()[1];
    ds_2(2) = SCDS()[2];
    ds.push_back(ds_1);
    ds.push_back(ds_2);
    // Use FCD to compute the acceleration in that direction to compute the
    // pulling force in CHARMm
    const bool FCD_Acceleration = FCDAcc();
    const double atomic_mass = AtomicMass(); // 18.2354; // amu; water
    const double Facc_scale = FaccScale(); // 1E26;
    // The basic coupling over the timescales; must be scaled to AKMA
    const bool movetime = true;
    const double time_to_AKMA = TimeAKMA();
    const double time_scale = TimeScale();
    // Use the hard coded charmm results (charmmfakeapi == true) or call charmm really (charmmfakeapi == false)
    const bool charmmhard = HARD();
    // Scale factor (by default c_CHARMm will be in N/m^2. This should be revised)
    const double c_scale = CScale(); //1E-9;
    //////////////////////////////////////////////////////////////////////////

    // Identity Matrix
    LINALG::Matrix < 3, 3 > I(true);
    for (int i = 0; i < 3; ++i) I(i, i) = 1.0;

    // Green-Lagrange Strain Tensor
    LINALG::Matrix < 3, 3 > E(false);
    E(0, 0) = (*glstrain)(0);
    E(1, 1) = (*glstrain)(1);
    E(2, 2) = (*glstrain)(2);
    E(0, 1) = 0.5 * (*glstrain)(3);
    E(1, 0) = 0.5 * (*glstrain)(3);
    E(1, 2) = 0.5 * (*glstrain)(4);
    E(2, 1) = 0.5 * (*glstrain)(4);
    E(0, 2) = 0.5 * (*glstrain)(5);
    E(2, 0) = 0.5 * (*glstrain)(5);

    // Right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix < 3, 3 > C(E);
    C.Scale(2.0);
    C += I;

    // Principal Invariants I1 = tr(C) and I3 = det(C)
    double I1 = C(0, 0) + C(1, 1) + C(2, 2);
    double I3 = C(0, 0) * C(1, 1) * C(2, 2) + C(0, 1) * C(1, 2) * C(2, 0)
	    + C(0, 2) * C(1, 0) * C(2, 1) - (C(0, 2) * C(1, 1) * C(2, 0)
	    + C(0, 1) * C(1, 0) * C(2, 2) + C(0, 0) * C(1, 2) * C(2, 1));

    // Calculation of C^-1 (Cinv)
    LINALG::Matrix < 3, 3 > Cinv(C);
    Cinv.Invert();

    /////////////////////////////////////////////////////////////////////CHARMm
    // CHARMm things come here
    if (gp == 0) {

	// Get the strains in the characteristic directions
	LINALG::Matrix < 3, 3 > V(C);
	LINALG::SerialDenseVector lambda(3);
	vector<LINALG::SerialDenseVector> dir_lambdas;
	vector<LINALG::Matrix < 3, 3 > > dir_eigenv;
	// go through number of directions
	for (int i = 0; i < (int) strain_type.size(); i++) {
	    if (strain_type[i].compare("principal") == 0) {
			V.SetCopy(C);
			EvalStrain(origin, xrefe, xcurr, V, lambda);
			// flip the unit vector in case it's showing not in the right direction.
			for (int j = 0; j < 3; j++) {
		    	if (ds[i](j) != 0 && ((ds[i](j) < 0 && V(j, 2) > 0) || (ds[i](j) > 0 && V(j, 2) < 0))) {
					V(0, 2) *= -1;
					V(1, 2) *= -1;
					V(2, 2) *= -1;
		    	}
			}
			dir_lambdas.push_back(lambda);
			dir_eigenv.push_back(V);
	    } else if (strain_type[i].compare("vector") == 0) {
			V.SetCopy(C);
			for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) V(k, l) = d[i](k) * V(k, l) * d[i](l);
			EvalStrain(origin, xrefe, xcurr, V, lambda);
			dir_lambdas.push_back(lambda);
			dir_eigenv.push_back(V);
	    } else if (strain_type[i].compare("none") == 0) {
			V.Clear();
			lambda.Zero();
			dir_lambdas.push_back(lambda);
			dir_eigenv.push_back(V);
	    } else {
			dserror("No valid strain type given for CHARMm!");
	    }
	}
	//cout << dir_lambdas[1](0) << " : " << dir_lambdas[1](1) << " : " << dir_lambdas[1](2) << endl;
	//cout << dir_eigenv[0] << endl;


	// Update and reconfigure history
	vector<double>* his;
	his = data_.GetMutable<vector<double> >("his_charmm");
	if ((*his)[0] < time) {
	    (*his)[1] = (*his)[0]; // time
	    for (int i = 0; i < (int) strain_type.size(); i++) {
		(*his)[3 + (i * 10)] = (*his)[2 + (i * 10)]; // lambda(0)
		(*his)[5 + (i * 10)] = (*his)[4 + (i * 10)]; // lambda(1)
		(*his)[7 + (i * 10)] = (*his)[6 + (i * 10)]; // lambda(2)
		(*his)[9 + (i * 10)] = (*his)[8 + (i * 10)]; // I1
		(*his)[11 + (i * 10)] = (*his)[10 + (i * 10)]; // v
	    }
	}
	(*his)[0] = time;
	for (int i = 0; i < (int) strain_type.size(); i++) {
	    (*his)[2 + (i * 10)] = dir_lambdas[i](0);
	    (*his)[4 + (i * 10)] = dir_lambdas[i](1);
	    (*his)[6 + (i * 10)] = dir_lambdas[i](2);
	    (*his)[8 + (i * 10)] = I1;
	}
	//cout  << (*his)[0] << " : " << (*his)[1] << " : ";
	//for(int i=0;i<(int)strain_type.size();i++) {
	//    cout    << (*his)[2+(i*8)] << " : " << (*his)[3+(i*8)] << " : "
	//            << (*his)[4+(i*8)] << " : " << (*his)[5+(i*8)] << " : "
	//            << (*his)[6+(i*8)] << " : " << (*his)[7+(i*8)] << " : "
	//            << (*his)[8+(i*8)] << " : " << (*his)[9+(i*8)] << " : "
	//	    << (*his)[10+(i*8)] << " : " << (*his)[11+(i*8)] << " : ";
	//}
	//cout <<  endl;

	// Prepare and call CHARMm in its beauty itself
	// get lambda t-dt information
	vector<double> lambda_his;
	for (int i = 0; i < (int) strain_type.size(); i++) {
	    lambda_his.push_back((*his)[7 + (i * 10)]);
	}

	// Data preparation for CHARMm
	// First charateristic direction (FCD)
	// calculate STARTD and ENDD for CHARMm (integrin)
	double FCD_STARTD = characteristic_length[0] * (lambda_his[0] - 1);
	double FCD_ENDD = characteristic_length[0] * (dir_lambdas[0](2) - 1); // Check for better way to choose!!!!
	// get direction for FCD (integrin)
	LINALG::SerialDenseVector FCD_direction(3);
	FCD_direction(0) = dir_eigenv[0](0, 2);
	FCD_direction(1) = dir_eigenv[0](1, 2);
	FCD_direction(2) = dir_eigenv[0](2, 2);
	//printf("%f",FCD_ENDD);
	//cout << "FCD: " << time << " STARTD: " << FCD_STARTD << " ENDD: " << FCD_ENDD << endl;

	// Compute the acceleration in FCD direction
	double FCD_v;
	double FCD_a;
	double FCD_Force;
	if (FCD_Acceleration) {
	    EvalAccForce(FCD_STARTD,FCD_ENDD,(*his)[1],time,(*his)[11],atomic_mass,Facc_scale,FCD_v,FCD_a,FCD_Force);
	    (*his)[10] = FCD_v;
	}

	// Second charateristic direction (SCD)
	// calculate STARTD and ENDD for CHARMm (collagen)
	double SCD_STARTD = characteristic_length[1] * (1 - lambda_his[1]);
	double SCD_ENDD = characteristic_length[1] * (1 - dir_lambdas[1](2));
	// get direction for SCD (collagen)
	LINALG::SerialDenseVector SCD_direction(3);
	SCD_direction(0) = dir_eigenv[1](0, 2);
	SCD_direction(1) = dir_eigenv[1](1, 2);
	SCD_direction(2) = dir_eigenv[1](2, 2);
	//cout << "SCD: " << "STARTD: " << SCD_STARTD << " ENDD: " << SCD_ENDD << SCD_direction << endl;

	// Compute the appropriate timestep for the moving time in MD
	double time_move = 0.0;
	if (movetime) time_move = time_scale * time_to_AKMA * (time - (*his)[1]);
	//cout << "Time: " << time_move << endl;

	// Check if results actually can be computed by CHARMm
	//if (STARTD != ENDD) dserror("STARTD and ENDD identical! CHARMm will not produce any results.");

	// Call API to CHARMM
	// Results vector: charmm_result
	// (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
	LINALG::SerialDenseVector direction(3);
	LINALG::SerialDenseVector charmm_result(6);
	map<string, double> CHARMmPar;
	CHARMmPar["FCD_STARTD"] = FCD_STARTD;
	CHARMmPar["FCD_ENDD"] = FCD_ENDD;
	CHARMmPar["FCD_dir_x"] = FCD_direction(0);
	CHARMmPar["FCD_dir_y"] = FCD_direction(1);
	CHARMmPar["FCD_dir_z"] = FCD_direction(2);
	CHARMmPar["FCD_Force"] = FCD_Force;
	CHARMmPar["SCD_STARTD"] = SCD_STARTD;
	CHARMmPar["SCD_ENDD"] = SCD_ENDD;
	CHARMmPar["SCD_dir_x"] = SCD_direction(0);
	CHARMmPar["SCD_dir_y"] = SCD_direction(1);
	CHARMmPar["SCD_dir_z"] = SCD_direction(2);
	CHARMmPar["movetime"] = time_move;
	if (charmmhard) {
	    // Just give the starting and ending strain in hard coded case
	    CHARMmfakeapi(FCD_STARTD, FCD_ENDD, charmm_result);
	} else {
	    if (FCD_STARTD != FCD_ENDD) {
		CHARMmfileapi(CHARMmPar, charmm_result);
	    }
	}

        vector<double>* his_mat;
        his_mat = data_.GetMutable<vector<double> >("his_mat");
        double I1_lastt = (*his)[9];
        if (FCD_STARTD == 0.0) (*his_mat)[1] = charmm_result[0];
        if (energy.compare("diff") == 0) {
            // Calculate new c (Neo-Hooke) parameter
            // c = E_FE / (I1(t_n) - I1(t_n-1)) [N/m^2]
            // E_FE = (E_MD(t_n) - E_MD(t_n-1)) * 1000 * 4.1868 * ( #Atoms / N_a )
            double E_MD = charmm_result[1] - charmm_result[0]; // kcal/mole
            double Volume = 0.5 * (charmm_result[4] + charmm_result[5]) * 1E-30; // A^3 *  (10^-10)^3
            double noAtoms = 0.5 * (charmm_result[3] + charmm_result[2]);

            double c;
            if (I1 != I1_lastt) {
                c = (1 / (I1 - I1_lastt)) * (1 / Volume) * E_MD * 1000 * 4.1868 * (noAtoms / 6.02214E23);
                //c = (1 / (I1 - 3)) * (1 / Volume) * E_MD * 1000 * 4.1868 * (noAtoms / 6.02214E23);
            } else { c = 0.0; }

            if (isnan(c)) c = 0;
            if (isinf(c)) c = 0;
            c = c_scale * c;

            if (FCD_STARTD != FCD_ENDD) {
                if (I1_lastt == 3) (*his_mat)[0] = c;
                else (*his_mat)[0] = c * ((I1 - I1_lastt) / (I1 - 3));
            } else {
                (*his_mat)[0] = 0.0;
            }

        } else if (energy.compare("store") == 0) {
            // Calculate new c (Neo-Hooke) parameter
            // c = E_FE / (I1(t_n) - 3) [N/m^2]
            // E_FE = (E_MD(t_n) - E_MD(t_0)) * 1000 * 4.1868 * ( #Atoms / N_a )

            if (FCD_STARTD == 0.0 || (*his_mat)[1] > charmm_result[0]) {
                (*his_mat)[1] = charmm_result[0]; // energy
                (*his_mat)[2] = charmm_result[2]; // noAtoms
                (*his_mat)[3] = charmm_result[4]; // volume
            }
            double E_MD = charmm_result[1] - (*his_mat)[1]; // kcal/mole
            double noAtoms = 0.5 * (charmm_result[3] + (*his_mat)[2]);
            double Volume = 0.5 * ( (*his_mat)[3] + charmm_result[5]) * 1E-30; // A^3 *  (10^-10)^3

            double c;
            if (I1 != I1_lastt) {
                c = (1 / (I1 - 3)) * (1 / Volume) * E_MD * 1000 * 4.1868 * (noAtoms / 6.02214E23);
            } else { c = 0.0; }

            if (isnan(c)) c = 0;
            if (isinf(c)) c = 0;
            c = c_scale * c;

            if (FCD_STARTD != FCD_ENDD) {
                (*his_mat)[0] = c;
            } else {
                (*his_mat)[0] = 0.0;
            }

        }
        if (FCD_STARTD != FCD_ENDD) cout << "ACE: c = " << (*his_mat)[0] << " I1(t_n-1) = " << I1_lastt << " I1(t_n) = " << I1<<  endl;

	//cout << "MD Result: " << charmm_result[0] << ":" << charmm_result[1] << " "
	//        << ( charmm_result[1] - charmm_result[0]) << " "
	//        << charmm_result[5] << endl;

    }

    //
    ///////////////////////////////////////////////////////////////////////////
    // Material Constants c1 and beta
    double c1(0.0);
    double nu = NUE(); // intermediate for testing purpose only
    double beta = nu / (1 - 2 * nu);
    if (time > 0.0) {
	vector<double>* his_mat = data_.GetMutable<vector<double> >("his_mat");
	if ((*his_mat)[0] != 0.0) c1 = (*his_mat)[0];
	//cout << time << " " << c1 << endl;
    } else {
	c1 = 1.0;
    }

    // Energy
    //double W = c1/beta * (pow(I3,-beta) - 1) + c1 * (I1-3);

    // PK2 Stresses
    LINALG::Matrix < 3, 3 > PK2(false);
    int i, j;
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++) {
	    PK2(i, j) = 2 * c1 * (I(i, j) - pow(I3, -beta) * Cinv(i, j));
	}

    // Transfer PK2 tensor to stress vector
    (*stress)(0) = PK2(0, 0);
    (*stress)(1) = PK2(1, 1);
    (*stress)(2) = PK2(2, 2);
    (*stress)(3) = PK2(0, 1);
    (*stress)(4) = PK2(1, 2);
    (*stress)(5) = PK2(0, 2);

    // Elasticity Tensor
    double delta6 = 4. * c1 * beta * pow(I3, -beta);
    double delta7 = 4. * c1 * pow(I3, -beta);

    int k, l;
    LINALG::Matrix < 9, 9 > ET(false);


    for (k = 0; k < 3; k++)
	for (l = 0; l < 3; l++) {
	    ET(k, l) = delta6 * (Cinv(0, 0) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(0, k) * Cinv(0, l) + Cinv(0, l) * Cinv(0, k));
	    ET(k + 3, l) = delta6 * (Cinv(1, 0) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(1, k) * Cinv(0, l) + Cinv(1, l) * Cinv(0, k));
	    ET(k + 3, l + 3) = delta6 * (Cinv(1, 1) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(1, k) * Cinv(1, l) + Cinv(1, l) * Cinv(1, k));
	    ET(k + 6, l) = delta6 * (Cinv(2, 0) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(2, k) * Cinv(0, l) + Cinv(2, l) * Cinv(0, k));
	    ET(k + 6, l + 3) = delta6 * (Cinv(2, 1) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(2, k) * Cinv(1, l) + Cinv(2, l) * Cinv(1, k));
	    ET(k + 6, l + 6) = delta6 * (Cinv(2, 2) * Cinv(k, l)) +
		    delta7 * 0.5 * (Cinv(2, k) * Cinv(2, l) + Cinv(2, l) * Cinv(2, k));
	}

    (*cmat)(0, 0) = ET(0, 0);
    (*cmat)(0, 1) = ET(1, 1);
    (*cmat)(0, 2) = ET(2, 2);
    (*cmat)(0, 3) = ET(1, 0);
    (*cmat)(0, 4) = ET(2, 1);
    (*cmat)(0, 5) = ET(2, 0);

    (*cmat)(1, 0) = ET(3, 3);
    (*cmat)(1, 1) = ET(4, 4);
    (*cmat)(1, 2) = ET(5, 5);
    (*cmat)(1, 3) = ET(4, 3);
    (*cmat)(1, 4) = ET(5, 4);
    (*cmat)(1, 5) = ET(5, 3);

    (*cmat)(2, 0) = ET(6, 6);
    (*cmat)(2, 1) = ET(7, 7);
    (*cmat)(2, 2) = ET(8, 8);
    (*cmat)(2, 3) = ET(7, 6);
    (*cmat)(2, 4) = ET(8, 7);
    (*cmat)(2, 5) = ET(8, 6);

    (*cmat)(3, 0) = ET(3, 0);
    (*cmat)(3, 1) = ET(4, 1);
    (*cmat)(3, 2) = ET(5, 2);
    (*cmat)(3, 3) = ET(4, 0);
    (*cmat)(3, 4) = ET(5, 1);
    (*cmat)(3, 5) = ET(5, 0);

    (*cmat)(4, 0) = ET(6, 3);
    (*cmat)(4, 1) = ET(7, 4);
    (*cmat)(4, 2) = ET(8, 5);
    (*cmat)(4, 3) = ET(7, 3);
    (*cmat)(4, 4) = ET(8, 4);
    (*cmat)(4, 5) = ET(8, 3);

    (*cmat)(5, 0) = ET(6, 0);
    (*cmat)(5, 1) = ET(7, 1);
    (*cmat)(5, 2) = ET(8, 2);
    (*cmat)(5, 3) = ET(7, 0);
    (*cmat)(5, 4) = ET(8, 1);
    (*cmat)(5, 5) = ET(8, 0);

    return;
}


/*----------------------------------------------------------------------*/
//! Evaluate strains in the charateristic directions
/*----------------------------------------------------------------------*/
void MAT::CHARMM::EvalStrain(const bool& origin,
	const LINALG::SerialDenseMatrix& xrefe,
	const LINALG::SerialDenseMatrix& xcurr,
	LINALG::Matrix < 3, 3 > & C,
	LINALG::SerialDenseVector& lambda) {
    LINALG::SerialDenseVector lambda2(3);
    LINALG::SerialDenseMatrix Ctmp(3, 3);
    if (origin) {
	// vector of dN/dxsi |r=s=t=0.0
	double dN0_vector[24] ={-0.125, -0.125, -0.125,
	    +0.125, -0.125, -0.125,
	    +0.125, +0.125, -0.125,
	    -0.125, +0.125, -0.125,
	    -0.125, -0.125, +0.125,
	    +0.125, -0.125, +0.125,
	    +0.125, +0.125, +0.125,
	    -0.125, +0.125, +0.125};

	// shape function derivatives, evaluated at origin (r=s=t=0.0)
	Epetra_DataAccess CV = Copy;
	Epetra_SerialDenseMatrix dN0(CV, dN0_vector, 3, 3, 8);

	// compute Jacobian, evaluated at element origin (r=s=t=0.0)
	LINALG::SerialDenseMatrix invJacobian0(3, 3);
	invJacobian0.Multiply('N', 'N', 1.0, dN0, xrefe, 0.0);
	const double detJacobian0 = LINALG::NonsymInverse3x3(invJacobian0);
	if (detJacobian0 < 0.0) dserror("Jacobian at origin negativ (CHARMMAPI)");

	//cout << invJacobian0 << endl;
	LINALG::SerialDenseMatrix N_XYZ(3, 8);
	//compute derivatives N_XYZ at gp w.r.t. material coordinates
	// by N_XYZ = J^-1 * N_rst
	N_XYZ.Multiply('N', 'N', 1.0, invJacobian0, dN0, 0.0);
	// (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
	LINALG::SerialDenseMatrix defgrd0(3, 3);
	defgrd0.Multiply('T', 'T', 1.0, xcurr, N_XYZ, 0.0);
	// Right Cauchy-Green tensor = F^T * F
	LINALG::SerialDenseMatrix C0(3, 3);
	C0.Multiply('T', 'N', 1.0, defgrd0, defgrd0, 0.0);

	// compute current eigenvalues of gaussian point C
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Ctmp(i, j) = C0(i, j);
	LINALG::SymmetricEigen(Ctmp, lambda2, 'V', false);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) C(i, j) = Ctmp(i, j);
    } else {
	// compute current eigenvalues of gaussian point C
	LINALG::SerialDenseMatrix Ctmp(3, 3);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Ctmp(i, j) = C(i, j);
	LINALG::SymmetricEigen(Ctmp, lambda2, 'V', false);
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) C(i, j) = Ctmp(i, j);
    }
    for (int i = 0; i < 3; i++) lambda(i) = sqrt(lambda2(i));
}


/*----------------------------------------------------------------------*/
//! Compute acceleration and force in characteristic direction
/*----------------------------------------------------------------------*/
void MAT::CHARMM::EvalAccForce(
	const double& STARTD,
	const double& ENDD,
	const double& time_STARTD,
	const double& time_ENDD,
	const double& v_his,
	const double& atomic_mass,
	const double& Facc_scale,
	double& v,
	double& a,
	double& Force) {

    double amu_to_kg = 1.66053886E-27;
    double v_0 = v_his;
    v = 0.0;
    // Compute velocity
    if (time_ENDD != time_STARTD) v = abs(ENDD - STARTD) / (time_ENDD - time_STARTD);
    else v = 0.0;
    // Round v and v_0 off, otherwise comparison is not working
    int d = 5;
    double n = v;
    v = floor(n * pow(10., d) + .5) / pow(10., d);
    n = v_0;
    v_0 = floor(n * pow(10., d) + .5) / pow(10., d);
    if (v == v_his) v_0 = 0.0; // Switch between tangent and secant?? Sure to do that??
    // Compute acceleration and Force
    if (time_ENDD != time_STARTD) a = (v - v_0) / (time_ENDD - time_STARTD);
    else a = 0.0;
    Force = atomic_mass * amu_to_kg * a * Facc_scale;

    //cout    << "ACC: " << a << " " << Force << " "
	//    << STARTD << " " << ENDD << " "
	//    << time_STARTD << " " << time_ENDD << " "
	//    << v << " " << v_his << endl;

}

/*----------------------------------------------------------------------*/
//! File based API to CHARMM
/*----------------------------------------------------------------------*/
void MAT::CHARMM::CHARMmfileapi(
	map<string, double>& CHARMmPar,
	LINALG::SerialDenseVector& charmm_result) {

    FILE* tty;
    ios_base::fmtflags flags = cout.flags(); // Save original flags
    ostringstream output(ios_base::out);
    ////////////////////////////////////////////////////////////////////////////
    // Variables needed for CHARMM and getting the results
    // Decide if parallel or seriell
    const bool use_old_results = Use_old_Results();
    const string serpar = Serpar(); // ser = seriell; par = mpirun; pbs = PBS Torque
    const string charmm = CHARMMEXE();
    const string input = INPUT();
    // FC6 setup
    //const char* path = "/home/metzke/ccarat.dev/codedev/charmm.fe.codedev/";
    //const char* path = "/home/metzke/projects/water/";
    // Mac setup
    //const char* path = "/Users/rmetzke/research/baci.dev/codedev/charmm.fe.codedev/";
    //const char* path = "/Users/rmetzke/research/projects/water/";
    //const char* charmm = "/Users/rmetzke/bin/charmm";
    //const char* mpicharmm = "/Users/rmetzke/bin/mpicharmm";
    //char* input = "1dzi_fem.inp";
    const string mdnature = "thermal"; // cold = minimization; thermal = fully dynamic with thermal energy; pert = pertubation
    output << "output/ACEcold_" << CHARMmPar["FCD_STARTD"] << "_" << CHARMmPar["FCD_ENDD"] << ".out";
    ////////////////////////////////////////////////////////////////////////////

    // Assemble all file and path names first
    ostringstream statusfile(ios_base::out);
    statusfile << Path() << "output/status_" << CHARMmPar["FCD_ENDD"] << ".out";

    // Print out the beginning of the CHARMM info line
    cout << std::setw(4) << left << "MD (" << showpoint << CHARMmPar["FCD_STARTD"] << std::setw(2) << "->" << CHARMmPar["FCD_ENDD"] << std::setw(3) << "): " << flush;

    // Check if the status file already exists
    // In that case skip the charmm call
    struct stat statusFileInfo;
    map<string, double> md_status;
    md_status.clear();
    if (stat(statusfile.str().c_str() ,&statusFileInfo) == 0) Reader(statusfile, md_status);
    if ( (stat(statusfile.str().c_str(),&statusFileInfo) != 0 && !md_status["CHARMMEND"]) || !use_old_results) {
	// Assemble the command line for charmm
	ostringstream command(ios_base::out);
	if (serpar.compare("ser") == 0) command << "cd " << Path() << " && " << charmm;
	else if (serpar.compare("par") == 0) command << "cd " << Path() << " && " << "openmpirun -np 2 " << charmm;
	else dserror("What you want now? Parallel or not!");
	command << " FCDSTARTD=" << CHARMmPar["FCD_STARTD"] << " FCDENDD=" << CHARMmPar["FCD_ENDD"]
		<< " FCDX=" << CHARMmPar["FCD_dir_x"] << " FCDY=" << CHARMmPar["FCD_dir_y"] << " FCDZ=" << CHARMmPar["FCD_dir_z"]
		<< " FCDForce=" << CHARMmPar["FCD_Force"]
		<< " SCDSTARTD=" << CHARMmPar["SCD_STARTD"] << " SCDENDD=" << CHARMmPar["SCD_ENDD"]
		<< " SCDX=" << CHARMmPar["SCD_dir_x"] << " SCDY=" << CHARMmPar["SCD_dir_y"] << " SCDZ=" << CHARMmPar["SCD_dir_z"]
		<< " MTIME=" << CHARMmPar["movetime"];
	if (serpar.compare("ser") == 0)
	    command << " < " << input << " > " << output.str();
	else if (serpar.compare("par") == 0) {
	    command << " INPUTFILE=" << input
		    << " < " << "stream.inp" << " > " << output.str();
	} else dserror("What you want now? Parallel or not!");
	cout << "0|" << flush;
	if ((tty = popen(command.str().c_str(), "r")) == NULL) dserror("CHARMM can not be started!");
	int runresult = pclose(tty);
	cout << runresult << "|";
    } else {
	cout << "-1|-1|" << flush;
    }


    // Read the results
    if (mdnature.compare("cold") == 0) {
        dserror("Cold currently not supported");
    } else if (mdnature.compare("thermal") == 0) {
	Readresults(CHARMmPar, charmm_result);
    } else {
	dserror("No included MD Simulation technique given!.");
    }
    //cout << endl;
    cout.flags(flags); // Set the flags to the way they were
}

//*----------------------------------------------------------------------*/
//! Read results from thermal CHARMm results files
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Readresults(map<string, double>& CHARMmPar,
	LINALG::SerialDenseVector& charmm_result) {

    // Check the status of CHARMm and if it run through
    ostringstream status(ios_base::out);
    status << Path() << "output/status_" << CHARMmPar["FCD_ENDD"] << ".out";
    map<string, double> md_status;
    md_status.clear();
    Reader(status, md_status);
    if (!md_status["CHARMMEND"]) dserror("CHARMm API: Run not successful till end.");
    else cout << std::setw(5) << left << "0" << flush;

    // Get the results at STARTD + ENDD
    ostringstream results(ios_base::out);
    map<string, double> md_STARTD;
    md_STARTD.clear();
    results << Path() << "output/results_" << CHARMmPar["FCD_STARTD"] << ".out";
    Reader(results, md_STARTD);
    results.str("");
    map<string, double> md_ENDD;
    md_ENDD.clear();
    results << Path() << "output/results_" << CHARMmPar["FCD_ENDD"] << ".out";
    Reader(results, md_ENDD);

    // Print out Information line
    cout << std::setw(4) << "dV:" << std::setw(15) << left << std::scientific << std::setprecision(6) << (md_ENDD["E"] - md_STARTD["E"]);
    cout << std::setw(8) << "#Atoms:" << std::setw(10) << left << fixed << std::setprecision(0) << md_ENDD["NUM"];
    cout << std::setw(8) << "Volume:" << std::setw(12) << left << std::setprecision(3) << md_ENDD["VOL"] << endl;

    ////////////////////////////////////////////////////////////////////////////
    // Results vector: charmm_result
    // (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
    charmm_result[0] = md_STARTD["E"];
    charmm_result[1] = md_ENDD["E"];
    charmm_result[2] = md_STARTD["NUM"];
    charmm_result[3] = md_ENDD["NUM"];
    charmm_result[4] = md_STARTD["VOL"];
    charmm_result[5] = md_ENDD["VOL"];
    ////////////////////////////////////////////////////////////////////////////

}


/*----------------------------------------------------------------------*/
//! Read status and results files and make results map
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Reader(const ostringstream& file,
	map<string, double>& content) {

    ifstream filestream(file.str().c_str());
    string line;
    if (filestream.is_open()) {
	while (!filestream.eof()) {
	    getline(filestream,line);
	    if (line.compare(0,1,"R") == 0) {
		vector<string> token;
		string buf;
		stringstream linestream(line);
		while (linestream >> buf)
		    token.push_back(buf);
		content[token[1].c_str()] = atof(token[2].c_str());
		token.clear();
	    }
	}
	filestream.close();
    } else dserror("CHARMM API READER: File cannot be opened: ");


    //cout << endl << "File Content: " << endl;
    //for( map<string, double>::iterator ii=content.begin(); ii!=content.end(); ++ii) {
    //	cout << (*ii).first << ": " << (*ii).second << endl;
    //}
    //cout << endl;

}


/*----------------------------------------------------------------------*/
//! Hard coupling without calling CHARMm
/*----------------------------------------------------------------------*/
void MAT::CHARMM::CHARMmfakeapi(const double STARTD,
	const double ENDD,
	LINALG::SerialDenseVector& charmm_result) {
    // Define the number n of steps / results from CHARMm (or any MD simluation)
    // If n=2 then it is assumes that always the same values should be used for all steps.
    const int n = 2;
    // Define roundoff for choosing in which step we are
    const double roundoff = 0.005;
    // Hard coded data structure (second variable):
    // (STARTD, Energy, # of Atoms, Volume)
    LINALG::SerialDenseMatrix MD(n, 4);

    ////////////////////////////////////////////////////////////////////////////
    // Hard coded results from MD
    //MD(0, 0) = 0.0;
    //MD(0, 1) = -330.912;
    //MD(0, 2) = 1202;
    //MD(0, 3) = 9954.29;

    //MD(1, 0) = -0.8125;
    //MD(1, 1) = -321.671;
    //MD(1, 2) = 1141;
    //MD(1, 3) = 9441.08;
    MD(0, 0) = 0.0;
    MD(0, 1) = -1.8;
    MD(0, 2) = 6;
    MD(0, 3) = 45.4763;

    MD(1, 0) = 0.0;
    MD(1, 1) = -0.0;
    MD(1, 2) = 6;
    MD(1, 3) = 46.3414;
    //
    ////////////////////////////////////////////////////////////////////////////

    // Compute the charmm_result vector
    // (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
    ios_base::fmtflags flags = cout.flags(); // Save original flags

    for (int i = n - 1; i >= 0; i--) {
	//cout << ENDD << " " << MD(i,0);
	if (abs(ENDD) == 0.0) { // start call at the beginning; just to give some information
	    i = 0;
	    cout << std::setw(4) << left << "MD (" << showpoint << STARTD << std::setw(2) << "->" << ENDD << std::setw(3) << "): " << flush;
	    charmm_result[0] = NAN;
	    charmm_result[1] = MD(i, 1);
	    charmm_result[2] = NAN;
	    charmm_result[3] = MD(i, 2);
	    charmm_result[4] = NAN;
	    charmm_result[5] = MD(i, 3);
	    cout << std::setw(4) << "V(0):" << std::setw(15) << left << std::scientific << std::setprecision(6) << (charmm_result[1]);
	    cout << std::setw(8) << "#Atoms:" << std::setw(10) << left << fixed << std::setprecision(0) << charmm_result[3] << std::setw(8) << "Volume:" << std::setw(12) << left << std::setprecision(2) << charmm_result[5] << endl;
	    i = -1; //break loop
	} else if (abs(ENDD) < (abs(MD(i, 0)) + roundoff) && abs(ENDD) > (abs(MD(i, 0)) - roundoff)) {
	    // main loop where basically at every step the data is given
	    cout << std::setw(4) << left << "MD (" << showpoint << STARTD << std::setw(2) << "->" << ENDD << std::setw(3) << "): " << flush;
	    charmm_result[0] = MD(i - 1, 1);
	    charmm_result[1] = MD(i, 1);
	    charmm_result[2] = MD(i - 1, 2);
	    charmm_result[3] = MD(i, 2);
	    charmm_result[4] = MD(i - 1, 3);
	    charmm_result[5] = MD(i, 3);
	    cout << std::setw(4) << "dV:" << std::setw(15) << left << std::scientific << std::setprecision(6) << (charmm_result[1] - charmm_result[0]);
	    cout << std::setw(8) << "#Atoms:" << std::setw(10) << left << fixed << std::setprecision(0) << charmm_result[3] << std::setw(8) << "Volume:" << std::setw(12) << left << std::setprecision(2) << charmm_result[5] << endl;
	    i = -1; //break loop
	} else {
	    // in case that only one dV is given, use it for all. If more then break.
	    if (n == 2) {
		cout << std::setw(4) << left << "MD (" << showpoint << STARTD << std::setw(2) << "->" << ENDD << std::setw(3) << "): " << flush;
		charmm_result[0] = MD(0, 1);
		charmm_result[1] = MD(1, 1);
		charmm_result[2] = MD(0, 2);
		charmm_result[3] = MD(1, 2);
		charmm_result[4] = MD(0, 3);
		charmm_result[5] = MD(1, 3);
		cout << std::setw(4) << "dV:" << std::setw(15) << left << std::scientific << std::setprecision(6) << (charmm_result[1] - charmm_result[0]);
		cout << std::setw(8) << "#Atoms:" << std::setw(10) << left << fixed << std::setprecision(0) << charmm_result[3] << std::setw(8) << "Volume:" << std::setw(12) << left << std::setprecision(2) << charmm_result[5] << endl;
		i = -1; //break loop
	    } else {
		dserror("No appropriate MD result found for ENDD");
	    }
	}
    }
    cout.flags(flags);
}


