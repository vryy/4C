/*!----------------------------------------------------------------------
\file proein.cpp
\brief CHARMm Interface to compute the mechanical properties of proteins

<pre>
Maintainer: Robert Metzke
	    metzke@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de
	    089 - 289-15244
</pre>
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "protein.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PROTEIN::PROTEIN(Teuchos::RCP<MAT::PAR::Material> matdata)
: Parameter(matdata),
density_(matdata->GetDouble("DENS")) {
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PROTEIN::PROTEIN()
: params_(NULL) {
}

MAT::PROTEIN::PROTEIN(MAT::PAR::PROTEIN* params)
: params_(params) {
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PROTEIN::Pack(vector<char>& data) const {
    data.resize(0);

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
void MAT::PROTEIN::Unpack(const vector<char>& data) {
    int position = 0;
    // extract type
    int type = 0;
    ExtractfromPack(position, data, type);
    if (type != UniqueParObjectId()) dserror("wrong instance type data");

    // matid
    int matid;
    ExtractfromPack(position, data, matid);
    // in post-process mode we do not have any instance of DRT::Problem
    if (DRT::Problem::NumInstances() > 0) {
	const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
	MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
	if (mat->Type() == MaterialType())
	    params_ = static_cast<MAT::PAR::PROTEIN*> (mat);
	else
	    dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    } else {
	params_ = NULL;
    }

    if (position != (int) data.size())
	dserror("Mismatch in size of data %d <-> %d", (int) data.size(), position);
}


/*----------------------------------------------------------------------*/
//! Compute second PK and constitutive tensor
/*----------------------------------------------------------------------*/
void MAT::PROTEIN::Evaluate(const LINALG::Matrix<NUM_STRESS_3D, 1 > * glstrain,
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

    dserror("Protein CHARMm API: not implemeted yet");
    return;
}


#endif

