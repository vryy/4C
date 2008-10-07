/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.H

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/rausch
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H"
#include "../drt_adapter/adapter_structure.H"
#include "../drt_structure/strugenalpha.H"
#include <math.h>

class Inv_analysis : public StruGenAlpha {

	private:
                int mp_;      // number of measurment points the
                               // paramters should be fittet to
		vector<double>   disp_;		                // displacment vector in all three direction
                vector<int>      bc_nodes_;
		vector<double>   load_;		                // displacment vector in all three direction
                int problem_type_;

                Epetra_SerialDenseVector final_value_;                                                // calculated displacment
                Epetra_SerialDenseVector final_value_o_;		// calculated displacment of the previous loop
		Epetra_SerialDenseVector measured_value_;	// measured displacment of the experiments (target value)
		Epetra_SerialDenseVector residual_value_;	// final_disp_ - measured_disp
		vector<double> storage_residual_value_;

                double mu_;
		double mu_plus_;
		double mu_minus_;
	        vector<double> storage_mu_;

                double error_;
		double error_o_;
		double tol_;

                int numb_run_;

                Epetra_SerialDenseVector p_;		        // parameter vector
		Epetra_SerialDenseVector p_o_;                       // old parameter vector
		vector<double> p_0_storage_;
		vector<double> p_1_storage_;
		vector<double> p_2_storage_;			// change of paramter to the next step
		bool negative_material_parameters_;

                double aim_function_;
		RefCountPtr<Epetra_CrsMatrix> jacobian_;
		vector<double> parameters_;
		RefCountPtr<Epetra_Vector>  ref_disp_;
		RefCountPtr<Epetra_Vector>  ref_load_;
		const MATERIAL& mat_;
		MATERIAL& nonconstmat_;

                vector<DRT::Condition*> surfneum_;
                vector<DRT::Condition*> surfdir_;


	public:

	  Inv_analysis(ParameterList& params,
	                          DRT::Discretization& dis,
	                          LINALG::Solver& solver,
	                          IO::DiscretizationWriter& output,
	                          bool init = true);
                void reset_parameters();
		void get_surfneum_nodes();
		void get_calculated_curve(const RefCountPtr<Epetra_Vector> value, int numb);
                void calculate_new_parameters();
		void evaluate();
		void Integrate();
                Teuchos::RCP<Epetra_Map> ConditionDofMap(const DRT::Discretization& dis, std::string condname);

		// virtual void Integrate();

};





#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
