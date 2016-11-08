/*!----------------------------------------------------------------------
\file solver_amgnxn_preconditioner.H

<pre>
\brief Declaration
\level 1
\maintainer Francesc Verdugo
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/
#ifndef AMGNXN_OPERATOR_H
#define AMGNXN_OPERATOR_H

#ifdef HAVE_MueLu

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_UseDefaultTypes.hpp>

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"
#include "solver_amgnxn_smoothers.H"
#include "solver_amgnxn_hierarchies.H"


namespace LINALG
{
  namespace SOLVER
  {

    class AMGnxn_Preconditioner : public PreconditionerType
    {
      public:
        AMGnxn_Preconditioner( FILE * outfile, Teuchos::ParameterList & params);

        virtual void Setup
          (
           bool create,
           Epetra_Operator * matrix,
           Epetra_MultiVector * x,
           Epetra_MultiVector * b
          );

        virtual void Setup(Teuchos::RCP<BlockSparseMatrixBase> A);


        virtual Epetra_Operator * PrecOperator() const ;

        /// linear operator used for preconditioning
        virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const;

        virtual void Print( std::ostream & stream ){stream << "AMGnxn";}
        virtual const std::string getParameterListName() const { return "AMGnxn Parameters"; }

      private:


        // Private variables
        Teuchos::RCP<Epetra_Operator> P_; // The underlying preconditioner object
        Teuchos::RCP<BlockSparseMatrixBase> A_; // A own copy of the system matrix
        Teuchos::ParameterList & params_;

    };//AMGnxn_Preconditioner

    class AMGnxn_Interface
    {
      public:
        AMGnxn_Interface(Teuchos::ParameterList& params, int NumBlocks);

        std::vector<std::string> GetMueLuXmlFiles(){ return xml_files_;}
        std::vector<int> GetNumPdes(){return num_pdes_;}
        std::vector<int> GetNullSpacesDim(){return null_spaces_dim_;}
        std::vector<Teuchos::RCP<std::vector<double> > > GetNullSpacesData(){return null_spaces_data_;}
        //int GetNumLevelAMG(){return NumLevelAMG_;}
        Teuchos::ParameterList GetPreconditionerParams(){return prec_params_;}
        Teuchos::ParameterList GetSmoothersParams(){return smoo_params_;}
        std::string GetPreconditionerType(){return prec_type_;}

      private:
        std::vector<std::string> xml_files_;
        std::vector<int> num_pdes_;
        std::vector<int> null_spaces_dim_;
        std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data_;
        //int NumLevelAMG_;
        Teuchos::ParameterList prec_params_;
        Teuchos::ParameterList smoo_params_;
        std::string prec_type_;

        void Params_TSI_AMG_BGS(Teuchos::ParameterList& params);

        // Helper function to convert int to std::string
        std::string ConvertInt(int number)
        {
          std::stringstream ss;
          ss << number;
          return ss.str();
        }
    };


    class AMGnxn_Operator : virtual public Epetra_Operator
    {
      public:
        AMGnxn_Operator(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            std::vector<int> num_pdes,
            std::vector<int> null_spaces_dim,
            std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
            const Teuchos::ParameterList& amgnxn_params,
            const Teuchos::ParameterList& smoothers_params,
            const Teuchos::ParameterList& muelu_params);


        // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
        virtual int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        virtual int   SetUseTranspose (bool UseTranspose)
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual double   NormInf () const
        {
          dserror("Function not implemented");
          return -1.0;
        }

        virtual const char *   Label () const
        {
          return "AMG(BlockSmoother)";
        }

        virtual bool   UseTranspose () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual bool   HasNormInf () const
        {
          dserror("Function not implemented");
          return false;
        }
        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Comm &   Comm () const {return A_->Comm();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorDomainMap () const {return A_->OperatorDomainMap();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorRangeMap () const {return A_->OperatorRangeMap();}

        void Setup();


      private:
        Teuchos::RCP<BlockSparseMatrixBase> A_;
        std::vector<Teuchos::ParameterList> muelu_lists_;
        std::vector<int> num_pdes_;
        std::vector<int> null_spaces_dim_;
        std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data_;
        Teuchos::ParameterList amgnxn_params_;
        Teuchos::ParameterList smoothers_params_;
        Teuchos::ParameterList muelu_params_;

        bool is_setup_flag_;

        Teuchos::RCP<AMGNXN::CoupledAmg> V_;


    }; // class AMGnxn_Operator

    class BlockSmoother_Operator : virtual public Epetra_Operator
    {
      public:
        BlockSmoother_Operator(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            std::vector<int> num_pdes,
            std::vector<int> null_spaces_dim,
            std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
            const Teuchos::ParameterList& amgnxn_params,
            const Teuchos::ParameterList& smoothers_params);


        // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
        virtual int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        virtual int   SetUseTranspose (bool UseTranspose)
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual double   NormInf () const
        {
          dserror("Function not implemented");
          return -1.0;
        }

        virtual const char *   Label () const
        {
          return "BlockSmoother(X)";
        }

        virtual bool   UseTranspose () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual bool   HasNormInf () const
        {
          dserror("Function not implemented");
          return false;
        }
        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Comm &   Comm () const {return A_->Comm();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorDomainMap () const {return A_->OperatorDomainMap();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorRangeMap () const {return A_->OperatorRangeMap();}

        void Setup(); //TODO

      private:
        Teuchos::RCP<BlockSparseMatrixBase> A_;
        std::vector<int> num_pdes_;
        std::vector<int> null_spaces_dim_;
        std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data_;
        Teuchos::ParameterList amgnxn_params_;
        Teuchos::ParameterList smoothers_params_;

        bool is_setup_flag_;
        Teuchos::RCP<AMGNXN::BlockedSmoother> S_;
        Teuchos::RCP<AMGNXN::GenericSmoother> Sbase_;

    }; // class BlockSmoother_Operator



    class Merged_Operator : virtual public Epetra_Operator
    {
      public:
        Merged_Operator(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            const Teuchos::ParameterList& amgnxn_params,
            const Teuchos::ParameterList& smoothers_params);


        // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
        virtual int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        virtual int   SetUseTranspose (bool UseTranspose)
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual double   NormInf () const
        {
          dserror("Function not implemented");
          return -1.0;
        }

        virtual const char *   Label () const
        {
          return "Merged matrix plus smoother";
        }

        virtual bool   UseTranspose () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual bool   HasNormInf () const
        {
          dserror("Function not implemented");
          return false;
        }
        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Comm &   Comm () const {return A_->Comm();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorDomainMap () const {return A_->OperatorDomainMap();}

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorRangeMap () const {return A_->OperatorRangeMap();}

        void Setup(); //TODO

      private:
        Teuchos::RCP<BlockSparseMatrixBase> A_;
        Teuchos::RCP<LINALG::SparseMatrix>  Asp_;
        Teuchos::ParameterList amgnxn_params_;
        Teuchos::ParameterList smoothers_params_;

        bool is_setup_flag_;
        Teuchos::RCP<AMGNXN::IfpackWrapper> S_;

    }; // class Merged_Operator



    void PrintMap(const Epetra_Map& Map,std::string prefix);


  } //namespace LINALG
} // namespace SOLVER

#endif // HAVE_MueLu
#endif // AMGNXN_OPERATOR_H
