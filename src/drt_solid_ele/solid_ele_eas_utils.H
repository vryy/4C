/*----------------------------------------------------------------------*/
/*! \file

\brief data structures for EAS related stuff in solid elements

\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_EAS_UTILS_H
#define SOLID_ELE_EAS_UTILS_H

#include "drt_dserror.H"
#include "linalg_serialdensematrix.H"
#include "linalg_fixedsizematrix.H"

namespace STR
{
  namespace ELEMENTS
  {
    enum class EASType  // with meaningfull value for matrix size info
    {
      soh8_easnone,
      eastype_h8_9,
      eastype_h8_21,
      eastype_sh8_7,
      eastype_sh18_9,
      eastype_undefined

    };


    template <int n_cond>
    class CondensatorBlockBase
    {
     public:
      virtual void Condense(LINALG::Matrix<n_cond, n_cond>& kbbInv, LINALG::Matrix<n_cond, 1>& rb,
          double* tar_matrix, double* tar_rhs)
      {
        return;
      }
      virtual double* GetKab()
      {
        dserror("stop");
        return NULL;
      }
      virtual double* GetKba()
      {
        dserror("stop");
        return NULL;
      }
      virtual void Recover(
          LINALG::Matrix<n_cond, n_cond>& kbbInv, double* a_incr, LINALG::Matrix<n_cond, 1>& b_incr)
      {
        dserror("stop");
        return;
      }
      virtual void Clear()
      {
        dserror("stop");
        return;
      }
    };

    template <int n_cond, int n_tar>
    class CondensatorBlock : public CondensatorBlockBase<n_cond>
    {
     public:
      virtual void Condense(LINALG::Matrix<n_cond, n_cond>& kbbInv, LINALG::Matrix<n_cond, 1>& rb,
          double* tar_matrix, double* tar_rhs)
      {
        static LINALG::Matrix<n_tar, n_cond> kabkbbinv(false);
        kabkbbinv.PutScalar(0.);
        kabkbbinv.Multiply(1., kab_, kbbInv, 0.);
        if (tar_matrix)
        {
          LINALG::Matrix<n_tar, n_tar> tarM(tar_matrix, true);
          tarM.Multiply(-1., kabkbbinv, kba_, 1.);
        }
        if (tar_rhs)
        {
          LINALG::Matrix<n_tar, 1> tarR(tar_rhs, true);
          tarR.Multiply(-1., kabkbbinv, rb, 1.);
        }
      }

      virtual void Recover(
          LINALG::Matrix<n_cond, n_cond>& kbbInv, double* a_incr, LINALG::Matrix<n_cond, 1>& b_incr)
      {
        LINALG::Matrix<n_tar, 1> a_incrM(a_incr, true);
        static LINALG::Matrix<n_cond, 1> tmp(false);
        tmp.PutScalar(0.);
        tmp.Multiply(1., kba_, a_incrM, 0.);
        b_incr.Multiply(-1., kbbInv, tmp, 1.);
        return;
      }

      virtual double* GetKab() { return kab_.A(); }
      virtual double* GetKba() { return kba_.A(); }

      virtual void Clear()
      {
        kab_.PutScalar(0.);
        kba_.PutScalar(0.);
      }
      CondensatorBlock() : kba_(true), kab_(true) {}

     protected:
      LINALG::Matrix<n_cond, n_tar> kba_;
      LINALG::Matrix<n_tar, n_cond> kab_;
    };

    class CondensatorBase
    {
     public:
      enum BlockName
      {
        block_eas,
        block_plast,
        block_temp,
        block_disp
      };
      virtual void Condense(double* tar_matrix, double* tar_rhs,
          const STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        dserror("stop");
        return;
      }
      virtual void Recover(std::map<STR::ELEMENTS::CondensatorBase::BlockName, double*> incr,
          const double step_length)
      {
        dserror("stop");
        return;
      }
      virtual double* GetKbbInv()
      {
        dserror("stop");
        return NULL;
      }
      virtual double* GetRb()
      {
        dserror("stop");
        return NULL;
      }
      virtual double* GetB()
      {
        dserror("stop");
        return NULL;
      }
      virtual double* GetBincr()
      {
        dserror("stop");
        return NULL;
      }
      virtual void ScaleStep(const double alpha)
      {
        dserror("stop");
        return;
      }
      virtual double* GetKab(STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        dserror("stop");
        return NULL;
      }
      virtual double* GetKba(STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        dserror("stop");
        return NULL;
      }
      virtual void AddBlock(STR::ELEMENTS::CondensatorBase::BlockName name, int n_tar)
      {
        dserror("stop");
        return;
      }
      virtual void Clear()
      {
        dserror("stop");
        return;
      }
      virtual void Invert()
      {
        dserror("stop");
        return;
      }
      virtual void ReduceStepLength(const double new_step_length)
      {
        dserror("stop");
        return;
      }
    };

    template <int n_cond>
    class Condensator : public CondensatorBase
    {
     public:
      virtual void AddBlock(STR::ELEMENTS::CondensatorBase::BlockName name, int n_tar)
      {
        switch (n_tar)
        {
          case 24:
            blocks_[name] = Teuchos::rcp(new CondensatorBlock<n_cond, 24>);
            break;
          default:
            dserror("block size not implemented");
        }
      }

      virtual void Clear()
      {
        kbbInv_.PutScalar(0.);
        rb_.PutScalar(0.);
        typedef typename std::map<STR::ELEMENTS::CondensatorBase::BlockName,
            Teuchos::RCP<CondensatorBlockBase<n_cond>>>::iterator it;
        for (it p = blocks_.begin(); p != blocks_.end(); ++p) p->second->Clear();
      }

      virtual double* GetB() { return b_.A(); }
      virtual double* GetRb() { return rb_.A(); }
      virtual double* GetBincr() { return b_iter_incr_.A(); }
      virtual double* GetKbbInv() { return kbbInv_.A(); }
      virtual double* GetKab(STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        if (blocks_.find(block) == blocks_.end()) dserror("block not set");
        return blocks_[block]->GetKab();
      }
      virtual double* GetKba(STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        if (blocks_.find(block) == blocks_.end()) dserror("block not set");
        return blocks_[block]->GetKba();
      }

      virtual void Invert()
      {
        const double norm1 = kbbInv_.Norm1();
        if (std::isnan(norm1) || std::isinf(norm1) || norm1 == 0. || norm1 > 1.e100)
        {
          for (std::size_t i = 0; i < kbbInv_.N(); ++i)
            for (std::size_t j = 0; j < kbbInv_.M(); ++j)
              kbbInv_(j, i) = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
          LINALG::FixedSizeSerialDenseSolver<n_cond, n_cond, 1> solve_Kbb;
          solve_Kbb.SetMatrix(kbbInv_);
          int err2 = solve_Kbb.Factor();
          int err = solve_Kbb.Invert();
          if ((err != 0) || (err2 != 0)) dserror("Inversion of EAS block failed");
        }
      }

      virtual void Condense(double* tar_matrix, double* tar_rhs,
          const STR::ELEMENTS::CondensatorBase::BlockName block)
      {
        if (blocks_.find(block) == blocks_.end()) dserror("block not set");
        blocks_[block]->Condense(kbbInv_, rb_, tar_matrix, tar_rhs);
      }

      virtual void Recover(std::map<STR::ELEMENTS::CondensatorBase::BlockName, double*> incr,
          const double step_length)
      {
        b_iter_incr_.PutScalar(0.);
        b_iter_incr_.Multiply(-1., kbbInv_, rb_, 0.);
        for (std::map<STR::ELEMENTS::CondensatorBase::BlockName, double*>::iterator p =
                 incr.begin();
             p != incr.end(); ++p)
        {
          if (blocks_.find(p->first) == blocks_.end()) dserror("block not set");
          blocks_[p->first]->Recover(kbbInv_, p->second, b_iter_incr_);
        }
        b_.Update(step_length, b_iter_incr_, 1.);
        old_step_length_ = step_length;
        return;
      }

      virtual void ReduceStepLength(const double new_step_length)
      {
        b_.Update(-old_step_length_, b_iter_incr_, 1.);
        b_.Update(new_step_length, b_iter_incr_, 1.);
        old_step_length_ = new_step_length;
      }

      //! ctor
      Condensator()
          : kbbInv_(true),
            rb_(true),
            b_(true),
            b_iter_incr_(true),
            b_incr_last_timestep_(true)
            //    ,blocks_(0)
            //
            ,
            old_step_length_(-1.)
      {
      }

     protected:
      LINALG::Matrix<n_cond, n_cond> kbbInv_;
      LINALG::Matrix<n_cond, 1> rb_;
      LINALG::Matrix<n_cond, 1> b_;
      LINALG::Matrix<n_cond, 1> b_iter_incr_;
      LINALG::Matrix<n_cond, 1> b_incr_last_timestep_;
      std::map<STR::ELEMENTS::CondensatorBase::BlockName,
          Teuchos::RCP<CondensatorBlockBase<n_cond>>>
          blocks_;
      double old_step_length_;
    };

  }  // namespace ELEMENTS

  namespace UTILS
  {
    template <::STR::ELEMENTS::EASType eastype>
    struct EasTypeToNumEas
    {
    };
    template <>
    struct EasTypeToNumEas<::STR::ELEMENTS::EASType::eastype_h8_9>
    {
      static const int neas = 9;
    };
    template <>
    struct EasTypeToNumEas<::STR::ELEMENTS::EASType::eastype_h8_21>
    {
      static const int neas = 21;
    };
    template <>
    struct EasTypeToNumEas<::STR::ELEMENTS::EASType::eastype_sh8_7>
    {
      static const int neas = 7;
    };
    template <>
    struct EasTypeToNumEas<::STR::ELEMENTS::EASType::eastype_sh18_9>
    {
      static const int neas = 9;
    };
    template <>
    struct EasTypeToNumEas<::STR::ELEMENTS::EASType::eastype_undefined>
    {
      static const int neas = -1000;
    };
  }  // namespace UTILS
}  // namespace STR

#endif  // SOLID_ELE_EAS_UTILS_H
