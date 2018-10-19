/*!----------------------------------------------------------------------
\file drt_fiber_node.H

\brief   This is basically a (3d-) node with an additional fiber direction.


\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264

\level 2


*----------------------------------------------------------------------*/
#ifndef FIBER_NODE_H
#define FIBER_NODE_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace DRT
{
  namespace FIBER
  {
    class FiberNodeType : public DRT::ParObjectType
    {
     public:
      std::string Name() const { return "FiberNodeType"; }

      static FiberNodeType& Instance() { return instance_; };

      virtual DRT::ParObject* Create(const std::vector<char>& data);

     private:
      static FiberNodeType instance_;
    };

    /*!
    \brief Node with additional fiber information

    */
    class FiberNode : public DRT::Node
    {
     public:
      //! @name Enums and Friends

      /*!
      \brief The Discretization is a friend of the fiber node
      */
      friend class Discretization;

      //@}

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id     (in): A globally unique fiber node
      \param coords (in): vector of nodal coordinates, length 3
      \param fiber  (in): fiber direction, length 3
      \param cosy   (in): circumferential and transverse direction, length 6
      \param angle  (in): helix and transverse angle, length 2
      \param owner  (in): Owner of this node.
      */
      FiberNode(int id, const double* coords, const double* fiber1, const double* fiber2,
          const double* cosy, const double* angle, const int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a fiber node

      \param old (in): The fiber node to copy

      */
      FiberNode(const DRT::FIBER::FiberNode& old);

      /*!
      \brief Deep copy the derived class and return
             pointer to it

      */
      virtual DRT::FIBER::FiberNode* Clone() const;

      /*!
      \brief Destructor

      */
      virtual ~FiberNode();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of drt_parobject.H.

      \return the parobject id
      */
      virtual int UniqueParObjectId() const
      {
        return FiberNodeType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this node

      \param data (in/out): a char vector to pack the data into

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this node

      \param data (in): a char vector to unpack the data from

      */
      virtual void Unpack(const std::vector<char>& data);

      //@}

      //! @name Additional access methods

      /*!
      \brief Return fiber direction
      */
      virtual inline const double* Fiber1() const { return &fiber1_[0]; }

      /*!
      \brief Return fiber direction
      */
      virtual inline const double* Fiber2() const { return &fiber2_[0]; }

      /*!
      \brief Return fiber direction
      */
      virtual inline const double* Cir() const { return &cir_[0]; }

      /*!
      \brief Return fiber direction
      */
      virtual inline const double* Tan() const { return &tan_[0]; }

      /*!
      \brief Return fiber direction
      */
      virtual inline double Helix() const { return helix_; }

      /*!
      \brief Return fiber direction
      */
      virtual inline double Transverse() const { return transverse_; }

      /*!
      \brief Print this node

      \param os ofstrem
      */
      virtual void Print(std::ostream& os) const;

      //@}

     protected:
      //! fiber direction 1
      std::vector<double> fiber1_;

      //! fiber direction 2
      std::vector<double> fiber2_;

      //! circumferential direction
      std::vector<double> cir_;

      //! tangential direction
      std::vector<double> tan_;

      //! helix angle
      double helix_;

      //! transverse angle
      double transverse_;

    };  // class FiberNode

  }  // namespace FIBER

}  // namespace DRT

#endif  // #ifndef FIBER_NODE_H
