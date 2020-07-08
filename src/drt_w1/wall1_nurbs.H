/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D solid-wall elements using NURBS shape functions.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef WALL1_NURBS_H
#define WALL1_NURBS_H

#include "wall1.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Wall1NurbsType : public DRT::ElementType
      {
       public:
        std::string Name() const { return "Wall1NurbsType"; }

        static Wall1NurbsType& Instance();

        virtual DRT::ParObject* Create(const std::vector<char>& data);

        virtual Teuchos::RCP<DRT::Element> Create(
            const std::string eletype, const std::string eledistype, const int id, const int owner);

        virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

        virtual void SetupElementDefinition(
            std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        {
          // do nothing. Definition inserted by normal wall element.
        }

        virtual void NodalBlockInformation(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

        virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
            const double* x0, int numdf, int dimns);

       private:
        static Wall1NurbsType instance_;
      };


      class Wall1Nurbs : public DRT::ELEMENTS::Wall1
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Wall1Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Wall1Nurbs(const Wall1Nurbs& old);

        /*!
        \brief Destructor

        */
        virtual ~Wall1Nurbs();

        /*!
        \brief Deep copy this instance of Wall1 and return pointer to it
        */

        virtual DRT::Element* Clone() const;


        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        virtual int UniqueParObjectId() const
        {
          return Wall1NurbsType::Instance().UniqueParObjectId();
        }


        /// Print this element
        void Print(std::ostream& os) const;

        virtual DRT::ElementType& ElementType() const { return Wall1NurbsType::Instance(); }

        /*!
        \brief Get shape type of element

        \return nurbs4 or nurbs9

        */
        virtual DiscretizationType Shape() const;


        /*!
        \brief Return number of lines of this element.
        */
        virtual int NumLine() const
        {
          if (NumNode() == 9 || NumNode() == 4)
          {
            return 4;
          }
          else
          {
            dserror("Could not determine number of lines");
            return -1;
          }
        }


        /*!
        \brief Get vector of Teuchos::RCPs to the lines of this element
        */
        virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();


        /*!
        \brief Get vector of Teuchos::RCPs to the surfaces of this element
        */
        virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();


       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef WALL1_NURBS_H
