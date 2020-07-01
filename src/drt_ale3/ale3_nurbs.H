/*----------------------------------------------------------------------------*/
/*! \file

\brief A nurbs implementation of the ale3 element

\level 2

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE3_NURBS_H
#define ALE3_NURBS_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ale3.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Ale3_NurbsType : public DRT::ElementType
      {
       public:
        std::string Name() const { return "Ale3_NurbsType"; }

        static Ale3_NurbsType& Instance();

        virtual DRT::ParObject* Create(const std::vector<char>& data);

        virtual Teuchos::RCP<DRT::Element> Create(
            const std::string eletype, const std::string eledistype, const int id, const int owner);

        virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

        virtual void NodalBlockInformation(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

        virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
            const double* x0, int numdf, int dimns);

        virtual void SetupElementDefinition(
            std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale3_NurbsType instance_;
      };

      /*----------------------------------------------------------------------------*/
      /*----------------------------------------------------------------------------*/
      class Ale3Nurbs : public DRT::ELEMENTS::Ale3
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Ale3Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Ale3Nurbs(const Ale3Nurbs& old);

        /*!
        \brief Destructor

        */
        virtual ~Ale3Nurbs();


        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        virtual int UniqueParObjectId() const
        {
          return Ale3_NurbsType::Instance().UniqueParObjectId();
        }

        virtual DRT::ElementType& ElementType() { return Ale3_NurbsType::Instance(); }

        /// Print this element
        void Print(std::ostream& os) const;


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
          return (0);
          /*
          if (NumNode()==27 || NumNode()==8)
          {
            return 12;
          }
          else {
            dserror("Could not determine number of lines");
            return -1;
          }
          */
        }


        /*!
        \brief Return number of surfaces of this element
        */
        virtual int NumSurface() const
        {
          if (NumNode() == 27 || NumNode() == 8)
          {
            return 6;
          }
          else
          {
            dserror("Could not determine number of surfaces");
            return -1;
          }
        }

       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef ALE3_NURBS_H
