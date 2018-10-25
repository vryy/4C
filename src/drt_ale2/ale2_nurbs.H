/*----------------------------------------------------------------------------*/
/*!
\file ale2_nurbs.H

\brief Nurbs verison of 2D ALE element

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE2_NURBS_H
#define ALE2_NURBS_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ale2.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Ale2_NurbsType : public Ale2Type
      {
       public:
        std::string Name() const { return "Ale2_NurbsType"; }

        static Ale2_NurbsType& Instance();

        virtual DRT::ParObject* Create(const std::vector<char>& data);

        virtual Teuchos::RCP<DRT::Element> Create(
            const std::string eletype, const std::string eledistype, const int id, const int owner);

        virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

        virtual void SetupElementDefinition(
            std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale2_NurbsType instance_;
      };


      class Ale2Nurbs : public DRT::ELEMENTS::Ale2
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Ale2Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Ale2Nurbs(const Ale2Nurbs& old);

        /*!
        \brief Destructor

        */
        virtual ~Ale2Nurbs();


        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        virtual int UniqueParObjectId() const
        {
          return Ale2_NurbsType::Instance().UniqueParObjectId();
        }


        /// Print this element
        void Print(std::ostream& os) const;

        virtual DRT::ElementType& ElementType() const { return Ale2_NurbsType::Instance(); }

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


       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef ALE2_NURBS_H
