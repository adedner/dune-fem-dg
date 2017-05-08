#ifndef FEMDG_EOCERRORLIST_HH
#define FEMDG_EOCERRORLIST_HH

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem-dg/misc/dgnorm.hh>
#include <dune/fem/misc/femeoc.hh>
namespace Dune
{
namespace Fem
{

  template< int i = 0 >
  class SubEOCErrorList
  {
  public:
    SubEOCErrorList()
    : names_( 0 ),
      ids_( 0 )
    {}


    template< class Measure, class Solution, class Model >
    static void setErrors ( Model& model, Solution &u )
    {
      const Measure measure( u );
      instance().addError( measure );
      int id = instance().id( measure );
      //only insert, if added in list
      if( id >= 0 )
        measure.add( model, u, id );
    }

  private:
    static SubEOCErrorList<i>& instance()
    {
      static SubEOCErrorList<i> list;
      return list;
    }

    template< class Measure >
    void addError( const Measure& measure )
    {
      auto it = std::find( std::begin(names_), std::end(names_), measure.name() );
      if( it == std::end(names_) )
      {
        names_.push_back( measure.name() );
        ids_.push_back( Dune::Fem::FemEoc::addEntry( names_.back() ) );
      }
    }

    template< class Measure >
    int id( const Measure& measure ) const
    {
      auto it = std::find( std::begin(names_), std::end(names_), measure.name() );
      return (it == std::end(names_)) ? -1 : std::distance( std::begin(names_), it);
    }

    std::vector< std::string > names_;
    std::vector< int > ids_;
  };

  using EOCErrorList = SubEOCErrorList<0>;

}
}
#endif
