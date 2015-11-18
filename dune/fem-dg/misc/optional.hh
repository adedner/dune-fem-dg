#ifndef DUNE_OPTIONAL_HH
#define DUNE_OPTIONAL_HH

#include <utility>

namespace Dune
{
namespace Fem
{

  template< class Obj >
  class OptionalObject
    : public Obj
  {
  public:
    typedef Obj BaseType;

    template< class... Args >
    OptionalObject( Args&&... args )
      : BaseType( std::forward<Args>( args )... )
    {}

    explicit operator bool() const
    {
      return true;
    }

    const Obj* value() const
    {
      return this;
    }

    Obj* value()
    {
      return this;
    }
  };

  template< class EmptyObj >
  class OptionalNullPtr
    : public EmptyObj
  {
    typedef EmptyObj BaseType;
  public:
    typedef void type;

    template< class... Args >
    OptionalNullPtr( Args&&... args)
      : BaseType( std::forward<Args>( args )... )
    {}

    explicit operator bool() const
    {
      return false;
    }

    const void* value() const
    {
      return nullptr;
    }

    void* value()
    {
      return nullptr;
    }

  };

}
}

#endif
