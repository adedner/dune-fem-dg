#ifndef DUNE_FEMDG_UNIQUEFUNCTIONNAME_HH
#define DUNE_FEMDG_UNIQUEFUNCTIONNAME_HH

#include <string>
#include <sstream>

namespace Dune
{
namespace Fem
{

  class FunctionIDGenerator
  {
    public:
      static FunctionIDGenerator& instance ()
      {
        static FunctionIDGenerator generator;
        return generator;
      }

      std::string nextId()
      {
        id_++;
        if( id_ == 0 )
          return "";
        std::stringstream s;
        s << "[" << id_ << "]";
        return s.str();
      }
      std::string id()
      {
        if( id_ == 0 )
          return "";
        std::stringstream s;
        s << "[" << id_ << "]";
        return s.str();
      }
    private:
      FunctionIDGenerator () : id_(-1) {}

      int id_;
  };



} // namespace Fem
} // namespace Dune

#endif
