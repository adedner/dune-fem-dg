#ifndef DUNE_FEMDG_UNIQUEFUNCTIONNAME_HH
#define DUNE_FEMDG_UNIQUEFUNCTIONNAME_HH

#include <string>
#include <sstream>

namespace Dune
{
namespace Fem
{

  /**
   * \brief This class generates a unique function name which can be used
   * for the
   *
   * Discrete Functions have to have a unique function name. If this is not
   * the case, errors might occur (paraview might crash).
   *
   */
  class FunctionIDGenerator
  {
    public:
      static FunctionIDGenerator& instance ()
      {
        static FunctionIDGenerator generator;
        return generator;
      }

      /**
       *  \brief return the next function name and make next function name
       *  the current one.
       */
      std::string nextId()
      {
        id_++;
        if( id_ == 0 )
          return "";
        std::stringstream s;
        s << "[" << id_ << "]";
        return s.str();
      }

      /**
       * \brief return the current function name.
       */
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
