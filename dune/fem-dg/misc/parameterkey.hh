#ifndef PARAMETERKEY_HH
#define PARAMETERKEY_HH

#include <string>

namespace Dune
{
  /**
   * \brief helper class to generate parameter keys with an additional prefix.
   */
  struct ParameterKey
  {
    static std::string generate( const std::string& prefix, const std::string& key )
    {
      return prefix + (prefix!=""? "." : "" ) + key;
    }
  };
}

#endif
