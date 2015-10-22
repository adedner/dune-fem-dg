#ifndef PARAMETERKEY_HH
#define PARAMETERKEY_HH

namespace Dune
{
  /**
   * \brief helper class to generate parameter keys with an additional prefix.
   *
   * This implemenation adds a prefix to the beginning of an key and adds
   * an separator ".", if the prefix is not empty.
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
