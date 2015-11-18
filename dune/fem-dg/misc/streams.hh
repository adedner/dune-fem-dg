#ifndef DUNE_FEM_DG_STREAMS_HH
#define DUNE_FEM_DG_STREAMS_HH

#if HAVE_SIONLIB
#if USE_SIONLIB
#warning "using SIONlib streams for I/O"

#include <dune/fem/io/streams/sionlibstreams.hh>

namespace Dune
{
namespace Fem
{

  struct PersistenceManagerTraits
  {
    typedef Fem :: SIONlibOutStream  BackupStreamType ;
    typedef Fem :: SIONlibInStream   RestoreStreamType ;
    static const bool singleBackupRestoreFile = true ;
  };

#define FEM_PERSISTENCEMANAGERSTREAMTRAITS  PersistenceManagerTraits
}
} // namespace Dune

#endif // #if USE_SIONLIB
#endif // #if HAVE_SIONLIB

#endif // #ifndef DUNE_FEM_DG_STREAMS_HH
