#ifndef DUNE_FEM_DG_STREAMS_HH
#define DUNE_FEM_DG_STREAMS_HH

#if HAVE_SIONLIB && USE_SIONLIB
#warning "using SIONlib streams for output"

#include <dune/fem/io/streams/sionlibstreams.hh>
namespace Dune {

  struct PersistenceManagerTraits
  {
    typedef Fem :: SIONlibOutStream  BackupStreamType ;
    typedef Fem :: SIONlibInStream   RestoreStreamType ;
    static const bool singleBackupRestoreFile = true ;
  };

#define FEM_PERSISTENCEMANAGERSTREAMTRAITS  PersistenceManagerTraits
}
#endif // #if HAVE_SIONLIB && USE_SIONLIB

#endif // #ifndef DUNE_FEM_DG_STREAMS_HH
