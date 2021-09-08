#ifndef DUNE_FEM_DG_SIMULATOR_HH
#define DUNE_FEM_DG_SIMULATOR_HH

#if defined GRIDDIM
#ifndef CODEDIM
#define CODEDIM GRIDDIM
#endif
#endif

// to make headercheck work
#ifndef POLORDER
#warning "Using default POLORDER = 1"
#define POLORDER 1
#endif

// in dbug mode also enable FieldVector checking and dune devel mode
#ifndef NDEBUG
#define DUNE_ISTL_WITH_CHECKING
#define DUNE_DEVEL_MODE
#endif

// -1 means higher order FV
#if POLORDER == -1
#define HIGHER_ORDER_FV
#undef POLORDER
#define POLORDER 0
#endif

// when POLORDER is defined this polynomial order is chosen to compile the program
// otherwise optional values from MIN_POLORD to MAX_POLORD are provided.
// This will increase compile time and program size and should only be chosen
// for special purpose
#if POLORDER
#define MIN_POLORD POLORDER
#define MAX_POLORD POLORDER
#else
#ifndef MIN_POLORD
#define MIN_POLORD 0
#endif
#ifndef MAX_POLORD
#define MAX_POLORD 4
#endif
#endif

#include <memory>
#include <dune/common/version.hh>

// streams for backup
#include <dune/fem-dg/misc/streams.hh>

#if defined USE_BASEFUNCTIONSET_CODEGEN || defined BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basisfunctionset/codegen.hh>
#endif

#ifdef USE_BASEFUNCTIONSET_CODEGEN
#warning "Using autogenerated code"
#include <autogeneratedcode.hh>
#endif

#include <dune/fem/misc/threads/threadpool.hh>
#if defined USE_SMP_PARALLEL
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem-dg/pass/threadpass.hh>

#if HAVE_DUNE_FEM_DG
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

#include <dune/fem-dg/algorithm/base.hh>

#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/flops.hh>
#include <dune/fem/common/forloop.hh>

// Dune includes
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basisfunctionset/codegen.hh>

std::string autoFilename(const int dim, const int polynomialOrder )
{
  std::stringstream name;
  name << "autogeneratedcode_" << dim << "d_k" << polynomialOrder << ".hh";
  return name.str();
}

void finalizeCodegen()
{
  std::string path = Dune::Fem::Parameter::commonInputPath() + "/autogeneratedcode";
  //////////////////////////////////////////////////
  //  write include header
  //////////////////////////////////////////////////
  std::ofstream file( path + ".hh" );

  if( file )
  {
    // we only cover dimension for 1 to 3
    for( int dim=1; dim<4; ++dim )
    {
      std::stringstream dimfilename;
      dimfilename      << path << "/autogeneratedcode_" << dim << "d.hh";
      std::stringstream dimfilenameflops;
      dimfilenameflops << path << "flops/autogeneratedcode_" << dim << "d.hh";
      file << "#if CODEDIM == " << dim << std::endl;
      file << "#ifdef COUNT_FLOPS" << std::endl;
      file << "#define DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC \"" << dimfilenameflops.str() << "\"" << std::endl;
      file << "#else" << std::endl;
      file << "#define DUNE_FEM_INCLUDE_AUTOGENERATEDCODE_FILENAME_SPEC \"" << dimfilename.str() << "\"" << std::endl;
      file << "#endif" << std::endl;
      file << "#endif" << std::endl;

#ifdef COUNT_FLOPS
      std::ofstream dimfile( dimfilenameflops.str().c_str() );
#else
      std::ofstream dimfile( dimfilename.str().c_str() );
#endif
      if( dimfile )
      {
        const int maxPolOrder = ( dim == 3 ) ? 4 : 8 ;
        // max polorder is 4 in 3d and 8 in 2d at the moment
        for( int i=0; i <= maxPolOrder; ++ i )
        {
          dimfile << "#if POLORDER == " << i << std::endl;
          dimfile << "#include \"" << autoFilename( dim, i ) << "\"" << std::endl;
          dimfile << "#endif" << std::endl;
        }
      }
    }
  }
  // dump all information
  std::cerr << "All automated code generated, bye, bye !! " << std::endl;
}
#endif

bool readParameters( int argc, char ** argv )
{
  // append all parameters from command line
  Dune::Fem::Parameter::append(argc,argv);
  // no parameters defined: search for parameter file
  if( argc <= 1 )
    Dune::Fem::Parameter::append("parameter");


  //get full test case path
  const std::string testCasePath = Dune::Fem::Parameter::getValue< std::string >( "testcase.path", "" );
  //get name for a test case
  const std::string testCaseName = Dune::Fem::Parameter::getValue< std::string >( "testcase", "" );

  std::string path( argv[0] );
  std::string targetName = path.substr(path.find_last_of("/") + 1);

  //append test case name, if given
  if( testCaseName != "" )
    targetName = targetName + "_" + testCaseName;
  //set important paths which are related to the test case path
  if( testCasePath != "" )
  {
    Dune::Fem::Parameter::append("fem.prefix.input",testCasePath);
    Dune::Fem::Parameter::append("fem.prefix",testCasePath + "data/" + targetName );
    Dune::Fem::Parameter::append("fem.eoc.outputpath",testCasePath + "data/" + targetName );
  }
  else
  {
    //fail: We need data to provide further information via parameter file
    if( !Dune::Fem::Parameter::exists( "fem.prefix" ) )
    {
      std::cout << "Please specify parameter 'fem.prefix'" << std::endl;
      return false;
    }
    if( !Dune::Fem::Parameter::exists( "fem.prefix.input" ) )
    {
      std::cout << "Please specify parameter 'fem.prefix.input'" << std::endl;
      return false;
    }
    if( !Dune::Fem::Parameter::exists( "fem.eoc.outputpath" ) )
    {
      std::cout << "Please specify parameter 'fem.eoc.outputpath'" << std::endl;
      return false;
    }
  }

  //read problem dependend file, if existent
  std::string fullpath = Dune::Fem::Parameter::commonInputPath() + "/parameters/" + targetName;
  if( std::ifstream( fullpath ) )
    Dune::Fem::Parameter::append( fullpath );
  else
  {
    std::cout << "Parameter file not found in " << fullpath << ". Aborting." << std::endl;
    return false;
  }

  // write parameters used (before simulation starts)
  Dune::Fem::Parameter::write("parameter.log");

  return true;
}



namespace Dune
{
namespace Fem
{

  struct FlopStartObject
  {
    FlopStartObject()
    {
      // initialize counters for master thread before all others
      (*this)();
    }

    void operator() () const
    {
      Dune::Fem::FlopCounter::start();
    }
  };

  struct FlopStopObject
  {
    void operator() () const
    {
      Dune::Fem::FlopCounter::stop();
    }
  };

  template <int polynomialOrder, class AlgorithmCreator >
  inline void simulate(const AlgorithmCreator& algCreator)
  {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
    //only one thread for codegen
    const int numThreads = 1;
#else
    // get number of desired threads (default is 1)
    const int numThreads = Dune::Fem::Parameter::getValue< int >("fem.parallel.numberofthreads", 1);
#endif
    // set maximal number of threads
    Dune :: Fem :: ThreadManager :: setMaxNumberThreads( numThreads );
    Dune :: Fem :: ThreadManager :: setNumThreads( numThreads );

    const bool countFlops = Dune::Fem::Parameter::getValue< bool >("femdg.flopcounter", false );

    // if flop count is enabled count floating point operations (PAPI needed)
    // start flop counters for all threads
    if( countFlops )
    {
      FlopStartObject startObj ;
      Dune::Fem::ThreadPool::run( startObj );
    }

    const auto& globalContainer = algCreator.template initContainer<polynomialOrder>();
    typedef typename AlgorithmCreator::template Algorithm< polynomialOrder > AlgorithmType;
    auto algorithm = std::make_unique< AlgorithmType >( algCreator.moduleName(), globalContainer );

    // run the algorithm
    compute( *algorithm );

    // stop flop counters for all threads
    if( countFlops )
    {
      FlopStopObject stopObj ;
      Dune::Fem::ThreadPool::run( stopObj );
      // print results
      Dune::Fem::FlopCounter::print( std::cout );
    }
  } // end simulate

  template <int polOrd>
  struct SimulatePolOrd
  {
    template< class AlgorithmCreator >
    static void apply( const AlgorithmCreator& algCreator, const int polynomialOrder, const bool computeAnyway )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      Dune::Fem::Codegen::CodegenInfo :: instance().setFileName( autoFilename( CODEDIM, polOrd ) );
#endif
      if( computeAnyway || polOrd == polynomialOrder )
      {
        if( Dune::Fem::Parameter::verbose() )
        {
          int pOrd = computeAnyway ? POLORDER : polOrd ;
          std::cout << "Simulator: run for polynomialOrder = " << pOrd << std::endl;
        }

        simulate< polOrd > ( algCreator );
      }
    }
  };

  struct Simulator
  {
    template <class AlgorithmCreator>
    static void run( const AlgorithmCreator& algCreator )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      using namespace Dune::Fem::Codegen;
#if COUNT_FLOPS
      CodegenInfo :: instance().setPath(
          Dune::Fem::Parameter::commonInputPath() + "/autogeneratedcodeflops");
#endif
      std::cout << "Generating Code \n";
      try
      {
#endif
        int polOrder = 1;
        polOrder = Dune::Fem::Parameter :: getValue("femdg.polynomialOrder", polOrder );

        // run through all available polynomial order and check with dynamic polOrder
        // when -DONLY_ONE_P was passed only POLORDER is used
        Dune::Fem::ForLoop< SimulatePolOrd, MIN_POLORD, MAX_POLORD > :: apply( algCreator, polOrder, bool(MIN_POLORD == MAX_POLORD) );
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      }
      catch (const CodegenInfoFinished& ) {}

      std::cerr << "Code for k="<< MAX_POLORD << " generated!! " << std::endl;
      CodegenInfo :: instance().dumpInfo();
      CodegenInfo :: instance().clear();
      std::remove( autoFilename( CODEDIM, MAX_POLORD ).c_str() );
      finalizeCodegen();
#endif
    }
  };

} // namespace Fem
} // namespace Dune
#endif // #ifndef DUNE_FEM_DG_SIMULATOR_HH
