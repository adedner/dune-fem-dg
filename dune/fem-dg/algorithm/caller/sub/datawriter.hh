#ifndef FEMDG_SUBDATAWRITER_HH
#define FEMDG_SUBDATAWRITER_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem-dg/misc/optional.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>


namespace Dune
{
namespace Fem
{

  template< class DiscreteFunctionImp >
  class Cons2PrimOutput
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;

    template< class... Args >
    Cons2PrimOutput( Args&&... args )
      : solution_( nullptr )
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ ) delete solution_;
      solution_ = new DiscreteFunctionType( alg->solution().name() + "[cons2prim]", alg->solution().space() );
    }

    /** \brief converts a discrete function of conservative variables to
     *    a discrete function of primitive variables for a visualization purpose only
     */
    template< class TimeProviderImp, class SubAlgImp >
    void prepare( const TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ == nullptr )
        init( alg );

      typedef typename SubAlgImp::DiscreteFunctionType                      InDiscreteFunctionType;
      typedef typename InDiscreteFunctionType::DiscreteFunctionSpaceType    InDiscreteFunctionSpaceType;
      typedef DiscreteFunctionImp                                           OutDiscreteFunctionType;
      typedef typename OutDiscreteFunctionType::DiscreteFunctionSpaceType   OutDiscreteFunctionSpaceType;

      typedef typename InDiscreteFunctionSpaceType::GridPartType  GridPartType;
      typedef typename InDiscreteFunctionSpaceType::RangeType     InRangeType;
      typedef typename OutDiscreteFunctionSpaceType::RangeType    OutRangeType;

      const auto& space =  alg->solution().space();
      solution_->clear();

      InRangeType cons(0.0);
      OutRangeType prim(0.0);

      for( const auto& entity : elements( space.gridPart() ) )
      {
        const auto& geo = entity.geometry();

        // get quadrature rule for L2 projection
        Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );

        typename InDiscreteFunctionType::LocalFunctionType consLF = alg->solution().localFunction( entity );
        typename OutDiscreteFunctionType::LocalFunctionType primLF = solution_->localFunction( entity );


        for( const auto qp : quad )
        {
          const auto& xgl = geo.global( qp.position() );
          consLF.evaluate( qp, cons );

          // it is useful to visualize better suited quantities
          alg->paraview( tp.time(), xgl, cons, prim );

          prim *=  qp.weight();
          primLF.axpy( qp, prim );
        }
      }
    }

    DiscreteFunctionType* data() const
    {
      return solution_;
    }

    ~Cons2PrimOutput()
    {
      if( solution_ ) delete solution_;
    }

  private:
    DiscreteFunctionType*   solution_;
  };


  template< class DiscreteFunctionImp >
  class ExactSolutionOutput
  {
  public:
    typedef DiscreteFunctionImp                        DiscreteFunctionType;

    template< class... Args >
    ExactSolutionOutput( Args&&... args )
      : solution_( nullptr )
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ ) delete solution_;
      solution_ = new DiscreteFunctionType( alg->solution().name() + "[exact]", alg->solution().space() );
    }


    template< class TimeProviderImp, class SubAlgImp >
    void prepare( TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ == nullptr )
        init( alg );

      auto ftf = alg->exactSolution( tp.time() );
      interpolate( gridFunctionAdapter( ftf, solution_->space().gridPart(), solution_->space().order()+2 ), *solution_ );
    }

    DiscreteFunctionType* data() const
    {
      return solution_;
    }

    ~ExactSolutionOutput()
    {
      if( solution_ ) delete solution_;
    }

  private:
    DiscreteFunctionType* solution_;
  };


  template< class DiscreteFunctionImp >
  class SolutionOutput
  {
  public:
    typedef DiscreteFunctionImp                         DiscreteFunctionType;

    template< class... Args >
    SolutionOutput( Args&&... args )
      : solution_( nullptr )
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {
      solution_ = &(alg->solution());
    }

    template< class TimeProviderImp, class SubAlgImp >
    void prepare( TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ == nullptr )
        init( alg );
    }

    DiscreteFunctionType* data() const
    {
      return solution_;
    }
  private:
    DiscreteFunctionType* solution_;
  };





  template< class... OutputImp >
  class SubDataWriter
  {
    static const int numAlgs = sizeof...( OutputImp );

    typedef std::make_index_sequence< numAlgs >                   IndexSequenceType;

    template< int i >
    struct Init
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get<i>( tuple ).init( std::forward<Args>(args)... );
      }
    };

    template< int i >
    struct Prepare
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get<i>( tuple ).prepare( std::forward<Args>(args)... );
      }
    };

    template< template< int > class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    typedef std::tuple< OutputImp... >                                 OutputTupleType;
    typedef std::tuple< typename OutputImp::DiscreteFunctionType*... > IOTupleType;

    SubDataWriter( const std::string keyPrefix = "" )
      : tuple_( outputTuple( IndexSequenceType()) ) //initialize with nullptr
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {
      ForLoopType< Init >::apply( tuple_, alg );
    }

    template< class TimeProviderImp, class SubAlgImp >
    void prepare( TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {
      ForLoopType< Prepare >::apply( tuple_, tp, alg );
    }

    IOTupleType dataTuple()
    {
      return dataTuple( tuple_, IndexSequenceType() );
    }

  private:
    template< std::size_t ... i >
    IOTupleType dataTuple ( const OutputTupleType &tuple, std::index_sequence< i ... > )
    {
      return std::make_tuple( (std::get< i >( tuple ).data() )... );
    }

    template< std::size_t ... i >
    OutputTupleType outputTuple ( std::index_sequence< i ... > )
    {
      return std::make_tuple( std::tuple_element_t<i,OutputTupleType>()... );
    }

    OutputTupleType tuple_;
  };


  template<>
  class SubDataWriter<>
  {
  public:
    typedef std::tuple<>                       IOTupleType;

    template< class... Args >
    SubDataWriter( Args&& ... )
    {}

    template <class ... Args>
    void init(Args&& ... ) const {}

    template <class ... Args>
    void prepare(Args&& ... ) const {}

    IOTupleType dataTuple()
    {
      return {};
    }
  private:
  };


  template< class Obj >
  class DataWriterOptional
    : public OptionalObject< Obj >
  {
    typedef OptionalObject< Obj >    BaseType;
  public:
    template< class... Args >
    DataWriterOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };

  template<>
  class DataWriterOptional< void >
    : public OptionalNullPtr< SubDataWriter<> >
  {
    typedef OptionalNullPtr< SubDataWriter<> >    BaseType;
  public:
    template< class... Args >
    DataWriterOptional( Args&&... args )
      : BaseType( std::forward<Args>(args)... )
    {}
  };
}
}

#endif
