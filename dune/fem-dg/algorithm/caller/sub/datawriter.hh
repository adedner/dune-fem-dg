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
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem-dg/misc/tupleutility.hh>

namespace Dune
{
namespace Fem
{

  template< class DiscreteFunctionImp >
  class Cons2PrimOutput
  {
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  public:
    typedef std::tuple< DiscreteFunctionType* >      DataType;

    template< class... Args >
    Cons2PrimOutput( Args&&... args )
      : solution_( nullptr ),
        space_( nullptr )
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {
      if( solution_ ) delete solution_;
      if( space_ ) delete space_;
      space_ = new DiscreteFunctionSpaceType( alg->solution().space().gridPart() );
      solution_ = new DiscreteFunctionType( alg->solution().name() + "[cons2prim]", *space_ );
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

      //typedef typename InDiscreteFunctionSpaceType::GridPartType  GridPartType;
      typedef typename InDiscreteFunctionSpaceType::RangeType     InRangeType;
      //typedef typename InDiscreteFunctionSpaceType::DomainType    InDomainType;
      typedef typename OutDiscreteFunctionSpaceType::RangeType    OutRangeType;
      typedef typename OutDiscreteFunctionSpaceType::DomainType   OutDomainType;

      //typedef Dune::Fem::CachingQuadrature< GridPartType, 0 >     QuadratureType;

      const auto& space =  alg->solution().space();
      solution_->clear();

      InRangeType cons(0.0);


      // create local function adapter
      typedef Dune::Fem::LocalAnalyticalFunctionBinder<OutDiscreteFunctionSpaceType> LocalAnalyticalFunctionType;
      LocalAnalyticalFunctionType localAnalyticalFunction(
          [&](const OutDomainType &xgl,double t,const typename OutDiscreteFunctionSpaceType::EntityType& entity)
          {
            OutRangeType prim(0);
            typename InDiscreteFunctionType::LocalFunctionType consLF = alg->solution().localFunction( entity );

            OutDomainType x = entity.geometry().local( xgl );

            consLF.evaluate( x, cons );

            //neglegt zero values
            if( cons[0] > 1e-15 )
              alg->paraview( t, xgl, cons, prim );

            return prim;
          });
      typedef Dune::Fem::LocalFunctionAdapter<LocalAnalyticalFunctionType> LocalAdaptedFunctionType;
      LocalAdaptedFunctionType localAdapted("local adapted function",std::move(localAnalyticalFunction),space_->gridPart(),2*space_->order()+3);

      // interpolate local adpated function over discrete function
      Dune::Fem::interpolate(localAdapted,*solution_);
    }

    template <class ... Args>
    void write(Args&& ... ) const {}

    DataType data() const
    {
      return std::make_tuple( solution_ );
    }

    ~Cons2PrimOutput()
    {
      if( solution_ ) delete solution_;
    }

  private:
    DiscreteFunctionType*       solution_;
    DiscreteFunctionSpaceType*  space_;
  };


  template< class DiscreteFunctionImp >
  class ExactSolutionOutput
  {
    typedef DiscreteFunctionImp                        DiscreteFunctionType;
  public:
    typedef std::tuple< DiscreteFunctionType* >        DataType;

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

    template <class ... Args>
    void write(Args&& ... ) const {}

    DataType data() const
    {
      return std::make_tuple( solution_ );
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
    typedef DiscreteFunctionImp                         DiscreteFunctionType;
  public:
    typedef std::tuple< DiscreteFunctionType* >         DataType;

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

    template <class ... Args>
    void write(Args&& ... ) const {}

    DataType data() const
    {
      return std::make_tuple( solution_ );
    }
  private:
    DiscreteFunctionType* solution_;
  };



  class VoidSolutionOutput
  {
  public:
    typedef std::tuple<>                              DataType;

    template< class... Args >
    VoidSolutionOutput( Args&&... args )
    {}

    template< class SubAlgImp >
    void init( const std::shared_ptr<SubAlgImp>& alg )
    {}

    template< class TimeProviderImp, class SubAlgImp >
    void prepare( TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {}

    template <class ... Args>
    void write(Args&& ... ) const
    {}

    DataType data() const
    {
      return std::make_tuple();
    }
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

    template< int i >
    struct Write
    {
      template< class Tuple, class ... Args >
      static void apply ( Tuple &tuple, Args && ... args )
      {
        std::get<i>( tuple ).write( std::forward<Args>(args)... );
      }
    };

    template< template< int > class Caller >
    using ForLoopType = ForLoop< Caller, 0, numAlgs - 1 >;

  public:

    typedef std::tuple< OutputImp... >                            OutputTupleType;
    typedef tuple_concat_t< typename OutputImp::DataType... >     IOTupleType;


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

    template< class TimeProviderImp, class SubAlgImp >
    void write( TimeProviderImp& tp, const std::shared_ptr<SubAlgImp>& alg )
    {
      ForLoopType< Write >::apply( tuple_, tp, alg );
    }

    IOTupleType dataTuple()
    {
      return dataTuple( tuple_, IndexSequenceType() );
    }

  private:
    template< std::size_t ... i >
    IOTupleType dataTuple ( const OutputTupleType& tuple, std::index_sequence< i ... > )
    {
      return std::tuple_cat( std::tuple_cat( std::get< i >( tuple ).data()...  ) );
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
