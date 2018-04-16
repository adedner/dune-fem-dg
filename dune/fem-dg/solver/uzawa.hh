#ifndef DUNE_UZAWASOLVER_HH
#define DUNE_UZAWASOLVER_HH


#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem-dg/misc/integral_constant.hh>

namespace Dune
{
namespace Fem
{
  /**
   * \brief Solves the Stokes equations with an UZAWA iteration
   *
   * \ingroup Solvers
   */
  template < class AssemblerType, class EllipticAlgorithmImp >
  class UzawaSolver : public Dune::Operator< typename AssemblerType::ContainerType::template Item1<1>::DiscreteFunction::DomainFieldType,
                                             typename AssemblerType::ContainerType::template Item1<1>::DiscreteFunction::RangeFieldType,
                                             typename AssemblerType::ContainerType::template Item1<1>::DiscreteFunction,
                                             typename AssemblerType::ContainerType::template Item1<1>::DiscreteFunction >
  {


    typedef typename AssemblerType::ContainerType                            ContainerType;

    typedef typename ContainerType::template Item1<0>::DiscreteFunction      DiscreteFunctionType;
    typedef typename ContainerType::template Item1<1>::DiscreteFunction      PressureDiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType         VeloSpaceType;
    typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureSpaceType;

    // for non istl version
    typedef typename ContainerType::template Item2<0,1>::Matrix              BOPType;
    typedef typename ContainerType::template Item2<1,0>::Matrix              BTOPType;
    typedef typename ContainerType::template Item2<1,1>::Matrix              COPType;

  public:
    /** \todo Please doc me! */
    //!Constructor:
    //!aufSolver is the InverseOperator for Solving the elliptic Problem A^-1
    template< class ContainerImp >
    UzawaSolver( ContainerImp& cont,
                 const EllipticAlgorithmImp& ellAlg,
                 const ParameterReader &parameter = Parameter::container()
               )
      : ellAlg_( ellAlg ),
        bop_( cont(_0,_1)->matrix() ),
        btop_( cont(_1,_0)->matrix() ),
        cop_( cont(_1,_1)->matrix() ),
        rhs1_( cont(_0)->rhs() ),
        rhs2_( cont(_1)->rhs() ),
        spc_( cont(_0)->solution()->space() ),
        pressurespc_( cont(_1)->solution()->space() ),
        velocity_( cont(_0)->solution() ),
        outer_absLimit_( parameter.getValue<double>("istl.absLimit",1.e-10) ),
        maxIter_( 3 * spc_.size() ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        iter_(0),
        linIter_(0)
    {}

    static_assert( (int)DiscreteFunctionType::DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == DiscreteFunctionType::GridType::dimension, "stokes assembler: velocity dimrange does not fit");
    static_assert( (int)PressureDiscreteFunctionType::DiscreteFunctionSpaceType::FunctionSpaceType::dimRange == 1 , "stokes assembler: pressure dimrange does not fit");

    void maxIterations( int maxIter )
    {
      maxIter_ = maxIter;
    }

    /** \todo Please doc me! */
    virtual void operator()(const PressureDiscreteFunctionType& arg,
                            PressureDiscreteFunctionType& pressure ) const
    {
      //get solver
      assert( ellAlg_.solver() );
      const auto& aufSolver = *ellAlg_.solver();

      Dune::Timer timer;
      Dune::Timer timer2;
      timer2.start();


      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType Field;
      Field delta, rho;

      DiscreteFunctionType f("f",spc_);
      // f := rhs1
      f.assign(*rhs1_);
      //DiscreteFunctionType velocity("velocity",spc_);
      velocity_->clear();
      DiscreteFunctionType tmp1("tmp1",spc_);

      tmp1.clear();
      DiscreteFunctionType xi("xi",spc_);
      xi.clear();
      PressureDiscreteFunctionType tmp2("tmp2",pressurespc_);
      tmp2.clear();

      //p<->d
      PressureDiscreteFunctionType d("d",pressurespc_);
      d.clear();
      PressureDiscreteFunctionType h("h",pressurespc_);
      h.clear();
      PressureDiscreteFunctionType g("g",pressurespc_);
      g.clear();
      PressureDiscreteFunctionType residuum("residuum",pressurespc_);


      // residuum = arg
      residuum.assign(arg);
      // B * pressure = tmp1
      bop_->apply(pressure,tmp1);
      // f -= tmp1
      f-=tmp1;

      timer2.stop();
      timer.start();
      // A^-1 * f = velocity
      aufSolver(f,*velocity_);
      timer.stop();
      timer2.start();
      // B^T * velocity = tmp2
      btop_->apply(*velocity_,tmp2);
      //=> tmp2 = B^T * A^-1 * ( F - B * d )

      // residuum -= tmp2
      residuum-=tmp2;
      //=> residuum = arg - B^T * A^-1 * ( F - B * d )
      tmp2.clear();

      // C * pressure = tmp2
      cop_->apply(pressure, tmp2);
      // residuum += tmp2
      residuum += tmp2;

      // d := residuum;
      d.assign(residuum);

      // save iteration number
      linIter_+=aufSolver.iterations();

      // delta = (residuum,residuum)
      delta = residuum.scalarProductDofs( residuum );
      while((delta > outer_absLimit_) && (iter_+=1 < maxIter_))
      {
        tmp1.clear();
        // B * d = tmp1
        bop_->apply(d,tmp1);

        timer2.stop();
        timer.start();
        // A^-1 * tmp1 = xi
        aufSolver(tmp1,xi);
        timer.stop();
        timer2.start();
        // B^T * xi = h
        btop_->apply(xi,h);
        // => h = B^T * A^-1 * B * d
        tmp2.clear();
        cop_->apply( d, tmp2 );
        h += tmp2;

        // rho = delta / d.scalarProductDofs( h );
        rho    = delta / d.scalarProductDofs( h );
        // pressure -= rho * d
        pressure.axpy( -rho, d );
        // velocity += rho * xi
        velocity_->axpy(rho,xi);
        // residuum -= rho * h
        residuum.axpy( -rho,h );

        double oldDelta = delta;

        // delta = (residuum,residuum)
        delta = residuum.scalarProductDofs( residuum );

        // save iteration number
        linIter_+=aufSolver.iterations();
        if( verbose_ > 0)
          std::cout << "SPcg-Iterationen " << iter_ << "   Residuum:"
                    << delta << "   lin. iter:" << aufSolver.iterations() <<std::endl;

        d *= delta / oldDelta;
        d += residuum;
      }
      std::cout << "UZAWA: SPcg-Iterationen " << iter_ << "   Residuum:"
                << delta << "   lin. iter:" << aufSolver.iterations() <<std::endl;
      std::cout << "UZAWA: Solving time (Poisson solves/total time SPcg): " << timer.elapsed() << " / " << timer2.elapsed() << std::endl;
      if( verbose_ > 0)
        std::cout << std::endl;
      //velocity_.assign(velocity);
    }

    int iterations() const {return iter_;}
    double averageLinIter() const {return linIter_/iter_;}

  private:
    // reference to operator which should be inverted

    //the CGSolver for A^-1
    const EllipticAlgorithmImp&                     ellAlg_;

    std::shared_ptr< BOPType >                      bop_;
    std::shared_ptr< BTOPType >                     btop_;
    std::shared_ptr< COPType >                      cop_;

    std::shared_ptr< DiscreteFunctionType >         rhs1_;
    std::shared_ptr< PressureDiscreteFunctionType > rhs2_;

    const VeloSpaceType&                            spc_;
    const PressureSpaceType&                        pressurespc_;
    std::shared_ptr< DiscreteFunctionType >         velocity_;

    // minial error to reach
    typename DiscreteFunctionType::RangeFieldType outer_absLimit_;

    // number of maximal iterations
    int maxIter_;

    // level of output
    int verbose_;

    mutable int iter_;
    mutable int linIter_;
  };

} // end namespace
} // end namespace Dune

#endif

