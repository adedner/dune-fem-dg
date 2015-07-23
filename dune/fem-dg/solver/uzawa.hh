#ifndef DUNE_UZAWASOLVER_HH
#define DUNE_UZAWASOLVER_HH


#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

//#include "../common/mymatrix.hh"
//#include "../common/spmatextra.hh"

//#include "matrixoperator.hh"

namespace Dune {
  //!CG Verfahren fuer Sattelpunkt Problem
  //!siehe NavSt Skript Siebert,Seite 36, Alg 3.34
  /** \brief Inversion operator using CG algorithm
   */


  template <class DiscreteFunctionType,class PressureDiscreteFunctionType,class AssemblerType,class InverseOperatorType>
  class UzawaSolver : public Operator<
    typename PressureDiscreteFunctionType::DomainFieldType,
    typename PressureDiscreteFunctionType::RangeFieldType,
    PressureDiscreteFunctionType,PressureDiscreteFunctionType>
  {

    // typedef Dune::SparseRowMatrixExtra< double > MatrixType ;

    typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureSpaceType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType VeloSpaceType;

    typedef AssemblerType MappingType;

    /**************** for non istl version*/
    typedef typename MappingType::PressureGradMatType BOPType;
    typedef typename MappingType::PressureDivMatType BTOPType;
    typedef typename MappingType::PressureStabMatType COPType;
    /*****************************************/

    // typedef typename OperatorType::BOPType BOPType;
    //     typedef typename OperatorType::BTOPType BTOPType;
    //     typedef typename OperatorType::COPType COPType;
  public:
    /** \todo Please doc me! */
    //!Constructor:
    //!aufSolver is the InverseOperator for Solving the elliptic Problem A^-1
    //!rhs1 is  stored as member,no better idea
    UzawaSolver(const MappingType& op,
                const InverseOperatorType& aufSolver,
                double redEps,
                double absLimit,
                int maxIter,
                int verbose=1
               )
      : op_(op), redEps_( redEps ), outer_absLimit_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_( verbose ),aufSolver_(aufSolver), bop_(op_.getBOP()),
        btop_(op_.getBTOP()),
        cop_(op_.getCOP()),
        rhs1_(aufSolver_.inv().affineShift()),
        rhs2_(op_.pressureRhs()),
        pressurespc_(op_.pressurespc()),
        spc_(op.spc()),
        velocity_("VELO",spc_),
        iter_(0),
        linIter_(0),
        tau_(0.2),
        inner_absLimit_( tau_*outer_absLimit_)
    {
    }


    UzawaSolver(const MappingType& op,
                const InverseOperatorType& aufSolver,
                const DiscreteFunctionType& rhs,
                double redEps,
                double absLimit,
                int maxIter,
                int verbose=1
                )
      : op_(op), redEps_( redEps ), outer_absLimit_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_( verbose ),aufSolver_(aufSolver), bop_(op_.getBOP()),
        btop_(op_.getBTOP()),
        cop_(op_.getCOP()),
        rhs1_(rhs),
        rhs2_(op_.pressureRhs()),
        pressurespc_(op_.pressurespc()),
        spc_(op.spc()),
        velocity_("VELO",spc_),
        iter_(0),
        linIter_(0),
        tau_(0.2),
        inner_absLimit_( tau_*outer_absLimit_)
    {
    }


    /** \todo Please doc me! */
    virtual void operator()(const PressureDiscreteFunctionType& arg,
                            PressureDiscreteFunctionType& pressure ) const
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType Field;
       Field gamma=0, delta, rho;


      DiscreteFunctionType f("f",spc_);
      // f := rhs1
      f.assign(rhs1_);
      DiscreteFunctionType velocity("velocity",spc_);
      velocity.clear();
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
      bop_.apply(pressure,tmp1);
      // f -= tmp1
      f-=tmp1;
#if 0
      aufSolver_.set( inner_absLimit_ );
#endif
      // A^-1 * f = velocity
      aufSolver_.inv()(f,velocity);
      // B^T * velocity = tmp2
      btop_.apply(velocity,tmp2);
      //=> tmp2 = B^T * A^-1 * ( F - B * d )

      // residuum -= tmp2
      residuum-=tmp2;
      //=> residuum = arg - B^T * A^-1 * ( F - B * d )
      tmp2.clear();

      // C * pressure = tmp2
      cop_.apply(pressure, tmp2);
      // residuum += tmp2
      residuum += tmp2;

      // d := residuum;
      d.assign(residuum);

      // save iteration number
      linIter_+=aufSolver_.inv().iterations();

      // delta = (residuum,residuum)
      delta = residuum.scalarProductDofs( residuum );
      while((delta > outer_absLimit_) && (iter_+=1 < maxIter_))
      {
        tmp1.clear();
        // B * d = tmp1
        bop_.apply(d,tmp1);
#if 0
        inner_absLimit_ = std::max( delta/redEps_, tau_*outer_absLimit_ );
        //inner_absLimit_ = tau_ * std::min( 1.0, outer_absLimit_ / std::min( delta, 1.0 ) );
        //inner_absLimit_ = tau_ * delta;
        aufSolver_.set( inner_absLimit_ );
#endif
        // A^-1 * tmp1 = xi
        aufSolver_.inv()(tmp1,xi);
        // B^T * xi = h
        btop_.apply(xi,h);
        // => h = B^T * A^-1 * B * d
        tmp2.clear();
        cop_.apply( d, tmp2 );
        h += tmp2;

        // rho = delta / d.scalarProductDofs( h );
        rho    = delta / d.scalarProductDofs( h );
        // pressure -= rho * d
        pressure.axpy( -rho, d );
        // velocity += rho * xi
        velocity.axpy(rho,xi);
        // residuum -= rho * h
        residuum.axpy( -rho,h );

        double oldDelta = delta;

        // delta = (residuum,residuum)
        delta = residuum.scalarProductDofs( residuum );

        // save iteration number
        linIter_+=aufSolver_.inv().iterations();
        if( verbose_ > 0)
          std::cout << "SPcg-Iterationen " << iter_ << "   Residuum:"
                    << delta << "   lin. iter:" << aufSolver_.inv().iterations() <<std::endl;

        d *= delta / oldDelta;
        d += residuum;

      }
      if( verbose_ > 0)
        std::cout << std::endl;
      velocity_.assign(velocity);
    }

    DiscreteFunctionType& velocity()     {
      return velocity_;
    }

    int iterations(){return iter_;}
    double averageLinIter(){return linIter_/iter_;}

  private:
    // reference to operator which should be inverted
    const AssemblerType & op_;

    // reduce error each step by
    double redEps_;

    // minial error to reach
    typename DiscreteFunctionType::RangeFieldType outer_absLimit_;

    // number of maximal iterations
    int maxIter_;

    // level of output
    int verbose_;

    //the CGSolver for A^-1
    const InverseOperatorType& aufSolver_;

    const BOPType& bop_;
    const BTOPType& btop_;
    const COPType& cop_;

    const DiscreteFunctionType& rhs1_;
    const PressureDiscreteFunctionType& rhs2_;
    const PressureSpaceType& pressurespc_;
    const VeloSpaceType& spc_;
    mutable DiscreteFunctionType velocity_;
    mutable int iter_;
    mutable int linIter_;
    const double tau_;
    mutable double inner_absLimit_;
  };


} // end namespace Dune

#endif

