#ifndef DUNE_UZAWASOLVER_HH
#define DUNE_UZAWASOLVER_HH


#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>


#include "../common/mymatrix.hh"
#include "../common/spmatextra.hh"

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
      : op_(op), _redEps ( redEps ), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , _verbose ( verbose ),aufSolver_(aufSolver), bop_(op_.getBOP()),
				btop_(op_.getBTOP()),
				cop_(op_.getCOP()),
				rhs1_(aufSolver_.affineShift()),
				rhs2_(op_.veloRhs()),
				pressurespc_(op_.pressurespc()),
				spc_(op.spc()),
				velocity_("VELO",spc_),
				iter_(0),
				linIter_(0)
		{

    }


    UzawaSolver(const MappingType& op,
								const InverseOperatorType& aufSolver,
								DiscreteFunctionType& rhs,
								double redEps,
								double absLimit,
								int maxIter,
								int verbose=1
								)
      : op_(op), _redEps ( redEps ), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , _verbose ( verbose ),aufSolver_(aufSolver), bop_(op_.getBOP()),
				btop_(op_.getBTOP()),
				cop_(op_.getCOP()),
				rhs1_(rhs),
				rhs2_(op_.veloRhs()),
				pressurespc_(op_.pressurespc()),
				spc_(op.spc()),
				velocity_("VELO",spc_),
				iter_(0),
				linIter_(0)
		{

    }


    /** \todo Please doc me! */
    virtual void operator()(const PressureDiscreteFunctionType& arg,
                            PressureDiscreteFunctionType& dest ) const
    {
			typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeFieldType Field;
      typedef typename PressureDiscreteFunctionType::DiscreteFunctionSpaceType PressureFunctionSpaceType;
       Field spa=0, spn, q, quad;


      DiscreteFunctionType f("f",spc_);
		  f.assign(rhs1_);
		  DiscreteFunctionType u("u",spc_);
      u.clear();
      DiscreteFunctionType tmp1("tmp1",spc_);

			tmp1.clear();
      DiscreteFunctionType xi("xi",spc_);
      xi.clear();
      PressureDiscreteFunctionType tmp2("tmp2",pressurespc_);
      tmp2.clear();



      //p<->d
      PressureDiscreteFunctionType p("p",pressurespc_);
      p.clear();
      PressureDiscreteFunctionType h("h",pressurespc_);
      h.clear();
      PressureDiscreteFunctionType g("g",pressurespc_);
      g.clear();
			PressureDiscreteFunctionType r("r",pressurespc_);
      r.assign(arg);

      bop_.apply(dest,tmp1);
			f-=tmp1;
			aufSolver_(f,u);
			linIter_+=aufSolver_.iterations();
		 	btop_.apply(u,tmp2);
		 	r-=tmp2;
    //   r.axpy(-1.,tmp2);
      tmp2.clear();

      p.assign(r);



      spn = r.scalarProductDofs( r );
      while((spn > epsilon_) && (iter_+=1 < maxIter_))
				{

					if(iter_ > 1)
						{
							const Field e = spn / spa;

							p *= e;
							p += r;
						}

					tmp1.clear();
					bop_.apply(p,tmp1);
					aufSolver_(tmp1,xi);
					linIter_+=aufSolver_.iterations();
					btop_.apply(xi,h);

					quad = p.scalarProductDofs( h );

					q    = spn / quad;

					dest.axpy( -q, p );
					u.axpy(q,xi);
					r.axpy( -q,h );

					spa = spn;


					spn = r.scalarProductDofs( r );

					if(_verbose > 0)
						std::cerr << " SPcg-Iterationen  " << iter_ << " Residuum:" << spn << "        \r";
				}
      if(_verbose > 0)
				std::cerr << "\n";
	    velocity_.assign(u);
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
    double _redEps;

    // minial error to reach
    typename DiscreteFunctionType::RangeFieldType epsilon_;

    // number of maximal iterations
    int maxIter_;

    // level of output
    int _verbose ;

    //the CGSolver for A^-1
    const InverseOperatorType& aufSolver_;

    const BOPType& bop_;
    const BTOPType& btop_;
    const COPType& cop_;

    const DiscreteFunctionType& rhs1_;
    const DiscreteFunctionType& rhs2_;
    const PressureSpaceType& pressurespc_;
    const VeloSpaceType& spc_;
    mutable DiscreteFunctionType velocity_;
    mutable int iter_;
		mutable int linIter_;
  };


} // end namespace Dune

#endif

