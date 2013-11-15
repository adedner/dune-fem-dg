#ifndef DUNE_EULERMODEL_HH
#define DUNE_EULERMODEL_HH

// system includes
#include <config.h>
#include <cmath>

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/space/common/functionspace.hh>

#include "../navierstokes/thermodynamics.hh"

#include <dune/fem-dg/models/defaultmodel.hh>
#include <dune/fem-dg/operator/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/rotator.hh>
#include <dune/fem-dg/operator/limiter/limitpass.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart >
    class EulerModelTraits 
    {
     public:
      typedef GridPart GridPartType;
      typedef typename GridPart::GridType                           GridType;
      enum{ dimDomain = GridType::dimensionworld };
      enum{ dimRange = dimDomain + 2 }; // the Euler equations
      enum{ dimGrad = dimRange * dimDomain };

      typedef typename GridType :: ctype  ctype ;

      typedef FunctionSpace< ctype, ctype, dimDomain, dimRange > FunctionSpaceType ;
      typedef typename ToNewDimRangeFunctionSpace< FunctionSpaceType, dimGrad > :: Type  GradientFunctionSpaceType;

      typedef typename FunctionSpaceType :: DomainType         DomainType;
      typedef typename FunctionSpaceType :: RangeType          RangeType;
      typedef typename FunctionSpaceType :: JacobianRangeType  JacobianRangeType;
      typedef JacobianRangeType                                FluxRangeType ; 
      typedef typename FunctionSpaceType :: RangeFieldType     FieldType ;
      typedef typename ToNewDimDomainFunctionSpace< 
        FunctionSpaceType, dimDomain-1 > :: Type :: DomainType    FaceDomainType ;

      typedef typename GradientFunctionSpaceType :: RangeType          GradientType;
      typedef GradientType                                        DiffusionType;
      typedef typename GradientFunctionSpaceType :: JacobianRangeType  JacobianFluxRangeType ;

      typedef typename GridPart::IntersectionIteratorType           IntersectionIteratorType;
      typedef typename IntersectionIteratorType::Intersection       IntersectionType;
      typedef typename GridType::template Codim<0>::Entity          EntityType;
      typedef typename GridType::template Codim<0>::EntityPointer   EntityPointerType;

      typedef Thermodynamics< dimDomain >                           ThermodynamicsType;

      typedef MinModLimiter< FieldType > LimiterFunctionType ;
      //typedef SuperBeeLimiter< FieldType > LimiterFunctionType ;
      //typedef VanLeerLimiter< FieldType > LimiterFunctionType ;
    };


    // Euler equations for dry atmosphere
    template< class GridPartType , class ProblemImp >
    class EulerModel :
      public DefaultModel< EulerModelTraits< GridPartType > >
    {
      typedef EulerModel< GridPartType, ProblemImp >            ThisType;
     public:
      typedef EulerModelTraits< GridPartType >                  Traits;
      typedef ProblemImp                                        ProblemType;

      enum { dimDomain = Traits :: dimDomain };
      enum { dimRange = Traits :: dimRange };

      typedef typename Traits :: GridType                       GridType;
      typedef typename Traits :: EntityType                     EntityType;
      typedef typename Traits :: EntityPointerType              EntityPointerType;
      typedef typename Traits :: IntersectionIteratorType       IntersectionIteratorType;
      typedef typename Traits :: IntersectionType               IntersectionType;
      typedef typename Traits :: FaceDomainType                 FaceDomainType;

      typedef typename Traits :: RangeType                      RangeType;
      typedef typename RangeType :: field_type                  FieldType ;
      typedef typename Traits :: DomainType                     DomainType;
      typedef typename Traits :: FluxRangeType                  FluxRangeType;
      typedef typename Traits :: GradientType              GradientType;
      typedef typename Traits :: JacobianRangeType              JacobianRangeType;
      typedef typename Traits :: JacobianFluxRangeType          JacobianFluxRangeType;

      typedef typename Traits::ThermodynamicsType               ThermodynamicsType;

     public:
      EulerModel( const ProblemType& problem )
        : gamma_( problem.gamma() )
        , problem_( problem )
        , fieldRotator_( 1 ) // insert number of fist velocity component   
      {}

      double gamma () const { return gamma_; }

      inline bool hasStiffSource() const { return false; }
      inline bool hasNonStiffSource() const { return false; }
      inline bool hasFlux() const { return true ; }
      inline double stiffSource( const EntityType& en,
                            const double time,
                            const DomainType& x,
                            const RangeType& u,
                            const GradientType& du,
                            RangeType & s) const
      {
        return stiffSource( en, time, x, u, s );
      }


      inline double stiffSource( const EntityType& en,
                            const double time,
                            const DomainType& x,
                            const RangeType& u,
                            const JacobianRangeType& jac,
                            RangeType & s) const
      {
        return stiffSource( en, time, x, u, s ); 
      }


      inline double stiffSource( const EntityType& en
                          , const double time
                          , const DomainType& x
                          , const RangeType& u
                          , RangeType& s ) const 
      {
        s = 0;
        return 0; 
      }


      inline double nonStiffSource( const EntityType& en,
                            const double time,
                            const DomainType& x,
                            const RangeType& u,
                            const GradientType& du,
                            RangeType & s) const
      {
        FieldMatrixConverter< GradientType, JacobianRangeType > jac( du );
        return nonStiffSource( en, time, x, u, jac, s );
      }



      template< class JacobianRangeTypeImp >
      inline double nonStiffSource( const EntityType& en,
                            const double time,
                            const DomainType& x,
                            const RangeType& u,
                            const JacobianRangeTypeImp& jac,
                            RangeType& s) const
      {
        s = 0;
        return 0; 
      }

      inline double pressure( const RangeType& u ) const 
      {
        return EulerAnalyticalFlux< dimDomain >().pressure( gamma_ , u );
      }

      inline void conservativeToPrimitive( const double time, 
                                           const DomainType& xgl,
                                           const RangeType& cons, 
                                           RangeType& prim,
                                           const bool ) const
      {
        problem_.evaluate( xgl, time, prim );
        //thermodynamics_.conservativeToPrimitiveEnergyForm( cons, prim );
      }

      inline void advection( const EntityType& en,
                             const double time,
                             const DomainType& x,
                             const RangeType& u,
                             FluxRangeType& f ) const 
      {
        EulerAnalyticalFlux<dimDomain>().analyticalFlux( gamma_ , u , f );
      }

      inline void eigenValues(const EntityType& en,
                              const double time,
                              const DomainType& x,
                              const RangeType& u,
                              RangeType& maxValue) const
      {
        std::cerr <<"eigenValues for problems/euler not implemented\n";
        abort();
      }

      inline double diffusionTimeStep( const IntersectionType& it,
                                       const double enVolume,
                                       const double circumEstimate,
                                       const double time,
                                       const FaceDomainType& x,
                                       const RangeType& u ) const
      {
        return 0;
      }

      // is not used
      inline  void jacobian( const EntityType& en,
                             const double time,
                             const DomainType& x,
                             const RangeType& u,
                             const FluxRangeType& du,
                             RangeType& A ) const
      {
        EulerAnalyticalFlux<dimDomain>().jacobian( gamma_ , u , du , A );
      }
      

      enum { Inflow = 1, Outflow = 2, Reflection = 3 , Slip = 4 };
      enum { MaxBnd = Slip };

      inline bool hasBoundaryValue( const IntersectionType& it,
                                    const double time,
                                    const FaceDomainType& x ) const
      { 
        const int bndId = problem_.boundaryId( it.boundaryId() );
        // on slip boundary we use boundaryFlux 
        return bndId != Slip;
      }

      // return iRight for insertion into the numerical flux 
      inline void boundaryValue( const IntersectionType& it,
                                 const double time,
                                 const FaceDomainType& x,
                                 const RangeType& uLeft,
                                 RangeType& uRight ) const
      {
        // Neumann boundary condition
        //uRight = uLeft;
        // 5 and 6 is also Reflection 
        //const int bndId = (it.boundaryId() > MaxBnd) ? MaxBnd : it.boundaryId();
        const int bndId = problem_.boundaryId( it.boundaryId() );

        assert( bndId > 0 );
        if( bndId == Inflow )
        {
          const DomainType xgl = it.geometry().global(x);
          problem_.evaluate(xgl, time, uRight);
          return ;
        }
        else if ( bndId == Outflow )
        {
          uRight = uLeft;
          return ;
        }
        else if ( bndId == Reflection )
        {
          uRight = uLeft;
          const DomainType unitNormal = it.unitOuterNormal( x );
          fieldRotator_.rotateForth( uRight , unitNormal );
          // Specific for euler: opposite sign for first component of momentum
          uRight[1] = -uRight[1];
          fieldRotator_.rotateBack( uRight, unitNormal );
          return ;
        }
        /*
        else if ( bndId == Slip )
        {
          RangeType tmp(uLeft);
          const DomainType unitNormal = it.unitOuterNormal(x);
          this->rot_.rotateForth( tmp , unitNormal );
          tmp[1] = 0.;
          tmp[2] /= tmp[0];
          tmp[3] = pressure(uLeft);
          prim2cons(tmp,uRight);
          this->rot_.rotateBack( uRight, unitNormal );
          return ;
        }
        */
        else
        {
          uRight = uLeft;
          return ;
          assert( false );
          abort();
        }
      }
     
      // boundary condition here is slip boundary cond. <u,n>=0
      // gLeft= p*[0 n(global(x)) 0]
      inline double boundaryFlux( const IntersectionType& it,
                                  const double time,
                                  const FaceDomainType& x,
                                  const RangeType& uLeft,
                                  RangeType& gLeft ) const
      {
        // Slip boundary condition 
        const DomainType normal = it.integrationOuterNormal( x ); 
        
        const double p = EulerAnalyticalFlux< dimDomain >().pressure( gamma_ , uLeft );
        gLeft = 0;
        for ( int i = 0 ; i < dimDomain ; ++i ) 
          gLeft[i+1] = normal[i] * p;
        return 0.;
      }

      void diffusion( const EntityType& en,
                      const double time,
                      const DomainType& x,
                      const RangeType& u,
                      const GradientType& v,
                      JacobianRangeType& diff ) const
      {
      }


      template <class JacobianRangeImp>
      void diffusion( const EntityType& en,
                      const double time,
                      const DomainType& x,
                      const RangeType& u,
                      const JacobianRangeImp& jac,
                      JacobianRangeType& diff ) const
      {
      }



      inline double diffusionBoundaryFlux( const IntersectionType& it,
                                           const double time,
                                           const FaceDomainType& x,
                                           const RangeType& uLeft,
                                           const GradientType& gradLeft,
                                           RangeType& gLeft ) const  
      {
        FieldMatrixConverter< GradientType, JacobianRangeType> jacLeft( gradLeft );
        return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
      }

      /** \brief boundary flux for the diffusion part
       */
      template< class JacobianRangeImp >
      inline double diffusionBoundaryFlux( const IntersectionType& it,
                                           const double time,
                                           const FaceDomainType& x,
                                           const RangeType& uLeft,
                                           const JacobianRangeImp& jacLeft,
                                           RangeType& gLeft ) const  
      {
        return 0.;
      }

      // here x is in global coordinates
      inline void maxSpeed( const EntityType& entity,
                            const double time,
                            const DomainType& x,
                            const DomainType& normal,
                            const RangeType& u,
                            double& advspeed,
                            double& totalspeed ) const
      {
        advspeed = EulerAnalyticalFlux< dimDomain >().maxSpeed( gamma_ , normal , u );
        totalspeed = advspeed;
      }

      inline const ProblemType& problem() const 
      {
        return problem_;
      }

      const double gamma_;

      /////////////////////////////////////////////////////////////////
      // Limiter section 
      ////////////////////////////////////////////////////////////////
      inline void velocity(
                 const EntityType& en,
                 const double time,
                 const DomainType& x,
                 const RangeType& u,
                 DomainType& velocity) const
      {
        for(int i=0; i<dimDomain; ++i)
        {
          // we store \rho u but do not need to divide by \rho here since only
          // sign is needed.
          velocity[i] = u[i+1];
        }
      }

      // we have physical check for this model 
      bool hasPhysical() const
      {
        return true;
      }

      // calculate jump between left and right value 
      inline bool physical(const EntityType& entity, 
                           const DomainType& xGlobal, 
                           const RangeType& u) const
      {
        if (u[0]<1e-8)
          return false;
        else
        {
          //std::cout << EulerAnalyticalFlux<dimDomain>().rhoeps(u) << std::endl;
          return (EulerAnalyticalFlux<dimDomain>().rhoeps(u) > 1e-8);
        }
      }

      // adjust average value if necessary 
      // (e.g. transform from conservative to primitive variables )
      void adjustAverageValue( const EntityType& entity,
                               const DomainType& xLocal,
                               RangeType& u ) const
      {
        // nothing to be done here for this test case 
      }

      // calculate jump between left and right value 
      inline void jump(const IntersectionType& it,
                       const double time, 
                       const FaceDomainType& x,
                       const RangeType& uLeft,
                       const RangeType& uRight,
                       RangeType& jump) const
      {
        // take pressure as shock detection values 
        const FieldType pl = pressure( uLeft );
        const FieldType pr = pressure( uRight );
        jump  = (pl-pr)/(0.5*(pl+pr));
      }

      // calculate jump between left and right value 
      inline void adaptationIndicator(
                       const IntersectionType& it,
                       const double time, 
                       const FaceDomainType& x,
                       const RangeType& uLeft,
                       const RangeType& uRight,
                       RangeType& indicator) const
      {
        // take density as shock detection values 
        indicator = (uLeft[0] - uRight[0])/(0.5 * (uLeft[0]+uRight[0]));

        const DomainType unitNormal = it.unitOuterNormal(x);
        RangeType ul(uLeft),ur(uRight);
        fieldRotator_.rotateForth( ul , unitNormal );
        fieldRotator_.rotateForth( ur , unitNormal );
        for (int i=1; i<dimDomain+1; ++i)
        {
          indicator += ((ul[i+1]/ul[0] - ur[i+1]/ur[0])/problem_.V());
        }
      }

     protected:
      const ThermodynamicsType thermodynamics_;
      const ProblemType& problem_;
      EulerFluxes::FieldRotator< DomainType, RangeType > fieldRotator_;
    };

  }

}

#endif
