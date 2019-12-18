#ifndef DUNE_FEMDG_MISC_TIMESTEPPINGPARAMETER_HH
#define DUNE_FEMDG_MISC_TIMESTEPPINGPARAMETER_HH

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

  /**
   *  \brief Parameter class for time stepping in instationary equations.
   *
   *  \ingroup ParameterClass
   */
  class TimeSteppingParameters
  : public Dune::Fem::LocalParameter< TimeSteppingParameters, TimeSteppingParameters >
  {
    protected:
    const std::string keyPrefix_;
    const Dune::Fem::ParameterReader parameter_;

    public:
    /**
     * \brief Constructor
     *
     * \param keyPrefix the key prefix for the parameter file.
     */
    TimeSteppingParameters( const std::string keyPrefix = "femdg.stepper.",
                            const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : keyPrefix_( keyPrefix ),
        parameter_( parameter )
    {}

    /**
     * \brief returns a fixed time step \f$ \Delta t=\text{const} \f$
     *
     * \note To choose the time step size adaptively (i.e. non fixed time steps),
     * set this value to \f$ \leq0 \f$.
     */
    virtual double fixedTimeStep() const
    {
      return parameter_.getValue< double >( keyPrefix_ + "fixedtimestep" , 0.0 );
    }

    /**
     * \brief return an additional scaling of fixedTimeStep() which is applied
     * for each eoc loop.
     *
     * Example: A value of \f$ 2 \f$ would half the fixed time step size for each eoc loop.
     */
    virtual double fixedTimeStepEocLoopFactor() const
    {
      return parameter_.getValue< double >( keyPrefix_ + "fixedtimestepeocloopfactor" , 1.0 );
    }

    /**
     * \brief returns the end time \f$ t_\text{start} \f$ of the time stepping.
     */
    virtual double startTime() const
    {
      return parameter_.getValue< double >( keyPrefix_ + "starttime" , 0.0 );
    }

    /**
     * \brief returns the end time \f$ t_\text{end} \f$ of the time stepping.
     */
    virtual double endTime() const
    {
      return parameter_.getValue< double >( keyPrefix_ + "endtime"/*, 1.0 */);
    }

    /**
     * \brief returns how often to print additional information about time stepping.
     *
     * \f$ 1 \f$ means every time step, \f$ 2 \f$ every second and so on.
     *
     * \note a value \f$ <1 \f$ disables printing.
     */
    virtual int printCount() const
    {
      return parameter_.getValue< int >( keyPrefix_ + "printcount" , -1 );
    }

    /**
     * \brief returns the maximal time step \f$ \Delta_\text{max}} t \f$ for the
     * time stepping.
     */
    virtual double maxTimeStep() const
    {
      return parameter_.getValue< double >( keyPrefix_ + "maxtimestep", std::numeric_limits<double>::max());
    }

    /**
     * \brief returns the maximal number of time steps.
     */
    virtual int maximalTimeSteps () const
    {
      return parameter_.getValue< int >(  keyPrefix_ + "maximaltimesteps", std::numeric_limits<int>::max());
    }

    /**
     * \brief returns true whether the last time step should reach \f$ t_\text{end} \f$ exactly or not.
     *
     * Usually, it is the case, that the last time step
     *
     * To avoid this inaccuracy (needed for eoc measurements etc.)
     * you can set this option to true.
     *
     * \note The last time steps may be time consuming because
     * the last time steps may become very small.
     */
    virtual bool stopAtEndTime() const
    {
      return parameter_.getValue< bool >( keyPrefix_ + "stopatendtime", bool(false) );
    }

  };
}}

#endif
