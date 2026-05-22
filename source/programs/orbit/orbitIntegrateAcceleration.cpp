/***********************************************/
/**
* @file orbitIntegrateAcceleration.cpp
*
* @brief Integrate an orbit from an acceleration time series.
*
* @author OpenAI
* @date 2026-05-17
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program integrates an \file{orbit}{instrument} from externally given
\file{accelerometer}{instrument} accelerations in the celestial reference frame (CRF).
The epochs and the initial position and velocity are taken from \configFile{inputfileOrbit}{instrument}.

This is intended for consistency checks such as
\[
  \M r(t) \rightarrow \ddot{\M r}(t) \rightarrow \M r(t),
\]
where the acceleration was derived from a precise orbit by numerical differentiation.
The acceleration is integrated by a moving interpolation polynomial of degree
\config{integrationDegree}.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Integrate an orbit from an acceleration time series.
* @ingroup programsGroup */
class OrbitIntegrateAcceleration
{
  UInt                integrationDegree;
  std::vector<Vector> coeff_g, coeff_tg;

  Matrix integrate2Position(Double deltaT, const_MatrixSliceRef g) const;
  Matrix integrate2Velocity(Double deltaT, const_MatrixSliceRef g) const;

public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(OrbitIntegrateAcceleration, PARALLEL, "integrate orbit from acceleration time series", Orbit, Instrument)

/***********************************************/

void OrbitIntegrateAcceleration::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut, fileNameOrbit, fileNameAccelerometer;

    readConfig(config, "outputfileOrbit",        fileNameOut,           Config::MUSTSET, "",  "integrated orbit");
    readConfig(config, "inputfileOrbit",         fileNameOrbit,         Config::MUSTSET, "",  "reference orbit; epochs and initial state are used");
    readConfig(config, "inputfileAccelerometer", fileNameAccelerometer, Config::MUSTSET, "",  "accelerations in CRF [m/s^2]");
    readConfig(config, "integrationDegree",      integrationDegree,     Config::DEFAULT, "7", "integration by polynomial approximation of degree n");
    if(isCreateSchema(config)) return;

    if(integrationDegree%2 == 0)
      throw(Exception("polynomial degree for integration must be odd."));

    // Init integration coefficients, same formulation as PreprocessingVariationalEquation.
    coeff_g.resize (integrationDegree);
    coeff_tg.resize(integrationDegree);
    for(UInt idx=0; idx<coeff_g.size(); idx++)
    {
      Matrix W(integrationDegree+1, integrationDegree+1);
      for(UInt i=0; i<W.columns(); i++)
      {
        W(i,0) = 1.;
        for(UInt n=1; n<W.rows(); n++)
          W(i,n) = (Double(i)-idx) * W(i,n-1);
      }
      inverse(W);

      coeff_g.at(idx)  = Vector(W.rows());
      coeff_tg.at(idx) = Vector(W.rows());
      for(UInt i=0; i<W.rows(); i++)
        for(UInt n=0; n<W.columns(); n++)
        {
          coeff_g.at(idx)(i)  += 1./(n+1.) * W(n,i);
          coeff_tg.at(idx)(i) += 1./(n+2.) * W(n,i);
        }
    }

    InstrumentFile orbitFile(fileNameOrbit);
    InstrumentFile accelerometerFile(fileNameAccelerometer);
    InstrumentFile::checkArcCount({orbitFile, accelerometerFile});

    std::vector<Arc> arcList(orbitFile.arcCount());
    logStatus<<"integrate acceleration arcs"<<Log::endl;
    Parallel::forEach(arcList, [&] (UInt arcNo)
    {
      OrbitArc         orbit         = orbitFile.readArc(arcNo);
      AccelerometerArc accelerometer = accelerometerFile.readArc(arcNo);
      Arc::checkSynchronized({orbit, accelerometer});

      if(orbit.size() <= integrationDegree)
        throw(Exception("Epoch count must be greater than integrationDegree."));

      const std::vector<Time> times      = orbit.times();
      const UInt              epochCount = orbit.size();
      const Double            deltaT     = (times.back()-times.front()).seconds()/(epochCount-1);

      Vector g(3*epochCount);
      for(UInt i=0; i<epochCount; i++)
      {
        g(3*i+0) = accelerometer.at(i).acceleration.x();
        g(3*i+1) = accelerometer.at(i).acceleration.y();
        g(3*i+2) = accelerometer.at(i).acceleration.z();
      }

      Matrix posInt = integrate2Position(deltaT, g);
      Matrix velInt = integrate2Velocity(deltaT, g);

      OrbitArc orbitOut;
      for(UInt i=0; i<epochCount; i++)
      {
        OrbitEpoch epoch;
        epoch.time         = times.at(i);
        epoch.position     = orbit.front().position + (times.at(i)-times.front()).seconds() * orbit.front().velocity + Vector3d(posInt.row(3*i, 3));
        epoch.velocity     = orbit.front().velocity + Vector3d(velInt.row(3*i, 3));
        epoch.acceleration = accelerometer.at(i).acceleration;
        orbitOut.push_back(epoch);
      }

      return orbitOut;
    }, comm);

    if(Parallel::isMaster(comm))
    {
      logStatus<<"write orbit data to file <"<<fileNameOut<<">"<<Log::endl;
      InstrumentFile::write(fileNameOut, arcList);
      Arc::printStatistics(arcList);
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix OrbitIntegrateAcceleration::integrate2Position(Double deltaT, const_MatrixSliceRef g) const
{
  try
  {
    const UInt epochCount = g.rows()/3;
    Matrix result(g.rows(), g.columns());  // integral_0^t (t-t') f(t') dt'

    UInt idx       = 0;
    UInt evalPoint = 0;
    Matrix int_g(3, g.columns());   // integral_0^t f(t') dt'
    Matrix int_tg(3, g.columns());  // integral_0^t t'f(t') dt'
    for(UInt i=0; i<epochCount-1; i++)
    {
      Matrix tmp(3, g.columns());
      for(UInt k=0; k<coeff_g.at(idx).rows(); k++)
        axpy(deltaT * coeff_g.at(idx)(k), g.row(3*(i+k-evalPoint),3), tmp);
      for(UInt k=0; k<coeff_tg.at(idx).rows(); k++)
        axpy(deltaT*deltaT * coeff_tg.at(idx)(k), g.row(3*(i+k-evalPoint),3), int_tg);
      axpy(i*deltaT, tmp, int_tg);
      int_g += tmp;

      if(i<integrationDegree/2)
        idx++, evalPoint++;
      if(i>=epochCount-1-integrationDegree/2-1)
        idx++, evalPoint++;

      axpy((i+1)*deltaT, int_g,  result.row(3*(i+1),3));
      axpy(-1.,          int_tg, result.row(3*(i+1),3));
    }

    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Matrix OrbitIntegrateAcceleration::integrate2Velocity(Double deltaT, const_MatrixSliceRef g) const
{
  try
  {
    const UInt epochCount = g.rows()/3;
    Matrix result(g.rows(), g.columns());

    UInt idx       = 0;
    UInt evalPoint = 0;
    for(UInt i=0; i<epochCount-1; i++)
    {
      copy(result.row(3*i,3), result.row(3*(i+1),3));
      for(UInt k=0; k<coeff_g.at(idx).rows(); k++)
        axpy(deltaT*coeff_g.at(idx)(k), g.row(3*(i+k-evalPoint),3), result.row(3*(i+1),3));

      if(i<integrationDegree/2)
        idx++, evalPoint++;
      if(i>=epochCount-1-integrationDegree/2-1)
        idx++, evalPoint++;
    }

    return result;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
