/***********************************************/
/**
* @file orbitAddVelocityAndAcceleration.cpp
*
* @brief Compute velocities and accelerations from a given orbit and save to one file.
*
* @author Norbert Zehentner
* @author Torsten Mayer-Guerr
* @date 2014-04-16
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program computes velocities and accelerations from a given \file{orbit}{instrument}
by differentiating a moving polynomial.
The values are saved in one output file which then contains orbit, velocity and acceleration.
)";

/***********************************************/

#include "programs/program.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Compute velocities and accelerations from a given orbit.
* @ingroup programsGroup */
class OrbitAddVelocityAndAcceleration
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(OrbitAddVelocityAndAcceleration, PARALLEL, "Compute velocities and accelerations from a given orbit", Orbit, Instrument)
GROOPS_RENAMED_PROGRAM(InstrumentOrbit2VelocityAcceleration, OrbitAddVelocityAndAcceleration, date2time(2020, 05, 25))

/***********************************************/

void OrbitAddVelocityAndAcceleration::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOut, fileNameIn;
    UInt     degree;

    readConfig(config, "outputfileOrbit",  fileNameOut, Config::MUSTSET, "", "");
    readConfig(config, "inputfileOrbit",   fileNameIn,  Config::MUSTSET, "", "");
    readConfig(config, "polynomialDegree", degree,      Config::DEFAULT, "8", "Polynomial degree, must be even!");
    if(isCreateSchema(config)) return;

    // ======================================================

    logStatus<<"read orbit data <"<<fileNameIn<<">"<<Log::endl;
    InstrumentFile orbitFile(fileNameIn);
    std::vector<Arc> arcList(orbitFile.arcCount());
    Parallel::forEach(arcList, [&] (UInt arcNo)
    {
      OrbitArc orbit = orbitFile.readArc(arcNo);
      UInt idx = 0;
      for(UInt idEpoch=0; idEpoch<orbit.size(); idEpoch++)
      {
        // find optimal interval
        // ---------------------
        while((idx+degree < orbit.size()) && (orbit.at(idx+degree).time < orbit.at(idEpoch).time))
          idx++;
        if(idx+degree >= orbit.size())
          break;

        UInt   idxOpt   = MAX_UINT;
        Double deltaOpt = 1e99;
        while((idx+degree < orbit.size()) && (orbit.at(idx).time <= orbit.at(idEpoch).time))
        {
          // interpolation point should be in the mid of the interval
          // => search minimum of the difference of the time before and after the interpolation point
          const Double delta = std::fabs(((orbit.at(idx+degree).time-orbit.at(idEpoch).time)-(orbit.at(idEpoch).time-orbit.at(idx).time)).seconds());
          if(delta <= deltaOpt)
          {
            idxOpt   = idx;
            deltaOpt = delta;
          }
          idx++;
        }
        idx = idxOpt;

        // polynomial interpolation——orbit
        // ------------------------
        // Matrix A(degree+1, degree+1);
        // for(UInt k=0; k<degree+1; k++)
        // {
        //   const Double factor = (orbit.at(idx+k).time-orbit.at(idEpoch).time).seconds();
        //   A(0,k) = 1.0;
        //   for(UInt n=1; n<=degree; n++)
        //     A(n,k) = factor * A(n-1,k);
        // }
        // Matrix coeff(degree+1, 2);
        // coeff(1, 0) = 1.; // velocity
        // coeff(2, 1) = 2.; // acceleration
        // solveInPlace(A, coeff);

        // orbit.at(idEpoch).velocity     = Vector3d();
        // orbit.at(idEpoch).acceleration = Vector3d();
        // for(UInt k=0; k<coeff.rows(); k++)
        // {
        //   orbit.at(idEpoch).velocity     += coeff(k,0) * orbit.at(idx+k).position;
        //   orbit.at(idEpoch).acceleration += coeff(k,1) * orbit.at(idx+k).position;
        // }

        // polynomial interpolation——velocity
        // ------------------------
        Matrix A(degree+1, degree+1);
        for(UInt k=0; k<degree+1; k++)
        {
          const Double factor = (orbit.at(idx+k).time-orbit.at(idEpoch).time).seconds();
          A(0,k) = 1.0;
          for(UInt n=1; n<=degree; n++)
            A(n,k) = factor * A(n-1,k);
        }
        Matrix coeff(degree+1, 1);
        coeff(1, 0) = 1.; // acceleration
        solveInPlace(A, coeff);

        orbit.at(idEpoch).acceleration = Vector3d();
        for(UInt k=0; k<coeff.rows(); k++)
        {
          orbit.at(idEpoch).acceleration += coeff(k,0) * orbit.at(idx+k).velocity;
        }

        // ----- �����������ղ�ֵ + ��ֵ΢�� -----
        // ׼�������ڵ�ʱ���λ��
        // const double dt = 0.025;  // half step for central difference (total step 0.05 s)
        // Time dt_time = seconds2time(dt);
        
        // std::vector<Time>   times(degree+1);
        // std::vector<Vector3d> pos(degree+1);
        // for(UInt k=0; k<=degree; ++k) {
        //   times[k] = orbit.at(idx+k).time;
        //   pos[k]   = orbit.at(idx+k).position;
        // }
        // // �����������ղ�ֵ�����������ڸ���ʱ�� t ��λ�ã�
        // auto lagrangeInterp = [&](const Time& t) -> Vector3d {
        //   Vector3d result(0,0,0);
        //   for(UInt i=0; i<=degree; ++i) {
        //     double L = 1.0;
        //     for(UInt j=0; j<=degree; ++j) {
        //       if(j==i) continue;
        //       L *= (t - times[j]).seconds() / (times[i] - times[j]).seconds();
        //     }
        //     result += L * pos[i];
        //   }
        //   return result;
        // };
        // // �������Ĳ�����������ʱ�̵�λ��
        // Time t_curr = orbit.at(idEpoch).time;
        // Time t_minus = t_curr - dt_time;
        // Time t_plus  = t_curr + dt_time;
        // Vector3d r_curr = lagrangeInterp(t_curr);
        // Vector3d r_minus = lagrangeInterp(t_minus);
        // Vector3d r_plus  = lagrangeInterp(t_plus);
        // // ���Ĳ�ּ����ٶȺͼ��ٶ�
        // orbit.at(idEpoch).velocity     = (r_plus - r_minus) / (2.0 * dt);
        // orbit.at(idEpoch).acceleration = (r_plus - 2.0 * r_curr + r_minus) / (dt * dt);
      }
      return orbit;
    }, comm);

    // write results
    // -------------
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
