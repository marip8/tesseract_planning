#ifndef TESSERACT_TIME_PARAMETERIZATION_CORE_TIME_PARAMETERIZATION_H
#define TESSERACT_TIME_PARAMETERIZATION_CORE_TIME_PARAMETERIZATION_H

#include <tesseract_time_parameterization/core/trajectory_container.h>

namespace tesseract_planning
{
/**
 * @brief Base class for performing time parameterization
 **/
struct TimeParameterization
{
  TimeParameterization() = default;
  virtual ~TimeParameterization() = default;

  /**
   * @brief Computes the time stamps for a trajectory
   * @return True if successful, otherwise false
   */
  virtual bool compute(TrajectoryContainer& trajectory) const = 0;
};

}  // namespace tesseract_planning

#endif  // TESSERACT_TIME_PARAMETERIZATION_CORE_TIME_PARAMETERIZATION_H
