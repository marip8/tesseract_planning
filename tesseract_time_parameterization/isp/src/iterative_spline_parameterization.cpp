/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2017, Ken Anderson
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Ken Anderson */

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <console_bridge/console.h>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_time_parameterization/isp/iterative_spline_parameterization.h>
#include "utils.hpp"

/** @brief Clamps values on [0.0, 1.0], replacing invalid values with 1.0 */
static Eigen::VectorXd clamp(const Eigen::VectorXd& in)
{
  Eigen::VectorXd out = Eigen::VectorXd::Ones(in.size());

  // Set scaling factors
  if ((in.array() > 0.0).all() && (in.array() <= 1.0).all())
    out = in;
  else
  {
    for (Eigen::Index idx = 0; idx < in.size(); idx++)
    {
      if (std::abs(in[idx]) < std::numeric_limits<double>::epsilon())
      {
        // NOLINTNEXTLINE
        CONSOLE_BRIDGE_logDebug("ISP: A factor of 0.0 was "
                                "specified at index %d, "
                                "defaulting to %f instead.",
                                idx,
                                out[idx]);
      }
      else
      {
        // NOLINTNEXTLINE
        CONSOLE_BRIDGE_logWarn("ISP: Invalid factor %f specified at "
                               "index %d, "
                               "defaulting to %f instead.",
                               in[idx],
                               idx,
                               out[idx]);
      }
    }
  }

  return out;
}

namespace tesseract_planning
{
IterativeSplineParameterization::IterativeSplineParameterization(
    const Eigen::Ref<const Eigen::VectorXd>& max_velocity,
    const Eigen::Ref<const Eigen::VectorXd>& max_acceleration,
    const Eigen::Ref<const Eigen::VectorXd>& max_velocity_scaling_factors,
    const Eigen::Ref<const Eigen::VectorXd>& max_acceleration_scaling_factors,
    bool add_points)
  : max_velocity_(max_velocity)
  , max_acceleration_(max_acceleration)
  , max_velocity_scaling_factors_(clamp(max_velocity_scaling_factors))
  , max_acceleration_scaling_factors_(clamp(max_acceleration_scaling_factors))
  , add_points_(add_points)
{
}

IterativeSplineParameterization::IterativeSplineParameterization(double max_velocity,
                                                                 double max_acceleration,
                                                                 double max_velocity_scaling_factor,
                                                                 double max_acceleration_scaling_factor,
                                                                 bool add_points)
  : IterativeSplineParameterization(Eigen::VectorXd::Constant(6, max_velocity),
                                    Eigen::VectorXd::Constant(6, max_acceleration),
                                    Eigen::VectorXd::Constant(6, max_velocity_scaling_factor),
                                    Eigen::VectorXd::Constant(6, max_acceleration_scaling_factor),
                                    add_points)
{
}

IterativeSplineParameterization::IterativeSplineParameterization(const std::vector<double>& max_velocity,
                                                                 const std::vector<double>& max_acceleration,
                                                                 double max_velocity_scaling_factor,
                                                                 double max_acceleration_scaling_factor,
                                                                 bool add_points)
  : IterativeSplineParameterization(
        Eigen::Map<const Eigen::VectorXd>(max_velocity.data(), static_cast<long>(max_velocity.size())),
        Eigen::Map<const Eigen::VectorXd>(max_acceleration.data(), static_cast<long>(max_acceleration.size())),
        Eigen::VectorXd::Constant(static_cast<long>(max_velocity.size()), max_velocity_scaling_factor),
        Eigen::VectorXd::Constant(static_cast<long>(max_acceleration.size()), max_acceleration_scaling_factor),
        add_points)
{
}

bool IterativeSplineParameterization::compute(TrajectoryContainer& trajectory) const
{
  if (trajectory.empty())
    return true;

  auto num_points = static_cast<std::size_t>(trajectory.size());

  if (max_velocity_.size() != trajectory.dof() || max_acceleration_.size() != trajectory.dof())
    return false;

  // JointTrajectory indexes in [point][joint] order.
  // We need [joint][point] order to solve efficiently,
  // so convert form here.

  std::vector<SingleJointTrajectory> t2(static_cast<std::size_t>(trajectory.dof()));

  const Eigen::VectorXd& start_vel = trajectory.getVelocity(0);
  const Eigen::VectorXd& last_vel = trajectory.getVelocity(static_cast<Eigen::Index>(num_points - 1));
  const Eigen::VectorXd& start_acc = trajectory.getAcceleration(0);
  const Eigen::VectorXd& last_acc = trajectory.getAcceleration(static_cast<Eigen::Index>(num_points - 1));

  for (std::size_t j = 0; j < static_cast<std::size_t>(trajectory.dof()); j++)
  {
    // Copy positions
    t2[j].positions_.resize(num_points, 0.0);
    for (std::size_t i = 0; i < num_points; i++)
      t2[j].positions_[i] = trajectory.getPosition(static_cast<Eigen::Index>(i))[static_cast<Eigen::Index>(j)];

    // Initialize velocities
    t2[j].velocities_.resize(num_points, 0.0);

    // Copy initial/final velocities if specified
    if (start_vel.size() > 0)
      t2[j].velocities_[0] = start_vel[static_cast<Eigen::Index>(j)];
    if (last_vel.size() > 0)
      t2[j].velocities_[num_points - 1] = last_vel[static_cast<Eigen::Index>(j)];

    // Initialize accelerations
    t2[j].accelerations_.resize(num_points, 0.0);
    t2[j].initial_acceleration_ = 0.0;
    t2[j].final_acceleration_ = 0.0;

    // Copy initial/final accelerations if specified
    if (start_acc.size() > 0)
      t2[j].initial_acceleration_ = start_acc[static_cast<Eigen::Index>(j)];
    t2[j].accelerations_[0] = t2[j].initial_acceleration_;
    if (last_acc.size() > 0)
      t2[j].final_acceleration_ = last_acc[static_cast<Eigen::Index>(j)];
    t2[j].accelerations_[num_points - 1] = t2[j].final_acceleration_;

    // Set bounds based on inputs
    t2[j].max_velocity_.resize(num_points);
    t2[j].min_velocity_.resize(num_points);
    t2[j].max_acceleration_.resize(num_points);
    t2[j].min_acceleration_.resize(num_points);
    Eigen::Map<Eigen::VectorXd> max_velocity_eigen(t2[j].max_velocity_.data(), static_cast<Eigen::Index>(num_points));
    Eigen::Map<Eigen::VectorXd> min_velocity_eigen(t2[j].min_velocity_.data(), static_cast<Eigen::Index>(num_points));
    Eigen::Map<Eigen::VectorXd> max_acceleration_eigen(t2[j].max_acceleration_.data(),
                                                       static_cast<Eigen::Index>(num_points));
    Eigen::Map<Eigen::VectorXd> min_acceleration_eigen(t2[j].min_acceleration_.data(),
                                                       static_cast<Eigen::Index>(num_points));

    max_velocity_eigen =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(num_points)) * max_velocity_[static_cast<Eigen::Index>(j)];
    min_velocity_eigen =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(num_points)) * -max_velocity_[static_cast<Eigen::Index>(j)];
    max_velocity_eigen.array() *= max_velocity_scaling_factors_.array();
    min_velocity_eigen.array() *= max_velocity_scaling_factors_.array();

    max_acceleration_eigen =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(num_points)) * max_acceleration_[static_cast<Eigen::Index>(j)];
    min_acceleration_eigen =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(num_points)) * -max_acceleration_[static_cast<Eigen::Index>(j)];
    max_acceleration_eigen.array() *= max_acceleration_scaling_factors_.array();
    min_acceleration_eigen.array() *= max_acceleration_scaling_factors_.array();

    // Error out if bounds don't make sense
    if ((max_velocity_eigen.array() <= 0.0).any() || (max_acceleration_eigen.array() <= 0.0).any())
    {
      CONSOLE_BRIDGE_logError("iterative_spline_parameterization: Joint %d max velocity %f and max acceleration %f "
                              "must be greater than zero or a solution won't be found.",
                              j,
                              t2[j].max_velocity_[0],
                              t2[j].max_acceleration_[0]);
      return false;
    }
    if ((min_velocity_eigen.array() >= 0.0).any() || (min_acceleration_eigen.array() >= 0.0).any())
    {
      CONSOLE_BRIDGE_logError("trajectory_processing.iterative_spline_parameterization: Joint %d min velocity %f and "
                              "min acceleration %f must be less than zero or a solution won't be found.",
                              j,
                              t2[j].min_velocity_[0],
                              t2[j].min_acceleration_[0]);
      return false;
    }
  }

  bool add_points = add_points_;
  if (num_points < 2)
    add_points = false;

  if (add_points)
  {
    // Insert 2nd and 2nd-last points
    // (required to force acceleration to specified values at endpoints)
    for (unsigned int j = 0; j < trajectory.dof(); j++)
    {
      double value = 0.9 * t2[j].positions_[0] + 0.1 * t2[j].positions_[1];
      t2[j].positions_.insert(t2[j].positions_.begin() + 1, value);

      value = t2[j].velocities_.front();
      t2[j].velocities_.insert(t2[j].velocities_.begin() + 1, value);

      value = t2[j].accelerations_.front();
      t2[j].accelerations_.insert(t2[j].accelerations_.begin() + 1, value);

      value = t2[j].max_velocity_.front();
      t2[j].max_velocity_.insert(t2[j].max_velocity_.begin() + 1, value);

      value = t2[j].min_velocity_.front();
      t2[j].min_velocity_.insert(t2[j].min_velocity_.begin() + 1, value);

      value = t2[j].max_acceleration_.front();
      t2[j].max_acceleration_.insert(t2[j].max_acceleration_.begin() + 1, value);

      value = t2[j].min_acceleration_.front();
      t2[j].min_acceleration_.insert(t2[j].min_acceleration_.begin() + 1, value);
    }
    num_points++;

    for (unsigned int j = 0; j < trajectory.dof(); j++)
    {
      double value = 0.1 * t2[j].positions_[num_points - 2] + 0.9 * t2[j].positions_[num_points - 1];
      t2[j].positions_.insert(t2[j].positions_.end() - 1, value);

      value = t2[j].velocities_.back();
      t2[j].velocities_.insert(t2[j].velocities_.end() - 1, value);

      value = t2[j].accelerations_.back();
      t2[j].accelerations_.insert(t2[j].accelerations_.end() - 1, value);

      value = t2[j].max_velocity_.back();
      t2[j].max_velocity_.insert(t2[j].max_velocity_.end() - 1, value);

      value = t2[j].min_velocity_.back();
      t2[j].min_velocity_.insert(t2[j].min_velocity_.end() - 1, value);

      value = t2[j].max_acceleration_.back();
      t2[j].max_acceleration_.insert(t2[j].max_acceleration_.end() - 1, value);

      value = t2[j].min_acceleration_.back();
      t2[j].min_acceleration_.insert(t2[j].min_acceleration_.end() - 1, value);
    }
    num_points++;
  }

  // Error check
  if (num_points < 4)
  {
    CONSOLE_BRIDGE_logError("iterative_spline_parameterization: number of waypoints %d, needs to be greater than 3.",
                            num_points);
    return false;
  }
  for (std::size_t j = 0; j < static_cast<std::size_t>(trajectory.dof()); j++)
  {
    if (t2[j].velocities_[0] > t2[j].max_velocity_[0] || t2[j].velocities_[0] < t2[j].min_velocity_[0])
    {
      CONSOLE_BRIDGE_logError("iterative_spline_parameterization: Initial velocity %f out of bounds.",
                              t2[j].velocities_[0]);
      return false;
    }

    if (t2[j].velocities_[num_points - 1] > t2[j].max_velocity_[num_points - 1] ||
        t2[j].velocities_[num_points - 1] < t2[j].min_velocity_[num_points - 1])
    {
      CONSOLE_BRIDGE_logError("iterative_spline_parameterization: Final velocity %f out of bounds.",
                              t2[j].velocities_[num_points - 1]);
      return false;
    }

    if (t2[j].accelerations_[0] > t2[j].max_acceleration_[0] || t2[j].accelerations_[0] < t2[j].min_acceleration_[0])
    {
      CONSOLE_BRIDGE_logError("iterative_spline_parameterization: Initial acceleration %f out of bounds\n",
                              t2[j].accelerations_[0]);
      return false;
    }

    if (t2[j].accelerations_[num_points - 1] > t2[j].max_acceleration_[num_points - 1] ||
        t2[j].accelerations_[num_points - 1] < t2[j].min_acceleration_[num_points - 1])
    {
      CONSOLE_BRIDGE_logError("iterative_spline_parameterization: Final acceleration %f out of bounds\n",
                              t2[j].accelerations_[num_points - 1]);
      return false;
    }
  }

  // Initialize times
  // start with valid velocities, then expand intervals
  // epsilon to prevent divide-by-zero
  std::vector<double> time_diff(static_cast<std::size_t>(num_points - 1), std::numeric_limits<double>::epsilon());
  for (unsigned int j = 0; j < trajectory.dof(); j++)
    init_times(static_cast<long>(num_points),
               time_diff.data(),
               t2[j].positions_.data(),
               t2[j].max_velocity_.data(),
               t2[j].min_velocity_.data());

  // Stretch intervals until close to the bounds
  while (true)
  {
    int loop = 0;

    // Calculate the interval stretches due to acceleration
    std::vector<double> time_factor(static_cast<std::size_t>(num_points - 1), 1.00);
    for (unsigned int j = 0; j < trajectory.dof(); j++)
    {
      // Move points to satisfy initial/final acceleration
      if (add_points)
      {
        adjust_two_positions(static_cast<int>(num_points),
                             time_diff.data(),
                             t2[j].positions_.data(),
                             t2[j].velocities_.data(),
                             t2[j].accelerations_.data(),
                             t2[j].initial_acceleration_,
                             t2[j].final_acceleration_);
      }

      fit_cubic_spline(static_cast<int>(num_points),
                       time_diff.data(),
                       t2[j].positions_.data(),
                       t2[j].velocities_.data(),
                       t2[j].accelerations_.data());
      for (unsigned i = 0; i < num_points; i++)
      {
        const double acc = t2[j].accelerations_[i];
        double atfactor = 1.0;
        if (acc > t2[j].max_acceleration_[i])
          atfactor = sqrt(acc / t2[j].max_acceleration_[i]);
        if (acc < t2[j].min_acceleration_[i])
          atfactor = sqrt(acc / t2[j].min_acceleration_[i]);
        if (atfactor > 1.01)  // within 1%
          loop = 1;
        atfactor = (atfactor - 1.0) / 16.0 + 1.0;  // 1/16th
        if (i > 0)
          time_factor[i - 1] = std::max(time_factor[i - 1], atfactor);
        if (i < num_points - 1)
          time_factor[i] = std::max(time_factor[i], atfactor);
      }
    }

    if (loop == 0)
      break;  // finished

    // Stretch
    for (unsigned i = 0; i < num_points - 1; i++)
      time_diff[i] *= time_factor[i];
  }

  // Final adjustment forces the trajectory within bounds
  globalAdjustment(t2, trajectory.dof(), static_cast<long>(num_points), time_diff);

  // Convert back to JointTrajectory form
  double time = 0;
  Eigen::Index idx = 0;
  for (unsigned int i = 0; i < num_points; i++)
  {
    Eigen::VectorXd uv(trajectory.dof());
    Eigen::VectorXd ua(trajectory.dof());
    for (unsigned int j = 0; j < trajectory.dof(); j++)
    {
      uv[j] = t2[j].velocities_[i];
      ua[j] = t2[j].accelerations_[i];
    }

    // Calculate time from start
    if (i > 0)
      time = time + time_diff[i - 1];

    // Do not process added points
    if (add_points && (i == 1 || i == num_points - 2))
    {
      time = time + time_diff[i - 1];
      continue;
    }

    trajectory.setData(idx++, uv, ua, time);
  }

  assert(trajectory.isTimeStrictlyIncreasing());
  return true;
}

}  // namespace tesseract_planning
