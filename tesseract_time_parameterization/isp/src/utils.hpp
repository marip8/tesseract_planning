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

#ifndef TESSERACT_TIME_PARAMETERIZATION_ISP_UTILS_HPP
#define TESSERACT_TIME_PARAMETERIZATION_ISP_UTILS_HPP

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <vector>
#include <limits>
#include <cmath>
TESSERACT_COMMON_IGNORE_WARNINGS_POP
#include <tesseract_common/utils.h>

namespace tesseract_planning
{
//////// Internal functions //////////////

/*
  Fit a 'clamped' cubic spline over a series of points.
  A cubic spline ensures continuous function across positions,
  1st derivative (velocities), and 2nd derivative (accelerations).
  'Clamped' means the first derivative at the endpoints is specified.

  Fitting a cubic spline involves solving a series of linear equations.
  The general form for each segment is:
    (tj-t_(j-1))*x"_(j-1) + 2*(t_(j+1)-t_(j-1))*x"j + (t_(j+1)-tj)*x"_j+1) =
          (x_(j+1)-xj)/(t_(j+1)-tj) - (xj-x_(j-1))/(tj-t_(j-1))
  And the first and last segment equations are clamped to specified values: x1_i and x1_f.

  Represented in matrix form:
  [ 2*(t1-t0)   (t1-t0)                              0              ][x0"]       [(x1-x0)/(t1-t0) - t1_i           ]
  [ t1-t0       2*(t2-t0)   t2-t1                                   ][x1"]       [(x2-x1)/(t2-t1) - (x1-x0)/(t1-t0)]
  [             t2-t1       2*(t3-t1)   t3-t2                       ][x2"] = 6 * [(x3-x2)/(t3/t2) - (x2-x1)/(t2-t1)]
  [                       ...         ...         ...               ][...]       [...                              ]
  [ 0                                    tN-t_(N-1)  2*(tN-t_(N-1)) ][xN"]       [t1_f - (xN-x_(N-1))/(tN-t_(N-1)) ]

  This matrix is tridiagonal, which can be solved solved in O(N) time
  using the tridiagonal algorithm.
  There is a forward propogation pass followed by a backsubstitution pass.

  n is the number of points
  dt contains the time difference between each point (size=n-1)
  x  contains the positions                          (size=n)
  x1 contains the 1st derivative (velocities)        (size=n)
     x1[0] and x1[n-1] MUST be specified.
  x2 contains the 2nd derivative (accelerations)     (size=n)
  x1 and x2 are filled in by the algorithm.
*/
// NOLINTNEXTLINE
static void fit_cubic_spline(const long n, const double dt[], const double x[], double x1[], double x2[])
{
  long i{ 0 };
  const double x1_i = x1[0], x1_f = x1[n - 1];

  // Tridiagonal alg - forward sweep
  // x1 and x2 used to store the temporary coefficients c and d
  // (will get overwritten during backsubstitution)
  double *c = x1, *d = x2;
  c[0] = 0.5;
  d[0] = 3.0 * ((x[1] - x[0]) / dt[0] - x1_i) / dt[0];
  for (i = 1; i <= n - 2; i++)
  {
    const double dt2 = dt[i - 1] + dt[i];
    const double a = dt[i - 1] / dt2;
    const double denom = 2.0 - a * c[i - 1];
    c[i] = (1.0 - a) / denom;
    d[i] = 6.0 * ((x[i + 1] - x[i]) / dt[i] - (x[i] - x[i - 1]) / dt[i - 1]) / dt2;
    d[i] = (d[i] - a * d[i - 1]) / denom;
  }
  const double denom = dt[n - 2] * (2.0 - c[n - 2]);
  d[n - 1] = 6.0 * (x1_f - (x[n - 1] - x[n - 2]) / dt[n - 2]);
  d[n - 1] = (d[n - 1] - dt[n - 2] * d[n - 2]) / denom;

  // Tridiagonal alg - backsubstitution sweep
  // 2nd derivative
  x2[n - 1] = d[n - 1];
  for (i = n - 2; i >= 0; i--)
    x2[i] = d[i] - c[i] * x2[i + 1];

  // 1st derivative
  x1[0] = x1_i;
  for (i = 1; i < n - 1; i++)
    x1[i] = (x[i + 1] - x[i]) / dt[i] - (2 * x2[i] + x2[i + 1]) * dt[i] / 6.0;
  x1[n - 1] = x1_f;
}

/*
  Modify the value of x[1] and x[N-2]
  so that 2nd derivative starts and ends at specified value.
  This involves fitting the spline twice,
  then solving for the specified value.

  x2_i and x2_f are the (initial and final) 2nd derivative at 0 and N-1
*/

static void adjust_two_positions(const long n,
                                 const double dt[],  // NOLINT
                                 double x[],         // NOLINT
                                 double x1[],        // NOLINT
                                 double x2[],        // NOLINT
                                 const double x2_i,
                                 const double x2_f)
{
  x[1] = x[0];
  x[n - 2] = x[n - 3];
  fit_cubic_spline(n, dt, x, x1, x2);
  double a0 = x2[0];
  double b0 = x2[n - 1];

  x[1] = x[2];
  x[n - 2] = x[n - 1];
  fit_cubic_spline(n, dt, x, x1, x2);
  double a2 = x2[0];
  double b2 = x2[n - 1];

  // we can solve this with linear equation (use two-point form)
  // if (a2 != a0)
  if (!tesseract_common::almostEqualRelativeAndAbs(a2, a0, 1e-5))
    x[1] = x[0] + ((x[2] - x[0]) / (a2 - a0)) * (x2_i - a0);

  // if (b2 != b0)
  if (!tesseract_common::almostEqualRelativeAndAbs(b2, b0, 1e-5))
    x[n - 2] = x[n - 3] + ((x[n - 1] - x[n - 3]) / (b2 - b0)) * (x2_f - b0);
}

/*
  Find time required to go max velocity on each segment.
  Increase a segment's time interval if the current time isn't long enough.
*/
// NOLINTNEXTLINE
static void init_times(long n, double dt[], const double x[], const double max_velocity[], double min_velocity[])
{
  for (long i = 0; i < n - 1; i++)
  {
    double time{ NAN };
    double dx = x[i + 1] - x[i];
    if (dx >= 0.0)
      time = (dx / max_velocity[i]);
    else
      time = (dx / min_velocity[i]);
    time += std::numeric_limits<double>::epsilon();  // prevent divide-by-zero

    if (dt[i] < time)
      dt[i] = time;
  }
}

/*
  Fit a spline, then check each interval to see if bounds are met.
  If all bounds met (no time adjustments made), return 0.
  If bounds not met (time adjustments made), slightly increase the
  surrounding time intervals and return 1.

  n is the number of points
  dt contains the time difference between each point (size=n-1)
  x  contains the positions                          (size=n)
  x1 contains the 1st derivative (velocities)        (size=n)
     x1[0] and x1[n-1] MUST be specified.
  x2 contains the 2nd derivative (accelerations)     (size=n)
  max_velocity is the max velocity for this joint.
  min_velocity is the min velocity for this joint.
  max_acceleration is the max acceleration for this joint.
  min_acceleration is the min acceleration for this joint.
  tfactor is the time adjustment (multiplication) factor.
  x1 and x2 are filled in by the algorithm.
*/

// static int fit_spline_and_adjust_times(const int n,
//                                       double dt[],
//                                       const double x[],
//                                       double x1[],
//                                       double x2[],
//                                       const double max_velocity,
//                                       const double min_velocity,
//                                       const double max_acceleration,
//                                       const double min_acceleration,
//                                       const double tfactor)
//{
//  int i, ret = 0;

//  fit_cubic_spline(n, dt, x, x1, x2);

//  // Instantaneous velocity is calculated at each point
//  for (i = 0; i < n - 1; i++)
//  {
//    const double vel = x1[i];
//    const double vel2 = x1[i + 1];
//    if (vel > max_velocity || vel < min_velocity || vel2 > max_velocity || vel2 < min_velocity)
//    {
//      dt[i] *= tfactor;
//      ret = 1;
//    }
//  }
//  // Instantaneous acceleration is calculated at each point
//  if (ret == 0)
//  {
//    for (i = 0; i < n - 1; i++)
//    {
//      const double acc = x2[i];
//      const double acc2 = x2[i + 1];
//      if (acc > max_acceleration || acc < min_acceleration || acc2 > max_acceleration || acc2 < min_acceleration)
//      {
//        dt[i] *= tfactor;
//        ret = 1;
//      }
//    }
//  }

//  return ret;
//}

// return global expansion multiplicative factor required
// to force within bounds.
// Assumes that the spline is already fit
// (fit_cubic_spline must have been called before this).
static double global_adjustment_factor(long n,
                                       double /*dt*/[],                  // NOLINT
                                       const double /*x*/[],             // NOLINT
                                       double x1[],                      // NOLINT
                                       double x2[],                      // NOLINT
                                       const double max_velocity[],      // NOLINT
                                       const double min_velocity[],      // NOLINT
                                       const double max_acceleration[],  // NOLINT
                                       const double min_acceleration[])  // NOLINT
{
  double tfactor2 = 1.00;

  // fit_cubic_spline(n, dt, x, x1, x2);

  for (long i = 0; i < n; i++)
  {
    double tfactor{ NAN };
    tfactor = x1[i] / max_velocity[i];
    if (tfactor2 < tfactor)
      tfactor2 = tfactor;
    tfactor = x1[i] / min_velocity[i];
    if (tfactor2 < tfactor)
      tfactor2 = tfactor;

    if (x2[i] >= 0)
    {
      tfactor = sqrt(fabs(x2[i] / max_acceleration[i]));
      if (tfactor2 < tfactor)
        tfactor2 = tfactor;
    }
    else
    {
      tfactor = sqrt(fabs(x2[i] / min_acceleration[i]));
      if (tfactor2 < tfactor)
        tfactor2 = tfactor;
    }
  }

  return tfactor2;
}

// The path of a single joint: positions, velocities, and accelerations
struct SingleJointTrajectory
{
  std::vector<double> positions_;  // joint's position at time[x]
  std::vector<double> velocities_;
  std::vector<double> accelerations_;
  double initial_acceleration_{ 0 };
  double final_acceleration_{ 0 };
  std::vector<double> min_velocity_;
  std::vector<double> max_velocity_;
  std::vector<double> min_acceleration_;
  std::vector<double> max_acceleration_;
};

// Expands the entire trajectory to fit exactly within bounds
static void globalAdjustment(std::vector<SingleJointTrajectory>& t2,
                             long num_joints,
                             long num_points,
                             std::vector<double>& time_diff)
{
  double gtfactor = 1.0;
  for (std::size_t j = 0; j < static_cast<std::size_t>(num_joints); j++)
  {
    double tfactor = global_adjustment_factor(num_points,
                                              time_diff.data(),
                                              t2[j].positions_.data(),
                                              t2[j].velocities_.data(),
                                              t2[j].accelerations_.data(),
                                              t2[j].max_velocity_.data(),
                                              t2[j].min_velocity_.data(),
                                              t2[j].max_acceleration_.data(),
                                              t2[j].min_acceleration_.data());
    if (tfactor > gtfactor)
      gtfactor = tfactor;
  }

  // printf("# Global adjustment: %0.4f%%\n", 100.0 * (gtfactor - 1.0));
  for (std::size_t i = 0; i < static_cast<std::size_t>(num_points - 1); i++)
    time_diff[i] *= gtfactor;

  for (std::size_t j = 0; j < static_cast<std::size_t>(num_joints); j++)
  {
    fit_cubic_spline(
        num_points, time_diff.data(), t2[j].positions_.data(), t2[j].velocities_.data(), t2[j].accelerations_.data());
  }
}

}  // namespace tesseract_planning

#endif  // TESSERACT_TIME_PARAMETERIZATION_ISP_UTILS_HPP
