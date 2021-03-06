/**
 * @file lvs_interpolation.h
 * @brief
 *
 * @author Levi Armstrong
 * @author Matthew Powelson
 * @date July 23, 2020
 * @version TODO
 * @bug No known bugs
 *
 * @copyright Copyright (c) 2020, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (Apache License)
 * @par
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * @par
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <vector>
#include <memory>
#include <algorithm>  // std::max
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_motion_planners/core/types.h>
#include <tesseract_motion_planners/core/utils.h>
#include <tesseract_motion_planners/simple/step_generators/lvs_interpolation.h>

namespace tesseract_planning
{
CompositeInstruction LVSInterpolateStateWaypoint(const JointWaypoint& start,
                                                 const JointWaypoint& end,
                                                 const PlanInstruction& base_instruction,
                                                 const PlannerRequest& request,
                                                 const ManipulatorInfo& manip_info,
                                                 double state_longest_valid_segment_length,
                                                 double translation_longest_valid_segment_length,
                                                 double rotation_longest_valid_segment_length,
                                                 int min_steps)
{
  // Joint waypoints should have joint names
  assert(static_cast<long>(start.joint_names.size()) == start.size());
  assert(static_cast<long>(end.joint_names.size()) == end.size());

  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  CompositeInstruction composite;

  // Initialize
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  auto world_to_base = request.env_state->link_transforms.at(fwd_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  // Calculate FK for start and end
  Eigen::Isometry3d p1 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), start));
  if (!fwd_kin->calcFwdKin(p1, start))
    throw std::runtime_error("LVSInterpolateStateWaypoint: failed to find forward kinematics solution!");
  p1 = world_to_base * p1 * tcp;

  Eigen::Isometry3d p2 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), end));
  if (!fwd_kin->calcFwdKin(p2, end))
    throw std::runtime_error("LVSInterpolateStateWaypoint: failed to find forward kinematics solution!");
  p2 = world_to_base * p2 * tcp;

  double trans_dist = (p2.translation() - p1.translation()).norm();
  // maybe rotation instead of linear
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  double joint_dist = (end - start).norm();

  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int joint_steps = int(joint_dist / state_longest_valid_segment_length) + 1;

  int steps = std::max(trans_steps, rot_steps);
  steps = std::max(steps, joint_steps);
  steps = std::max(steps, min_steps);

  // Linearly interpolate in joint space
  Eigen::MatrixXd states = interpolate(start, end, steps);

  // Get move type base on base instruction type
  MoveInstructionType move_type;
  if (base_instruction.isLinear())
    move_type = MoveInstructionType::LINEAR;
  else if (base_instruction.isFreespace())
    move_type = MoveInstructionType::FREESPACE;
  else
    throw std::runtime_error("LVS Interpolation: Unsupported Move Instruction Type!");

  // Convert to MoveInstructions
  for (long i = 1; i < states.cols(); ++i)
  {
    MoveInstruction move_instruction(StateWaypoint(fwd_kin->getJointNames(), states.col(i)), move_type);
    move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
    move_instruction.setDescription(base_instruction.getDescription());
    move_instruction.setProfile(base_instruction.getProfile());
    composite.push_back(move_instruction);
  }

  return composite;
}

CompositeInstruction LVSInterpolateStateWaypoint(const JointWaypoint& start,
                                                 const CartesianWaypoint& end,
                                                 const PlanInstruction& base_instruction,
                                                 const PlannerRequest& request,
                                                 const ManipulatorInfo& manip_info,
                                                 double state_longest_valid_segment_length,
                                                 double translation_longest_valid_segment_length,
                                                 double rotation_longest_valid_segment_length,
                                                 int min_steps)
{
  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  // Joint waypoints should have joint names
  assert(static_cast<long>(start.joint_names.size()) == start.size());

  // Initialize
  auto inv_kin = request.env->getManipulatorManager()->getInvKinematicSolver(mi.manipulator);
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  if (inv_kin->getJointNames() != fwd_kin->getJointNames())
    throw std::runtime_error("Forward and Inverse Kinematic objects joints are not ordered the same!");

  auto world_to_base = request.env_state->link_transforms.at(inv_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  assert(checkJointPositionFormat(fwd_kin->getJointNames(), start));
  Eigen::VectorXd j1 = start;

  // Calculate p2 in kinematics base frame without tcp for accurate comparison with p1
  Eigen::Isometry3d p2 = end * tcp.inverse();
  p2 = world_to_base.inverse() * p2;

  // Calculate FK for start
  Eigen::Isometry3d p1 = Eigen::Isometry3d::Identity();
  if (!fwd_kin->calcFwdKin(p1, j1))
    throw std::runtime_error("LVSInterpolateStateWaypoint: failed to find forward kinematics solution!");

  // Calculate steps based on cartesian information
  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);

  Eigen::VectorXd j2, j2_final;
  CompositeInstruction composite;

  bool j2_found = inv_kin->calcInvKin(j2, p2, j1);
  if (j2_found)
  {
    // Find closest solution to the start state
    double dist = std::numeric_limits<double>::max();
    const auto dof = inv_kin->numJoints();
    long num_solutions = j2.size() / dof;
    j2_final = j2.middleRows(0, dof);
    for (long i = 0; i < num_solutions; ++i)
    {
      /// @todo: May be nice to add contact checking to find best solution, but may not be neccessary because this is
      /// used to generate the seed.
      auto solution = j2.middleRows(i * dof, dof);
      double d = (solution - j1).norm();
      if (d < dist)
      {
        j2_final = solution;
        dist = d;
      }
    }
    double joint_dist = (j2_final - j1).norm();
    int state_steps = int(joint_dist / state_longest_valid_segment_length) + 1;
    steps = std::max(steps, state_steps);
  }

  // Check min steps requirement
  steps = std::max(steps, min_steps);

  // Get move type base on base instruction type
  MoveInstructionType move_type;
  if (base_instruction.isLinear())
    move_type = MoveInstructionType::LINEAR;
  else if (base_instruction.isFreespace())
    move_type = MoveInstructionType::FREESPACE;
  else
    throw std::runtime_error("LVS Interpolation: Unsupported Move Instruction Type!");

  if (j2_found)
  {
    // Linearly interpolate in joint space
    Eigen::MatrixXd states = interpolate(j1, j2_final, steps);

    // Convert to MoveInstructions
    for (long i = 1; i < states.cols(); ++i)
    {
      MoveInstruction move_instruction(StateWaypoint(fwd_kin->getJointNames(), states.col(i)), move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }
  else
  {
    // Convert to MoveInstructions
    StateWaypoint swp(fwd_kin->getJointNames(), j1);
    for (long i = 1; i < steps; ++i)
    {
      MoveInstruction move_instruction(swp, move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }

  return composite;
}

CompositeInstruction LVSInterpolateStateWaypoint(const CartesianWaypoint& start,
                                                 const JointWaypoint& end,
                                                 const PlanInstruction& base_instruction,
                                                 const PlannerRequest& request,
                                                 const ManipulatorInfo& manip_info,
                                                 double state_longest_valid_segment_length,
                                                 double translation_longest_valid_segment_length,
                                                 double rotation_longest_valid_segment_length,
                                                 int min_steps)
{
  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  // Joint waypoints should have joint names
  assert(static_cast<long>(end.joint_names.size()) == end.size());

  // Initialize
  auto inv_kin = request.env->getManipulatorManager()->getInvKinematicSolver(mi.manipulator);
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  if (inv_kin->getJointNames() != fwd_kin->getJointNames())
    throw std::runtime_error("Forward and Inverse Kinematic objects joints are not ordered the same!");

  auto world_to_base = request.env_state->link_transforms.at(inv_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  // Calculate p1 in kinematics base frame without tcp
  Eigen::Isometry3d p1 = start * tcp.inverse();
  p1 = world_to_base.inverse() * p1;

  // Calculate FK for end state
  Eigen::Isometry3d p2 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), end));
  Eigen::VectorXd j2 = end;
  if (!fwd_kin->calcFwdKin(p2, j2))
    throw std::runtime_error("LVSInterpolateStateWaypoint: failed to find forward kinematics solution!");

  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);

  Eigen::VectorXd j1, j1_final;
  CompositeInstruction composite;

  bool j1_found = inv_kin->calcInvKin(j1, p1, j2);
  if (j1_found)
  {
    // Find closest solution to the start state
    double dist = std::numeric_limits<double>::max();
    const auto dof = inv_kin->numJoints();
    long num_solutions = j1.size() / dof;
    j1_final = j1.middleRows(0, dof);
    for (long i = 0; i < num_solutions; ++i)
    {
      /// @todo: May be nice to add contact checking to find best solution, but may not be neccessary because this is
      /// used to generate the seed.
      auto solution = j1.middleRows(i * dof, dof);
      double d = (solution - j2).norm();
      if (d < dist)
      {
        j1_final = solution;
        dist = d;
      }
    }
    double joint_dist = (j1_final - j2).norm();
    int state_steps = int(joint_dist / state_longest_valid_segment_length) + 1;
    steps = std::max(steps, state_steps);
  }

  // Check min steps requirement
  steps = std::max(steps, min_steps);

  // Get move type base on base instruction type
  MoveInstructionType move_type;
  if (base_instruction.isLinear())
    move_type = MoveInstructionType::LINEAR;
  else if (base_instruction.isFreespace())
    move_type = MoveInstructionType::FREESPACE;
  else
    throw std::runtime_error("LVS Interpolation: Unsupported Move Instruction Type!");

  if (j1_found)
  {
    // Linearly interpolate in joint space
    Eigen::MatrixXd states = interpolate(j1_final, j2, steps);

    // Convert to MoveInstructions
    for (long i = 1; i < states.cols(); ++i)
    {
      MoveInstruction move_instruction(StateWaypoint(fwd_kin->getJointNames(), states.col(i)), move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }
  else
  {
    // Convert to MoveInstructions
    StateWaypoint swp(fwd_kin->getJointNames(), j2);
    for (long i = 1; i < steps; ++i)
    {
      MoveInstruction move_instruction(swp, move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }

  return composite;
}

CompositeInstruction LVSInterpolateStateWaypoint(const CartesianWaypoint& start,
                                                 const CartesianWaypoint& end,
                                                 const PlanInstruction& base_instruction,
                                                 const PlannerRequest& request,
                                                 const ManipulatorInfo& manip_info,
                                                 double state_longest_valid_segment_length,
                                                 double translation_longest_valid_segment_length,
                                                 double rotation_longest_valid_segment_length,
                                                 int min_steps)
{
  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  auto inv_kin = request.env->getManipulatorManager()->getInvKinematicSolver(mi.manipulator);
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  if (inv_kin->getJointNames() != fwd_kin->getJointNames())
    throw std::runtime_error("Forward and Inverse Kinematic objects joints are not ordered the same!");

  auto world_to_base = request.env_state->link_transforms.at(inv_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  // Get IK seed
  Eigen::VectorXd seed = request.env_state->getJointValues(inv_kin->getJointNames());

  // Calculate IK for start and end
  Eigen::Isometry3d p1 = start * tcp.inverse();
  p1 = world_to_base.inverse() * p1;
  Eigen::VectorXd j1, j1_final;
  bool found_j1 = inv_kin->calcInvKin(j1, p1, seed);

  Eigen::Isometry3d p2 = end * tcp.inverse();
  p2 = world_to_base.inverse() * p2;
  Eigen::VectorXd j2, j2_final;
  bool found_j2 = inv_kin->calcInvKin(j2, p2, seed);

  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);

  if (found_j1 && found_j2)
  {
    // Find closest solution to the end state
    double dist = std::numeric_limits<double>::max();
    const auto dof = inv_kin->numJoints();
    long j1_num_solutions = j1.size() / dof;
    long j2_num_solutions = j2.size() / dof;
    j1_final = j1.middleRows(0, dof);
    j2_final = j2.middleRows(0, dof);
    for (long i = 0; i < j1_num_solutions; ++i)
    {
      auto j1_solution = j1.middleRows(i * dof, dof);
      for (long j = 0; j < j2_num_solutions; ++j)
      {
        /// @todo: May be nice to add contact checking to find best solution, but may not be neccessary because this is
        /// used to generate the seed.
        auto j2_solution = j2.middleRows(j * dof, dof);
        double d = (j2_solution - j1_solution).norm();
        if (d < dist)
        {
          j1_final = j1_solution;
          j2_final = j2_solution;
          dist = d;
        }
      }
    }
    double joint_dist = (j1_final - j2_final).norm();
    int state_steps = int(joint_dist / state_longest_valid_segment_length) + 1;
    steps = std::max(steps, state_steps);
  }
  else if (found_j1)
  {
    double dist = std::numeric_limits<double>::max();
    const auto dof = inv_kin->numJoints();
    long j1_num_solutions = j1.size() / dof;
    j1_final = j1.middleRows(0, dof);
    for (long i = 0; i < j1_num_solutions; ++i)
    {
      auto j1_solution = j1.middleRows(i * dof, dof);
      double d = (seed - j1_solution).norm();
      if (d < dist)
      {
        j1_final = j1_solution;
        dist = d;
      }
    }
  }
  else if (found_j2)
  {
    double dist = std::numeric_limits<double>::max();
    const auto dof = inv_kin->numJoints();
    long j2_num_solutions = j2.size() / dof;
    j2_final = j2.middleRows(0, dof);
    for (long i = 0; i < j2_num_solutions; ++i)
    {
      auto j2_solution = j2.middleRows(i * dof, dof);
      double d = (seed - j2_solution).norm();
      if (d < dist)
      {
        j2_final = j2_solution;
        dist = d;
      }
    }
  }

  // Check min steps requirement
  steps = std::max(steps, min_steps);

  // Get move type base on base instruction type
  MoveInstructionType move_type;
  if (base_instruction.isLinear())
    move_type = MoveInstructionType::LINEAR;
  else if (base_instruction.isFreespace())
    move_type = MoveInstructionType::FREESPACE;
  else
    throw std::runtime_error("LVS Interpolation: Unsupported Move Instruction Type!");

  CompositeInstruction composite;

  if (found_j1 && found_j2)
  {
    // Linearly interpolate in joint space
    Eigen::MatrixXd states = interpolate(j1_final, j2_final, steps);

    // Convert to MoveInstructions
    for (long i = 1; i < states.cols(); ++i)
    {
      MoveInstruction move_instruction(StateWaypoint(inv_kin->getJointNames(), states.col(i)), move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }
  else if (found_j1)
  {
    // Convert to MoveInstructions
    StateWaypoint swp(inv_kin->getJointNames(), j1_final);
    for (long i = 1; i < steps; ++i)
    {
      MoveInstruction move_instruction(swp, move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }
  else if (found_j2)
  {
    // Convert to MoveInstructions
    StateWaypoint swp(inv_kin->getJointNames(), j2_final);
    for (long i = 1; i < steps; ++i)
    {
      MoveInstruction move_instruction(swp, move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }
  else
  {
    // Convert to MoveInstructions
    StateWaypoint swp(inv_kin->getJointNames(), seed);
    for (long i = 1; i < steps; ++i)
    {
      MoveInstruction move_instruction(swp, move_type);
      move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
      move_instruction.setDescription(base_instruction.getDescription());
      move_instruction.setProfile(base_instruction.getProfile());
      composite.push_back(move_instruction);
    }
  }

  return composite;
}

CompositeInstruction LVSInterpolateCartStateWaypoint(const JointWaypoint& start,
                                                     const JointWaypoint& end,
                                                     const PlanInstruction& base_instruction,
                                                     const PlannerRequest& request,
                                                     const ManipulatorInfo& manip_info,
                                                     double /*state_longest_valid_segment_length*/,
                                                     double translation_longest_valid_segment_length,
                                                     double rotation_longest_valid_segment_length,
                                                     int min_steps)
{
  /// @todo: Need to create a cartesian state waypoint and update the code below
  throw std::runtime_error("Not implemented, PR's are welcome!");

  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  // Initialize
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  auto world_to_base = request.env_state->link_transforms.at(fwd_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  CompositeInstruction composite;

  // Calculate FK for start and end
  Eigen::Isometry3d p1 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), start));
  if (!fwd_kin->calcFwdKin(p1, start))
    throw std::runtime_error("fixedSizeLinearInterpolation: failed to find forward kinematics solution!");
  p1 = world_to_base * p1 * tcp;

  Eigen::Isometry3d p2 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), end));
  if (!fwd_kin->calcFwdKin(p2, end))
    throw std::runtime_error("fixedSizeLinearInterpolation: failed to find forward kinematics solution!");
  p2 = world_to_base * p2 * tcp;

  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);
  steps = std::max(steps, min_steps);

  // Linear interpolation in cartesian space
  tesseract_common::VectorIsometry3d poses = interpolate(p1, p2, steps);

  // Convert to MoveInstructions
  for (const auto& pose : poses)
  {
    tesseract_planning::MoveInstruction move_instruction(CartesianWaypoint(pose), MoveInstructionType::LINEAR);
    move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
    move_instruction.setDescription(base_instruction.getDescription());
    move_instruction.setProfile(base_instruction.getProfile());
    composite.push_back(move_instruction);
  }

  return composite;
}

CompositeInstruction LVSInterpolateCartStateWaypoint(const JointWaypoint& start,
                                                     const CartesianWaypoint& end,
                                                     const PlanInstruction& base_instruction,
                                                     const PlannerRequest& request,
                                                     const ManipulatorInfo& manip_info,
                                                     double /*state_longest_valid_segment_length*/,
                                                     double translation_longest_valid_segment_length,
                                                     double rotation_longest_valid_segment_length,
                                                     int min_steps)
{
  /// @todo: Need to create a cartesian state waypoint and update the code below
  throw std::runtime_error("Not implemented, PR's are welcome!");

  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  // Initialize
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  auto world_to_base = request.env_state->link_transforms.at(fwd_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  CompositeInstruction composite;

  // Calculate FK for start and end
  Eigen::Isometry3d p1 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), start));
  if (!fwd_kin->calcFwdKin(p1, start))
    throw std::runtime_error("fixedSizeLinearInterpolation: failed to find forward kinematics solution!");
  p1 = world_to_base * p1 * tcp;

  Eigen::Isometry3d p2 = end;

  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);
  steps = std::max(steps, min_steps);

  // Linear interpolation in cartesian space
  tesseract_common::VectorIsometry3d poses = interpolate(p1, p2, steps);

  // Convert to MoveInstructions
  for (const auto& pose : poses)
  {
    tesseract_planning::MoveInstruction move_instruction(CartesianWaypoint(pose), MoveInstructionType::LINEAR);
    move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
    move_instruction.setDescription(base_instruction.getDescription());
    move_instruction.setProfile(base_instruction.getProfile());
    composite.push_back(move_instruction);
  }

  return composite;
}

CompositeInstruction LVSInterpolateCartStateWaypoint(const CartesianWaypoint& start,
                                                     const JointWaypoint& end,
                                                     const PlanInstruction& base_instruction,
                                                     const PlannerRequest& request,
                                                     const ManipulatorInfo& manip_info,
                                                     double /*state_longest_valid_segment_length*/,
                                                     double translation_longest_valid_segment_length,
                                                     double rotation_longest_valid_segment_length,
                                                     int min_steps)
{
  /// @todo: Need to create a cartesian state waypoint and update the code below
  throw std::runtime_error("Not implemented, PR's are welcome!");

  assert(!(manip_info.empty() && base_instruction.getManipulatorInfo().empty()));
  ManipulatorInfo mi = manip_info.getCombined(base_instruction.getManipulatorInfo());

  // Initialize
  auto fwd_kin = request.env->getManipulatorManager()->getFwdKinematicSolver(mi.manipulator);
  auto world_to_base = request.env_state->link_transforms.at(fwd_kin->getBaseLinkName());
  const Eigen::Isometry3d& tcp = request.env->findTCP(mi);

  CompositeInstruction composite;

  // Calculate FK for start and end
  Eigen::Isometry3d p1 = start;

  Eigen::Isometry3d p2 = Eigen::Isometry3d::Identity();
  assert(checkJointPositionFormat(fwd_kin->getJointNames(), end));
  if (!fwd_kin->calcFwdKin(p2, end))
    throw std::runtime_error("fixedSizeLinearInterpolation: failed to find forward kinematics solution!");
  p2 = world_to_base * p2 * tcp;

  double trans_dist = (p2.translation() - p1.translation()).norm();
  double rot_dist = Eigen::Quaterniond(p1.linear()).angularDistance(Eigen::Quaterniond(p2.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);
  steps = std::max(steps, min_steps);

  // Linear interpolation in cartesian space
  tesseract_common::VectorIsometry3d poses = interpolate(p1, p2, steps);

  // Convert to MoveInstructions
  for (const auto& pose : poses)
  {
    tesseract_planning::MoveInstruction move_instruction(CartesianWaypoint(pose), MoveInstructionType::LINEAR);
    move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
    move_instruction.setDescription(base_instruction.getDescription());
    move_instruction.setProfile(base_instruction.getProfile());
    composite.push_back(move_instruction);
  }

  return composite;
}

CompositeInstruction LVSInterpolateCartStateWaypoint(const CartesianWaypoint& start,
                                                     const CartesianWaypoint& end,
                                                     const PlanInstruction& base_instruction,
                                                     const PlannerRequest& /*request*/,
                                                     const ManipulatorInfo& /*manip_info*/,
                                                     double /*state_longest_valid_segment_length*/,
                                                     double translation_longest_valid_segment_length,
                                                     double rotation_longest_valid_segment_length,
                                                     int min_steps)
{
  /// @todo: Need to create a cartesian state waypoint and update the code below
  throw std::runtime_error("Not implemented, PR's are welcome!");

  CompositeInstruction composite;

  double trans_dist = (end.translation() - start.translation()).norm();
  double rot_dist = Eigen::Quaterniond(start.linear()).angularDistance(Eigen::Quaterniond(end.linear()));
  int trans_steps = int(trans_dist / translation_longest_valid_segment_length) + 1;
  int rot_steps = int(rot_dist / rotation_longest_valid_segment_length) + 1;
  int steps = std::max(trans_steps, rot_steps);
  steps = std::max(steps, min_steps);

  // Linear interpolation in cartesian space
  tesseract_common::VectorIsometry3d poses = interpolate(start, end, steps);

  // Convert to MoveInstructions
  for (const auto& pose : poses)
  {
    tesseract_planning::MoveInstruction move_instruction(CartesianWaypoint(pose), MoveInstructionType::LINEAR);
    move_instruction.setManipulatorInfo(base_instruction.getManipulatorInfo());
    move_instruction.setDescription(base_instruction.getDescription());
    move_instruction.setProfile(base_instruction.getProfile());
    composite.push_back(move_instruction);
  }
  return composite;
}

}  // namespace tesseract_planning
