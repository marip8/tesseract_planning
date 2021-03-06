/**
 * @file fixed_size_interpolation.h
 * @brief Provides interpolators where the number of steps are specified
 *
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
#ifndef TESSERACT_MOTION_PLANNERS_FIXED_SIZE_INTERPOLATION_H
#define TESSERACT_MOTION_PLANNERS_FIXED_SIZE_INTERPOLATION_H

#include <tesseract_common/macros.h>
TESSERACT_COMMON_IGNORE_WARNINGS_PUSH
#include <vector>
#include <memory>
TESSERACT_COMMON_IGNORE_WARNINGS_POP

#include <tesseract_command_language/command_language.h>
#include <tesseract_motion_planners/core/types.h>

namespace tesseract_planning
{
CompositeInstruction fixedSizeInterpolateStateWaypoint(const JointWaypoint& start,
                                                       const JointWaypoint& end,
                                                       const PlanInstruction& base_instruction,
                                                       const PlannerRequest& request,
                                                       const ManipulatorInfo& manip_info,
                                                       int steps);

CompositeInstruction fixedSizeInterpolateStateWaypoint(const JointWaypoint& start,
                                                       const CartesianWaypoint& end,
                                                       const PlanInstruction& base_instruction,
                                                       const PlannerRequest& request,
                                                       const ManipulatorInfo& manip_info,
                                                       int steps);

CompositeInstruction fixedSizeInterpolateStateWaypoint(const CartesianWaypoint& start,
                                                       const JointWaypoint& end,
                                                       const PlanInstruction& base_instruction,
                                                       const PlannerRequest& request,
                                                       const ManipulatorInfo& manip_info,
                                                       int steps);

CompositeInstruction fixedSizeInterpolateStateWaypoint(const CartesianWaypoint& start,
                                                       const CartesianWaypoint& end,
                                                       const PlanInstruction& base_instruction,
                                                       const PlannerRequest& request,
                                                       const ManipulatorInfo& manip_info,
                                                       int steps);

CompositeInstruction fixedSizeInterpolateCartStateWaypoint(const JointWaypoint& start,
                                                           const JointWaypoint&,
                                                           const PlanInstruction& base_instruction,
                                                           const PlannerRequest& request,
                                                           const ManipulatorInfo& manip_info,
                                                           int steps);

CompositeInstruction fixedSizeInterpolateCartStateWaypoint(const JointWaypoint& start,
                                                           const CartesianWaypoint& end,
                                                           const PlanInstruction& base_instruction,
                                                           const PlannerRequest& request,
                                                           const ManipulatorInfo& manip_info,
                                                           int steps);

CompositeInstruction fixedSizeInterpolateCartStateWaypoint(const CartesianWaypoint& start,
                                                           const JointWaypoint& end,
                                                           const PlanInstruction& base_instruction,
                                                           const PlannerRequest& request,
                                                           const ManipulatorInfo& manip_info,
                                                           int steps);

CompositeInstruction fixedSizeInterpolateCartStateWaypoint(const CartesianWaypoint& start,
                                                           const CartesianWaypoint& end,
                                                           const PlanInstruction& base_instruction,
                                                           const PlannerRequest& request,
                                                           const ManipulatorInfo& manip_info,
                                                           int steps);

}  // namespace tesseract_planning

#endif  // TESSERACT_MOTION_PLANNERS_FIXED_SIZE_INTERPOLATION_H
