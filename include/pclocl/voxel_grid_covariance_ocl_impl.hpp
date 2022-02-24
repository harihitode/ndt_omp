/*
 * Copyright 2021 Tier IV inc. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * This file includes works by PCL.
 *
 * ======== ORIGINAL LICENSE AND COPYRIGHTS BELOW ========
 *
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *
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
 *   * Neither the name of the copyright holder(s) nor the names of its
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
 *
 */

#ifndef PCL_VOXEL_GRID_COVARIANCE_IMPL_OCL_H_
#define PCL_VOXEL_GRID_COVARIANCE_IMPL_OCL_H_

#include <pcl/common/common.h>
#include <pcl/filters/boost.h>
#include "voxel_grid_covariance_ocl.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#define OCL_CREATE_BUFFER_CHECK(ret_mem, context, flags, size, host_ptr, errcode_ret) \
  ret_mem = clCreateBuffer(context, flags, size, host_ptr, &errcode_ret); \
if (errcode_ret != CL_SUCCESS) { \
  return -1; \
 }

#define OCL_ALLOC_SVM_CHECK(ret_mem, mem_type, context, flags, size, alignment) \
  ret_mem = static_cast<mem_type>(clSVMAlloc(context, flags, size, alignment)); \
  if (ret_mem == NULL) {                                                \
    return -1;                                                          \
  }

#define OCL_RELEASE_MEMORY_CHECK(mem)           \
  if (mem != NULL) {                            \
    clReleaseMemObject(mem);                    \
    mem = NULL;                                 \
  }

#define OCL_FREE_SVM_CHECK(context, svm_pointer)  \
  if (svm_pointer != NULL) {                      \
    clSVMFree(context, svm_pointer);              \
    svm_pointer = NULL;                           \
  }

//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> void
pclocl::VoxelGridCovariance<PointT>::applyFilter (PointCloud &output)
{
  voxel_centroids_leaf_indices_.clear ();

  // Has the input dataset been set already?
  if (!input_)
  {
    PCL_WARN ("[pcl::%s::applyFilter] No input dataset given!\n", getClassName ().c_str ());
    output.width = output.height = 0;
    output.points.clear ();
    return;
  }

  // Copy the header (and thus the frame_id) + allocate enough space for points
  output.height = 1;                          // downsampling breaks the organized structure
  output.is_dense = true;                     // we filter out invalid points
  output.points.clear ();

  Eigen::Vector4f min_p, max_p;
  // Get the minimum and maximum dimensions
  if (!filter_field_name_.empty ()) // If we don't want to process the entire cloud...
      pcl::getMinMax3D<PointT> (input_, filter_field_name_, static_cast<float> (filter_limit_min_), static_cast<float> (filter_limit_max_), min_p, max_p, filter_limit_negative_);
  else
	  pcl::getMinMax3D<PointT> (*input_, min_p, max_p);

  // Check that the leaf size is not too small, given the size of the data
  int64_t dx = static_cast<int64_t>((max_p[0] - min_p[0]) * inverse_leaf_size_[0])+1;
  int64_t dy = static_cast<int64_t>((max_p[1] - min_p[1]) * inverse_leaf_size_[1])+1;
  int64_t dz = static_cast<int64_t>((max_p[2] - min_p[2]) * inverse_leaf_size_[2])+1;

  if((dx*dy*dz) > std::numeric_limits<int32_t>::max())
  {
    PCL_WARN("[pcl::%s::applyFilter] Leaf size is too small for the input dataset. Integer indices would overflow.", getClassName().c_str());
    output.clear();
    return;
  }

  // Compute the minimum and maximum bounding box values
  min_b_[0] = static_cast<int> (floor (min_p[0] * inverse_leaf_size_[0]));
  max_b_[0] = static_cast<int> (floor (max_p[0] * inverse_leaf_size_[0]));
  min_b_[1] = static_cast<int> (floor (min_p[1] * inverse_leaf_size_[1]));
  max_b_[1] = static_cast<int> (floor (max_p[1] * inverse_leaf_size_[1]));
  min_b_[2] = static_cast<int> (floor (min_p[2] * inverse_leaf_size_[2]));
  max_b_[2] = static_cast<int> (floor (max_p[2] * inverse_leaf_size_[2]));

  // Compute the number of divisions needed along all axis
  div_b_ = max_b_ - min_b_ + Eigen::Vector4i::Ones ();
  div_b_[3] = 0;

  // Clear the leaves
  leaves_.clear ();
//  leaves_.reserve(8192);

  // Set up the division multiplier
  divb_mul_ = Eigen::Vector4i (1, div_b_[0], div_b_[0] * div_b_[1], 0);

  int centroid_size = 4;

  if (downsample_all_data_)
    centroid_size = boost::mpl::size<FieldList>::value;

  // ---[ RGB special case
  std::vector<pcl::PCLPointField> fields;
  int rgba_index = -1;
  rgba_index = pcl::getFieldIndex<PointT> ("rgb", fields);
  if (rgba_index == -1)
    rgba_index = pcl::getFieldIndex<PointT> ("rgba", fields);
  if (rgba_index >= 0)
  {
    rgba_index = fields[rgba_index].offset;
    centroid_size += 3;
  }

  // If we don't want to process the entire cloud, but rather filter points far away from the viewpoint first...
  if (!filter_field_name_.empty ())
  {
    // Get the distance field index
    std::vector<pcl::PCLPointField> fields;
    int distance_idx = pcl::getFieldIndex<PointT> (filter_field_name_, fields);
    if (distance_idx == -1)
      PCL_WARN ("[pcl::%s::applyFilter] Invalid filter field name. Index is %d.\n", getClassName ().c_str (), distance_idx);

    // First pass: go over all points and insert them into the right leaf
    for (size_t cp = 0; cp < input_->points.size (); ++cp)
    {
      if (!input_->is_dense)
        // Check if the point is invalid
        if (!std::isfinite (input_->points[cp].x) ||
            !std::isfinite (input_->points[cp].y) ||
            !std::isfinite (input_->points[cp].z))
          continue;

      // Get the distance value
      const uint8_t* pt_data = reinterpret_cast<const uint8_t*> (&input_->points[cp]);
      float distance_value = 0;
      memcpy (&distance_value, pt_data + fields[distance_idx].offset, sizeof (float));

      if (filter_limit_negative_)
      {
        // Use a threshold for cutting out points which inside the interval
        if ((distance_value < filter_limit_max_) && (distance_value > filter_limit_min_))
          continue;
      }
      else
      {
        // Use a threshold for cutting out points which are too close/far away
        if ((distance_value > filter_limit_max_) || (distance_value < filter_limit_min_))
          continue;
      }

      int ijk0 = static_cast<int> (floor (input_->points[cp].x * inverse_leaf_size_[0]) - static_cast<float> (min_b_[0]));
      int ijk1 = static_cast<int> (floor (input_->points[cp].y * inverse_leaf_size_[1]) - static_cast<float> (min_b_[1]));
      int ijk2 = static_cast<int> (floor (input_->points[cp].z * inverse_leaf_size_[2]) - static_cast<float> (min_b_[2]));

      // Compute the centroid leaf index
      int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

      Leaf& leaf = leaves_[idx];
      if (leaf.nr_points == 0)
      {
        leaf.centroid.resize (centroid_size);
        leaf.centroid.setZero ();
      }

      Eigen::Vector3d pt3d (input_->points[cp].x, input_->points[cp].y, input_->points[cp].z);
      // Accumulate point sum for centroid calculation
      leaf.mean_ += pt3d;
      // Accumulate x*xT for single pass covariance calculation
      leaf.cov_ += pt3d * pt3d.transpose ();

      // Do we need to process all the fields?
      if (!downsample_all_data_)
      {
        Eigen::Vector4f pt (input_->points[cp].x, input_->points[cp].y, input_->points[cp].z, 0);
        leaf.centroid.template head<4> () += pt;
      }
      else
      {
        // Copy all the fields
        Eigen::VectorXf centroid = Eigen::VectorXf::Zero (centroid_size);
        // ---[ RGB special case
        if (rgba_index >= 0)
        {
          // fill r/g/b data
          int rgb;
          memcpy (&rgb, reinterpret_cast<const char*> (&input_->points[cp]) + rgba_index, sizeof (int));
          centroid[centroid_size - 3] = static_cast<float> ((rgb >> 16) & 0x0000ff);
          centroid[centroid_size - 2] = static_cast<float> ((rgb >> 8) & 0x0000ff);
          centroid[centroid_size - 1] = static_cast<float> ((rgb) & 0x0000ff);
        }
        pcl::for_each_type<FieldList> (pcl::NdCopyPointEigenFunctor<PointT> (input_->points[cp], centroid));
        leaf.centroid += centroid;
      }
      ++leaf.nr_points;
    }
  }
  // No distance filtering, process all data
  else
  {
    // First pass: go over all points and insert them into the right leaf
    for (size_t cp = 0; cp < input_->points.size (); ++cp)
    {
      if (!input_->is_dense)
        // Check if the point is invalid
        if (!std::isfinite (input_->points[cp].x) ||
            !std::isfinite (input_->points[cp].y) ||
            !std::isfinite (input_->points[cp].z))
          continue;

      int ijk0 = static_cast<int> (floor (input_->points[cp].x * inverse_leaf_size_[0]) - static_cast<float> (min_b_[0]));
      int ijk1 = static_cast<int> (floor (input_->points[cp].y * inverse_leaf_size_[1]) - static_cast<float> (min_b_[1]));
      int ijk2 = static_cast<int> (floor (input_->points[cp].z * inverse_leaf_size_[2]) - static_cast<float> (min_b_[2]));

      // Compute the centroid leaf index
      int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

      //int idx = (((input_->points[cp].getArray4fmap () * inverse_leaf_size_).template cast<int> ()).matrix () - min_b_).dot (divb_mul_);
      Leaf& leaf = leaves_[idx];
      if (leaf.nr_points == 0)
      {
        leaf.centroid.resize (centroid_size);
        leaf.centroid.setZero ();
      }

      Eigen::Vector3d pt3d (input_->points[cp].x, input_->points[cp].y, input_->points[cp].z);
      // Accumulate point sum for centroid calculation
      leaf.mean_ += pt3d;
      // Accumulate x*xT for single pass covariance calculation
      leaf.cov_ += pt3d * pt3d.transpose ();

      // Do we need to process all the fields?
      if (!downsample_all_data_)
      {
        Eigen::Vector4f pt (input_->points[cp].x, input_->points[cp].y, input_->points[cp].z, 0);
        leaf.centroid.template head<4> () += pt;
      }
      else
      {
        // Copy all the fields
        Eigen::VectorXf centroid = Eigen::VectorXf::Zero (centroid_size);
        // ---[ RGB special case
        if (rgba_index >= 0)
        {
          // Fill r/g/b data, assuming that the order is BGRA
          int rgb;
          memcpy (&rgb, reinterpret_cast<const char*> (&input_->points[cp]) + rgba_index, sizeof (int));
          centroid[centroid_size - 3] = static_cast<float> ((rgb >> 16) & 0x0000ff);
          centroid[centroid_size - 2] = static_cast<float> ((rgb >> 8) & 0x0000ff);
          centroid[centroid_size - 1] = static_cast<float> ((rgb) & 0x0000ff);
        }
        pcl::for_each_type<FieldList> (pcl::NdCopyPointEigenFunctor<PointT> (input_->points[cp], centroid));
        leaf.centroid += centroid;
      }
      ++leaf.nr_points;
    }
  }

  // Second pass: go over all leaves and compute centroids and covariance matrices
  output.points.reserve (leaves_.size ());
  if (searchable_)
    voxel_centroids_leaf_indices_.reserve (leaves_.size ());
  int cp = 0;
  if (save_leaf_layout_)
    leaf_layout_.resize (div_b_[0] * div_b_[1] * div_b_[2], -1);

  // Eigen values and vectors calculated to prevent near singular matrices
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
  Eigen::Matrix3d eigen_val;
  Eigen::Vector3d pt_sum;

  // Eigen values less than a threshold of max eigen value are inflated to a set fraction of the max eigen value.
  double min_covar_eigvalue;

  for (auto it = leaves_.begin (); it != leaves_.end (); ++it)
  {

    // Normalize the centroid
    Leaf& leaf = it->second;

    // Normalize the centroid
    leaf.centroid /= static_cast<float> (leaf.nr_points);
    // Point sum used for single pass covariance calculation
    pt_sum = leaf.mean_;
    // Normalize mean
    leaf.mean_ /= leaf.nr_points;

    // If the voxel contains sufficient points, its covariance is calculated and is added to the voxel centroids and output clouds.
    // Points with less than the minimum points will have a can not be accurately approximated using a normal distribution.
    if (leaf.nr_points >= min_points_per_voxel_)
    {
      if (save_leaf_layout_)
        leaf_layout_[it->first] = cp++;

      output.push_back (PointT ());

      // Do we need to process all the fields?
      if (!downsample_all_data_)
      {
        output.points.back ().x = leaf.centroid[0];
        output.points.back ().y = leaf.centroid[1];
        output.points.back ().z = leaf.centroid[2];
      }
      else
      {
        pcl::for_each_type<FieldList> (pcl::NdCopyEigenPointFunctor<PointT> (leaf.centroid, output.back ()));
        // ---[ RGB special case
        if (rgba_index >= 0)
        {
          // pack r/g/b into rgb
          float r = leaf.centroid[centroid_size - 3], g = leaf.centroid[centroid_size - 2], b = leaf.centroid[centroid_size - 1];
          int rgb = (static_cast<int> (r)) << 16 | (static_cast<int> (g)) << 8 | (static_cast<int> (b));
          memcpy (reinterpret_cast<char*> (&output.points.back ()) + rgba_index, &rgb, sizeof (float));
        }
      }

      // Stores the voxel indices for fast access searching
      if (searchable_)
        voxel_centroids_leaf_indices_.push_back (static_cast<int> (it->first));

      // Single pass covariance calculation
      leaf.cov_ = (leaf.cov_ - 2 * (pt_sum * leaf.mean_.transpose ())) / leaf.nr_points + leaf.mean_ * leaf.mean_.transpose ();
      leaf.cov_ *= (leaf.nr_points - 1.0) / leaf.nr_points;

      //Normalize Eigen Val such that max no more than 100x min.
      eigensolver.compute (leaf.cov_);
      eigen_val = eigensolver.eigenvalues ().asDiagonal ();
      leaf.evecs_ = eigensolver.eigenvectors ();

      if (eigen_val (0, 0) < 0 || eigen_val (1, 1) < 0 || eigen_val (2, 2) <= 0)
      {
        leaf.nr_points = -1;
        continue;
      }

      // Avoids matrices near singularities (eq 6.11)[Magnusson 2009]

      min_covar_eigvalue = min_covar_eigvalue_mult_ * eigen_val (2, 2);
      if (eigen_val (0, 0) < min_covar_eigvalue)
      {
        eigen_val (0, 0) = min_covar_eigvalue;

        if (eigen_val (1, 1) < min_covar_eigvalue)
        {
          eigen_val (1, 1) = min_covar_eigvalue;
        }

        leaf.cov_ = leaf.evecs_ * eigen_val * leaf.evecs_.inverse ();
      }
      leaf.evals_ = eigen_val.diagonal ();

      leaf.icov_ = leaf.cov_.inverse ();
      if (leaf.icov_.maxCoeff () == std::numeric_limits<float>::infinity ( )
          || leaf.icov_.minCoeff () == -std::numeric_limits<float>::infinity ( ) )
      {
        leaf.nr_points = -1;
      }

    }
  }

  output.width = static_cast<uint32_t> (output.points.size ());
}

//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> int
pclocl::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint(const Eigen::MatrixXi& relative_coordinates, const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const
{
	neighbors.clear();

	// Find displacement coordinates
	Eigen::Vector4i ijk(static_cast<int> (floor(reference_point.x / leaf_size_[0])),
		static_cast<int> (floor(reference_point.y / leaf_size_[1])),
		static_cast<int> (floor(reference_point.z / leaf_size_[2])), 0);
	Eigen::Array4i diff2min = min_b_ - ijk;
	Eigen::Array4i diff2max = max_b_ - ijk;
	neighbors.reserve(relative_coordinates.cols());

	// Check each neighbor to see if it is occupied and contains sufficient points
	// Slower than radius search because needs to check 26 indices
	for (int ni = 0; ni < relative_coordinates.cols(); ni++)
	{
		Eigen::Vector4i displacement = (Eigen::Vector4i() << relative_coordinates.col(ni), 0).finished();
		// Checking if the specified cell is in the grid
		if ((diff2min <= displacement.array()).all() && (diff2max >= displacement.array()).all())
		{
			auto leaf_iter = leaves_.find(((ijk + displacement - min_b_).dot(divb_mul_)));
			if (leaf_iter != leaves_.end() && leaf_iter->second.nr_points >= min_points_per_voxel_)
			{
				LeafConstPtr leaf = &(leaf_iter->second);
				neighbors.push_back(leaf);
			}
		}
	}

	return (static_cast<int> (neighbors.size()));
}

//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> int
pclocl::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const
{
	neighbors.clear();

	// Find displacement coordinates
	Eigen::MatrixXi relative_coordinates = pcl::getAllNeighborCellIndices();
	return getNeighborhoodAtPoint(relative_coordinates, reference_point, neighbors);
}

//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> int
pclocl::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint7(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const
{
	neighbors.clear();

	Eigen::MatrixXi relative_coordinates(3, 7);
	relative_coordinates.setZero();
	relative_coordinates(0, 1) = 1;
	relative_coordinates(0, 2) = -1;
	relative_coordinates(1, 3) = 1;
	relative_coordinates(1, 4) = -1;
	relative_coordinates(2, 5) = 1;
	relative_coordinates(2, 6) = -1;

	return getNeighborhoodAtPoint(relative_coordinates, reference_point, neighbors);
}


//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> int
pclocl::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint1(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const
{
	neighbors.clear();
	return getNeighborhoodAtPoint(Eigen::MatrixXi::Zero(3,1), reference_point, neighbors);
}



//////////////////////////////////////////////////////////////////////////////////////////
template<typename PointT> void
pclocl::VoxelGridCovariance<PointT>::getDisplayCloud (pcl::PointCloud<pcl::PointXYZ>& cell_cloud)
{
  cell_cloud.clear ();

  int pnt_per_cell = 1000;
  boost::mt19937 rng;
  boost::normal_distribution<> nd (0.0, leaf_size_.head (3).norm ());
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor (rng, nd);

  Eigen::LLT<Eigen::Matrix3d> llt_of_cov;
  Eigen::Matrix3d cholesky_decomp;
  Eigen::Vector3d cell_mean;
  Eigen::Vector3d rand_point;
  Eigen::Vector3d dist_point;

  // Generate points for each occupied voxel with sufficient points.
  for (auto it = leaves_.begin (); it != leaves_.end (); ++it)
  {
    Leaf& leaf = it->second;

    if (leaf.nr_points >= min_points_per_voxel_)
    {
      cell_mean = leaf.mean_;
      llt_of_cov.compute (leaf.cov_);
      cholesky_decomp = llt_of_cov.matrixL ();

      // Random points generated by sampling the normal distribution given by voxel mean and covariance matrix
      for (int i = 0; i < pnt_per_cell; i++)
      {
        rand_point = Eigen::Vector3d (var_nor (), var_nor (), var_nor ());
        dist_point = cell_mean + cholesky_decomp * rand_point;
        cell_cloud.push_back (pcl::PointXYZ (static_cast<float> (dist_point (0)), static_cast<float> (dist_point (1)), static_cast<float> (dist_point (2))));
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
typename pclocl::VoxelGridCovariance<PointT>::kdtree_node_petit * pclocl::VoxelGridCovariance<PointT>::kdtree(
                                                                                                                kdtree_node_petit * root_node, kdtree_node_petit * parent_node, unsigned * alloc_pointer,
                                                                                                                point * pointlist, int left, int right, unsigned depth,
                                                                                                                int * node_indexes)
{
  point new_location;
  if (right < 0) {
    // not allocate
    return NULL;
  }
  // allocate
  kdtree_node_petit * new_node = &root_node[*alloc_pointer];
  *alloc_pointer = *alloc_pointer + 1; // count up

  // set parent
  new_node->parent = parent_node;

  // index list included in this node
  int new_left_index = left;
  int new_right_index = left + right + 1;
  new_node->left_right_index = ((new_left_index & 0x0000FFFF) << 16) | (new_right_index & 0x0000FFFF);
  // tree's depth
  // new_node->depth = depth;
  int median = right / 2;

  if (right == 0 || depth > MAX_DEPTH) {
    // new_node is leaf
    // new_node->location = pointlist[median];
    new_location = pointlist[median];
    // no child node
    new_node->child1 = NULL;
    new_node->child2 = NULL;

    qsort(pointlist, right + 1, sizeof(point), comparex);
    // new_node->leftmost = pointlist[0].x;
    // new_node->rightmost = pointlist[right].x;

    qsort(pointlist, right + 1, sizeof(point), comparey);
    // new_node->downmost = pointlist[0].y;
    // new_node->upmost = pointlist[right].y;

    qsort(pointlist, right + 1, sizeof(point), comparez);
    // new_node->zlowmost = pointlist[0].z;
    // new_node->zupmost = pointlist[right].z;

    for (int i = 0; i < right + 1; i++) {
      node_indexes[new_left_index + i] = pointlist[i].id;
    }

    return new_node;
  }

  // change the dividing dirction
  new_node->axis = depth % 3;

  // sorting by space (1 direction) for binary division
  if (new_node->axis == XAxis) {
    // ascending sorting by point.x
    qsort(pointlist, right + 1, sizeof(point), comparex);
    // new_node->leftmost = pointlist[0].x;
    // new_node->rightmost = pointlist[right].x;
  } else if (new_node->axis == YAxis) {
    // by point.y
    qsort(pointlist, right + 1, sizeof(point), comparey);
    // new_node->downmost = pointlist[0].y;
    // new_node->upmost = pointlist[right].y;
  } else if (new_node->axis == ZAxis) {
    // by point.z
    qsort(pointlist, right + 1, sizeof(point), comparez);
    // new_node->zlowmost = pointlist[0].z;
    // new_node->zupmost = pointlist[right].z;
  }


  node_indexes[new_left_index + median] = pointlist[median].id;
  // printf("node->left_index+median: %d\n", new_node->left_index+median);
  // new_node->location = pointlist[median];
  new_location = pointlist[median];
  new_node->child2 = kdtree(root_node, new_node, alloc_pointer,
                            pointlist + median + 1,
                            /*left  =*/ left + median + 1,
                            /*right =*/ right - (median + 1),
                            /*depth =*/ depth + 1, node_indexes);

  new_node->child1 = kdtree(root_node, new_node, alloc_pointer,
                            pointlist,
                            /*left  =*/ left,
                            /*right =*/ median - 1,
                            /*depth =*/ depth + 1, node_indexes);

  // calculate bounding box
  switch (new_node->axis) {
  case XAxis:
    new_node->axis_val = new_location.x;
    // if (new_node->child2 != NULL && new_node->child1 != NULL) {
    //   new_node->upmost = std::max(new_node->child2->upmost,
    //                               new_node->child1->upmost);
    //   new_node->downmost = std::min(new_node->child2->downmost,
    //                                 new_node->child1->downmost);
    //   new_node->zupmost = std::max(new_node->child2->zupmost,
    //                                new_node->child1->zupmost);
    //   new_node->zlowmost = std::min(new_node->child2->zlowmost,
    //                                 new_node->child1->zlowmost);
    // } else if (new_node->child2 != NULL) {
    //   new_node->upmost = new_node->child2->upmost;
    //   new_node->downmost = new_node->child2->downmost;
    //   new_node->zupmost = new_node->child2->zupmost;
    //   new_node->zlowmost = new_node->child2->zlowmost;
    // } else if (new_node->child1 != NULL) {
    //   new_node->upmost = new_node->child1->upmost;
    //   new_node->downmost = new_node->child1->downmost;
    //   new_node->zupmost = new_node->child1->zupmost;
    //   new_node->zlowmost = new_node->child1->zlowmost;
    // } else {
    //   new_node->upmost = new_node->location.y;
    //   new_node->downmost = new_node->location.y;
    //   new_node->zupmost = new_node->location.z;
    //   new_node->zlowmost =  new_node->location.z;
    // }
    break;
  case YAxis:
    new_node->axis_val = new_location.y;
    // if (new_node->child2 != NULL && new_node->child1 != NULL) {
    //   new_node->rightmost = std::max(new_node->child2->rightmost,
    //                                  new_node->child1->rightmost);
    //   new_node->leftmost = std::min(new_node->child2->leftmost,
    //                                 new_node->child1->leftmost);
    //   new_node->zupmost = std::max(new_node->child2->zupmost,
    //                                new_node->child1->zupmost);
    //   new_node->zlowmost = std::min(new_node->child2->zlowmost,
    //                                 new_node->child1->zlowmost);
    // } else if (new_node->child2 != NULL) {
    //   new_node->rightmost = new_node->child2->rightmost;
    //   new_node->leftmost = new_node->child2->leftmost;
    //   new_node->zupmost = new_node->child2->zupmost;
    //   new_node->zlowmost = new_node->child2->zlowmost;
    // } else if (new_node->child1 != NULL) {
    //   new_node->rightmost = new_node->child1->rightmost;
    //   new_node->leftmost = new_node->child1->leftmost;
    //   new_node->zupmost = new_node->child1->zupmost;
    //   new_node->zlowmost = new_node->child1->zlowmost;
    // } else {
    //   new_node->rightmost = new_node->location.x;
    //   new_node->leftmost = new_node->location.x;
    //   new_node->zupmost = new_node->location.z;
    //   new_node->zlowmost = new_node->location.z;
    // }
    break;
  case ZAxis:
  default:
    new_node->axis_val = new_location.z;
    // if (new_node->child2 != NULL && new_node->child1 != NULL) {
    //   new_node->rightmost = std::max(new_node->child2->rightmost,
    //                                  new_node->child1->rightmost);
    //   new_node->leftmost = std::min(new_node->child2->leftmost,
    //                                 new_node->child1->leftmost);
    //   new_node->upmost = std::max(new_node->child2->upmost,
    //                               new_node->child1->upmost);
    //   new_node->downmost = std::min(new_node->child2->downmost,
    //                                 new_node->child1->downmost);
    // } else if (new_node->child2 != NULL) {
    //   new_node->rightmost = new_node->child2->rightmost;
    //   new_node->leftmost = new_node->child2->leftmost;
    //   new_node->upmost = new_node->child2->upmost;
    //   new_node->downmost = new_node->child2->downmost;
    // } else if (new_node->child1 != NULL) {
    //   new_node->rightmost = new_node->child1->rightmost;
    //   new_node->leftmost = new_node->child1->leftmost;
    //   new_node->upmost = new_node->child1->upmost;
    //   new_node->downmost = new_node->child1->downmost;
    // } else {
    //   new_node->rightmost = new_node->location.x;
    //   new_node->leftmost = new_node->location.x;
    //   new_node->upmost = new_node->location.y;
    //   new_node->downmost = new_node->location.y;
    // }
    break;
  }
  return new_node;
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
int pclocl::VoxelGridCovariance<PointT>::createOCLMemoryObjects()
{
  size_t map_points_size, kdtree_size, map_indexes_size, map_mean_size, map_inverse_cov_size;
  cl_uint num_platforms, num_devices;
  cl_int ret;

  if (num_centroids_ > MAX_PCL_MAP_NUM) {
    return -1;
  }

  // memory object
  map_points_size = num_centroids_ * sizeof(float);
  map_indexes_size = num_centroids_ * sizeof(int);
  map_inverse_cov_size = num_centroids_ * 3 * 3 * sizeof(float);
  map_mean_size = num_centroids_ * 3 * sizeof(float);
  kdtree_size = num_centroids_ * sizeof(kdtree_node_petit);
  // should be have:
  // parent[pointer], child1[pointer], child2[pointer]: 12 bytes
  // left_index, right_index: (2^16 - 2^16): 4 bytes
  // axis_val: 4 bytes
  // axis: 4 bytes
  // total 24 bytes
  size_t kdtree_size_least = num_centroids_ * 24;
  printf("MAP_INFO: #VOXELS %lu\n", num_centroids_);
  printf("MAP_INFO: map_tree_size %lu [bytes]\n", kdtree_size_least);
  printf("MAP_INFO: map_points_size (x, y, z): %lu [bytes]\n", map_points_size * 3);
  printf("MAP_INFO: map_mean_size: %lu [bytes]\n", map_mean_size);
  printf("MAP_INFO: map_inverse_cov_size: %lu [bytes]\n", map_inverse_cov_size);
  printf("MAP_INFO: map_indexes_size: %lu [bytes]\n", map_indexes_size);
  printf("MAP_INFO: TOTAL memory: %lu [bytes]\n", kdtree_size_least + map_points_size * 3 + map_mean_size + map_inverse_cov_size + map_indexes_size);

  OCL_CREATE_BUFFER_CHECK(d_map_points_x_, context_, CL_MEM_READ_WRITE, map_points_size, NULL, ret);
  OCL_CREATE_BUFFER_CHECK(d_map_points_y_, context_, CL_MEM_READ_WRITE, map_points_size, NULL, ret);
  OCL_CREATE_BUFFER_CHECK(d_map_points_z_, context_, CL_MEM_READ_WRITE, map_points_size, NULL, ret);
  OCL_CREATE_BUFFER_CHECK(d_map_mean_, context_, CL_MEM_READ_WRITE, map_mean_size, NULL, ret);
  OCL_CREATE_BUFFER_CHECK(d_map_inverse_cov_, context_, CL_MEM_READ_WRITE, map_inverse_cov_size, NULL, ret);
  OCL_CREATE_BUFFER_CHECK(d_node_indexes_, context_, CL_MEM_READ_WRITE, map_indexes_size, NULL, ret);
  // allocation of shared virtual memory
  OCL_ALLOC_SVM_CHECK(kdtree_root_, kdtree_node_petit*, context_, CL_MEM_READ_WRITE, kdtree_size, 0);
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
void pclocl::VoxelGridCovariance<PointT>::releaseOCLMemoryObjects()
{
  OCL_RELEASE_MEMORY_CHECK(d_map_points_x_);
  OCL_RELEASE_MEMORY_CHECK(d_map_points_y_);
  OCL_RELEASE_MEMORY_CHECK(d_map_points_z_);
  OCL_RELEASE_MEMORY_CHECK(d_node_indexes_);
  OCL_RELEASE_MEMORY_CHECK(d_map_mean_);
  OCL_RELEASE_MEMORY_CHECK(d_map_inverse_cov_);
  OCL_FREE_SVM_CHECK(context_, kdtree_root_);
}

#define PCL_INSTANTIATE_VoxelGridCovariance(T) template class PCL_EXPORTS pcl::VoxelGridCovariance<T>;

#endif
