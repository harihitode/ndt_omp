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

#ifndef PCL_VOXEL_GRID_COVARIANCE_OCL_H_
#define PCL_VOXEL_GRID_COVARIANCE_OCL_H_

#include <pcl/pcl_macros.h>
#include <pcl/filters/boost.h>
#include <pcl/filters/voxel_grid.h>
#include <map>
#include <unordered_map>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_HPP_TARGET_OPENCL_VERSION 210
#define CL_HPP_ENABLE_EXCEPTIONS

#include <CL/cl2.hpp>

#define MAX_PCL_MAP_NUM (20000)
#define XAxis 0
#define YAxis 1
#define ZAxis 2
#define MAX_DEPTH 9

namespace pclocl
{
  /** \brief A searchable voxel structure containing the mean and covariance of the data.
    * \note For more information please see
    * <b>Magnusson, M. (2009). The Three-Dimensional Normal-Distributions Transform —
    * an Efficient Representation for Registration, Surface Analysis, and Loop Detection.
    * PhD thesis, Orebro University. Orebro Studies in Technology 36</b>
    * \author Brian Okorn (Space and Naval Warfare Systems Center Pacific)
    */
  template<typename PointT>
  class VoxelGridCovariance : public pcl::VoxelGrid<PointT>
  {
    protected:
      using pcl::VoxelGrid<PointT>::filter_name_;
      using pcl::VoxelGrid<PointT>::getClassName;
      using pcl::VoxelGrid<PointT>::input_;
      using pcl::VoxelGrid<PointT>::indices_;
      using pcl::VoxelGrid<PointT>::filter_limit_negative_;
      using pcl::VoxelGrid<PointT>::filter_limit_min_;
      using pcl::VoxelGrid<PointT>::filter_limit_max_;
      using pcl::VoxelGrid<PointT>::filter_field_name_;

      using pcl::VoxelGrid<PointT>::downsample_all_data_;
      using pcl::VoxelGrid<PointT>::leaf_layout_;
      using pcl::VoxelGrid<PointT>::save_leaf_layout_;
      using pcl::VoxelGrid<PointT>::leaf_size_;
      using pcl::VoxelGrid<PointT>::min_b_;
      using pcl::VoxelGrid<PointT>::max_b_;
      using pcl::VoxelGrid<PointT>::inverse_leaf_size_;
      using pcl::VoxelGrid<PointT>::div_b_;
      using pcl::VoxelGrid<PointT>::divb_mul_;

      typedef typename pcl::traits::fieldList<PointT>::type FieldList;
      typedef typename pcl::Filter<PointT>::PointCloud PointCloud;
      typedef typename PointCloud::Ptr PointCloudPtr;
      typedef typename PointCloud::ConstPtr PointCloudConstPtr;

    public:

#if PCL_VERSION >= PCL_VERSION_CALC(1, 10, 0)
      typedef pcl::shared_ptr< pcl::VoxelGrid<PointT> > Ptr;
      typedef pcl::shared_ptr< const pcl::VoxelGrid<PointT> > ConstPtr;
#else
      typedef boost::shared_ptr< pcl::VoxelGrid<PointT> > Ptr;
      typedef boost::shared_ptr< const pcl::VoxelGrid<PointT> > ConstPtr;
#endif

      /** \brief Simple structure to hold a centroid, covariance and the number of points in a leaf.
        * Inverse covariance, eigen vectors and eigen values are precomputed. */
      struct Leaf
      {
        /** \brief Constructor.
         * Sets \ref nr_points, \ref icov_, \ref mean_ and \ref evals_ to 0 and \ref cov_ and \ref evecs_ to the identity matrix
         */
        Leaf () :
          nr_points (0),
          mean_ (Eigen::Vector3d::Zero ()),
          centroid (),
          cov_ (Eigen::Matrix3d::Identity ()),
          icov_ (Eigen::Matrix3d::Zero ()),
          evecs_ (Eigen::Matrix3d::Identity ()),
          evals_ (Eigen::Vector3d::Zero ())
        {
        }

        /** \brief Get the voxel covariance.
          * \return covariance matrix
          */
        Eigen::Matrix3d
        getCov () const
        {
          return (cov_);
        }

        /** \brief Get the inverse of the voxel covariance.
          * \return inverse covariance matrix
          */
        Eigen::Matrix3d
        getInverseCov () const
        {
          return (icov_);
        }

        /** \brief Get the voxel centroid.
          * \return centroid
          */
        Eigen::Vector3d
        getMean () const
        {
          return (mean_);
        }

        /** \brief Get the eigen vectors of the voxel covariance.
          * \note Order corresponds with \ref getEvals
          * \return matrix whose columns contain eigen vectors
          */
        Eigen::Matrix3d
        getEvecs () const
        {
          return (evecs_);
        }

        /** \brief Get the eigen values of the voxel covariance.
          * \note Order corresponds with \ref getEvecs
          * \return vector of eigen values
          */
        Eigen::Vector3d
        getEvals () const
        {
          return (evals_);
        }

        /** \brief Get the number of points contained by this voxel.
          * \return number of points
          */
        int
        getPointCount () const
        {
          return (nr_points);
        }

        /** \brief Number of points contained by voxel */
        int nr_points;

        /** \brief 3D voxel centroid */
        Eigen::Vector3d mean_;

        /** \brief Nd voxel centroid
         * \note Differs from \ref mean_ when color data is used
         */
        Eigen::VectorXf centroid;

        /** \brief Voxel covariance matrix */
        Eigen::Matrix3d cov_;

        /** \brief Inverse of voxel covariance matrix */
        Eigen::Matrix3d icov_;

        /** \brief Eigen vectors of voxel covariance matrix */
        Eigen::Matrix3d evecs_;

        /** \brief Eigen values of voxel covariance matrix */
        Eigen::Vector3d evals_;

      };

      /** \brief Simple structure to hold 3-D position and id */
      typedef struct point {
        float x, y, z;
        int id;
      } point;

      /** \brief Simple structure of a kdtree node. */
      typedef struct tag_kdtree_node {
        int depth; // level of the node
        int left_index, right_index; // index of lower nodes than this node
        float rightmost, leftmost, upmost, downmost, zlowmost, zupmost; // minmax of the lower nodes than this node
        point location; // axis location of voxel division
        int axis; // axis of voxel division
        float axis_val; // value of axis delegate point
        struct tag_kdtree_node * parent;
        struct tag_kdtree_node * child1;
        struct tag_kdtree_node * child2;
      } kdtree_node;

      // 24 bytes in accelerator
      typedef struct tag_kdtree_node_petit {
        uint32_t left_right_index; // min&max indices of points included below this node
        uint32_t axis; // axis of division
        float axis_val; // value of division
        struct tag_kdtree_node_petit * parent;
        struct tag_kdtree_node_petit * child1;
        struct tag_kdtree_node_petit * child2;
      } kdtree_node_petit;

      /** \brief Pointer to VoxelGridCovariance leaf structure */
      typedef Leaf* LeafPtr;

      /** \brief Const pointer to VoxelGridCovariance leaf structure */
      typedef const Leaf* LeafConstPtr;

      typedef std::map<size_t, Leaf> Map;

      /** \brief KdTree generated using \ref voxel_centroids_ (used for searching). */
      kdtree_node_petit * kdtree_root_;

      /** \brief array of points of \ref voxel_centroids_. */
      point datapoints_[MAX_PCL_MAP_NUM];
      float map_points_x_[MAX_PCL_MAP_NUM];
      float map_points_y_[MAX_PCL_MAP_NUM];
      float map_points_z_[MAX_PCL_MAP_NUM];

      /** \brief array of Mean of voxel covariance matrix. */
      float map_mean_[MAX_PCL_MAP_NUM][3];

      /** \brief array of Inverse of voxel covariance matrix. */
      float map_inverse_cov_[MAX_PCL_MAP_NUM][3][3];

      /** \brief KdTree indexes. */
      int node_indexes_[MAX_PCL_MAP_NUM];

      /** \brief num of points of \ref voxel_centroids_. */
      size_t num_centroids_;

      /** \brief OpenCL context. */
      cl_context context_;

      /** \brief OpenCL memory objects. */
      cl_mem d_map_points_x_, d_map_points_y_, d_map_points_z_;
      cl_mem d_map_mean_, d_map_inverse_cov_;
      cl_mem d_node_indexes_;

    public:

      /** \brief Constructor.
       * Sets \ref leaf_size_ to 0 and \ref searchable_ to false.
       */
      VoxelGridCovariance () :
        searchable_ (true),
        min_points_per_voxel_ (6),
        min_covar_eigvalue_mult_ (0.01),
        leaves_ (),
        voxel_centroids_ (),
        voxel_centroids_leaf_indices_ (),
        kdtree_ ()
      {
        downsample_all_data_ = false;
        save_leaf_layout_ = false;
        leaf_size_.setZero ();
        min_b_.setZero ();
        max_b_.setZero ();
        filter_name_ = "VoxelGridCovariance";
      }

      /** \brief Finalize.
       * release OpenCL memory objects.
       */
      inline void finalize()
      {
        releaseOCLMemoryObjects();
      }

      inline void setContext(cl_context context)
      {
        if (context_ != NULL) {
          std::cerr << "error : OpenCL context is already set" << std::endl;
          return;
        }
        context_ = context;
      }

      /** \brief Set the minimum number of points required for a cell to be used (must be 3 or greater for covariance calculation).
        * \param[in] min_points_per_voxel the minimum number of points for required for a voxel to be used
        */
      inline void
      setMinPointPerVoxel (int min_points_per_voxel)
      {
        if(min_points_per_voxel > 2)
        {
          min_points_per_voxel_ = min_points_per_voxel;
        }
        else
        {
          PCL_WARN ("%s: Covariance calculation requires at least 3 points, setting Min Point per Voxel to 3 ", this->getClassName ().c_str ());
          min_points_per_voxel_ = 3;
        }
      }

      /** \brief Get the minimum number of points required for a cell to be used.
        * \return the minimum number of points for required for a voxel to be used
        */
      inline int
      getMinPointPerVoxel ()
      {
        return min_points_per_voxel_;
      }

      /** \brief Set the minimum allowable ratio between eigenvalues to prevent singular covariance matrices.
        * \param[in] min_covar_eigvalue_mult the minimum allowable ratio between eigenvalues
        */
      inline void
      setCovEigValueInflationRatio (double min_covar_eigvalue_mult)
      {
        min_covar_eigvalue_mult_ = min_covar_eigvalue_mult;
      }

      /** \brief Get the minimum allowable ratio between eigenvalues to prevent singular covariance matrices.
        * \return the minimum allowable ratio between eigenvalues
        */
      inline double
      getCovEigValueInflationRatio ()
      {
        return min_covar_eigvalue_mult_;
      }

      /** \brief Filter cloud and initializes voxel structure.
       * \param[out] output cloud containing centroids of voxels containing a sufficient number of points
       * \param[in] searchable flag if voxel structure is searchable, if true then kdtree is built
       */
      inline void
      filter (PointCloud &output, bool searchable = false)
      {
        searchable_ = searchable;
        applyFilter (output);

        voxel_centroids_ = PointCloudPtr (new PointCloud (output));

        constructKdTree();
      }

      /** \brief Initializes voxel structure.
       * \param[in] searchable flag if voxel structure is searchable, if true then kdtree is built
       */
      inline void
      filter (bool searchable = false)
      {
        searchable_ = searchable;
        voxel_centroids_ = PointCloudPtr (new PointCloud);
        applyFilter (*voxel_centroids_);

        constructKdTree();
      }

      /** \brief Construct KdTree from voxel structure.
       */
      inline void constructKdTree()
      {
        if (!searchable_ || num_centroids_ == 0 || num_centroids_ > MAX_PCL_MAP_NUM) {
          return;
        }
        // Initiates kdtree of the centroids of voxels containing a sufficient number of points
        for (int i = 0; i < num_centroids_; i++) {
          datapoints_[i].x = map_points_x_[i] = voxel_centroids_->points[i].x;
          datapoints_[i].y = map_points_y_[i] = voxel_centroids_->points[i].y;
          datapoints_[i].z = map_points_z_[i] = voxel_centroids_->points[i].z;
          datapoints_[i].id = i;
        }
        int ret = createOCLMemoryObjects();
        if (ret != 0) {
          std::cerr << "error : failed to create OpenCL memory objects" << std::endl;
          return;
        }
        // dump tree
        unsigned alloc_pointer = 0;
        if (kdtree(kdtree_root_, NULL, &alloc_pointer, datapoints_, 0, num_centroids_ - 1, 0, node_indexes_) == NULL) {
          std::cerr << "error : failed to construct kdtree" << std::endl;
          return;
        }
        {
          const uint32_t root_address = 0x20120000; // SPM1
          const uint32_t tree_node_size = 4 + 4 + 4 + (4 * 3);
          FILE * tree_file = fopen("sample_tree.txt", "w");
          fprintf(tree_file, "# SPM1 (note; big-endian); %08x\n", root_address);
          fprintf(tree_file, "# tree size; %08x\n", tree_node_size);
          fprintf(tree_file, "kdtree_root:\n");
          for (int i = 0; i < num_centroids_; i++) {
            fprintf(tree_file, "%08x\n", kdtree_root_[i].left_right_index);
            fprintf(tree_file, "%08x\n", kdtree_root_[i].axis);
            union f_and_i {
              uint32_t i;
              float f;
            } axis_val;
            axis_val.f = kdtree_root_[i].axis_val;
            fprintf(tree_file, "%08x\n", axis_val.i);
            if (kdtree_root_[i].parent != NULL) {
              fprintf(tree_file, "%08lx\n", root_address + (kdtree_root_[i].parent - kdtree_root_) * tree_node_size);
            } else {
              fprintf(tree_file, "%08x\n", 0);
            }
            if (kdtree_root_[i].child1 != NULL) {
              fprintf(tree_file, "%08lx\n", root_address + (kdtree_root_[i].child1 - kdtree_root_) * tree_node_size);
            } else {
              fprintf(tree_file, "%08x\n", 0);
            }
            if (kdtree_root_[i].child2 != NULL) {
              fprintf(tree_file, "%08lx\n", root_address + (kdtree_root_[i].child2 - kdtree_root_) * tree_node_size);
            } else {
              fprintf(tree_file, "%08x\n", 0);
            }
          }
          fclose(tree_file);
        }
        // dump infos:
        union f_and_i {
          uint32_t i;
          float f;
        } tmp;
        {
          FILE * fp = fopen("map_points_x.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "map_points_x:\n");
          for (int i = 0; i < num_centroids_; i++) {
            tmp.f = map_points_x_[i];
            fprintf(fp, "%08x\n", tmp.i);
          }
          fclose(fp);
        }
        {
          FILE * fp = fopen("map_points_y.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "map_points_y:\n");
          for (int i = 0; i < num_centroids_; i++) {
            tmp.f = map_points_y_[i];
            fprintf(fp, "%08x\n", tmp.i);
          }
          fclose(fp);
        }
        {
          FILE * fp = fopen("map_points_z.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "map_points_z:\n");
          for (int i = 0; i < num_centroids_; i++) {
            tmp.f = map_points_z_[i];
            fprintf(fp, "%08x\n", tmp.i);
          }
          fclose(fp);
        }
        {
          FILE * fp = fopen("map_mean.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "map_mean:\n");
          for (int i = 0; i < num_centroids_; i++) {
            for (int j = 0; j < 3; j++) {
              tmp.f = map_mean_[i][j];
              fprintf(fp, "%08x\n", tmp.i);
            }
          }
          fclose(fp);
        }
        {
          FILE * fp = fopen("map_inverse_cov.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "map_inverse_cov:\n");
          for (int i = 0; i < num_centroids_; i++) {
            for (int j = 0; j < 3; j++) {
              for (int k = 0; k < 3; k++) {
                tmp.f = map_inverse_cov_[i][j][k];
                fprintf(fp, "%08x\n", tmp.i);
              }
            }
          }
          fclose(fp);
        }
        {
          FILE * fp = fopen("node_indices.txt", "w");
          fprintf(fp, "# big endian\n");
          fprintf(fp, "node_indices:\n");
          for (int i = 0; i < num_centroids_; i++) {
            fprintf(fp, "%08x\n", node_indexes_[i]);
          }
          fclose(fp);
        }
      }

      /** \brief Get the voxel containing point p.
       * \param[in] index the index of the leaf structure node
       * \return const pointer to leaf structure
       */
      inline LeafConstPtr
      getLeaf (int index)
      {
        auto leaf_iter = leaves_.find (index);
        if (leaf_iter != leaves_.end ())
        {
          LeafConstPtr ret (&(leaf_iter->second));
          return ret;
        }
        else
          return NULL;
      }

      /** \brief Get the voxel containing point p.
       * \param[in] p the point to get the leaf structure at
       * \return const pointer to leaf structure
       */
      inline LeafConstPtr
      getLeaf (PointT &p)
      {
        // Generate index associated with p
        int ijk0 = static_cast<int> (floor (p.x * inverse_leaf_size_[0]) - min_b_[0]);
        int ijk1 = static_cast<int> (floor (p.y * inverse_leaf_size_[1]) - min_b_[1]);
        int ijk2 = static_cast<int> (floor (p.z * inverse_leaf_size_[2]) - min_b_[2]);

        // Compute the centroid leaf index
        int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

        // Find leaf associated with index
        auto leaf_iter = leaves_.find (idx);
        if (leaf_iter != leaves_.end ())
        {
          // If such a leaf exists return the pointer to the leaf structure
          LeafConstPtr ret (&(leaf_iter->second));
          return ret;
        }
        else
          return NULL;
      }

      /** \brief Get the voxel containing point p.
       * \param[in] p the point to get the leaf structure at
       * \return const pointer to leaf structure
       */
      inline LeafConstPtr
      getLeaf (Eigen::Vector3f &p)
      {
        // Generate index associated with p
        int ijk0 = static_cast<int> (floor (p[0] * inverse_leaf_size_[0]) - min_b_[0]);
        int ijk1 = static_cast<int> (floor (p[1] * inverse_leaf_size_[1]) - min_b_[1]);
        int ijk2 = static_cast<int> (floor (p[2] * inverse_leaf_size_[2]) - min_b_[2]);

        // Compute the centroid leaf index
        int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

        // Find leaf associated with index
        auto leaf_iter = leaves_.find (idx);
        if (leaf_iter != leaves_.end ())
        {
          // If such a leaf exists return the pointer to the leaf structure
          LeafConstPtr ret (&(leaf_iter->second));
          return ret;
        }
        else
          return NULL;

      }

      /** \brief Get the voxels surrounding point p, not including the voxel containing point p.
       * \note Only voxels containing a sufficient number of points are used (slower than radius search in practice).
       * \param[in] reference_point the point to get the leaf structure at
       * \param[out] neighbors
       * \return number of neighbors found
       */
	  int getNeighborhoodAtPoint(const Eigen::MatrixXi&, const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const ;
	  int getNeighborhoodAtPoint(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const ;
	  int getNeighborhoodAtPoint7(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const ;
	  int getNeighborhoodAtPoint1(const PointT& reference_point, std::vector<LeafConstPtr> &neighbors) const ;

      /** \brief Get the leaf structure map
       * \return a map containing all leaves
       */
      inline const Map&
      getLeaves ()
      {
        return leaves_;
      }

      /** \brief Get a pointcloud containing the voxel centroids
       * \note Only voxels containing a sufficient number of points are used.
       * \return a map containing all leaves
       */
      inline PointCloudPtr
      getCentroids ()
      {
        return voxel_centroids_;
      }


      /** \brief Get a cloud to visualize each voxels normal distribution.
       * \param[out] cell_cloud a cloud created by sampling the normal distributions of each voxel
       */
      void
      getDisplayCloud (pcl::PointCloud<pcl::PointXYZ>& cell_cloud);

      /** \brief Search for the k-nearest occupied voxels for the given query point.
       * \note Only voxels containing a sufficient number of points are used.
       * \param[in] point the given query point
       * \param[in] k the number of neighbors to search for
       * \param[out] k_leaves the resultant leaves of the neighboring points
       * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
       * \return number of neighbors found
       */
      int
      nearestKSearch (const PointT &point, int k,
                      std::vector<LeafConstPtr> &k_leaves, std::vector<float> &k_sqr_distances)
      {
        // k_leaves.clear ();

        // // Check if kdtree has been built
        // if (!searchable_)
        // {
        //   PCL_WARN ("%s: Not Searchable", this->getClassName ().c_str ());
        //   return 0;
        // }

        // // Find k-nearest neighbors in the occupied voxel centroid cloud
        // std::vector<int> k_indices;
        // k = kdtree_.nearestKSearch (point, k, k_indices, k_sqr_distances);

        // // Find leaves corresponding to neighbors
        // k_leaves.reserve (k);
        // for (std::vector<int>::iterator iter = k_indices.begin (); iter != k_indices.end (); iter++)
        // {
        //   k_leaves.push_back (&leaves_[voxel_centroids_leaf_indices_[*iter]]);
        // }
        // return k;
        return 0;
      }

      /** \brief Search for the k-nearest occupied voxels for the given query point.
       * \note Only voxels containing a sufficient number of points are used.
       * \param[in] cloud the given query point
       * \param[in] index the index
       * \param[in] k the number of neighbors to search for
       * \param[out] k_leaves the resultant leaves of the neighboring points
       * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
       * \return number of neighbors found
       */
      inline int
      nearestKSearch (const PointCloud &cloud, int index, int k,
                      std::vector<LeafConstPtr> &k_leaves, std::vector<float> &k_sqr_distances)
      {
        if (index >= static_cast<int> (cloud.points.size ()) || index < 0)
          return (0);
        return (nearestKSearch (cloud.points[index], k, k_leaves, k_sqr_distances));
      }


      /** \brief Search for all the nearest occupied voxels of the query point in a given radius.
       * \note Only voxels containing a sufficient number of points are used.
       * \param[in] point the given query point
       * \param[in] radius the radius of the sphere bounding all of p_q's neighbors
       * \param[out] k_leaves the resultant leaves of the neighboring points
       * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
       * \param[in] max_nn
       * \return number of neighbors found
       */
      int
      radiusSearch (const PointT &point, double radius, std::vector<LeafConstPtr> &k_leaves,
                    std::vector<float> &k_sqr_distances, unsigned int max_nn = 0) const
      {
      //   k_leaves.clear ();

      //   // Check if kdtree has been built
      //   if (!searchable_)
      //   {
      //     PCL_WARN ("%s: Not Searchable", this->getClassName ().c_str ());
      //     return 0;
      //   }

      //   // Find neighbors within radius in the occupied voxel centroid cloud
      //   std::vector<int> k_indices;
      //   int k = kdtree_.radiusSearch (point, radius, k_indices, k_sqr_distances, max_nn);

      //   // Find leaves corresponding to neighbors
      //   k_leaves.reserve (k);
      //   for (std::vector<int>::iterator iter = k_indices.begin (); iter != k_indices.end (); iter++)
      //   {
		  // auto leaf = leaves_.find(voxel_centroids_leaf_indices_[*iter]);
		  // if (leaf == leaves_.end()) {
			//   std::cerr << "error : could not find the leaf corresponding to the voxel" << std::endl;
			//   std::cin.ignore(1);
		  // }
      //     k_leaves.push_back (&(leaf->second));
      //   }
      //   return k;
        return 0;
      }

      /** \brief Search for all the nearest occupied voxels of the query point in a given radius.
       * \note Only voxels containing a sufficient number of points are used.
       * \param[in] cloud the given query point
       * \param[in] index a valid index in cloud representing a valid (i.e., finite) query point
       * \param[in] radius the radius of the sphere bounding all of p_q's neighbors
       * \param[out] k_leaves the resultant leaves of the neighboring points
       * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
       * \param[in] max_nn
       * \return number of neighbors found
       */
      inline int
      radiusSearch (const PointCloud &cloud, int index, double radius,
                    std::vector<LeafConstPtr> &k_leaves, std::vector<float> &k_sqr_distances,
                    unsigned int max_nn = 0) const
      {
        if (index >= static_cast<int> (cloud.points.size ()) || index < 0)
          return (0);
        return (radiusSearch (cloud.points[index], radius, k_leaves, k_sqr_distances, max_nn));
      }

    protected:

      /** \brief Filter cloud and initializes voxel structure.
       * \param[out] output cloud containing centroids of voxels containing a sufficient number of points
       */
      void applyFilter (PointCloud &output);

      /** \brief Construct KdTree recursively.
       * \param[out] root_node pointer of root node.
       * \param[in] parent_node pointer of parent node.
       * \param[in] alloc_pointer index of node to be allocated.
       * \param[in] pointlist points of nodes.
       * \param[in] left left index of nodes.
       * \param[in] right right index of nodes.
       * \param[in] depth depth of the node.
       * \param[out] node_indexes indexes of nodes.
       */
      kdtree_node_petit * kdtree (kdtree_node_petit * root_node,
                                  kdtree_node_petit * parent_node,
                                  unsigned * alloc_pointer, point * pointlist,
                                  int left, int right, unsigned depth, int * node_indexes);

      /** \brief Create OpenCL memory objects.
       * \return error code of creation.
       */
      int createOCLMemoryObjects();

      /** \brief release OpenCL memory objects.
       */
      void releaseOCLMemoryObjects();

      /** \brief Return whether a > 0 or not.
       * \param[in] a target value.
       * \return positive value if a > 0, negative value if not.
       */
      static int sign (float a) {
        return (a > 0) - (a < 0);
      }

      /** \brief Return whether a.x > b.x or not.
       * \param[in] a target point.
       * \param[in] b target point.
       * \return positive value if a.x > b.x, negative value if not.
       */
      static int comparex (const void * a, const void * b) {
        return sign(((point *)a)->x - ((point *)b)->x);
      }

      /** \brief Return whether a.y > b.y or not.
       * \param[in] a target point.
       * \param[in] b target point.
       * \return positive value if a.y > b.y, negative value if not.
       */
      static int comparey (const void * a, const void * b) {
        return sign(((point *)a)->y - ((point *)b)->y);
      }

      /** \brief Return whether a.z > b.z or not.
       * \param[in] a target point.
       * \param[in] b target point.
       * \return positive value if a.z > b.z, negative value if not.
       */
      static int comparez (const void * a, const void * b) {
        return sign(((point *)a)->z - ((point *)b)->z);
      }

      /** \brief Flag to determine if voxel structure is searchable. */
      bool searchable_;

      /** \brief Minimum points contained with in a voxel to allow it to be usable. */
      int min_points_per_voxel_;

      /** \brief Minimum allowable ratio between eigenvalues to prevent singular covariance matrices. */
      double min_covar_eigvalue_mult_;

      /** \brief Voxel structure containing all leaf nodes (includes voxels with less than a sufficient number of points). */
	  Map leaves_;

      /** \brief Point cloud containing centroids of voxels containing at least minimum number of points. */
      PointCloudPtr voxel_centroids_;

      /** \brief Indices of leaf structures associated with each point in \ref voxel_centroids_ (used for searching). */
      std::vector<int> voxel_centroids_leaf_indices_;

      /** \brief KdTree generated using \ref voxel_centroids_ (used for searching). */
      pcl::KdTreeFLANN<PointT> kdtree_;
  };
}

#endif  //#ifndef PCL_VOXEL_GRID_COVARIANCE_H_
