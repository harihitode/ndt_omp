/*
 *  Software License Agreement (BSD License)
 *
 *  OpenCL Kernel version.
 *  Copyright (c) 2021, Tier IV, Inc.
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
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
 * $Id$
 *
 */

/** \brief Simple structure to hold 3-D position and id */
typedef struct point {
  float x, y, z;
  int id;
} point;

/** \brief Simple structure of a kdtree node. */
typedef struct tag_kdtree_node {
  int depth;
  int left_index, right_index;
  float rightmost, leftmost, upmost, downmost, zlowmost, zupmost;
  point location;
  int axis;
  float axis_val;
  global struct tag_kdtree_node * parent;
  global struct tag_kdtree_node * child1;
  global struct tag_kdtree_node * child2;
} kdtree_node;

/** \brief Simple structure of a kdtree node. */
typedef struct tag_kdtree_node_petit {
  int left_right_index; // [upper 16bit is left, lower 16bit is right]
  int axis;
  float axis_val;
  // global struct tag_kdtree_node_petit * parent;
  // global struct tag_kdtree_node_petit * child1;
  // global struct tag_kdtree_node_petit * child2;
  int parent;
  int child1;
  int child2;
} kdtree_node_petit;

global kdtree_node_petit * get_node(global kdtree_node_petit * root, int addr) {
  const int base = 0x20120000;
  if (addr == 0) {
    return NULL;
  }
  int index = (addr - base) / 24;
  return &(root[index]);
}

/** \brief Compute point derivatives.
 * \note Equation 6.18-21 [Magnusson 2009].
 * \param[in] x point from the input cloud
 * \param[in] compute_hessian flag to calculate hessian, unnessissary for step calculation.
 * \param[in] j_ang precomputed angular gradient
 * \param[in] h_ang precomputed angular hessian
 */
void computePointDerivatives(const float x[4], float point_gradient_[4][6], float point_hessian_[24][6], global const float j_ang[8][4], global const float h_ang[16][4])
{
  float x4[4];
  x4[0] = x[0];
  x4[1] = x[1];
  x4[2] = x[2];
  x4[3] = 0.0f;

  // Calculate first derivative of Transformation Equation 6.17 w.r.t. transform vector p.
  // Derivative w.r.t. ith element of transform vector corresponds to column i, Equation 6.18 and 6.19 [Magnusson 2009]
  float x_j_ang[8];
  {
    x_j_ang[0] =
      j_ang[0][0] * x4[0] + j_ang[0][1] * x4[1] + j_ang[0][2] * x4[2] + j_ang[0][3] * x4[3];
    x_j_ang[1] =
      j_ang[1][0] * x4[0] + j_ang[1][1] * x4[1] + j_ang[1][2] * x4[2] + j_ang[1][3] * x4[3];
    x_j_ang[2] =
      j_ang[2][0] * x4[0] + j_ang[2][1] * x4[1] + j_ang[2][2] * x4[2] + j_ang[2][3] * x4[3];
    x_j_ang[3] =
      j_ang[3][0] * x4[0] + j_ang[3][1] * x4[1] + j_ang[3][2] * x4[2] + j_ang[3][3] * x4[3];
    x_j_ang[4] =
      j_ang[4][0] * x4[0] + j_ang[4][1] * x4[1] + j_ang[4][2] * x4[2] + j_ang[4][3] * x4[3];
    x_j_ang[5] =
      j_ang[5][0] * x4[0] + j_ang[5][1] * x4[1] + j_ang[5][2] * x4[2] + j_ang[5][3] * x4[3];
    x_j_ang[6] =
      j_ang[6][0] * x4[0] + j_ang[6][1] * x4[1] + j_ang[6][2] * x4[2] + j_ang[6][3] * x4[3];
    x_j_ang[7] =
      j_ang[7][0] * x4[0] + j_ang[7][1] * x4[1] + j_ang[7][2] * x4[2] + j_ang[7][3] * x4[3];
  }

  point_gradient_[1][3] = x_j_ang[0];
  point_gradient_[2][3] = x_j_ang[1];
  point_gradient_[0][4] = x_j_ang[2];
  point_gradient_[1][4] = x_j_ang[3];
  point_gradient_[2][4] = x_j_ang[4];
  point_gradient_[0][5] = x_j_ang[5];
  point_gradient_[1][5] = x_j_ang[6];
  point_gradient_[2][5] = x_j_ang[7];

  float x_h_ang[16];
  {
    x_h_ang[0] =
      h_ang[0][0] * x4[0] + h_ang[0][1] * x4[1] + h_ang[0][2] * x4[2] + h_ang[0][3] * x4[3];
    x_h_ang[1] =
      h_ang[1][0] * x4[0] + h_ang[1][1] * x4[1] + h_ang[1][2] * x4[2] + h_ang[1][3] * x4[3];
    x_h_ang[2] =
      h_ang[2][0] * x4[0] + h_ang[2][1] * x4[1] + h_ang[2][2] * x4[2] + h_ang[2][3] * x4[3];
    x_h_ang[3] =
      h_ang[3][0] * x4[0] + h_ang[3][1] * x4[1] + h_ang[3][2] * x4[2] + h_ang[3][3] * x4[3];
    x_h_ang[4] =
      h_ang[4][0] * x4[0] + h_ang[4][1] * x4[1] + h_ang[4][2] * x4[2] + h_ang[4][3] * x4[3];
    x_h_ang[5] =
      h_ang[5][0] * x4[0] + h_ang[5][1] * x4[1] + h_ang[5][2] * x4[2] + h_ang[5][3] * x4[3];
    x_h_ang[6] =
      h_ang[6][0] * x4[0] + h_ang[6][1] * x4[1] + h_ang[6][2] * x4[2] + h_ang[6][3] * x4[3];
    x_h_ang[7] =
      h_ang[7][0] * x4[0] + h_ang[7][1] * x4[1] + h_ang[7][2] * x4[2] + h_ang[7][3] * x4[3];
    x_h_ang[8] =
      h_ang[8][0] * x4[0] + h_ang[8][1] * x4[1] + h_ang[8][2] * x4[2] + h_ang[8][3] * x4[3];
    x_h_ang[9] =
      h_ang[9][0] * x4[0] + h_ang[9][1] * x4[1] + h_ang[9][2] * x4[2] + h_ang[9][3] * x4[3];
    x_h_ang[10] =
      h_ang[10][0] * x4[0] + h_ang[10][1] * x4[1] + h_ang[10][2] * x4[2] + h_ang[10][3] * x4[3];
    x_h_ang[11] =
      h_ang[11][0] * x4[0] + h_ang[11][1] * x4[1] + h_ang[11][2] * x4[2] + h_ang[11][3] * x4[3];
    x_h_ang[12] =
      h_ang[12][0] * x4[0] + h_ang[12][1] * x4[1] + h_ang[12][2] * x4[2] + h_ang[12][3] * x4[3];
    x_h_ang[13] =
      h_ang[13][0] * x4[0] + h_ang[13][1] * x4[1] + h_ang[13][2] * x4[2] + h_ang[13][3] * x4[3];
    x_h_ang[14] =
      h_ang[14][0] * x4[0] + h_ang[14][1] * x4[1] + h_ang[14][2] * x4[2] + h_ang[14][3] * x4[3];
    x_h_ang[15] =
      h_ang[15][0] * x4[0] + h_ang[15][1] * x4[1] + h_ang[15][2] * x4[2] + h_ang[15][3] * x4[3];
  }

  // Vectors from Equation 6.21 [Magnusson 2009]
  float a[4];
  float b[4];
  float c[4];
  float d[4];
  float e[4];
  float f[4];
  a[0] = 0.0f;
  a[1] = x_h_ang[0];
  a[2] = x_h_ang[1];
  a[3] = 0.0f;
  b[0] = 0.0f;
  b[1] = x_h_ang[2];
  b[2] = x_h_ang[3];
  b[3] = 0.0f;
  c[0] = 0.0f;
  c[1] = x_h_ang[4];
  c[2] = x_h_ang[5];
  c[3] = 0.0f;
  d[0] = x_h_ang[6];
  d[1] = x_h_ang[7];
  d[2] = x_h_ang[8];
  d[3] = 0.0f;
  e[0] = x_h_ang[9];
  e[1] = x_h_ang[10];
  e[2] = x_h_ang[11];
  e[3] = 0.0f;
  f[0] = x_h_ang[12];
  f[1] = x_h_ang[13];
  f[2] = x_h_ang[14];
  f[3] = x_h_ang[15];

  // Calculate second derivative of Transformation Equation 6.17 w.r.t. transform vector p.
  // Derivative w.r.t. ith and jth elements of transform vector corresponds to the 3x1 block matrix starting at
  // (3i,j), Equation 6.20 and 6.21 [Magnusson 2009]
  point_hessian_[12][3] = a[0];
  point_hessian_[13][3] = a[1];
  point_hessian_[14][3] = a[2];
  point_hessian_[15][3] = a[3];

  point_hessian_[16][3] = b[0];
  point_hessian_[17][3] = b[1];
  point_hessian_[18][3] = b[2];
  point_hessian_[19][3] = b[3];

  point_hessian_[20][3] = c[0];
  point_hessian_[21][3] = c[1];
  point_hessian_[22][3] = c[2];
  point_hessian_[23][3] = c[3];

  point_hessian_[12][4] = b[0];
  point_hessian_[13][4] = b[1];
  point_hessian_[14][4] = b[2];
  point_hessian_[15][4] = b[3];

  point_hessian_[16][4] = d[0];
  point_hessian_[17][4] = d[1];
  point_hessian_[18][4] = d[2];
  point_hessian_[19][4] = d[3];

  point_hessian_[20][4] = e[0];
  point_hessian_[21][4] = e[1];
  point_hessian_[22][4] = e[2];
  point_hessian_[23][4] = e[3];

  point_hessian_[12][5] = c[0];
  point_hessian_[13][5] = c[1];
  point_hessian_[14][5] = c[2];
  point_hessian_[15][5] = c[3];

  point_hessian_[16][5] = e[0];
  point_hessian_[17][5] = e[1];
  point_hessian_[18][5] = e[2];
  point_hessian_[19][5] = e[3];

  point_hessian_[20][5] = f[0];
  point_hessian_[21][5] = f[1];
  point_hessian_[22][5] = f[2];
  point_hessian_[23][5] = f[3];
}

float updateDerivatives(
                        float score_gradient[6], float hessian[6][6], const float point_gradient4[4][6],
                        const float point_hessian_[24][6], const float x_trans[3], global const float c_inv[3][3], const float gauss_d1_,
                        const float gauss_d2_, const float gauss_d3_)
{
  float x_trans4[4];
  x_trans4[0] = x_trans[0];
  x_trans4[1] = x_trans[1];
  x_trans4[2] = x_trans[2];
  x_trans4[3] = 0.0f;

  float c_inv4[4][4];
  c_inv4[0][0] = c_inv[0][0];
  c_inv4[0][1] = c_inv[0][1];
  c_inv4[0][2] = c_inv[0][2];
  c_inv4[0][3] = 0.0f;
  c_inv4[1][0] = c_inv[1][0];
  c_inv4[1][1] = c_inv[1][1];
  c_inv4[1][2] = c_inv[1][2];
  c_inv4[1][3] = 0.0f;
  c_inv4[2][0] = c_inv[2][0];
  c_inv4[2][1] = c_inv[2][1];
  c_inv4[2][2] = c_inv[2][2];
  c_inv4[2][3] = 0.0f;
  c_inv4[3][0] = 0.0f;
  c_inv4[3][1] = 0.0f;
  c_inv4[3][2] = 0.0f;
  c_inv4[3][3] = 0.0f;

  float x_trans4_x_c_inv4[4];
  x_trans4_x_c_inv4[0] = x_trans4[0] * c_inv4[0][0] + x_trans4[1] * c_inv4[0][1] +
    x_trans4[2] * c_inv4[0][2] + x_trans4[3] * c_inv4[0][3];

  x_trans4_x_c_inv4[1] = x_trans4[0] * c_inv4[1][0] + x_trans4[1] * c_inv4[1][1] +
    x_trans4[2] * c_inv4[1][2] + x_trans4[3] * c_inv4[1][3];

  x_trans4_x_c_inv4[2] = x_trans4[0] * c_inv4[2][0] + x_trans4[1] * c_inv4[2][1] +
    x_trans4[2] * c_inv4[2][2] + x_trans4[3] * c_inv4[2][3];

  x_trans4_x_c_inv4[3] = x_trans4[0] * c_inv4[3][0] + x_trans4[1] * c_inv4[3][1] +
    x_trans4[2] * c_inv4[3][2] + x_trans4[3] * c_inv4[3][3];

  float x_trans4_dot_x_trans4_x_c_inv4 =
    x_trans4[0] * x_trans4_x_c_inv4[0] + x_trans4[1] * x_trans4_x_c_inv4[1] +
    x_trans4[2] * x_trans4_x_c_inv4[2] + x_trans4[3] * x_trans4_x_c_inv4[3];
  // e^(-d_2/2 * (x_k - mu_k)^T Sigma_k^-1 (x_k - mu_k)) Equation 6.9 [Magnusson 2009]
  // float e_x_cov_x = exp(-gauss_d2_ * x_trans4_dot_x_trans4_x_c_inv4 * 0.5f);
  float e_x_cov_x = exp(-gauss_d2_ * x_trans4_dot_x_trans4_x_c_inv4 * 0.5f);
  // Calculate probability of transtormed points existance, Equation 6.9 [Magnusson 2009]
  float score_inc = -gauss_d1_ * e_x_cov_x;

  e_x_cov_x = gauss_d2_ * e_x_cov_x;

  // Error checking for invalid values.
  if (e_x_cov_x > 1.0 || e_x_cov_x < 0.0 || e_x_cov_x != e_x_cov_x) return 0;

  // Reusable portion of Equation 6.12 and 6.13 [Magnusson 2009]
  e_x_cov_x *= gauss_d1_;

  float c_inv4_x_point_gradient4[4][6];
  c_inv4_x_point_gradient4[0][0] =
    c_inv4[0][0] * point_gradient4[0][0] + c_inv4[0][1] * point_gradient4[1][0] +
    c_inv4[0][2] * point_gradient4[2][0] + c_inv4[0][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[0][1] =
    c_inv4[0][0] * point_gradient4[0][1] + c_inv4[0][1] * point_gradient4[1][1] +
    c_inv4[0][2] * point_gradient4[2][1] + c_inv4[0][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[0][2] =
    c_inv4[0][0] * point_gradient4[0][2] + c_inv4[0][1] * point_gradient4[1][2] +
    c_inv4[0][2] * point_gradient4[2][2] + c_inv4[0][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[0][3] =
    c_inv4[0][0] * point_gradient4[0][3] + c_inv4[0][1] * point_gradient4[1][3] +
    c_inv4[0][2] * point_gradient4[2][3] + c_inv4[0][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[0][4] =
    c_inv4[0][0] * point_gradient4[0][4] + c_inv4[0][1] * point_gradient4[1][4] +
    c_inv4[0][2] * point_gradient4[2][4] + c_inv4[0][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[0][5] =
    c_inv4[0][0] * point_gradient4[0][5] + c_inv4[0][1] * point_gradient4[1][5] +
    c_inv4[0][2] * point_gradient4[2][5] + c_inv4[0][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[1][0] =
    c_inv4[1][0] * point_gradient4[0][0] + c_inv4[1][1] * point_gradient4[1][0] +
    c_inv4[1][2] * point_gradient4[2][0] + c_inv4[1][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[1][1] =
    c_inv4[1][0] * point_gradient4[0][1] + c_inv4[1][1] * point_gradient4[1][1] +
    c_inv4[1][2] * point_gradient4[2][1] + c_inv4[1][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[1][2] =
    c_inv4[1][0] * point_gradient4[0][2] + c_inv4[1][1] * point_gradient4[1][2] +
    c_inv4[1][2] * point_gradient4[2][2] + c_inv4[1][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[1][3] =
    c_inv4[1][0] * point_gradient4[0][3] + c_inv4[1][1] * point_gradient4[1][3] +
    c_inv4[1][2] * point_gradient4[2][3] + c_inv4[1][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[1][4] =
    c_inv4[1][0] * point_gradient4[0][4] + c_inv4[1][1] * point_gradient4[1][4] +
    c_inv4[1][2] * point_gradient4[2][4] + c_inv4[1][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[1][5] =
    c_inv4[1][0] * point_gradient4[0][5] + c_inv4[1][1] * point_gradient4[1][5] +
    c_inv4[1][2] * point_gradient4[2][5] + c_inv4[1][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[2][0] =
    c_inv4[2][0] * point_gradient4[0][0] + c_inv4[2][1] * point_gradient4[1][0] +
    c_inv4[2][2] * point_gradient4[2][0] + c_inv4[2][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[2][1] =
    c_inv4[2][0] * point_gradient4[0][1] + c_inv4[2][1] * point_gradient4[1][1] +
    c_inv4[2][2] * point_gradient4[2][1] + c_inv4[2][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[2][2] =
    c_inv4[2][0] * point_gradient4[0][2] + c_inv4[2][1] * point_gradient4[1][2] +
    c_inv4[2][2] * point_gradient4[2][2] + c_inv4[2][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[2][3] =
    c_inv4[2][0] * point_gradient4[0][3] + c_inv4[2][1] * point_gradient4[1][3] +
    c_inv4[2][2] * point_gradient4[2][3] + c_inv4[2][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[2][4] =
    c_inv4[2][0] * point_gradient4[0][4] + c_inv4[2][1] * point_gradient4[1][4] +
    c_inv4[2][2] * point_gradient4[2][4] + c_inv4[2][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[2][5] =
    c_inv4[2][0] * point_gradient4[0][5] + c_inv4[2][1] * point_gradient4[1][5] +
    c_inv4[2][2] * point_gradient4[2][5] + c_inv4[2][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[3][0] =
    c_inv4[3][0] * point_gradient4[0][0] + c_inv4[3][1] * point_gradient4[1][0] +
    c_inv4[3][2] * point_gradient4[2][0] + c_inv4[3][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[3][1] =
    c_inv4[3][0] * point_gradient4[0][1] + c_inv4[3][1] * point_gradient4[1][1] +
    c_inv4[3][2] * point_gradient4[2][1] + c_inv4[3][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[3][2] =
    c_inv4[3][0] * point_gradient4[0][2] + c_inv4[3][1] * point_gradient4[1][2] +
    c_inv4[3][2] * point_gradient4[2][2] + c_inv4[3][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[3][3] =
    c_inv4[3][0] * point_gradient4[0][3] + c_inv4[3][1] * point_gradient4[1][3] +
    c_inv4[3][2] * point_gradient4[2][3] + c_inv4[3][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[3][4] =
    c_inv4[3][0] * point_gradient4[0][4] + c_inv4[3][1] * point_gradient4[1][4] +
    c_inv4[3][2] * point_gradient4[2][4] + c_inv4[3][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[3][5] =
    c_inv4[3][0] * point_gradient4[0][5] + c_inv4[3][1] * point_gradient4[1][5] +
    c_inv4[3][2] * point_gradient4[2][5] + c_inv4[3][3] * point_gradient4[3][5];

  ////

  float x_trans4_dot_c_inv4_x_point_gradient4[6];
  x_trans4_dot_c_inv4_x_point_gradient4[0] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][0] + x_trans4[1] * c_inv4_x_point_gradient4[1][0] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][0] + x_trans4[3] * c_inv4_x_point_gradient4[3][0];

  x_trans4_dot_c_inv4_x_point_gradient4[1] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][1] + x_trans4[1] * c_inv4_x_point_gradient4[1][1] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][1] + x_trans4[3] * c_inv4_x_point_gradient4[3][1];

  x_trans4_dot_c_inv4_x_point_gradient4[2] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][2] + x_trans4[1] * c_inv4_x_point_gradient4[1][2] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][2] + x_trans4[3] * c_inv4_x_point_gradient4[3][2];

  x_trans4_dot_c_inv4_x_point_gradient4[3] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][3] + x_trans4[1] * c_inv4_x_point_gradient4[1][3] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][3] + x_trans4[3] * c_inv4_x_point_gradient4[3][3];

  x_trans4_dot_c_inv4_x_point_gradient4[4] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][4] + x_trans4[1] * c_inv4_x_point_gradient4[1][4] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][4] + x_trans4[3] * c_inv4_x_point_gradient4[3][4];

  x_trans4_dot_c_inv4_x_point_gradient4[5] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][5] + x_trans4[1] * c_inv4_x_point_gradient4[1][5] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][5] + x_trans4[3] * c_inv4_x_point_gradient4[3][5];

  score_gradient[0] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[0];
  score_gradient[1] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[1];
  score_gradient[2] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[2];
  score_gradient[3] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[3];
  score_gradient[4] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[4];
  score_gradient[5] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[5];

  ///
  float point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[6][6];
  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][0] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][1] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][2] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][3] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][4] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][5] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][0] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][1] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][2] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][3] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][4] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][5] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][0] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][1] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][2] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][3] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][4] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][5] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][0] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][1] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][2] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][3] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][4] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][5] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][0] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][1] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][2] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][3] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][4] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][5] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][0] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][1] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][2] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][3] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][4] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][5] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][5];

  float x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[6];

  for (int i = 0; i < 6; i++) {
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    for (int j = 0; j < 6; j++) {
      hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                    point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    }
  }

  return score_inc;
}

/** \brief Return the distance between two points.
 * \param[in] point1 input point 1
 * \param[in] point2 input point 2
 * \return the distance between two points
 */
float calculateDistance(const float point1[4], const float point2[4])
{
  float dist_square[4];
  for (int d = 0; d < 4; d++) {
    dist_square[d] = (point1[d] - point2[d]) * (point1[d] - point2[d]);
  }
  return sqrt(dist_square[0] + dist_square[1] + dist_square[2] + dist_square[3]);
}

/** \brief Return the indexes of n neighbor voxels of the reference point.
 * \param[in] lidar_point_x x value of input point (n elements)
 * \param[in] lidar_point_y y value of input point (n elements)
 * \param[in] lidar_point_z z value of input point (n elements)
 * \param[in] map_points_x 1-D array of x value of map points (? elements)
 * \param[in] map_points_y 1-D array of y value of map points (? elements)
 * \param[in] map_points_z 1-D array of z value of map points (? elements)
 * \param[in] node_indexes 1-D array of indexes of map points (? elements)
 * \param[in] root_node kdtree root node
 * \param[in] n number of input points
 * \param[in] limit limit of neighbors
 * \param[in] radius search radius
 * \param[out] neighbor_candidate_indexes 1-D array of indexes of neighbor candidates (? elements)
 * \param[out] neighbor_candidate_dists 1-D array of distances of neighbor candidates (? elements)
 */
void radiusSearch(
                  const global float * lidar_point_x, const global float * lidar_point_y, const global float * lidar_point_z,
                  const global float * map_points_x, const global float * map_points_y, const global float * map_points_z,
                  const global int * node_indexes, global kdtree_node_petit * root_node, const int n, const int limit, const float radius,
                  global int * neighbor_candidate_indexes, global float * neighbor_candidate_dists)
{
  float reference_point[4];
  reference_point[0] = *lidar_point_x;
  reference_point[1] = *lidar_point_y;
  reference_point[2] = *lidar_point_z;
  reference_point[3] = 0.0;

  global kdtree_node_petit  *previous_node = get_node(root_node, root_node->parent);
  global kdtree_node_petit  *current_node = root_node;
  int neighbors_count = 0;

  float dist;
  // radius search
  while (current_node) {
    int left_index = (current_node->left_right_index >> 16) & 0x0000FFFF;
    int right_index = (current_node->left_right_index) & 0x0000FFFF;
    if ((!current_node->child1) && (!current_node->child2)) {
      // reached to leaf node
      for (int i = left_index; i < right_index; ++i) {
        int index = node_indexes[i];
        float map_point[4];
        map_point[0] = map_points_x[index];
        map_point[1] = map_points_y[index];
        map_point[2] = map_points_z[index];
        map_point[3] = 0.0;
        dist = calculateDistance(reference_point, map_point);
        if (dist < radius) {
          neighbor_candidate_indexes[neighbors_count] = index;
          neighbor_candidate_dists[neighbors_count] = dist;
          neighbors_count++;
        }
        if (neighbors_count >= limit) {
          break;
        }
      }
      // get back to parent node
      previous_node = current_node;
      current_node = get_node(root_node, current_node->parent);
    } else {
      // select best_node, which to be searched first
      float val = reference_point[current_node->axis];
      float diff = val - current_node->axis_val;

      global kdtree_node_petit * best_node;
      global kdtree_node_petit * other_node;

      // not to select NULL child
      if (!current_node->child1) {
        best_node = get_node(root_node, current_node->child2);
        other_node = get_node(root_node, current_node->child1);
      } else if (!current_node->child2) {
        best_node = get_node(root_node, current_node->child1);
        other_node = get_node(root_node, current_node->child2);
      } else {
        // both nodes are not NULL
        if (diff < 0) {
          best_node = get_node(root_node, current_node->child1);
          other_node = get_node(root_node, current_node->child2);
        } else {
          best_node = get_node(root_node, current_node->child2);
          other_node = get_node(root_node, current_node->child1);
        }
      }

      // calc distance, if this is the first time to visit
      if (previous_node == get_node(root_node, current_node->parent)) {
        int index = node_indexes[left_index + (right_index - left_index-1) / 2];
        float map_point[4];
        map_point[0] = map_points_x[index];
        map_point[1] = map_points_y[index];
        map_point[2] = map_points_z[index];
        map_point[3] = 0.0;

        dist = calculateDistance(reference_point, map_point);
        if (dist < radius) {
          neighbor_candidate_indexes[neighbors_count] = index;
          neighbor_candidate_dists[neighbors_count] = dist;
          neighbors_count++;
        }
        if (neighbors_count >= limit) {
          break;
        }

        // move to the child node to be searched first
        previous_node = current_node;
        current_node = best_node;
        continue;
      }
      // return from the search of first child
      if (previous_node == best_node) {
        // test if there are neighbors below the other subtree
        if (fabs(diff) > radius) {
          // if nothing, back to parent node
          previous_node = current_node;
          current_node = get_node(root_node, current_node->parent);
          continue;
        }

        // go to another subtree
        if (other_node) {
          previous_node = current_node;
          current_node = other_node;
          continue;
        } else {
          // enable checked flag, if the other subtree is NULL
          previous_node = other_node;
        }
      }

      // search for both subtrees is completed
      if (previous_node == other_node) {
        // back to parent node
        previous_node = current_node;
        current_node = get_node(root_node, current_node->parent);
        continue;
      }

      // error
      // printf("radiusSearch error\n");
      break;
    }
  }
}

/** \brief Compute derivatives of probability function w.r.t. the transformation vector in OpenCL.
 * \param[in] lidar_points_x 1-D array of x value of input points (n elements)
 * \param[in] lidar_points_y 1-D array of y value of input points (n elements)
 * \param[in] lidar_points_z 1-D array of z value of input points (n elements)
 * \param[in] map_points_x 1-D array of x value of map points (? elements)
 * \param[in] map_points_y 1-D array of y value of map points (? elements)
 * \param[in] map_points_z 1-D array of z value of map points (? elements)
 * \param[in] node_indexes 1-D array of indexes of map points (? elements)
 * \param[in] root_node kdtree root node
 * \param[in] n number of input points
 * \param[in] limit limit of neighbors
 * \param[in] radius search radius
 * \param[out] neighbor_candidate_indexes 1-D array of indexes of neighbor candidates (? elements)
 * \param[out] neighbor_candidate_dists 1-D array of distances of neighbor candidates (? elements)
 * \param[in] map_mean array of mean of map points (? elements)
 * \param[in] map_inverse_cov array of inversed covariance of map points (? elements)
 * \param[in] input_points_x 1-D array of x value of NormalDistributionsTransform.input_ points (n elements)
 * \param[in] input_points_y 1-D array of y value of NormalDistributionsTransform.input_ points (n elements)
 * \param[in] input_points_z 1-D array of z value of NormalDistributionsTransform.input_ points (n elements)
 * \param[in] j_ang 2-D matrix of NormalDistributionsTransform.j_ang
 * \param[in] h_ang 2-D matrix of NormalDistributionsTransform.h_ang
 * \param[out] scores 1-D array of scores of input points (n elements)
 * \param[out] score_gradients 1-D array of gradient scores of input points (n elements)
 * \param[out] hessians 1-D array of hessians of input points (n elements)
 * \param[in] gauss_d1 value of NormalDistributionsTransform.gauss_d1_
 * \param[in] gauss_d2 value of NormalDistributionsTransform.gauss_d2_
 * \param[in] gauss_d3 value of NormalDistributionsTransform.gauss_d3_
 */
kernel void computeDerivativesCL(const global float * lidar_points_x,
                                 const global float * lidar_points_y,
                                 const global float * lidar_points_z,
                                 const global float * map_points_x,
                                 const global float * map_points_y,
                                 const global float * map_points_z,
                                 const global int * node_indexes,
                                 global kdtree_node_petit * root_node,
                                 const int n,
                                 const int limit,
                                 const float radius,
                                 global int * neighbor_candidate_indexes,
                                 global float * neighbor_candidate_dists,
                                 const global float map_mean[][3],
                                 const global float map_inverse_cov[][3][3],
                                 const global float * input_points_x,
                                 const global float * input_points_y,
                                 const global float * input_points_z,
                                 const global float j_ang[8][4],
                                 const global float h_ang[16][4],
                                 global float * scores,
                                 global float score_gradients[][6],
                                 global float hessians[][6][6],
                                 const float gauss_d1,
                                 const float gauss_d2,
                                 const float gauss_d3
                                 )
{
  // 1d range kernel for the point cloud
  int item_index = get_global_id(0);
  if (item_index >= n)
    return;

  // initialize candidates
  for (int i = 0; i < limit; i++) {
    neighbor_candidate_indexes[limit*item_index + i] = -1;
    neighbor_candidate_dists[limit*item_index + i] = 0.0f;
  }

  // Original Point and Transformed Point (for math)
  float x[4];
  float x_trans[3];

  // Initialize Point Gradient and Hessian
  float point_gradient_[4][6];
  float point_hessian_[24][6];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 6; j++) {
      point_gradient_[i][j] = 0.0f;
    }
  }
  for (int i = 0; i < 3; i++) {
    point_gradient_[i][i] = 1.0f;
  }
  for (int i = 0; i < 24; i++) {
    for (int j = 0; j < 6; j++) {
      point_hessian_[i][j] = 0.0f;
    }
  }

  // Find nieghbors (Radius search has been experimentally faster than direct neighbor checking.
  radiusSearch(&lidar_points_x[item_index], &lidar_points_y[item_index], &lidar_points_z[item_index],
               map_points_x, map_points_y, map_points_z, node_indexes, root_node, n, limit, radius,
               (neighbor_candidate_indexes + limit*item_index), (neighbor_candidate_dists + limit*item_index));

  float score_pt = 0;
  float score_gradient_pt[6];
  float hessian_pt[6][6];
  for (int i = 0; i < 6; i++) {
    score_gradient_pt[i] = 0.0f;
  }
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      hessian_pt[i][j] = 0.0f;
    }
  }

  for (int i = 0; i < limit; i++) {
    int candidate_index = neighbor_candidate_indexes[limit*item_index + i];
    if (candidate_index != -1) {
      x[0] = input_points_x[item_index];
      x[1] = input_points_y[item_index];
      x[2] = input_points_z[item_index];
      x[3] = 0.0f;

      // Denorm point, x_k' in Equations 6.12 and 6.13 [Magnusson 2009]
      x_trans[0] = lidar_points_x[item_index] - map_mean[candidate_index][0];
      x_trans[1] = lidar_points_y[item_index] - map_mean[candidate_index][1];
      x_trans[2] = lidar_points_z[item_index] - map_mean[candidate_index][2];

      // Compute derivative of transform function w.r.t. transform vector, J_E and H_E in Equations 6.18 and 6.20
      // [Magnusson 2009]
      computePointDerivatives(x, point_gradient_, point_hessian_, j_ang, h_ang);
      // Update score, gradient and hessian, lines 19-21 in Algorithm 2, according to Equations 6.10, 6.12 and 6.13,
      // respectively [Magnusson 2009]
      score_pt += updateDerivatives(
                                    score_gradient_pt, hessian_pt, point_gradient_, point_hessian_, x_trans, map_inverse_cov[candidate_index],
                                    gauss_d1, gauss_d2, gauss_d3);
    }
  }
  scores[item_index] += score_pt;
  for (int i = 0; i < 6; i++) {
    score_gradients[item_index][i] += score_gradient_pt[i];
  }
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      hessians[item_index][i][j] += hessian_pt[i][j];
    }
  }
  return;
}
