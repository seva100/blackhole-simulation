#pragma once

#include "glm/glm.hpp"
#include "vector"

struct Params {
	int xRes;
	int yRes;
	double black_hole_mass;
	// Let's use left vector set (x ->, y from us, z down to up)
	glm::dvec3 camera_pos;
	glm::dvec3 up;
	glm::dvec3 right;
	glm::dvec3 view_dir;
	glm::dvec2 view_angle;
	double disk_bh_rad_ratio;
	//glm::dvec3 disk_normal;
	bool planet_enable[2];
	glm::dvec3 planet_center[2];
	double planet_rad[2];
	glm::vec3 planet_color[2];
	int alpha_blending_enable;
	int antialiasing_rays;
};

struct Segment {
	glm::dvec3 start;
	glm::dvec3 end;
};

struct Sphere {
	glm::dvec3 center;
	double radius;
	double mass;	// useful for black hole
};

struct Disk {
	glm::dvec3 center;
	double in_rad, out_rad;
	glm::dvec3 normal;
};

struct SRay
{
  SRay()  {}
  SRay(glm::dvec3 start, glm::dvec3 dir) : m_start(start), m_dir(dir) {};
  glm::dvec3 m_start;
  glm::dvec3 m_dir;
};

struct SCamera
{
  glm::dvec3 m_pos;          // Camera position and orientation
  glm::dvec3 m_forward;      // Orthonormal basis
  glm::dvec3 m_right;
  glm::dvec3 m_up;

  glm::dvec2 m_viewAngle;    // View angles, rad
  glm::uvec2 m_resolution;  // Image resolution: w, h

  std::vector<glm::vec3> m_pixels;  // Pixel array
};

struct SMesh
{
  std::vector<glm::dvec3> m_vertices;  // vertex positions
  std::vector<glm::uvec3> m_triangles;  // vetrex indices
};