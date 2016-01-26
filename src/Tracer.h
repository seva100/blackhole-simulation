#pragma once

#include "glm/glm.hpp"
#include "Types.h"
#include "Scene.h"

#include "string"
#include "atlimage.h"
#include <vector>
#include <utility>

class CTracer
{
public:
  SRay MakeRay(glm::uvec2 pixel_pos, double x_shift, double y_shift);  // Create ray for specified pixel
  glm::vec3 TraceRay(SRay ray); // Trace ray, compute its color
  void RenderImage(Params &params);
  void SaveImageToFile(std::string fileName);
  CImage* LoadImageFromFile(std::string fileName);
  std::pair<int, glm::vec4> intersection(SRay &ray);
  std::pair<int, glm::vec4> CTracer::intersection(Segment &segm);

public:
  SCamera m_camera;
  CScene* m_pScene;
  // objects
  Sphere black_hole;
  Disk disk;
  bool planet_enable[2];
  Sphere planets[2];
  glm::vec3 planet_colors[2];
  // properties
  bool alpha_blending_enable;
  int antialiasing_rays;
};