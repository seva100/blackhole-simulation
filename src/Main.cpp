#include "Tracer.h"
#include "Constants.h"
#include "stdio.h"
#include "inih/cpp/INIReader.h"
#include <cmath>
#include <iostream>

#define PI 3.1415926536

void main(int argc, char** argv)
{
	CTracer tracer;
	CScene scene;
	Params params;

	if (argc == 2) // There is input file in parameters
	{
		INIReader reader(argv[1]);
		if (reader.ParseError() < 0) {
			std::cout << "Can't load configuration file.\n";
			return;
		}

		params.xRes = reader.GetInteger("screen", "xRes", 1024);
		params.yRes = reader.GetInteger("screen", "yRes", 768);
		params.black_hole_mass = reader.GetReal("black hole", "mass", 1e37);
		double black_hole_radius = 2 * grav_const * params.black_hole_mass / light_speed / light_speed;
		// Black hole lies in (0, 0, 0)
		params.camera_pos = glm::dvec3(reader.GetReal("camera pos", "x", 8e10),
			reader.GetReal("camera pos", "y", 0), 
			reader.GetReal("camera pos", "z", 1e10));
		params.up = glm::dvec3(reader.GetReal("up", "x", -1), 
			reader.GetReal("up", "y", 0), 
			reader.GetReal("up", "z", 6));
		params.right = glm::dvec3(reader.GetReal("right", "x", 0),
			reader.GetReal("right", "y", 1),
			reader.GetReal("right", "z", 0));
		params.view_dir = glm::dvec3(reader.GetReal("view_dir", "x", -6), 
			reader.GetReal("view_dir", "y", 0), 
			reader.GetReal("view_dir", "z", -1));
		params.view_angle = glm::dvec2(reader.GetReal("view_angle", "x", PI / 2.0), 0);
		params.view_angle.y = 2 * atan(double(params.yRes) / double(params.xRes) * tan(params.view_angle.x / 2));
		params.disk_bh_rad_ratio = reader.GetInteger("disk", "radius ratio", 2);
		//params.disk_normal = glm::dvec3(0, 0, 1);
		params.planet_enable[0] = reader.Get("planet 0", "enable", "yes") == "yes";
		params.planet_enable[0] = reader.Get("planet `", "enable", "yes") == "yes";
		params.planet_center[0] = glm::dvec3(reader.GetReal("planet 0", "center x", 0),
			reader.GetReal("planet 0", "center y", 0),
			reader.GetReal("planet 0", "center z", 0));
		params.planet_center[1] = glm::dvec3(reader.GetReal("planet 1", "center x", 0),
			reader.GetReal("planet 1", "center y", 0),
			reader.GetReal("planet 1", "center z", 0));
		params.planet_rad[0] = reader.GetReal("planet 0", "radius", 0);
		params.planet_rad[1] = reader.GetReal("planet 1", "radius", 0);
		params.planet_color[0] = glm::vec3(reader.GetReal("planet 0", "color r", 0),
			reader.GetReal("planet 0", "color b", 0),
			reader.GetReal("planet 0", "color g", 0)) / 255.0f;
		params.planet_color[1] = glm::vec3(reader.GetReal("planet 1", "color r", 0),
			reader.GetReal("planet 1", "color b", 0),
			reader.GetReal("planet 1", "color g", 0)) / 255.0f;
		params.alpha_blending_enable = int(reader.Get("alpha blending", "enable", "yes") == "yes");
		params.antialiasing_rays = reader.GetInteger("antialiasing", "rays through pixel", 1);
	}
	else {
		std::cout << "No config specified. Using default values." << std::endl;
		params.xRes = 800;  // Default resolution
		params.yRes = 600;
		params.black_hole_mass = 1e37; // Default mass
		double black_hole_radius = 2 * grav_const * params.black_hole_mass / light_speed / light_speed;
		// Black hole lies in (0, 0, 0)
		params.camera_pos = glm::vec3(8e10, 0, 4e9);		// Default camera position
		params.up = glm::vec3(-1, 0, 6);					// Default vectors
		params.right = glm::vec3(0, 1, 0);
		params.view_dir = glm::vec3(-6, 0, -1);
		params.view_angle = glm::vec2(PI / 2.0, 0.0);		// Default angle
		params.view_angle.y = 2 * atan(double(params.yRes) / double(params.xRes) * tan(params.view_angle.x / 2));
		params.disk_bh_rad_ratio = 4;
		//params.disk_normal = glm::dvec3(0, 0, 1);
		params.planet_enable[0] = params.planet_enable[1] = true;
		params.planet_center[0] = glm::dvec3(2e10, black_hole_radius * 2, -black_hole_radius);
		params.planet_center[1] = glm::dvec3(3e10, -black_hole_radius * 1.5, -black_hole_radius);
		params.planet_rad[0] = black_hole_radius / 3;
		params.planet_rad[1] = black_hole_radius / 5;
		params.planet_color[0] = glm::vec3(222, 155, 40) / 255.0f;
		params.planet_color[1] = glm::vec3(100, 150, 100) / 255.0f;
		params.alpha_blending_enable = 1;
		params.antialiasing_rays = 2;
	}
	tracer.m_pScene = &scene;
	tracer.RenderImage(params);
	tracer.SaveImageToFile("Result.png");
}