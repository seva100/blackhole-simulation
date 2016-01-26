#include "Tracer.h"
#include "Constants.h"
#include "atlimage.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#define PI 3.1415926536
enum { NO_INTERSN, DISK_INTERSN, SPHERE_INTERSN, OTHER_INTERSN };

//using namespace glm;
std::vector<CImage *> saved_images;
std::vector<std::pair<unsigned, unsigned> > image_shapes;

SRay CTracer::MakeRay(glm::uvec2 pixel_pos, double x_shift=0.5, double y_shift=0.5)
{
	SRay ray;
	ray.m_start = m_camera.m_pos;
	// m_dir = ViewDir + ((cv + 0.5) / W - 0.5) * Right + ((cy + 0.5) / H - 0.5) * Up
	ray.m_dir = m_camera.m_forward + ((pixel_pos.x + x_shift) / m_camera.m_resolution.x - 0.5) * m_camera.m_right + 
		((pixel_pos.y + y_shift) / m_camera.m_resolution.y - 0.5) * m_camera.m_up;
	ray.m_dir /= glm::length(ray.m_dir);
	return ray;
}

template <class T1, class T2>
inline bool is_close(T1 num, T2 to, double tol = 1e-7) {
	return (((num - to) <= tol) && ((num - to) >= -tol));
}

bool crosses_sphere(const Segment &sgm, const Sphere &sph) {
	// Checks if segment intersects sphere surface.
	// x^2+y^2+z^2 (</>)? r^2
	double side1 = glm::length(sgm.start - sph.center);
	double side2 = glm::length(sgm.end - sph.center);
	bool check = ((side1 - sph.radius) * (side2 - sph.radius) <= 0);		// segment start and end are on different sides of sphere
	return check;
}

std::pair<bool, double>    // [intersects or not], line parameter of intersection pnt
crosses_sphere(const SRay &ray, const Sphere &sph) {
	// Checks if ray will ever cross sphere surface.
	// ray: a = a0 + dt, where a, a0, d are vectors
	// sphere: center s0, radius R
	// condition: ||d||^2 t^2 + 2t (d, a0 - s0) + (||a0 - s0||^2 - R^2) = 0;   t > 0

	// coefficients of quadratic equation:
	double a = glm::length(ray.m_dir);
	a *= a;
	double b = 2 * glm::dot(ray.m_dir, ray.m_start - sph.center);
	double c = glm::length(ray.m_start - sph.center);
	c = c * c - sph.radius * sph.radius;

	if (is_close(a, 0)) {
		return std::make_pair(false, 0);
	}
	double sqrtD = sqrt(b * b - 4 * a * c);
	double t[2];
	t[0] = (-b + sqrtD) / (2 * a);
	t[1] = (-b - sqrtD) / (2 * a);
	if ((t[0] > 0) || (t[1] > 0)) {
		return std::make_pair(true, std::min<double>(t[0], t[1]));
	}
	return std::make_pair(false, 0);
}

std::pair<bool, double>	  // [intersects or not], line parameter of intersection pnt
crosses_disk(const SRay &ray, const Disk &disk) {
	// let w be disk normal, a0 - ray start, d - ray dir
	// t = -<w,a0 - s0> / <w, d> 
	double tol = 1e-7;

	double t = 0;
	if (is_close(glm::dot(disk.normal, ray.m_dir), 0)) {
		return std::make_pair(false, t);
	}
	t = -glm::dot(disk.normal, ray.m_start - disk.center) / glm::dot(disk.normal, ray.m_dir);
	if (t <= tol) {
		return std::make_pair(false, t);
	}
	double dist_to_center = glm::length(disk.center - (ray.m_start + t * ray.m_dir));
	if (disk.in_rad < dist_to_center && dist_to_center <= disk.out_rad) {
		return std::make_pair(true, t);
	}
	return std::make_pair(false, t);
}

std::pair<unsigned, unsigned> img_shape(int img_no) {
	CImage* pImage = saved_images[img_no];
	return std::make_pair(pImage->GetWidth(), pImage->GetHeight());
}

double img_get_pxl(int img_no, int i, int j, int channel) {
	// Reading input texture sample
	CImage* pImage = saved_images[img_no];
	unsigned char *pData = (unsigned char*)pImage->GetBits();
	int pitch = pImage->GetPitch();
	int byte_pp = pImage->GetBPP() / 8;
	//fout << "Initial sizeof: " << sizeof(pData) << std::endl;
	if (pitch < 0)
	{
		//fout << "NEGATIVE; pitch = " << pitch << std::endl;
		pData += pitch * (image_shapes[img_no].second - 1);
		//pitch = -pitch;
	}
	/*for (int i = 0; i < pImage->GetHeight(); i++) // Image lines
	{
		for (int j = 0; j < pImage->GetWidth(); j++) // Pixels in line
		{
			unsigned char b = pCurrentLine[j * 4];
			unsigned char g = pCurrentLine[j * 4 + 1];
			unsigned char r = pCurrentLine[j * 4 + 2];
			unsigned char alpha = pCurrentLine[j * 4 + 3];
		}
		pCurrentLine += pitch;
	}*/
	//fout << i << ' ' << pitch << ' ' << j << ' ' << byte_pp << ' ' << channel << std::endl << i * pitch + j * byte_pp + channel <<
	//	' ' << sizeof(pData) << std::endl;
	unsigned char *pxl_addr = (unsigned char *)pImage->GetPixelAddress(i, j); // j, i?
	double pxl = double(*(pxl_addr + channel));
	return pxl;
}

glm::vec4 img_get_pxl_rgba(int img_no, int i, int j) {
	int chnls[] = { 2, 1, 0, 3 };
	float color_ch[4];
	for (int k = 0; k < 4; ++k) {
		color_ch[k] = float(img_get_pxl(img_no, i, j, chnls[k]));
	}
	glm::vec4 ans(color_ch[0], color_ch[1], color_ch[2], color_ch[3]);
	ans /= 255.0;
	return ans;
}

glm::vec3 rgb_cut(const glm::vec4 color) {
	return glm::vec3(color.r, color.g, color.b);
}

glm::vec4 disk_pnt_color(const Disk &disk, const glm::dvec3 &pnt) {
	glm::vec4 color;
	glm::dvec3 diff = pnt - disk.center;
	std::pair<unsigned, unsigned> texture_shape = img_shape(0);
	double ratio = disk.out_rad / disk.in_rad;
	double pic_ratio = 1020.0 / 200.0;		// measured empirically
	diff = diff / (disk.in_rad * pic_ratio);
	//diff = diff / (disk.in_rad * ratio);
	diff.x *= double(texture_shape.first) / 2.0;
	diff.y *= double(texture_shape.second) / 2.0;
	glm::ivec2 pxl_pos(texture_shape.second / 2.0 + diff.x, texture_shape.first / 2.0 + diff.y);
	color = img_get_pxl_rgba(0, pxl_pos.x, pxl_pos.y);
	//fout << color.x << ' ' << color.y << ' ' << color.z << ' ' << color.w << std::endl;
	//return glm::vec3(0, 1, 0);
	return color;
}

std::pair<int, glm::vec4> CTracer::intersection(SRay &ray) {
	glm::vec4 color(0, 0, 1, 0);
	// checking if ray crosses sphere (black hole)
	std::pair<bool, double> sphere_intersn = crosses_sphere(ray, black_hole);
	std::pair<bool, double> disk_intersn = crosses_disk(ray, disk);
	glm::dvec3 intersn_pnt;
	glm::vec4 disk_intersn_color;
	if (disk_intersn.first) {
		intersn_pnt = ray.m_start + ray.m_dir * disk_intersn.second;
		disk_intersn_color = disk_pnt_color(disk, intersn_pnt);
	}
	if (sphere_intersn.first) {
		if (!disk_intersn.first) {		// checking alpha-channel
			//std::cout << "intersects" << std::endl;
			//fout << "ONLY SPHERE" << std::endl;
			color = glm::vec4(0, 0, 0, 0);
			return std::make_pair(SPHERE_INTERSN, color);
		}
		else {
			if (sphere_intersn.second < disk_intersn.second) {
				//fout << "BOTH, sphere first" << std::endl;
				color = glm::vec4(0, 0, 0, 0);
				return std::make_pair(SPHERE_INTERSN, color);
			}
			else {
				//fout << "BOTH, sphere second" << std::endl;
				// Alpha-mixing
				float alpha = disk_intersn_color.w;
				return std::make_pair(DISK_INTERSN, alpha * disk_intersn_color +
					(1 - alpha) * glm::vec4(0, 0, 0, 0));
			}
		}
	}
	for (int i = 0; i < 2; ++i) {
		sphere_intersn = crosses_sphere(ray, planets[i]);
		if (sphere_intersn.first) {
			if (!disk_intersn.first) {		// checking alpha-channel
				//std::cout << "intersects" << std::endl;
				//fout << "ONLY SPHERE" << std::endl;
				color = glm::vec4(0, planet_colors[i]);
				return std::make_pair(SPHERE_INTERSN, color);
			}
			else {
				if (sphere_intersn.second < disk_intersn.second) {
					//fout << "BOTH, sphere first" << std::endl;
					color = glm::vec4(0, planet_colors[i]);
					return std::make_pair(SPHERE_INTERSN, color);
				}
				else {
					//fout << "BOTH, sphere second" << std::endl;
					float alpha = disk_intersn_color.w;
					return std::make_pair(DISK_INTERSN, alpha * disk_intersn_color +
						(1 - alpha) * glm::vec4(0, planet_colors[i]));
				}
			}
		}
	}
	if (disk_intersn.first && disk_intersn_color.w != 0) {		// checking alpha-channel
		//fout << "ONLY DISK" << std::endl;
		return std::make_pair(DISK_INTERSN, disk_intersn_color);
	}
	//fout << "NONE" << std::endl;
	return std::make_pair(NO_INTERSN, color);
}

inline std::pair<int, glm::vec4> CTracer::intersection(Segment &segm) {
	// checking if ray crosses sphere (black hole)
	if (crosses_sphere(segm, black_hole)) {
		return std::make_pair(SPHERE_INTERSN, glm::vec4(0, 0, 0, 1));
	}
	for (int i = 0; i < 2; ++i) {
		if (planet_enable[i] && crosses_sphere(segm, planets[i])) {
			return std::make_pair(SPHERE_INTERSN, glm::vec4(planet_colors[i], 1));
		}
	}
	glm::dvec3 dir = segm.end - segm.start;
	std::pair<bool, double> disk_intersn = crosses_disk(SRay(segm.start, dir), disk);
	if (disk_intersn.first) {
		double t_max = 1;
		// segment is a part of ray SRay(segm.start, dir) with condition 0 <= t <= t_max
		if (disk_intersn.second > t_max) {
			return std::make_pair(NO_INTERSN, glm::vec4(0, 0, 0, 1));
		}
		glm::dvec3 intersn_pnt = segm.start + dir * disk_intersn.second;
		glm::vec4 disk_intersn_color = disk_pnt_color(disk, intersn_pnt);
		return std::make_pair(DISK_INTERSN, disk_intersn_color);
	}
	return std::make_pair(NO_INTERSN, glm::vec4(0, 0, 0, 0));
}

glm::vec4 alpha_blend(glm::vec4 &first, glm::vec4 &second) {
	const double ztol = 1e-7;

	float alpha1 = first.w;
	float alpha2 = second.w;
	float result_alpha = alpha1 + alpha2 * (1 - alpha1);
	glm::vec3 result_rgb;
	if (result_alpha < ztol) {
		result_rgb = glm::vec3(0, 0, 0);
	}
	else {
		result_rgb = (alpha1 * rgb_cut(first) + alpha2 * (1 - alpha1) * rgb_cut(second)) / result_alpha;
	}
	return glm::vec4(result_rgb, result_alpha);
}

glm::vec3 CTracer::TraceRay(SRay ray)
{
	int max_iter = 400;
	int straight_iters_num = 7;
	double time_step = 4;	// in seconds
	double eps = 1e-5;
	double ztol = 1e-7;

	glm::vec4 disk_color(-1, -1, -1, -1);
	bool crossed_disk = false;
	bool crossed_other = false;
	glm::vec3 color(0, 0, 1);
	glm::vec4 other_color;
	
	glm::dvec3 cur_pos(m_camera.m_pos), cur_speed(ray.m_dir / glm::length(ray.m_dir) * light_speed);
	int iter = 0;
	glm::dvec3 cur_dir(0.0, 0.0, 0.0);
	glm::dvec3 prev_dir(0.0, 0.0, 0.0);
	Segment segm;
	SRay changed;
	double cam_to_bh = glm::length(m_camera.m_pos - black_hole.center);
	glm::dvec3 accel;
	int last_iter_straight = -1;
	double cur_pos_abs_sq;
	//while (iter < max_iter) {
	while (true) {
		// simple:
		//cur_pos += ray.m_dir * light_speed * time_step;
		//cur_speed = ray.m_dir / glm::length(ray.m_dir) * light_speed;
		// with acceleration:
		cur_pos_abs_sq = glm::dot(cur_pos, cur_pos);
		accel = -(cur_pos / glm::length(cur_pos)) * grav_const * black_hole.mass / cur_pos_abs_sq;		// assuming coordinates center is in center of black hole.
		//fout << "accel: " << accel.x << ' ' << accel.y << ' ' << accel.z << ", module: " << glm::length(accel) << std::endl;
		segm.start = cur_pos;
		//changed.m_start = cur_pos;
		cur_pos += (cur_speed + accel * time_step / 2.0) * time_step;
		cur_speed += accel * time_step;
		cur_speed *= light_speed / glm::length(cur_speed);
		// choosing non-uniform time step
		//time_step = 2 / glm::length(accel);
	
		// Checking intersection for segment.
		segm.end = cur_pos;
		auto intersn = intersection(segm);
		if (intersn.first == SPHERE_INTERSN) {
			other_color = intersn.second;
			crossed_other = true;
			if (!alpha_blending_enable) {
				return rgb_cut(other_color);
			}
			break;
		}
		else if (intersn.first == DISK_INTERSN
			&& intersn.second.w > ztol) {
			if (!crossed_disk) {
				disk_color = intersn.second;
				crossed_disk = true;
				if (!alpha_blending_enable) {
					return rgb_cut(disk_color);
				}
			}
			else {
				// encountered disk previously. 
				// It is impossible to come to this point if alpha blending is disabled.
				disk_color = alpha_blend(disk_color, intersn.second);
			}
		}
		
		// another way is to compute cross product
		//if (glm::dot(cur_pos, cur_pos) > 3 * (cam_to_bh + black_hole.radius) * (cam_to_bh + black_hole.radius) ||
		//if	(glm::dot(cur_speed, prev_dir) < (1 - eps) * light_speed * light_speed) {		// speeds must be of length light_speed
									    														// has put light_speed^2 in the right side of equation
		if (glm::length(glm::cross(cur_speed, prev_dir)) < eps * light_speed * light_speed &&
			glm::dot(cur_pos, cur_pos) > 2 * (cam_to_bh + black_hole.radius) * (cam_to_bh + black_hole.radius)) {
			if (last_iter_straight == -1) {
				last_iter_straight = iter;
			}
			else if (iter - last_iter_straight > straight_iters_num) {
				glm::dvec2 spher(0, 0);
				glm::dvec3 dir_normd = cur_speed / glm::length(cur_speed);
				spher.x = atan2(dir_normd.x, dir_normd.y) + PI;		// (at least in MS library) atan2 returns value in [-pi, pi]
				spher.y = asin(dir_normd.z);
				//fout << "final iter no: " << iter << std::endl;
				//fout << "spher coords: " << spher.y << ' ' << spher.x << std::endl;
				std::pair<unsigned, unsigned> stars_shape = img_shape(1);
				double H = stars_shape.second;
				spher.x *= H / PI;	// transforming [0, 2pi) to [0, 2H)
				spher.y = (spher.y + PI / 2.0) * H / PI;		// transforming [-pi/2, pi/2] to [0, H]
				spher.x = glm::clamp(spher.x, 0.0, 2 * H - 1);
				spher.y = glm::clamp(spher.y, 0.0, H - 1);
				//fout << "ray start: " << ray.m_start.x << ' ' << ray.m_start.y << ' ' << ray.m_start.z << std::endl;
				//fout << "ray dir: " << ray.m_dir.x << ' ' << ray.m_dir.y << ' ' << ray.m_dir.z << std::endl;
				//fout << "spher coords: " << spher.y << ' ' << spher.x << std::endl;
				
				other_color = img_get_pxl_rgba(1, int(spher.x), int(spher.y));
				other_color.w = 1;
				crossed_other = true;
				break;
			}
		}
		prev_dir = cur_speed;
		
		iter++;
	}
	//fout << "min dist to black hole: " << *min_element(dist_to_bh.begin(), dist_to_bh.end()) << std::endl;
	//std::cout << alpha_blending_enable << std::endl;
	if (!alpha_blending_enable) {
		color = rgb_cut(other_color);
	}
	else {
		if (crossed_disk) {
			if (!crossed_other) {
				color = rgb_cut(disk_color);
			}
			else {
				// Alpha-blending
				color = rgb_cut(alpha_blend(disk_color, other_color));
			}
		}
		else {
			color = rgb_cut(other_color);
		}
	}
	
	return color;
}

void CTracer::RenderImage(Params & params)
{
	// Reading input texture sample
	CImage* pImage = LoadImageFromFile("data/disk_32.png");
	if (!pImage) {
		std::cout << "Received NULL pointer when loading texture. Probably wrong file." << std::endl;
		return;
	}
	saved_images.push_back(pImage);
	image_shapes.push_back(img_shape(0));
	// Reading background (stars) texture
	CImage* stars = LoadImageFromFile("data/stars.jpg");
	if (!stars) {
		std::cout << "Received NULL pointer when loading texture. Probably wrong file." << std::endl;
		return;
	}
	saved_images.push_back(stars);
	image_shapes.push_back(img_shape(1));

	// Filling in properties
	m_camera.m_resolution = glm::uvec2(params.xRes, params.yRes);
	m_camera.m_pixels.resize(params.xRes * params.yRes);
	m_camera.m_pos = params.camera_pos;
	params.up /= glm::length(params.up);
	params.up *= params.yRes;
	m_camera.m_up = params.up;
	params.right /= glm::length(params.right);
	params.right *= params.xRes;
	m_camera.m_right = params.right;
	params.view_dir /= glm::length(params.view_dir);
	params.view_dir *= params.yRes / (2.0 * tan(params.view_angle.y / 2.0));
	m_camera.m_forward = params.view_dir;
	m_camera.m_viewAngle = params.view_angle;
	double black_hole_radius = 2 * grav_const * params.black_hole_mass / light_speed / light_speed;
	black_hole.center = glm::dvec3(0, 0, 0);
	black_hole.radius = black_hole_radius;
	black_hole.mass = params.black_hole_mass;
	// setting disk params
	disk.center = glm::dvec3(0, 0, 0);
	disk.in_rad = black_hole_radius;
	disk.out_rad = black_hole_radius * params.disk_bh_rad_ratio;
	disk.normal = glm::vec3(0, 0, 1);
	for (int i = 0; i < 2; ++i) {
		planet_enable[i] = params.planet_enable[i];
		planets[i].center = params.planet_center[i];
		planets[i].radius = params.planet_rad[i];
		planet_colors[i] = params.planet_color[i];
	}
	alpha_blending_enable = params.alpha_blending_enable;
	antialiasing_rays = params.antialiasing_rays;
	
	// Rendering
	double gap, x_shift, y_shift;
	glm::vec3 rays_accum(0, 0, 0);
	SRay ray;
	for (int i = 0; i < params.yRes; i++) {
		for (int j = 0; j < params.xRes; j++)
		{
			rays_accum.r = rays_accum.b = rays_accum.g = 0;
			x_shift = y_shift = gap = 1.0 / (antialiasing_rays + 1);
			for (int x_rays = 0; x_rays < antialiasing_rays; ++x_rays, x_shift += gap) {
				for (int y_rays = 0; y_rays < antialiasing_rays; ++y_rays, y_shift += gap) {
					ray = MakeRay(glm::uvec2(j, i), x_shift, y_shift);
					rays_accum += TraceRay(ray);
				}
			}
			rays_accum /= (antialiasing_rays * antialiasing_rays);
			m_camera.m_pixels[i * params.xRes + j] = rays_accum;
			//m_camera.m_pixels[i * params.xRes + j] = rgb_cut(img_get_pxl_rgba(1, i, j)) / 255.0f;
		}
	}
}

void CTracer::SaveImageToFile(std::string fileName)
{
	CImage image;

	int width = m_camera.m_resolution.x;
	int height = m_camera.m_resolution.y;

	image.Create(width, height, 24);

	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch = -pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			glm::vec3 color = m_camera.m_pixels[textureDisplacement + j];

			// clamp(val, minVal, maxVal) = min(max(x, minVal), maxVal)
			// shows that val is in [minVal, maxVal]
			imageBuffer[imageDisplacement + j * 3] = unsigned(glm::clamp(color.b, 0.0f, 1.0f) * 255.0f);
			imageBuffer[imageDisplacement + j * 3 + 1] = unsigned(glm::clamp(color.g, 0.0f, 1.0f) * 255.0f);
			imageBuffer[imageDisplacement + j * 3 + 2] = unsigned(glm::clamp(color.r, 0.0f, 1.0f) * 255.0f);
		}

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

	image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
	CImage* pImage = new CImage;

	if (SUCCEEDED(pImage->Load(fileName.c_str())))
		return pImage;
	else
	{
		delete pImage;
		return NULL;
	}
}