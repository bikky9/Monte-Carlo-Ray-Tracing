#include <light.hpp>
#include <math.h>
#include <iostream>
using namespace rt;

#define max(a,b) ((a) > (b))?(a):(b)

light_t::light_t() {}
light_t::~light_t() { }


point_light_t::point_light_t(const Vector3d& _pos, const Vector3d& _col, const double _ka): pos(_pos), col(_col), ka(_ka) 
{ }

point_light_t::~point_light_t()
{ }
bool point_light_t::intersect(hit_t& result, const ray_t& _ray) const{
	return false;
}
color_t point_light_t::direct(const Vector3d& hitpt, const Vector3d& normal, const material_t* mat, const scene_t* scn,ray_t& ray) const
{
	Eigen::Vector3d litvec = (pos - hitpt).normalized() , Kd = mat->get_diffuse(),Ks = mat->get_specular();
	Eigen::Vector3d Ambient_comp = ka*col;
	Eigen::Vector3d Local_illumination = Ambient_comp;

	bool shadowed = false;
	std::vector<object_t*>::const_iterator oit;
	hit_t hit;
	for (oit=scn->objs.begin(); oit!=scn->objs.end(); oit++)
	{	
		ray_t shadow_ray(hitpt,litvec);
		if ((*oit)->intersect(hit,shadow_ray))
		{
			shadow_ray.maxt = hit.second;
			if(shadow_ray.maxt * shadow_ray.maxt < (pos - hitpt).dot(pos - hitpt)){
				shadowed = true;
				break;
			}	
		}
	}
	if(!shadowed){

		Eigen::Vector3d temp(Kd[0]*col[0],Kd[1]*col[1],Kd[2]*col[2]);
		double tempk = max(0,litvec.dot(normal));
		Eigen::Vector3d Diffuse_comp = temp * tempk/(litvec.norm()*normal.norm());

		Eigen::Vector3d refected_ray = 2*(litvec.dot(normal))*normal - litvec;
		
		Eigen::Vector3d view_dir = ray.origin - hitpt;
		double s = mat->get_shininess() , dotp = refected_ray.dot(view_dir)/(refected_ray.norm()*view_dir.norm()) ;
		double tempk2 = max(0,dotp);

		Eigen::Vector3d temp2(Ks[0]*col[0],Ks[1]*col[1],Ks[2]*col[2]);
		Eigen::Vector3d Specular_comp = tempk * tempk2 * temp2 * pow(dotp,s-1)/litvec.dot(normal);
		Local_illumination += Diffuse_comp + Specular_comp;
	}

	color_t colo(Local_illumination[0],Local_illumination[1],Local_illumination[2]);
	return colo;
}
		

void point_light_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Light Properties: -------------------------------"<<std::endl;
	stream<<"Type: Point Light"<<std::endl;
	stream<<"Position: "<<pos.format(CommaInitFmt)<<std::endl;
	stream<<"Color: "<<col.format(CommaInitFmt)<<std::endl;
	stream<<"Ambient Coefficient: "<<ka<<std::endl<<std::endl;
}

area_light_t::area_light_t(const Vector3d &_pos, const Vector3d &_col, const double _ka, double _radius) : pos(_pos), col(_col), ka(_ka), radius(_radius)
{
}

area_light_t::~area_light_t()
{
}
bool area_light_t::intersect(hit_t &result, const ray_t &_ray)const
{	
	std::string name = "light";
	color_t col(0);
	simplemat_t mat(name,col,col,col,col,0.0,0.0,false,false);
	material_t* matptr = &mat;
	sphere_t lightobj(matptr,pos,radius);
	return lightobj.intersect(result, _ray);
}
std::vector<Eigen::Vector3d> area_light_t::sample_from_light(int samples) const{
	std::vector<Eigen::Vector3d> result;
	for(int i=0;i<samples;i++){
		double u = drand48();
		double v = drand48();
		Eigen::Vector3d relative_ray = SampleFromHemisphere(u, v);
		relative_ray[2] *= -1;
		result.push_back(pos+relative_ray*radius);
	}
	return result;
}
color_t area_light_t::direct(const Vector3d &hitpt, const Vector3d &normal, const material_t *mat, const scene_t *scn, ray_t &ray) const
{	
	int samples = 1;
	std::vector<Eigen::Vector3d> sample_lights = area_light_t::sample_from_light(samples);
	color_t ans(0.0);
	for(int i=0;i<samples;i++){
		Eigen::Vector3d litpos = sample_lights[i];
		Eigen::Vector3d litvec = (litpos - hitpt).normalized(), Kd = mat->get_diffuse(), Ks = mat->get_specular();
		Eigen::Vector3d Ambient_comp(ka, ka, ka);
		Eigen::Vector3d Local_illumination = Ambient_comp;

		bool shadowed = false;
		std::vector<object_t *>::const_iterator oit;
		hit_t hit;
		for (oit = scn->objs.begin(); oit != scn->objs.end(); oit++)
		{
			ray_t shadow_ray(hitpt, litvec);
			if ((*oit)->intersect(hit, shadow_ray))
			{
				shadow_ray.maxt = hit.second;
				if (shadow_ray.maxt * shadow_ray.maxt < (litpos - hitpt).dot(litpos - hitpt))
				{
					shadowed = true;
					break;
				}
			}
		}
		if (!shadowed)
		{

			Eigen::Vector3d temp(Kd[0] * col[0], Kd[1] * col[1], Kd[2] * col[2]);
			double tempk = max(0, litvec.dot(normal));
			Eigen::Vector3d Diffuse_comp = temp * tempk / (litvec.norm() * normal.norm());

			Eigen::Vector3d refected_ray = 2 * (litvec.dot(normal)) * normal - litvec;

			Eigen::Vector3d view_dir = ray.origin - hitpt;
			double s = mat->get_shininess(), dotp = refected_ray.dot(view_dir) / (refected_ray.norm() * view_dir.norm());
			double tempk2 = max(0, dotp);

			Eigen::Vector3d temp2(Ks[0] * col[0], Ks[1] * col[1], Ks[2] * col[2]);
			Eigen::Vector3d Specular_comp = tempk * tempk2 * temp2 * pow(dotp, s - 1) / litvec.dot(normal);
			Local_illumination += Diffuse_comp + Specular_comp;
		}

		color_t colo(col[0] * Local_illumination[0], col[1] * Local_illumination[1], col[2] * Local_illumination[2]);
		ans+=colo;
	}
	return ans/samples;
}

void area_light_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

	stream << "Light Properties: -------------------------------" << std::endl;
	stream << "Type: Area Light" << std::endl;
	stream << "Position: " << pos.format(CommaInitFmt) << std::endl;
	stream << "Color: " << col.format(CommaInitFmt) << std::endl;
	stream << "Ambient Coefficient: " << ka << std::endl;
	stream << "Radius: " << radius << std::endl
		   << std::endl;
}