#include <iostream>
#include <integrator.hpp>
#include <math.h>
#define EPSILON 0.00001
using namespace rt;

color_t whitted_integrator_t::radiance(const scene_t* _scn, ray_t& _ray, int d) const
{
	bool found_intersection=false;
	color_t d_col(0.0f);
	if(d >= _scn->intg->depth) return _scn->img->get_bgcolor();


	std::vector<object_t*>::const_iterator oit;
	hit_t hit, minhit;
	Eigen::Vector3d hitpt, normal;

	for (oit=_scn->objs.begin(); oit!=_scn->objs.end(); oit++)
	{
		if ((*oit)->intersect(hit, _ray))
		{
		  _ray.maxt = hit.second;
		  minhit = hit;
		  
		  hitpt = _ray.origin+_ray.maxt*_ray.direction;
		  normal = (*oit)->get_normal(hitpt);

		  found_intersection=true;
		}
		
	}
	
	
	if(found_intersection)
	{
		std::list<light_t*>::const_iterator lit;
		for (lit = _scn->lits.begin(); lit != _scn->lits.end(); lit++)
		{
			d_col += (*lit)->direct(hitpt, normal, minhit.first->get_material(), _scn, _ray);
		}

		double eta = minhit.first->get_material()->get_eta();

		bool can_transmit = minhit.first->get_material()->get_is_transmit(), can_reflect = minhit.first->get_material()->get_is_reflect();
		if(_ray.direction.dot(normal) > 0){
			normal = -normal;
		}else{
			eta = 1.0 / eta;
		}

		if(can_reflect && can_transmit){
			Eigen::Vector3d ref_ray = (_ray.direction - 2 * (_ray.direction.dot(normal)) * normal).normalized();

			double TIR_check = 1 - pow(eta, 2) * (1 - pow(normal.dot(_ray.direction), 2));
			Eigen::Vector3d transmit_ray = (TIR_check > 0) ? (eta * _ray.direction + (-eta * normal.dot(_ray.direction) - sqrt(TIR_check)) * normal).normalized() : ref_ray;

			ray_t reflected_ray(hitpt + normal * EPSILON, ref_ray), transmitted_ray(hitpt - normal * EPSILON, transmit_ray);
			d_col += minhit.first->get_material()->get_transmit() * radiance(_scn, transmitted_ray, d + 1);
		}
		else if(can_reflect){
			Eigen::Vector3d ref_ray = (_ray.direction - 2 * (_ray.direction.dot(normal)) * normal).normalized();
			ray_t reflected_ray(hitpt + normal * EPSILON ,ref_ray);
			d_col += minhit.first->get_material()->get_reflect() * radiance(_scn,reflected_ray,d+1);
		}
	}

	return d_col;
}

color_t path_integrator_t::radiance(const scene_t *_scn, ray_t &_ray, int d) const{
	bool found_intersection = false;
	color_t d_col(0.0f);
	if (d >= _scn->intg->depth)
		return _scn->img->get_bgcolor();

	std::vector<object_t *>::const_iterator oit;
	hit_t hit, minhit;
	Eigen::Vector3d hitpt, normal;

	for (oit = _scn->objs.begin(); oit != _scn->objs.end(); oit++)
	{
		if ((*oit)->intersect(hit, _ray))
		{
			_ray.maxt = hit.second;
			minhit = hit;

			hitpt = _ray.origin + _ray.maxt * _ray.direction;
			normal = (*oit)->get_normal(hitpt);

			found_intersection = true;
		}
	}
	hit_t hit1;

	std::list<light_t *>::const_iterator lit;
	for(lit = _scn->lits.begin(); lit != _scn->lits.end(); lit++)
	{
		if ((*lit)->intersect(hit1, _ray))
		{
			return (*lit)->get_color();
		}
	}
	
	if (found_intersection)
	{
		color_t f = minhit.first->get_material()->get_diffuse();
		double p = f.x() > f.y() && f.x() > f.z() ? f.x() : f.y() > f.z() ? f.y() : f.z();
		if (d > 5)
		{
			if (drand48() < p)
				f = f * (1 / p);
			else
				return color_t(0.0);
		}
		for (lit = _scn->lits.begin(); lit != _scn->lits.end(); lit++)
		{
			d_col += (*lit)->direct(hitpt, normal, minhit.first->get_material(), _scn, _ray);
		}

		double eta = minhit.first->get_material()->get_eta();

		bool can_transmit = minhit.first->get_material()->get_is_transmit(), can_reflect = minhit.first->get_material()->get_is_reflect();
		if (_ray.direction.dot(normal) > 0)
		{
			normal = -normal;
		}
		else
		{
			eta = 1.0 / eta;
		}
		if (!can_reflect && !can_transmit)
		{
			double u = drand48();
			double v = drand48();
			Eigen::Vector3d relative_ray = SampleFromHemisphere(u,v);
			Eigen::Vector3d Nt,Nb;
			createCoordinateSystem(normal, Nt, Nb);
			Eigen::Vector3d sample_direction(
				relative_ray[0] * Nb[0] + relative_ray[1] * normal[0] + relative_ray[2] * Nt[0],
				relative_ray[0] * Nb[1] + relative_ray[1] * normal[1] + relative_ray[2] * Nt[1],
				relative_ray[0] * Nb[2] + relative_ray[1] * normal[2] + relative_ray[2] * Nt[2]
			);
			ray_t sample_ray(hitpt+normal*EPSILON, sample_direction);
			d_col += f*radiance(_scn,sample_ray, d+1)*M_1_PI*abs(sample_ray.direction.dot(normal))*2*M_PI;
		}
		else if (can_reflect && !can_transmit)
		{
			Eigen::Vector3d ref_ray = (_ray.direction - 2 * (_ray.direction.dot(normal)) * normal).normalized();
			ray_t reflected_ray(hitpt + normal * EPSILON, ref_ray);
			d_col += minhit.first->get_material()->get_reflect() * radiance(_scn, reflected_ray, d + 1);
			return d_col;
		}
		else if (can_transmit){//dielectric like glass using frenel reflectance
			Eigen::Vector3d ref_ray = (_ray.direction - 2 * (_ray.direction.dot(normal)) * normal).normalized();

			double TIR_check = 1 - pow(eta, 2) * (1 - pow(normal.dot(_ray.direction), 2));
			Eigen::Vector3d transmit_ray = (TIR_check > 0) ? (eta * _ray.direction + (-eta * normal.dot(_ray.direction) - sqrt(TIR_check)) * normal).normalized() : ref_ray;

			ray_t reflected_ray(hitpt + normal * EPSILON, ref_ray), transmitted_ray(hitpt - normal * EPSILON, transmit_ray);
			Eigen::Vector3d n = minhit.first->get_normal(hitpt);
			Eigen::Vector3d nl = (n.dot(_ray.direction) < 0) ? n : -1 * n;
			bool into = n.dot(nl) > 0; // Ray from outside going in?
			double nc = 1, nt = minhit.first->get_material()->get_eta();
			double ddn = _ray.direction.dot(nl);

			if (TIR_check<=0)
			{
				d_col += f * radiance(_scn, transmitted_ray, d + 1);
				return d_col;
			}

			Vector3d tdir = transmitted_ray.direction;

			double a = nt - nc, b = nt + nc;
			double R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
			double Re = R0 + (1 - R0) * c * c * c * c * c;
			double Tr = 1 - Re, P = .25 + .5 * Re;
			double RP = Re / P, TP = Tr / (1 - P);

			if (d > 2)
			{
				if (drand48() < P)
					d_col += f * RP * radiance(_scn, reflected_ray, d + 1);
				else
					d_col += f * TP * radiance(_scn, transmitted_ray, d + 1);
			}
			else
			{
				d_col += f * (Re * radiance(_scn, reflected_ray, d + 1) +
							  Tr * radiance(_scn, transmitted_ray, d + 1));
			}
		}
	}

	return d_col;
}