#include <resolver.hpp>

#include <rt.hpp>


using namespace rt;

void rt::render(const scene_t* scn)
{
	unsigned int w=scn->img->get_width();
	unsigned int h=scn->img->get_height();
	unsigned int samples = scn->img->get_no_of_samples();
#pragma omp parallel for
	for (unsigned int i=0; i<w; i++)
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples, 100. * i / (w - 1));
#pragma omp parallel for
		for (unsigned int j=0; j<h; j++)
		{
			ray_t ray;
			int d=0;
			std::vector<Eigen::Vector2f> psample=scn->img->sample_pixels(i,j);
			color_t col(0.0,0.0,0.0);
#pragma omp parallel for reduction(+:col)

			for(unsigned int k=0;k<samples;k++){
				scn->cam->sample_ray(ray, psample[k]);
				col += scn->intg->radiance(scn, ray, d);
			}
			col /= samples;
			// col *= scn->intg->radiance(scn, ray, d);

			scn->img->set_pixel(i,j,col);
		}
	}
	printf("\n");
}

int main(int argc, char **argv)
{
	if (argc != 2) 
		{
        	std::cerr << "Syntax: " << argv[0] << " <scene.xml>" << std::endl;
        	return -1;
		}

	filesystem::path path(argv[1]);

	try
	{
		if (path.extension() == "xml") 
		{
			std::string scene_filename(argv[1]);
  			rt::scene_t scn(scene_filename);

  			rt::render(&scn);

			std::string img_filename = scene_filename;
    		size_t lastdot = img_filename.find_last_of(".");
    		if (lastdot != std::string::npos)
     			img_filename.erase(lastdot, std::string::npos);
   			img_filename += ".ppm";
  	
  			scn.img->write(img_filename);
  		}
  		else
  		{
  			std::cerr<<"Error: Unknown file type."<<std::endl;
  			return -1;
  		}
  } 
  catch (const std::exception &e)
  {
  	std::cerr<<"Error: "<<e.what()<<std::endl;
  	return -1;
  }
  	
  return 0;
}
