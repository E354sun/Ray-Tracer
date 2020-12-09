// Winter 2019

#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include "Ray.hpp"
#include "Mesh.hpp"
#include "GeometryNode.hpp"
#include "PhongMaterial.hpp"
#include "glassMaterial.hpp"
#include "A4.hpp"
#include "OPT_A5.hpp"
#include <cmath>
#include <cstdlib>

using namespace std;

static void printMat(glm::mat4 &M) {

    cout << "Matrix begin" << endl;
    for(int i=0; i<4; i++) {
    	cout << glm::to_string(M[i]) << endl;
    }

    cout << "Matrix end" << endl; 
}  	

glm::mat4 DCtoWorld (
	const double fovy,
	const double width,
	const double height,
	const glm::vec3 &eye,
	const glm::vec3 &view,
	const glm::vec3 &up)
{
    const double d = glm::distance(eye, view);
    // const double d = 1;
    const double Proj_h = 2 * d * glm::tan(glm::radians(fovy/2));
    const double Proj_w = (width / height) * Proj_h;

    // my implement of Lookat 
    glm::mat4 VtoW = glm::inverse(glm::lookAt(eye, view, up));

    cout << "Good: " << endl;
    printMat(VtoW);

    glm::vec3 v = glm::normalize(view);
    glm::vec3 u = glm::cross(glm::normalize(up), v);
    glm::vec3 w = glm::cross(u, v);

    glm::mat4 VtoW_bad;
    VtoW_bad[0] = glm::vec4(u, 0);
    VtoW_bad[1] = glm::vec4(w, 0);
    VtoW_bad[2] = glm::vec4(v, 0);
    VtoW_bad[3] = glm::vec4(eye, 1);

    cout << "Bad: " << endl;
    printMat(VtoW_bad);
    
    // from screen 
    glm::mat4 StoV;
    StoV = glm::translate(StoV, glm::vec3(-width/2, -height/2, -d));
    StoV = glm::scale(glm::mat4(), 
    		glm::vec3(Proj_w/width, -Proj_h/height, 1)) * StoV;

    return VtoW * StoV;
}

// Reference: https://stackoverflow.com/questions/15846867/glossy-reflection-in-ray-tracing
glm::vec3 pur(glm::vec3 ray, float exponent){	
	assert(exponent >= 1);
	
	glm::vec3 e_0 {0, 1, 0};
	glm::vec3 e_1 {0, 0, 1};
	
	glm::vec3 w = glm::normalize(ray);
	
	glm::vec3 u = glm::cross(w, e_0);
	if (glm::length(u) < 0.1) {
		u = glm::cross(w, e_1);
	}
	glm::vec3 v = glm::cross(w, u);
	
	// Reference: https://stackoverflow.com/questions/33986375/generate-random-between-0-1
	std::random_device generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	
	float phi = distribution(generator) * 2 * glm::pi<double>();
    float cosPhi = glm::cos(phi);
	float sinPhi = glm::sin(phi);

	float cosTheta = glm::pow(distribution(generator), 1.0f/(exponent+1));
	float sinTheta = glm::sqrt(1 - std::pow(cosTheta, 2));
	
	glm::vec3 A = w*cosTheta + u*cosPhi*sinTheta + v*sinPhi*sinTheta;

	return glm::normalize(A);
}

std::vector<glm::vec3> purturbedRays(glm::vec3 R, glm::vec3 normal, float exponent, int size)
{
	assert(size > 0);
	std::vector<glm::vec3> rays;
	
	while (rays.size() < size) {
		glm::vec3 new_R = pur(R, exponent);
		if (glm::dot(normal, new_R) <= 0) continue;
        rays.push_back(new_R);
	}

	return rays;
}

glm::vec3 trace_colour (const std::list<Light *> &lights, 
			const glm::vec3 & ambient, 
			SceneNode * root,
                        Ray &r, 
			const int depth = 2,
			const int in_side = 0) {
    glm::vec3 phongLight = glm::vec3(0.0f);
    glm::vec3 reflectedColor;
    glm::vec3 refractedColor;
    glm::vec3 motionColor = glm::vec3(0.0f);
    glm::vec3 final_colour = glm::vec3(0.0f);

    Intersection intersection = root->Intersect(r);
    if (intersection.hit) {
	    phongLight = rayTracing(lights, ambient, root, r, intersection);
	    final_colour += phongLight;

	if (depth > 0) {
#ifdef IF_REFLECTION
	    Ray reflected_r = getReflect(intersection.normal, intersection.pos, r.direction);
	    reflectedColor = trace_colour(lights, ambient, root, reflected_r, depth-1);
  #ifdef IF_GLOSSY_REFLECTION
		float glossy_coeff = 0.5f;
		std::vector<glm::vec3> purtubed_reflected_rays
		  = purturbedRays(reflected_r.direction, intersection.normal, 4.0, 10);

		for (glm::vec3 dir : purtubed_reflected_rays) {
		  Ray glossy_r = Ray(intersection.pos+0.01f*dir, dir);
		  reflectedColor += 0.1 * trace_colour(lights, ambient, root, glossy_r, 0, in_side);
		}
  #endif
#endif		

#ifdef IF_REFRACTION
	    	double n_i = 1.52f;
	   	double n_t = 1.0f;

	    	if (in_side % 2 == 1) {
			std::swap(n_i, n_t);
	    	}
		
	    	Ray refracted_r = getRefract(intersection.normal, intersection.pos, r.direction, n_i, n_t);
	    
	    	if (refracted_r.direction == glm::vec3(0, 0, 0)) {
		     return final_colour;
	    	}
		
	   	refractedColor = trace_colour(lights, ambient, root, refracted_r, depth-1, in_side+1);
		
  #ifdef IF_GLOSSY_REFRACTION
		std::vector<glm::vec3> purtubed_refracted_rays
		 = purturbedRays(refracted_r.direction, intersection.normal, 4.0, 10);
		for (glm::vec3 dir : purtubed_refracted_rays) {
		    Ray glossy_r = Ray(intersection.pos+0.01f*dir, dir);
		    refractedColor = += 0.1 * trace_colour(lights, ambient, root, glossy_r, 0, in_side);
		}
  #endif
		double refract_coef = ((PhongMaterial *)intersection.colour)->m_refract_coef;
		double reflect_coef = (1 - refract_coef) / 4;
	    	final_colour = (1-reflect_coef-refract_coef)* final_colour + 
				reflect_coef * reflectedColor + refract_coef * refractedColor;
#endif

	}
    } else {
		glm::vec3 colour_d = glm::normalize(r.direction);
		phongLight = glm::vec3(0.0f, 0.5f, 0.5f) + colour_d.y * glm::vec3(1.0f, 1.0f, 0.0f);
		final_colour += phongLight;
    }
    return final_colour;
}

#ifdef IF_TEXTURE
glm::vec3 rayTracing (const std::list<Light *> &lights, const glm::vec3 & ambient, SceneNode * root,
			Ray &r, Intersection &intersection) {
    glm::vec3 kd = intersection.texture->getColour(intersection.u, intersection.v);

	glm::vec3 ks = glm::vec3(0.5, 0.7, 0.5);
	double shininess = 20;
    glm::vec3 phongLight = glm::vec3(0.0f);

	glm::vec3 v = glm::normalize(-r.direction);
    phongLight = ambient * kd;
    glm::vec3 intersection_pos = intersection.pos+0.01*intersection.normal;
	glm::vec3 inter_normal = intersection.normal;
        
    for (Light* light : lights) {
        Ray shadow = Ray(intersection_pos, light->position-intersection_pos);
	    glm::vec3 l = glm::normalize(shadow.direction);
		glm::vec3 diffuse = glm::vec3(0.0f);
		glm::vec3 specular = glm::vec3(0.0f);
	     				
		Intersection s_inter = root->Intersect(shadow);
		if (!s_inter.hit) {
	        	float nl = max((float)glm::dot(inter_normal, l), 0.0f);
			diffuse = kd * nl;
			if (nl > 0.0) {
		        glm::vec3 s = glm::normalize(v+l);
				float ns = max((float)glm::dot(inter_normal, s), 0.0f);
				
				specular = ks * std::pow(ns, shininess);
			}
        }
				
		phongLight += light->colour * (diffuse + specular);
	}

    return phongLight;
}

#else 
glm::vec3 rayTracing (const std::list<Light *> &lights, const glm::vec3 & ambient, SceneNode * root,
			Ray &r, Intersection &intersection) {
    
    glm::vec3 phongLight = glm::vec3(0.0f);

    if (intersection.colour->m_matType == MaterialType::Phong) {
	if (PhongMaterial* colour = dynamic_cast<PhongMaterial*>(intersection.colour)) {
	    glm::vec3 kd = colour->m_kd;
	    glm::vec3 ks = colour->m_ks;
	    double shininess = colour->m_shininess;
	    glm::vec3 v = glm::normalize(-r.direction);

	    phongLight = ambient * kd;
	    glm::vec3 intersection_pos = intersection.pos+0.01*intersection.normal;
	    glm::vec3 inter_normal = intersection.normal;

	    for (Light* light : lights) {
		glm::vec3 diffuse = glm::vec3(0.0f);
		glm::vec3 specular = glm::vec3(0.0f);

#ifdef IF_SOFT_SHADOW
		double shadow_blur = 40;
		for (int i=0; i < shadow_blur; i++) {
		    std::random_device generator;
		    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
		    double p1 = distribution(generator);
		    double p2 = distribution(generator);
		    glm::vec3 area_light = light->position + glm::vec3(p1*20, p2*20, 0);
#ifdef IF_MOTION_BLUR
            Ray shadow = Ray(intersection_pos, glm::normalize(area_light-intersection_pos), r.time);
#else
		    Ray shadow = Ray(intersection_pos, glm::normalize(area_light-intersection_pos));
#endif
		    double drop = length(shadow.direction);
		    double attenuation = 1.0 / (light->falloff[0] + light->falloff[1] * drop 
				+ light->falloff[2] * drop * drop);
		    glm::vec3 l = glm::normalize(shadow.direction);
		    
		    Intersection s_inter = root->Intersect(shadow);
		    if (!s_inter.hit) {
			    float nl = max((float)glm::dot(inter_normal, l), 0.0f);
			    diffuse += kd * nl * attenuation;
			    if (nl > 0.0) {
				glm::vec3 s = glm::normalize(v+l);
				float ns = max((float)glm::dot(inter_normal, s), 0.0f);

				specular += ks * std::pow(ns, shininess) * attenuation;
			    }
		    }
		}
		
		phongLight += light->colour * (diffuse + specular)/shadow_blur;
#else 

#ifdef IF_MOTION_BLUR
		Ray shadow = Ray(intersection_pos, glm::normalize(light->position-intersection_pos), r.time);
#else
		Ray shadow = Ray(intersection_pos, glm::normalize(light->position-intersection_pos));
#endif

				double drop = length(shadow.direction);     
                double attenuation = 1.0 / (light->falloff[0] + light->falloff[1] * drop
                                + light->falloff[2] * drop * drop);     
                glm::vec3 l = glm::normalize(shadow.direction);     
                    
                Intersection s_inter = root->Intersect(shadow);     
                if (!s_inter.hit) {
		    float nl = max((float)glm::dot(inter_normal, l), 0.0f);
		    diffuse += kd * nl * attenuation;
		    if (nl > 0.0) {
			glm::vec3 s = glm::normalize(v+l);
                        float ns = max((float)glm::dot(inter_normal, s), 0.0f);
			specular += ks * std::pow(ns, shininess) * attenuation;
		    }
		}
		phongLight += light->colour * (diffuse + specular);
#endif
	    }
	} else {
	    cout << "ERROR: Not Phong Material" << endl;  
	}
    } 

    return phongLight;
}

#endif

Ray getReflect (glm::vec3 normal, glm::vec3 position, glm::vec3 direction) {
    glm::vec3 ref_dir = glm::normalize(direction - 2 * normal * glm::dot(direction, normal));
    Ray r = Ray(position+0.01*normal, ref_dir);

    return r;
}

Ray getRefract (glm::vec3 normal, glm::vec3 position, glm::vec3 direction, double n_i, double n_t) {
    double nr = n_i/n_t;
    double c1 = glm::dot(normal, glm::normalize(direction));
    double s1 = 1 - c1 * c1;
    double c2 = 1 - nr * nr * s1;

    if (c2 <= 0) {
// std::cout << "generate fails" << endl;
	    return Ray(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0));
    }

    glm::vec3 dir = glm::normalize(nr * direction + (nr*c1-sqrt(c2)) * normal); 
    Ray r = Ray(position-0.01*normal, dir);
    return r;
}


void A4_Render(
		// What to render  
		SceneNode * root,

		// Image to write to, set to a given width and height  
		Image & image,

		// Viewing parameters  
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters  
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
) {
  std::cout << "Calling A4_Render(\n" <<
		  "\t" << *root <<
          "\t" << "Image(width:" << image.width() << ", height:" << image.height() << ")\n"
          "\t" << "eye:  " << glm::to_string(eye) << std::endl <<
		  "\t" << "view: " << glm::to_string(view) << std::endl <<
		  "\t" << "up:   " << glm::to_string(up) << std::endl <<
		  "\t" << "fovy: " << fovy << std::endl <<
          "\t" << "ambient: " << glm::to_string(ambient) << std::endl <<
		  "\t" << "lights{" << std::endl;

	for(const Light * light : lights) {
		std::cout << "\t\t" <<  *light << std::endl;
	}
	std::cout << "\t}" << std::endl;
	std:: cout <<")" << std::endl;

	size_t h = image.height();
	size_t w = image.width();

	glm::mat4 DCtoW = DCtoWorld(fovy, w, h, eye, view, up);

#ifdef IF_ANTI_ALIASING
	glm::vec3 final_colour_left_pixel = glm::vec3(0.0f);
	std::vector<glm::vec3> last_line_colours;
	while (last_line_colours.size() < w) {
	    last_line_colours.push_back(glm::vec3(0.0, 0.0, 0.0));
	}
#endif

	for (uint y = 0; y < h; ++y) {
	    for (uint x = 0; x < w; ++x) {
		glm::vec3 final_colour = glm::vec3(0.0f);

		glm::vec3 pixel = glm::vec3(DCtoW * glm::vec4(x, y, 1, 1));
	    glm::vec3 origin = eye;

#ifdef IF_MOTION_BLUR
        for (double i=1; i<=11; i+=1) {
            // double time = sin(i/10*(M_PI/2));
            if (i == 11) i=10;
            double time = sqrt(i/10);

            Ray r = Ray(origin, glm::normalize(pixel-origin), time);
#else 
		Ray r = Ray(origin, glm::normalize(pixel-origin));
#endif

        glm::vec3 direction = pixel-origin;
        // Ray r = Ray(origin, glm::normalize(pixel-origin));

#ifdef IF_DEPTH_OF_FIELD
	        float focal_plane = -500.0f; 
		for (int i = 0; i < 41; i++) {
		  std::random_device generator;
		  std::uniform_real_distribution<double> distribution(0.0,1.0);
		  double p1 = distribution(generator);
		  double p2 = distribution(generator);
		  glm::vec3 move = glm::vec3(p1*10, p2*10, 0);
		  glm::vec3 eye_pos = eye + move;
		  float ratio = (direction.z-focal_plane) / direction.z;
		  glm::vec3 focal_dir = direction * ratio;
		  focal_dir = focal_dir - move;

		  Ray new_r = Ray(eye_pos, glm::vec3(focal_dir));
		  final_colour += 0.025 * trace_colour(lights, ambient, root, new_r);
		}
#else 	
    #ifdef IF_MOTION_BLUR
        final_colour += trace_colour(lights, ambient, root, r, 0);
    #else
		final_colour += trace_colour(lights, ambient, root, r);
    #endif
#endif

#ifdef IF_ANTI_ALIASING
		glm::vec3 diff1 = final_colour_left_pixel - final_colour;
		glm::vec3 diff2 = last_line_colours.at(x) - final_colour;

		glm::vec3 aa_final_colour = glm::vec3(0.0);
		uint aa = 10;
		if (glm::length(diff1) > 0.003 || glm::length(diff2) > 0.003) {		    
		    for (double i=0; i<aa; i+=1) {
			for (double j=0; j<aa; j+=1) {
			    glm::vec3 new_pixel = glm::vec3(DCtoW * glm::vec4(x+i/aa, y+i/aa, 1, 1));
			    r = Ray(origin, glm::normalize(new_pixel-origin));  
			    aa_final_colour += (trace_colour(lights, ambient, root, r) / std::pow(aa, 2));
			}
		    }

		    final_colour = 0.6*aa_final_colour + 0.2*(final_colour_left_pixel + last_line_colours.at(x));
		}

		final_colour_left_pixel = final_colour;
		last_line_colours.at(x) = final_colour;
#endif

#ifdef IF_MOTION_BLUR
        }
        final_colour = final_colour/11;
#endif

        image(x, y, 0) = final_colour.r;
	    image(x, y, 1) = final_colour.g;
  		image(x, y, 2) = final_colour.b;
	    }
	}
}
