// Winter 2019
#include <algorithm>
#include "Primitive.hpp"
#include "polyroots.hpp"
#include "Ray.hpp"

#include <iostream>
#include <cmath>
#include <limits>
#include <cfloat>
#include "OPT_A5.hpp"


using namespace glm;

int hit_times_bottom = 0;
int hit_times_top = 0;

Primitive::~Primitive()
{
}

Intersection Primitive::Intersect(Ray &r) {
    Intersection i;
    return i;
}

Sphere::Sphere() {
#ifdef IF_MOTION_BLUR
    m_sphere = new NonhierSphere(glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 0.0), 1.0);
#else
    m_sphere = new NonhierSphere(glm::vec3(0.0, 0.0, 0.0), 1.0);
#endif
}


Sphere::~Sphere()
{
	delete m_sphere;
}

Intersection Sphere::Intersect(Ray &r) {
    Intersection hit = m_sphere->Intersect(r);
    return hit;
}

Cube::Cube() {
    m_box = new NonhierBox(glm::vec3(0.0, 0.0, 0.0), 1.0);
}

Cube::~Cube()
{
	delete m_box;
}

Intersection Cube::Intersect(Ray &r) {
    Intersection hit = m_box->Intersect(r);
    return hit;
}

NonhierSphere::~NonhierSphere()
{
}

#ifdef IF_MOTION_BLUR
Intersection NonhierSphere::Intersect(Ray &r) 
{
        // glm::vec3 m_pos = glm::mix(m_pos1, m_pos2, r.time);
        glm::vec3 m_pos = m_pos1 + (float)r.time * (m_pos2 - m_pos1);     

		double roots[2];
  		const double A = dot(r.direction, r.direction);
  		const double B = 2 * dot(r.direction, r.origin-m_pos);
  		const double C = dot(r.origin-m_pos, r.origin-m_pos) 
  							- m_radius * m_radius;

  		size_t numRoots = quadraticRoots(A, B, C, roots);
  		
  		Intersection i;

  		if (numRoots == 0 ||
  			(numRoots == 1 && roots[0] < 0) ||
  			(numRoots == 2 && roots[0] < 0 && roots[1] < 0)) {
  			return i;
  		} 
  		
  		else if (numRoots == 1) {		// tangent line
		    i.t = roots[0];
		    i.hit = true;
		    i.pos = r.origin + (r.direction*i.t);
		    i.normal = glm::normalize(i.pos-m_pos);
  		} 
  		else {
		    if (roots[0] < 0) {
			    i.t = roots[1];
			}
		    else if (roots[1] < 0) {
			    i.t = roots[0];
			}	    
		    else {
		    	i.t = std::min(roots[0], roots[1]);
		    }
		    
		    i.hit = true;
		    i.pos = r.origin + (r.direction*i.t);
		    i.normal = glm::normalize(i.pos-m_pos);
  		}
  		return i;
}

#else
Intersection NonhierSphere::Intersect(Ray &r) 
{
		double roots[2];
  		const double A = dot(r.direction, r.direction);
  		const double B = 2 * dot(r.direction, r.origin-m_pos);
  		const double C = dot(r.origin-m_pos, r.origin-m_pos) 
  							- m_radius * m_radius;

  		size_t numRoots = quadraticRoots(A, B, C, roots);
  		
  		Intersection i;

  		if (numRoots == 0 ||
  			(numRoots == 1 && roots[0] < 0) ||
  			(numRoots == 2 && roots[0] < 0 && roots[1] < 0)) {
  			return i;
  		} 
  		
  		else if (numRoots == 1) {		// tangent line
		    i.t = roots[0];
		    i.hit = true;
		    i.pos = r.origin + (r.direction*i.t);
		    i.normal = glm::normalize(i.pos-m_pos);
  		} 
  		else {
		    if (roots[0] < 0) {
			    i.t = roots[1];
			}
		    else if (roots[1] < 0) {
			    i.t = roots[0];
			}	    
		    else {
		    	i.t = std::min(roots[0], roots[1]);
		    }
		    
		    i.hit = true;
		    i.pos = r.origin + (r.direction*i.t);
		    i.normal = glm::normalize(i.pos-m_pos);
  		}

#ifdef IF_TEXTURE
        if (i.hit) {
            glm::vec3 pos = i.pos-m_pos;
            float phi = atan2(pos.z, pos.x);
            float theta = asin(pos.y);
            i.u = 1-(phi+M_PI) / (2*M_PI);
            i.v = (theta+M_PI/2) / M_PI; 
        }
#endif
  		return i;
}
#endif

NonhierBox::~NonhierBox()
{

}

Intersection NonhierBox::Intersect(Ray &r) 
{
		glm::vec3 direction = r.direction;
		glm::vec3 origin = r.origin;
		
		glm::vec3 n = glm::vec3(0.0f);
		glm::vec3 normal = glm::vec3(0.0f);

		Intersection i;
		// i.normal = glm::vec3(0.0f);

		float t = std::numeric_limits<float>::infinity();

		std::vector<glm::vec3> normals;
		normals.push_back(glm::vec3(0, 0, 1)); // 
		normals.push_back(glm::vec3(0, 1, 0));
		normals.push_back(glm::vec3(1, 0, 0));
		normals.push_back(glm::vec3(0, 0, -1));
		normals.push_back(glm::vec3(0, -1, 0));
		normals.push_back(glm::vec3(-1, 0, 0));

		std::vector<glm::vec3> corners;
		corners.push_back(glm::vec3(0, 0, 1)); //
		corners.push_back(glm::vec3(0, 1, 0));
		corners.push_back(glm::vec3(0, 0, 0));

		corners.push_back(glm::vec3(0, 0, 1)); //
		corners.push_back(glm::vec3(0, 1, 0));
		corners.push_back(glm::vec3(0, 1, 1));

		corners.push_back(glm::vec3(0, 0, 0));
		corners.push_back(glm::vec3(1, 0, 0));
		corners.push_back(glm::vec3(0, 0, 1));

		corners.push_back(glm::vec3(1, 0, 1));
		corners.push_back(glm::vec3(1, 0, 0));
		corners.push_back(glm::vec3(0, 0, 1));

		corners.push_back(glm::vec3(0, 1, 0));
		corners.push_back(glm::vec3(1, 1, 0));
		corners.push_back(glm::vec3(0, 1, 1));

		corners.push_back(glm::vec3(1, 1, 1));
		corners.push_back(glm::vec3(1, 1, 0));
		corners.push_back(glm::vec3(0, 1, 1));

		corners.push_back(glm::vec3(1, 0, 1)); //
		corners.push_back(glm::vec3(1, 1, 0));
		corners.push_back(glm::vec3(1, 0, 0));

		corners.push_back(glm::vec3(1, 0, 1)); //
		corners.push_back(glm::vec3(1, 1, 0));
		corners.push_back(glm::vec3(1, 1, 1));

		corners.push_back(glm::vec3(1, 1, 0)); //
		corners.push_back(glm::vec3(0, 1, 0));
		corners.push_back(glm::vec3(0, 0, 0));

		corners.push_back(glm::vec3(0, 0, 0)); //
		corners.push_back(glm::vec3(1, 0, 0));
		corners.push_back(glm::vec3(1, 1, 0));

		corners.push_back(glm::vec3(1, 1, 1)); //
		corners.push_back(glm::vec3(0, 1, 1));
		corners.push_back(glm::vec3(0, 0, 1));

		corners.push_back(glm::vec3(0, 0, 1)); //
		corners.push_back(glm::vec3(1, 0, 1));
		corners.push_back(glm::vec3(1, 1, 1));

		for (int i=0; i < 36; i ++) {
		    corners.at(i) *= m_size;
		    corners.at(i) += m_pos;
		}

		for (int i=0; i < 36; i += 3) {
                const glm::vec3 a = corners.at(i);
                const glm::vec3 b = corners.at(i+1);
                const glm::vec3 c = corners.at(i+2);

                glm::mat3 A;
                A[0] = a-b;
                A[1] = a-c;
                A[2] = direction;

                if (glm::determinant(A) < 0.001f &&
                    glm::determinant(A) > 0.001f ) {
                    continue;
                }

                glm::vec3 T = glm::inverse(A) * (a-origin);

                if (T.x > -0.00001 && T.y > -0.00001 && T.x+T.y < 1.00001) {
                    if (T.z < t && T.z > -0.00001) {
                        t = T.z;
                    }
                }

                n = glm::normalize(glm::cross(b-a, c-a));
                if (glm::dot(n, origin+(direction*t)-a) < 0.0001f &&
                    glm::dot(n, origin+(direction*t)-a) > -0.0001f) {
                    normal = n;
                }
        }

        if (std::isinf(t)) {
	    	i.hit = false;
            return i;
        }

        i.t = t;
		i.hit = true;
		i.normal = normal;
		i.pos = r.origin + t*r.direction;

#ifdef IF_TEXTURE
        // if (i.normal == glm::vec3(0, 0, -1)) {
            glm::vec3 p = i.pos - m_pos;
            i.u = p.x;
            i.v = p.y;

            // std::cout << "u, v: " << i.u << ", " << i.v << std::endl;
        // }
#endif

        return i;
}

Cylinder::Cylinder() {
}

Cylinder::~Cylinder()
{

}

Intersection Cylinder::Intersect(Ray &r) 
{
		double roots[2];
		glm::vec3 d = r.direction;
		glm::vec3 o = r.origin;
		
  		const double A = d.x*d.x + d.z*d.z;
  		const double B = 2 * (o.x*d.x + o.z*d.z);
  		const double C = o.x*o.x + o.z*o.z - 1;

  		size_t numRoots = quadraticRoots(A, B, C, roots);
  		
  		Intersection i;

  		if (numRoots == 0 ||
  			(numRoots == 1 && roots[0] < 0) ||
  			(numRoots == 2 && roots[0] < 0 && roots[1] < 0)) {
  			return i;
  		}
  		else if (numRoots == 1 || (numRoots == 2 && roots[1] < 0)) {
		    i.t = roots[0];	    
		    i.pos = r.origin + (r.direction*i.t);
		    
		    if (i.pos.y <= 1 && i.pos.y >= -1) {
			i.normal = glm::normalize(glm::vec3(i.pos.x, 0, i.pos.z));
			i.hit = true;
		    }
  		}
		else if (numRoots == 2 && roots[0] < 0) {
		    i.t = roots[1];
		    i.pos = r.origin + (r.direction*i.t);

		    if (i.pos.y <= 1 && i.pos.y >= -1) {
			i.normal = glm::normalize(glm::vec3(i.pos.x, 0, i.pos.z));
			 i.hit = true;
		    }
		}
  		else {
		    float t1 = std::min(roots[0], roots[1]);
		    float t2 = std::max(roots[0], roots[1]);

		    glm::vec3 p1 = r.origin + t1*r.direction;
		    glm::vec3 p2 = r.origin + t2*r.direction;

		    if (p1.y <= 1 && p1.y >= -1) {
			i.normal = glm::normalize(glm::vec3(p1.x, 0, p1.z));
			i.t = t1;
			i.pos = p1;
			i.hit = true;
		    }
		    else if (p2.y <= 1 && p2.y >= -1) {
			i.normal = glm::normalize(glm::vec3(p2.x, 0, p2.z));
			i.t = t2;
			i.pos = p2;
			i.hit = true;
		    }
  		}

		glm::vec3 top_n = glm::vec3(0, 1, 0);
		glm::vec3 top_origin = r.origin-glm::vec3(0, 1, 0);
		float top_t = -glm::dot(top_n, top_origin) / glm::dot(top_n, r.direction);
		glm::vec3 top_pos = r.origin + top_t * r.direction;

		glm::vec3 bottom_n = glm::vec3(0, -1, 0);
		glm::vec3 bottom_origin = r.origin-glm::vec3(0, -1, 0);
		float bottom_t = -glm::dot(bottom_n, bottom_origin) / glm::dot(bottom_n, r.direction);
		glm::vec3 bottom_pos = r.origin + bottom_t * r.direction;

		if (top_pos.x * top_pos.x + top_pos.z * top_pos.z <= 1) {
// hit_times_top ++;
// std::cout << "hit top circle: " << hit_times_top << std::endl;
		  if (top_t > 0 && top_t < i.t) {
		    i.hit = true;
		    i.t = top_t;
		    i.pos = top_pos;
		    i.normal = top_n;
// hit_times_top ++;
// std::cout << "hit top circle: " << hit_times_top << std::endl;
		  }
		}
		if (bottom_pos.x * bottom_pos.x + bottom_pos.z * bottom_pos.z <= 1) {
		  // if (bottom_t > 0 && bottom_t < i.t) {
// hit_times_bottom ++;
// std::cout << "hit bottom circle: " << hit_times_bottom << std::endl; 
		  if (bottom_t > 0 && bottom_t < i.t) {
		    i.hit = true;
		    i.t = bottom_t;
		    i.pos = bottom_pos;
		    i.normal = bottom_n;
		  }
		}

  		return i;
}

Cone::Cone() {
}

Cone::~Cone()
{

}


Intersection Cone::Intersect(Ray &r) 
{
    /* cone */
		double roots[2];
		
		glm::vec3 d = r.direction;
		glm::vec3 o = r.origin;
  		const double A = d.x*d.x + d.z*d.z - d.y*d.y;
  		const double B = 2 * (o.x*d.x + o.z*d.z - o.y*d.y + d.y);
  		const double C = o.x*o.x + o.z*o.z - o.y*o.y + 2*o.y - 1;
  		size_t numRoots = quadraticRoots(A, B, C, roots);
  		Intersection i;

  		if (numRoots == 0 ||
  			(numRoots == 1 && roots[0] < 0) ||
  			(numRoots == 2 && roots[0] < 0 && roots[1] < 0)) {
  			i.hit = false;
  		}
	       
  		else if (numRoots == 1 || (numRoots == 2 && roots[1] < 0)) {
		    i.t = roots[0];	    
		    i.pos = r.origin + (r.direction*i.t);
		    
		    if (i.pos.y <= 1 && i.pos.y >= 0) {
				double length = sqrt(i.pos.x*i.pos.x + i.pos.z*i.pos.z);
				glm::vec3 X = glm::vec3(1.0, i.pos.x/length, 0.0);
				glm::vec3 Z = glm::vec3(0.0, i.pos.z/length, 1.0);
				i.normal = glm::normalize(glm::cross(X, Z));
				i.hit = true;
		    }
  		}
	       	else if (numRoots == 2 && roots[0] < 0) {
		    i.t = roots[1];
		    i.pos = r.origin + (r.direction*i.t);

		    if (i.pos.y <= 1 && i.pos.y >= 0) {
			double length = sqrt(i.pos.x*i.pos.x + i.pos.z*i.pos.z);
			glm::vec3 X = glm::vec3(1.0, i.pos.x/length, 0.0);
			glm::vec3 Z = glm::vec3(0.0, i.pos.z/length, 1.0);
			i.normal = glm::normalize(glm::cross(X, Z));
			 i.hit = true;
		    }
		}
		else {
		    float t1 = std::min(roots[0], roots[1]);
		    float t2 = std::max(roots[0], roots[1]);
		    glm::vec3 p1 = r.origin + t1*r.direction;
		    glm::vec3 p2 = r.origin + t2*r.direction;

		    if (p1.y <= 1 && p1.y >= 0) {
			double length = sqrt(p1.x*p1.x + p1.z*p1.z);
			glm::vec3 X = glm::vec3(1.0, p1.x/length, 0.0);
			glm::vec3 Z = glm::vec3(0.0, p1.z/length, 1.0);
			i.normal = glm::normalize(glm::cross(X, Z));
			i.t = t1;
			i.pos = p1;
			i.hit = true;
		    }
		    else if (p2.y <= 1 && p2.y >= 0) {
			double length = sqrt(p2.x*p2.x + p2.z*p2.z);
			glm::vec3 X = glm::vec3(1.0, p2.x/length, 0.0);
			glm::vec3 Z = glm::vec3(0.0, p2.z/length, 1.0);
			i.normal = glm::normalize(glm::cross(X, Z));
			i.t = t2;
			i.pos = p2;
			i.hit = true;
		    } 
		}

   /* circle plane */
		glm::vec3 circle_n = glm::vec3(0, -1, 0);
		float circle_t = -glm::dot(circle_n, r.origin) / glm::dot(circle_n, r.direction);
		glm::vec3 circle_pos = r.origin + circle_t * r.direction;

		if (circle_pos.x * circle_pos.x + circle_pos.z * circle_pos.z <= 1) {
		  if (circle_t > 0 && circle_t < i.t) {
		    i.hit = true;
		    i.t = circle_t;
		    i.pos = circle_pos;
		    i.normal = circle_n;
		  }
		}

		return i;
}
