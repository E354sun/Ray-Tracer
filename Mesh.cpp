// Winter 2019

#include <iostream>
#include <fstream>
#include <algorithm>

#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <cmath>

#include "cs488-framework/ObjFileDecoder.hpp"
#include "Mesh.hpp"
#include "OPT_A5.hpp"

using namespace std;

Mesh::Mesh( const std::string& fname )
	: m_vertices()
	, m_faces()
	, boundingVolume(NULL)
{
	std::string code;
	double vx, vy, vz;
	size_t s1, s2, s3;

	std::ifstream ifs( fname.c_str() );
	while( ifs >> code ) {
		if( code == "v" ) {
			ifs >> vx >> vy >> vz;
			m_vertices.push_back( glm::vec3( vx, vy, vz ) );
		} else if( code == "f" ) {
			ifs >> s1 >> s2 >> s3;
			m_faces.push_back( Triangle( s1 - 1, s2 - 1, s3 - 1 ) );
		}
	}

  	float min_x = m_vertices.at(0).x;
  	float min_y = m_vertices.at(0).y;
  	float min_z = m_vertices.at(0).z;

  	float max_x = m_vertices.at(0).x;
  	float max_y = m_vertices.at(0).y;
  	float max_z = m_vertices.at(0).z;

  	for (const glm::vec3& v: m_vertices) {
    	min_x = std::min(min_x, v.x);
    	min_y = std::min(min_y, v.y);
    	min_z = std::min(min_z, v.z);

    	max_x = std::max(max_x, v.x);
    	max_y = std::max(max_y, v.y);
    	max_z = std::max(max_z, v.z);
  	}

  	glm::vec3 center = glm::vec3(min_x+max_x, min_y+max_y, min_z+max_z) / 2;
  	double radius = glm::length(glm::vec3(max_x-min_x, max_y-min_y, max_z-min_z)) / 2;

#ifdef IF_MOTION_BLUR
    boundingVolume = new NonhierSphere(center, center, radius);
#else
  	boundingVolume = new NonhierSphere(center, radius);
#endif
}

Intersection Mesh::Intersect(Ray &r) {
	Ray r_copy = Ray(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0));
	r_copy.origin = r.origin;
	r_copy.direction = r.direction;
	Intersection i;
	
	i = boundingVolume->Intersect(r_copy);
	
	if (!i.hit) {
		return i;
	}
	

	double t = std::numeric_limits<double>::infinity();
	glm::vec3 n = glm::vec3(0.0f);
	glm::vec3 direction = glm::vec3(r.direction);
	glm::vec3 origin = glm::vec3(r.origin);

	for (const Triangle &triangle : m_faces) {
    		const glm::vec3 a = m_vertices.at(triangle.v1);
    		const glm::vec3 b = m_vertices.at(triangle.v2);
    		const glm::vec3 c = m_vertices.at(triangle.v3);

		glm::mat3 A;
		A[0] = a-b;
		A[1] = a-c;
		A[2] = direction;

		if (glm::determinant(A) < 0.001f &&
		    glm::determinant(A) > 0.001f ) {
		    continue;
		}

		glm::vec3 T = glm::inverse(A) * (a-origin);

		if (T.x > 0 && T.y > 0 && T.x+T.y < 1) {
		    if (T.z < t && T.z > 0) {
			t = T.z;
		    }
		}
	
		n = glm::normalize(glm::cross(b-a, c-a));
		if (glm::dot(n, origin+direction*t-a) < 0.01f &&
		    glm::dot(n, origin+direction*t-a) > -0.01f) {
		    i.normal = n;
		}
	}
	
	if (t == std::numeric_limits<double>::infinity()) {
	    i.hit = false;
		return i;
	}

	i.t = t;
	i.hit = true;
	i.pos = r.origin + t*r.direction;

	return i;
}


std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh {";
  /*
  
  for( size_t idx = 0; idx < mesh.m_verts.size(); ++idx ) {
  	const MeshVertex& v = mesh.m_verts[idx];
  	out << glm::to_string( v.m_position );
	if( mesh.m_have_norm ) {
  	  out << " / " << glm::to_string( v.m_normal );
	}
	if( mesh.m_have_uv ) {
  	  out << " / " << glm::to_string( v.m_uv );
	}
  }

*/
  out << "}";
  return out;
}
