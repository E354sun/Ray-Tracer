#include "glassMaterial.hpp"
#include <glm/ext.hpp>
#include <iostream>

using namespace std;

glassMaterial::glassMaterial() : Material() {}

glassMaterial::glassMaterial(double transmittance, double reflectivity, double refractive_idx)
 : Material(refractive_idx), transmittance(transmittance), reflectivity(reflectivity) 
  {
      m_matType = MaterialType::Glass;
  }

glassMaterial::~glassMaterial()
{}

/*
glm::vec3 glassMaterial::getColor(glm::vec3 pHit, glm::vec3 pNormal,
				  Light *light, glm::mat4 inv, Material *lastMat) {

  return vec3(0);
}
*/
