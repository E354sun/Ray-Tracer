#pragma once

#include <glm/glm.hpp>
#include "Material.hpp"

class glassMaterial : public Material {
public:
  glassMaterial();
  glassMaterial(double transmittance, double reflectivity, double refractive_idx);

  virtual ~glassMaterial();
  
  double reflectivity;
  double transmittance;
  // double refrective_idx;

  /*
  glm::vec3 getColor(glm::vec3 pHit, 
    glm::vec3 pNormal, 
    Light *light,
    glm::mat4 inv,
    Material *lastMat);
  */
};
