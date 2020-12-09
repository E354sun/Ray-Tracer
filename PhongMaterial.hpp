// Winter 2019

#pragma once

#include <glm/glm.hpp>

#include "Material.hpp"

class PhongMaterial : public Material {
public:
  // PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess);
  PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess, double refract_coef); 

  virtual ~PhongMaterial();

  glm::vec3 m_kd;
  glm::vec3 m_ks;
  double m_shininess;

  double m_refract_coef;
};
