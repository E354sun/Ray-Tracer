// Winter 2019

#pragma once

enum class MaterialType {
    Phong,
    Glass
};

class Material {
public:
  virtual ~Material();
  Material();
  Material(double refract_idx);

  // virtual glm::vec3 getColor(glm::vec3 pHit, glm::vec3 pNormal, Light *light, glm::mat4 inv, Material *lastMat);
  double refractive_idx;
  MaterialType m_matType;
// protected:
// Material();
};
