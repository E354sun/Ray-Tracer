// Fall 2018

#pragma once

#include <glm/glm.hpp>
#include "OPT_A5.hpp"
#include "Ray.hpp"

class Primitive {
public:
  virtual Intersection Intersect(Ray &r); 
  virtual ~Primitive();
};

#ifdef IF_MOTION_BLUR
class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos1, const glm::vec3& pos2, double radius)
    : m_pos1(pos1), m_pos2(pos2), m_radius(radius)
  {
  }

  virtual ~NonhierSphere();

  Intersection Intersect(Ray &r) override;

private:
  glm::vec3 m_pos1;
  glm::vec3 m_pos2;

  double m_radius;
};

#else 
class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }

  virtual ~NonhierSphere();

  Intersection Intersect(Ray &r) override;

private:
  glm::vec3 m_pos;

  double m_radius;
};
#endif

class NonhierBox : public Primitive {
public:
  NonhierBox(const glm::vec3& pos, double size)
    : m_pos(pos), m_size(size)
  {
  }

  Intersection Intersect(Ray &r) override;
  virtual ~NonhierBox();

private:
  glm::vec3 m_pos;
  double m_size;
};

class Sphere : public Primitive {
  NonhierSphere *m_sphere;
public:
  virtual Intersection Intersect(Ray &r);
  virtual ~Sphere();
  Sphere();
};

class Cube : public Primitive {
  NonhierBox *m_box;
public:
  virtual Intersection Intersect(Ray &r);
  virtual ~Cube();
  Cube();
};

class Cone : public Primitive {

    public:
    Cone();
    virtual ~Cone();
    Intersection Intersect(Ray &r);
};


class Cylinder : public Primitive {

    public:
    virtual ~Cylinder();
    Intersection Intersect(Ray &r);
    Cylinder();
};

/*
class NonhierCone : public Primitive {
public:
  NonhierCone(const glm::vec3& pos, double radius, double height)
    : m_pos(pos), m_radius(radius), m_height(height)
  {
  }

  Intersection Intersect(Ray &r) override;
  virtual ~NonhierCone();

private:
  glm::vec3 m_pos;
  double m_radius;
  double m_height;
};

class NonhierCylinder : public Primitive {
public:
  NonhierCylinder(const glm::vec3& pos, double radius, double height)
    : m_pos(pos), m_radius(radius), m_height(height)
  {
  }

  Intersection Intersect(Ray &r) override;
  virtual ~NonhierCylinder();

private:
  glm::vec3 m_pos;
  double m_radius;
  double m_height;
};
*/
