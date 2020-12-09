// Winter 2019

#include "PhongMaterial.hpp"

PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess, double refract_coef )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
	, m_refract_coef(refract_coef)
{
    m_matType = MaterialType::Phong;
}
/*
PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess, double refractIdx )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
	, m_refractIdx(refractIdx)
{}
*/
PhongMaterial::~PhongMaterial()
{}
