// Winter 2019

#pragma once

#include "SceneNode.hpp"
#include "Primitive.hpp"
#include "Material.hpp"
#include "texture.hpp"
#include "OPT_A5.hpp"

#ifdef IF_TEXTURE
class GeometryNode : public SceneNode {
public:
	GeometryNode( const std::string & name, Primitive *prim, 
		Texture *tex = nullptr );

	void setMaterial( Texture *tex );

	Intersection Intersect(Ray &r);

	Texture *m_texture;
    // Material *m_material;
	Primitive *m_primitive;
};

#else 
class GeometryNode : public SceneNode {
public:
	GeometryNode( const std::string & name, Primitive *prim, 
		Material *mat = nullptr );

	void setMaterial( Material *material );

	Intersection Intersect(Ray &r);

	Material *m_material;
	Primitive *m_primitive;
};
#endif
