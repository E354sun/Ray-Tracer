// Winter 2019

#include "GeometryNode.hpp"
#include <glm/glm.hpp>
#include <algorithm>
#include <glm/gtx/string_cast.hpp>

using namespace std;

//---------------------------------------------------------------------------------------
#ifdef IF_TEXTURE
GeometryNode::GeometryNode(
	const std::string & name, Primitive *prim, Texture *tex )
	: SceneNode( name )
	, m_texture( tex )
	, m_primitive( prim )
{
	m_nodeType = NodeType::GeometryNode;
}

void GeometryNode::setMaterial(Texture *tex )
{
	// Obviously, there's a potential memory leak here.  A good solution
	// would be to use some kind of reference counting, as in the 
	// C++ shared_ptr.  But I'm going to punt on that problem here.
	// Why?  Two reasons:
	// (a) In practice we expect the scene to be constructed exactly
	//     once.  There's no reason to believe that materials will be
	//     repeatedly overwritten in a GeometryNode.
	// (b) A ray tracer is a program in which you compute once, and 
	//     throw away all your data.  A memory leak won't build up and
	//     crash the program.

	m_texture = tex;
}

#else 
GeometryNode::GeometryNode(
	const std::string & name, Primitive *prim, Material *mat )
	: SceneNode( name )
	, m_material( mat )
	, m_primitive( prim )
{
	m_nodeType = NodeType::GeometryNode;
}

void GeometryNode::setMaterial( Material *mat )
{
	// Obviously, there's a potential memory leak here.  A good solution
	// would be to use some kind of reference counting, as in the 
	// C++ shared_ptr.  But I'm going to punt on that problem here.
	// Why?  Two reasons:
	// (a) In practice we expect the scene to be constructed exactly
	//     once.  There's no reason to believe that materials will be
	//     repeatedly overwritten in a GeometryNode.
	// (b) A ray tracer is a program in which you compute once, and 
	//     throw away all your data.  A memory leak won't build up and
	//     crash the program.

	m_material = mat;
}
#endif


Intersection GeometryNode::Intersect(Ray &r) {
	Intersection i;

    r.origin = glm::vec3(invtrans * glm::vec4(r.origin, 1));
	r.direction = glm::vec3(invtrans * glm::vec4(r.direction, 0));
	i = m_primitive->Intersect(r);

	if (i.hit) {
#ifdef IF_TEXTURE
        i.texture = m_texture;
#else 
        i.colour = m_material;
#endif

	    i.pos = glm::vec3(trans * glm::vec4(i.pos, 1));
	    i.normal = glm::normalize(glm::vec3(glm::transpose(invtrans) * glm::vec4(i.normal, 0)));
	}
	
	r.direction = glm::vec3(trans * glm::vec4(r.direction, 0));
	r.origin = glm::vec3(trans * glm::vec4(r.origin, 1));

	return i;
}
