// Winter 2019

#include "SceneNode.hpp"

#include "cs488-framework/MathUtils.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
using namespace std;

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/transform.hpp>

// using namespace glm;


// Static class variable
unsigned int SceneNode::nodeInstanceCount = 0;


//---------------------------------------------------------------------------------------
SceneNode::SceneNode(const std::string& name)
  : m_name(name),
	m_nodeType(NodeType::SceneNode),
	trans(glm::mat4()),
	invtrans(glm::mat4()),
	m_nodeId(nodeInstanceCount++)
{

}

//---------------------------------------------------------------------------------------
// Deep copy
SceneNode::SceneNode(const SceneNode & other)
	: m_nodeType(other.m_nodeType),
	  m_name(other.m_name),
	  trans(other.trans),
	  invtrans(other.invtrans)
{
	for(SceneNode * child : other.children) {
		this->children.push_front(new SceneNode(*child));
	}
}

//---------------------------------------------------------------------------------------
SceneNode::~SceneNode() {
	for(SceneNode * child : children) {
		delete child;
	}
}

//---------------------------------------------------------------------------------------
void SceneNode::set_transform(const glm::mat4& m) {
	trans = m;
	invtrans = glm::inverse(m);
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_transform() const {
	return trans;
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_inverse() const {
	return invtrans;
}

//---------------------------------------------------------------------------------------
void SceneNode::add_child(SceneNode* child) {
	children.push_back(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::remove_child(SceneNode* child) {
	children.remove(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::rotate(char axis, float angle) {
	glm::vec3 rot_axis;

	switch (axis) {
		case 'x':
			rot_axis = glm::vec3(1,0,0);
			break;
		case 'y':
			rot_axis = glm::vec3(0,1,0);
	        break;
		case 'z':
			rot_axis = glm::vec3(0,0,1);
	        break;
		default:
			break;
	}
	glm::mat4 rot_matrix = glm::rotate(degreesToRadians(angle), rot_axis);
	set_transform( rot_matrix * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::scale(const glm::vec3 & amount) {
	set_transform( glm::scale(amount) * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::translate(const glm::vec3& amount) {
	set_transform( glm::translate(amount) * trans );
}


//---------------------------------------------------------------------------------------
int SceneNode::totalSceneNodes() const {
	return nodeInstanceCount;
}

Intersection SceneNode::Intersect(Ray &r){
	Intersection i;
	
	r.origin = glm::vec3(invtrans * glm::vec4(r.origin, 1));
	r.direction = glm::vec3(invtrans * glm::vec4(r.direction, 0));

	for (SceneNode * child : children) {
		Intersection temp = child->Intersect(r);
		if (temp.hit && temp.t < i.t) {
		    i = temp;
		}
	}
	
    	r.origin = glm::vec3(trans * glm::vec4(r.origin, 1));
	r.direction = glm::vec3(trans * glm::vec4(r.direction, 0));
	
	i.pos = glm::vec3(trans * glm::vec4(i.pos, 1));
	i.normal = glm::normalize(glm::vec3(glm::transpose(invtrans) * glm::vec4(i.normal, 1)));
	 
	return i;
}

//---------------------------------------------------------------------------------------
std::ostream & operator << (std::ostream & os, const SceneNode & node) {

	//os << "SceneNode:[NodeType: ___, name: ____, id: ____, isSelected: ____, transform: ____"
	switch (node.m_nodeType) {
		case NodeType::SceneNode:
			os << "SceneNode";
			break;
		case NodeType::GeometryNode:
			os << "GeometryNode";
			break;
		case NodeType::JointNode:
			os << "JointNode";
			break;
	}
	os << ":[";

	os << "name:" << node.m_name << ", ";
	os << "id:" << node.m_nodeId;

	os << "]\n";
	return os;
}
