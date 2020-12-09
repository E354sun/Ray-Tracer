#pragma once

#include <glm/glm.hpp>
#include "PhongMaterial.hpp"
#include <algorithm>
#include <iostream>
#include "OPT_A5.hpp"
#include "texture.hpp"

class Ray {				
public:
    glm::vec3 origin;
    glm::vec3 direction;
    double time;// dest - origin 
    Ray(const glm::vec3 &o, const glm::vec3 &d) 
    	: origin(o), direction(d) {
            time = 0.0;        
        };
    Ray(const glm::vec3 &o, const glm::vec3 &d, double t)
        : origin(o), direction(d), time(t) {};
};

class Intersection {
public:
    glm::vec3 pos;
    glm::vec3 normal;
    float t;
    bool hit;

#ifdef IF_TEXTURE
    Texture* texture;
    float u;
    float v;
#else 
    Material* colour;
#endif

    // default constructor
    Intersection() {
        pos = glm::vec3(0.0f);
	    normal = glm::vec3(0.0f);
	    t = std::numeric_limits<float>::infinity();
	    hit = false;
    }
};
