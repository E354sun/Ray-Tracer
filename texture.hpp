#pragma once
#include <glm/glm.hpp>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

class Texture {
public:
    Texture(std::string filename);
    glm::vec3 getColour(float u, float v);
                    
    std::vector<unsigned char> image;
    unsigned width;
    unsigned height;
};
