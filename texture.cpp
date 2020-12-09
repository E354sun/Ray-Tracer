#include "texture.hpp"
#include <lodepng/lodepng.h>
#include <iostream>

Texture::Texture (std::string filename)
{
    unsigned error = lodepng::decode(image, width, height, filename);
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

glm::vec3 Texture::getColour(float u, float v)
{    
    int i = u * width;
    int j = (1-v) * height -0.001;

    i = std::max(i, 0);
    j = std::max(j, 0);
    i = std::min((unsigned int)i, width-1);
    j = std::min((unsigned int)j, height-1);

    int pos = (i+j*width) * 4;
    // int pos = (i+j*width);

    return glm::vec3(image[pos]/255.0, image[pos+1]/255.0, image[pos+2]/255.0);
              
}
