// Winter 2019
#pragma once

#include <glm/glm.hpp>

#include "SceneNode.hpp"
#include "Light.hpp"
#include "Image.hpp"

glm::mat4 DCtoWorld (
        const double fovy,
        const double width,
        const double height,
        const glm::vec3 &eye,
        const glm::vec3 &view,
        const glm::vec3 &up);

Ray getRefract (glm::vec3 normal, glm::vec3 position, glm::vec3 direction, double n_i, double n_t);

Ray getReflect (glm::vec3 normal, glm::vec3 position, glm::vec3 direction);

glm::vec3 rayTracing (const std::list<Light *> &lights, const glm::vec3 & ambient, SceneNode * root,
                        Ray &r, Intersection &intersection);

void A4_Render(
		// What to render
		SceneNode * root,

		// Image to write to, set to a given width and height
		Image & image,

		// Viewing parameters
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
);

/*
void trace_colour (const std::list<Light *> &lights,
const glm::vec3 & ambient,
SceneNode * root,
Ray &r,
glm::vec3 &final_colour);
// const int depth = 2,
// const int in_side = 0);
*/
