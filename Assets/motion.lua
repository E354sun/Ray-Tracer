-- A simple scene with some miscellaneous geometry.
-- This file is very similar to nonhier.lua, but interposes
-- an additional transformation on the root node.  
-- The translation moves the scene, and the position of the camera
-- and lights have been modified accordingly.

mat1 = gr.material({0.7, 1.0, 0.7}, {0.5, 0.7, 0.5}, 25, 0)
mat2 = gr.material({0.5, 0.5, 0.5}, {0.5, 0.7, 0.5}, 25, 0)
mat3 = gr.material({1.0, 0.6, 0.1}, {0.5, 0.7, 0.5}, 25, 0)
mat4 = gr.material({0.7, 0.6, 1.0}, {0.5, 0.4, 0.8}, 25, 0)
white = gr.material({0.7, 0.7, 0.7}, {0.5, 0.4, 0.8}, 25, 0)

scene = gr.node( 'scene' )
-- scene:translate(0, 0, -800)
--[[
s1 = gr.nh_sphere('s1', {50, -50, 600}, 20)
scene:add_child(s1)
s1:set_material(mat1)

s2 = gr.nh_sphere('s2', {0, -25, 450}, 20)
scene:add_child(s2)
s2:set_material(mat1)
]]--

s3 = gr.nh_sphere('s3', {-50, 25, 300}, {-50, -25, 300}, 35)
scene:add_child(s3)
s3:set_material(mat2)

s4 = gr.nh_sphere('s3', {50, 0, 300}, {50, 0, 300}, 35)
scene:add_child(s4)
s4:set_material(mat2)

plane = gr.cube('plane')
plane:scale(1500, 200, 1000)
plane:translate(-600, -250, 200)
scene:add_child(plane)
plane:set_material(white)

-- b1 = gr.nh_box('b1', {-200, -125, 0}, 175)
-- scene:add_child(b1)
-- b1:set_material(mat4)
--
--[[
s4 = gr.nh_sphere('s4', {-100, 25, 150}, 20)
scene:add_child(s4)
s4:set_material(mat3)

s5 = gr.nh_sphere('s5', {-150, 40, 10}, 20)
scene:add_child(s5)
s5:set_material(mat1)
]]--
-- A small stellated dodecahedron.

-- steldodec = gr.mesh( 'dodec', 'smstdodeca.obj' )
-- steldodec:set_material(mat3)
-- scene:add_child(steldodec)

white_light = gr.light({-400.0, 300.0, 600.0}, {0.7, 0.7, 0.7}, {1, 0, 0})
yellow_light = gr.light({150.0, 300.0, -600.0}, {0.7, 0.7, 0.5}, {1, 0, 0})
-- magenta_light = gr.light({400.0, 100.0, -650.0}, {0.7, 0.0, 0.7}, {1, 0, 0})

gr.render(scene, 'motion.png', 512, 512,
	  {0, 0, 800}, {0, 0, -1}, {0, 1, 0}, 50,
	  {0.3, 0.3, 0.3}, {white_light, yellow_light})
