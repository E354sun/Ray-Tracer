-- A simple scene with five spheres

mat1 = gr.material({0.7, 1.0, 0.7}, {0.5, 0.7, 0.5}, 25, 0)
mat2 = gr.material({0.5, 0.5, 0.5}, {0.5, 0.7, 0.5}, 25, 0)
mat3 = gr.material({1.0, 0.6, 0.1}, {0.5, 0.7, 0.5}, 25, 0)
white = gr.material({1.0, 1.0, 1.0}, {0.5, 0.5, 0.5}, 25, 0)

scene_root = gr.node('root')
-- scene_one = gr.node('one')
-- scene_root:add_child(scene_one)

--[[
s1 = gr.nh_sphere('s1', {0, 0, -400}, 100)
scene_root:add_child(s1)
scene_one:add_child(s1)
s1:set_material(mat1)

s2 = gr.nh_sphere('s2', {200, 50, -100}, 150)
s1:add_child(s2)
s2:set_material(mat1)

s3 = gr.nh_sphere('s3', {0, -1200, -500}, 1000)
scene_root:add_child(s3)
s3:set_material(mat2)

s4 = gr.nh_sphere('s4', {-100, 25, -300}, 50)
scene_root:add_child(s4)
s4:set_material(mat3)

s5 = gr.nh_sphere('s5', {0, 100, -250}, 25)
scene_root:add_child(s5)
s5:set_material(mat1)
]]-- 

plane = gr.cube('plane')
plane:scale(1200, 200, 1000)
plane:translate(-500, -250, -500)
scene_root:add_child(plane)
plane:set_material(white)

s4 = gr.cone('s4')
scene_root:add_child(s4)
s4:set_material(mat3)
-- s4:rotate('Z', 90)
s4:scale(50, 100, 50)
-- s4:rotate('X', -30)
-- s4:rotate('Z', 50)
s4:translate(-100, 0, 200)

s5 = gr.cylinder('s4')
scene_root:add_child(s5)
s5:set_material(mat1)
--s5:rotate('X', 30)
--s5:rotate('Z', 90)
s5:scale(50, 50, 50)
s5:translate(100, 0, 200)

white_light = gr.light({-100.0, 200.0, 500.0}, {0.9, 0.9, 0.9}, {1, 0, 0})
-- magenta_light = gr.light({400.0, 100.0, 150.0}, {0.7, 0.0, 0.7}, {1, 0, 0})

gr.render(scene_root, 'simple.png', 512, 512,
	  {0, 100, 800}, {0, 0, -800}, {0, 1, 0}, 50,
	  {0.3, 0.3, 0.3}, {white_light, magenta_light})
