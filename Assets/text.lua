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

s2 = gr.sphere('s2')
s2:scale(40, 40, 40)
s2:rotate('Y', -90)
s2:translate(160, -85,300)
scene:add_child(s2)
s2:set_material('jupiter.png')

s1 = gr.sphere('s1')
s1:scale(30, 30, 30)
s1:rotate('Y', -90)
s1:translate(-5, -80, 300)
scene:add_child(s1)
s1:set_material('mars.png')

s3 = gr.sphere('s3')
s3:scale(25, 22, 25)
s3:rotate('Y', -90)
s3:translate(50, 75,300)
scene:add_child(s3)
s3:set_material('uranus.png')

s4 = gr.sphere('s4')
s4:scale(15, 15, 15)
s4:rotate('Y', -90)
s4:translate(-25, 105, 300)
scene:add_child(s4)
s4:set_material('Mercury.png')

s5 = gr.sphere('s5')
s5:scale(15, 15, 15)
s5:rotate('Y', -135)
s5:translate(125, 30, 500)
scene:add_child(s5)
s5:set_material('me.png')

-- s5 = gr.nh_sphere('s5', {105, 105, 500}, 50)
-- scene:add_child(s5)
-- s5:set_material('me.png')

plane = gr.cube('plane')
plane:scale(512, 512, 512)
plane:translate(-256, -256, -200)
scene:add_child(plane)
plane:set_material('solar_system.png')

white_light = gr.light({-600.0, -30.0, 600.0}, {0.7, 0.7, 0.7}, {1, 0, 0})
yellow_light = gr.light({150.0, 0.0, -800.0}, {0.7, 0.7, 0.5}, {1, 0, 0})

gr.render(scene, 'text.png', 512, 512,
	  {0, 0, 800}, {0, 0, -1}, {0, 1, 0}, 50,
	  {0.3, 0.3, 0.3}, {white_light, yellow_light})
