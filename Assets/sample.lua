-- A simple scene with some miscellaneous geometry.

mat1 = gr.material({0.7, 1.0, 0.7}, {0.5, 0.7, 0.5}, 25)
mat2 = gr.material({0.5, 0.5, 0.5}, {0.5, 0.7, 0.5}, 25)
mat3 = gr.material({1.0, 0.6, 0.1}, {0.5, 0.7, 0.5}, 25)
mat4 = gr.material({0.7, 0.6, 1.0}, {0.5, 0.4, 0.8}, 25)
mat5 = gr.material({0.9, 0.4, 0.7}, {0.5, 0.2, 0.9}, 25)

scene_root = gr.node('root')

s1 = gr.sphere('s1')
scene_root:add_child(s1)
s1:set_material(mat4)
s1:scale(100, 20,100)
s1:translate(200,0, 50)

s2 = gr.cube('s2')
scene_root:add_child(s2)
s2:set_material(mat2)
s2:scale(100,100,500)
s2:rotate('z',-45)
s2:rotate('y',-45)
s2:translate(0,100,0)

-- A small stellated dodecahedron.
steldodec = gr.mesh( 'dodec', 'smstdodeca.obj' )
steldodec:set_material(mat1)
scene_root:add_child(steldodec)
steldodec:translate(200,100,0)

nhs = gr.nh_sphere('nhs', {-100, -300, -400}, 350)
scene_root:add_child(nhs)
nhs:set_material(mat5)

nhc = gr.nh_box('nhc', {250, -500, -300}, 200)
scene_root:add_child(nhc)
nhc:set_material(mat1)


white_light = gr.light({-100.0, 150.0, 400.0}, {0.9, 0.9, 0.9}, {1, 0, 0})
magenta_light = gr.light({400.0, 100.0, 150.0}, {0.7, 0.4, 0.4}, {1, 0, 0})

gr.render(scene_root, 'sample.png', 500, 500,
	  {0, 0, 800}, {0, 0, -1}, {0, 1, 0}, 50,
	  {0.3, 0.3, 0.3}, {white_light, magenta_light})
