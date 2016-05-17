# Recorded script from Mayavi2

from numpy import array
from mayavi import mlab
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
scene = engine.scenes[0]
vtk_file_reader = engine.open(u'out.vtk', scene)
from mayavi.modules.surface import Surface
surface = Surface()
engine.add_filter(surface, vtk_file_reader)
from mayavi.modules.iso_surface import IsoSurface
iso_surface1 = IsoSurface()
vtk_file_reader = engine.scenes[0].children[0]
engine.add_filter(iso_surface1, vtk_file_reader)
scene.scene.x_minus_view()
hand=mlab.gcf(engine)
mlab.figure(hand,bgcolor=(1,1,1))
surface.actor.property.opacity = 0.78
scene.scene.save(u'out.png')
