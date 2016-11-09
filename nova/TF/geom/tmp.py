from OCC.Display.SimpleGui import init_display
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox

from OCC.Display.WebGl import x3dom_renderer

display, start_display, add_menu, add_function_to_menu = init_display(backend_str='qt-pyqt4')
my_box = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()

#x3d = Construct.compound([TFcage['wp'],TFcage['case'],PFcage,GScage]) 
my_renderer = x3dom_renderer.X3DomRenderer()
my_renderer.DisplayShape(my_box)

display.DisplayShape(my_box, update=True)
start_display()
