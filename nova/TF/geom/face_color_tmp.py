from __future__ import print_function

from random import random

from OCC.AIS import AIS_ColoredShape
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Display.OCCViewer import color
from OCC.Display.SimpleGui import init_display

from OCCUtils import Topo

display, start_display, add_menu, add_function_to_menu = init_display()

my_box = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()

ais = AIS_ColoredShape(my_box)


for fc in Topo(my_box).faces():
    # set a custom color per-face
    ais.SetCustomColor(fc, color(random(), random(), random()))


display.Context.Display(ais.GetHandle())
display.FitAll()

start_display()