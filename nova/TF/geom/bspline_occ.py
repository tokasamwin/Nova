#!/usr/bin/env python

"""
Simple example of bspline surface created directly (not using Bezier surface)
"""

from OCC.gp import *
from OCC.Geom import *
from OCC.TColGeom import *
from OCC.TColgp import * 
from OCC.TColStd import *
from OCC.Display.SimpleGui import init_display

def bspline_surface():
    """
    Try to create B-spline surface directly
    """

    # Set U and V degree to 2
    udeg = 2
    vdeg = 2

    # Non-periodic surface
    uperiod = False
    vperiod = False

    # Create 2D array of poles (control points)
    poles = TColgp_Array2OfPnt(1, udeg + 1, 1, vdeg + 1)
    poles.SetValue(1, 1, gp_Pnt(1, 1, 1))
    poles.SetValue(1, 2, gp_Pnt(2, 1, 2))
    poles.SetValue(1, 3, gp_Pnt(3, 1, 1))
    poles.SetValue(2, 1, gp_Pnt(1, 2, 1))
    poles.SetValue(2, 2, gp_Pnt(2, 2, 2))
    poles.SetValue(2, 3, gp_Pnt(3, 2, 0))
    poles.SetValue(3, 1, gp_Pnt(1, 3, 2))
    poles.SetValue(3, 2, gp_Pnt(2, 3, 1))
    poles.SetValue(3, 3, gp_Pnt(3, 3, 0))

    # Length of uknots and umult has to be same
    # Same rule is for vknots and vmult
    uknot_len = umult_len = 2
    vknot_len = vmult_len = 2

    # Knots for U and V direction
    uknots = TColStd_Array1OfReal(1, uknot_len)
    vknots = TColStd_Array1OfReal(1, vknot_len)

    # Main curves begins and ends at first and last points
    uknots.SetValue(1, 0.0)
    uknots.SetValue(2, 1.0)
    vknots.SetValue(1, 0.0)
    vknots.SetValue(2, 1.0)

    # Multiplicities for U and V direction
    umult = TColStd_Array1OfInteger(1, umult_len)
    vmult = TColStd_Array1OfInteger(1, vmult_len)

    # First and last multiplicities are set to udeg + 1 (vdeg respectively),
    # because we want main curves to start and finish on the first and
    # the last points
    umult.SetValue(1, udeg + 1)
    umult.SetValue(2, udeg + 1)
    vmult.SetValue(1, vdeg + 1)
    vmult.SetValue(2, vdeg + 1)

    # Some other rules, that has to hold:
    # poles.ColLength == sum(umult(i)) - udeg - 1 (e.g.: 3 == 6 - 2 - 1)

    # Try to create surface
    BSPLSURF = Geom_BSplineSurface(poles, uknots, vknots, umult, vmult, udeg, vdeg, uperiod, vperiod)

    # Display surface
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.EraseAll()
    display.DisplayShape(BSPLSURF.GetHandle(), update=True)
    start_display()

if __name__ == '__main__':
    bspline_surface()