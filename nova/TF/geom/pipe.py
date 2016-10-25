# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.8.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'D:/Code/Nova/nova/TF/geom')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
geomObj_2 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
geomObj_3 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
sk = geompy.Sketcher2D()
sk.addPoint(0.000000, 0.000000)
sk.addSegmentRelative(50.000000, 0.000000)
sk.addArcRadiusLength(100.000000, 90.000000)
Sketch_1 = sk.wire(geomObj_3)
[Edge_1,Edge_2] = geompy.ExtractShapes(Sketch_1, geompy.ShapeType["EDGE"], True)
[Vertex_5,Vertex_6] = geompy.ExtractShapes(Edge_1, geompy.ShapeType["VERTEX"], True)
[Vertex_3,Vertex_4] = geompy.ExtractShapes(Edge_2, geompy.ShapeType["VERTEX"], True)
Edge_3 = geompy.MakeTangentOnCurve(Edge_2, 1)
[Vertex_1,Vertex_2] = geompy.ExtractShapes(Edge_3, geompy.ShapeType["VERTEX"], True)
Wire_1 = geompy.MakeWire([Edge_3, Edge_1, Edge_2], 1e-007)
Wire_1_vertex_3 = geompy.GetSubShape(Wire_1, [3])
Wire_1_edge_2 = geompy.GetSubShape(Wire_1, [2])
Circle_1 = geompy.MakeCircle(Wire_1_vertex_3, Wire_1_edge_2, 30)
Wire_1_vertex_4 = geompy.GetSubShape(Wire_1, [4])
Wire_1_edge_5 = geompy.GetSubShape(Wire_1, [5])
Circle_2 = geompy.MakeCircle(Wire_1_vertex_4, Wire_1_edge_2, 100)
Wire_1_vertex_6 = geompy.GetSubShape(Wire_1, [6])
Wire_1_edge_7 = geompy.GetSubShape(Wire_1, [7])
Circle_3 = geompy.MakeCircle(Wire_1_vertex_6, Wire_1_edge_7, 10)
Wire_1_vertex_8 = geompy.GetSubShape(Wire_1, [8])
Circle_4 = geompy.MakeCircle(Wire_1_vertex_8, Wire_1_edge_7, 50)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Sketch_1, 'Sketch_1' )
geompy.addToStudyInFather( Sketch_1, Edge_1, 'Edge_1' )
geompy.addToStudyInFather( Sketch_1, Edge_2, 'Edge_2' )
geompy.addToStudy( Edge_3, 'Edge_3' )
geompy.addToStudyInFather( Edge_3, Vertex_1, 'Vertex_1' )
geompy.addToStudyInFather( Edge_3, Vertex_2, 'Vertex_2' )
geompy.addToStudy( Wire_1, 'Wire_1' )
geompy.addToStudyInFather( Wire_1, Wire_1_vertex_3, 'Wire_1:vertex_3' )
geompy.addToStudyInFather( Wire_1, Wire_1_edge_2, 'Wire_1:edge_2' )
geompy.addToStudy( Circle_1, 'Circle_1' )
geompy.addToStudyInFather( Wire_1, Wire_1_vertex_4, 'Wire_1:vertex_4' )
geompy.addToStudyInFather( Wire_1, Wire_1_edge_5, 'Wire_1:edge_5' )
geompy.addToStudy( Circle_2, 'Circle_2' )
geompy.addToStudyInFather( Wire_1, Wire_1_vertex_6, 'Wire_1:vertex_6' )
geompy.addToStudyInFather( Wire_1, Wire_1_edge_7, 'Wire_1:edge_7' )
geompy.addToStudy( Circle_3, 'Circle_3' )
geompy.addToStudyInFather( Wire_1, Wire_1_vertex_8, 'Wire_1:vertex_8' )
geompy.addToStudy( Circle_4, 'Circle_4' )
geompy.addToStudyInFather( Edge_1, Vertex_5, 'Vertex_5' )
geompy.addToStudyInFather( Edge_2, Vertex_3, 'Vertex_3' )
geompy.addToStudyInFather( Edge_2, Vertex_4, 'Vertex_4' )
geompy.addToStudyInFather( Edge_1, Vertex_6, 'Vertex_6' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
