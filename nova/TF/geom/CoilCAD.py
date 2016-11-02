import sys
import json
import math
import numpy as np

import seaborn as sns
colors = sns.color_palette('Set2',12)
rgba = list(colors[3])
rgba.append(0)  # add alpha

from OCC.gp import gp_Pnt, gp_OX, gp_Vec, gp_Trsf, gp_DZ, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Pnt2d, gp_Dir2d, gp_Ax2d, gp_Dir
from OCC.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.GCE2d import GCE2d_MakeSegment
from OCC.Geom import Geom_Plane, Geom_CylindricalSurface, Handle_Geom_Plane 
from OCC.Geom import Handle_Geom_Surface, Geom_BezierCurve
from OCC.GeomFill import GeomFill_Pipe
from OCC.BRepFill import BRepFill_PipeShell
from OCC.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Geom2d import Geom2d_Ellipse, Geom2d_TrimmedCurve, Geom2d_BezierCurve
from OCC.Geom2d import Handle_Geom2d_Ellipse, Handle_Geom2d_Curve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
from OCC.BRepPrimAPI import BRepPrimAPI_MakePrism, BRepPrimAPI_MakeCylinder
from OCC.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid 
from OCC.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCC.BRepLib import breplib
from OCC.BRep import BRep_Tool_Surface, BRep_Builder
from OCC.TopoDS import topods, TopoDS_Edge, TopoDS_Compound, TopoDS_Solid
from OCC.TopoDS import TopoDS_Builder
from OCC.TopExp import TopExp_Explorer
from OCC.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCC.TopTools import TopTools_ListOfShape
import OCC.STEPControl as STEPControl
from OCC.GProp import GProp_GProps
from OCC.BRepGProp import (brepgprop_LinearProperties,
                           brepgprop_SurfaceProperties,
                           brepgprop_VolumeProperties)

from OCC.Display.SimpleGui import *
from OCC.TColgp import TColgp_Array1OfPnt
import OCCUtils
from OCCUtils.Common import GpropsFromShape
from OCCUtils import Construct

sys.path.insert(0,r'D:/Code/Nova/nova/TF/geom')

def add_line(point):
    p = []
    for i in range(2):
        p.append(gp_Pnt(point[i]['r'],0,point[i]['z']))
    edge = BRepBuilderAPI_MakeEdge(p[0],p[1]).Edge()
    return edge

with open('D:/Code/Nova/Data/salome_input.json','r') as f:  # _S16L
    data = json.load(f)
loop,cs,pf,nTF = data['p'],data['section'],data['pf'],data['nTF']
color = data['color']
PFsupport,CSsupport = data['PFsupport'],data['CSsupport']

w = []
for segment in loop:
    Pnt = TColgp_Array1OfPnt(1,4)
    for i in range(4):
        point = segment['p{:d}'.format(i)]
        Pnt.SetValue(i+1,gp_Pnt(point['r'],0,point['z']))
    curve = Geom_BezierCurve(Pnt)
    w.append(BRepBuilderAPI_MakeEdge(curve.GetHandle()).Edge())
    
inner_edge = add_line([loop[3]['p3'],loop[0]['p0']])
inner_curve = BRepBuilderAPI_MakeWire(inner_edge)
upper_curve = BRepBuilderAPI_MakeWire(w[0])
lower_curve = BRepBuilderAPI_MakeWire(w[3])

 # outboard loop
dflat = (loop[1]['p3']['r']-loop[2]['p0']['r'])**2+\
(loop[1]['p3']['z']-loop[2]['p0']['z'])**2
if dflat > 0:
    outboard_flat = add_line([loop[1]['p3'],loop[2]['p0']])
    outer_curve = BRepBuilderAPI_MakeWire(*[w[1],outboard_flat,w[2]])
    full_curve = BRepBuilderAPI_MakeWire(*[w[0],w[1],outboard_flat,w[2],w[3],
                                           inner_edge])
else:
    outer_curve = BRepBuilderAPI_MakeWire(*[w[1],w[2]])
    full_curve = BRepBuilderAPI_MakeWire(*[w[0],w[1],w[2],w[3]])
    full_curve = BRepBuilderAPI_MakeWire(full_curve.Wire(),inner_edge)
   
# extract winding pack dimesions
ro,zo = loop[0]['p0']['r'],loop[0]['p0']['z']
rtop,ztop = loop[0]['p3']['r'],loop[0]['p3']['z']
width,depth = cs['winding_pack']['width'],cs['winding_pack']['depth']
inboard,outboard = cs['case']['inboard'],cs['case']['outboard']
external = cs['case']['external']
side,nose = cs['case']['side'],cs['case']['nose']

# create winding pack upper
wp_upper = BRepBuilderAPI_MakePolygon() 
for dr,dy in zip([0,0,-1,-1],[-1,1,1,-1]):
    wp_upper.Add(gp_Pnt(ro+dr*width-inboard,dy*depth/2,zo))
wp_upper.Close()
wp_upper = wp_upper.Wire()

# create winding pack top
wp_top = BRepBuilderAPI_MakePolygon() 
for dy,dz in zip([1,1,-1,-1],[1,0,0,1]):
    wp_top.Add(gp_Pnt(rtop,dy*depth/2,ztop+outboard+dz*width))
wp_top.Close()
wp_top = wp_top.Wire()

# create case inboard
case_upper = BRepBuilderAPI_MakePolygon() 
theta = math.pi/nTF
rsep = (depth/2+side)/math.tan(theta)
rnose = ro-(width+inboard+nose)
if rsep <= rnose:
    ynose = depth/2+side
else:
    ynose = rnose*math.tan(theta)
case_upper.Add(gp_Pnt(ro,depth/2+side,zo))
if rsep > rnose:
    case_upper.Add(gp_Pnt(rsep,depth/2+side,zo))
case_upper.Add(gp_Pnt(rnose,ynose,zo))
case_upper.Add(gp_Pnt(rnose,-ynose,zo))
if rsep > rnose:
    case_upper.Add(gp_Pnt(rsep,-depth/2-side,zo))
case_upper.Add(gp_Pnt(ro,-depth/2-side,zo))
case_upper.Close()
case_upper = case_upper.Wire()

case_top = BRepBuilderAPI_MakePolygon() # create case outboard
for dy,dz in zip([1,1,-1,-1],[1,0,0,1]):
    case_top.Add(gp_Pnt(rtop,dy*(depth/2+side),
                        ztop+dz*(width+outboard+external)))
case_top.Close() 
case_top = case_top.Wire()

# transform inner sections
trans = gp_Trsf()
trans.SetTranslation(gp_Pnt(loop[0]['p0']['r'],0,loop[0]['p0']['z']),
                     gp_Pnt(loop[3]['p3']['r'],0,loop[3]['p3']['z']))
case_lower = BRepBuilderAPI_Transform(case_upper,trans).Shape()
case_lower = topods.Wire(case_lower)  # wire
wp_lower = BRepBuilderAPI_Transform(wp_upper,trans).Shape()
wp_lower = topods.Wire(wp_lower)  # wire

# transform mid sections
trans = gp_Trsf()
r1,z1 = loop[0]['p0']['r'],loop[0]['p0']['z']
r2,z2 = loop[3]['p3']['r'],loop[3]['p3']['z']
rm,zm = r1+(r2-r1)/2,z1+(z2-z1)/2
trans.SetTranslation(gp_Pnt(r1,0,z1),gp_Pnt(rm,0,zm))
case_mid = BRepBuilderAPI_Transform(case_upper,trans).Shape()
case_mid = topods.Wire(case_mid)  # wire
wp_mid = BRepBuilderAPI_Transform(wp_upper,trans).Shape()
wp_mid = topods.Wire(wp_mid)  # wire

# mirror outer sections
mirror = gp_Trsf()
mirror.SetMirror(gp_Pnt(loop[0]['p3']['r'],0,loop[0]['p3']['z'])) 
wp_bottom = BRepBuilderAPI_Transform(wp_top,mirror).Shape()
wp_bottom = topods.Wire(wp_bottom)  # wire
case_bottom = BRepBuilderAPI_Transform(case_top,mirror).Shape()
case_bottom = topods.Wire(case_bottom)  # wire

# translate outer sections
trans = gp_Trsf()
trans.SetTranslation(gp_Pnt(loop[0]['p3']['r'],0,loop[0]['p3']['z']),
                     gp_Pnt(loop[3]['p0']['r'],0,loop[3]['p0']['z']))
wp_bottom = BRepBuilderAPI_Transform(wp_bottom,trans).Shape()
wp_bottom = topods.Wire(wp_bottom)  # wire
case_bottom = BRepBuilderAPI_Transform(case_bottom,trans).Shape()
case_bottom = topods.Wire(case_bottom)  # wire

# rotate outboard
rotate = gp_Trsf()
axis = gp_Ax1(gp_Pnt(loop[0]['p3']['r'],0,loop[0]['p3']['z']),gp_Dir(0,1,0))
angle = -loop[2]['p0']['t']
rotate.SetRotation(axis,angle)
case_outer = BRepBuilderAPI_Transform(case_top,rotate).Shape()
wp_outer = BRepBuilderAPI_Transform(wp_top,rotate).Shape()

# translate outboard
trans = gp_Trsf()
trans.SetTranslation(gp_Pnt(loop[0]['p3']['r'],0,loop[0]['p3']['z']),
                     gp_Pnt(loop[1]['p3']['r'],0,loop[1]['p3']['z']))
case_outer = BRepBuilderAPI_Transform(case_outer,trans).Shape()
case_outer = topods.Wire(case_outer)  # wire
wp_outer = BRepBuilderAPI_Transform(wp_outer,trans).Shape()
wp_outer= topods.Wire(wp_outer)  # wire

TFprofile = {}
TFprofile['case'] = {'upper':[case_upper,case_top],
                     'outer':[case_top,case_outer,case_bottom],
                     'lower':[case_bottom,case_lower],
                     'inner':[case_lower,case_mid,case_upper]}
TFprofile['wp'] = {'upper':[wp_upper,wp_top],
                     'outer':[wp_top,wp_outer,wp_bottom],
                     'lower':[wp_bottom,wp_lower],
                     'inner':[wp_lower,wp_mid,wp_upper]}
TFcurve = {'upper':upper_curve,'outer':outer_curve,'lower':lower_curve,
           'inner':inner_curve}
           
#solid = {}   
builder = TopoDS_Builder()
TF = {}
for part in TFprofile:
    TF[part] = TopoDS_Compound()
    builder.MakeCompound(TF[part])
    for curve in TFprofile[part]:
        pipe = BRepFill_PipeShell(TFcurve[curve].Wire())
        for profile in TFprofile[part][curve]:
            pipe.Add(profile)
        pipe.Build()
        pipe.MakeSolid()
        builder.Add(TF[part],pipe.Shape())
  
TFcage = {}
rotate = gp_Trsf()
axis = gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1))
angle = 2*np.pi/nTF
rotate.SetRotation(axis,angle)

for part in TFprofile:
    TFcage[part] = TopoDS_Compound()
    builder.MakeCompound(TFcage[part])
    stamp = TF[part]
    for i in range(nTF):
        stamp = BRepBuilderAPI_Transform(stamp,rotate).Shape()
        builder.Add(TFcage[part],stamp)
        
        
#a = Construct.make_polygon([[0,0],[0,5],[5,5],[5,0]],closed=True)


# Export to STEP
step_export = STEPControl.STEPControl_Writer()
step_export.Transfer(TFcage['case'],STEPControl.STEPControl_AsIs)
step_export.Write('TF.stp')

prop = GpropsFromShape(TF['case'])
print('gprop',prop.volume().Mass())

prop = GProp_GProps()
brepgprop_VolumeProperties(TF['case'],prop,1e-6)
print('case volume',prop.Mass())

prop = GProp_GProps()
brepgprop_VolumeProperties(TF['wp'],prop,1e-6)
print('wp volume',prop.Mass())

prop = GProp_GProps()
brepgprop_VolumeProperties(TFcage['wp'],prop,1e-6)
print('cold mass volume',prop.Mass())



display,start_display = init_display()[:2]
display.DisplayColoredShape(TFcage['case'],'BLUE')
display.DisplayColoredShape(TFcage['wp'],'RED')
display.DisplayColoredShape(case_mid,'BLUE')
display.DisplayColoredShape(full_curve.Wire(),'BLUE')
#display.DisplayColoredShape(shell['case']['comp'].Shape(),'RED')


#rgba = list(colors[0])
#rgba.append(0)  # add alpha
#display.DisplayColoredShape(myBody,Quantity_Color(*rgba))
start_display()

'''


OY = geompy.MakeVectorDXDYDZ(0,1,0) 
OZ = geompy.MakeVectorDXDYDZ(0,0,1)   
ztop,zo,dt = CSsupport['ztop'],CSsupport['zo'],CSsupport['dt']
rnode,ynode = [],[]
for point in ['wp','nose']:
    rnode.append(CSsupport['r{}'.format(point)])
    ynode.append(CSsupport['y{}'.format(point)])
vseat = []
for sign in [-1,1]:  
    v = []
    v.append(geompy.MakeVertex(rnode[0],sign*ynode[0],ztop))
    v.append(geompy.MakeVertex(rnode[1],sign*ynode[1],ztop))
    vseat.append(v[::sign])
    v.append(geompy.MakeVertex(rnode[1],sign*ynode[1]-sign*dt,ztop))
    v.append(geompy.MakeVertex(rnode[0],sign*ynode[0]-sign*dt,ztop))
    v.append(v[0])
    support_loop = geompy.MakePolyline(v)
    face = geompy.MakeFaceWires([support_loop],1)
    solid = geompy.MakePrismVecH(face,OZ,-0.999*(ztop-zo))
    shell['case']['solid'] = geompy.MakeFuseList([shell['case']['solid'],
                                                  solid],True,True)
vseat.append(vseat[0])
seat_loop = geompy.MakePolyline(vseat)
face = geompy.MakeFaceWires([seat_loop],1)
solid = geompy.MakePrismVecH(face,OZ,-dt)
shell['case']['solid'] = geompy.MakeFuseList([shell['case']['solid'],solid],
                                             True,True)

PFsup = {'dt':side,'n':2}
for name in PFsupport:
    v = []
    for node in PFsupport[name]:
        vertex = geompy.MakeVertex(node[0],0,node[1])
        v.append(vertex)
    v.append(v[0])
    support_loop = geompy.MakePolyline(v)
    face = geompy.MakeFaceWires([support_loop],1)
    solid = geompy.MakePrismVecH2Ways(face,OY,PFsup['dt']/2)
    if PFsup['n'] > 1:
        rnode = PFsupport[name][0][0]
        if rnode < loop[0]['p3']['r']:
            shift = ynose+(depth/2+side-ynose)*\
            math.acos((loop[0]['p3']['r']-rnode)/
                        (loop[0]['p3']['r']-loop[0]['p0']['r']))/(math.pi/2)
        else:
            shift = depth/2+side
        solid = geompy.MakeTranslation(solid,0,-shift+PFsup['dt']/2,0)
        space = (2*shift-PFsup['dt'])/(PFsup['n']-1)
        solid = geompy.MakeMultiTranslation1D(solid,OY,space,PFsup['n'])                                               
        shell['case']['solid'] = geompy.MakeFuseList([shell['case']['solid'],
                                                     solid],True,True)

    v = []
    for node,so in zip(PFsupport[name][:2],[1,-1]):
        for s1 in [1,-1]:
            vertex = geompy.MakeVertex(node[0],so*s1*shift,node[1])
            v.append(vertex)
    v.append(v[0])
    support_loop = geompy.MakePolyline(v)
    face = geompy.MakeFaceWires([support_loop],1)
    if PFsupport[name][0][1] > PFsupport[name][2][1]:
        sign = -1
    else:
        sign = 1
    solid = geompy.MakePrismVecH(face,OZ,sign*PFsup['dt'])
    shell['case']['solid'] = geompy.MakeFuseList([shell['case']['solid'],
                                                     solid],True,True)

#shell['case']['solid'] = geompy.MakeCutList(shell['case']['solid'],\
#[shell['wp']['solid']],True)

for p in profiles:
    solid_id = geompy.addToStudy(shell[p]['solid'],'{}_solid'.format(p))
    gg.createAndDisplayGO(solid_id)
   
mag_part = geompy.MakeCompound([shell['wp']['solid'],shell['case']['solid']])
mag_part_id = geompy.addToStudy(mag_part,'TF')
gg.createAndDisplayGO(mag_part_id)
    
TF_cage = geompy.MultiRotate1DNbTimes(mag_part,None,nTF)
c = color[0]
TF_cage.SetColor(SALOMEDS.Color(c[0],c[1],c[2]))
rotate_id = geompy.addToStudy(TF_cage,'TF_cage')
gg.createAndDisplayGO(rotate_id)

coil_set = [TF_cage]
for i in range(len(pf)):
    c = color[i+1]
    coil = 'Coil{:d}'.format(i)
    coil_face = geompy.MakeFaceHW(pf[coil]['dz'],pf[coil]['dr'],3)
    coil_face = geompy.MakeTranslation(coil_face,
                                       pf[coil]['r'],0,pf[coil]['z'])
    coil = geompy.MakeRevolution(coil_face,OZ,360*math.pi/180.0)
    coil.SetColor(SALOMEDS.Color(c[0],c[1],c[2]))
    coil_set.append(coil)
    coil_id = geompy.addToStudy(coil,'coil_{:d}'.format(i))
    gg.createAndDisplayGO(coil_id)
 
coil_set = geompy.MakeCompound(coil_set)
geompy.addToStudy(coil_set,'coil_set')
geompy.ExportSTEP(coil_set,'D:/Code/Nova/nova/TF/geom/coil_set.step',
                  GEOM.LU_METER)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
'''