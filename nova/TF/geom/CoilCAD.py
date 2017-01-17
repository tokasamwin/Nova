import sys
import json
import math
import numpy as np
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
import seaborn as sns
colors = sns.color_palette('Paired',12)
QC = [[] for i in range(len(colors))]
for i,c in enumerate(colors):
    QC[i] = Quantity_Color(*c,Quantity_TOC_RGB)
    QC[i].ChangeIntensity(-50)
    
from collections import OrderedDict
from OCC.Display.WebGl import x3dom_renderer
from OCC.gp import gp_Pnt, gp_Trsf, gp_Ax1, gp_Ax2
from OCC.gp import gp_Dir
from OCC.GC import GC_MakeArcOfCircle
from OCC.Geom import Geom_BezierCurve
from OCC.BRepFill import BRepFill_PipeShell
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.BRepPrimAPI import BRepPrimAPI_MakeRevol
from OCC.TopoDS import topods
from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop_VolumeProperties
from OCC.Display.SimpleGui import init_display
from OCC.TColgp import TColgp_Array1OfPnt
from OCCUtils.Common import GpropsFromShape
from OCCUtils import Construct
from OCC import UnitsAPI
UnitsAPI.unitsapi_SetCurrentUnit('LENGTH','meter')
from amigo.IO import trim_dir
                                 
sys.path.insert(0,r'D:/Code/Nova/nova/TF/geom')

def add_line(point):
    p = []
    for i in range(2):
        p.append(gp_Pnt(point[i]['r'],0,point[i]['z']))
    edge = BRepBuilderAPI_MakeEdge(p[0],p[1]).Edge()
    return edge
    
datadir = trim_dir('../../../Data/') 
with open(datadir+'occ_input.json','r') as f:
    data = json.load(f)
loop,cs,pf,nTF = data['p'],data['section'],data['pf'],data['nTF']
color = data['color']
PFsupport,CSsupport = data['PFsupport'],data['CSsupport']
Gsupport,OISsupport = data['Gsupport'], data['OISsupport']

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
if dflat > 1e-3:
    outboard_flat = add_line([loop[1]['p3'],loop[2]['p0']])
    outer_curve = BRepBuilderAPI_MakeWire(*[w[1],outboard_flat,w[2]])
    full_curve = BRepBuilderAPI_MakeWire(*[w[0],w[1],outboard_flat,w[2]])
    full_curve = BRepBuilderAPI_MakeWire(full_curve.Wire(),w[3])
    full_curve = BRepBuilderAPI_MakeWire(full_curve.Wire(),inner_edge)
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
TFprofile['case'] = OrderedDict()
TFprofile['case']['upper'] = [case_upper,case_top]
TFprofile['case']['outer'] = [case_top,case_outer,case_bottom]
TFprofile['case']['lower'] = [case_bottom,case_lower]
TFprofile['case']['inner'] = [case_lower,case_mid,case_upper]

TFprofile['wp'] = OrderedDict()
TFprofile['wp']['upper'] = [wp_upper,wp_top]
TFprofile['wp']['outer'] = [wp_top,wp_outer,wp_bottom]
TFprofile['wp']['lower'] = [wp_bottom,wp_lower]
TFprofile['wp']['inner'] = [wp_lower,wp_mid,wp_upper]
TFcurve = {'upper':upper_curve,'outer':outer_curve,'lower':lower_curve,
           'inner':inner_curve}
           
TF = {}
for part in TFprofile:
    segments = []
    for curve in TFprofile[part]:
        pipe = BRepFill_PipeShell(TFcurve[curve].Wire())
        for profile in TFprofile[part][curve]:
            pipe.Add(profile)
        pipe.Build()
        segments.append(pipe.Shape())
    quilt = Construct.sew_shapes(segments)
    TF[part] = Construct.make_solid(quilt)

# outer inter-coil supports
OISloop = [[] for i in range(2)]
pnt = [[] for i in range(2)]
OISface = [[] for i in range(2)] 
yo = depth/2+side  # coil half thickness
ro = yo/np.arctan(theta)  # radial offset

for name in OISsupport:
    OISloop[0] = BRepBuilderAPI_MakePolygon()  # inter-coil support
    for node in OISsupport[name]:
        pnt[0] = gp_Pnt(node[0],yo,node[1])
        OISloop[0].Add(pnt[0])
    OISloop[0].Close()
    OISloop[0] = OISloop[0].Wire()  
    OISloop[1] = BRepBuilderAPI_MakePolygon()  # inter-coil support
    for node in OISsupport[name]:
        r1 = node[0]-ro
        dy = r1*np.sin(theta)
        dr = dy*np.tan(theta)
        pnt[1] = gp_Pnt(node[0]-dr,yo+dy,node[1])
        OISloop[1].Add(pnt[1])
    OISloop[1].Close()
    OISloop[1] = OISloop[1].Wire()  

    path = Construct.make_line(pnt[0],pnt[1])
    path = Construct.make_wire(path)
    shell = BRepFill_PipeShell(path)
    shell.Add(OISloop[0])
    shell.Add(OISloop[1])
    shell.Build()
    for i in range(2):
        OISface[i] = Construct.make_face(OISloop[i])
    quilt = Construct.sew_shapes([shell.Shape(),OISface[0],OISface[1]])  # 
    Lplate = Construct.make_solid(quilt)
    
    node = OISsupport[name][0]
    Rplate = Construct.mirror_axe2(Lplate,gp_Ax2(gp_Pnt(node[0],0,node[1]),
                                                 gp_Dir(0,1,0),gp_Dir(0,0,1)))
    TF['case'] = Construct.boolean_fuse(TF['case'],Lplate) 
    TF['case'] = Construct.boolean_fuse(TF['case'],Rplate)
          
# CS seat
ztop,zo,dt = CSsupport['ztop'],CSsupport['zo'],CSsupport['dt']
yfactor = 0.8
rnode,ynode = [],[]
for point in ['wp','nose']:
    rnode.append(CSsupport['r{}'.format(point)])
    ynode.append(CSsupport['y{}'.format(point)])
    
CSloop = BRepBuilderAPI_MakePolygon() # create CStop
CSloop.Add(gp_Pnt(rnode[0],yfactor*ynode[0],ztop))
CSloop.Add(gp_Pnt(rnode[1],yfactor*ynode[1],ztop))
CSloop.Add(gp_Pnt(rnode[1],-yfactor*ynode[1],ztop))
CSloop.Add(gp_Pnt(rnode[0],-yfactor*ynode[0],ztop))
CSloop.Close()
CSloop = CSloop.Wire()
CSface = Construct.make_face(CSloop)
CSbody = Construct.make_prism(CSface,Construct.gp_Vec(0,0,zo-ztop))
CSbody = topods.Solid(CSbody)
TF['case'] = Construct.boolean_fuse(TF['case'],CSbody)  # join seat to body

# GS base
node = Gsupport['base']
rGS = np.mean([node[0][0],node[1][0]])
GSwidth = node[1][0]-node[0][0]
zGS = node[1][1]
Gloop = BRepBuilderAPI_MakePolygon()
for n in node:
    Gloop.Add(gp_Pnt(n[0],-depth/2-side,n[1]))
Gloop.Close()
Gloop = Gloop.Wire()
Gface = Construct.make_face(Gloop)
Gbody = Construct.make_prism(Gface,Construct.gp_Vec(0,depth+2*side,0))
Gbody = topods.Solid(Gbody)
TF['case'] = Construct.boolean_fuse(TF['case'],Gbody)  # join seat to body

Gax = Construct.gp_Ax1(gp_Pnt(rGS,0,zGS-GSwidth/2),gp_Dir(1,0,0))
pin = Construct.gp_Circ()
pin.SetRadius(GSwidth/4)
pin.SetAxis(Gax)
pin = Construct.make_edge(pin)
pin = Construct.make_wire(pin)
pin = Construct.make_face(pin)
pin = Construct.make_prism(pin,Construct.gp_Vec(GSwidth,0,0))
pin = Construct.translate_topods_from_vector(pin,
                                             Construct.gp_Vec(-GSwidth/2,0,0))

arc = GC_MakeArcOfCircle(gp_Pnt(rGS,-GSwidth/2,zGS-GSwidth/2),
                         gp_Pnt(rGS,0,zGS-GSwidth),
                         gp_Pnt(rGS,GSwidth/2,zGS-GSwidth/2))
arc = Construct.make_edge(arc.Value())
base = BRepBuilderAPI_MakePolygon()
base.Add(gp_Pnt(rGS,-GSwidth/2,zGS-GSwidth/2))
base.Add(gp_Pnt(rGS,-GSwidth/2,zGS))
base.Add(gp_Pnt(rGS,GSwidth/2,zGS))
base.Add(gp_Pnt(rGS,GSwidth/2,zGS-GSwidth/2))
Gloop = Construct.make_wire([arc,base.Wire()])
Gface = Construct.make_face(Gloop)
Gbody = Construct.make_prism(Gface,Construct.gp_Vec(GSwidth/3,0,0))
Gbody = topods.Solid(Gbody)
Gtag = Construct.translate_topods_from_vector(Gbody,
                                             Construct.gp_Vec(-GSwidth/6,0,0))
Gbase = Construct.translate_topods_from_vector(Gtag,\
        Construct.gp_Vec(-GSwidth/2+GSwidth/6,0,0))
Gbase = Construct.boolean_fuse(Gbase,\
                               Construct.translate_topods_from_vector(Gtag,\
                               Construct.gp_Vec(+GSwidth/2-GSwidth/6,0,0)))
Gbase = Construct.boolean_fuse(Gbase,pin)
TF['case'] = Construct.boolean_fuse(TF['case'],Gbase)  # join GSbase to case

Gtag = Construct.rotate(Gtag,Gax,180)   
r2c = 1.5*GSwidth
rtube,ttube = GSwidth/3,GSwidth/9  # GS leg radius, thickness
rect = BRepBuilderAPI_MakePolygon()
rect.Add(gp_Pnt(rGS-GSwidth/6,-GSwidth/2,zGS-GSwidth))
rect.Add(gp_Pnt(rGS-GSwidth/6,GSwidth/2,zGS-GSwidth))
rect.Add(gp_Pnt(rGS+GSwidth/6,GSwidth/2,zGS-GSwidth))
rect.Add(gp_Pnt(rGS+GSwidth/6,-GSwidth/2,zGS-GSwidth))
rect.Close()
rect = rect.Wire()
rect_face = Construct.make_face(rect)
rpath = Construct.make_line(gp_Pnt(rGS,0,zGS-GSwidth),
                           gp_Pnt(rGS,0,zGS-r2c))
rpath = Construct.make_wire(rpath)
circ = Construct.gp_Circ()
circ.SetRadius(rtube)
circ.SetAxis(Construct.gp_Ax1(gp_Pnt(rGS,0,zGS-r2c),gp_Dir(0,0,1)))
circ = Construct.make_edge(circ)
circ = Construct.make_wire(circ)
circ_face = Construct.make_face(circ)

circ_cut = Construct.gp_Circ()
circ_cut.SetRadius(rtube-ttube)
circ_cut.SetAxis(Construct.gp_Ax1(gp_Pnt(rGS,0,zGS-r2c),gp_Dir(0,0,1)))
circ_cut = Construct.make_edge(circ_cut)
circ_cut = Construct.make_wire(circ_cut)
circ_face_cut = Construct.make_face(circ_cut)

pipe = BRepFill_PipeShell(rpath)
pipe.Add(rect)
pipe.Add(circ)
pipe.Build()
quilt = Construct.sew_shapes([pipe.Shape(),rect_face,circ_face])  # 
GStrans = Construct.make_solid(quilt)
GStrans = Construct.boolean_fuse(GStrans,Gtag)
GSt_upper = Construct.boolean_cut(GStrans,pin)  # make upper hole

alpha = np.arctan((Gsupport['radius']*np.tan(theta)-GSwidth)/
                     (Gsupport['zbase']-Gsupport['zfloor']-GSwidth/2))
Pin2pin = np.sqrt((Gsupport['radius']*np.tan(theta)-GSwidth)**2+
                     (Gsupport['zbase']-Gsupport['zfloor']-GSwidth/2)**2)

tube_body = Construct.make_prism(circ_face,\
            Construct.gp_Vec(0,0,-(Pin2pin+GSwidth-2*r2c)))
tube_body = topods.Solid(tube_body)

tube_body_cut = Construct.make_prism(circ_face_cut,\
                Construct.gp_Vec(0,0,-(Pin2pin+GSwidth-2*r2c)))
tube_body_cut = topods.Solid(tube_body_cut)
tube_body = Construct.boolean_cut(tube_body,tube_body_cut)

GSt_lower = Construct.rotate(GSt_upper,Gax,180)
GSt_lower = Construct.translate_topods_from_vector(GSt_lower,\
            Construct.gp_Vec(0,0,-Pin2pin))
GSt_lower = Construct.rotate(GSt_lower,\
            Construct.gp_Ax1(gp_Pnt(rGS,0,zGS),gp_Dir(0,0,1)),90)

leg = Construct.boolean_fuse(GSt_upper,tube_body)
leg = Construct.boolean_fuse(leg,GSt_lower)
leftleg = Construct.rotate(leg,Gax,-alpha*180/np.pi)
rightleg = Construct.rotate(leg,Gax,alpha*180/np.pi)
GSstrut = Construct.boolean_fuse(leftleg,rightleg)

# lower pin y-axis
circ = Construct.gp_Circ()
circ.SetRadius(rtube-ttube)
circ.SetAxis(Construct.gp_Ax1(gp_Pnt(rGS,-Gsupport['radius']*np.tan(theta),
                                     Gsupport['zfloor']),gp_Dir(0,1,0)))
circ = Construct.make_edge(circ)
circ = Construct.make_wire(circ)
circ_face = Construct.make_face(circ)
lower_pin = Construct.make_prism(circ_face,\
            Construct.gp_Vec(0,2*Gsupport['radius']*np.tan(theta),0))
lower_pin = topods.Solid(lower_pin)


PFsup = {'dt':side,'n':3}
for name in PFsupport:
    PFSloop = BRepBuilderAPI_MakePolygon()  # PF support
    for node in PFsupport[name]:
        PFSloop.Add(gp_Pnt(node[0],0,node[1]))
    PFSloop.Close()
    PFSloop = PFSloop.Wire()
    PFSface = Construct.make_face(PFSloop)
    PFSbody = Construct.make_prism(PFSface,Construct.gp_Vec(0,PFsup['dt'],0))
    PFSbody = Construct.translate_topods_from_vector(PFSbody,\
    Construct.gp_Vec(0,-PFsup['dt']/2,0))
    if PFsup['n'] > 1:
        rnode = PFsupport[name][0][0]
        if rnode < loop[0]['p3']['r']:
            shift = ynose+(depth/2+side-ynose)*\
            math.acos((loop[0]['p3']['r']-rnode)/
                        (loop[0]['p3']['r']-loop[0]['p0']['r']
                         +width+nose+inboard))/(math.pi/2)
        else:
            shift = depth/2+side
        PFSbody = Construct.translate_topods_from_vector(PFSbody,\
        Construct.gp_Vec(0,-shift+PFsup['dt']/2,0))    
        space = (2*shift-PFsup['dt'])/(PFsup['n']-1)
        TF['case'] = Construct.boolean_fuse(TF['case'],PFSbody)
        for i in np.arange(1,PFsup['n']):
            PFSbody = Construct.translate_topods_from_vector(PFSbody,\
            Construct.gp_Vec(0,space,0))
            TF['case'] = Construct.boolean_fuse(TF['case'],PFSbody)
    PFSloop = BRepBuilderAPI_MakePolygon()  # PF support
    for node,so in zip(PFsupport[name][:2],[1,-1]):
        for s1 in [1,-1]:
            PFSloop.Add(gp_Pnt(node[0],so*s1*shift,node[1]))
    PFSloop.Close()
    PFSloop = PFSloop.Wire()
    PFSface = Construct.make_face(PFSloop)
    if PFsupport[name][0][1] > PFsupport[name][2][1]:
        sign = -1
    else:
        sign = 1
    vector = np.array(PFsupport[name][2])-np.array(PFsupport[name][1])
    vector /= np.linalg.norm(vector) 
    vector *= PFsup['dt']  # sign*
    PFSbody = Construct.make_prism(PFSface,
                                   Construct.gp_Vec(vector[0],0,vector[1]))
    TF['case'] = Construct.boolean_fuse(TF['case'],PFSbody)  # join seat to body 
TF['case'] = Construct.boolean_cut(TF['case'],TF['wp'])

TFcage = {}
rotate = gp_Trsf()
axis = gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1))
angle = 2*np.pi/nTF
rotate.SetRotation(axis,angle)

for part in TFprofile:
    stamp,compound = TF[part],[]
    for i in range(nTF):
        stamp = BRepBuilderAPI_Transform(stamp,rotate).Shape()
        compound.append(stamp)
    TFcage[part] = Construct.compound(compound) 
        
stamp,compound = GSstrut,[]    
for i in range(nTF):
    stamp = BRepBuilderAPI_Transform(stamp,rotate).Shape()
    compound.append(stamp)
GScage = Construct.compound(compound) 
    
PFcoil = []
for i in range(len(pf)):
    c = color[i+1]
    PFloop = BRepBuilderAPI_MakePolygon() # create CStop
    coil = 'Coil{:d}'.format(i)
    r,z = pf[coil]['r'],pf[coil]['z']
    dr,dz = pf[coil]['dr'],pf[coil]['dz']
    for sr,sz in zip([-1,1,1,-1],[-1,-1,1,1]):
        PFloop.Add(gp_Pnt(sr/2*dr,0,sz/2*dz))
    PFloop.Close()
    PFloop = Construct.translate_topods_from_vector(PFloop.Wire(),
                                                    Construct.gp_Vec(r,0,z))
    PFface = Construct.make_face(PFloop)
    ax = gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1))
    PFcoil.append(BRepPrimAPI_MakeRevol(PFface,ax).Shape())
PFcage = Construct.compound(PFcoil) 

x3d = Construct.compound([TFcage['wp'],TFcage['case'],PFcage,GScage]) 
my_renderer = x3dom_renderer.X3DomRenderer()
my_renderer.DisplayShape(x3d)
       


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
display.DisplayColoredShape(TFcage['case'],QC[0])
display.DisplayColoredShape(GScage,QC[1])
display.DisplayColoredShape(PFcage,QC[5])
display.FitAll()
start_display()

# Export to STEP
scale = gp_Trsf()  # to mm
scale.SetScaleFactor(1000)
TF['wp'] = BRepBuilderAPI_Transform(TF['wp'],scale).Shape()
TF['case'] = BRepBuilderAPI_Transform(TF['case'],scale).Shape()
GSstrut = BRepBuilderAPI_Transform(GSstrut,scale).Shape()

from aocxchange import step_ocaf
export = step_ocaf.StepOcafExporter('./TF_{:d}.stp'.format(nTF))
export.add_shape(TF['wp'],color=colors[0],layer='winding_pack')
export.add_shape(topods.Compound(TF['case']),color=colors[1],layer='case')
export.add_shape(GSstrut,color=colors[2],layer='gravity_support')
export.write_file()