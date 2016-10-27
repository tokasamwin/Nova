import sys
import salome
import json
import math
gg = salome.ImportComponentGUI("GEOM")
from itertools import count
vindex = count(0)
iindex = count(0)
findex = count(0)
sindex = count(0)
cindex = count(0)

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0,r'D:/Code/Nova/nova/TF/geom')

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

def add_line(point):
    v = []
    for i in range(2):
        vertex = geompy.MakeVertex(point[i]['r'],0,point[i]['z'])
        v.append(vertex)
    line = geompy.MakeLineTwoPnt(v[0],v[1])
    return line

with open('D:/Code/Nova/Data/salome_input.json','r') as f:  # _S16L
    data = json.load(f)
loop,cs,pf,nTF = data['p'],data['section'],data['pf'],data['nTF']
color = data['color']
PFsupport = data['PFsupport']

geompy = geomBuilder.New(theStudy)
w = []
for segment in loop:
    v = []
    for i in range(4):
        point = segment['p{:d}'.format(i)]
        vertex = geompy.MakeVertex(point['r'],0,point['z'])
        v.append(vertex)
    curve = geompy.MakeBezier(v,False)
    w.append(curve)
    
inboard_flat = add_line([loop[3]['p3'],loop[0]['p0']])
dflat = (loop[1]['p3']['r']-loop[2]['p0']['r'])**2+\
(loop[1]['p3']['z']-loop[2]['p0']['z'])**2
if dflat > 0:
    outboard_flat = add_line([loop[1]['p3'],loop[2]['p0']])
    outer_loop = geompy.MakeWire([w[1],outboard_flat,w[2]])
else:
    outer_loop = geompy.MakeWire([w[1],w[2]])

ro,zo = loop[0]['p0']['r'],loop[0]['p0']['z']
rtop,ztop = loop[0]['p3']['r'],loop[0]['p3']['z']
width,depth = cs['winding_pack']['width'],cs['winding_pack']['depth']
inboard,outboard = cs['case']['inboard'],cs['case']['outboard']
external = cs['case']['external']
side,nose = cs['case']['side'],cs['case']['nose']

v = []  # create winding pack upper
for dr,dy in zip([0,0,-1,-1],[1,-1,-1,1]):
    vertex = geompy.MakeVertex(ro+dr*width-inboard,dy*depth/2,zo)
    v.append(vertex)
v.append(v[0])
wp_upper = geompy.MakePolyline(v)

v = []  # create winding pack top
for dy,dz in zip([1,1,-1,-1],[0,1,1,0]):
    vertex = geompy.MakeVertex(rtop,dy*depth/2,ztop+outboard+dz*width)
    v.append(vertex)
v.append(v[0])
wp_top = geompy.MakePolyline(v)
    
v = []  # create case inboard
theta = math.pi/nTF
rsep = (depth/2+side)/math.tan(theta)

rnose = ro-(width+inboard+nose)
if rsep <= rnose:
    ynose = depth/2+side
else:
    ynose = rnose*math.tan(theta)
v.append(geompy.MakeVertex(ro,depth/2+side,zo))
if rsep > rnose:
    v.append(geompy.MakeVertex(rsep,depth/2+side,zo))    
v.append(geompy.MakeVertex(rnose,ynose,zo)) 
v.append(geompy.MakeVertex(rnose,-ynose,zo)) 
if rsep > rnose:
    v.append(geompy.MakeVertex(rsep,-depth/2-side,zo)) 
v.append(geompy.MakeVertex(ro,-depth/2-side,zo))
v.append(v[0])
case_upper = geompy.MakePolyline(v)

v = []  # create case outboard
for dy,dz in zip([1,1,-1,-1],[0,1,1,0]):
    vertex = geompy.MakeVertex(rtop,dy*(depth/2+side),
                               ztop+dz*(width+outboard+external))
    v.append(vertex)
v.append(v[0])
case_top = geompy.MakePolyline(v)  

#upper and lower vertex
vupper = geompy.MakeVertex(loop[0]['p0']['r'],0,loop[0]['p0']['z'])
vlower = geompy.MakeVertex(loop[3]['p3']['r'],0,loop[3]['p3']['z'])
case_lower = geompy.MakeTranslationTwoPoints(case_upper,vupper,vlower)
wp_lower = geompy.MakeTranslationTwoPoints(wp_upper,vupper,vlower)

#top and bottom vertex
vtop = geompy.MakeVertex(loop[0]['p3']['r'],0,loop[0]['p3']['z'])
vbottom = geompy.MakeVertex(loop[3]['p0']['r'],0,loop[3]['p0']['z'])
wp_bottom = geompy.MakeTranslationTwoPoints(wp_top,vtop,vbottom,False)
wp_bottom = geompy.MakeMirrorByPoint(wp_bottom,vbottom)
case_bottom = geompy.MakeTranslationTwoPoints(case_top,vtop,vbottom,False)
case_bottom = geompy.MakeMirrorByPoint(case_bottom,vbottom)

# make shells
shell = {}
profiles = {'case':{'upper':case_upper,'top':case_top,
                    'bottom':case_bottom,'lower':case_lower},
            'wp':{'upper':wp_upper,'top':wp_top,
                    'bottom':wp_bottom,'lower':wp_lower}}
                    
for p in profiles:
    shell[p] = {}
    shell[p]['upper'] = geompy.MakePipeWithDifferentSections(\
    [profiles[p]['upper'],profiles[p]['top']],[vupper,vtop],w[0],\
    theWithContact=0,theWithCorrection=0)
    shell[p]['lower'] = geompy.MakePipeWithDifferentSections(\
    [profiles[p]['bottom'],profiles[p]['lower']],[vbottom,vlower],w[-1],\
    theWithContact=0,theWithCorrection=0)
    shell[p]['outer'] = geompy.MakePipe(profiles[p]['top'],outer_loop)
    shell[p]['inner'] = geompy.MakePipe(profiles[p]['upper'],inboard_flat)
    shell[p]['shell'] = geompy.MakeShell([shell[p]['upper'],shell[p]['lower'],
                                shell[p]['outer'],shell[p]['inner']])
    shell[p]['solid'] = geompy.MakeSolid([shell[p]['shell']])

OY = geompy.MakeVectorDXDYDZ(0,1,0) 
OZ = geompy.MakeVectorDXDYDZ(0,0,1)   
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
#gg.createAndDisplayGO(coil_set_id)

geompy.ExportSTEP(coil_set,'D:/Code/Nova/nova/TF/geom/coil_set.step',
                  GEOM.LU_METER)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
