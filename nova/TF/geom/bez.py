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
import math
import SALOMEDS

def add_line(point):
    v = []
    for i in range(2):
        vertex = geompy.MakeVertex(point[i]['r'],0,point[i]['z'])
        v.append(vertex)
        geompy.addToStudy(vertex,'vertex_{}'.format(next(vindex)))
    line = geompy.MakeLineTwoPnt(v[0],v[1])
    line_id = geompy.addToStudy(line, 'line_{:d}'.format(next(iindex)))
    gg.createAndDisplayGO(line_id)
    return line

with open('D:/Code/Nova/Data/test_sal.vi','r') as f:
    loop = json.load(f)

with open('D:/Code/Nova/Data/CS_sal.vi','r') as f:
    cs = json.load(f)
    
nTF = 18  # store in json

geompy = geomBuilder.New(theStudy)
w = []
for segment in loop:
    v = []
    for i in range(4):
        point = segment['p{:d}'.format(i)]
        vertex = geompy.MakeVertex(point['r'],0,point['z'])
        v.append(vertex)
        geompy.addToStudy(vertex,'vertex_{}'.format(next(vindex)))
    curve = geompy.MakeBezier(v,False)
    w.append(curve)
    curve_id = geompy.addToStudy(curve, 'Curve_{:d}'.format(next(cindex)))
    gg.createAndDisplayGO(curve_id)
    
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
    geompy.addToStudy(vertex,'vertex_{:d}'.format(next(vindex)))
    v.append(vertex)
v.append(v[0])
wp_upper = geompy.MakePolyline(v)
wp_upper_id = geompy.addToStudy(wp_upper,'wp_upper')
gg.createAndDisplayGO(wp_upper_id)  

v = []  # create winding pack top
for dy,dz in zip([1,1,-1,-1],[0,1,1,0]):
    vertex = geompy.MakeVertex(rtop,dy*depth/2,ztop+outboard+dz*width)
    geompy.addToStudy(vertex,'vertex_{:d}'.format(next(vindex)))
    v.append(vertex)
v.append(v[0])
wp_top = geompy.MakePolyline(v)
wp_top_id = geompy.addToStudy(wp_top,'wp_top')
gg.createAndDisplayGO(wp_top_id)  
    
v = []  # create case inboard
theta = math.pi/nTF
rsep = (depth/2+side)/math.tan(theta)
rnose = ro-(width+inboard+nose)
ynose = rnose*math.tan(theta)
v.append(geompy.MakeVertex(ro,depth/2+side,zo))
v.append(geompy.MakeVertex(rsep,depth/2+side,zo))    
v.append(geompy.MakeVertex(rnose,ynose,zo)) 
v.append(geompy.MakeVertex(rnose,-ynose,zo)) 
v.append(geompy.MakeVertex(rsep,-depth/2-side,zo)) 
v.append(geompy.MakeVertex(ro,-depth/2-side,zo))
for vertex in v:
    geompy.addToStudy(vertex,'vertex_{:d}'.format(next(vindex)))
v.append(v[0])
case_upper = geompy.MakePolyline(v)
case_upper_id = geompy.addToStudy(case_upper,'case_upper')
gg.createAndDisplayGO(case_upper_id)  

v = []  # create case outboard
for dy,dz in zip([1,1,-1,-1],[0,1,1,0]):
    vertex = geompy.MakeVertex(rtop,dy*(depth/2+side),
                               ztop+dz*(width+outboard+external))
    geompy.addToStudy(vertex,'vertex_{:d}'.format(next(vindex)))
    v.append(vertex)
v.append(v[0])
case_top = geompy.MakePolyline(v)
case_top_id = geompy.addToStudy(case_top,'case_top')
gg.createAndDisplayGO(case_top_id)  

#upper and lower vertex
vupper = geompy.MakeVertex(loop[0]['p0']['r'],0,loop[0]['p0']['z'])
geompy.addToStudy(vupper,'vertex_{:d}'.format(next(vindex)))
vlower = geompy.MakeVertex(loop[3]['p3']['r'],0,loop[3]['p3']['z'])
geompy.addToStudy(vlower,'vertex_{:d}'.format(next(vindex)))

case_lower = geompy.MakeTranslationTwoPoints(case_upper,vupper,vlower)
case_lower_id = geompy.addToStudy(case_lower,'case_lower')
gg.createAndDisplayGO(case_lower_id)

wp_lower = geompy.MakeTranslationTwoPoints(wp_upper,vupper,vlower)
wp_lower_id = geompy.addToStudy(wp_lower,'wp_lower')
gg.createAndDisplayGO(wp_lower_id)

#top and bottom vertex
vtop = geompy.MakeVertex(loop[0]['p3']['r'],0,loop[0]['p3']['z'])
geompy.addToStudy(vtop,'vertex_{:d}'.format(next(vindex)))
vbottom = geompy.MakeVertex(loop[3]['p0']['r'],0,loop[3]['p0']['z'])
geompy.addToStudy(vbottom,'vertex_{:d}'.format(next(vindex)))

wp_bottom = geompy.MakeTranslationTwoPoints(wp_top,vtop,vbottom,False)
wp_bottom = geompy.MakeMirrorByPoint(wp_bottom,vbottom)
wp_bottom_id = geompy.addToStudy(wp_bottom,'case_bottom')
gg.createAndDisplayGO(wp_bottom_id)

case_bottom = geompy.MakeTranslationTwoPoints(case_top,vtop,vbottom,False)
case_bottom = geompy.MakeMirrorByPoint(case_bottom,vbottom)
case_bottom_id = geompy.addToStudy(case_bottom,'case_bottom')
gg.createAndDisplayGO(case_bottom_id)

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
    shell[p]['upper_id'] = geompy.addToStudy(shell[p]['upper'],\
    'shell+{:d}'.format(next(sindex)))
    gg.createAndDisplayGO(shell[p]['upper_id'])
    
    shell[p]['lower'] = geompy.MakePipeWithDifferentSections(\
    [profiles[p]['bottom'],profiles[p]['lower']],[vbottom,vlower],w[-1],\
    theWithContact=0,theWithCorrection=0)
    shell[p]['lower_id'] = geompy.addToStudy(shell[p]['lower'],\
    'shell+{:d}'.format(next(sindex)))
    gg.createAndDisplayGO(shell[p]['lower_id'])

    shell[p]['outer'] = geompy.MakePipe(profiles[p]['top'],outer_loop)
    shell[p]['outer_id'] = geompy.addToStudy(shell[p]['outer'],\
    'shell+{:d}'.format(next(sindex)))
    gg.createAndDisplayGO(shell[p]['outer_id'])

    shell[p]['inner'] = geompy.MakePipe(profiles[p]['upper'],inboard_flat)
    shell[p]['inner_id'] = geompy.addToStudy(shell[p]['inner'],\
    'shell+{:d}'.format(next(sindex)))
    gg.createAndDisplayGO(shell[p]['inner_id'])

#,case_bottom,case_lower
#,vbottom,vlower

#if salome.sg.hasDesktop():
#  salome.sg.updateObjBrowser(1)
