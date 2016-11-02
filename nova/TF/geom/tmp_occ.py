
import OCC

from OCC.Display.SimpleGui import init_display

from OCC.gp import gp_Pnt,gp_Vec
from OCC.GC import GC_MakeArcOfCircle,GC_MakeSegment
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge as MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire as MakeWire
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace as MakeFace
from OCC.BRepPrimAPI import BRepPrimAPI_MakePrism as MakePrism
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder as MakeCylinder
from OCC.BRepFilletAPI import BRepFilletAPI_MakeFillet as MakeFillet
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform as Transform
from OCC.TopExp import TopExp_Explorer as Explorer
from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse as Fuse
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid as MakeThickSolid
from OCC.Quantity import Quantity_Color
from OCC.BRep import BRep_Tool
from OCC.Standard import Standard_Type

import seaborn as sns
colors = sns.color_palette('Set2',12)

'''
#from OCC.gp import *
#from OCC.GC import *
#from OCC.BRepBuilderAPI import *
#from OCC.TopoDS import *
#from OCC.TopExp import *
#from OCC.TopAbs import *
#from OCC.BRepAlgoAPI import *
#from OCC.BRepFilletAPI import *
#from OCC.BRepPrimAPI import *
from OCC.Geom import *
from OCC.Geom2d import *
from OCC.GCE2d import *
#from OCC.BRepOffsetAPI import *
from OCC.IGESControl import *
from OCC.TopTools import *
import OCC.BRepLib as BRepLib
import OCC.BRep as BRep
import OCC.STEPControl as STEPControl
import math
'''
display,start_display = init_display()[:2]

myWidth = 50
myHeight = 70
myThickness = 30

# Profile : Define Support Points
aPnt1 = gp_Pnt(-myWidth/2,0,0)
aPnt2 = gp_Pnt(-myWidth/2,-myThickness/4,0)
aPnt3 = gp_Pnt(0,-myThickness/2,0)
aPnt4 = gp_Pnt(myWidth/2,-myThickness/4,0)
aPnt5 = gp_Pnt(myWidth/2,0,0)

# Profile : Define the Geometry
aArcOfCircle = GC_MakeArcOfCircle(aPnt2,aPnt3,aPnt4)
aSegment1 = GC_MakeSegment(aPnt1,aPnt2)
aSegment2 = GC_MakeSegment(aPnt4,aPnt5)

# Profile : Define the Topology
aEdge1 = MakeEdge(aSegment1.Value())
aEdge2 = MakeEdge(aArcOfCircle.Value())
aEdge3 = MakeEdge(aSegment2.Value())
aWire  = MakeWire(aEdge1.Edge(),aEdge2.Edge(),aEdge3.Edge())

# Complete Profile
xAxis = OCC.gp.gp_OX()
aTrsf = OCC.gp.gp_Trsf()
aTrsf.SetMirror(xAxis)
aBRepTrsf = Transform(aWire.Wire(),aTrsf)
aMirroredShape = aBRepTrsf.Shape()
aMirroredWire = OCC.TopoDS.TopoDS_Shape(aMirroredShape)
aMirroredWire = OCC.TopoDS.topods().Wire(aMirroredWire)

mkWire = MakeWire()
mkWire.Add(aWire.Wire())
mkWire.Add(aMirroredWire)
myWireProfile = mkWire.Wire()

# Body : Prism the Profile
myFaceProfile = MakeFace(myWireProfile)
if myFaceProfile.IsDone():
    bottomFace = myFaceProfile.Face()
aPrismVec = gp_Vec(0,0,myHeight)
myBody = MakePrism(myFaceProfile.Shape(),aPrismVec)

# Body : Apply Fillets
mkFillet = MakeFillet(myBody.Shape())
aEdgeExplorer = Explorer(myBody.Shape(),OCC.TopAbs.TopAbs_EDGE)
while aEdgeExplorer.More():
    aEdge = OCC.TopoDS.topods_Edge(aEdgeExplorer.Current())
    mkFillet.Add(myThickness/12,aEdge)
    aEdgeExplorer.Next()
myBody = mkFillet.Shape()

# Body : Add the Neck   
neckLocation = OCC.gp.gp_Pnt(0,0,myHeight)
neckNormal = OCC.gp.gp_DZ()
neckAx2 = OCC.gp.gp_Ax2(neckLocation,neckNormal)
myNeckRadius = myThickness/4
myNeckHeight = myHeight/10
MKCylinder = MakeCylinder(neckAx2,myNeckRadius,myNeckHeight)
myNeck = MKCylinder.Shape()
myBody = Fuse(myBody,myNeck).Shape()

# Body : Create a Hollowed Solid
faceToRemove = OCC.TopoDS.TopoDS_Face()
aFaceExplorer = OCC.TopExp.TopExp_Explorer(myBody,OCC.TopAbs.TopAbs_FACE)
facesToRemove = OCC.TopTools.TopTools_ListOfShape()
i = 0
#while aFaceExplorer.More():
for i in range(31):
    i += 1
    aFace = OCC.TopoDS.topods_Face(aFaceExplorer.Current())
    aSurface = BRep_Tool.Surface(aFace)
    #print(aSurface.GetObject().DynamicType(),OCC.Geom.Geom_Plane)
    print(aFace)
    aFaceExplorer.Next()
#facesToRemove.Append(aFace) 

myBody = MakeThickSolid(myBody,facesToRemove,-myThickness/50,1e-3).Shape()

#display.DisplayColoredShape(myWireProfile,'BLUE')
rgba = list(colors[0])
rgba.append(0)  # add alpha
display.DisplayColoredShape(myBody,Quantity_Color(*rgba))
start_display()
 
'''
# Threading : Create Surfaces
neckAx2_bis = gp_Ax3(neckLocation , neckNormal)
aCyl1 = Geom_CylindricalSurface(neckAx2_bis , myNeckRadius * 0.99)
aCyl2 = Geom_CylindricalSurface(neckAx2_bis , myNeckRadius * 1.05)
aCyl1_handle = aCyl1.GetHandle()
aCyl2_handle = aCyl2.GetHandle()

# Threading : Define 2D Curves
aPnt = gp_Pnt2d(2. * math.pi , myNeckHeight / 2.)
aDir = gp_Dir2d(2. * math.pi , myNeckHeight / 4.)
aAx2d = gp_Ax2d(aPnt , aDir)
aMajor = 2. * math.pi
aMinor = myNeckHeight / 10.
anEllipse1 = Geom2d_Ellipse(aAx2d , aMajor , aMinor)
anEllipse2 = Geom2d_Ellipse(aAx2d , aMajor , aMinor / 4.)
anEllipse1_handle = anEllipse1.GetHandle()
anEllipse2_handle = anEllipse2.GetHandle()
aArc1 = Geom2d_TrimmedCurve(anEllipse1_handle, 0 , math.pi)
aArc2 = Geom2d_TrimmedCurve(anEllipse2_handle, 0 , math.pi)
aArc1_handle = aArc1.GetHandle()
aArc2_handle = aArc2.GetHandle()
anEllipsePnt1 = anEllipse1.Value(0)
anEllipsePnt2 = anEllipse1.Value(math.pi)
aSegment = GCE2d_MakeSegment(anEllipsePnt1 , anEllipsePnt2)

# Threading : Build Edges and Wires
aEdge1OnSurf1 = BRepBuilderAPI_MakeEdge( aArc1_handle , aCyl1_handle)
aEdge2OnSurf1 = BRepBuilderAPI_MakeEdge( aSegment.Value() , aCyl1_handle)
aEdge1OnSurf2 = BRepBuilderAPI_MakeEdge( aArc2_handle , aCyl2_handle)
aEdge2OnSurf2 = BRepBuilderAPI_MakeEdge( aSegment.Value() , aCyl2_handle)

threadingWire1 = BRepBuilderAPI_MakeWire(aEdge1OnSurf1.Edge() , aEdge2OnSurf1.Edge())
threadingWire2 = BRepBuilderAPI_MakeWire(aEdge1OnSurf2.Edge() , aEdge2OnSurf2.Edge())
BRepLib.BRepLib().BuildCurves3d(threadingWire1.Shape())
BRepLib.BRepLib().BuildCurves3d(threadingWire2.Shape())

# Create Threading
aTool = BRepOffsetAPI_ThruSections(True)
aTool.AddWire(threadingWire1.Wire())
aTool.AddWire(threadingWire2.Wire())
aTool.CheckCompatibility(False)
myThreading = aTool.Shape()

#Building the resulting compound
aRes = TopoDS_Compound()
aBuilder = BRep.BRep_Builder()
aBuilder.MakeCompound(aRes)
aBuilder.Add(aRes, myBody)
aBuilder.Add(aRes, myThreading)

# Export to STEP
step_export = STEPControl.STEPControl_Writer()
step_export.Transfer(aRes,STEPControl.STEPControl_AsIs)
step_export.Write('bottle.stp') # Export bottle and threads

display.DisplayColoredShape(myBody,'RED') # Display bottle but not threads
display.DisplayColoredShape(bottomFace,'YELLOW')
display.DisplayColoredShape(myWireProfile,'BLUE')
start_display()

'''