#author Rohit
#HET018 random sphere microstructre plate
# units in mm and N
# ignore this code
from abaqus import *
from abaqusConstants import *
import __main__


import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import math as m
import random as rm
import os

#changing working directory

os.chdir(r"C:\Users\Rohit\Desktop\IITB\DPS HER\HET018")

#function to check if point is in reqd region

def dist(x1,y1,z1=0,x2=0,y2=0,z2=0):
    d = m.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
    return d

#function to find the volume of sphere

def sphereVol(rad):
    vol = 4*m.pi*rad*rad*rad/3
    return vol

#input plate parameters

plateOuterR = 0.1250 #mm
plateInnerR = 0.025 #mm
plateThick = 0.005 #mm
plateVol = m.pi*plateThick*(plateOuterR*plateOuterR-plateInnerR*plateInnerR)

rqdvf = 0.001    #required volume fraction

#creating list of coordinates

x = []
y = []
z = []
r = []
totalVolM = 0
ns = 0
vf = 0

while vf<rqdvf:
    rp = rm.uniform(0.0015,0.0025)
    xp = rm.uniform(-0.125 + rp, 0.125 - rp)
    yp = rm.uniform(-0.125 + rp, 0.125 - rp)
    zp = rm.uniform(0 + rp,0.005 - rp)
    d0 = dist(xp, yp)
    inter = 0
    if d0 > 0.025 + rp and d0 < 0.125 - rp:
        for j in range(len(x)):
            d = dist(x[j], y[j], z[j], xp, yp, zp)
            if d < rp + r[j]:
                inter = 1
                break
        if inter == 0:
            r.append(rp)
            x.append(xp)
            y.append(yp)
            z.append(zp)
            print('point created', rp, xp, yp, zp)
            totalVolM = totalVolM + sphereVol(rp)
            print('volume of sphere = ',sphereVol(rp))
            print('total volume of sphere = ',totalVolM)
            plateVol = plateVol - sphereVol(rp)
            print('volume of disc = ',plateVol)
            vf = totalVolM/plateVol
            print('volume fraction =',vf)
            ns = ns + 1
            print('number of sphere = ',ns)

print('radius and center generated')


#creating materials

mdb.models['Model-1'].Material(name='Ferrite')
mdb.models['Model-1'].materials['Ferrite'].Elastic(table=((210000,0.3),))
mdb.models['Model-1'].materials['Ferrite'].Plastic(table=((302.0, 0.0), (304.0, 
        0.0001), (305.0, 0.0002), (307.0, 0.0004), (310.0, 0.0008), (317.0, 
        0.0016), (328.0, 0.0032), (344.0, 0.0064), (366.0, 0.0125), (398.0, 
        0.025), (443.0, 0.05), (499.0, 0.1), (536.0, 0.15), (562.0, 0.2) ))

mdb.models['Model-1'].Material(name='Martensite')
mdb.models['Model-1'].materials['Martensite'].Elastic(table=((225000,0.3),))
mdb.models['Model-1'].materials['Martensite'].Plastic(table=((1510.0, 0.0), 
        (1532.0, 0.0001), (1553.0, 0.0002), (1590.0, 0.0004), (1650.0, 0.0008), 
        (1744.0, 0.0016), (1872.0, 0.0032), (2024.0, 0.0064), (2162.0, 0.0125), 
        (2247.0, 0.025), (2268.0, 0.05), (2269.0, 0.1), (2269.0, 0.2)))

print("created materials")

#creating sections

mdb.models['Model-1'].HomogeneousSolidSection(name='FerriteSec',material='Ferrite', thickness=None)

mdb.models['Model-1'].HomogeneousSolidSection(name='MartensiteSec',material='Martensite', thickness=None)

print("created sections")

#creating parts

#punch

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.3)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=3)
s.ConstructionLine(point1=(0.0, -0.15), point2=(0.0, 0.15))
s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, 0.05), point2=(0.05, 0.0), direction=CLOCKWISE)
s.Line(point1=(0.05, 0.0), point2=(0.05, -0.1))
p = mdb.models['Model-1'].Part(name='Punch', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['Punch']
p.AnalyticRigidSurfRevolve(sketch=s)
del mdb.models['Model-1'].sketches['__profile__']
p.ReferencePoint(point=p.vertices.findAt((0.0,0.05,0.0)))

print("created punch")

#disc

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.6)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.25, 0.0))
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.025, 0.0))
p = mdb.models['Model-1'].Part(name='Disc', dimensionality=THREE_D,type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Disc']
p.BaseSolidExtrude(sketch=s, depth=0.005)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
f, e, d = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f.findAt((0.125,0.125,0.005)), sketchUpEdge=e.findAt((0.025,0.0,0.005),),sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.005))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=1.41,gridSpacing=0.03, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.125, 0.0))
f = p.faces
pickedFaces = f.findAt((0.15,0.15,0.005))
e, d = p.edges, p.datums
e1 = p.edges.findAt((0.025,0.0,0.005),)
p.PartitionFaceBySketch(sketchUpEdge= e1, faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
pickedCells = p.cells.findAt(((0.03,0.03,0.003),))
pickedEdges =(e.findAt((0.125,0.0,0.005),), )
p.PartitionCellByExtrudeEdge(line = d[3],cells=pickedCells,edges = pickedEdges,sense = REVERSE)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
c = p.cells.getByBoundingBox(-0.25,-0.25,-0.005,0.25,0.25,0.005)
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[6], cells=c)
c = p.cells.getByBoundingBox(-0.25,-0.25,-0.005,0.25,0.25,0.005)
p.PartitionCellByDatumPlane(datumPlane=d[5], cells=c)

print("created disc")

#creating assembly

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Disc']
a.Instance(name='Disc-1', part=p, dependent=ON)

#creating island
i = 0
while i < ns:
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=0.006)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.sketchOptions.setValues(decimalPlaces=4)
    s.ConstructionLine(point1=(0.0, -0.003), point2=(0.0, 0.003))
    s.FixedConstraint(entity=g[2])
    s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, r[i]), point2=(0.0, -1*r[i]), direction=CLOCKWISE)
    s.Line(point1=(0.0, -1*r[i]), point2=(0.0, r[i]))
    s.VerticalConstraint(entity=g[4], addUndoState=False)
    p = mdb.models['Model-1'].Part(name='Island-%d'%(i+1), dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Island-%d'%(i+1)]
    p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
    p = mdb.models['Model-1'].parts['Island-%d'%(i+1)]
    del mdb.models['Model-1'].sketches['__profile__']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    c = p.cells.findAt(((0.0,0.0,0.0),))
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[2], cells=c)
    c = p.cells.getByBoundingBox(-r[i], -r[i], -r[i], r[i], r[i], r[i])
    p.PartitionCellByDatumPlane(datumPlane=d[3], cells=c)
    c = p.cells.getByBoundingBox(-r[i],-r[i],-r[i],r[i],r[i],r[i])
    p.PartitionCellByDatumPlane(datumPlane=d[4], cells=c)
    print('created island-%d'%(i+1))
    a.Instance(name='Island-%d'%(i+1), part=p, dependent=ON)
    a.translate(instanceList=('Island-%d'%(i+1), ), vector=(x[i], y[i], z[i]))
    if i == 0:
        selectedInstances = (a.instances['Island-%d'%(i+1)],)
    else:
        sphereInstances = a.instances['Island-%d'%(i+1)]
        selectedInstances = selectedInstances + (sphereInstances,)
    i+= 1

selectedInstances = selectedInstances + (a.instances['Disc-1'],)

#merging instances
a.InstanceFromBooleanMerge(name='Plate', instances=selectedInstances,keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
print('merged islands')
a.regenerate()

p=mdb.models['Model-1'].parts['Punch']
a.Instance(name='Punch-1',part=p,dependent=ON)
a.rotate(instanceList=('Punch-1',),axisPoint=(0.25,0.0,0.0),axisDirection=(-0.5,0.0,0.0),angle=270.0)
a.translate(instanceList=('Punch-1',),vector=(0.0,0.0,-0.05))
a.regenerate()

print("created assembly")

i = 0
while i < ns:
	del mdb.models['Model-1'].parts['Island-%d'%(i+1)]
	i += 1

#assigning section

#creating Ferrite set
fp = ((0.025,0.025,0.003),(-0.025,0.025,0.003),(-0.025,-0.025,0.003),(0.025,-0.025,0.003),
      (0.125,0.125,0.003),(-0.125,0.125,0.003),(-0.125,-0.125,0.003),(0.125,-0.125,0.003),)
p = mdb.models['Model-1'].parts['Plate']
Fcell = p.cells.findAt((fp[0],))
i = 0
while i < 8:
    Fcell = Fcell + p.cells.findAt((fp[i],))
    i = i + 1
p.Set(cells=Fcell,name='Set-F')
region = p.sets['Set-F']
p.SectionAssignment(region=p.sets['Set-F'],sectionName='FerriteSec',offset=0.0,offsetType =MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)

#creating martensite set

mpl = []
i = 0
while i < ns:
    mpl.append((x[i] + r[i] / 3, y[i] + r[i] / 3, z[i] + r[i] / 3))
    mpl.append((x[i] - r[i] / 3, y[i] + r[i] / 3, z[i] + r[i] / 3))
    mpl.append((x[i] + r[i] / 3, y[i] - r[i] / 3, z[i] + r[i] / 3))
    mpl.append((x[i] - r[i] / 3, y[i] - r[i] / 3, z[i] + r[i] / 3))
    mpl.append((x[i] + r[i] / 3, y[i] + r[i] / 3, z[i] - r[i] / 3))
    mpl.append((x[i] - r[i] / 3, y[i] + r[i] / 3, z[i] - r[i] / 3))
    mpl.append((x[i] + r[i] / 3, y[i] - r[i] / 3, z[i] - r[i] / 3))
    mpl.append((x[i] - r[i] / 3, y[i] - r[i] / 3, z[i] - r[i] / 3))
    i = i + 1
mp = tuple(mpl)

i = 0
Mcell=p.cells.findAt((mp[i],))
while i < len(mp) :
    Mcell=Mcell+p.cells.findAt((mp[i],))
    i = i +1
p.Set(cells=Mcell, name='Set-M')
p.SectionAssignment(region = p.sets['Set-M'],sectionName='MartensiteSec',offset=0.0,offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)

print("assigned section")

a.regenerate()

#creating step

mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial',maxNumInc=1000, initialInc=0.001, minInc=1e-09)

print('created step')

#creating interaction

mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD,
        allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
        table=((0.2,),), shearStressLimit=None, maximumElasticSlip=FRACTION,
        fraction=0.005, elasticSlipStiffness=None)

f = a.instances['Plate-1'].faces
sur=f.findAt(((0.026,0.026,0.0),)) + f.findAt(((-0.026,0.026,0.0),)) + f.findAt(((-0.026,-0.026,0.0),)) + f.findAt(((0.026,-0.026,0.0),)) + f.findAt(((0.0176776699,0.0176776699,0.003),)) + f.findAt(((-0.0176776699,0.0176776699,0.003),))+ f.findAt(((-0.0176776699,-0.0176776699,0.003),))+ f.findAt(((0.0176776699,-0.0176776699,0.003),))
slaveSurf=a.Surface(side1Faces=sur,name='PlateSurf')

f=a.instances['Punch-1'].faces
side1Faces1=f.findAt(((0.0,0.0,0.0),))
masterSurf=a.Surface(side1Faces=side1Faces1, name='PunchSurf')

mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1',createStepName='Step-1',
                                                 master=masterSurf,slave=slaveSurf,
                                                 sliding=FINITE,thickness=ON,
                                                 interactionProperty='IntProp-1',
                                                 adjustMethod=NONE, initialClearance=OMIT,
                                                 datumAxis=None, clearanceRegion=None)

print("interaction created")

#creating boundary conditions

f=a.instances['Plate-1'].faces
encastface = f.findAt(((0.15,0.15,0.0),)) + f.findAt(((-0.15,0.15,0.0),)) + f.findAt(((-0.15,-0.15,0.0),)) + f.findAt(((0.15,-0.15,0.0),))
BCset1 = a.Set(faces = encastface, name = 'BC1')
mdb.models['Model-1'].EncastreBC(name='BC-1',createStepName='Initial',region = BCset1,localCsys=None)

r = a.instances['Punch-1'].referencePoints
refPoints = r.findAt((0.0,0.0,0.0))
region = a.Set(referencePoints = (refPoints,),name='BC2')
mdb.models['Model-1'].DisplacementBC(name = 'BC-2',createStepName = 'Step-1',region = region, u1 =0.0, u2= 0.0, u3 = 0.1,
                                     ur1 = 0.0, ur2 = 0.0, ur3 =0.0, amplitude = UNSET, fixed = OFF, distributionType = UNIFORM,
                                     fieldName ='', localCsys=None)

print("Boundary conditions created")

#creating mesh

p = mdb.models['Model-1'].parts['Plate']
e = p.edges
c = p.cells
pickedEdges = e.getByBoundingBox(0.15,0.15,-0.001,0.3,0.3,0.006)
p.seedEdgeByNumber(edges=pickedEdges, number=10, constraint=FINER)
pickedRegions = c.getByBoundingBox(-0.25,-0.25,-0.006,0.25,0.25,0.006)#SequenceFromMask(mask=('[#ffffffff:2 #f ]', ), )
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
cells = c.getByBoundingBox(-0.25,-0.25,-0.006,0.25,0.25,0.006)#getSequenceFromMask(mask=('[#ffffffff:2 #f ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,
       elemType3))
p.seedPart(size=0.007, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

print("created mesh")

a.regenerate()

#creating job

mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)

print("Job created")
session.viewports['Viewport: 1'].setValues(displayedObject = a)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap=session.viewports['Viewport: 1'].colorMappings['Material']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()

print('working')
