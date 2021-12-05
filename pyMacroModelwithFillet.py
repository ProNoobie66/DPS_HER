# author Rohit
# Macro-model for HER with code
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

Mdb()

ConeAngle = 30.0
#creating materials

mdb.models['Model-1'].Material(name='DPSteel')
mdb.models['Model-1'].materials['DPSteel'].Elastic(table=((201596.0, 0.3), ))
mdb.models['Model-1'].materials['DPSteel'].Plastic(table=((391.09558, 0.0), (449.3053137, 0.003045919), (494.4013464, 0.00930571),
														  (531.8412594, 0.018590825), (565.7182048, 0.028843932), (594.3309203, 0.040275795),
														  (619.3229891, 0.052320065), (641.6317208, 0.064745475), (660.7841464,0.077497753),
														  (678.3416248, 0.090470826), (694.0209775, 0.103655957), (707.5533536, 0.116944257),
														  (720.3426705, 0.130388436), (731.2873431, 0.143904754), (734.9915566, 0.152185816),
														  (740.5952913, 0.157490893), (747.8894289, 0.171112221), (752.1646528, 0.184914886),
														  (753.0648037, 0.198742115), (753.2674725, 0.202191939)))
print("Material DP Steel created")

#creating sections

mdb.models['Model-1'].HomogeneousSolidSection(name='DPSteelSection', material='DPSteel', thickness=None)
print("Section DP Steel created")

#creating plate

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=100.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(50.0, 0.0))
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(5.0, 0.0))
p = mdb.models['Model-1'].Part(name='Plate', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Plate']
p.BaseSolidExtrude(sketch=s, depth=1.0)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
e = p.edges
p.Round(radius=0.5, edgeList=(e.findAt((5.0,0.0,0.0),),))

print("Plate with fillet created")

#partitioning the plate

p = mdb.models['Model-1'].parts['Plate']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f.findAt((35,0,1),), sketchUpEdge=e.findAt((5,0,1),), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 1.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=282.07, gridSpacing=7.05, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(25.0, 0.0))
pickedFaces = f.findAt((35,0,1),)
e1 = p.edges.findAt((5,0.0,1),)
p.PartitionFaceBySketch(sketchUpEdge= e1, faces=pickedFaces, sketch=s1)
#p.PartitionFaceBySketch(sketchUpEdge=e.findAt((25.0,0.0,1.0),), faces=pickedFaces, sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
c = p.cells
pickedCells = c.getByBoundingBox(-50,-50,-1,50,50,2)#getSequenceFromMask(mask=('[#1 ]', ), )
pickedEdges =(e.findAt((25,0,1),), )
p.PartitionCellByExtrudeEdge(line=d1[4], cells=pickedCells, edges=pickedEdges, sense=REVERSE)
print("Plate circular partition complete")

#creating conical punch

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=100.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -50.0), point2=(0.0, 50.0))
s.FixedConstraint(entity=g.findAt((0.0,1.0),))#[2])
s.Line(point1=(0.0, 0.0), point2=(20.0, -25.0))
s.FixedConstraint(entity=v.findAt((0.0,0.0),))
s.AngularDimension(line1=g.findAt((20.0,-25.0),), line2=g.findAt((0.0,1.0),), textPoint=(5, -15), value=ConeAngle)
#g[3],g[2]
p = mdb.models['Model-1'].Part(name='Punch', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['Punch']
p.AnalyticRigidSurfRevolve(sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
v1, e1, d1, n = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v1.findAt((0.0,0.0,0.0),))
print("Punch created")

#assigning sections
p = mdb.models['Model-1'].parts['Plate']
c = p.cells
cells = c.getByBoundingBox(-50,-50,-1,50,50,2)#SequenceFromMask(mask=('[#3 ]',), )
region = regionToolset.Region(cells=cells)
p = mdb.models['Model-1'].parts['Plate']
p.SectionAssignment(region=region, sectionName='DPSteelSection', offset=0.0,
					offsetType=MIDDLE_SURFACE, offsetField='',
					thicknessAssignment=FROM_SECTION)
print("section assigned")

#creating assembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate']
a.Instance(name='Plate-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Punch']
a.Instance(name='Punch-1', part=p, dependent=ON)
a.rotate(instanceList=('Punch-1',), axisPoint=(50.0, 0.0, 0.0),
		 axisDirection=(-100.0, 0.0, 0.0), angle=-90.0)
a = mdb.models['Model-1'].rootAssembly
s = a.instances['Punch-1'].faces
side1Faces1 = s.findAt(((0.0,0.0,0.0),))#getSequenceFromMask(mask=('[#1 ]',), )
a.Surface(side1Faces=side1Faces1, name='MasSur')
s = a.instances['Plate-1'].faces
side1Faces1 = s.findAt(((5,0,0.9),)) + s.findAt(((6.0,0.0,0.0),)) + s.findAt(((5.1464466,0.0,0.1464466),))#getSequenceFromMask(mask=('[#38 ]',), )
a.Surface(side1Faces=side1Faces1, name='SlavSur')
print("assembly and contact surface created")

#creating step
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', maxNumInc=1000, initialInc=0.001, minInc=1e-08, nlgeom=ON)
print("Step created")

#creating interaction property
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
	formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
	pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
	table=((0.2,),), shearStressLimit=None, maximumElasticSlip=FRACTION,
	fraction=0.005, elasticSlipStiffness=None)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
	pressureOverclosure=HARD, allowSeparation=ON,
	constraintEnforcementMethod=DEFAULT)

#creating interactions
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1',
	createStepName='Step-1', master=a.surfaces['MasSur'],
	slave=a.surfaces['SlavSur'], sliding=FINITE, thickness=ON,
	interactionProperty='IntProp-1', adjustMethod=NONE,
	initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
print("Interactions created")

#creating boundry conditions
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Plate-1'].faces
faces1 = f1.findAt(((45,0,1),))+f1.findAt(((45,0,0),)) #getSequenceFromMask(mask=('[#6 ]',), )
region = a.Set(faces=faces1, name='BCSet-1')
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial',
								 region=region, localCsys=None)
r1 = a.instances['Punch-1'].referencePoints
refPoints = r1.findAt((0,0,0))#[2],)
region = a.Set(referencePoints=(refPoints,), name='BCSet-2')
mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1',
									 region=region, u1=0.0, u2=0.0, u3=25.0, ur1=0.0, ur2=0.0, ur3=0.0,
									 amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='',
									 localCsys=None)
print("Boundry conditions created")

#creating mesh
# p = mdb.models['Model-1'].parts['Plate']
# p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
# p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
# c = p.cells
# pickedCells = c.getByBoundingBox(-150,-150,-1.5,150,150,1.5)#SequenceFromMask(mask=('[#3 ]', ), )
# d = p.datums
# p.PartitionCellByDatumPlane(datumPlane=d[7], cells=pickedCells)
# c = p.cells
# pickedCells = c.getByBoundingBox(-150,-150,-1.5,150,150,1.5)
# p.PartitionCellByDatumPlane(datumPlane=d[8], cells=pickedCells)
# e = p.edges
# pickedEdges = e.findAt(((40,0,0),)) + e.findAt(((40,0,1),)) + e.findAt(((0,40,0),)) + e.findAt(((0,40,1),)) + e.findAt(((-40,0,0),)) + e.findAt(((-40,0,1),)) + e.findAt(((0,-40,0),)) + e.findAt(((0,-40,1),))
# p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
# pickedEdges = e.findAt(((35.355339,35.355339,0),)) + e.findAt(((35.355339,35.355339,1),)) + e.findAt(((-35.355339,35.355339,0),)) + e.findAt(((-35.355339,35.355339,1),)) + e.findAt(((-35.355339,-35.355339,0),)) + e.findAt(((-35.355339,-35.355339,1),)) + e.findAt(((35.355339,-35.355339,0),)) + e.findAt(((35.355339,-35.355339,1),)) #e.getSequenceFromMask(mask=('[#41000000 #412a001 ]', ), )
# p.seedEdgeByNumber(edges=pickedEdges, number=10, constraint=FINER)
# pickedEdges1 = e.getSequenceFromMask(mask=('[#20020240 ]', ), )
# pickedEdges2 = e.getSequenceFromMask(mask=('[#880 #500 ]', ), )
# p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1,
#         end2Edges=pickedEdges2, minSize=0.05, maxSize=1.0, constraint=FINER)
# pickedEdges = e.getSequenceFromMask(mask=('[#8000401 #800 ]', ), )
# p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
# p.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
# p.generateMesh()
p = mdb.models['Model-1'].parts['Plate']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[7], cells=pickedCells)
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#f ]', ), )
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[8], cells=pickedCells)
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#200a00a #14002 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
pickedEdges = e.getSequenceFromMask(mask=('[#41000000 #412a001 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=10, constraint=FINER)
pickedEdges1 = e.getSequenceFromMask(mask=('[#20020240 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#880 #500 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1,
        end2Edges=pickedEdges2, minSize=0.05, maxSize=1.0, constraint=FINER)
pickedEdges = e.getSequenceFromMask(mask=('[#8000401 #800 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
p.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
print("mesh generated")

#creating job
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)
print("Job created")


print("working")
