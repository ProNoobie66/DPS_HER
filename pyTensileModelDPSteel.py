# python script to create tensile model for DPSteel
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
import numpy as np


# function to find distance
def dist(x1, y1, z1=0, x2=0, y2=0, z2=0):
	d = m.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))
	return d


# function to find volume of sphere
def sphereVol(rad):
	vol = 4 * m.pi * rad * rad * rad / 3
	return vol


# gauge area parameter
L = 50.0  # um
B = 20.0  # um
H = 8.0  # um
gaugeVol = L * B * H
rqdvf = 0.1  # required volume fraction

# creating list of coordinates

x = []
y = []
z = []
r = []
totalVolM = 0
noOfParticle = 0
vf = 0

while vf < rqdvf:
	rp = rm.uniform(1.5, 2.5)
	xp = rm.uniform(-B / 2 + rp, B / 2 - rp)
	yp = rm.uniform(-L / 2 + rp, L / 2 - rp)
	zp = rm.uniform(-H / 2 + rp, H / 2 - rp)
	inter = 0
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
		print('volume of sphere = ', sphereVol(rp))
		print('total volume of sphere = ', totalVolM)
		gaugeVol = gaugeVol - sphereVol(rp)
		print('volume of disc = ', gaugeVol)
		vf = totalVolM / gaugeVol
		print('volume fraction =', vf)
		noOfParticle = noOfParticle + 1
		print('number of sphere = ', noOfParticle)

print('radius and center generated')
# creating a model
Mdb()
model = mdb.Model(name='Model-1', modelType=STANDARD_EXPLICIT)

# creating materials
model.Material(name='DP Steel')
model.materials['DP Steel'].Elastic(table=((0.2, 0.3),))
model.materials['DP Steel'].Plastic(table=((391.09558e-6, 0.0), (
	449.3053137e-6, 0.003045919), (494.4013464e-6, 0.00930571), (531.8412594e-6,
																 0.018590825), (565.7182048e-6, 0.028843932),
										   (594.3309203e-6, 0.040275795),
										   (619.3229891e-6, 0.052320065), (641.6317208e-6, 0.064745475),
										   (660.7841464e-6, 0.077497753), (678.3416248e-6, 0.090470826),
										   (694.0209775e-6, 0.103655957),
										   (707.5533536e-6, 0.116944257), (720.3426705e-6, 0.130388436),
										   (731.2873431e-6,
											0.143904754), (740.5952913e-6, 0.157490893), (747.8894289e-6, 0.171112221),
										   (752.1646528e-6, 0.184914886), (753.0648037e-6, 0.198742115),
										   (753.2674725e-6, 0.202191939)))

model.Material(name='Martensite')
model.materials['Martensite'].Elastic(table=((0.225, 0.3),))
model.materials['Martensite'].Plastic(table=((1506.83447, 0.0),
											 (1507.5, 1.7e-05), (1583.749217, 0.000142766), (1659.93414, 0.000308327),
											 (1736.151251, 0.000626444), (1812.368335, 0.001134407),
											 (1888.585418, 0.001891694),
											 (1964.802502, 0.002999662), (2041.019585, 0.004650164),
											 (2117.236669, 0.007276208),
											 (2193.453752, 0.012263619), (2269.57, 0.18943869)))
model.Material(name='Ferrite')
model.materials['Ferrite'].Elastic(table=((0.21, 0.3),))
model.materials['Ferrite'].Plastic(
	table=((303.6851816, 0.0), (322.1037187, 0.002096385), (340.5222556, 0.005395182), (358.9407922,
																						0.009946118),
		   (377.3593285, 0.015820343), (395.7778642, 0.023113982), (414.1963992, 0.031953744),
		   (432.6149334, 0.042505035),
		   (451.0334663, 0.054983778), (469.4519976, 0.069673965), (487.8705266, 0.0869544),
		   (506.2890524, 0.107340882), (524.7075735, 0.131555765), (543.1260875, 0.160649168),
		   (561.5445901, 0.196225817)))

print('Created Material')

# creating sections

model.HomogeneousSolidSection(name='DPSteelSection', material='DP Steel', thickness=None)
model.HomogeneousSolidSection(name='MartensiteSection', material='Martensite', thickness=None)
model.HomogeneousSolidSection(name='FerriteSection', material='Ferrite', thickness=None)

print('created section')

# creating tensile part sketch

s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.Line(point1=(-35.0, 0.0), point2=(35.0, 0.0))
s.Line(point1=(35.0, 0.0), point2=(35.0, 30.0))
c1 = s.Line(point1=(35.0, 30.0), point2=(10.0, 56.65))
idc1 = c1.id
c2 = s.Line(point1=(10.0, 56.65), point2=(10.0, 106.65))
idc2 = c2.id
c3 = s.Line(point1=(10.0, 106.65), point2=(35.0, 133.3))
idc3 = c3.id
s.Line(point1=(35.0, 133.3), point2=(35.0, 163.3))
s.Line(point1=(35.0, 163.3), point2=(-35.0, 163.3))
s.Line(point1=(-35.0, 163.3), point2=(-35.0, 133.3))
c4 = s.Line(point1=(-35.0, 133.3), point2=(-10.0, 106.65))
idc4 = c4.id
c5 = s.Line(point1=(-10.0, 106.65), point2=(-10.0, 56.65))
idc5 = c5.id
c6 = s.Line(point1=(-10.0, 56.65), point2=(-35.0, 30.0))
idc6 = c6.id
s.Line(point1=(-35.0, 30.0), point2=(-35.0, 0.0))

# creating fillet

s.FilletByRadius(radius=20.0, curve1=g[idc1], nearPoint1=(21.8456115722656, 46.6164588928223), curve2=g[idc2],
				 nearPoint2=(10.4230804443359, 63.1363372802734))
s.FilletByRadius(radius=20.0, curve1=g[idc2], nearPoint1=(10.4230804443359, 97.8850631713867), curve2=g[idc3],
				 nearPoint2=(14.9920959472656, 110.132568359375))
s.FilletByRadius(radius=20.0, curve1=g[idc4], nearPoint1=(-18.9899597167969, 113.835289001465), curve2=g[idc5],
				 nearPoint2=(-10.4230499267578, 101.01815032959))
s.FilletByRadius(radius=20.0, curve1=g[idc5], nearPoint1=(-10.9941864013672, 63.1363372802734), curve2=g[idc6],
				 nearPoint2=(-14.4209442138672, 52.8826217651367))

tensilePart = model.Part(name='TensileMatrix', dimensionality=THREE_D, type=DEFORMABLE_BODY)
tensilePart = model.parts['TensileMatrix']
tensilePart.BaseSolidExtrude(sketch=s, depth=8.0)

del model.sketches['__profile__']

# creating datum planes
tp1 = tensilePart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=30.0)
idtp1 = tp1.id
tp2 = tensilePart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=56.5)
idtp2 = tp2.id
tp3 = tensilePart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=106.5)
idtp3 = tp3.id
tp4 = tensilePart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=133.3)
idtp4 = tp4.id

# creating partitions

c = tensilePart.cells
d = tensilePart.datums
pickedCells = c.findAt((0.0, 15.0, 4.0))
tensilePart.PartitionCellByDatumPlane(datumPlane=d[idtp1], cells=pickedCells)
pickedCells = c.findAt((0.0, 50.0, 4.0))
tensilePart.PartitionCellByDatumPlane(datumPlane=d[idtp2], cells=pickedCells)
pickedCells = c.findAt((0.0, 70.0, 4.0))
tensilePart.PartitionCellByDatumPlane(datumPlane=d[idtp3], cells=pickedCells)
pickedCells = c.findAt((0.0, 120.0, 4.0))
tensilePart.PartitionCellByDatumPlane(datumPlane=d[idtp4], cells=pickedCells)

# creating particles and assembly

assembly = model.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name='TensileMatrix-1', part=tensilePart, dependent=ON)
assembly.translate(instanceList=('TensileMatrix-1',), vector=(0.0, -81.65, -4.0))

# creating martensite Islands

i = 0
while i < noOfParticle:
	s = model.ConstrainedSketch(name='__profile__', sheetSize=5)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.sketchOptions.setValues(decimalPlaces=4)
	centerConstructLine = s.ConstructionLine(point1=(0.0, -3), point2=(0.0, 3))
	centerConstructLineID = centerConstructLine.id
	s.FixedConstraint(entity=g[centerConstructLineID])
	s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, r[i]), point2=(0.0, -r[i]), direction=CLOCKWISE)
	s.Line(point1=(0.0, -r[i]), point2=(0.0, r[i]))
	p = model.Part(name='Particle-{}'.format(i + 1), dimensionality=THREE_D, type=DEFORMABLE_BODY)
	p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
	p = model.parts['Particle-{}'.format(i + 1)]
	del model.sketches['__profile__']
	XYP = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
	XYPID = XYP.id
	XZP = p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
	XZPID = XZP.id
	YZP = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
	YZPID = YZP.id
	c = p.cells.findAt(((0.0, 0.0, 0.0),))
	d = p.datums
	p.PartitionCellByDatumPlane(datumPlane=d[XYPID], cells=c)
	c = p.cells.getByBoundingBox(-r[i], -r[i], -r[i], r[i], r[i], r[i])
	p.PartitionCellByDatumPlane(datumPlane=d[XZPID], cells=c)
	c = p.cells.getByBoundingBox(-r[i], -r[i], -r[i], r[i], r[i], r[i])
	p.PartitionCellByDatumPlane(datumPlane=d[YZPID], cells=c)
	print('Created Particle-{}'.format(i + 1))
	# creating assembly
	assembly.Instance(name='Particle-{}'.format(i + 1), part=p, dependent=ON)
	assembly.translate(instanceList=('Particle-{}'.format(i + 1),), vector=(x[i], y[i], z[i]))
	if i == 0:
		selectedInstances = (assembly.instances['Particle-{}'.format(i + 1)],)
	else:
		sphereInstances = assembly.instances['Particle-{}'.format(i + 1)]
		selectedInstances = selectedInstances + (sphereInstances,)
	i += 1

selectedInstances = selectedInstances + (assembly.instances['TensileMatrix-1'],)
# particles placed on their coordinates

# merging instances
assembly.InstanceFromBooleanMerge(name='Sample', instances=selectedInstances, keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
print('merged particles and matrix')
assembly.regenerate()

i = 0
while i < noOfParticle:
	del model.parts['Particle-{}'.format(i+1)]
	i += 1

# assigning section

# creating ferrite section
sample = model.parts['Sample']
ferriteCell = sample.cells.findAt(((-11.0,-25.0,0.0),))
ferriteCell += sample.cells.findAt(((0.0,-70.0,0.0),)) + sample.cells.findAt(((0.0,-40.0,0.0),))
ferriteCell += sample.cells.findAt(((0.0,70.0,0.0),)) + sample.cells.findAt(((0.0,40.0,0.0),))
sample.Set(cells=ferriteCell,name='Set-Ferrite')
region = sample.sets['Set-Ferrite']
sample.SectionAssignment(region=region,sectionName='FerriteSection',offset=0.0,offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)
print('Ferrite Section Assigned')

mpl = []  #list of coordinates of martensite region
i = 0
while i < noOfParticle:
	mpl.append((x[i] + r[i] / 2, y[i] + r[i] / 2, z[i] + r[i] / 2))
	mpl.append((x[i] - r[i] / 2, y[i] + r[i] / 2, z[i] + r[i] / 2))
	mpl.append((x[i] + r[i] / 2, y[i] - r[i] / 2, z[i] + r[i] / 2))
	mpl.append((x[i] - r[i] / 2, y[i] - r[i] / 2, z[i] + r[i] / 2))
	mpl.append((x[i] + r[i] / 2, y[i] + r[i] / 2, z[i] - r[i] / 2))
	mpl.append((x[i] - r[i] / 2, y[i] + r[i] / 2, z[i] - r[i] / 2))
	mpl.append((x[i] + r[i] / 2, y[i] - r[i] / 2, z[i] - r[i] / 2))
	mpl.append((x[i] - r[i] / 2, y[i] - r[i] / 2, z[i] - r[i] / 2))
	i = i + 1
mp = tuple(mpl)  # tuple list of all coordinates of martensite cells
i = 0
martensiteCell= sample.cells.findAt((mp[i],))
while i < len(mp):
	martensiteCell = martensiteCell + sample.cells.findAt((mp[i],))
	i = i +1
sample.Set(cells=martensiteCell, name='Set-Martensite')
sample.SectionAssignment(region = sample.sets['Set-Martensite'],sectionName='MartensiteSection',offset=0.0,offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)

print("assigned section")

assembly.regenerate()

# creating step

model.StaticStep(name='Step-1',previous='Initial',maxNumInc=1000,initialInc=0.001,minInc=1e-8, maxInc = 0.05, nlgeom=ON)
print('Created Step')

# creating boundary conditions
c = assembly.instances['Sample-1'].cells
cells = c.findAt(((0.0,-70.0,0.0),))
region = assembly.Set(cells=cells,name='BCSet1')
model.EncastreBC(name='BC-1',createStepName='Initial', region=region, localCsys=None)
cells = c.findAt(((0.0,70.0,0.0),))
region = assembly.Set(cells=cells,name='BCSet2')
model.DisplacementBC(name='BC-2',createStepName='Step-1',region=region, u1=0.0,u2=12.0,u3=0.0,ur1=0.0,ur2=0.0,ur3=0.0,
					 amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='',localCsys=None)
