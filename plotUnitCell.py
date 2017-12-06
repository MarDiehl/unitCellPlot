#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-^

import numpy as np
import os
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options ', description = """
Creates unit cell including slip planes and directions for visualization

""", version = scriptID)

parser.add_option('-r', '--radians',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees [%default]')
parser.add_option('-e', '--eulers',
                  dest='eulers',
                  type = 'float', nargs = 3, metavar = ' '.join(['float']*3),
                  help = 'crystallographic orientation as euler angles')

parser.set_defaults(eulers  = (0.,0.,0.),                                                   # no rotation about 1,1,1
                    radians = False)

(options, filenames) = parser.parse_args()

eulers=np.array(options.eulers)
if not options.radians: eulers*=(np.pi/180.)

ori = damask.Orientation(matrix=np.dot(np.array([[0.0,-1.0,0.0],[-1.0,0.0,0.0],[0.0,0.0,-1.0]]),
                                       damask.Orientation(Eulers=eulers).asMatrix()))
ori=damask.Orientation()
coordinates=[]
for x in np.linspace(-0.5,0.5,3):
  for y in np.linspace(-0.5,0.5,3):
    for z in np.linspace(-0.5,0.5,3):
      coordinates.append([x,y,z])
coordinates=np.array(coordinates)
for i,c in enumerate(coordinates):
    coordinates[i,:] = np.dot(ori.asMatrix(),c)

Miller110=np.array(
  [[ 0, 1, 1],
   [ 0,-1, 1],
   [ 1, 0, 1],
   [-1, 0, 1],
   [ 1, 1, 0],
   [-1, 1, 0]],dtype='float')  
for i,p in enumerate(Miller110):
    Miller110[i,:] = np.dot(ori.asMatrix(),p)

Miller111=np.array(
  [[ 1,-1, 1],
   [-1,-1, 1],
   [ 1, 1, 1],
   [-1, 1, 1]],dtype='float')
for i,d in enumerate(Miller111):
    Miller111[i,:] = np.dot(ori.asMatrix(),d)

Miller112=np.array(
  [[ 2, 1, 1],
   [-2, 1, 1],
   [ 2,-1, 1],
   [ 2, 1,-1],
   [ 1, 2, 1],
   [-1, 2, 1],
   [ 1,-2, 1],
   [ 1, 2,-1],
   [ 1, 1, 2],
   [-1, 1, 2],
   [ 1,-1, 2],
   [ 1, 1,-2]],dtype='float')
for i,p in enumerate(Miller112):
    Miller112[i,:] = np.dot(ori.asMatrix(),p)

XDMF_file=open('slipSystems.xdmf','w')
XDMF_file.write(
'<?xml version="1.0" ?>\n'
'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
'<Xdmf Version="2.0" xmlns:xi="[http://www.w3.org/2001/XInclude]">\n'
'  <Domain>\n'
'    <Grid Name="Unit cell">\n'
'      <Topology TopologyType="Hexahedron" NumberOfElements="1">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="8 1">\n'
'          0 6 8 2 18 24 26 20\n'
'        </DataItem>\n'
'      </Topology>\n'
'      <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Geometry>\n'
'    </Grid>\n'
'    <Grid Name="{110}-planes (bcc)">\n'
'      <Topology TopologyType="Quadrilateral" NumberOfElements="12">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="4 12">\n'
'          2 6 24 20\n'
'          2 6 24 20\n'
'          0 18 26 8\n'
'          0 18 26 8\n'
'          2 8 24 18\n'
'          2 8 24 18\n'
'          6 26 20 0\n'
'          6 26 20 0\n'
'          6 8 20 18\n'
'          6 8 20 18\n'
'          26 2 0 24\n'
'          26 2 0 24\n'
'        </DataItem>\n'
'      </Topology>\n'
'      <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Geometry>\n'
'      <Attribute Name="System Index" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="12" DataType="Int">\n'
'          1 2 3 4 5 6 7 8 9 10 11 12\n'
'        </DataItem>\n'
'      </Attribute>\n'
'      <Attribute Name="Plane Normal" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="3 12" DataType="Float">\n')
for p in Miller110:
  for m in range(2):
    XDMF_file.write('          {} {} {}\n'.format(p[0],p[1],p[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Attribute>\n'
'    </Grid>\n'
'    <Grid Name="{112}-planes (bcc)">\n'
'      <Topology TopologyType="Quadrilateral" NumberOfElements="12">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="4 12">\n'
'          18 15 08 11\n'
'          00 15 26 11\n'
'          02 17 24 09\n'
'          20 17 06 09\n'
'          06 21 20 05\n'
'          03 24 23 02\n'
'          00 21 26 05\n'
'          08 03 18 23\n'
'          02 07 24 19\n'
'          06 25 20 01\n'
'          08 25 18 01\n'
'          07 26 19 00\n'
'        </DataItem>\n'
'      </Topology>\n'
'       <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'<!--      \n'
'          0.0 0.0 0.0\n'
'          0.0 0.0 0.5\n'
'          0.0 0.0 1.0\n'
'        3 \n'
'          0.0 0.5 0.0\n'
'          0.0 0.5 0.5\n'
'          0.0 0.5 1.0\n'
'        6\n'
'          0.0 1.0 0.0\n'
'          0.0 1.0 0.5\n'
'          0.0 1.0 1.0\n'
'        9\n'
'          0.5 0.0 0.0\n'
'          0.5 0.0 0.5\n'
'          0.5 0.0 1.0\n'
'       12\n'
'          0.5 0.5 0.0\n'
'          0.5 0.5 0.5\n'
'          0.5 0.5 1.0\n'
'       15\n'
'          0.5 1.0 0.0\n'
'          0.5 1.0 0.5\n'
'          0.5 1.0 1.0\n'
'       18\n'
'          1.0 0.0 0.0\n'
'          1.0 0.0 0.5\n'
'          1.0 0.0 1.0\n'
'       21\n'
'          1.0 0.5 0.0\n'
'          1.0 0.5 0.5\n'
'          1.0 0.5 1.0\n'
'       24   \n'
'          1.0 1.0 0.0\n'
'          1.0 1.0 0.5\n'
'          1.0 1.0 1.0\n'
'-->\n'
'      </Geometry>\n'
'      <Attribute Name="System Index" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="12" DataType="Int">\n'
'          13 14 15 16 17 18 19 20 21 22 23 24\n'
'        </DataItem>\n'
'      </Attribute>\n'
'      <Attribute Name="Plane Normal" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="3 12" DataType="Float">\n')
for p in Miller112:
  XDMF_file.write('          {} {} {}\n'.format(p[0],p[1],p[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Attribute>\n'
'    </Grid>\n'
'    <Grid Name="{111}-planes (fcc)">\n'
'      <Topology TopologyType="Triangle" NumberOfElements="12">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="3 12">\n'
'          02 06 18\n'
'          02 06 18\n'
'          02 06 18\n'
'          06 18 26\n'
'          06 18 26\n'
'          06 18 26\n'
'          02 26 06\n'
'          02 26 06\n'
'          02 26 06\n'
'          18 26 02\n'
'          18 26 02\n'
'          18 26 02\n'
'        </DataItem>\n'
'      </Topology>\n'
'      <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Geometry>\n'
'      <Attribute Name="System Index" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="12" DataType="Int">\n'
'          1 2 3 4 5 6 7 8 9 10 11 12\n'
'        </DataItem>\n'
'      </Attribute>\n'
'      <Attribute Name="Plane Normal" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="3 12" DataType="Float">\n')
for i in range(3):
  XDMF_file.write('          {} {} {}\n'.format(Miller111[2,0],Miller111[2,1],Miller111[2,2]))
for i in range(3):
  XDMF_file.write('          {} {} {}\n'.format(Miller111[1,0],Miller111[1,1],Miller111[1,2]))
for i in range(3):
  XDMF_file.write('          {} {} {}\n'.format(Miller111[3,0],Miller111[3,1],Miller111[3,2]))
for i in range(3):
  XDMF_file.write('          {} {} {}\n'.format(Miller111[0,0],Miller111[0,1],Miller111[0,2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Attribute>\n'
'    </Grid>\n'
'    <Grid Name="[111]-directions (bcc)">\n'
'      <Topology TopologyType="Polyline" NumberOfElements="24" NodesPerElement="2">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="2 24">\n'
'          6 20\n'
'          6 20\n'
'          6 20\n'
'          6 20\n'
'          6 20\n'
'          6 20\n'
'          2 24\n'
'          2 24\n'
'          2 24\n'
'          2 24\n'
'          2 24\n'
'          2 24\n'
'          0 26\n'
'          0 26\n'
'          0 26\n'
'          0 26\n'
'          0 26\n'
'          0 26\n'
'          8 18\n'
'          8 18\n'
'          8 18\n'
'          8 18\n'
'          8 18\n'
'          8 18\n'
'        </DataItem>\n'
'      </Topology>\n'
'      <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Geometry>\n'
'      <Attribute Name="System Index" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="24" DataType="Int">\n'
'          1 8 10 16 17 22 2 6 12 15 18 21 3 7 11 14 19 24 4 5 9 13 20 23\n'
'        </DataItem>\n'
'      </Attribute>\n'
'      <Attribute Name="Slip Direction" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="3 24" DataType="Float">\n')
for d in Miller111:
  for m in range(6):
    XDMF_file.write('          {} {} {}\n'.format(d[0],d[1],d[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Attribute>\n'
'    </Grid>\n'
'    <Grid Name="[110]-directions (fcc)">\n'
'      <Topology TopologyType="Polyline" NumberOfElements="12" NodesPerElement="2">\n'
'        <DataItem Format="XML" DataType="Int" Dimensions="2 12">\n'
'          02 06\n'
'          18 02\n'
'          06 18\n'
'          18 26\n'
'          26 06\n'
'          18 06\n'
'          06 02\n'
'          26 06\n'
'          26 02\n'
'          26 18\n'
'          18 02\n'
'          26 02\n'
'        </DataItem>\n'
'      </Topology>\n'
'      <Geometry GeometryType="XYZ">\n'
'        <DataItem Name="Points" Format="XML" Dimensions="27 3" DataType="Float">\n')
for c in coordinates:
  XDMF_file.write('          {} {} {}\n'.format(c[0],c[1],c[2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Geometry>\n'
'      <Attribute Name="System Index" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="12" DataType="Int">\n'
'          1 2 3 4 5 6 7 8 9 10 11 12\n'
'        </DataItem>\n'
'      </Attribute>\n'
'      <Attribute Name="Slip Direction" Center="Cell">\n'
'        <DataItem Format="XML" Dimensions="3 12" DataType="Float">\n')
for i in [1,3,5,0,2,5,1,2,4,0,3,4]:
  XDMF_file.write('          {} {} {}\n'.format(Miller110[i,0],Miller110[i,1],Miller110[i,2]))
XDMF_file.write(\
'        </DataItem>\n'
'      </Attribute>\n'
'    </Grid>\n'
'  </Domain>\n'
'</Xdmf>\n'
)
