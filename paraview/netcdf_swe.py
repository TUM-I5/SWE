#!/usr/bin/python

try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# Select the files with PyQt4
from PyQt4 import QtGui, QtCore
def nullMessageOutput(type, msg):
	pass
QtCore.qInstallMsgHandler(nullMessageOutput)
files = QtGui.QFileDialog.getOpenFileNames(None, 'Select SWE output files ...', QtCore.QString(), 'NetCDF (*.nc)')

# List of all sources we load
sources = []

for file in files:
	# Create NetCDF reader
	reader = NetCDFReader( FileName=[str(file)] )
	reader.Dimensions = '(y, x)'

	sources.append(reader)

# Group Data if mor then one
if len(sources) > 1:
	
	group = GroupDatasets( Input=sources )
else:
	group = sources[0]

# Create Calculator, to show "min(h,h+b)"
calc = Calculator( Input=group )
calc.AttributeMode = 'point_data'
calc.Function = 'min(h,h+b)'

# Specify the data representation
representation = Show()
representation.Representation = 'Surface'

# Create nice color table form -1 to 2
table = GetLookupTableForArray( "Result", 1, RGBPoints=[-1.0, 0.23, 0.299, 0.754, 2.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, LockScalarRange=1 )

# Color result
representation.ColorArrayName = 'Result'
representation.LookupTable = table

GetRenderView().ResetCamera()
Render()
