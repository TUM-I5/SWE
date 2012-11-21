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

out_00_nc = NetCDFReader( FileName=[str(files[0])] )

#AnimationScene1 = GetAnimationScene()
out_00_nc.Dimensions = '(y, x)'

#AnimationScene1.EndTime = 15.003315925598145
#AnimationScene1.PlayMode = 'Snap To TimeSteps'

#a2DRenderView1 = Create2DRenderView()
#a2DRenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
#a2DRenderView1.RemoteRenderThreshold = 3.0
#a2DRenderView1.ViewTime = 0.0
#a2DRenderView1.LODResolution = 50.0
#a2DRenderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
#a2DRenderView1.LODThreshold = 5.0

#AnimationScene1.ViewModules = a2DRenderView1

#DataRepresentation1 = Show()
#DataRepresentation1.ColorArrayName = 'b'
#DataRepresentation1.UseXYPlane = 0

Calculator1 = Calculator()

#a1_b_PVLookupTable = GetLookupTableForArray( "b", 1, NanColor=[0.25, 0.0, 0.0], RGBPoints=[-260.0, 0.23, 0.299, 0.754, -255.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

#a1_b_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

#a2DRenderView1.CameraPosition = [500.0, 500.0, 2724.048059039453]
#a2DRenderView1.CameraFocalPoint = [500.0, 500.0, 0.0]
#a2DRenderView1.CameraClippingRange = [2696.8075784490584, 2764.9087799250456]
#a2DRenderView1.CameraParallelScale = 705.0355174539663

#DataRepresentation1.LookupTable = a1_b_PVLookupTable

Calculator1.AttributeMode = 'point_data'

Calculator1.Function = 'min(h,h+b)'

DataRepresentation2 = Show()
#DataRepresentation2.ColorArrayName = 'b'
DataRepresentation2.Representation = 'Slice'
#DataRepresentation2.UseXYPlane = 0

#DataRepresentation1.Visibility = 0

a1_Result_PVLookupTable = GetLookupTableForArray( "Result", 1, RGBPoints=[-1.0, 0.23, 0.299, 0.754, 2.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, LockScalarRange=1 )

#a1_Result_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

DataRepresentation2.ColorArrayName = 'Result'
DataRepresentation2.LookupTable = a1_Result_PVLookupTable

GetRenderView.ResetCamera()
Render()
