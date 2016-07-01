''' This script to be run through VisIt. Python 2.7. e.g.:
visit -nowin -cli -s this_script.py'''

from __future__ import division # make float division the default
import sys
# visit_utils has several nice helper functions, among them "query" which just
#   returns the numerical output of Query.
from visit_utils import *


path_to_dumps = "C:\Users\Arviragus\Dropbox\MacrophyteProject\Visit3D\sample_viz_IB3d\dumps.visit"


### This bit may not be necessary - it just sets a default view ###
# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (0, 1, 0, 1)
View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 0
SetView2D(View2DAtts)
# End spontaneous state

### Open data, plot x-magnitude, add slice ###
OpenDatabase("localhost:"+path_to_dumps, 0) #opens data to first time step
last_time = TimeSliderGetNStates() - 1 #returns the number of time slider states
                                       #   subtract 1 for base 0
SetTimeSliderState(last_time)
# Define some new variables
DefineScalarExpression("U_x_abs", "abs(<U_x>)") # abs of U in x direction
DefineScalarExpression("U_y_abs", "abs(<U_y>)")
DefineScalarExpression("U_z_abs", "abs(<U_z>)")
AddPlot("Pseudocolor", "U_x_abs")
AddOperator("Slice", 0) # the 0 here means "apply operator only to this plot"
SetActivePlots(1) # there are two plots now, since levels is auto added

### Set slice attributes ###
SliceAtts = SliceAttributes()
# Specify the origin via intercept
SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
SliceAtts.originPoint = (0, 0, 0)
SliceAtts.originIntercept = -0.45 # This is where on the Z Axis we are!
SliceAtts.originPercent = 0
SliceAtts.originZone = 0
SliceAtts.originNode = 0
SliceAtts.normal = (0, 0, 1)
# Set plane orthogonal to the Z Axis
SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
SliceAtts.upAxis = (0, 1, 0)
SliceAtts.project2d = 0 # Turn off 2D projection!
SliceAtts.interactive = 1
SliceAtts.flip = 0
SliceAtts.originZoneDomain = 0
SliceAtts.originNodeDomain = 0
SliceAtts.meshName = "amr_mesh"
SliceAtts.theta = 0
SliceAtts.phi = 90
SetOperatorOptions(SliceAtts, 0) # Apply Slice settings
DrawPlots() # Draw the plot

### Get the average value and store it ###
x_abs_avgs = []
x_abs_avgs.append(query("Average Value"))

### Move the slice and get another average ###
SliceAtts.originIntercept = -0.4
SetOperatorOptions(SliceAtts, 0) # Operator number 0
DrawPlots()
x_abs_avgs.append(query("Average Value"))

print x_abs_avgs