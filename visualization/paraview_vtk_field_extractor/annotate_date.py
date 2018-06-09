"""
A macro for annotating with YYYY-MM-DD based on time values that are assumed
to be years since 0000-01-01.

Add to ParaView via Macros -> Add new macro..

Apply by selecting the imported pdv file in the pipeline browser, then running
Macros -> annotate_date

Xylar Asay-Davis
24-Aug-2017
"""

# import the simple module from the paraview
import paraview.simple

# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
source = paraview.simple.GetActiveSource()

# create a new 'Programmable Filter'
programmableFilter1 = \
    paraview.simple.ProgrammableFilter(Input=source)
programmableFilter1.OutputDataSetType = 'vtkTable'
programmableFilter1.Script = \
    'from paraview.simple import GetAnimationScene\n' \
    'import numpy\n' \
    '\n' \
    'daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]\n' \
    'cumSumDays = numpy.append([0], numpy.cumsum(daysInMonth))\n' \
    '\n' \
    't = GetAnimationScene().TimeKeeper.Time\n' \
    'days = int(365.*t+0.5)\n' \
    'year = int(days/365.)\n' \
    'days -= 365*year\n' \
    'index = numpy.nonzero(days < cumSumDays)[0][0]\n' \
    'month = index\n' \
    'days = days - cumSumDays[index-1]\n' \
    'day = days+1\n' \
    '\n' \
    'datetime = \'{:04d}-{:02d}-{:02d}\'.format(year, month, day)\n' \
    '\n' \
    'outputarray = vtk.vtkStringArray()\n' \
    'outputarray.SetName("datetime")\n' \
    'outputarray.SetNumberOfTuples(1)\n' \
    'outputarray.SetValue(0, "{}".format(datetime))\n' \
    'output.RowData.AddArray(outputarray)'
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# create a new 'Python Annotation'
pythonAnnotation1 = paraview.simple.PythonAnnotation(Input=programmableFilter1)
pythonAnnotation1.ArrayAssociation = 'Row Data'
pythonAnnotation1.Expression = \
    '"{}".format(input.RowData["datetime"].GetValue(0))'

# find view
renderView1 = paraview.simple.GetActiveView()
# show data in view
pythonAnnotation1Display = paraview.simple.Show(pythonAnnotation1, renderView1)
pythonAnnotation1Display.FontSize = 12
# update the view to ensure updated data information
renderView1.Update()

