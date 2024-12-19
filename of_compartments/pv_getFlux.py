# trace generated using paraview version 5.10.1
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import argparse
import pathlib
import os
import subprocess

parser = argparse.ArgumentParser(description="Process a bioreactor sim.")
parser.add_argument(
    "--case_dir",
    "-c",
    type=pathlib.Path,
    help="the path to the case folder",
)
parser.add_argument(
    "--phase",
    "-p",
    choices=["water", "air"],
    help="the path to the case folder",
)

args = parser.parse_args()
case_dir = args.case_dir
phase = args.phase
output_dir = case_dir


case_name = os.path.basename(case_dir)
foam_path = (
    subprocess.check_output(f"ls {os.path.join(case_dir, '*.foam')}", shell=True)
    .strip()
    .decode("utf-8")
)

# phase = "air"
# output_dir = "."


def getCellCenters(source):
    cellCenters1 = CellCenters(registrationName="CellCenters1", Input=source)
    cellCenters1.VertexCells = 1

    calculator_centers = Calculator(
        registrationName="Calculator_center", Input=cellCenters1
    )
    calculator_centers.ResultArrayName = "Center"
    calculator_centers.Function = "coords"
    pointDatatoCellData1 = PointDatatoCellData(
        registrationName="PointDatatoCellData1", Input=calculator_centers
    )
    pointDatatoCellData1.PointDataArraytoprocess = ["Center"]
    passArrays1 = PassArrays(registrationName="PassArrays1", Input=pointDatatoCellData1)
    passArrays1.CellDataArrays = ["Center"]
    return passArrays1


def getSurfaceNormals(source):
    centers = getCellCenters(source)
    generateSurfaceNormals1 = GenerateSurfaceNormals(
        registrationName="GenerateSurfaceNormals1", Input=source
    )
    generateSurfaceNormals1.FlipNormals = 1
    generateSurfaceNormals1.ComputeCellNormals = 1
    programmableFilter1 = ProgrammableFilter(
        registrationName="ProgrammableFilter1", Input=[generateSurfaceNormals1, centers]
    )
    programmableFilter1.Script = """import numpy as np
normals = inputs[0].CellData["Normals"]
centers = inputs[1].CellData["Center"]
negation = -1 * (dot(normals, centers) < 0) + 1 * (dot(normals, centers) >= 0)
output.CellData.append(negation * normals, "normals_corrected")
    """
    return programmableFilter1


def getFlux(source, phase):
    calculator1 = Calculator(registrationName="Calculator1", Input=source)
    calculator1.AttributeType = "Cell Data"
    calculator1.ResultArrayName = "flux"
    if phase == "water":
        calculator1.Function = (
            'dot(normals_corrected, "UMean.water")*(1-"alphaMean.air")'
        )
    else:
        calculator1.Function = 'dot(normals_corrected, "UMean.air")*"alphaMean.air"'
    return calculator1


def getFluxPerBoundary(source):
    programmableFilter1 = ProgrammableFilter(
        registrationName="ProgrammableFilter1", Input=source
    )
    programmableFilter1.OutputDataSetType = "vtkTable"
    programmableFilter1.Script = """import numpy as np
from vtk.util import numpy_support

input0 = inputs[0]

dataArray = input0.CellData["flux"]


nBlocks = len(dataArray.Arrays)

flux_pos = np.zeros(nBlocks)
flux_neg = np.zeros(nBlocks)
names = np.empty(nBlocks, dtype=object)
for i in range(nBlocks):
    names[i] = ""
    a = dataArray.Arrays[i]
    if type(a) == vtk.numpy_interface.dataset_adapter.VTKNoneArray:
        continue
    areas = input0.CellData["Area"].Arrays[i]
    names[i] = inputs[0].GetBlock(0).GetBlock(0).GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
    pos = a > 0
    flux_pos[i] = np.sum(np.dot(a * pos,  areas))

    neg = a < 0
    flux_neg[i] = np.sum(np.dot(a * neg, areas))


vtkarr = numpy_support.numpy_to_vtk( flux_pos, deep=True, array_type=vtk.VTK_DOUBLE )
vtkarr.SetNumberOfComponents( 1 )
vtkarr.SetNumberOfTuples( len(flux_pos) )
vtkarr.SetName("flux_pos")
output.GetRowData().AddArray(vtkarr)

vtkarr = numpy_support.numpy_to_vtk( flux_neg, deep=True, array_type=vtk.VTK_DOUBLE )
vtkarr.SetNumberOfComponents( 1 )
vtkarr.SetNumberOfTuples( len(flux_neg) )
vtkarr.SetName("flux_neg")
output.GetRowData().AddArray(vtkarr)

arr = vtk.vtkStringArray()
arr.SetName("name")
arr.SetNumberOfComponents(1)
for name in names:
    arr.InsertNextValue(name)
output.AddColumn(arr)

"""


case = OpenDataFile(foam_path)
case.CaseType = "Reconstructed Case"
case.Readzones = 1
case.Copydatatocellzones = 1
# case = GetActiveSource()

materialLibrary1 = GetMaterialLibrary()
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

extractBlock1 = ExtractBlock(registrationName="ExtractBlock1", Input=case)
extractBlock1.Selectors = ["/Root/zones/faceZones"]

normals = getSurfaceNormals(extractBlock1)
extractBlock1 = AppendAttributes(
    registrationName="AppendAttributes1", Input=[normals, extractBlock1]
)

flux = getFlux(extractBlock1, phase)

cellSize1 = CellSize(registrationName="CellSize1", Input=flux)
cellSize1.ComputeVertexCount = 0
cellSize1.ComputeLength = 0
cellSize1.ComputeVolume = 0

result = getFluxPerBoundary(cellSize1)

spreadSheetView1 = CreateView("SpreadSheetView")
spreadSheetView1.ColumnToSort = ""
spreadSheetView1.BlockSize = 1024
resultDisplay = Show(result, spreadSheetView1, "SpreadSheetRepresentation")

layout1 = GetLayoutByName("Layout #1")

AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

spreadSheetView1.Update()

ExportView(
    os.path.join(output_dir, f"flux_data_{phase}.csv"),
    view=spreadSheetView1,
)
