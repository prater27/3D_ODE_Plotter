#include <vtkInteractorStyleUnicam.h>
#include <vtkInteractorStyleJoystickCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <chrono>
#include <vtkMultiThreader.h>

#include <iostream>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSphereSource.h>
#include <vtkElevationFilter.h>
#include <vtkVectorText.h>
#include <vtkCommand.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkQuaternion.h>
#include <vtkTransform.h>
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

//CoordinateSystem
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkNamedColors.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>

//Spline
#include <vtkKochanekSpline.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkProperty.h>
#include <vtkGlyph3DMapper.h>
#include <vtkNamedColors.h>

//Plane
#include <vtkPlaneSource.h>
#include <vtkLegendBoxActor.h>
#include <vtkTransformPolyDataFilter.h>

//Animation
#include <vtkProgrammableFilter.h>
#include <vtkCommand.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkDataSetAttributes.h>

#include <string>
#include <array>

#include "vtkPolyDataAlgorithm.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkVertex.h>

#include <cmath>

#include <vtkHedgeHog.h>

#include <vtkGlyph3D.h>

#include <vtkArrowSource.h>
#include <vtkLookupTable.h>


int timePointsProcessedCounter=1;
bool animStart=false;

void KeypressCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
    vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

    std::string ret = "Return";
    std::string pushed(iren->GetKeySym());
  // Handle an arrow key
  if(pushed == ret && animStart==false)
  {
        //startTime = std::chrono::steady_clock::now();
    //    iren->CreateOneShotTimer(1);
        iren->CreateRepeatingTimer(50);
        animStart=true;
  }
}

struct MeshParameters{
    int nx, ny, nz, nt;
    double xMax, yMax, zMax, xMin,yMin, zMin, tMin, tMax;
};

void generatePoints(const MeshParameters& meshParam,
                    const vtkSmartPointer<vtkPoints>& points,
                    const vtkSmartPointer<vtkCellArray>& verts,
                    const vtkSmartPointer<vtkPolyData>& polyDataPoints,
                    std::vector<double>& timeVector){
    for(int x=0; x<meshParam.nx;x++){
        for(int y=0; y<meshParam.ny;y++){
            for(int z=0; z<meshParam.nz;z++){
                double point[3] = {meshParam.xMin + x*((meshParam.xMax-meshParam.xMin)/meshParam.nx),
                                   meshParam.yMin + y*((meshParam.yMax-meshParam.yMin)/meshParam.ny),
                                   meshParam.zMin + z*(meshParam.zMax-meshParam.zMin)/meshParam.nz};
                vtkIdType id = points->InsertNextPoint(point);
                verts->InsertNextCell(1, &id);            }
        }
    }
    polyDataPoints->SetPoints(points);
    polyDataPoints->SetVerts(verts);

    for(int j=0; j<meshParam.nt; j++){
        double  taux = meshParam.tMin+j*((meshParam.tMax-meshParam.tMin)/meshParam.nt);
        timeVector.push_back(taux);
    }
}

void evaluateFunction(const double t, const double (&position)[3], double (&result)[3]){
    result[0] = sin((1+1*t/10)*(position[1]-position[0]));
    result[1] = cos((1+1*t/10)*position[0]*(50-position[2])-position[1]);
    result[2] = (1+1*t/10)*position[0]*position[1]-75*position[2];
}

void generateVectorFieldAtT(const vtkSmartPointer<vtkPolyData>& polyDataPoints,
                            const vtkSmartPointer<vtkDoubleArray>& vectorField, const double& t,
                            const vtkSmartPointer<vtkLookupTable>& colorLookupTable){

       vtkSmartPointer<vtkDoubleArray> magnitudes = vtkSmartPointer<vtkDoubleArray>::New();

        vectorField->Reset();
        polyDataPoints->GetPointData()->Reset();

        for (int i=0; i<polyDataPoints->GetNumberOfPoints(); i++){
            double result[3];
            double position[3];

            polyDataPoints->GetPoint(i, position);
            evaluateFunction(t, position, result);
            vectorField->InsertTuple3(i, result[0], result[1], result[2]);
            double scalarResult = sqrt((pow(result[0],2)+pow(result[1],2)+pow(result[2],2)));
            magnitudes->InsertValue(i,scalarResult);
        }
        double scalarMagnitudes[2];
        polyDataPoints->GetPointData()->SetVectors(vectorField);
        polyDataPoints->GetPointData()->SetScalars(magnitudes);
        polyDataPoints->GetScalarRange(scalarMagnitudes);

// Create the color map
        colorLookupTable->ResetAnnotations();
        colorLookupTable->SetTableRange(scalarMagnitudes[0], scalarMagnitudes[1]);
        colorLookupTable->Build();

        polyDataPoints->GetPointData()->GetScalars()->SetName("Magnitudes");
        polyDataPoints->GetPointData()->Update();
}

struct CallBackParameters{
    vtkSmartPointer<vtkPolyData> polyDataPoints;
    vtkSmartPointer<vtkDoubleArray> vectorField;
    std::vector<double> timeVector;
    vtkSmartPointer<vtkPolyDataMapper> polyDataMapper;

    CallBackParameters(const vtkSmartPointer<vtkPolyData>& polyDataPoints_, const vtkSmartPointer<vtkDoubleArray>& vectorField_,
                       const std::vector<double>& timeVector_, const vtkSmartPointer<vtkPolyDataMapper>& polyDataMapper_):
        polyDataPoints(polyDataPoints_),
        vectorField(vectorField_),
        timeVector(timeVector_),
        polyDataMapper(polyDataMapper_)
        {}
};

void TimerCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData) )
{
    CallBackParameters* params = static_cast<CallBackParameters*>(clientData);
    vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

if(timePointsProcessedCounter<params->timeVector.size()){
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    generateVectorFieldAtT(params->polyDataPoints, params->vectorField,
                           params->timeVector[timePointsProcessedCounter], colorLookupTable);

    //params->polyDataPoints->GetPointData()->SetVectors(params->vectorField);
   // params->polyDataPoints->GetPointData()->SetActiveScalars("Magnitudes");
   // params->polyDataMapper->SetLookupTable(colorLookupTable);
  //  params->polyDataMapper->SetScalarRange(params->polyDataPoints->GetScalarRange());

    iren->Render();
    timePointsProcessedCounter++;
 //   iren->CreateOneShotTimer(1);
}
}



int main()
{

//Generation of the base mesh
vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
vtkSmartPointer<vtkPolyData> polyDataPoints = vtkSmartPointer<vtkPolyData>::New();
MeshParameters meshParams;
std::vector<double> timeVector;

meshParams.nx = 10;
meshParams.ny = 10;
meshParams.nz = 10;
meshParams.nt = 500;

meshParams.xMin = -100;
meshParams.yMin = -100;
meshParams.zMin = -100;

meshParams.xMax = 100;
meshParams.yMax = 100;
meshParams.zMax = 100;

meshParams.tMin=0;
meshParams.tMax=100;

//It generates the mesh of spatial and temporal points
generatePoints(meshParams, points, verts, polyDataPoints, timeVector);

//Colors
  vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
  vtkColor3d backgroundColor = colors->GetColor3d("SlateGray");
  vtkColor3d legendBackgroundColor = colors->GetColor3d("Black");

//STATIC ELEMENTS

//Plane
  // Create planes xy
    double planeSize = 100.0;
    vtkSmartPointer<vtkPlaneSource> planeSourceXY = vtkSmartPointer<vtkPlaneSource>::New();
    planeSourceXY->SetCenter(0.0, 0.0, 0.0);
    planeSourceXY->SetPoint1(planeSize/2, 0.0, 0.0);
    planeSourceXY->SetPoint2(0.0, planeSize/2, 0.0);
    planeSourceXY->Update();

    vtkPolyData* planeXY = planeSourceXY->GetOutput();

    // Create a mapper and actor for the plane xy
    vtkSmartPointer<vtkPolyDataMapper> mapperPlaneXY = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPlaneXY->SetInputData(planeXY);

    vtkSmartPointer<vtkActor> actorPlaneXY = vtkSmartPointer<vtkActor>::New();
    actorPlaneXY->SetMapper(mapperPlaneXY);
    actorPlaneXY->GetProperty()->SetOpacity(0.5);


    // Create planes xz
    vtkSmartPointer<vtkPlaneSource> planeSourceXZ = vtkSmartPointer<vtkPlaneSource>::New();
    planeSourceXZ->SetCenter(0.0, 0.0, 0.0);
    planeSourceXZ->SetPoint1(planeSize/2, 0.0, 0.0);
    planeSourceXZ->SetPoint2(0.0, 0.0, planeSize/2);
    planeSourceXZ->Update();

    vtkPolyData* planeXZ = planeSourceXZ->GetOutput();

    // Create a mapper and actor for the plane xz
    vtkSmartPointer<vtkPolyDataMapper> mapperPlaneXZ = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPlaneXZ->SetInputData(planeXZ);

    vtkSmartPointer<vtkActor> actorPlaneXZ = vtkSmartPointer<vtkActor>::New();
    actorPlaneXZ->SetMapper(mapperPlaneXZ);
    actorPlaneXZ->GetProperty()->SetOpacity(0.5);


    // Create planes yz
    vtkSmartPointer<vtkPlaneSource> planeSourceYZ = vtkSmartPointer<vtkPlaneSource>::New();
    planeSourceYZ->SetCenter(0.0, 0.0, 0.0);
    planeSourceYZ->SetPoint1(0.0, planeSize/2, 0.0);
    planeSourceYZ->SetPoint2(0.0, 0.0, planeSize/2);
    planeSourceYZ->Update();

    vtkPolyData* planeYZ = planeSourceYZ->GetOutput();

    // Create a mapper and actor for the plane yz
    vtkSmartPointer<vtkPolyDataMapper> mapperPlaneYZ = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperPlaneYZ->SetInputData(planeYZ);

    vtkSmartPointer<vtkActor> actorPlaneYZ = vtkSmartPointer<vtkActor>::New();
    actorPlaneYZ->SetMapper(mapperPlaneYZ);
    actorPlaneYZ->GetProperty()->SetOpacity(0.5);


//Coordinate System Origin
    //Create actor coordinate system axes
    vtkSmartPointer<vtkAxesActor> actorAxesOrigin = vtkSmartPointer<vtkAxesActor>::New();
    //Change the size of the axes, change label size, change label text, change label color...// axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1,0,0);
    actorAxesOrigin->SetTotalLength(planeSize/3, planeSize/3, planeSize/3);
    actorAxesOrigin->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(1);
    actorAxesOrigin->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(1);
    actorAxesOrigin->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(1);
    actorAxesOrigin->SetXAxisLabelText("");
    actorAxesOrigin->SetYAxisLabelText("");
    actorAxesOrigin->SetZAxisLabelText("");


//Generate Vector field
    vtkSmartPointer<vtkDoubleArray> vectorField = vtkSmartPointer<vtkDoubleArray>::New();
        vectorField->SetName("Vector_Field");
        vectorField->SetNumberOfComponents(3);
        vectorField->SetNumberOfTuples(polyDataPoints->GetNumberOfPoints());

//DYNAMIC ELEMENTS

//Points
    //Allocate memory for the  Mapper for the Points PolyData
    vtkSmartPointer<vtkPolyDataMapper> mapperPoints = vtkSmartPointer<vtkPolyDataMapper>::New();
    //Allocate memory for the Points actor
    vtkSmartPointer<vtkActor> actorPoints = vtkSmartPointer<vtkActor>::New();

/**Spline
    //Allocate memory for all the splines and related functions
    vtkSmartPointer<vtkKochanekSpline> xSpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkKochanekSpline> ySpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkKochanekSpline> zSpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    //Allocates memory for the Spline's actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> mapperSpline = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actorSpline = vtkSmartPointer<vtkActor>::New();*/

//Renderer and windows
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

//VISUALIZATION AND FIRST FRAME CREATION
    //Create Points PolyData Mapper
    mapperPoints->SetInputData(polyDataPoints);
    //Create Points PolyData actor
    actorPoints->SetMapper(mapperPoints);
    actorPoints->GetProperty()->SetPointSize(2);

/**Spline - In the first frame there is not passed path, only initial position
  spline->SetXSpline(xSpline);
  spline->SetYSpline(ySpline);
  spline->SetZSpline(zSpline);
  spline->SetPoints(points);

  functionSource->SetParametricFunction(spline);
  functionSource->SetUResolution(50 * points->GetNumberOfPoints());
  functionSource->SetVResolution(50 * points->GetNumberOfPoints());
  functionSource->SetWResolution(50 * points->GetNumberOfPoints());
  functionSource->Update();

  mapperSpline->SetInputConnection(functionSource->GetOutputPort());
  actorSpline->SetMapper(mapperSpline);
  actorSpline->GetProperty()->SetColor(colors->GetColor3d("DarkSlateGrey").GetData());
  actorSpline->GetProperty()->SetLineWidth(3.0);


  //Creates the initial structures that will be passed to the function called by the timer callback,
  //containing the objects to be updated
  UpdateParametersArguments updateParametersArguments(spline, functionSource, points, results, polyDataPoints,
                                                      verts, actorAxes, transform);
  void* arguments = static_cast<void*>(&updateParametersArguments);
*/


  //Add camera interactor
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  style->SetCurrentRenderer(renderer);


  vtkSmartPointer<vtkHedgeHog> hhog = vtkSmartPointer<vtkHedgeHog>::New();
  vtkSmartPointer<vtkPolyDataMapper> hhogMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();

  vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
  vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
  glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
  vtkSmartPointer<vtkPolyDataMapper> glyph3DMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
  glyph3DMapper->SetInputConnection(glyph3D->GetOutputPort());

  //glyph3D->SetVectorModeToUseNormal();
  vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
  generateVectorFieldAtT(polyDataPoints, vectorField, timeVector.at(0), colorLookupTable);

  polyDataPoints->GetPointData()->SetVectors(vectorField);
  polyDataPoints->GetPointData()->SetActiveScalars("Magnitudes");
  glyph3DMapper->SetLookupTable(colorLookupTable);
  glyph3DMapper->SetScalarRange(polyDataPoints->GetScalarRange());
  hhog->SetInputData(polyDataPoints);
  hhog->SetScaleFactor(0.001);
  hhog->Update();

  hhogMapper->SetInputConnection(hhog->GetOutputPort());

  vtkSmartPointer<vtkActor> hhogActor = vtkSmartPointer<vtkActor>::New();
  hhogActor->SetMapper(hhogMapper);

  glyph3D->SetInputData(polyDataPoints);
  glyph3D->SetScaleFactor(0.001);
  glyph3D->SetColorModeToColorByVector();
  glyph3D->SetScaleModeToScaleByVector();
  glyph3D->OrientOn();

  glyph3D->Update();

  vtkSmartPointer<vtkActor> glyph3DActor = vtkSmartPointer<vtkActor>::New();
  glyph3DActor->SetMapper(glyph3DMapper);

  glyph3DMapper->SetScalarModeToUsePointFieldData();
  glyph3DMapper->SetColorModeToMapScalars();
  glyph3DMapper->ScalarVisibilityOn();
  glyph3DMapper->SelectColorArray("Magnitudes");
  glyph3DMapper->Update();

//Render
    //Add actors
//    renderer->AddActor(actorPoints);
  //  renderer->AddActor(actorSpline);
 //   renderer->AddActor(actorPlaneXY);
 //   renderer->AddActor(actorPlaneXZ);
  //  renderer->AddActor(actorPlaneYZ);
    renderer->AddActor(actorAxesOrigin);
 //   renderer->AddActor(hhogActor);
    renderer->AddActor(glyph3DActor);
    //Set Background color
    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
    //Sets the window rendering process with the set renderer
    renderWindow->AddRenderer(renderer);
    //Sets a windown interactor with the rendering window
    renderWindowInteractor->SetRenderWindow(renderWindow);

    //Callback for the keypress event - Press enter and the animation starts
    vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    keypressCallback->SetCallback ( KeypressCallbackFunction );
    renderWindowInteractor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

    CallBackParameters callbackParams(polyDataPoints, vectorField, timeVector, glyph3DMapper);
    void* callBackArgums = static_cast<void*>(&callbackParams);

    //Initialize must be called prior to creating timer events.
    renderWindowInteractor->Initialize();
    vtkSmartPointer<vtkCallbackCommand> timerCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    timerCallback->SetCallback(TimerCallbackFunction);
    timerCallback->SetClientData(callBackArgums);
    renderWindowInteractor->AddObserver ( vtkCommand::TimerEvent, timerCallback );

  //Renders the window
  renderWindow->Render();
  //Starts the window interactor
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
