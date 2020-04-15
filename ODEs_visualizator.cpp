#include "resultsreader.h"
#include <algorithm>

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
    //    iren->CreateOneShotTimer(1);
        iren->CreateRepeatingTimer(20);
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

struct LocalCubeBounds{
    double xMin, xMax, yMin, yMax, zMin, zMax;
    LocalCubeBounds(const MeshParameters& meshParams, const double (&pos)[3]){
        int surroundingSize = 8;
        xMin = pos[0] - (1.01*(meshParams.nx/surroundingSize))*((meshParams.xMax-meshParams.xMin)/meshParams.nx);
        xMax = pos[0] + (1.01*(meshParams.nx/surroundingSize))*((meshParams.xMax-meshParams.xMin)/meshParams.nx);
        yMin = pos[1] - (1.01*(meshParams.ny/surroundingSize))*((meshParams.yMax-meshParams.yMin)/meshParams.ny);
        yMax = pos[1] + (1.01*(meshParams.ny/surroundingSize))*((meshParams.yMax-meshParams.yMin)/meshParams.ny);
        zMin = pos[2] - (1.01*(meshParams.nz/surroundingSize))*((meshParams.zMax-meshParams.zMin)/meshParams.nz);
        zMax = pos[2] + (1.01*(meshParams.nz/surroundingSize))*((meshParams.zMax-meshParams.zMin)/meshParams.nz);
    }
};

vtkNew<vtkPolyData> generateCube(LocalCubeBounds localBounds){
    std::array<std::array<double, 3>, 8> pts = {{{{localBounds.xMin, localBounds.yMin, localBounds.zMin}},
                                                 {{localBounds.xMin, localBounds.yMin, localBounds.zMax}},
                                                 {{localBounds.xMin, localBounds.yMax, localBounds.zMin}},
                                                 {{localBounds.xMin, localBounds.yMax, localBounds.zMax}},
                                                 {{localBounds.xMax, localBounds.yMin, localBounds.zMin}},
                                                 {{localBounds.xMax, localBounds.yMin, localBounds.zMax}},
                                                 {{localBounds.xMax, localBounds.yMax, localBounds.zMin}},
                                                 {{localBounds.xMax, localBounds.yMax, localBounds.zMax}}}};
    // The ordering of the corner points on each face.
    std::array<std::array<vtkIdType, 4>, 6> ordering = {{{{0, 1, 2, 3}},
                                                         {{4, 5, 6, 7}},
                                                         {{0, 1, 5, 4}},
                                                         {{1, 2, 6, 5}},
                                                         {{2, 3, 7, 6}},
                                                         {{3, 0, 4, 7}}}};
    vtkNew<vtkPolyData> cube;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;

    // Load the point, cell, and data attributes.
    for (auto i = 0ul; i < pts.size(); ++i)
    {
      points->InsertPoint(i, pts[i].data());
    }
    for (auto&& i : ordering)
    {
      polys->InsertNextCell(vtkIdType(i.size()), i.data());
    }

    // We now assign the pieces to the vtkPolyData.
    cube->SetPoints(points);
    cube->SetPolys(polys);

    return cube;
}

void evaluateFunction(const double t, const double (&position)[3], double (&result)[3]){
    result[0] = 10.0*(position[1] - position[0])*cos(t)*(t+1)/(t+1);
    result[1] = (28.0*position[0] - position[1]-position[0]*position[2])*sin(t)*(t+1)/(t+1);
    result[2] = (-(8.00/3.00)*position[2] + position[0]*position[1])*log(log(t+1)+1)*(t+1)/(t+1);
}


void getPositionsAroundSpline(const double (&posSpline)[3], const vtkSmartPointer<vtkPolyData>& polyDataPoints,
                            const vtkSmartPointer<vtkPolyData>& polyDataPointsLocal,
                            const MeshParameters& meshParams){
    vtkSmartPointer<vtkExtractEnclosedPoints> extractor = vtkSmartPointer<vtkExtractEnclosedPoints>::New();
    vtkNew<vtkCubeSource> cube;

    LocalCubeBounds localBounds(meshParams, posSpline);

    extractor->SetInputData(polyDataPoints);
    extractor->SetSurfaceData(generateCube(localBounds));
    extractor->SetOutput(polyDataPointsLocal);
    extractor->SetTolerance(.0001);
    extractor->CheckSurfaceOn();
    extractor->Update();
}


void generateLocalVectorFieldAtT(const vtkSmartPointer<vtkPolyData>& polyDataPoints,
                            const vtkSmartPointer<vtkDoubleArray>& vectorField, const double& t,
                            const vtkSmartPointer<vtkLookupTable>& colorLookupTable,
                            const vtkSmartPointer<vtkPolyData>& polyDataPointsLocal,
                            const double (&posSpline)[3], const MeshParameters& meshParams){

        vtkSmartPointer<vtkDoubleArray> magnitudes = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> magnitudesLocal = vtkSmartPointer<vtkDoubleArray>::New();

        polyDataPointsLocal->Reset();//Modificacion para Vector tangente solo
        polyDataPointsLocal->GetPoints()->InsertNextPoint(posSpline);//Modificacion para Vector tangente solo

        vtkSmartPointer<vtkDoubleArray> vectorFieldLocal = vtkSmartPointer<vtkDoubleArray>::New();
        vectorFieldLocal->SetName("Vector_Field_Local");
        vectorFieldLocal->SetNumberOfComponents(3);
        vectorFieldLocal->SetNumberOfTuples(polyDataPointsLocal->GetNumberOfPoints());

        vectorField->Reset();
        polyDataPoints->GetPointData()->Reset();
        polyDataPointsLocal->GetPointData()->Reset();

        for (int i=0; i<polyDataPoints->GetNumberOfPoints(); i++){
            double result[3];
            double position[3];

            polyDataPoints->GetPoint(i, position);

            evaluateFunction(t, position, result);
            vectorField->InsertTuple3(i, result[0], result[1], result[2]);
            double scalarResult = sqrt((pow(result[0],2)+pow(result[1],2)+pow(result[2],2)));
            magnitudes->InsertValue(i,scalarResult);
        }

        //Code for just drawing the tangent vector
        double result[3];
        evaluateFunction(t, posSpline, result);
        vectorFieldLocal->InsertTuple3(0, result[0], result[1], result[2]);
        double scalarResult = sqrt((pow(result[0],2)+pow(result[1],2)+pow(result[2],2)));
        magnitudesLocal->InsertValue(0,scalarResult);

        double scalarMagnitudes[2];
        polyDataPoints->GetPointData()->SetVectors(vectorField);
        polyDataPoints->GetPointData()->SetScalars(magnitudes);
        polyDataPoints->GetScalarRange(scalarMagnitudes);

        polyDataPointsLocal->GetPointData()->SetVectors(vectorFieldLocal);
        polyDataPointsLocal->GetPointData()->SetScalars(magnitudesLocal);

// Create the color map
        colorLookupTable->ResetAnnotations();
        colorLookupTable->SetTableRange(scalarMagnitudes[0], scalarMagnitudes[1]);
        colorLookupTable->Build();

        polyDataPoints->GetPointData()->GetScalars()->SetName("Magnitudes");
        polyDataPoints->GetPointData()->Update();

        polyDataPointsLocal->GetPointData()->GetScalars()->SetName("Magnitudes");
        polyDataPointsLocal->GetPointData()->Update();
}


//Dibuja todo el campo vectorial
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


void updateSpline(const double (&pos)[3], const vtkSmartPointer<vtkParametricSpline>& spline,
                  const vtkSmartPointer<vtkParametricFunctionSource>& functionSource,
                  const vtkSmartPointer<vtkPoints>& points)
{
    points->InsertNextPoint(pos);

    spline->SetNumberOfPoints(points->GetNumberOfPoints());
    spline->SetPoints(points);

    functionSource->SetParametricFunction(spline);
    functionSource->SetUResolution(50 * points->GetNumberOfPoints());
    functionSource->SetVResolution(50 * points->GetNumberOfPoints());
    functionSource->SetWResolution(50 * points->GetNumberOfPoints());
    functionSource->Update();
}


struct CallBackParameters{
    vtkSmartPointer<vtkPolyData> polyDataPoints;
    vtkSmartPointer<vtkDoubleArray> vectorField;
    std::vector<double> timeVector;

    vtkSmartPointer<vtkParametricSpline> spline;
    vtkSmartPointer<vtkParametricFunctionSource> functionSource;
    vtkSmartPointer<vtkPoints> pointsSpline;
    std::vector<ResultsLine> results;

    vtkSmartPointer<vtkPolyData> polyDataPointsLocal;
    MeshParameters meshParams;
    vtkSmartPointer<vtkCamera> camera;


    CallBackParameters(const vtkSmartPointer<vtkPolyData>& polyDataPoints_, const vtkSmartPointer<vtkDoubleArray>& vectorField_,
                       const std::vector<double>& timeVector_, const     vtkSmartPointer<vtkParametricSpline>& spline_,
                       const vtkSmartPointer<vtkParametricFunctionSource>& functionSource_,
                       const vtkSmartPointer<vtkPoints>& pointsSpline_, const std::vector<ResultsLine>& results_,
                       const vtkSmartPointer<vtkPolyData>& polyDataPointsLocal_, const MeshParameters& meshParams_,
                       const vtkSmartPointer<vtkCamera> camera_):
        polyDataPoints(polyDataPoints_),
        vectorField(vectorField_),
        timeVector(timeVector_),
        spline(spline_),
        functionSource(functionSource_),
        pointsSpline(pointsSpline_),
        results(results_),
        polyDataPointsLocal(polyDataPointsLocal_),
        meshParams(meshParams_),
        camera(camera_)
        {}
};

void TimerCallbackFunction ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData) )
{
    CallBackParameters* params = static_cast<CallBackParameters*>(clientData);
    vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

if(timePointsProcessedCounter<params->timeVector.size()){
//   draws the spline
   updateSpline(params->results[timePointsProcessedCounter].pos, params->spline, params->functionSource,
   params->pointsSpline);

    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
//This is the original vector field
//    generateVectorFieldAtT(params->polyDataPoints, params->vectorField,
//                           params->timeVector[timePointsProcessedCounter], colorLookupTable);

    //This part generates the "Tangent vector" highlighted effect
      generateLocalVectorFieldAtT(params->polyDataPoints, params->vectorField,
                             params->timeVector[timePointsProcessedCounter], colorLookupTable, params->polyDataPointsLocal,
                             params->results[timePointsProcessedCounter].pos, params->meshParams);

  //  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform> ::New();
 //   transform->RotateZ(0.5);
 //   params->camera->ApplyTransform(transform);

    iren->Render();
    timePointsProcessedCounter++;
}
}

int main(int argc, char* argv[])
{

    // Verify input arguments
    if ( argc != 2 )
    {
      std::cout << "Usage: " << argv[0]
                << " Filename(.xyz)" << std::endl;
      return EXIT_FAILURE;
    }

  // Read the file
    ResultsReader reader;
    std::vector<ResultsLine> results;
    //Reads the data, and exits in case there were any problem while reading points
    if (reader.read(std::string(argv[1]), results)==0){return 0;};


//Points Spline
        // Allocate objects to hold points and vertex cells.
        vtkSmartPointer<vtkPoints> pointsSpline = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> vertsSpline = vtkSmartPointer<vtkCellArray>::New();
        //Allocate memory for the  Mapper for the Spline
        vtkSmartPointer<vtkPolyDataMapper> mapperPointsSpline = vtkSmartPointer<vtkPolyDataMapper>::New();
        //Allocate memory for the Spline actor
        vtkSmartPointer<vtkActor> actorPointsSpline = vtkSmartPointer<vtkActor>::New();


    //VISUALIZATION AND FIRST FRAME CREATION
    //Points Spline
        vtkIdType id = pointsSpline->InsertNextPoint(results.begin()->pos);
        vertsSpline->InsertNextCell(1, &id);
        //Create Points PolyData actor
        actorPointsSpline->SetMapper(mapperPointsSpline);
        actorPointsSpline->GetProperty()->SetPointSize(2);


//Generation of the base mesh (positions where the vectors will be drawen
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkPolyData> polyDataPoints = vtkSmartPointer<vtkPolyData>::New();
        MeshParameters meshParams;
        std::vector<double> timeVector;

        double xMinMax[2];
        double yMinMax[2];
        double zMinMax[2];

        for (uint i=0; i<results.size(); i++){
            if(i == 0){
                xMinMax[0] = results[i].pos[0];
                xMinMax[1] = results[i].pos[0];
            }
            else if(results[i].pos[0] < xMinMax[0]){
                xMinMax[0] = results[i].pos[0];
            }
            else if(results[i].pos[0] > xMinMax[1]){
                xMinMax[1] = results[i].pos[0];
            }
        }

        for (uint i=0; i<results.size(); i++){
            if(i == 0){
                yMinMax[0] = results[i].pos[1];
                yMinMax[1] = results[i].pos[1];
            }
            else if(results[i].pos[1] < yMinMax[0]){
                yMinMax[0] = results[i].pos[1];
            }
            else if(results[i].pos[1] > yMinMax[1]){
                yMinMax[1] = results[i].pos[1];
            }
        }

        for (uint i=0; i<results.size(); i++){
            if(i == 0){
                zMinMax[0] = results[i].pos[2];
                zMinMax[1] = results[i].pos[2];
            }
            else if(results[i].pos[2] < zMinMax[0]){
                zMinMax[0] = results[i].pos[2];
            }
            else if(results[i].pos[2] > zMinMax[1]){
                zMinMax[1] = results[i].pos[2];
            }
        }

    meshParams.nx = 7;
    meshParams.ny = 7;
    meshParams.nz = 7;
    meshParams.nt = results.size();

    meshParams.xMin = 1.1*xMinMax[0];
    meshParams.yMin = 1.1*yMinMax[0];
    meshParams.zMin = 1.1*zMinMax[0];

    meshParams.xMax = 1.1*xMinMax[1];
    meshParams.yMax = 1.1*yMinMax[1];
    meshParams.zMax = 1.1*zMinMax[1];


    meshParams.tMin=0;
    meshParams.tMax=100;
    //meshParams.tMax=results.end()->t;

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
    actorAxesOrigin->SetXAxisLabelText("x");
    actorAxesOrigin->SetYAxisLabelText("y");
    actorAxesOrigin->SetZAxisLabelText("z");


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

//Spline setting
    //Allocate memory for all the splines and related functions
    vtkSmartPointer<vtkKochanekSpline> xSpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkKochanekSpline> ySpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkKochanekSpline> zSpline = vtkSmartPointer<vtkKochanekSpline>::New();
    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    //Allocates memory for the Spline's actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> mapperSpline = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actorSpline = vtkSmartPointer<vtkActor>::New();//*/

    //Spline - In the first frame there is not passed path, only initial position
      spline->SetXSpline(xSpline);
      spline->SetYSpline(ySpline);
      spline->SetZSpline(zSpline);
      spline->SetPoints(pointsSpline);

      functionSource->SetParametricFunction(spline);
      functionSource->SetUResolution(50 * points->GetNumberOfPoints());
      functionSource->SetVResolution(50 * points->GetNumberOfPoints());
      functionSource->SetWResolution(50 * points->GetNumberOfPoints());
      functionSource->Update();

      mapperSpline->SetInputConnection(functionSource->GetOutputPort());
      actorSpline->SetMapper(mapperSpline);
      actorSpline->GetProperty()->SetColor(colors->GetColor3d("Goldenrod").GetData());
      actorSpline->GetProperty()->SetLineWidth(3.0);

//LookUp Table (colors)
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();

//Generate and set vector field
    generateVectorFieldAtT(polyDataPoints, vectorField, timeVector.at(0), colorLookupTable);

    polyDataPoints->GetPointData()->SetVectors(vectorField);
    polyDataPoints->GetPointData()->SetActiveScalars("Magnitudes");
/**
//HedgeHog
    vtkSmartPointer<vtkHedgeHog> hhog = vtkSmartPointer<vtkHedgeHog>::New();
    vtkSmartPointer<vtkPolyDataMapper> hhogMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> hhogActor = vtkSmartPointer<vtkActor>::New();

    hhog->SetInputData(polyDataPoints);
    hhog->SetScaleFactor(0.01);
    hhog->SetVectorModeToUseVector();
    hhog->Update();

    hhogMapper->SetInputConnection(hhog->GetOutputPort());
    hhogMapper->SetLookupTable(colorLookupTable);
    hhogMapper->SetScalarRange(polyDataPoints->GetScalarRange());
    hhogMapper->SetScalarModeToUseFieldData();
    hhogMapper->SetColorModeToMapScalars();
    hhogMapper->ScalarVisibilityOn();
    hhogMapper->SelectColorArray("Magnitudes");
    hhogMapper->Update();

    hhogActor->SetMapper(hhogMapper);
*/
//Glyph3D
    vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkPolyDataMapper> glyph3DMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> glyph3DActor = vtkSmartPointer<vtkActor>::New();

    glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
    glyph3D->SetInputData(polyDataPoints);
    glyph3D->SetScaleFactor(0.025);
    glyph3D->SetColorModeToColorByVector();
    glyph3D->SetScaleModeToScaleByVector();
    glyph3D->OrientOn();
    glyph3D->Update();

    glyph3DMapper->SetInputConnection(glyph3D->GetOutputPort());
    glyph3DMapper->SetLookupTable(colorLookupTable);
    glyph3DMapper->SetScalarRange(polyDataPoints->GetScalarRange());
    glyph3DMapper->SetScalarModeToUsePointFieldData();
    glyph3DMapper->SetColorModeToMapScalars();
    glyph3DMapper->ScalarVisibilityOn();
    glyph3DMapper->SelectColorArray("Magnitudes");
    glyph3DMapper->Update();

    glyph3DActor->SetMapper(glyph3DMapper);
    glyph3DActor->GetProperty()->SetOpacity(0.025);


//Glyph3D Local points
     vtkSmartPointer<vtkPoints> pointsLocal = vtkSmartPointer<vtkPoints>::New();
     vtkSmartPointer<vtkCellArray> vertsLocal = vtkSmartPointer<vtkCellArray>::New();
     vtkSmartPointer<vtkPolyData> polyDataPointsLocal = vtkSmartPointer<vtkPolyData>::New();
     polyDataPointsLocal->SetPoints(pointsLocal);
     polyDataPointsLocal->SetVerts(vertsLocal);

//Glyph3D Local
    vtkSmartPointer<vtkGlyph3D> glyph3DLocal = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkPolyDataMapper> glyph3DLocalMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> glyph3DLocalActor = vtkSmartPointer<vtkActor>::New();

    glyph3DLocal->SetSourceConnection(arrowSource->GetOutputPort());
    glyph3DLocal->SetInputData(polyDataPointsLocal);
    glyph3DLocal->SetScaleFactor(0.05);
    glyph3DLocal->SetColorModeToColorByVector();
    glyph3DLocal->SetScaleModeToScaleByVector();
    glyph3DLocal->OrientOn();
    glyph3DLocal->Update();

    glyph3DLocalMapper->SetInputConnection(glyph3DLocal->GetOutputPort());
    glyph3DLocalMapper->SetLookupTable(colorLookupTable);
    glyph3DLocalMapper->SetScalarRange(polyDataPoints->GetScalarRange());
    glyph3DLocalMapper->SetScalarModeToUsePointFieldData();
    glyph3DLocalMapper->SetColorModeToMapScalars();
    glyph3DLocalMapper->ScalarVisibilityOn();
    glyph3DLocalMapper->SelectColorArray("Magnitudes");
    glyph3DLocalMapper->Update();

    glyph3DLocalActor->SetMapper(glyph3DLocalMapper);
    glyph3DLocalActor->GetProperty()->SetOpacity(1.0);


//Firs frame creation data creation
    mapperPoints->SetInputData(polyDataPoints);
    actorPoints->SetMapper(mapperPoints);
    actorPoints->GetProperty()->SetPointSize(2);

//Render
    //Renderer and windows
     vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
     vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
     vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();


//Camera
     vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
     camera->SetPosition(0, 120, 1);
     camera->SetFocalPoint(0, 0, 25);

     renderer->SetActiveCamera(camera);


    //Add camera interactor
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);
    style->SetCurrentRenderer(renderer);



    //Add actors
    //renderer->AddActor(actorPoints);
    //renderer->AddActor(actorPlaneXY);
    //renderer->AddActor(actorPlaneXZ);
    //renderer->AddActor(actorPlaneYZ);
    //renderer->AddActor(hhogActor);
    renderer->AddActor(glyph3DActor);
    renderer->AddActor(glyph3DLocalActor);
    //renderer->AddActor(actorAxesOrigin);
    renderer->AddActor(actorSpline);
    renderWindow->AddRenderer(renderer);
    renderWindowInteractor->SetRenderWindow(renderWindow);
    //Set Background color
    renderer->SetBackground(colors->GetColor3d("Black").GetData());
    //Sets the window rendering process with the set renderer
    //Sets a windown interactor with the rendering window

//Callback for the keypress event - Press enter and the animation starts
    vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    keypressCallback->SetCallback ( KeypressCallbackFunction );
    renderWindowInteractor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

//Callback for the timer event (refreshing the image with time)
//    CallBackParameters callbackParams(polyDataPoints, vectorField, timeVector, glyph3DMapper,
//                                      spline, functionSource, pointsSpline, results);

    CallBackParameters callbackParams(polyDataPoints, vectorField, timeVector, 
					spline, functionSource, pointsSpline, results, polyDataPointsLocal, meshParams, camera);
    void* callBackArgums = static_cast<void*>(&callbackParams);

    //Initialize must be called prior to creating timer events
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
