#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkColorTransferFunction.h>
#include <vtkDecimatePro.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <QFileInfo>

#include "cad_viewer.h"
#include <cstdio>
#include <QMouseEvent>
#include <gmsh.h>
#include <QWheelEvent>

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkActor2D.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include <vtkHardwareSelector.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkInformation.h>
#include <vtkAxesActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkLine.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <cstdio>
#include <vtkMapper.h>
#include <vtkFeatureEdges.h> 
#include <vtkCellPicker.h>
#include <vtkCoordinate.h>
#include <vtkMath.h>
#include <cmath>

#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>

// ====================================================================
// 【引擎 1】：视线深度投影锚定 
// ====================================================================
void AnchorPivot(vtkRenderer* renderer, int mouseX, int mouseY) {
    if (!renderer) return;
    vtkCamera* cam = renderer->GetActiveCamera();
    if (!cam) return;

    double pos[3], fp[3], dir[3];
    cam->GetPosition(pos);
    cam->GetFocalPoint(fp);

    dir[0] = fp[0] - pos[0]; dir[1] = fp[1] - pos[1]; dir[2] = fp[2] - pos[2];
    if(vtkMath::Norm(dir) < 1e-6) return;
    vtkMath::Normalize(dir);

    double bounds[6];
    renderer->ComputeVisiblePropBounds(bounds);
    double diag = 1.0;
    if (vtkMath::AreBoundsInitialized(bounds)) {
        diag = std::sqrt(std::pow(bounds[1]-bounds[0],2) + std::pow(bounds[3]-bounds[2],2) + std::pow(bounds[5]-bounds[4],2));
        if (diag < 1e-6) diag = 1.0;
    }

    double origin[3] = {0.0, 0.0, 0.0};
    double vecOrigin[3] = { origin[0]-pos[0], origin[1]-pos[1], origin[2]-pos[2] };
    double depthMacro = vtkMath::Dot(vecOrigin, dir); 
    if (depthMacro < diag * 0.01) depthMacro = diag * 0.5;
    double targetDepth = depthMacro;

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.001);
    if (picker->Pick(mouseX, mouseY, 0.0, renderer)) {
        double pickPos[3];
        picker->GetPickPosition(pickPos);
        double distToPick = std::sqrt(vtkMath::Distance2BetweenPoints(pos, pickPos));
        
        if (distToPick < diag * 0.5) {
            double vecPick[3] = { pickPos[0]-pos[0], pickPos[1]-pos[1], pickPos[2]-pos[2] };
            targetDepth = vtkMath::Dot(vecPick, dir); 
        }
    }

    double newFp[3] = { pos[0] + dir[0]*targetDepth, pos[1] + dir[1]*targetDepth, pos[2] + dir[2]*targetDepth };
    cam->SetFocalPoint(newFp);
    renderer->ResetCameraClippingRange();
}

// ====================================================================
// 【引擎 2】：指向性无级平滑缩放与空气墙防穿模
// ====================================================================
void performSmartZoom(vtkRenderer* renderer, vtkRenderWindowInteractor* interactor, int mouseX, int mouseY, bool zoomIn) {
    if (!renderer || !interactor) return;
    vtkCamera* cam = renderer->GetActiveCamera();
    if (!cam) return;

    double bounds[6];
    renderer->ComputeVisiblePropBounds(bounds);
    double diag = 1.0;
    if (vtkMath::AreBoundsInitialized(bounds)) {
        diag = std::sqrt(std::pow(bounds[1]-bounds[0],2) + std::pow(bounds[3]-bounds[2],2) + std::pow(bounds[5]-bounds[4],2));
        if (diag < 1e-6) diag = 1.0;
    }
    double minDist = diag * 0.0005; 

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.001);
    bool hit = picker->Pick(mouseX, mouseY, 0.0, renderer);
    
    double targetPos[3];
    if (hit) {
        picker->GetPickPosition(targetPos);
    } else {
        double fp[4]; cam->GetFocalPoint(fp); fp[3] = 1.0;
        renderer->SetWorldPoint(fp); renderer->WorldToDisplay();
        double* displayFp = renderer->GetDisplayPoint();
        renderer->SetDisplayPoint(mouseX, mouseY, displayFp[2]); renderer->DisplayToWorld();
        double* worldPos = renderer->GetWorldPoint();
        for(int i=0; i<3; ++i) targetPos[i] = worldPos[i]/worldPos[3];
    }

    double camPos[3], camFp[3];
    cam->GetPosition(camPos); cam->GetFocalPoint(camFp);
    double dir[3] = { targetPos[0] - camPos[0], targetPos[1] - camPos[1], targetPos[2] - camPos[2] };
    double distToTarget = vtkMath::Norm(dir);
    
    if (distToTarget < 1e-6) return;
    vtkMath::Normalize(dir);

    double step = 0.0;
    if (zoomIn) {
        step = distToTarget * 0.15; 
        if (hit && (distToTarget - step < minDist)) {
            step = distToTarget - minDist; 
        }
        if (step <= 1e-6) return; 
    } else {
        step = -distToTarget * 0.20; 
    }

    for (int i = 0; i < 3; ++i) {
        camPos[i] += dir[i] * step;
        camFp[i]  += dir[i] * step;
    }
    cam->SetPosition(camPos); cam->SetFocalPoint(camFp);

    renderer->ResetCameraClippingRange();
    double clip[2]; cam->GetClippingRange(clip);
    double newDist = distToTarget - step;
    double clipNear = newDist * 0.01; 
    if (clipNear < 1e-6) clipNear = 1e-6;
    cam->SetClippingRange(clipNear, clip[1]);
    
    interactor->GetRenderWindow()->Render();
}

// --------------------------------------------------------------------
// 边高亮辅助 (前向声明)
vtkSmartPointer<vtkActor> makeEdgeActor(const TopoDS_Edge& E, double r, double g, double b, double w);

// 商业级 Trackball 交互器
// --------------------------------------------------------------------
class CustomTrackballStyle : public vtkInteractorStyleTrackballCamera {
public:
    static CustomTrackballStyle* New();
    vtkTypeMacro(CustomTrackballStyle, vtkInteractorStyleTrackballCamera);
    CADViewer* Viewer = nullptr;
    int StartPos[2];

    void OnLeftButtonDown() override {
        int* clickPos = this->Interactor->GetEventPosition();
        this->FindPokedRenderer(clickPos[0], clickPos[1]);
        if (this->CurrentRenderer) AnchorPivot(this->CurrentRenderer, clickPos[0], clickPos[1]);
        this->Interactor->GetEventPosition(this->StartPos);
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown(); 
    }

    void OnMouseMove() override {
        vtkInteractorStyleTrackballCamera::OnMouseMove();
        if (Viewer && Viewer->m_edgeSelectMode) {
            int* pos = this->Interactor->GetEventPosition();
            this->FindPokedRenderer(pos[0], pos[1]);
            int oid, eidx; TopoDS_Edge E; gp_Pnt pt;
            Viewer->findClosestEdge(pos[0], pos[1], oid, eidx, E, pt);
            Viewer->clearEdgePreview();
            if (oid >= 0 && eidx >= 0) {
                CADViewer::EdgeKey key{oid, eidx};
                bool already = Viewer->m_selectedEdgeMap.count(key);
                // 已选中→绿色, 未选中→红色 (线宽5px覆盖下方蓝色4px)
                auto a = makeEdgeActor(E, already?0.2:1.0, already?0.9:0.3, already?0.2:0.2, 5.0);
                Viewer->m_renderer->AddActor(a);
                Viewer->m_previewEdge = a;
                Viewer->GetRenderWindow()->Render();
            }
        }
    }

    void OnLeftButtonUp() override {
        int endPos[2]; this->Interactor->GetEventPosition(endPos);
        int dx = endPos[0] - StartPos[0]; int dy = endPos[1] - StartPos[1];
        vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
        if (std::abs(dx) < 5 && std::abs(dy) < 5 && Viewer) {
            Viewer->processSingleClickSelection(endPos[0], endPos[1]);
        }
    }

    void OnMouseWheelForward() override {
        int* pos = this->Interactor->GetEventPosition();
        this->FindPokedRenderer(pos[0], pos[1]); 
        performSmartZoom(this->CurrentRenderer, this->Interactor, pos[0], pos[1], true);
    }
    void OnMouseWheelBackward() override {
        int* pos = this->Interactor->GetEventPosition();
        this->FindPokedRenderer(pos[0], pos[1]); 
        performSmartZoom(this->CurrentRenderer, this->Interactor, pos[0], pos[1], false);
    }
};
vtkStandardNewMacro(CustomTrackballStyle);

// --------------------------------------------------------------------
// 自定义 RubberBand (框选) 交互器
// --------------------------------------------------------------------
class CustomRubberBandStyle : public vtkInteractorStyleTrackballCamera {
public:
    static CustomRubberBandStyle* New();
    vtkTypeMacro(CustomRubberBandStyle, vtkInteractorStyleTrackballCamera);
    CADViewer* Viewer = nullptr;
    bool IsSelecting = false; int StartPos[2];
    vtkSmartPointer<vtkPolyData> BoxPolyData; vtkSmartPointer<vtkActor2D> BoxActor;

    CustomRubberBandStyle() {
        BoxPolyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New(); points->SetNumberOfPoints(4);
        BoxPolyData->SetPoints(points);
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New(); lines->InsertNextCell(5);
        lines->InsertCellPoint(0); lines->InsertCellPoint(1); lines->InsertCellPoint(2); lines->InsertCellPoint(3); lines->InsertCellPoint(0); 
        BoxPolyData->SetLines(lines);
        vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New(); mapper->SetInputData(BoxPolyData);
        BoxActor = vtkSmartPointer<vtkActor2D>::New(); BoxActor->SetMapper(mapper);
        BoxActor->GetProperty()->SetColor(0.2, 0.6, 1.0); BoxActor->GetProperty()->SetLineWidth(1.5);
    }
    void OnLeftButtonDown() override {
        this->IsSelecting = true; this->Interactor->GetEventPosition(this->StartPos);
        this->FindPokedRenderer(this->StartPos[0], this->StartPos[1]); 
        if (this->CurrentRenderer) { this->CurrentRenderer->AddActor2D(BoxActor); UpdateBox(this->StartPos); }
    }
    void OnMouseMove() override {
        if (this->IsSelecting && this->CurrentRenderer) {
            int endPos[2]; this->Interactor->GetEventPosition(endPos); UpdateBox(endPos); 
        } else { vtkInteractorStyleTrackballCamera::OnMouseMove(); }
    }
    void OnLeftButtonUp() override {
        if (!this->IsSelecting) { vtkInteractorStyleTrackballCamera::OnLeftButtonUp(); return; }
        this->IsSelecting = false; int endPos[2]; this->Interactor->GetEventPosition(endPos);
        if (this->CurrentRenderer) { this->CurrentRenderer->RemoveActor2D(BoxActor); this->Interactor->GetRenderWindow()->Render(); }
        if (Viewer) {
            int x0 = StartPos[0]; int y0 = StartPos[1]; int x1 = endPos[0]; int y1 = endPos[1];
            if (std::abs(x1 - x0) < 5 && std::abs(y1 - y0) < 5) Viewer->processSingleClickSelection(x1, y1);
            else Viewer->processBoxSelection(std::min(x0, x1), std::min(y0, y1), std::max(x0, x1), std::max(y0, y1));
            emit Viewer->boxSelectionCompleted(); 
        }
    }
    void UpdateBox(int* endPos) {
        vtkPoints* pts = BoxPolyData->GetPoints();
        pts->SetPoint(0, StartPos[0], StartPos[1], 0); pts->SetPoint(1, endPos[0], StartPos[1], 0);
        pts->SetPoint(2, endPos[0], endPos[1], 0); pts->SetPoint(3, StartPos[0], endPos[1], 0);
        pts->Modified(); this->Interactor->GetRenderWindow()->Render();
    }
    
    void OnMouseWheelForward() override { 
        int* pos = this->Interactor->GetEventPosition();
        this->FindPokedRenderer(pos[0], pos[1]);
        performSmartZoom(this->CurrentRenderer, this->Interactor, pos[0], pos[1], true);
    }
    void OnMouseWheelBackward() override { 
        int* pos = this->Interactor->GetEventPosition();
        this->FindPokedRenderer(pos[0], pos[1]);
        performSmartZoom(this->CurrentRenderer, this->Interactor, pos[0], pos[1], false);
    }
};
vtkStandardNewMacro(CustomRubberBandStyle);

// --------------------------------------------------------------------
// CADViewer 主体
// --------------------------------------------------------------------
CADViewer::CADViewer(QWidget *parent) : QVTKWidget(parent) {
    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    setTheme(true);
    this->GetRenderWindow()->AddRenderer(m_renderer);

    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
    vtkMapper::SetResolveCoincidentTopologyPolygonOffsetParameters(1.0, 1.0);
    vtkMapper::SetResolveCoincidentTopologyLineOffsetParameters(-2.0, -2.0);  

    m_picker = vtkSmartPointer<vtkPropPicker>::New();

    vtkCamera* camera = m_renderer->GetActiveCamera();
    camera->SetPosition(0, 0, 1); camera->SetFocalPoint(0, 0, 0); camera->SetViewUp(0, 1, 0);
    m_renderer->ResetCamera(-1, 1, -1, 1, -1, 1); 

    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
    axes->SetXAxisLabelText("x"); axes->SetYAxisLabelText("y"); axes->SetZAxisLabelText("z");
    axes->SetTotalLength(1.0, 1.0, 1.0);
    m_axesWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    m_axesWidget->SetOrientationMarker(axes);
    m_axesWidget->SetInteractor(this->GetRenderWindow()->GetInteractor());
    m_axesWidget->SetViewport(0.0, 0.0, 0.2, 0.2); 
    this->GetRenderWindow()->GetInteractor()->Initialize();
    m_axesWidget->SetEnabled(1); m_axesWidget->InteractiveOff();

    // ---- COMSOL 风格栅格 (三个面贴着几何包围盒) ----
    m_gridPlanes.resize(3);
    for (auto& gp : m_gridPlanes) {
        gp.majorActor = vtkSmartPointer<vtkActor>::New();
        gp.majorActor->GetProperty()->SetColor(0.55, 0.55, 0.55);
        gp.majorActor->GetProperty()->SetLineWidth(1.2);
        gp.majorActor->SetPickable(false);
        gp.majorActor->SetVisibility(false);
        gp.minorActor = vtkSmartPointer<vtkActor>::New();
        gp.minorActor->GetProperty()->SetColor(0.75, 0.75, 0.75);
        gp.minorActor->GetProperty()->SetLineWidth(1.0);
        gp.minorActor->GetProperty()->SetOpacity(0.8);
        gp.minorActor->SetPickable(false);
        gp.minorActor->SetVisibility(false);
        m_renderer->AddActor(gp.majorActor);
        m_renderer->AddActor(gp.minorActor);
    }

    setInteractionModes(false, false, false);
}

CADViewer::~CADViewer() {}

void CADViewer::setGeometryTransparent(bool isTransparent) {
    if (m_isTransparent != isTransparent) {
        m_isTransparent = isTransparent;
        updateFaceColors();
    }
}

void CADViewer::setTheme(bool isLight) {
    if (isLight) {
        m_renderer->SetBackground(1.0, 1.0, 1.0); m_renderer->SetBackground2(0.8, 0.85, 0.95);
        m_renderer->SetGradientBackground(true);
    } else {
        m_renderer->SetGradientBackground(false); m_renderer->SetBackground(0.15, 0.15, 0.2);
    }
    if(this->GetRenderWindow()) this->GetRenderWindow()->Render();
}

void CADViewer::clearScene() {
    for (auto const& [actor, id] : m_actorToFaceId) m_renderer->RemoveActor(actor);
    for (auto const& [id, actor] : m_idToEdgeActor) m_renderer->RemoveActor(actor);

    m_actorToFaceId.clear(); m_idToActor.clear(); m_idToEdgeActor.clear();
    m_selectedFaceIds.clear(); m_hiddenFaceIds.clear();
    m_faceToObject.clear();
    m_occFaceIdCounter = 100000;
    m_highlightedEdge = nullptr;
    m_viewMode = VIEW_NORMAL;

    clearAllEdgeHighlights();
    clearMesh();
}

void CADViewer::loadShape(const TopoDS_Shape& shape) {
    BRepMesh_IncrementalMesh(const_cast<TopoDS_Shape&>(shape), 0.1);
    int faceIdCounter = 1;
    for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        TopLoc_Location location;
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        if (triangulation.IsNull()) { faceIdCounter++; continue; }

        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        auto points = vtkSmartPointer<vtkPoints>::New();
        auto triangles = vtkSmartPointer<vtkCellArray>::New();
        for (int i = 1; i <= triangulation->NbNodes(); ++i) {
            gp_Pnt pt = triangulation->Node(i); pt.Transform(location.Transformation());
            points->InsertNextPoint(pt.X(), pt.Y(), pt.Z());
        }
        for (int i = 1; i <= triangulation->NbTriangles(); ++i) {
            int n1, n2, n3; triangulation->Triangle(i).Get(n1, n2, n3);
            vtkIdType ids[3] = {n1-1, n2-1, n3-1};
            triangles->InsertNextCell(3, ids);
        }
        polyData->SetPoints(points); polyData->SetPolys(triangles);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.85, 0.88, 0.9);
        m_renderer->AddActor(actor);

        auto featEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featEdges->SetInputData(polyData);
        featEdges->BoundaryEdgesOn(); featEdges->FeatureEdgesOn();
        featEdges->SetFeatureAngle(30); featEdges->ManifoldEdgesOff();
        featEdges->NonManifoldEdgesOff(); featEdges->Update();
        auto edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        edgeMapper->SetInputData(featEdges->GetOutput());
        edgeMapper->ScalarVisibilityOff();
        auto edgeActor = vtkSmartPointer<vtkActor>::New();
        edgeActor->SetMapper(edgeMapper);
        edgeActor->GetProperty()->SetColor(0, 0, 0);
        edgeActor->GetProperty()->SetLineWidth(2.0);
        edgeActor->GetProperty()->SetLighting(false);
        m_renderer->AddActor(edgeActor);

        m_actorToFaceId[actor.Get()] = faceIdCounter;
        m_idToActor[faceIdCounter] = actor.Get();
        m_idToEdgeActor[faceIdCounter] = edgeActor;
        faceIdCounter++;
    }
    updateGridExtent();
    m_renderer->ResetCamera(); this->GetRenderWindow()->Render();
}

bool CADViewer::loadStepFile(const QString& filePath) {
    clearScene();
    STEPControl_Reader reader;
    if (reader.ReadFile(filePath.toUtf8().constData()) != IFSelect_RetDone) return false;
    reader.TransferRoots();
    loadShape(reader.OneShape());
    return true;
}

bool CADViewer::loadIgesFile(const QString& filePath) {
    clearScene();
    IGESControl_Reader reader;
    if (reader.ReadFile(filePath.toUtf8().constData()) != IFSelect_RetDone) return false;
    reader.TransferRoots();
    loadShape(reader.OneShape());
    return true;
}

bool CADViewer::loadBrepFile(const QString& filePath) {
    clearScene();
    TopoDS_Shape shape;
    BRep_Builder B;
    if (!BRepTools::Read(shape, filePath.toUtf8().constData(), B)) return false;
    loadShape(shape);
    return true;
}

// ====================================================================
// 将 OCC Shape 注册为可拾取的面 (与 STEP 导入面一致)
// ====================================================================
int CADViewer::addOCCShape(const TopoDS_Shape& shape, double r, double g, double b) {
    fprintf(stderr, "[DEBUG addOCCShape] shape.IsNull=%d, faceCount via explorer...\n", shape.IsNull());
    if (shape.IsNull()) return -1;
    BRepMesh_IncrementalMesh(const_cast<TopoDS_Shape&>(shape), 0.5);
    int firstId = m_occFaceIdCounter;
    int nFaces = 0;

    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        TopLoc_Location loc;
        Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
        if (tri.IsNull()) continue;

        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        auto pts = vtkSmartPointer<vtkPoints>::New();
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for (int i = 1; i <= tri->NbNodes(); i++) {
            gp_Pnt p = tri->Node(i).Transformed(loc.Transformation());
            pts->InsertNextPoint(p.X(), p.Y(), p.Z());
        }
        for (int i = 1; i <= tri->NbTriangles(); i++) {
            int n1,n2,n3; tri->Triangle(i).Get(n1,n2,n3);
            vtkIdType ids[3]={n1-1,n2-1,n3-1};
            cells->InsertNextCell(3,ids);
        }
        polyData->SetPoints(pts); polyData->SetPolys(cells);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(r, g, b);
        m_renderer->AddActor(actor);

        auto featEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featEdges->SetInputData(polyData);
        featEdges->BoundaryEdgesOn(); featEdges->FeatureEdgesOn();
        featEdges->SetFeatureAngle(30); featEdges->ManifoldEdgesOff();
        featEdges->NonManifoldEdgesOff(); featEdges->Update();
        auto edgeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        edgeMapper->SetInputData(featEdges->GetOutput());
        edgeMapper->ScalarVisibilityOff();
        auto edgeActor = vtkSmartPointer<vtkActor>::New();
        edgeActor->SetMapper(edgeMapper);
        edgeActor->GetProperty()->SetColor(0,0,0);
        edgeActor->GetProperty()->SetLineWidth(2.0);
        edgeActor->GetProperty()->SetLighting(false);
        m_renderer->AddActor(edgeActor);

        m_actorToFaceId[actor.Get()] = m_occFaceIdCounter;
        m_idToActor[m_occFaceIdCounter] = actor.Get();
        m_idToEdgeActor[m_occFaceIdCounter] = edgeActor;
        m_faceToObject[m_occFaceIdCounter] = firstId;
        m_faceIdToShape[m_occFaceIdCounter] = shape;
        m_occFaceIdCounter++;
        nFaces++;
    }
    fprintf(stderr, "[DEBUG addOCCShape] firstId=%d, nFaces=%d\n", firstId, nFaces);
    updateGridExtent();   // 内部创建的体素也更新栅格范围
    m_renderer->ResetCamera();
    GetRenderWindow()->Render();
    return firstId;
}

// ====================================================================
// 【核心管线控制】：统一管理透明、网格遮挡与选取逻辑
// ====================================================================
void CADViewer::updateFaceColors() {
    for (auto const& [id, actor] : m_idToActor) {
        if (m_geometryHidden) {
            actor->SetVisibility(false);
            if (m_idToEdgeActor.count(id)) m_idToEdgeActor[id]->SetVisibility(false);
            if (m_idToMeshActor.count(id)) m_idToMeshActor[id]->SetVisibility(false);
            continue;
        }
        bool shouldShow = true;
        double r, g, b;

        // 1. 判断颜色态
        if (m_viewMode == VIEW_HIDDEN_ONLY) {
            shouldShow = m_hiddenFaceIds.count(id);
            r = 0.9; g = 0.3; b = 0.3; 
        } else if (m_viewMode == VIEW_ALL) {
            shouldShow = true;
            if (m_selectedFaceIds.count(id)) { r = 0.0; g = 0.8; b = 1.0; }
            else { r = 0.85; g = 0.88; b = 0.9; }
        } else {
            shouldShow = (m_hiddenFaceIds.count(id) == 0); 
            if (m_selectedFaceIds.count(id)) { r = 0.0; g = 0.8; b = 1.0; }
            else { r = 0.85; g = 0.88; b = 0.9; }
        }

        actor->GetProperty()->SetColor(r, g, b);

        // 2. 网格审查模式
        if (m_isTransparent) {
            // 【终极方案】：如果该面有网格，完全隐藏底层 CAD 实体，用网格实体面代替，杜绝Z-Fighting
            if (m_idToMeshActor.count(id)) {
                actor->SetVisibility(false);
                
                m_idToMeshActor[id]->SetVisibility(shouldShow);
                m_idToMeshActor[id]->GetProperty()->SetColor(r, g, b); // 颜色与CAD保持一致
                m_idToMeshActor[id]->GetProperty()->SetOpacity(1.0);   // 网格面完全不透明
            } else {
                // 没有网格的面，作为背景半透明显示
                actor->SetVisibility(shouldShow);
                actor->GetProperty()->SetOpacity(0.20); 
            }
            
            // 始终保留 CAD 黑色特征边界轮廓
            if (m_idToEdgeActor.count(id)) {
                m_idToEdgeActor[id]->SetVisibility(shouldShow); 
                m_idToEdgeActor[id]->GetProperty()->SetOpacity(1.0);
            }
        } else {
            // 普通建模模式
            actor->SetVisibility(shouldShow);
            actor->GetProperty()->SetOpacity(1.0); 
            
            if (m_idToEdgeActor.count(id)) {
                m_idToEdgeActor[id]->SetVisibility(shouldShow); 
                m_idToEdgeActor[id]->GetProperty()->SetOpacity(1.0);
            }
            if (m_idToMeshActor.count(id)) m_idToMeshActor[id]->SetVisibility(false);
        }
    }

    // 处理无对应CAD面的网格面(内部几何: CAD ID 100000+ vs Gmsh tag 1,2...)
    for (auto& [faceId, meshActor] : m_idToMeshActor) {
        if (m_idToActor.count(faceId)) continue; // 已在上面循环处理
        if (m_geometryHidden) { meshActor->SetVisibility(m_hiddenFaceIds.count(faceId) == 0); continue; }
        bool vis = true;
        if (m_viewMode == VIEW_HIDDEN_ONLY) vis = m_hiddenFaceIds.count(faceId);
        else if (m_viewMode == VIEW_ALL) vis = true;
        else vis = (m_hiddenFaceIds.count(faceId) == 0);
        meshActor->SetVisibility(vis);
        if (vis) { meshActor->GetProperty()->SetColor(0.85, 0.88, 0.9); meshActor->GetProperty()->SetOpacity(1.0); }
    }

    this->GetRenderWindow()->Render();
}

void CADViewer::setSelectedFaces(const std::set<int>& faceIds) { m_selectedFaceIds = faceIds; updateFaceColors(); }
void CADViewer::toggleFaceSelection(int faceId) { if (m_selectedFaceIds.count(faceId)) m_selectedFaceIds.erase(faceId); else m_selectedFaceIds.insert(faceId); updateFaceColors(); emit faceSelectionChanged(m_selectedFaceIds); }
void CADViewer::clearSelection() { m_selectedFaceIds.clear(); clearAllEdgeHighlights(); updateFaceColors(); emit faceSelectionChanged(m_selectedFaceIds); }

void CADViewer::setInteractionModes(bool boxSelect, bool encloseMode, bool hideMode) {
    m_boxSelectMode = boxSelect; m_encloseMode = encloseMode; m_hideMode = hideMode;
    if (m_boxSelectMode) {
        vtkSmartPointer<CustomRubberBandStyle> style = vtkSmartPointer<CustomRubberBandStyle>::New();
        style->Viewer = this; this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
    } else {
        vtkSmartPointer<CustomTrackballStyle> style = vtkSmartPointer<CustomTrackballStyle>::New();
        style->Viewer = this; this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
    }
}

void CADViewer::processSingleClickSelection(int x,int y) {
    vtkSmartPointer<vtkCellPicker> cellPicker = vtkSmartPointer<vtkCellPicker>::New();
    cellPicker->SetTolerance(0.005);
    bool hit = cellPicker->Pick(x, y, 0, m_renderer);
    vtkActor* pickedActor = hit ? vtkActor::SafeDownCast(cellPicker->GetViewProp()) : nullptr;
    if (!pickedActor) { m_picker->Pick(x, y, 0, m_renderer); pickedActor = vtkActor::SafeDownCast(m_picker->GetViewProp()); }
    if (pickedActor && m_actorToFaceId.count(pickedActor)) {
        int faceId = m_actorToFaceId[pickedActor];
        if (m_edgeSelectMode) {
            int oid, eidx; TopoDS_Edge E; gp_Pnt pt;
            findClosestEdge(x, y, oid, eidx, E, pt);
            clearEdgePreview();
            if (oid >= 0 && eidx >= 0) {
                // 统一使用IsSame几何匹配: 遍历所有条目移除重复项
                bool removed = false;
                for (auto it = m_selectedEdgeMap.begin(); it != m_selectedEdgeMap.end(); ) {
                    if (E.IsSame(it->second.second)) {
                        m_renderer->RemoveActor(it->second.first);
                        it = m_selectedEdgeMap.erase(it);
                        removed = true;
                    } else {
                        ++it;
                    }
                }
                if (!removed) {
                    auto a = makeEdgeActor(E, 0.2, 0.8, 1.0, 4.0);
                    m_renderer->AddActor(a);
                    m_selectedEdgeMap[{oid, eidx}] = {a, E};
                }
                m_lastEdgeIdx = eidx; m_lastObjId = oid;
                GetRenderWindow()->Render();
                emit faceSelectionChanged({faceId});
            }
            return;
        }
        if (m_objectSelectMode && m_faceToObject.count(faceId)) {
                // 对象选择模式: 选中整个对象的所有面
                int objId = m_faceToObject[faceId];
                if (m_hideMode) {
                    for (auto& [fid, oid] : m_faceToObject)
                        if (oid == objId) m_hiddenFaceIds.insert(fid);
                    emit faceSelectionChanged(m_selectedFaceIds);
                } else {
                    std::set<int> objFaces;
                    for (auto& [fid, oid] : m_faceToObject)
                        if (oid == objId) objFaces.insert(fid);
                    // toggle: 若全部已选则取消, 否则选中
                    bool allSel = true;
                    for (int fid : objFaces) if (!m_selectedFaceIds.count(fid)) { allSel = false; break; }
                    if (allSel) for (int fid : objFaces) m_selectedFaceIds.erase(fid);
                    else for (int fid : objFaces) m_selectedFaceIds.insert(fid);
                    emit faceSelectionChanged(m_selectedFaceIds);
                }
            } else {
                // 边界选择模式
                if (m_hideMode) { m_hiddenFaceIds.insert(faceId); }
                else { toggleFaceSelection(faceId); }
            }
            updateFaceColors();
        }
}

void CADViewer::processBoxSelection(int xmin, int ymin, int xmax, int ymax) {
    // 边选择模式: 直接遍历所有边投影到屏幕检查, 不依赖面拾取
    if (m_edgeSelectMode) {
        vtkSmartPointer<vtkCoordinate> coord = vtkSmartPointer<vtkCoordinate>::New();
        coord->SetViewport(m_renderer);
        std::set<EdgeKey> newEdges;
        std::set<int> processedObjIds;
        for (auto& [fid, oid] : m_faceToObject) {
            if (processedObjIds.count(oid)) continue;
            processedObjIds.insert(oid);
            if (!m_faceIdToShape.count(fid)) continue;
            const TopoDS_Shape& sh = m_faceIdToShape[fid];
            int ei = 0;
            for (TopExp_Explorer ex(sh, TopAbs_EDGE); ex.More(); ex.Next(), ei++) {
                TopoDS_Edge E = TopoDS::Edge(ex.Current());
                try {
                    BRepAdaptor_Curve adapt(E);
                    double u1 = adapt.FirstParameter(), u2 = adapt.LastParameter();
                    int insideCount = 0, nSamples = 9;
                    for (int s = 0; s <= nSamples; s++) {
                        double u = u1 + (u2 - u1) * s / (double)nSamples;
                        gp_Pnt pt; adapt.D0(u, pt);
                        coord->SetValue(pt.X(), pt.Y(), pt.Z());
                        int* dp = coord->GetComputedDisplayValue(m_renderer);
                        if (dp && dp[0] >= xmin && dp[0] <= xmax && dp[1] >= ymin && dp[1] <= ymax)
                            insideCount++;
                    }
                    if ((!m_encloseMode && insideCount > 0) || (m_encloseMode && insideCount > nSamples)) {
                        EdgeKey key{oid, ei};
                        if (!m_selectedEdgeMap.count(key)) {
                            auto a = makeEdgeActor(E, 0.2, 0.8, 1.0, 4.0);
                            m_renderer->AddActor(a);
                            m_selectedEdgeMap[key] = {a, E};
                        }
                        newEdges.insert(key);
                    }
                } catch (...) {}
            }
        }
        for (auto it = m_selectedEdgeMap.begin(); it != m_selectedEdgeMap.end(); ) {
            if (!newEdges.count(it->first)) {
                m_renderer->RemoveActor(it->second.first);
                it = m_selectedEdgeMap.erase(it);
            } else ++it;
        }
        GetRenderWindow()->Render();
        m_lastObjId = -1; m_lastEdgeIdx = -1;
        if (!m_selectedEdgeMap.empty())
            emit faceSelectionChanged({m_selectedEdgeMap.begin()->first.objId});
        return;
    }

    // 面/对象模式: 先拾取面
    std::set<int> pickedFaceIds;
    if (m_encloseMode) {
        for (auto const& [actor, id] : m_actorToFaceId) {
            if (actor->GetVisibility() == false) continue;
            vtkPolyData* pd = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInput());
            if (!pd) continue;
            bool isFullyEnclosed = true;
            for (vtkIdType i = 0; i < pd->GetNumberOfPoints(); ++i) {
                double p[3]; pd->GetPoint(i, p);
                m_renderer->SetWorldPoint(p[0], p[1], p[2], 1.0); m_renderer->WorldToDisplay();
                double* dpt = m_renderer->GetDisplayPoint();
                if (dpt[0] < xmin || dpt[0] > xmax || dpt[1] < ymin || dpt[1] > ymax) { isFullyEnclosed = false; break; }
            }
            if (isFullyEnclosed) pickedFaceIds.insert(id);
        }
    } else {
        vtkSmartPointer<vtkHardwareSelector> selector = vtkSmartPointer<vtkHardwareSelector>::New();
        selector->SetRenderer(m_renderer); selector->SetArea(xmin, ymin, xmax, ymax);
        selector->SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
        vtkSmartPointer<vtkSelection> selection; selection.TakeReference(selector->Select());
        if (selection) {
            for (unsigned int i = 0; i < selection->GetNumberOfNodes(); i++) {
                vtkSelectionNode* node = selection->GetNode(i);
                vtkActor* actor = vtkActor::SafeDownCast(node->GetProperties()->Get(vtkSelectionNode::PROP()));
                if (actor && m_actorToFaceId.count(actor)) pickedFaceIds.insert(m_actorToFaceId[actor]);
            }
        }
    }
    if (pickedFaceIds.empty()) return;

    if (m_objectSelectMode) {
        // 对象选择模式: 将面按对象分组
        std::set<int> objIds;
        for (int fid : pickedFaceIds) {
            auto it = m_faceToObject.find(fid);
            if (it != m_faceToObject.end()) objIds.insert(it->second);
        }
        std::set<int> allObjFaces;
        for (auto& [fid, oid] : m_faceToObject)
            if (objIds.count(oid)) allObjFaces.insert(fid);
        if (m_hideMode) m_hiddenFaceIds.insert(allObjFaces.begin(), allObjFaces.end());
        else m_selectedFaceIds.insert(allObjFaces.begin(), allObjFaces.end());
        updateFaceColors();
        emit faceSelectionChanged(m_selectedFaceIds);
        return;
    }

    // 默认面选择模式
    if (m_hideMode) m_hiddenFaceIds.insert(pickedFaceIds.begin(), pickedFaceIds.end());
    else m_selectedFaceIds.insert(pickedFaceIds.begin(), pickedFaceIds.end());
    updateFaceColors();
    emit faceSelectionChanged(m_selectedFaceIds);
}

void CADViewer::mousePressEvent(QMouseEvent *event) { QVTKWidget::mousePressEvent(event); }
void CADViewer::wheelEvent(QWheelEvent *event) { QVTKWidget::wheelEvent(event); }

void CADViewer::setViewXY() { m_renderer->GetActiveCamera()->SetFocalPoint(0,0,0); m_renderer->GetActiveCamera()->SetPosition(0,0,1); m_renderer->GetActiveCamera()->SetViewUp(0,1,0); m_renderer->ResetCamera(); m_renderer->ResetCameraClippingRange(); this->GetRenderWindow()->Render(); }
void CADViewer::setViewYZ() { m_renderer->GetActiveCamera()->SetFocalPoint(0,0,0); m_renderer->GetActiveCamera()->SetPosition(1,0,0); m_renderer->GetActiveCamera()->SetViewUp(0,0,1); m_renderer->ResetCamera(); m_renderer->ResetCameraClippingRange(); this->GetRenderWindow()->Render(); }
void CADViewer::setViewZX() { m_renderer->GetActiveCamera()->SetFocalPoint(0,0,0); m_renderer->GetActiveCamera()->SetPosition(0,-1,0); m_renderer->GetActiveCamera()->SetViewUp(0,0,1); m_renderer->ResetCamera(); m_renderer->ResetCameraClippingRange(); this->GetRenderWindow()->Render(); }
void CADViewer::resetCameraView() {
    m_renderer->ResetCamera(); m_renderer->ResetCameraClippingRange();
    updateGridExtent();
    this->GetRenderWindow()->Render();
}
void CADViewer::hideSelected() { m_hiddenFaceIds.insert(m_selectedFaceIds.begin(), m_selectedFaceIds.end()); updateFaceColors(); }
void CADViewer::showUnhidden()   { m_viewMode = VIEW_NORMAL;      updateFaceColors(); }
void CADViewer::showOnlyHidden() { m_viewMode = VIEW_HIDDEN_ONLY; updateFaceColors(); }
void CADViewer::showAll()        { m_viewMode = VIEW_ALL;         updateFaceColors(); }
void CADViewer::resetHidden() { m_hiddenFaceIds.clear(); updateFaceColors(); }

void CADViewer::clearMesh() {
    for (auto const& [id, actor] : m_idToMeshActor) {
        m_renderer->RemoveActor(actor);
        // 【关键】：从 Picker 映射中注销
        m_actorToFaceId.erase(actor);
    }
    m_idToMeshActor.clear();
    this->GetRenderWindow()->Render();
}

// ====================================================================
// 【核心修改】：二维网格以“Surface + 黑色黑框”形式渲染
// ====================================================================
void CADViewer::loadMesh(const FaceMeshMap& meshData) {
    clearMesh(); 

    for (const auto& pair : meshData) {
        int faceId = pair.first;
        const RenderableMesh& rMesh = pair.second;

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (size_t i = 0; i < rMesh.vertices.size(); i += 3) {
            points->InsertNextPoint(rMesh.vertices[i], rMesh.vertices[i + 1], rMesh.vertices[i + 2]);
        }

        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
        for (size_t i = 0; i < rMesh.indices.size(); i += 3) {
            vtkIdType pts[3] = { static_cast<vtkIdType>(rMesh.indices[i]), static_cast<vtkIdType>(rMesh.indices[i + 1]), static_cast<vtkIdType>(rMesh.indices[i + 2]) };
            triangles->InsertNextCell(3, pts);
        }

        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points); polyData->SetPolys(triangles);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);
        mapper->ScalarVisibilityOff(); 

        vtkSmartPointer<vtkActor> meshActor = vtkSmartPointer<vtkActor>::New();
        meshActor->SetMapper(mapper);
        
        // 【专业级渲染设置】：渲染为表面实体，并开启边线显示
        meshActor->GetProperty()->SetRepresentationToSurface();
        meshActor->GetProperty()->EdgeVisibilityOn(); 
        meshActor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); // 黑色网格边框
        meshActor->GetProperty()->SetLineWidth(1.0); 
        meshActor->GetProperty()->SetLighting(true);           // 开启曲面光照立体感

        m_renderer->AddActor(meshActor);
        m_idToMeshActor[faceId] = meshActor;
        
        // 【关键】：将网格 Actor 注册到映射表！这样网格不透明覆盖原 CAD 时，鼠标依然能拾取面！
        m_actorToFaceId[meshActor.Get()] = faceId;
    }

    updateGridExtent();

    updateFaceColors();
}

// ====================================================================
// 将 VTK 高阶单元 (Lagrange/Quadratic) 转化为标准线性单元,
// 避免 VTK 7.x 的 vtkUnstructuredGridGeometryFilter 报 UnknownClass 错误
// ====================================================================
static vtkSmartPointer<vtkUnstructuredGrid> ConvertHighOrderCellsToLinear(vtkUnstructuredGrid* grid) {
    bool needsConversion = false;
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++) {
        int t = grid->GetCellType(i);
        if (t >= 68 && t <= 74) { needsConversion = true; break; }
    }
    if (!needsConversion) return nullptr;

    vtkSmartPointer<vtkUnstructuredGrid> result = vtkSmartPointer<vtkUnstructuredGrid>::New();
    result->SetPoints(grid->GetPoints());
    result->Allocate(grid->GetNumberOfCells());

    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++) {
        int cellType = grid->GetCellType(i);
        vtkIdType npts;
        vtkIdType* pts;
        grid->GetCellPoints(i, npts, pts);

        switch (cellType) {
            case 68: result->InsertNextCell(VTK_LINE,     2, pts); break; // LAGRANGE_CURVE
            case 69: result->InsertNextCell(VTK_TRIANGLE, 3, pts); break; // LAGRANGE_TRIANGLE
            case 70: result->InsertNextCell(VTK_QUAD,     4, pts); break; // LAGRANGE_QUADRILATERAL
            case 71: result->InsertNextCell(VTK_TETRA,    4, pts); break; // LAGRANGE_TETRAHEDRON
            case 72: result->InsertNextCell(VTK_HEXAHEDRON,8, pts); break; // LAGRANGE_HEXAHEDRON
            case 73: result->InsertNextCell(VTK_WEDGE,    6, pts); break; // LAGRANGE_WEDGE
            case 74: result->InsertNextCell(VTK_PYRAMID,  5, pts); break; // LAGRANGE_PYRAMID
            default: result->InsertNextCell(cellType, npts, pts); break;
        }
    }

    result->GetPointData()->ShallowCopy(grid->GetPointData());
    result->GetCellData()->ShallowCopy(grid->GetCellData());
    return result;
}

// ====================================================================
// COMSOL 风格表面云图渲染 (Quality / Smoothing / 边界选择 可调)
//   quality  : 1.0 = 最高分辨率,  0.0 = 最低 (控制 DecimatePro)
//   smoothing: 0.0 = 不光滑,      1.0 = 最光滑 (控制 WindowedSinc)
// ====================================================================
void CADViewer::loadVtkContourSurface(const QString& vtkFilePath, const QString& fieldName,
                                       const QString& mshFilePath,
                                       const std::set<int>& selectedGroupTags,
                                       double quality, double smoothing,
                                       const std::set<int>& directEntityTags) {
    if (m_vtkResultActor) { m_renderer->RemoveActor(m_vtkResultActor); m_vtkResultActor = nullptr; }
    if (m_scalarBar)      { m_renderer->RemoveActor2D(m_scalarBar); m_scalarBar = nullptr; }

    // --------------- 1. 读取 VTK 体积网格 (获取节点坐标 + 标量数据) ---------------
    vtkSmartPointer<vtkAlgorithm> reader;
    if (vtkFilePath.endsWith(".pvtu", Qt::CaseInsensitive)) {
        auto r = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
        r->SetFileName(vtkFilePath.toStdString().c_str());
        r->Update();
        reader = r;
    } else if (vtkFilePath.endsWith(".vtu", Qt::CaseInsensitive)) {
        auto r = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        r->SetFileName(vtkFilePath.toStdString().c_str());
        r->Update();
        reader = r;
    } else {
        return;
    }

    vtkUnstructuredGrid* volGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputDataObject(0));
    if (!volGrid || volGrid->GetNumberOfPoints() == 0) return;
    auto linearGrid = ConvertHighOrderCellsToLinear(volGrid);
    if (linearGrid) volGrid = linearGrid;

    // 标量
    std::string fn = fieldName.toStdString();
    vtkDataArray* volScalars = volGrid->GetPointData()->GetArray(fn.c_str());
    if (!volScalars) {
        qWarning("loadVtkContourSurface: field '%s' not in VTK", qPrintable(fieldName));
        return;
    }

    // 用 VTK 节点建 KdTree (用于将网格面节点坐标映射到 VTK 标量)
    auto vtkPts = vtkSmartPointer<vtkPolyData>::New();
    vtkPts->SetPoints(volGrid->GetPoints());
    auto ptLoc = vtkSmartPointer<vtkKdTreePointLocator>::New();
    ptLoc->SetDataSet(vtkPts);
    ptLoc->BuildLocator();

    // --------------- 2. 从 .msh 文件读取选中物理组的 2D 边界单元 ---------------
    vtkSmartPointer<vtkPolyData> workingSurface;
    if (!mshFilePath.isEmpty() && !selectedGroupTags.empty()) {
        auto renderPts = vtkSmartPointer<vtkPoints>::New();
        auto renderCells = vtkSmartPointer<vtkCellArray>::New();
        auto renderScalars = vtkSmartPointer<vtkFloatArray>::New();
        renderScalars->SetName(fn.c_str());
        renderScalars->SetNumberOfComponents(1);

        vtkIdType nVtkPts = volGrid->GetNumberOfPoints();
        vtkIdType ptOff = 0;
        bool gmshOk = false;

        try {
            try { gmsh::initialize(); } catch (...) {}
            gmsh::open(mshFilePath.toStdString());

            // 一次性获取整个网格的节点坐标 (tag → xyz)
            std::vector<std::size_t> allNodeTags;
            std::vector<double> allCoords, allParamCoord;
            gmsh::model::mesh::getNodes(allNodeTags, allCoords, allParamCoord);
            std::map<std::size_t, int> tagToCoordIdx;
            for (size_t ni = 0; ni < allNodeTags.size(); ni++)
                tagToCoordIdx[allNodeTags[ni]] = static_cast<int>(ni);

            // 逐个物理组收集 2D 边界单元
            std::map<std::size_t, vtkIdType> tagToLocal; // global node tag → render point ID
            // 确定要查询的 entity tags: 优先用 directEntityTags, 否则通过 PG 查询
            std::set<int> entityTagsToRender;
            if (!directEntityTags.empty()) {
                entityTagsToRender = directEntityTags;
            } else {
                for (int groupTag : selectedGroupTags) {
                    std::vector<int> physEntityTags;
                    gmsh::model::getEntitiesForPhysicalGroup(2, groupTag, physEntityTags);
                    for (int et : physEntityTags) entityTagsToRender.insert(et);
                }
            }
            for (int entityTag : entityTagsToRender) {
                    std::vector<int> elemTypes;
                    std::vector<std::vector<std::size_t>> elemTags, nodeTagsPerElem;
                    gmsh::model::mesh::getElements(elemTypes, elemTags, nodeTagsPerElem,
                                                   2, entityTag);
                for (size_t ei = 0; ei < elemTypes.size(); ei++) {
                    int etype = elemTypes[ei];
                    int nNodesPerElem = (etype == 2) ? 3 : (etype == 3) ? 4 : 0;
                    if (nNodesPerElem == 0) continue;

                    const auto& nodes = nodeTagsPerElem[ei];
                    int nElems = static_cast<int>(nodes.size() / nNodesPerElem);

                    // 为每个单元顶点分配渲染点 (首次遇到时从全局坐标表获取 + KdTree 标量查值)
                    for (int k = 0; k < nElems * nNodesPerElem; k++) {
                        std::size_t tag = nodes[k];
                        if (tagToLocal.count(tag)) continue;

                        auto ci = tagToCoordIdx.find(tag);
                        if (ci == tagToCoordIdx.end()) continue;
                        int idx = ci->second;
                        double p[3] = {allCoords[3*idx], allCoords[3*idx+1], allCoords[3*idx+2]};
                        renderPts->InsertNextPoint(p);
                        vtkIdType nearest = ptLoc->FindClosestPoint(p);
                        float val = 0.0f;
                        if (nearest >= 0 && nearest < nVtkPts)
                            val = static_cast<float>(volScalars->GetComponent(nearest, 0));
                        renderScalars->InsertNextValue(val);
                        tagToLocal[tag] = ptOff++;
                    }
                    // 构建单元
                    for (int k = 0; k < nElems; k++) {
                        vtkIdType ids[4];
                        bool valid = true;
                        for (int j = 0; j < nNodesPerElem; j++) {
                            auto it = tagToLocal.find(nodes[k * nNodesPerElem + j]);
                            if (it == tagToLocal.end()) { valid = false; break; }
                            ids[j] = it->second;
                        }
                        if (!valid) continue;
                        if (nNodesPerElem == 3)
                            renderCells->InsertNextCell(3, ids);
                        else
                            renderCells->InsertNextCell(4, ids);
                    }
                }
            }
            gmsh::clear();
            gmshOk = true;
        } catch (std::exception& e) {
            qWarning("loadVtkContourSurface: Gmsh error: %s", e.what());
            try { gmsh::clear(); } catch (...) {}
        }

        if (gmshOk && renderCells->GetNumberOfCells() > 0) {
            auto pd = vtkSmartPointer<vtkPolyData>::New();
            pd->SetPoints(renderPts);
            pd->SetPolys(renderCells);
            pd->GetPointData()->AddArray(renderScalars);
            pd->GetPointData()->SetActiveScalars(fn.c_str());

            auto faceNorms = vtkSmartPointer<vtkPolyDataNormals>::New();
            faceNorms->SetInputData(pd);
            faceNorms->SplittingOff();
            faceNorms->ConsistencyOn();
            faceNorms->ComputePointNormalsOn();
            faceNorms->Update();
            workingSurface = faceNorms->GetOutput();
            if (!workingSurface || workingSurface->GetNumberOfCells() == 0)
                workingSurface = pd;
        }
    }

    // 无边界选择: 提取 VTK 体网格外表面
    if (!workingSurface || workingSurface->GetNumberOfCells() == 0) {
        auto surfFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfFilter->SetInputData(volGrid);
        surfFilter->Update();
        if (surfFilter->GetOutput() && surfFilter->GetOutput()->GetNumberOfCells() > 0) {
            auto norms = vtkSmartPointer<vtkPolyDataNormals>::New();
            norms->SetInputConnection(surfFilter->GetOutputPort());
            norms->SplittingOff();
            norms->ConsistencyOn();
            norms->ComputePointNormalsOn();
            norms->Update();
            workingSurface = norms->GetOutput();
        }
    }
    if (!workingSurface || workingSurface->GetNumberOfCells() == 0) return;

    // --------------- 3. 分辨率控制 ---------------
    if (quality < 0.99) {
        auto dec = vtkSmartPointer<vtkDecimatePro>::New();
        dec->SetInputData(workingSurface);
        dec->SetTargetReduction(1.0 - quality);
        dec->SetPreserveTopology(true);
        dec->SetFeatureAngle(60.0);
        dec->Update();
        workingSurface = dec->GetOutput();
    }

    // --------------- 4. 平滑控制 ---------------
    if (smoothing > 0.01) {
        auto smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
        smoother->SetInputData(workingSurface);
        int iter = static_cast<int>(smoothing * 40.0);
        double passBand = 0.1 - (smoothing * 0.099);
        if (passBand < 0.001) passBand = 0.001;
        smoother->SetNumberOfIterations(iter);
        smoother->SetPassBand(passBand);
        smoother->SetFeatureAngle(60.0);
        smoother->SetEdgeAngle(15.0);
        smoother->BoundarySmoothingOff();
        smoother->Update();
        workingSurface = smoother->GetOutput();
    }

    // --------------- 5. 标量场验证 + 值域 ---------------
    vtkPointData* pd = workingSurface->GetPointData();
    vtkDataArray* scalars = pd->GetArray(fn.c_str());
    if (!scalars) {
        scalars = workingSurface->GetCellData()->GetArray(fn.c_str());
        if (scalars) workingSurface->GetCellData()->SetActiveScalars(fn.c_str());
    }
    if (!scalars) {
        QStringList avail;
        for (int i = 0; i < pd->GetNumberOfArrays(); i++)
            avail << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "(null)");
        qWarning("loadVtkContourSurface: field '%s' missing. Available: %s",
                 qPrintable(fieldName), qPrintable(avail.join(", ")));
        return;
    }
    pd->SetActiveScalars(fn.c_str());

    double range[2];
    scalars->GetRange(range);
    if (range[1] - range[0] < 1e-12) {
        double c = range[0];
        range[0] = c - 0.5;
        range[1] = c + 0.5;
    }

    // --------------- 6. 256 阶 COMSOL 彩虹 LUT ---------------
    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToHSV();
    for (int i = 0; i <= 256; i++) {
        float t = i / 256.0f;
        double val = range[0] + t * (range[1] - range[0]);
        float r, g, b;
        if      (t < 0.125f) { float s = t / 0.125f;       r=0.231f*(1-s)+0.0f*s;   g=0.298f*(1-s)+0.502f*s; b=0.753f*(1-s)+0.753f*s; }
        else if (t < 0.375f) { float s=(t-0.125f)/0.25f;   r=0.0f;     g=0.502f*(1-s)+0.706f*s; b=0.753f*(1-s)+0.804f*s; }
        else if (t < 0.500f) { float s=(t-0.375f)/0.125f;  r=0.0f;     g=0.706f*(1-s)+0.812f*s; b=0.804f*(1-s)+0.0f*s;   }
        else if (t < 0.625f) { float s=(t-0.500f)/0.125f;  r=0.0f*(1-s)+0.655f*s; g=0.812f*(1-s)+0.89f*s; b=0.0f; }
        else if (t < 0.875f) { float s=(t-0.625f)/0.25f;   r=0.655f*(1-s)+1.0f*s; g=0.89f*(1-s)+0.706f*s;  b=0.0f; }
        else                 { float s=(t-0.875f)/0.125f;   r=1.0f*(1-s)+0.906f*s; g=0.706f*(1-s)+0.0f*s;   b=0.0f; }
        ctf->AddRGBPoint(val, r, g, b);
    }

    // --------------- 7. Mapper + Actor ---------------
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(workingSurface);
    mapper->SetScalarRange(range[0], range[1]);
    mapper->SetLookupTable(ctf);
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(fn.c_str());
    mapper->ScalarVisibilityOn();
    mapper->InterpolateScalarsBeforeMappingOn();

    m_vtkResultActor = vtkSmartPointer<vtkActor>::New();
    m_vtkResultActor->SetMapper(mapper);
    m_vtkResultActor->GetProperty()->SetLighting(false);
    m_vtkResultActor->GetProperty()->SetEdgeVisibility(false);
    m_vtkResultActor->SetVisibility(true);
    m_renderer->AddActor(m_vtkResultActor);

    // --------------- 8. 标量栏 ---------------
    // 标题带单位
    std::string title;
    if (fn == "rho") title = "rho (C/m^3/ε₀)";
    else if (fn == "E_scalar") title = "|E| (V/m)";
    else if (fn == "phi") title = "φ (V)";
    else title = fn;

    m_scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    m_scalarBar->SetLookupTable(ctf);
    m_scalarBar->SetTitle(title.c_str());
    m_scalarBar->SetNumberOfLabels(5);
    m_scalarBar->SetBarRatio(0.15);
    m_scalarBar->SetPosition(0.85, 0.1);
    m_scalarBar->SetWidth(0.1);
    m_scalarBar->SetHeight(0.6);
    // 标题: 黑色, 小字体
    m_scalarBar->GetTitleTextProperty()->SetColor(0, 0, 0);
    m_scalarBar->GetTitleTextProperty()->SetFontSize(10);
    m_scalarBar->GetTitleTextProperty()->BoldOff();
    m_scalarBar->GetTitleTextProperty()->ShadowOff();
    // 刻度标签: 黑色, 小字体
    m_scalarBar->GetLabelTextProperty()->SetColor(0, 0, 0);
    m_scalarBar->GetLabelTextProperty()->SetFontSize(6);
    m_scalarBar->GetLabelTextProperty()->ShadowOff();
    m_renderer->AddActor2D(m_scalarBar);

    m_geometryHidden = true;
    m_resultVisible = true;
    updateFaceColors();
    resetCameraView();
}

// ====================================================================
// 三维切面渲染: vtkCutter + vtkPlane 对体网格切平面, 插值标量后在片段着色器查色
// ====================================================================
void CADViewer::loadVtkSliceView(const QString& vtkFilePath, const QString& fieldName,
                                  bool sx, double px, bool sy, double py, bool sz, double pz) {
    if (!sx && !sy && !sz) return;
    if (m_vtkResultActor) { m_renderer->RemoveActor(m_vtkResultActor); m_vtkResultActor = nullptr; }
    if (m_scalarBar)      { m_renderer->RemoveActor2D(m_scalarBar); m_scalarBar = nullptr; }

    // 1. 读取 VTK 体网格
    vtkSmartPointer<vtkAlgorithm> reader;
    if (vtkFilePath.endsWith(".pvtu", Qt::CaseInsensitive)) {
        auto r = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
        r->SetFileName(vtkFilePath.toStdString().c_str());
        r->Update();
        reader = r;
    } else if (vtkFilePath.endsWith(".vtu", Qt::CaseInsensitive)) {
        auto r = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
         r->SetFileName(vtkFilePath.toStdString().c_str());
        r->Update();
        reader = r;
    } else return;

    vtkUnstructuredGrid* volGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutputDataObject(0));
    if (!volGrid || volGrid->GetNumberOfPoints() == 0) return;
    auto linearGrid = ConvertHighOrderCellsToLinear(volGrid);
    if (linearGrid) volGrid = linearGrid;

    // 2. 对各启用轴用 vtkCutter 切平面
    auto appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    auto doSlice = [&](double origin[3], double normal[3]) {
        auto plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetOrigin(origin);
        plane->SetNormal(normal);
        auto cutter = vtkSmartPointer<vtkCutter>::New();
        cutter->SetCutFunction(plane);
        cutter->SetInputData(volGrid);
        cutter->Update();
        if (cutter->GetOutput() && cutter->GetOutput()->GetNumberOfCells() > 0)
            appendFilter->AddInputData(cutter->GetOutput());
    };

    if (sx) { double o[3]={px,0,0}, n[3]={1,0,0}; doSlice(o,n); }
    if (sy) { double o[3]={0,py,0}, n[3]={0,1,0}; doSlice(o,n); }
    if (sz) { double o[3]={0,0,pz}, n[3]={0,0,1}; doSlice(o,n); }

    appendFilter->Update();
    if (!appendFilter->GetOutput() || appendFilter->GetOutput()->GetNumberOfCells() == 0) {
        qWarning("loadVtkSliceView: 切面与网格无交集, 请调整切面位置坐标");
        return;
    }

    // 3. 法向
    auto norms = vtkSmartPointer<vtkPolyDataNormals>::New();
    norms->SetInputConnection(appendFilter->GetOutputPort());
    norms->SplittingOff();
    norms->ConsistencyOn();
    norms->ComputePointNormalsOn();
    norms->Update();
    vtkPolyData* sliceSurface = norms->GetOutput();

    // 4. 标量场
    std::string fn = fieldName.toStdString();
    vtkDataArray* scalars = sliceSurface->GetPointData()->GetArray(fn.c_str());
    if (!scalars) {
        scalars = sliceSurface->GetCellData()->GetArray(fn.c_str());
        if (scalars) sliceSurface->GetCellData()->SetActiveScalars(fn.c_str());
    }
    if (!scalars) {
        QStringList avail;
        auto* pd = sliceSurface->GetPointData();
        for (int i = 0; i < pd->GetNumberOfArrays(); i++)
            avail << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "(null)");
        qWarning("loadVtkSliceView: field '%s' missing. Available: %s",
                 qPrintable(fieldName), qPrintable(avail.join(", ")));
        return;
    }
    sliceSurface->GetPointData()->SetActiveScalars(fn.c_str());

    double range[2]; scalars->GetRange(range);
    if (range[1] - range[0] < 1e-12) { double c = range[0]; range[0] = c - 0.5; range[1] = c + 0.5; }

    // 5. 256 阶彩虹 LUT
    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToHSV();
    for (int i = 0; i <= 256; i++) {
        float t = i / 256.0f; double val = range[0] + t * (range[1] - range[0]);
        float r,g,b;
        if      (t < 0.125f) { float s=t/0.125f;       r=0.231f*(1-s)+0.0f*s;   g=0.298f*(1-s)+0.502f*s; b=0.753f*(1-s)+0.753f*s; }
        else if (t < 0.375f) { float s=(t-0.125f)/0.25f;r=0.0f;     g=0.502f*(1-s)+0.706f*s; b=0.753f*(1-s)+0.804f*s; }
        else if (t < 0.500f) { float s=(t-0.375f)/0.125f;r=0.0f;   g=0.706f*(1-s)+0.812f*s; b=0.804f*(1-s)+0.0f*s; }
        else if (t < 0.625f) { float s=(t-0.500f)/0.125f;r=0.0f*(1-s)+0.655f*s; g=0.812f*(1-s)+0.89f*s; b=0.0f; }
        else if (t < 0.875f) { float s=(t-0.625f)/0.25f; r=0.655f*(1-s)+1.0f*s; g=0.89f*(1-s)+0.706f*s; b=0.0f; }
        else                 { float s=(t-0.875f)/0.125f;r=1.0f*(1-s)+0.906f*s; g=0.706f*(1-s)+0.0f*s;   b=0.0f; }
        ctf->AddRGBPoint(val, r, g, b);
    }

    // 6. Mapper + Actor
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(sliceSurface);
    mapper->SetScalarRange(range[0], range[1]);
    mapper->SetLookupTable(ctf);
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(fn.c_str());
    mapper->ScalarVisibilityOn();
    mapper->InterpolateScalarsBeforeMappingOn();

    m_vtkResultActor = vtkSmartPointer<vtkActor>::New();
    m_vtkResultActor->SetMapper(mapper);
    m_vtkResultActor->GetProperty()->SetLighting(false);
    m_vtkResultActor->SetVisibility(true);
    m_renderer->AddActor(m_vtkResultActor);

    // 7. 标量栏
    std::string title;
    if (fn == "rho") title = "rho (C/m^3/ε₀)";
    else if (fn == "E_scalar") title = "|E| (V/m)";
    else if (fn == "phi") title = "φ (V)";
    else title = fn;
    m_scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    m_scalarBar->SetLookupTable(ctf);
    m_scalarBar->SetTitle(title.c_str());
    m_scalarBar->SetNumberOfLabels(5);
    m_scalarBar->SetBarRatio(0.15);
    m_scalarBar->SetPosition(0.85, 0.1);
    m_scalarBar->SetWidth(0.1);
    m_scalarBar->SetHeight(0.6);
    m_scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
    m_scalarBar->GetTitleTextProperty()->SetFontSize(10);
    m_scalarBar->GetTitleTextProperty()->BoldOff();
    m_scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
    m_scalarBar->GetLabelTextProperty()->SetFontSize(6);
    m_renderer->AddActor2D(m_scalarBar);

    updateGridExtent();

    m_geometryHidden = true;
    m_resultVisible = true;
    updateFaceColors();
    resetCameraView();
}

void CADViewer::findClosestEdge(int mx, int my, int& objId, int& edgeIdx, TopoDS_Edge& edge, gp_Pnt& pickPt) {
    objId = -1; edgeIdx = -1;
    fprintf(stderr, "[DEBUG findClosestEdge] searching at (%d,%d), m_faceToObject.size()=%zu\n", mx, my, m_faceToObject.size());
    // 屏幕空间距离: 将各边投影到屏幕, 计算到鼠标的 2D 像素距离
    vtkSmartPointer<vtkCoordinate> coord = vtkSmartPointer<vtkCoordinate>::New();
    coord->SetViewport(m_renderer);
    double bestPxDist = 12.0;  // 12 像素内有效
    for (auto& [fid, oid] : m_faceToObject) {
        if (!m_faceIdToShape.count(fid)) continue;
        const TopoDS_Shape& sh = m_faceIdToShape[fid];
        int ei = 0;
        for (TopExp_Explorer ex(sh, TopAbs_EDGE); ex.More(); ex.Next(), ei++) {
            TopoDS_Edge E = TopoDS::Edge(ex.Current());
            try {
                double u1, u2;
                Handle(Geom_Curve) curve = BRep_Tool::Curve(E, u1, u2);
                BRepAdaptor_Curve adapt(E);
                double lu1 = adapt.FirstParameter(), lu2 = adapt.LastParameter();
                // 沿曲线采样, 投影到屏幕, 计算像素距离
                double minPx = 1e30;
                int nSamples = 64;
                for (int s = 0; s <= nSamples; s++) {
                    double u = lu1 + (lu2 - lu1) * s / (double)nSamples;
                    gp_Pnt pt; adapt.D0(u, pt);
                    coord->SetValue(pt.X(), pt.Y(), pt.Z());
                    int* dp = coord->GetComputedDisplayValue(m_renderer);
                    if (!dp) continue;
                    double dx = dp[0] - mx, dy = dp[1] - my;
                    double pxDist = std::sqrt(dx*dx + dy*dy);
                    if (pxDist < minPx) minPx = pxDist;
                }
                if (minPx < bestPxDist) {
                    bestPxDist = minPx;
                    objId = oid; edgeIdx = ei; edge = E;
                    // 近似 3D 点用于其余逻辑
                    double uMid = lu1 + (lu2 - lu1) * 0.5;
                    gp_Pnt pm; adapt.D0(uMid, pm); pickPt = pm;
                }
            } catch (...) {}
        }
    }
}
std::vector<std::pair<int,int>> CADViewer::getSelectedEdges() const {
    std::vector<std::pair<int,int>> result;
    for (auto& [key, val] : m_selectedEdgeMap)
        result.push_back({key.objId, key.edgeIdx});
    return result;
}

void CADViewer::clearEdgePreview() { if (m_previewEdge) { m_renderer->RemoveActor(m_previewEdge); m_previewEdge=nullptr; } }
void CADViewer::clearAllEdgeHighlights() {
    clearEdgePreview();
    for (auto& [k, val] : m_selectedEdgeMap) m_renderer->RemoveActor(val.first);
    m_selectedEdgeMap.clear();
    if (m_highlightedEdge) { m_renderer->RemoveActor(m_highlightedEdge); m_highlightedEdge=nullptr; }
    GetRenderWindow()->Render();
}

void CADViewer::setEdgeHighlightsVisible(bool visible) {
    for (auto& [key, val] : m_selectedEdgeMap)
        val.first->SetVisibility(visible ? 1 : 0);
    GetRenderWindow()->Render();
}

void CADViewer::highlightEdgeInShape(const TopoDS_Shape& shape, int edgeIdx, int objId) {
    if (shape.IsNull()) return;
    int cnt = 0;
    for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next(), cnt++) {
        if (cnt == edgeIdx) {
            TopoDS_Edge E = TopoDS::Edge(ex.Current());
            EdgeKey key{objId, edgeIdx};
            if (m_selectedEdgeMap.count(key)) return; // already highlighted
            auto a = makeEdgeActor(E, 0.2, 0.8, 1.0, 4.0);
            m_renderer->AddActor(a);
            m_selectedEdgeMap[key] = {a, E};
            return;
        }
    }
}

// 边高亮辅助
vtkSmartPointer<vtkActor> makeEdgeActor(const TopoDS_Edge& E, double r, double g, double b, double w) {
    auto pts=vtkSmartPointer<vtkPoints>::New(); auto cells=vtkSmartPointer<vtkCellArray>::New();
    BRepAdaptor_Curve curve(E);
    double u1=curve.FirstParameter(), u2=curve.LastParameter();
    int n=21; cells->InsertNextCell(n);
    for (int i=0; i<n; i++) {
        double u=u1+(u2-u1)*i/(n-1.0); gp_Pnt pt; curve.D0(u,pt);
        pts->InsertNextPoint(pt.X(),pt.Y(),pt.Z()); cells->InsertCellPoint(i);
    }
    auto pd=vtkSmartPointer<vtkPolyData>::New(); pd->SetPoints(pts); pd->SetLines(cells);
    auto m=vtkSmartPointer<vtkPolyDataMapper>::New(); m->SetInputData(pd);
    auto a=vtkSmartPointer<vtkActor>::New(); a->SetMapper(m);
    a->GetProperty()->SetColor(r,g,b); a->GetProperty()->SetLineWidth(w); a->GetProperty()->SetLighting(false);
    a->PickableOff();  // 边高亮不参与拾取, 点击穿透到下面的面
    return a;
}

// ====================================================================
// COMSOL 风格栅格 (工作平面参考线)
// ====================================================================
void CADViewer::setGridVisible(bool v) {
    m_gridVisible = v;
    if (v) rebuildGridActor();   // 确保 mapper 已构建 (首次点击/参数修改后)
    for (auto& gp : m_gridPlanes) {
        gp.majorActor->SetVisibility(v ? 1 : 0);
        gp.minorActor->SetVisibility(v ? 1 : 0);
        for (auto& la : gp.labels) la->SetVisibility(v ? 1 : 0);
    }
    GetRenderWindow()->Render();
}

void CADViewer::setGridSpacing(double major) {
    if (major <= 0) return;
    m_gridSpacing = major;
    if (m_gridVisible) rebuildGridActor();
}

void CADViewer::setGridSubdivisions(int n) {
    if (n < 1 || n > 10) return;
    m_gridSubdiv = n;
    if (m_gridVisible) rebuildGridActor();
}

void CADViewer::setGridAutoExtent(bool on) {
    m_gridAutoExtent = on;
    if (m_gridVisible) { updateGridExtent(); rebuildGridActor(); }
}

// 从可见包围盒计算三个面的栅格范围 (自动模式)
// 平面 0 (XY): z=zmin 底面, 平面 1 (XZ): y=ymin 背面, 平面 2 (YZ): x=xmin 侧面
// 显式遍历几何/网格/结果 actors 计算组合包围盒
// (不使用 ComputeVisiblePropBounds: 它受渲染器可见性状态影响, 且会把
//  栅格 actor 算进去, 导致范围停留在初始化值或无限膨胀)
bool CADViewer::computeSceneBounds(double b[6]) {
    bool first = true;
    auto merge = [&](vtkActor* actor) {
        if (!actor || !actor->GetVisibility()) return;
        double ab[6];
        actor->GetBounds(ab);
        if (ab[0] > ab[1] || ab[2] > ab[3] || ab[4] > ab[5]) return;  // 未初始化
        if (first) {
            for (int i = 0; i < 6; i++) b[i] = ab[i];
            first = false;
        } else {
            b[0] = std::min(b[0], ab[0]); b[1] = std::max(b[1], ab[1]);
            b[2] = std::min(b[2], ab[2]); b[3] = std::max(b[3], ab[3]);
            b[4] = std::min(b[4], ab[4]); b[5] = std::max(b[5], ab[5]);
        }
    };
    for (auto& [fid, actor] : m_idToActor)      merge(actor);  // 几何面
    for (auto& [fid, actor] : m_idToMeshActor)  merge(actor);  // 网格面
    if (m_vtkResultActor) merge(m_vtkResultActor.Get());       // 结果
    return !first;
}

void CADViewer::updateGridExtent() {
    if (!m_gridAutoExtent) return;   // 手动模式: 不自动更新

    double b[6];
    if (!computeSceneBounds(b)) {
        fprintf(stderr, "[Grid] 无可见几何, 保留当前范围\n"); fflush(stderr);
        return;   // 无几何时保留默认范围
    }

    const double pad = 0.0;          // 不外扩: 刻度首尾精确 = 几何 min/max
    double sx = b[1] - b[0], sy = b[3] - b[2], sz = b[5] - b[4];
    if (sx < 1e-12) sx = 1; if (sy < 1e-12) sy = 1; if (sz < 1e-12) sz = 1;
    double offset = std::max(sx, std::max(sy, sz)) / 10.0;   // 悬浮距离 = 包围盒最大边长 / 10

    // 三个互不相邻的面, 避免交叉:
    //   XY 面: z=zmin (底面)     u=x, v=y
    //   XZ 面: y=ymax (顶背面)   u=x, v=z
    //   YZ 面: x=xmax (右侧面)   u=y, v=z
    m_gridRange[0][0] = b[0] - pad * sx;  m_gridRange[0][1] = b[1] + pad * sx;
    m_gridRange[0][2] = b[2] - pad * sy;  m_gridRange[0][3] = b[3] + pad * sy;
    m_gridPlanePos[0] = b[4]-offset;   // zmin

    m_gridRange[1][0] = b[0] - pad * sx;  m_gridRange[1][1] = b[1] + pad * sx;
    m_gridRange[1][2] = b[4] - pad * sz;  m_gridRange[1][3] = b[5] + pad * sz;
    m_gridPlanePos[1] = b[3]+offset;   // ymax

    m_gridRange[2][0] = b[2] - pad * sy;  m_gridRange[2][1] = b[3] + pad * sy;
    m_gridRange[2][2] = b[4] - pad * sz;  m_gridRange[2][3] = b[5] + pad * sz;
    m_gridPlanePos[2] = b[1]+offset;   // xmax

    fprintf(stderr, "[Grid] bounds=(%.2f,%.2f, %.2f,%.2f, %.2f,%.2f) XYrange=(%.1f,%.1f)\n",
            b[0], b[1], b[2], b[3], b[4], b[5], m_gridRange[0][0], m_gridRange[0][1]);
    fflush(stderr);

    if (m_gridVisible) rebuildGridActor();
}

// 重建三个包围盒面的主/次栅格 polydata + 刻度标签
void CADViewer::rebuildGridActor() {
    // 数字格式: 按主间距数量级自适应小数位
    std::string fmt;
    double a = std::fabs(m_gridSpacing);
    if (a >= 100.0)      fmt = "%.0f";
    else if (a >= 1.0)   fmt = "%.1f";
    else if (a >= 0.01)  fmt = "%.2f";
    else                 fmt = "%.4g";

    for (int p = 0; p < 3; p++) {
        double u0 = m_gridRange[p][0], u1 = m_gridRange[p][1];
        double v0 = m_gridRange[p][2], v1 = m_gridRange[p][3];
        double w  = m_gridPlanePos[p];
        double major = m_gridSpacing;
        double minor = major / std::max(1, m_gridSubdiv);

        // 线数上限保护 (防止超大网格导致卡顿)
        auto nLines = [&](double step) {
            int nu = (int)std::floor((u1 - u0) / step + 1e-9) + 1;
            int nv = (int)std::floor((v1 - v0) / step + 1e-9) + 1;
            return nu + nv;
        };
        if (nLines(major) > 500) {   // 间距过密 → 按 500 条钳制
            double span = std::max(u1 - u0, v1 - v0);
            major = span / 250.0;
            minor = major / std::max(1, m_gridSubdiv);
        }
        if (nLines(minor) > 2000) minor = std::max(major / 2, (u1 - u0) / 1000.0);

        // 3D 坐标组装: 平面 (u,v) → (x,y,z)
        // p=0 (XY): (u,v,w); p=1 (XZ): (u,w,v); p=2 (YZ): (w,u,v)
        auto makePt = [&](double u, double v) -> std::array<double,3> {
            if (p == 0) return {u, v, w};
            if (p == 1) return {u, w, v};
            return {w, u, v};
        };

        GridPlaneSet& gp = m_gridPlanes[p];

        // ---- 主栅格线 ----
        auto ptsM = vtkSmartPointer<vtkPoints>::New();
        auto linesM = vtkSmartPointer<vtkCellArray>::New();
        vtkIdType nid = 0;
        auto addLine = [&](double ua, double va, double ub, double vb) {
            auto pa = makePt(ua, va), pb = makePt(ub, vb);
            ptsM->InsertNextPoint(pa[0], pa[1], pa[2]);
            ptsM->InsertNextPoint(pb[0], pb[1], pb[2]);
            vtkIdType ids[2] = {nid, nid + 1}; nid += 2;
            linesM->InsertNextCell(2, ids);
        };
        for (double u = u0; u <= u1 + 1e-9; u += major) addLine(u, v0, u, v1);
        for (double v = v0; v <= v1 + 1e-9; v += major) addLine(u0, v, u1, v);
        auto pdM = vtkSmartPointer<vtkPolyData>::New();
        pdM->SetPoints(ptsM); pdM->SetLines(linesM);
        auto mapperM = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperM->SetInputData(pdM);
        gp.majorActor->SetMapper(mapperM);

        // ---- 次栅格线 (细分) ----
        auto ptsS = vtkSmartPointer<vtkPoints>::New();
        auto linesS = vtkSmartPointer<vtkCellArray>::New();
        nid = 0;
        auto addLineS = [&](double ua, double va, double ub, double vb) {
            auto pa = makePt(ua, va), pb = makePt(ub, vb);
            ptsS->InsertNextPoint(pa[0], pa[1], pa[2]);
            ptsS->InsertNextPoint(pb[0], pb[1], pb[2]);
            vtkIdType ids[2] = {nid, nid + 1}; nid += 2;
            linesS->InsertNextCell(2, ids);
        };
        // 次栅格线避开主栅格线位置 (偏移半步长)
        double half = minor * 0.5;
        for (double u = u0 + half; u <= u1 + 1e-9; u += minor)
            if (std::abs(std::fmod((u - u0), major)) > 1e-9) addLineS(u, v0, u, v1);
        for (double v = v0 + half; v <= v1 + 1e-9; v += minor)
            if (std::abs(std::fmod((v - v0), major)) > 1e-9) addLineS(u0, v, u1, v);
        auto pdS = vtkSmartPointer<vtkPolyData>::New();
        pdS->SetPoints(ptsS); pdS->SetLines(linesS);
        auto mapperS = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperS->SetInputData(pdS);
        gp.minorActor->SetMapper(mapperS);

        // ---- 主刻度坐标标签 (栅格边缘外侧, 3D 文字参与深度测试, 被几何遮挡) ----
        for (auto& la : gp.labels) m_renderer->RemoveActor(la);
        gp.labels.clear();

        double off = (u1 - u0 + v1 - v0) * 0.015;   // 标签离边缘的偏移
        double fudge = off * 0.5;                    // 沿轴方向微移, 避免文字压线
        double th = std::max(u1 - u0, v1 - v0) * 0.03;   // 字高 (世界单位, VectorText 默认高 1.0)
        int nu = (int)std::floor((u1 - u0) / major + 1e-9) + 1;
        int nv = (int)std::floor((v1 - v0) / major + 1e-9) + 1;
        // 密度保护: 主刻度线过密时跳过标签
        if (nu + nv <= 80) {
            auto addLabel = [&](double u, double v, const char* text) {
                auto pa = makePt(u, v);
                // vtkVectorText + vtkFollower: 3D 文字, 面向相机, 参与深度测试
                // (几何后面的数字会被遮挡, 不再透过几何可见)
                auto textSrc = vtkSmartPointer<vtkVectorText>::New();
                textSrc->SetText(text);
                auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                mapper->SetInputConnection(textSrc->GetOutputPort());
                auto label = vtkSmartPointer<vtkFollower>::New();
                label->SetMapper(mapper);
                label->SetPosition(pa[0], pa[1], pa[2]);
                label->SetCamera(m_renderer->GetActiveCamera());
                label->GetProperty()->SetColor(0.12, 0.12, 0.12);
                label->GetProperty()->SetLineWidth(1.0);
                label->SetScale(th);
                label->SetPickable(false);
                label->SetVisibility(m_gridVisible ? 1 : 0);  // 跟随栅格状态, 重建后不消失
                m_renderer->AddActor(label);
                gp.labels.push_back(label);
            };

            char buf[64];
            // u 轴刻度标签: 常规刻度 (不含末端), 循环后强制补 u1(max); u0(min) 含在循环中
            for (double u = u0; u < u1 - 1e-9; u += major) {
                std::snprintf(buf, sizeof(buf), fmt.c_str(), u);
                addLabel(u + fudge, v0 - off, buf);
            }
            std::snprintf(buf, sizeof(buf), fmt.c_str(), u1);
            addLabel(u1 + fudge, v0 - off, buf);   // 补 max
            // v 轴刻度标签: 同理
            for (double v = v0; v < v1 - 1e-9; v += major) {
                std::snprintf(buf, sizeof(buf), fmt.c_str(), v);
                addLabel(u0 - off, v + fudge, buf);
            }
            std::snprintf(buf, sizeof(buf), fmt.c_str(), v1);
            addLabel(u0 - off, v1 + fudge, buf);   // 补 max
        }
    }

    GetRenderWindow()->Render();
}

void CADViewer::clearFaceGroup(int startFid) {
    // 只删除与 startFid 同属一个对象的面 (通过 m_faceToObject 判断)
    int targetObjId = -1;
    auto it = m_faceToObject.find(startFid);
    if (it != m_faceToObject.end()) targetObjId = it->second;
    if (targetObjId < 0) return;

    for (auto& [fid, objId] : m_faceToObject) {
        if (objId != targetObjId) continue;
        if (m_idToActor.count(fid)) {
            m_renderer->RemoveActor(m_idToActor[fid]);
            m_actorToFaceId.erase(m_idToActor[fid]);
            m_idToActor.erase(fid);
        }
        if (m_idToEdgeActor.count(fid)) {
            m_renderer->RemoveActor(m_idToEdgeActor[fid]);
            m_idToEdgeActor.erase(fid);
        }
        m_selectedFaceIds.erase(fid);
        m_hiddenFaceIds.erase(fid);
    }
    // 清除已遍历的面 (不能在遍历中直接erase)
    for (auto it2 = m_faceToObject.begin(); it2 != m_faceToObject.end(); ) {
        if (it2->second == targetObjId) it2 = m_faceToObject.erase(it2);
        else ++it2;
    }
    GetRenderWindow()->Render();
}

void CADViewer::setGeometryHidden(bool hidden) {
    if (m_geometryHidden != hidden) {
        m_geometryHidden = hidden;
        updateFaceColors();
    }
}

void CADViewer::setResultVisible(bool visible) {
    m_resultVisible = visible;
    if (m_vtkResultActor)  m_vtkResultActor->SetVisibility(visible ? 1 : 0);
    if (m_scalarBar)       m_scalarBar->SetVisibility(visible ? 1 : 0);
    this->GetRenderWindow()->Render();
}