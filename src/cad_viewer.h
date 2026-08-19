#pragma once

#include <QVTKWidget.h> 
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPropPicker.h> 
#include <vtkOrientationMarkerWidget.h>
#include <set>
#include <map>
#include <vector>

#include "mesh_data.h"

class vtkScalarBarActor;
class vtkFollower;
#include <TopoDS_Edge.hxx>
class TopoDS_Shape;
class gp_Pnt;

enum ViewMode { VIEW_NORMAL = 0, VIEW_HIDDEN_ONLY, VIEW_ALL };

class CADViewer : public QVTKWidget {
    Q_OBJECT
    friend class CustomTrackballStyle;
public:
    explicit CADViewer(QWidget *parent = nullptr);
    ~CADViewer();

    bool loadStepFile(const QString& filePath);
    bool loadIgesFile(const QString& filePath);
    bool loadBrepFile(const QString& filePath);
    int addOCCShape(const TopoDS_Shape& shape, double r=0.85, double g=0.88, double b=0.9);
    void setSelectedFaces(const std::set<int>& faceIds);
    void clearSelection();

    void setViewXY(); void setViewYZ(); void setViewZX(); void resetCameraView();

    void setInteractionModes(bool boxPick, bool encloseMode, bool hideMode);
    void resetHidden(); void showUnhidden(); void showOnlyHidden(); void showAll();
    void hideSelected();

    void loadMesh(const FaceMeshMap& meshData);
    void clearMesh();

    void clearScene();

    void loadVtkContourSurface(const QString& vtkFilePath, const QString& fieldName,
                               const QString& mshFilePath,
                               const std::set<int>& selectedGroupTags = {},
                               double quality = 1.0, double smoothing = 0.0,
                               const std::set<int>& directEntityTags = {});
    void loadVtkSliceView(const QString& vtkFilePath, const QString& fieldName,
                          bool sx, double px, bool sy, double py, bool sz, double pz);
    bool hasResultActor() const { return m_vtkResultActor != nullptr; }

    vtkRenderer* GetRenderer() { return m_renderer; }

    // ---- COMSOL 风格栅格 (三个面贴着几何包围盒) ----
    void setGridVisible(bool v);
    void setGridSpacing(double major);     // 主栅格间距
    void setGridSubdivisions(int n);       // 每主间距细分 (1~10)
    void setGridAutoExtent(bool on);       // 范围自动适配几何
    bool gridVisible() const { return m_gridVisible; }
    double gridSpacing() const { return m_gridSpacing; }
    int gridSubdivisions() const { return m_gridSubdiv; }

    void setObjectSelectMode(bool on) { m_objectSelectMode = on; }
    void setEdgeSelectMode(bool on) { m_edgeSelectMode = on; if (!on) m_highlightedEdge = nullptr; }
    const std::map<int,int>& faceObjectMap() const { return m_faceToObject; }
    void clearFaceGroup(int startFid);
    void setGeometryHidden(bool hidden);
    void setResultVisible(bool visible);

    void setGeometryTransparent(bool isTransparent);

    void processSingleClickSelection(int x, int y);
    void processBoxSelection(int xmin, int ymin, int xmax, int ymax);

signals:
    void faceSelectionChanged(const std::set<int>& selectedIds);
    void boxSelectionCompleted();

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;

private:
    void setTheme(bool isLight);
    void loadShape(const TopoDS_Shape& shape);
    void updateFaceColors();
    void toggleFaceSelection(int faceId);

    void updateGridExtent();               // 从可见包围盒计算栅格范围 (自动模式)
    void rebuildGridActor();               // 重建主/次栅格 polydata
    // 显式遍历几何/网格/结果 actors 计算组合包围盒 (不依赖渲染器可见性状态)
    bool computeSceneBounds(double b[6]);

    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkPropPicker> m_picker;
    vtkSmartPointer<vtkOrientationMarkerWidget> m_axesWidget;
    vtkSmartPointer<vtkActor> m_vtkResultActor;
    vtkSmartPointer<vtkScalarBarActor> m_scalarBar;

    // ---- COMSOL 风格栅格成员 (三面贴包围盒) ----
    struct GridPlaneSet {
        vtkSmartPointer<vtkActor> majorActor;               // 主栅格线
        vtkSmartPointer<vtkActor> minorActor;               // 次栅格线
        std::vector<vtkSmartPointer<vtkFollower>> labels;   // 刻度数字 (3D, 参与深度测试)
    };
    std::vector<GridPlaneSet> m_gridPlanes;  // [0]=XY(z=zmin) [1]=XZ(y=ymax) [2]=YZ(x=xmax) 互不相邻不交叉
    bool m_gridVisible = false;
    double m_gridSpacing = 10.0;
    int  m_gridSubdiv = 5;
    bool m_gridAutoExtent = true;
    double m_gridRange[3][4] = {{-50, 50, -50, 50}, {-50, 50, -50, 50}, {-50, 50, -50, 50}};  // 每平面两轴范围
    double m_gridPlanePos[3] = {0.0, 0.0, 0.0};   // 每平面固定轴位置 (zmin/ymin/xmin)

    std::map<vtkActor*, int> m_actorToFaceId;
    std::map<int, vtkActor*> m_idToActor;
    std::map<int, vtkSmartPointer<vtkActor>> m_idToEdgeActor; 
    std::map<int, vtkSmartPointer<vtkActor>> m_idToMeshActor; 

    std::set<int> m_selectedFaceIds;
    std::set<int> m_hiddenFaceIds;

    bool m_boxSelectMode = false;
    bool m_encloseMode = false;
    bool m_hideMode = false;
    
    bool m_isTransparent = false;
    bool m_geometryHidden = false;
    bool m_resultVisible = false;
    bool m_objectSelectMode = false;
    bool m_edgeSelectMode = false;
    vtkSmartPointer<vtkActor> m_highlightedEdge;
    vtkSmartPointer<vtkActor> m_previewEdge;
    struct EdgeKey { int objId; int edgeIdx; bool operator<(const EdgeKey& o) const { return objId!=o.objId?objId<o.objId:edgeIdx<o.edgeIdx; } };
    std::map<EdgeKey, std::pair<vtkSmartPointer<vtkActor>, TopoDS_Edge>> m_selectedEdgeMap;
    int m_lastObjId=-1, m_lastEdgeIdx=-1; 
  public:
    int lastPickedObjId() const { return m_lastObjId; } int lastPickedEdgeIdx() const { return m_lastEdgeIdx; }
    bool isEdgeSelected(int objId, int edgeIdx) const { return m_selectedEdgeMap.count(EdgeKey{objId, edgeIdx}); }
    void setEdgeHighlightsVisible(bool visible);
    void highlightEdgeInShape(const TopoDS_Shape& shape, int edgeIdx, int objId);
    void clearAllEdgeHighlights();
    std::vector<std::pair<int,int>> getSelectedEdges() const;
  private:
    void findClosestEdge(int mx,int my,int&objId,int&edgeIdx,TopoDS_Edge&edge,gp_Pnt&pt);
    void clearEdgePreview();

    int m_occFaceIdCounter = 100000;
    std::map<int, TopoDS_Shape> m_faceIdToShape;
    int m_lastEdgeNode = -1, m_lastEdgeLocalIdx = -1;
    std::map<int, int> m_faceToObject;  // faceId → objectId (用于对象选择)
    int m_viewMode = VIEW_NORMAL;
};