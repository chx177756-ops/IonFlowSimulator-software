#pragma once
#include <QMainWindow>
#include <map>
#include <set>
#include <QString>
#include <QTreeWidgetItem>
#include <QProcess>
#include <QTableWidget>
#include <QTextEdit>
#include <QComboBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QDoubleSpinBox>
#include <QProgressBar>
#include <QPushButton>
#include <QTimer>
#include <QDockWidget>
#include <QStatusBar>
#include <QListWidget>
#include <QStackedWidget>
#include <QLabel>
#include <QTabWidget>
#include <QToolButton>
#include <QLineEdit>
#include <atomic>
#include <thread>
#include <TopoDS_Shape.hxx>
#include "mesh_data.h"
#include "project_config.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void setProjectPath(const QString& path);

private slots:
    void importGeometry();
    void exportGeometry();
    void importExternalMesh();

    void showTreeContextMenu(const QPoint& pos);
    
    void onTreeSelectionChanged(); 
    void onFaceSelectionChanged(const std::set<int>& selectedIds);
    void onGroupNameEdited(const QString& newName); 
    
    void onMeshGroupTargetChanged(int index);
    void onEnableLocalFieldToggled(bool checked); 
    void onMeshParamsChanged(); 

    void generateMesh();
    void exportMesh();
    void exportMeshToFile(const QString& path);

private:
    void updateInteractionModes(); 
    bool runMeshEngine(const QString& savePath); 

    void setupGeometryToolbar();
    void setupSolverUI();
    void setupPostProcessUI();
    void setupBoundaryUI();
    void setupParamUI();
    void showCaseLibrary();
    void evalAllParams();
    void refreshParamTable(int pageIdx);
    void refreshBoundaryUI();
    void refreshSolverMeshCombo();
    void refreshPostProcessGroups(); // 【新增】：刷新后处理界面的物理组列表
    void generateJsonAndCalculate();
    // 对导入网格: 基于 m_groupData 同步物理组, 写出统一缓存路径 .msh + .pgroups.json
    // 返回同步后的 .msh 路径; 失败返回空串
    QString syncImportedMeshToCache();

    void logMessage(const QString& msg, const QString& level = "info");

    // 几何参数页
    QWidget* m_pageGeomParam = nullptr;
    // 求解器 UI 组件
    QWidget* m_pageSolver;
    QTextEdit* m_textConsole;
    QComboBox* m_comboSolverMesh;
    QDoubleSpinBox *m_spinE0 = nullptr, *m_spinRho = nullptr, *m_spinK = nullptr;
    QDoubleSpinBox *m_spinWindX = nullptr, *m_spinWindY = nullptr, *m_spinWindZ = nullptr;
    QDoubleSpinBox *m_spinWRho = nullptr, *m_spinGoalConv = nullptr;
    QDoubleSpinBox *m_spinTolE = nullptr, *m_spinTolRho = nullptr;
    QSpinBox* m_spinCores = nullptr;
    QComboBox* m_comboSolverType = nullptr;

    // 边界条件 UI 组件
    QComboBox* m_comboBoundaryGroup;
    QCheckBox* m_chkApplyBoundary;
    QDoubleSpinBox* m_spinBdrVoltage;
    QComboBox* m_cmbBdrMode = nullptr;
    QStackedWidget* m_stackBdrInput = nullptr;
    QLineEdit* m_editBdrExpr = nullptr;
    QCheckBox* m_chkBdrCorona;

    // 参数系统 UI 组件
    QWidget* m_pageParams = nullptr;
    QWidget* m_pageBoundary = nullptr;
    QTabWidget* m_paramTabs = nullptr;
    QLabel* m_paramStatus = nullptr;
    std::map<QDoubleSpinBox*, std::string> m_paramBindings;
    std::map<QDoubleSpinBox*, std::string> m_spinboxNames;  // spinbox→名称(持久化)
    std::vector<QDoubleSpinBox*> m_xfmSpins;  // 变换页面的局部spinbox
    std::map<QToolButton*, QDoubleSpinBox*> m_bindBtns;
    void installBindButtons();
    void showBindDialog(QDoubleSpinBox* spin, QToolButton* btn);
    void pushParamsToBindings();

    // 【修改】：后处理 UI 组件
    QWidget* m_pagePostProcess;
    QComboBox* m_comboPostField;      
    QRadioButton* m_radioSurface;     
    QRadioButton* m_radioVolume;      
    QListWidget* m_listPostGroups;    // 用于多选物理组
    QDoubleSpinBox* m_spinSliceX;
    QDoubleSpinBox* m_spinSliceY;
    QDoubleSpinBox* m_spinSliceZ;
    QCheckBox* m_chkSliceX;
    QCheckBox* m_chkSliceY;
    QCheckBox* m_chkSliceZ;
    QCheckBox* m_chkAutoRange;
    QDoubleSpinBox* m_spinRangeMin;
    QDoubleSpinBox* m_spinRangeMax;
    QPushButton* m_btnRenderPost;

    QString m_lastVtkFilePath;
    QString m_lastFieldName;
    QString m_projectFilePath;
    bool m_projectModified = false;

    void saveProject(const QString& path);
    void loadProject(const QString& path);
    bool maybeSave();
    void closeEvent(QCloseEvent* event) override;

    QProcess* m_solverProcess = nullptr; 
    ProjectConfig m_config;              

    Ui::MainWindow *ui;
    
    QTreeWidgetItem* m_rootGeomNode;
    QComboBox* m_geomSelectMode = nullptr;

    // 几何序列节点 (COMSOL 风格)
    struct GeomNode {
        QString type;
        QTreeWidgetItem* treeItem = nullptr;
        bool built = false;
        // 体素参数
        double w=100, h=100, d=100, e=0, px=0, py=0, pz=0, rotAng=0; int rotAxis=0;
        bool asSolid = true;  // 实体/表面
        // 布尔参数
        std::vector<int> inputIndices, toolIndices;
        bool keepInterior = false;
        bool keepInputs = true;
        bool keepTools = true;
        // 倒角: 存储选中的边 (nodeIndex, edgeLocalIdx)
        std::vector<std::pair<int,int>> edgeSelections;
        // 扫掠: 截面面的局部索引 (在builtShape的TopAbs_FACE遍历中的位置)
        int sweepProfileFaceIdx = 0;
        // 节点级参数绑定 (spinbox名称 → 参数名)
        std::map<std::string, std::string> nodeParamBindings;
        // 构建结果
        TopoDS_Shape builtShape;
        int builtFaceId = -1;
    };
    std::vector<GeomNode> m_geomNodes;
    void buildGeometrySequence(int upToIndex,
                               int forceIdx = -1,
                               const std::vector<int>* forceInputs = nullptr);
    void loadGeomNodeToUI(int idx);
    void saveGeomNodeFromUI(int idx);
    // 保存指定节点 UI 参数并链式构建 (btnBuild "构建/更新几何体" 与工具栏"全部构建"共用)
    void buildGeometryFromUI(int idx);

    // 几何设置面板控件 (用于 load/save)
    QDoubleSpinBox *m_geomW, *m_geomH, *m_geomD, *m_geomE, *m_geomPX, *m_geomPY, *m_geomPZ, *m_geomRotAng;
    QComboBox* m_geomRotAxis;
    QCheckBox *m_chkKeepInteriorU, *m_chkKeepInteriorS;
    QCheckBox *m_chkKeepInputsU, *m_chkKeepInputsS, *m_chkKeepToolsS;
    QComboBox* m_geomSolid;
    QListWidget *m_listBoolUnion, *m_listBoolSubIn, *m_listBoolSubTool;
    QStackedWidget* m_geomStackPages;
    QStackedWidget* m_geomPrimStack;
    QPushButton *m_btnAddUnion, *m_btnAddSubIn, *m_btnAddSubTool;
    QListWidget *m_listXfmT, *m_listXfmR, *m_listXfmM, *m_listXfmS, *m_listXfmO;
    QListWidget *m_listFilletEdges, *m_listChamfEdges, *m_listSwProf, *m_listSwPath;
    QListWidget *m_listArrLin, *m_listArrCir;
    QDoubleSpinBox *m_spinFilletR, *m_spinChamfD;
    QComboBox *m_cmbArrDir; QDoubleSpinBox *m_spArrSpc; QSpinBox *m_spArrCnt;
    QComboBox *m_cmbCirAx; QDoubleSpinBox *m_spCirCx, *m_spCirCy, *m_spCirCz, *m_spCirAng; QSpinBox *m_spCirCnt;
    int m_lastGeomIdx = -1;
    bool m_sweepSwitching = false;  // 扫掠面板按钮互斥切换标志

    QTreeWidgetItem* m_rootGroupNode;
    QTreeWidgetItem* m_rootMeshNode;
    QTreeWidgetItem* m_rootSolverNode;
    QTreeWidgetItem* m_rootPostNode;
    QTreeWidgetItem* m_rootParamNode;
    QTreeWidgetItem* m_rootBoundaryNode;
    
    QString m_currentStepFile;  

    struct BoundaryParams {
        bool apply = true;
        bool useFunction = false;
        double voltage = 0.0;
        std::string expression;
        bool is_corona = false;
    };

    std::map<QTreeWidgetItem*, QString> m_meshFiles;         
    std::map<QTreeWidgetItem*, bool>    m_isImportedMesh;    
    std::map<QTreeWidgetItem*, FaceMeshMap> m_meshDataMap;   

    std::map<QTreeWidgetItem*, std::set<int>> m_groupData;
    std::map<QTreeWidgetItem*, int> m_groupTags;  // Gmsh physical group tag per tree item
    std::map<int, int> m_elemToEntity;  // mesh element tag → Gmsh entity tag (imported meshes)
    std::map<int, int> m_faceToEntity;  // viewer face ID → Gmsh entity tag (imported meshes)
    std::map<int, std::set<int>> m_entityToFaces; // entity tag → viewer face IDs
    // 物理组元数据 (仅导入网格存储)
    struct GroupMeta { QString name; int tag; std::set<int> faceIds; };
    std::map<QTreeWidgetItem*, std::vector<GroupMeta>> m_importedGroupMeta;
    std::map<QTreeWidgetItem*, std::map<std::string, BoundaryParams>> m_importedBoundaryConfigs;
    QTreeWidgetItem* m_activeMeshNode = nullptr;
    std::map<QTreeWidgetItem*, std::map<QTreeWidgetItem*, MeshGroupParams>> m_meshConfigs;
    std::map<QTreeWidgetItem*, BoundaryParams> m_boundaryConfigs;

    bool m_isUpdatingUI = false;

    QProgressBar* m_progressBar;
    QPushButton* m_btnCancelMesh;
    QPushButton* m_btnCancelCompute;
    QPushButton* m_btnPartition = nullptr;
    QSpinBox* m_spinPartitions = nullptr;
    QTimer* m_meshLogTimer;
    std::atomic<bool> m_isMeshing;
    std::thread m_meshThread;
    FaceMeshMap m_lastMeshData;
    std::string m_lastMeshLog;
    QString m_meshSavePath;
    QString m_meshCachePath;
    QTreeWidgetItem* m_meshingNode;
    double m_simulatedProgress;

private slots:
    void onBoundaryGroupTargetChanged(int index);
    void onBoundaryParamsChanged();
    void onMeshFinished();
};