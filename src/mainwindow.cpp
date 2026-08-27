#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "cad_viewer.h"
#include "param_parser.h"
#include "geom_icons.h"
#include <QTabWidget>
#include <QTableWidget>
#include <QHeaderView>
#include <QInputDialog>
#include <QBoxLayout>
#include <QLabel>
#include <QDialogButtonBox>
#include "mesh_engine.h"
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkFeatureEdges.h>
#include <gmsh.h>
#include <QFileDialog>
#include <QMessageBox>
#include <QMenu>
#include <QScrollArea>
#include <QGroupBox>
#include <QFormLayout>
#include <QHeaderView>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QToolButton>
#include <QGridLayout>
#include <thread>
#include <QDir>
#include <QScrollBar>
#include <QSlider>
#include <QDebug>
#include <fstream>
#include <array>
#include "nlohmann/json.hpp" 
#include <QStandardItemModel>
#include <QSplitter>
#include <QLabel>
#include <QDateTime>
#include <QCloseEvent>
#include <QInputDialog>
#include <QFile>
#include <QFrame>
#include <QInputDialog>

// OpenCASCADE geometry creation
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakeWedge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_GTrsf.hxx>
#include <Standard_Failure.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <BRep_Builder.hxx>
#include <BRepLib.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <Geom_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomConvert_ApproxSurface.hxx>
#include <TopoDS_Compound.hxx>
#include <BOPAlgo_Builder.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopLoc_Location.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Shell.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <cmath>

using json = nlohmann::json;

// ====================================================================
// 流式布局: 一行放不下时自动换行 (参考 Qt 官方 Custom Layout Example)
// 用于: 几何操作面板 + 画布上方视图/选择功能栏, 画布变窄时自动换行成多行
// ====================================================================
class FlowLayout : public QLayout {
public:
    explicit FlowLayout(QWidget* parent = nullptr, int hSpacing = 6, int vSpacing = 6)
        : QLayout(parent), m_hSpace(hSpacing), m_vSpace(vSpacing) {}
    ~FlowLayout() { while (QLayoutItem* item = takeAt(0)) delete item; }

    void addItem(QLayoutItem* item) override { itemList.append(item); }
    int count() const override { return itemList.size(); }
    QLayoutItem* itemAt(int index) const override { return itemList.value(index); }
    QLayoutItem* takeAt(int index) override {
        if (index >= 0 && index < itemList.size()) return itemList.takeAt(index);
        return nullptr;
    }
    Qt::Orientations expandingDirections() const override { return {}; }
    bool hasHeightForWidth() const override { return true; }
    int heightForWidth(int width) const override { return doLayout(QRect(0, 0, width, 0), true); }
    QSize sizeHint() const override { return minimumSize(); }
    QSize minimumSize() const override {
        QSize size;
        for (const QLayoutItem* item : itemList) size = size.expandedTo(item->minimumSize());
        const QMargins m = contentsMargins();
        size += QSize(m.left() + m.right(), m.top() + m.bottom());
        return size;
    }
    void setGeometry(const QRect& rect) override {
        QLayout::setGeometry(rect);
        doLayout(rect, false);
    }

private:
    int doLayout(const QRect& rect, bool testOnly) const {
        const QMargins m = contentsMargins();
        int x = rect.x() + m.left();
        int y = rect.y() + m.top();
        int right = rect.right() - m.right();
        int lineHeight = 0;
        for (QLayoutItem* item : itemList) {
            int w = item->sizeHint().width();
            int nextX = x + w + m_hSpace;
            if (nextX - m_hSpace > right && lineHeight > 0) {   // 放不下 → 换行
                x = rect.x() + m.left();
                y += lineHeight + m_vSpace;
                nextX = x + w + m_hSpace;
                lineHeight = 0;
            }
            if (!testOnly) item->setGeometry(QRect(QPoint(x, y), item->sizeHint()));
            x = nextX;
            lineHeight = qMax(lineHeight, item->sizeHint().height());
        }
        return y + lineHeight - rect.y() + m.bottom();
    }

    QList<QLayoutItem*> itemList;
    int m_hSpace, m_vSpace;
};

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    this->setStyleSheet(R"(
        QMainWindow { background-color: #F3F5F8; }
        QDockWidget { border: 1px solid #D1D5DB; font-size: 16px; }
        QDockWidget::title { background-color: #E5E7EB; padding: 6px; font-weight: bold; color: #374151; }
        QGroupBox { border: 1px solid #CED4DA; border-radius: 6px; margin-top: 14px; background-color: #FAFAFA; }
        QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; left: 10px; padding: 0 5px; color: #2B5797; font-weight: bold; font-size: 13px; }
        QTreeWidget, QTableWidget, QListWidget { border: 1px solid #CED4DA; border-radius: 4px; background-color: #FFFFFF; alternate-background-color: #F8F9FA; font-size: 13px; }
        QHeaderView::section { background-color: #E9ECEF; padding: 4px; border: 1px solid #DEE2E6; font-weight: bold; font-size: 13px; }
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox { border: 1px solid #CED4DA; border-radius: 4px; padding: 4px 8px; background-color: #FFFFFF; selection-background-color: #007BFF; font-size: 13px; }
        QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border: 1px solid #80BDFF; }
        QPushButton { background-color: #F8F9FA; border: 1px solid #CED4DA; border-radius: 4px; padding: 6px 12px; color: #212529; font-weight: 500; font-size: 13px; }
        QPushButton:hover { background-color: #E2E6EA; border-color: #DAE0E5; }
        QPushButton:pressed { background-color: #DAE0E5; border-color: #D3D9DF; }
        QPushButton:disabled { background-color: #E9ECEF; color: #ADB5BD; }
        QPushButton#btnGenerateMesh, QPushButton#btnCalculate, QPushButton#btnRenderPost { background-color: #2B5797; color: white; border: none; }
        QPushButton#btnGenerateMesh:hover, QPushButton#btnCalculate:hover, QPushButton#btnRenderPost:hover { background-color: #1E3F70; }
        QPushButton#btnGenerateMesh:pressed, QPushButton#btnCalculate:pressed, QPushButton#btnRenderPost:pressed { background-color: #152C4D; }
        QLabel, QCheckBox, QToolButton, QTabBar::tab { font-size: 13px; }
        QMenuBar { font-size: 13px; }
        QMenuBar::item { font-size: 13px; }
        QMenu { font-size: 13px; }
        QMenu::item { font-size: 13px; }
    )");

    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) maxThreads = 4; 
    ui->spinMeshThreads->setRange(1, 128);
    ui->spinMeshThreads->setValue(maxThreads); 

    ui->treeGroups->setHeaderHidden(true);

    // 几何创建模块 (最上方)
    m_rootGeomNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootGeomNode->setText(0, "创建几何");
    m_rootGeomNode->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    m_rootGeomNode->setExpanded(true);

    m_rootGroupNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootGroupNode->setText(0, "物理组 (边界选择)");
    m_rootGroupNode->setExpanded(true);

    m_rootMeshNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootMeshNode->setText(0, "网格生成方案");
    m_rootMeshNode->setExpanded(true);

    m_rootSolverNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootSolverNode->setText(0, "求解器与参数设置");
    m_rootSolverNode->setExpanded(true);

    m_rootPostNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootPostNode->setText(0, "后处理与可视化");
    m_rootPostNode->setExpanded(true);

    m_rootParamNode = new QTreeWidgetItem(ui->treeGroups);
    m_rootParamNode->setText(0, "全局参数");
    m_rootParamNode->setExpanded(false);

    // 树节点重排: 全局参数置顶, 边界条件插入网格与求解器之间
    int paramIdx = ui->treeGroups->indexOfTopLevelItem(m_rootParamNode);
    ui->treeGroups->insertTopLevelItem(0, ui->treeGroups->takeTopLevelItem(paramIdx));
    int solverIdx = ui->treeGroups->indexOfTopLevelItem(m_rootSolverNode);
    m_rootBoundaryNode = new QTreeWidgetItem(); m_rootBoundaryNode->setText(0, "边界条件");
    ui->treeGroups->insertTopLevelItem(solverIdx, m_rootBoundaryNode);

    ui->stackedWidgetSettings->setCurrentIndex(3);

    ui->dockModelBuilder->setFeatures(QDockWidget::DockWidgetMovable);
    ui->dockSettings->setFeatures(QDockWidget::DockWidgetMovable);
    this->splitDockWidget(ui->dockModelBuilder, ui->dockSettings, Qt::Horizontal);

    QDockWidget* dockConsole = new QDockWidget("消息与系统日志", this);
    dockConsole->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetClosable);
    dockConsole->setAllowedAreas(Qt::BottomDockWidgetArea | Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);

    // 停靠窗标题字号: 统一放大为 16px (QSS 的 ::title 不支持 font, 需直接设置 QFont; 子控件有 13px 规则覆盖不受影响)
    QFont dockTitleFont = this->font();
    dockTitleFont.setPixelSize(16);
    ui->dockModelBuilder->setFont(dockTitleFont);
    ui->dockSettings->setFont(dockTitleFont);
    dockConsole->setFont(dockTitleFont);

    m_textConsole = new QTextEdit();
    m_textConsole->setReadOnly(true);
    m_textConsole->setStyleSheet(R"(
        QTextEdit { background-color: #1E1E1E; color: #CCCCCC; font-family: 'Consolas', 'Courier New', monospace; font-size: 13px; border: none; padding: 8px; line-height: 1.5; }
    )");
    dockConsole->setWidget(m_textConsole);
    this->setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
    this->setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);
    this->addDockWidget(Qt::BottomDockWidgetArea, dockConsole);

    // 设置停靠窗内容包一层滚动区: 使宽度可自由伸缩
    // (否则 QStackedWidget 的 minimumSizeHint 取所有页面最大值, 会把 dock 最小宽度撑死)
    {
        QWidget* oldContent = ui->dockSettings->widget();   // dockWidgetContents_2 (含 stackedWidgetSettings)
        auto* settingsScroll = new QScrollArea(ui->dockSettings);
        settingsScroll->setWidgetResizable(true);
        settingsScroll->setFrameShape(QFrame::NoFrame);
        settingsScroll->setWidget(oldContent);
        ui->dockSettings->setWidget(settingsScroll);
    }

    QStatusBar* statBar = this->statusBar();
    m_progressBar = new QProgressBar();
    m_progressBar->setRange(0, 100);
    m_progressBar->setValue(0);
    m_progressBar->setMinimumWidth(200);
    m_progressBar->setMaximumWidth(250);
    m_progressBar->setMaximumHeight(16);
    m_progressBar->setTextVisible(true);
    m_progressBar->setAlignment(Qt::AlignCenter);
    m_progressBar->setStyleSheet(
        "QProgressBar { border: 1px solid #CED4DA; border-radius: 2px; background-color: #E9ECEF; text-align: center; color: #212529; font-size: 10px; }"
        "QProgressBar::chunk { background-color: #28A745; border-radius: 1px; }"
    );
    m_progressBar->setVisible(false);

    m_btnCancelMesh = new QPushButton("终止操作");
    m_btnCancelMesh->setMaximumHeight(20);
    m_btnCancelMesh->setStyleSheet(
        "QPushButton { background-color: #DC3545; color: white; border: none; border-radius: 2px; padding: 0px 10px; font-size: 11px; font-weight: bold; }"
        "QPushButton:hover { background-color: #C82333; }"
        "QPushButton:pressed { background-color: #BD2130; }"
        "QPushButton:disabled { background-color: #E9ECEF; color: #6C757D; }"
    );
    m_btnCancelMesh->setVisible(false);

    m_btnCancelCompute = new QPushButton("终止计算");
    m_btnCancelCompute->setMaximumHeight(20);
    m_btnCancelCompute->setStyleSheet(
        "QPushButton { background-color: #DC3545; color: white; border: none; border-radius: 2px; padding: 0px 10px; font-size: 11px; font-weight: bold; }"
        "QPushButton:hover { background-color: #C82333; }"
        "QPushButton:pressed { background-color: #BD2130; }"
        "QPushButton:disabled { background-color: #E9ECEF; color: #6C757D; }"
    );
    m_btnCancelCompute->setVisible(false);
    connect(m_btnCancelCompute, &QPushButton::clicked, this, [=](){
        if (m_solverProcess && m_solverProcess->state() == QProcess::Running) {
            logMessage("正在终止计算进程...", "warning");
            m_solverProcess->kill();
            m_btnCancelCompute->setEnabled(false);
            m_btnCancelCompute->setText("正在终止...");
        }
    });

    statBar->addPermanentWidget(m_progressBar);
    statBar->addPermanentWidget(m_btnCancelMesh);
    statBar->addPermanentWidget(m_btnCancelCompute);

    m_meshLogTimer = new QTimer(this);
    m_isMeshing = false;
    connect(m_meshLogTimer, &QTimer::timeout, this, [=](){
        if (m_isMeshing) {
            m_simulatedProgress += (99.0 - m_simulatedProgress) * 0.03; 
            m_progressBar->setValue(static_cast<int>(m_simulatedProgress));
        }
    });

    connect(m_btnCancelMesh, &QPushButton::clicked, this, [=](){
        if (m_isMeshing) {
            logMessage("正在请求终止网格生成，请稍候...", "warning");
            try { gmsh::option::setNumber("General.Cancel", 1); } catch (...) {}
            m_btnCancelMesh->setEnabled(false);
            m_btnCancelMesh->setText("正在终止...");
        }
    });

    ui->mainToolBar->setVisible(false);  // 隐藏工具栏, 导入功能归入 File 菜单
    setupGeometryToolbar();

    // ========== 几何参数页: 体素参数 + 布尔参数 (QStackedWidget 切换) ==========
    m_pageGeomParam = new QWidget();
    auto* geomOuterLayout = new QVBoxLayout(m_pageGeomParam);

    // 公共构建按钮
    QPushButton* btnBuild = new QPushButton("构建/更新几何体");
    btnBuild->setMinimumHeight(35);
    btnBuild->setStyleSheet("QPushButton { background-color: #2B5797; color: white; border: none; }");
    geomOuterLayout->addWidget(btnBuild);

    m_geomStackPages = new QStackedWidget();
    geomOuterLayout->addWidget(m_geomStackPages);

    // -- Page 0: 体素参数 (动态标签) --
    m_geomW = new QDoubleSpinBox(); m_geomW->setRange(1e-9,1e6); m_geomW->setDecimals(6); m_geomW->setValue(100); m_geomW->setSuffix(" m");
    m_geomH = new QDoubleSpinBox(); m_geomH->setRange(1e-9,1e6); m_geomH->setDecimals(6); m_geomH->setValue(100); m_geomH->setSuffix(" m");
    m_geomD = new QDoubleSpinBox(); m_geomD->setRange(1e-9,1e6); m_geomD->setDecimals(6); m_geomD->setValue(100); m_geomD->setSuffix(" m");
    m_geomE = new QDoubleSpinBox(); m_geomE->setRange(-1e6,1e6); m_geomE->setValue(3); m_geomE->setDecimals(0);
    QLabel* lblDim1 = new QLabel("长 (W):");
    QLabel* lblDim2 = new QLabel("宽 (H):");
    QLabel* lblDim3 = new QLabel("高 (D):");
    QLabel* lblDim4 = new QLabel("");
    m_geomPX = new QDoubleSpinBox(); m_geomPX->setRange(-1e6,1e6); m_geomPX->setDecimals(6); m_geomPX->setValue(0); m_geomPX->setSuffix(" m");
    m_geomPY = new QDoubleSpinBox(); m_geomPY->setRange(-1e6,1e6); m_geomPY->setDecimals(6); m_geomPY->setValue(0); m_geomPY->setSuffix(" m");
    m_geomPZ = new QDoubleSpinBox(); m_geomPZ->setRange(-1e6,1e6); m_geomPZ->setDecimals(6); m_geomPZ->setValue(0); m_geomPZ->setSuffix(" m");
    m_geomRotAng = new QDoubleSpinBox(); m_geomRotAng->setRange(-360,360); m_geomRotAng->setDecimals(4); m_geomRotAng->setValue(0); m_geomRotAng->setSuffix(" deg");
    m_geomRotAxis = new QComboBox(); m_geomRotAxis->addItems({"X","Y","Z"});

    auto* primPage = new QWidget();
    auto* primLayout = new QFormLayout(primPage);
    primLayout->addRow(lblDim1, m_geomW);
    primLayout->addRow(lblDim2, m_geomH);
    primLayout->addRow(lblDim3, m_geomD);
    primLayout->addRow(lblDim4, m_geomE);
    primLayout->addRow("位置 X:", m_geomPX); primLayout->addRow("位置 Y:", m_geomPY); primLayout->addRow("位置 Z:", m_geomPZ);
    primLayout->addRow("旋转轴:", m_geomRotAxis); primLayout->addRow("旋转角度:", m_geomRotAng);
    m_geomStackPages->addWidget(primPage);

    // 动态标签更新 (loadGeomNodeToUI 中调用)
    auto updatePrimLabels = [=](const QString& t) {
        if (t=="立方体")      { lblDim1->setText("长 (W):"); lblDim2->setText("宽 (H):"); lblDim3->setText("高 (D):");
                                m_geomW->show(); m_geomH->show(); m_geomD->show(); }
        else if (t=="球体")   { lblDim1->setText("半径 (R):"); lblDim2->setText(""); lblDim3->setText("");
                                m_geomW->show(); m_geomH->hide(); m_geomD->hide(); }
        else if (t=="圆柱")   { lblDim1->setText("半径 (R):"); lblDim2->setText("高度 (H):"); lblDim3->setText("");
                                m_geomW->show(); m_geomH->show(); m_geomD->hide(); }
        else if (t=="圆锥")   { lblDim1->setText("底半径 (R1):"); lblDim2->setText("顶半径 (R2):"); lblDim3->setText("高度 (H):");
                                m_geomW->show(); m_geomH->show(); m_geomD->show(); }
        else if (t=="圆环")   { lblDim1->setText("主半径 (R1):"); lblDim2->setText("副半径 (R2):"); lblDim3->setText("");
                                m_geomW->show(); m_geomH->show(); m_geomD->hide(); }
        else                  { lblDim1->setText("长 (W):"); lblDim2->setText("宽 (H):"); lblDim3->setText("高 (D):");
                                m_geomW->show(); m_geomH->show(); m_geomD->show(); }
    };
    m_geomPrimStack = nullptr; // 不再使用栈, 改用动态标签

    // -- Page 1: 布尔 - 并集/交集 --
    auto* boolUnionPage = new QWidget();
    auto* boolULayout = new QVBoxLayout(boolUnionPage);
    m_btnAddUnion = new QPushButton("+ 加入对象");
    boolULayout->addWidget(m_btnAddUnion);
    m_listBoolUnion = new QListWidget(); m_listBoolUnion->setSelectionMode(QAbstractItemView::MultiSelection);
    boolULayout->addWidget(m_listBoolUnion);
    m_geomStackPages->addWidget(boolUnionPage);

    // -- Page 2: 布尔 - 差集 --
    auto* boolSubPage = new QWidget();
    auto* boolSLayout = new QVBoxLayout(boolSubPage);
    boolSLayout->addWidget(new QLabel("要添加的对象:"));
    m_btnAddSubIn = new QPushButton("+ 加入对象"); boolSLayout->addWidget(m_btnAddSubIn);
    m_listBoolSubIn = new QListWidget(); m_listBoolSubIn->setSelectionMode(QAbstractItemView::MultiSelection);
    boolSLayout->addWidget(m_listBoolSubIn);
    boolSLayout->addWidget(new QLabel("要减去的对象:"));
    m_btnAddSubTool = new QPushButton("+ 加入对象"); boolSLayout->addWidget(m_btnAddSubTool);
    m_listBoolSubTool = new QListWidget(); m_listBoolSubTool->setSelectionMode(QAbstractItemView::MultiSelection);
    boolSLayout->addWidget(m_listBoolSubTool);
    m_geomStackPages->addWidget(boolSubPage);

    ui->stackedWidgetSettings->addWidget(m_pageGeomParam);

    // --- 保留内部边界复选框 (每个布尔页面独立) ---
    m_chkKeepInteriorU = new QCheckBox("保留内部边界");
    m_chkKeepInteriorS = new QCheckBox("保留内部边界");
    // 加入对象开关按钮 (需先定义, 供布尔和变换页面共用)
    // mode: 0=对象选择, 1=边选择(添加边), 2=边选择(添加节点), 3=面选择(添加节点-扫掠截面)
    auto setupToggleBtn = [this](QPushButton* btn, QListWidget* list, int mode=0) {
        bool isEdgeMode = (mode == 1 || mode == 2);
        bool isFaceMode = (mode == 3);
        btn->setCheckable(true);
        QString label, onLabel;
        if (mode == 0)      { label = "关闭: 加入对象";    onLabel = "开启: 点击对象添加"; }
        else if (mode == 1) { label = "关闭: 加入边";      onLabel = "开启: 点击边添加"; }
        else if (mode == 2) { label = "关闭: 加入路径边";  onLabel = "开启: 点击边添加路径"; }
        else                { label = "关闭: 加入面";      onLabel = "开启: 点击面添加截面"; }
        btn->setText(label);
        btn->setStyleSheet("QPushButton { background-color: #E0E0E0; color: #333; }");
        QMetaObject::Connection* conn = new QMetaObject::Connection();
        connect(btn, &QPushButton::toggled, [=](bool on) mutable {
            if (on) {
                btn->setText(onLabel);
                btn->setStyleSheet("QPushButton { background-color: #DC3545; color: white; }");
                int curIdx = m_lastGeomIdx;
                if (curIdx >= 0 && curIdx < (int)m_geomNodes.size() && m_geomNodes[curIdx].built) {
                    buildGeometrySequence(curIdx - 1);
                    m_geomNodes[curIdx].built = false;
                    m_geomNodes[curIdx].treeItem->setForeground(0, QBrush(QColor(128,128,128)));
                }
                if (isEdgeMode) {
                    ui->cadViewerWidget->setEdgeSelectMode(true);
                    m_geomSelectMode->setCurrentIndex(2);
                    ui->cadViewerWidget->clearAllEdgeHighlights();
                    if (mode == 1) {
                        for (int r = 0; r < list->count(); r++) {
                            int ni = list->item(r)->data(Qt::UserRole).toInt();
                            int ei = list->item(r)->data(Qt::UserRole+1).toInt();
                            if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].built)
                                ui->cadViewerWidget->highlightEdgeInShape(m_geomNodes[ni].builtShape, ei, m_geomNodes[ni].builtFaceId);
                        }
                    }
                } else if (isFaceMode) {
                    // 面选择模式: 恢复列表中面的高亮
                    ui->cadViewerWidget->setObjectSelectMode(false);
                    ui->cadViewerWidget->setEdgeSelectMode(false);
                    m_geomSelectMode->setCurrentIndex(0);
                    std::set<int> selFaces;
                    for (int r = 0; r < list->count(); r++) {
                        int ni = list->item(r)->data(Qt::UserRole).toInt();
                        int localIdx = list->item(r)->data(Qt::UserRole+1).toInt();
                        if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].built)
                            selFaces.insert(m_geomNodes[ni].builtFaceId + localIdx);
                    }
                    ui->cadViewerWidget->setSelectedFaces(selFaces);
                } else {
                    // 对象选择模式: 恢复列表中对象所有面的高亮
                    ui->cadViewerWidget->setObjectSelectMode(true);
                    m_geomSelectMode->setCurrentIndex(1);
                    std::set<int> selFaces;
                    const auto& fom = ui->cadViewerWidget->faceObjectMap();
                    for (int r = 0; r < list->count(); r++) {
                        int ni = list->item(r)->data(Qt::UserRole).toInt();
                        if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].built) {
                            int objId = m_geomNodes[ni].builtFaceId;
                            for (auto& [fid, oid] : fom)
                                if (oid == objId) selFaces.insert(fid);
                        }
                    }
                    ui->cadViewerWidget->setSelectedFaces(selFaces);
                }
                // 信号: 根据mode处理不同的列表绑定
                *conn = connect(ui->cadViewerWidget, &CADViewer::faceSelectionChanged, [=](const std::set<int>& fids) mutable {
                    if (fids.empty()) { list->clear(); return; }
                    if (mode == 1) {
                        // 边模式(添加边): 双向同步 (ni, ei) 对
                        int pickedObj = ui->cadViewerWidget->lastPickedObjId();
                        int pickedEdge = ui->cadViewerWidget->lastPickedEdgeIdx();
                        if (pickedObj < 0 || pickedEdge < 0) {
                            // 框选边: 用viewer中所有选中边重建列表
                            if (!fids.empty()) {
                                auto edges = ui->cadViewerWidget->getSelectedEdges();
                                list->clear();
                                for (auto& [objId, edgeIdx] : edges) {
                                    for (int ni = 0; ni < (int)m_geomNodes.size(); ni++) {
                                        if (m_geomNodes[ni].built && m_geomNodes[ni].builtFaceId == objId) {
                                            auto* li = new QListWidgetItem(QString("%1 - 边%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(edgeIdx));
                                            li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, edgeIdx);
                                            list->addItem(li);
                                            break;
                                        }
                                    }
                                }
                            }
                            return;
                        }
                        for (int ni = 0; ni < (int)m_geomNodes.size(); ni++) {
                            if (!m_geomNodes[ni].built) continue;
                            if (m_geomNodes[ni].builtFaceId == pickedObj) {
                                bool inList = false; int lr = -1;
                                for (int r=0; r<list->count(); r++)
                                    if (list->item(r)->data(Qt::UserRole).toInt()==ni && list->item(r)->data(Qt::UserRole+1).toInt()==pickedEdge) { inList=true; lr=r; break; }
                                bool edgeSel = ui->cadViewerWidget->isEdgeSelected(pickedObj, pickedEdge);
                                if (edgeSel && !inList) {
                                    auto* li = new QListWidgetItem(QString("%1 - 边%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(pickedEdge));
                                    li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, pickedEdge);
                                    list->addItem(li);
                                } else if (!edgeSel && inList) {
                                    delete list->takeItem(lr);
                                }
                                break;
                            }
                        }
                    } else if (mode == 2) {
                        // 边模式(添加节点): 点击边→添加/移除其所属节点
                        int pickedObj = ui->cadViewerWidget->lastPickedObjId();
                        if (pickedObj < 0) return;
                        for (int ni = 0; ni < (int)m_geomNodes.size(); ni++) {
                            if (!m_geomNodes[ni].built) continue;
                            if (m_geomNodes[ni].builtFaceId == pickedObj) {
                                bool inList = false; int lr = -1;
                                for (int r=0; r<list->count(); r++)
                                    if (list->item(r)->data(Qt::UserRole).toInt()==ni) { inList=true; lr=r; break; }
                                if (!inList) {
                                    auto* li = new QListWidgetItem(m_geomNodes[ni].treeItem->text(0));
                                    li->setData(Qt::UserRole, ni);
                                    list->addItem(li);
                                } else {
                                    delete list->takeItem(lr);
                                }
                                break;
                            }
                        }
                    } else if (mode == 3) {
                        // mode 3: 面选择→存储具体面索引 (ni, localFaceIdx)
                        // 清空列表, 用当前选中的面重建
                        list->clear();
                        const auto& fom = ui->cadViewerWidget->faceObjectMap();
                        for (int fid : fids) {
                            auto it = fom.find(fid);
                            if (it == fom.end()) continue;
                            int objId = it->second;
                            for (int ni = 0; ni < (int)m_geomNodes.size(); ni++) {
                                if (!m_geomNodes[ni].built) continue;
                                if (m_geomNodes[ni].builtFaceId == objId) {
                                    int localIdx = fid - objId;
                                    auto* li = new QListWidgetItem(QString("%1 - 面%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(localIdx));
                                    li->setData(Qt::UserRole, ni);
                                    li->setData(Qt::UserRole+1, localIdx);
                                    list->addItem(li);
                                    break;
                                }
                            }
                        }
                    } else {
                        // mode 0: 对象模式→双向同步节点
                        for (int ni = 0; ni < (int)m_geomNodes.size(); ni++) {
                            if (!m_geomNodes[ni].built) continue;
                            const auto& fom = ui->cadViewerWidget->faceObjectMap();
                            auto it = fom.find(m_geomNodes[ni].builtFaceId);
                            int objId = it != fom.end() ? it->second : -1;
                            if (objId < 0) continue;
                            bool sel = false;
                            for (auto& [fid,oid] : fom) if (oid==objId && fids.count(fid)) { sel=true; break; }
                            bool inList = false; int lr = -1;
                            for (int r=0; r<list->count(); r++) if (list->item(r)->data(Qt::UserRole).toInt()==ni) { inList=true; lr=r; break; }
                            if (sel && !inList) { auto* li = new QListWidgetItem(m_geomNodes[ni].treeItem->text(0)); li->setData(Qt::UserRole,ni); list->addItem(li); }
                            else if (!sel && inList) delete list->takeItem(lr);
                        }
                    }
                });
            } else {
                btn->setText(label); btn->setStyleSheet("QPushButton { background-color: #E0E0E0; color: #333; }");
                ui->cadViewerWidget->setSelectedFaces(std::set<int>());
                // 互斥切换时不改变viewer模式(由新按钮的ON处理器设置)
                if (!m_sweepSwitching) {
                    ui->cadViewerWidget->setObjectSelectMode(false); ui->cadViewerWidget->setEdgeSelectMode(false);
                    m_geomSelectMode->setCurrentIndex(0);
                }
                if (isEdgeMode) ui->cadViewerWidget->setEdgeHighlightsVisible(false);
                QObject::disconnect(*conn);
            }
        });
    };

    // Page 3-7: 变换 (平移, 旋转, 镜像, 缩放, 偏移)
    auto makeXfmPage = [&](const QString& type) -> QPair<QPushButton*, QListWidget*> {
        auto* pg = new QWidget(); auto* lay = new QVBoxLayout(pg);
        auto* btn = new QPushButton("+ 加入对象"); btn->setCheckable(true);
        lay->addWidget(btn);
        auto* list = new QListWidget(); list->setSelectionMode(QAbstractItemView::MultiSelection);
        lay->addWidget(list);
        QCheckBox* chk = new QCheckBox("保留输入对象"); chk->setChecked(true);
        lay->addWidget(chk);
        // 特定参数
        if (type=="平移") {
            auto* fl = new QFormLayout(); lay->addLayout(fl);
            auto* dx = new QDoubleSpinBox(); dx->setRange(-1e6,1e6); dx->setValue(0); dx->setDecimals(6); dx->setSuffix(" m");
            auto* dy = new QDoubleSpinBox(); dy->setRange(-1e6,1e6); dy->setValue(0); dy->setDecimals(6); dy->setSuffix(" m");
            auto* dz = new QDoubleSpinBox(); dz->setRange(-1e6,1e6); dz->setValue(0); dz->setDecimals(6); dz->setSuffix(" m");
            m_xfmSpins.push_back(dx); m_xfmSpins.push_back(dy); m_xfmSpins.push_back(dz);
            fl->addRow("X 位移:", dx); fl->addRow("Y 位移:", dy); fl->addRow("Z 位移:", dz);
        } else if (type=="旋转") {
            auto* fl = new QFormLayout(); lay->addLayout(fl);
            auto* ax = new QComboBox(); ax->addItems({"X","Y","Z"});
            auto* ang = new QDoubleSpinBox(); ang->setRange(-360,360); ang->setValue(0); ang->setDecimals(4); ang->setSuffix(" deg");
            auto* cx = new QDoubleSpinBox(); cx->setRange(-1e6,1e6); cx->setValue(0); cx->setDecimals(6); cx->setSuffix(" m");
            auto* cy = new QDoubleSpinBox(); cy->setRange(-1e6,1e6); cy->setValue(0); cy->setDecimals(6); cy->setSuffix(" m");
            auto* cz = new QDoubleSpinBox(); cz->setRange(-1e6,1e6); cz->setValue(0); cz->setDecimals(6); cz->setSuffix(" m");
            m_xfmSpins.push_back(ang); m_xfmSpins.push_back(cx); m_xfmSpins.push_back(cy); m_xfmSpins.push_back(cz);
            fl->addRow("旋转轴:", ax); fl->addRow("角度:", ang);
            fl->addRow("中心 X:", cx); fl->addRow("中心 Y:", cy); fl->addRow("中心 Z:", cz);
        } else if (type=="镜像") {
            auto* fl = new QFormLayout(); lay->addLayout(fl);
            auto* pl = new QComboBox(); pl->addItems({"XY平面","YZ平面","ZX平面"});
            auto* pp = new QDoubleSpinBox(); pp->setRange(-1e6,1e6); pp->setValue(0); pp->setDecimals(6); pp->setSuffix(" m");
            m_xfmSpins.push_back(pp);
            fl->addRow("镜像平面:", pl);
            fl->addRow("平面位置:", pp);
        } else if (type=="缩放") {
            auto* fl = new QFormLayout(); lay->addLayout(fl);
            auto* fac = new QDoubleSpinBox(); fac->setRange(1e-9,100); fac->setDecimals(6); fac->setValue(1.0); fac->setSingleStep(0.1);
            auto* cx = new QDoubleSpinBox(); cx->setRange(-1e6,1e6); cx->setValue(0); cx->setDecimals(6); cx->setSuffix(" m");
            auto* cy = new QDoubleSpinBox(); cy->setRange(-1e6,1e6); cy->setValue(0); cy->setDecimals(6); cy->setSuffix(" m");
            auto* cz = new QDoubleSpinBox(); cz->setRange(-1e6,1e6); cz->setValue(0); cz->setDecimals(6); cz->setSuffix(" m");
            m_xfmSpins.push_back(fac); m_xfmSpins.push_back(cx); m_xfmSpins.push_back(cy); m_xfmSpins.push_back(cz);
            fl->addRow("缩放因子:", fac); fl->addRow("中心 X:", cx); fl->addRow("中心 Y:", cy); fl->addRow("中心 Z:", cz);
        } else if (type=="偏移") {
            auto* fl = new QFormLayout(); lay->addLayout(fl);
            auto* dist = new QDoubleSpinBox(); dist->setRange(-1e6,1e6); dist->setValue(1); dist->setDecimals(6); dist->setSuffix(" m");
            m_xfmSpins.push_back(dist);
            fl->addRow("偏移距离:", dist);
        }
        m_geomStackPages->addWidget(pg);
        return {btn, list};
    };
    auto pT = makeXfmPage("平移"), pR = makeXfmPage("旋转"), pM = makeXfmPage("镜像"),
         pS = makeXfmPage("缩放"), pO = makeXfmPage("偏移");
    m_listXfmT=pT.second; m_listXfmR=pR.second; m_listXfmM=pM.second;
    m_listXfmS=pS.second; m_listXfmO=pO.second;
    setupToggleBtn(pT.first, pT.second); setupToggleBtn(pR.first, pR.second);
    setupToggleBtn(pM.first, pM.second); setupToggleBtn(pS.first, pS.second);
    setupToggleBtn(pO.first, pO.second);

    // Page 8: 倒圆角 — 边列表 + 半径
    auto makeFilletPage = [&]() {
        auto* pg = new QWidget(); auto* lay = new QVBoxLayout(pg);
        auto* btn = new QPushButton("+ 加入边"); btn->setCheckable(true);
        lay->addWidget(btn); m_listFilletEdges = new QListWidget(); lay->addWidget(m_listFilletEdges);
        auto* fl = new QFormLayout(); lay->addLayout(fl);
        m_spinFilletR = new QDoubleSpinBox(); m_spinFilletR->setRange(1e-9,1e6); m_spinFilletR->setDecimals(6); m_spinFilletR->setValue(5); m_spinFilletR->setSuffix(" m");
        fl->addRow("圆角半径:", m_spinFilletR);
        m_geomStackPages->addWidget(pg);
        setupToggleBtn(btn, m_listFilletEdges, 1);
    };
    makeFilletPage();

    // Page 9: 倒斜角 — 边列表 + 距离
    auto* chamfPage = new QWidget(); auto* chamfLay = new QVBoxLayout(chamfPage);
    auto* btnChamf = new QPushButton("+ 加入边"); btnChamf->setCheckable(true);
    chamfLay->addWidget(btnChamf); m_listChamfEdges = new QListWidget(); chamfLay->addWidget(m_listChamfEdges);
    auto* chamfFl = new QFormLayout(); chamfLay->addLayout(chamfFl);
    m_spinChamfD = new QDoubleSpinBox(); m_spinChamfD->setRange(1e-9,1e6); m_spinChamfD->setDecimals(6); m_spinChamfD->setValue(5); m_spinChamfD->setSuffix(" m");
    chamfFl->addRow("斜角距离:", m_spinChamfD);
    m_geomStackPages->addWidget(chamfPage);
    setupToggleBtn(btnChamf, m_listChamfEdges, 1);

    // Page 10: 扫掠 — 截面列表(面选择) + 路径列表(边选择), 互斥
    auto* sweepPage = new QWidget(); auto* sweepLay = new QVBoxLayout(sweepPage);
    sweepLay->addWidget(new QLabel("截面轮廓 (面):"));
    auto* btnSwProf = new QPushButton("+ 加入截面"); btnSwProf->setCheckable(true);
    sweepLay->addWidget(btnSwProf); m_listSwProf = new QListWidget(); sweepLay->addWidget(m_listSwProf);
    sweepLay->addWidget(new QLabel("扫掠路径 (边):"));
    auto* btnSwPath = new QPushButton("+ 加入路径"); btnSwPath->setCheckable(true);
    sweepLay->addWidget(btnSwPath); m_listSwPath = new QListWidget(); sweepLay->addWidget(m_listSwPath);
    m_geomStackPages->addWidget(sweepPage);
    // 截面: 模式3=面选择添加节点, 路径: 模式1=边选择(选择具体边作为路径)
    setupToggleBtn(btnSwProf, m_listSwProf, 3);
    setupToggleBtn(btnSwPath, m_listSwPath, 1);
    // 互斥: 开启一个自动关闭另一个
    connect(btnSwProf, &QPushButton::toggled, [=](bool on) { if (on) { m_sweepSwitching = true; btnSwPath->setChecked(false); m_sweepSwitching = false; } });
    connect(btnSwPath, &QPushButton::toggled, [=](bool on) { if (on) { m_sweepSwitching = true; btnSwProf->setChecked(false); m_sweepSwitching = false; } });

    // Page 11: 线性阵列 — 对象列表 + 方向/间距/数量
    auto* arrLinPage = new QWidget(); auto* arrLinLay = new QVBoxLayout(arrLinPage);
    auto* btnArrLin = new QPushButton("+ 加入对象"); btnArrLin->setCheckable(true);
    arrLinLay->addWidget(btnArrLin); m_listArrLin = new QListWidget(); arrLinLay->addWidget(m_listArrLin);
    auto* arrLinFl = new QFormLayout(); arrLinLay->addLayout(arrLinFl);
    m_cmbArrDir = new QComboBox(); m_cmbArrDir->addItems({"X","Y","Z"});
    m_spArrSpc = new QDoubleSpinBox(); m_spArrSpc->setRange(1e-9,1e6); m_spArrSpc->setDecimals(6); m_spArrSpc->setValue(100); m_spArrSpc->setSuffix(" m");
    m_spArrCnt = new QSpinBox(); m_spArrCnt->setRange(2,100); m_spArrCnt->setValue(3);
    arrLinFl->addRow("方向:", m_cmbArrDir); arrLinFl->addRow("间距:", m_spArrSpc); arrLinFl->addRow("数量:", m_spArrCnt);
    auto* chkArrLinKeep = new QCheckBox("保留输入对象"); chkArrLinKeep->setChecked(true);
    arrLinLay->addWidget(chkArrLinKeep);
    m_geomStackPages->addWidget(arrLinPage);
    setupToggleBtn(btnArrLin, m_listArrLin);

    // Page 12: 圆形阵列 — 对象列表 + 轴/中心/角度/数量
    auto* arrCirPage = new QWidget(); auto* arrCirLay = new QVBoxLayout(arrCirPage);
    auto* btnArrCir = new QPushButton("+ 加入对象"); btnArrCir->setCheckable(true);
    arrCirLay->addWidget(btnArrCir); m_listArrCir = new QListWidget(); arrCirLay->addWidget(m_listArrCir);
    auto* arrCirFl = new QFormLayout(); arrCirLay->addLayout(arrCirFl);
    m_cmbCirAx = new QComboBox(); m_cmbCirAx->addItems({"X","Y","Z"});
    m_spCirCx = new QDoubleSpinBox(); m_spCirCx->setRange(-1e6,1e6); m_spCirCx->setDecimals(6); m_spCirCx->setValue(0); m_spCirCx->setSuffix(" m");
    m_spCirCy = new QDoubleSpinBox(); m_spCirCy->setRange(-1e6,1e6); m_spCirCy->setDecimals(6); m_spCirCy->setValue(0); m_spCirCy->setSuffix(" m");
    m_spCirCz = new QDoubleSpinBox(); m_spCirCz->setRange(-1e6,1e6); m_spCirCz->setDecimals(6); m_spCirCz->setValue(0); m_spCirCz->setSuffix(" m");
    m_spCirAng = new QDoubleSpinBox(); m_spCirAng->setRange(0,360); m_spCirAng->setDecimals(4); m_spCirAng->setValue(90); m_spCirAng->setSuffix(" deg");
    m_spCirCnt = new QSpinBox(); m_spCirCnt->setRange(2,100); m_spCirCnt->setValue(4);
    arrCirFl->addRow("旋转轴:", m_cmbCirAx); arrCirFl->addRow("中心 X:", m_spCirCx); arrCirFl->addRow("中心 Y:", m_spCirCy); arrCirFl->addRow("中心 Z:", m_spCirCz);
    arrCirFl->addRow("角度间距:", m_spCirAng); arrCirFl->addRow("数量:", m_spCirCnt);
    auto* chkArrCirKeep = new QCheckBox("保留输入对象"); chkArrCirKeep->setChecked(true);
    arrCirLay->addWidget(chkArrCirKeep);
    m_geomStackPages->addWidget(arrCirPage);
    setupToggleBtn(btnArrCir, m_listArrCir);

    m_chkKeepInputsU = new QCheckBox("保留输入对象"); m_chkKeepInputsU->setChecked(true);
    m_chkKeepInputsS = new QCheckBox("保留要添加的对象"); m_chkKeepInputsS->setChecked(true);
    m_chkKeepToolsS = new QCheckBox("保留要减去的对象"); m_chkKeepToolsS->setChecked(true);
    boolULayout->addWidget(m_chkKeepInteriorU);
    boolULayout->addWidget(m_chkKeepInputsU);
    boolSLayout->addWidget(m_chkKeepInteriorS);
    boolSLayout->addWidget(m_chkKeepInputsS);
    boolSLayout->addWidget(m_chkKeepToolsS);

    setupToggleBtn(m_btnAddUnion, m_listBoolUnion);
    setupToggleBtn(m_btnAddSubIn, m_listBoolSubIn);
    setupToggleBtn(m_btnAddSubTool, m_listBoolSubTool);

    // ---- 构建按钮 (逻辑在 buildGeometryFromUI, 与工具栏"全部构建"共用) ----
    connect(btnBuild, &QPushButton::clicked, this, [=]() {
        // 使用已记录的选中节点索引, 而非 currentItem (避免竞态)
        buildGeometryFromUI(m_lastGeomIdx);
    });

    setupBoundaryUI();
    setupSolverUI();
    setupPostProcessUI();
    setupParamUI();
    installBindButtons();

    // 分区控件: 插入到 pageMesh 布局中 btnExportMesh 之后
    {
        QWidget* partContainer = new QWidget();
        QVBoxLayout* partOuter = new QVBoxLayout(partContainer);
        partOuter->setContentsMargins(0, 4, 0, 0);
        QGroupBox* partGroup = new QGroupBox("并行网格分区");
        QVBoxLayout* partLay = new QVBoxLayout(partGroup);
        QWidget* partRow = new QWidget();
        QHBoxLayout* partRowLay = new QHBoxLayout(partRow);
        partRowLay->setContentsMargins(0, 0, 0, 0);
        partRowLay->addWidget(new QLabel("分区数:"));
        m_spinPartitions = new QSpinBox();
        m_spinPartitions->setRange(1, 128);
        int physCores = std::max(1, (int)std::thread::hardware_concurrency() / 2);
        m_spinPartitions->setValue(physCores);
        m_spinPartitions->setToolTip("与 MPI 并行进程数一致可最大化并行效率");
        partRowLay->addWidget(m_spinPartitions);
        m_btnPartition = new QPushButton("预分区网格 (并行加速)");
        m_btnPartition->setEnabled(false);
        m_btnPartition->setMinimumHeight(32);
        partRowLay->addWidget(m_btnPartition);
        partLay->addWidget(partRow);
        partOuter->addWidget(partGroup);

        QVBoxLayout* pageLayout = qobject_cast<QVBoxLayout*>(ui->pageMesh->layout());
        if (pageLayout) {
            int count = pageLayout->count();
            // 插入到倒数第2个位置 (spacer 之前)
            pageLayout->insertWidget(count - 1, partContainer);
        }
    }

    QWidget* emptyPage = new QWidget();
    ui->stackedWidgetSettings->addWidget(emptyPage);
    ui->stackedWidgetSettings->setCurrentIndex(ui->stackedWidgetSettings->count() - 1);

    connect(ui->cadViewerWidget, &CADViewer::faceSelectionChanged, this, &MainWindow::onFaceSelectionChanged);
    connect(ui->cadViewerWidget, &CADViewer::boxSelectionCompleted, this, [=](){ ui->btnBoxPick->setChecked(false); });
    // File menu actions
    connect(ui->actionImport, &QAction::triggered, this, &MainWindow::importGeometry);
    connect(ui->actionExportGeometry, &QAction::triggered, this, &MainWindow::exportGeometry);
    connect(ui->actionExportMesh, &QAction::triggered, this, &MainWindow::exportMesh);
    if (ui->actionImportMesh) connect(ui->actionImportMesh, &QAction::triggered, this, &MainWindow::importExternalMesh);

    ui->menuFile->addSeparator();
    QAction* actCaseLib = ui->menuFile->addAction("案例库...");
    connect(actCaseLib, &QAction::triggered, this, &MainWindow::showCaseLibrary);

    QAction* actSave = ui->menuFile->addAction("保存项目");
    actSave->setShortcut(QKeySequence::Save);
    connect(actSave, &QAction::triggered, this, [=](){
        if (m_projectFilePath.isEmpty()) {
            QString path = QFileDialog::getSaveFileName(this, "保存项目", "", "IonFlow Project (*.ion)");
            if (path.isEmpty()) return;
            m_projectFilePath = path;
        }
        saveProject(m_projectFilePath);
    });

    QAction* actSaveAs = ui->menuFile->addAction("另存为...");
    actSaveAs->setShortcut(QKeySequence("Ctrl+Shift+S"));
    connect(actSaveAs, &QAction::triggered, this, [=](){
        QString path = QFileDialog::getSaveFileName(this, "另存为", m_projectFilePath, "IonFlow Project (*.ion)");
        if (path.isEmpty()) return;
        m_projectFilePath = path;
        saveProject(path);
    });
    connect(ui->btnClearSelection, &QToolButton::clicked, ui->cadViewerWidget, &CADViewer::clearSelection);
    connect(ui->btnXY, &QToolButton::clicked, ui->cadViewerWidget, &CADViewer::setViewXY);
    connect(ui->btnYZ, &QToolButton::clicked, ui->cadViewerWidget, &CADViewer::setViewYZ);
    connect(ui->btnZX, &QToolButton::clicked, ui->cadViewerWidget, &CADViewer::setViewZX);
    connect(ui->btnResetHidden, &QToolButton::clicked, ui->cadViewerWidget, &CADViewer::resetHidden);

    // COMSOL 风格栅格: 显示切换按钮 + 栅格设置按钮 (插在 spacer 之前, 贴着左侧按钮)
    QToolButton* btnGrid = new QToolButton();
    btnGrid->setText("栅格");
    btnGrid->setCheckable(true);
    btnGrid->setChecked(false);
    btnGrid->setToolTip("显示/隐藏工作平面栅格");
    QToolButton* btnGridSettings = new QToolButton();
    btnGridSettings->setText("栅格设置");
    btnGridSettings->setToolTip("设置栅格平面/间距/细分");
    QHBoxLayout* toolsLayout = qobject_cast<QHBoxLayout*>(ui->frameGraphicsTools->layout());
    if (toolsLayout) {
        int spacerIdx = toolsLayout->count();
        for (int i = 0; i < toolsLayout->count(); i++)
            if (toolsLayout->itemAt(i)->spacerItem()) { spacerIdx = i; break; }
        toolsLayout->insertWidget(spacerIdx, btnGridSettings);
        toolsLayout->insertWidget(spacerIdx, btnGrid);
    }
    connect(btnGrid, &QToolButton::toggled, [=](bool on){ ui->cadViewerWidget->setGridVisible(on); });
    connect(btnGridSettings, &QToolButton::clicked, [=]() {
        QDialog dlg(this);
        dlg.setWindowTitle("栅格设置");
        auto* lay = new QVBoxLayout(&dlg);
        auto* form = new QFormLayout();

        auto* spinSpacing = new QDoubleSpinBox();
        spinSpacing->setRange(1e-6, 1e9);
        spinSpacing->setDecimals(6);
        spinSpacing->setValue(ui->cadViewerWidget->gridSpacing());

        auto* spinSubdiv = new QSpinBox();
        spinSubdiv->setRange(1, 10);
        spinSubdiv->setValue(ui->cadViewerWidget->gridSubdivisions());

        form->addRow("主栅格间距:", spinSpacing);
        form->addRow("每主间距细分:", spinSubdiv);
        lay->addLayout(form);

        auto* btnBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        connect(btnBox, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
        connect(btnBox, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
        lay->addWidget(btnBox);

        if (dlg.exec() == QDialog::Accepted) {
            ui->cadViewerWidget->setGridSpacing(spinSpacing->value());
            ui->cadViewerWidget->setGridSubdivisions(spinSubdiv->value());
        }
    });

    // 选择模式: 边界选择 / 对象选择
    m_geomSelectMode = new QComboBox();
    m_geomSelectMode->addItem("边界选择"); m_geomSelectMode->addItem("对象选择"); m_geomSelectMode->addItem("边选择");
    m_geomSelectMode->setMaximumWidth(100);
    if (toolsLayout) {
        int si = toolsLayout->count();
        for (int i = 0; i < toolsLayout->count(); i++)
            if (toolsLayout->itemAt(i)->spacerItem()) { si = i; break; }
        toolsLayout->insertWidget(si, m_geomSelectMode);
    }

    // ============================================================
    // 视图/选择功能栏: HBox → FlowLayout
    // 画布变窄时按钮自动换行成多行, 不再撑大 centralwidget 最小宽度
    // (解除画布左右缩放限制; 所有控件已就位后再替换布局)
    // ============================================================
    {
        QLayout* oldLay = ui->frameGraphicsTools->layout();
        QList<QWidget*> widgets;
        while (oldLay && oldLay->count() > 0) {
            QLayoutItem* it = oldLay->takeAt(0);
            if (QWidget* w = it->widget()) widgets.append(w);
            delete it;   // 丢弃尾部 spacer: FlowLayout 无需右端撑开
        }
        delete oldLay;
        auto* toolsFlow = new FlowLayout(ui->frameGraphicsTools, 6, 6);
        toolsFlow->setContentsMargins(8, 4, 8, 4);   // 与原 .ui margin 一致
        for (QWidget* w : widgets) toolsFlow->addWidget(w);
    }

    connect(m_geomSelectMode, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int idx){
        ui->cadViewerWidget->setObjectSelectMode(idx == 1);
        ui->cadViewerWidget->setEdgeSelectMode(idx == 2);
    });

    connect(ui->btnBoxPick, &QToolButton::toggled, this, &MainWindow::updateInteractionModes);
    connect(ui->btnHideMode, &QToolButton::toggled, this, &MainWindow::updateInteractionModes);
    // 边选择模式下点击隐藏自动切换为面选择模式 (边不可隐藏)
    connect(ui->btnHideMode, &QToolButton::toggled, [=](bool on){
        if (on && m_geomSelectMode && m_geomSelectMode->currentIndex() == 2)
            m_geomSelectMode->setCurrentIndex(0);
    });
    connect(ui->comboSelectMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::updateInteractionModes);
    connect(ui->comboViewMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [=](int index){
        if (index == 0) ui->cadViewerWidget->showUnhidden();
        else if (index == 1) ui->cadViewerWidget->showOnlyHidden();
        else if (index == 2) ui->cadViewerWidget->showAll();
    });
    connect(ui->treeGroups, &QTreeWidget::customContextMenuRequested, this, &MainWindow::showTreeContextMenu);
    connect(ui->treeGroups, &QTreeWidget::currentItemChanged, this, &MainWindow::onTreeSelectionChanged);
    connect(ui->treeGroups, &QTreeWidget::itemChanged, this, [this](QTreeWidgetItem* item, int col) {
        if (col != 0 || !item || item->parent() != m_rootGroupNode) return;
        // 重命名持久化: 更新 m_importedGroupMeta 中存储的组名
        if (m_activeMeshNode && m_importedGroupMeta.count(m_activeMeshNode)) {
            QString newName = item->text(0);
            for (auto& gm : m_importedGroupMeta[m_activeMeshNode])
                if (m_groupTags.count(item) && gm.tag == m_groupTags[item])
                    { gm.name = newName; break; }
            m_projectModified = true;
        }
    });
    connect(ui->btnClear, &QPushButton::clicked, ui->cadViewerWidget, &CADViewer::clearSelection);
    connect(ui->editGroupName, &QLineEdit::textEdited, this, &MainWindow::onGroupNameEdited);
    connect(ui->btnGenerateMesh, &QPushButton::clicked, this, &MainWindow::generateMesh);
    connect(ui->btnExportMesh, &QPushButton::clicked, this, &MainWindow::exportMesh);
    connect(ui->comboMeshGroupTarget, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onMeshGroupTargetChanged);
    connect(ui->chkEnableLocalField, &QCheckBox::toggled, this, &MainWindow::onEnableLocalFieldToggled);
    connect(ui->spinSizeMin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onMeshParamsChanged);
    connect(ui->spinSizeMax, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onMeshParamsChanged);
    connect(ui->spinDistMin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onMeshParamsChanged);
    connect(ui->spinDistMax, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onMeshParamsChanged);

    // 分区按钮
    connect(m_btnPartition, &QPushButton::clicked, [this]() {
        if (!m_activeMeshNode || !m_meshFiles.count(m_activeMeshNode)) return;

        // 分区源: 导入网格先同步到缓存 (物理组一致 + 路径与计算 mesh_file 一致),
        // 内部网格直接用缓存文件
        QString mshPath;
        bool isImported = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode)
                          && m_isImportedMesh[m_activeMeshNode];
        if (isImported) {
            mshPath = syncImportedMeshToCache();
            if (mshPath.isEmpty()) {
                logMessage("同步导入网格失败，无法分区", "error");
                return;
            }
        } else {
            mshPath = m_meshFiles[m_activeMeshNode];
        }

        int nParts = m_spinPartitions->value();
        QFileInfo fi(mshPath);
        QString prefix = fi.absolutePath() + "/" + fi.completeBaseName() + "_part";

        m_btnPartition->setEnabled(false);
        logMessage(QString("Gmsh 分区中 (%1 分区)...").arg(nParts), "info");

        // Step 1: gmsh CLI 分区 (等效 CLI: gmsh mesh.msh -save -format msh2
        //          -o mesh_part.msh -part N -part_split -part_physicals)
        auto* step1 = new QProcess(this);
        step1->setProcessChannelMode(QProcess::MergedChannels);
        connect(step1, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
                [=](int code, QProcess::ExitStatus) {
            if (code != 0) {
                QString err = QString::fromLocal8Bit(step1->readAll());
                logMessage("Gmsh 分区失败: " + err.left(200), "error");
                m_btnPartition->setEnabled(true);
                step1->deleteLater();
                return;
            }
            logMessage("Gmsh 分区完成，转换 MFEM 并行格式...", "info");

            // Step 2: gmsh_partition_to_mfem 转换
            QString converter = QApplication::applicationDirPath()
                                + "/gmsh_partition_to_mfem";
            if (!QFile::exists(converter)) {
                logMessage("未找到 gmsh_partition_to_mfem 工具，仅完成 Gmsh 分区", "warning");
                logMessage(QString("分区文件: %1_1.msh ~ %1_%2.msh")
                           .arg(prefix).arg(nParts), "info");
                m_btnPartition->setEnabled(true);
                step1->deleteLater();
                return;
            }

            QString mfemPrefix = mshPath + ".part";
            QStringList args;
            args << "-np" << QString::number(nParts) << converter;
            for (int i = 1; i <= nParts; i++)
                args << (prefix + "_" + QString::number(i) + ".msh");
            args << "-o" << mfemPrefix;

            auto* step2 = new QProcess(this);
            connect(step2, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
                    [=](int code2, QProcess::ExitStatus) {
                if (code2 == 0)
                    logMessage(QString("预分区完成 → %1.part.*").arg(mshPath), "success");
                else {
                    QString err2 = QString::fromLocal8Bit(step2->readAll());
                    logMessage("MFEM 格式转换失败: " + err2.left(200), "error");
                }
                m_btnPartition->setEnabled(true);
                step2->deleteLater();
            });
            step2->start("mpirun", args);
            step1->deleteLater();
        });
        step1->start("gmsh", QStringList() << mshPath
                     << "-save" << "-format" << "msh2"
                     << "-o" << (prefix + ".msh")
                     << "-part" << QString::number(nParts)
                     << "-part_split" << "-part_physicals");
    });

    logMessage("系统初始化完成。欢迎使用离子流场模拟器。", "info");
    QString geomIconDir = QApplication::applicationDirPath() + "/icons";
    logMessage(QString("几何图标目录: %1 (可放入同名 .png/.svg/.jpg 自定义图标)").arg(geomIconDir), "info");
}

void MainWindow::setProjectPath(const QString& path) {
    if (path == "__untitled__") {
        m_projectFilePath.clear();  // 未保存, 保存时才选路径
        m_projectModified = false;
        setWindowTitle("离子流场模拟器 - 未命名项目");
        return;
    }
    m_projectFilePath = path;
    QFileInfo fi(path);
    setWindowTitle(QString("离子流场模拟器 - %1").arg(fi.fileName()));
    if (fi.exists() && fi.size() > 0) {
        loadProject(path);
    } else {
        m_projectModified = false;
    }
}

MainWindow::~MainWindow() {
    if (m_solverProcess && m_solverProcess->state() == QProcess::Running) {
        m_solverProcess->kill();
        m_solverProcess->waitForFinished();
    }
    delete ui;
}

void MainWindow::logMessage(const QString& msg, const QString& level) {
    QString timeStr = QDateTime::currentDateTime().toString("HH:mm:ss");
    QString color = "#CCCCCC"; 
    
    if (level == "success") color = "#4CAF50"; 
    else if (level == "error") color = "#F44336"; 
    else if (level == "warning") color = "#FFC107"; 
    else if (level == "highlight") color = "#64B5F6"; 
    
    QString html = QString("<span style='color:#888888;'>[%1]</span> <span style='color:%2;'>%3</span>")
                    .arg(timeStr).arg(color).arg(msg.toHtmlEscaped().replace("\n", "<br>"));
    m_textConsole->append(html);
}

// ====================================================================
// 几何创建工具栏 (COMSOL 风格: 构建 | 导入/导出 | 体素 | 操作 | 变换)
// ====================================================================
void MainWindow::setupGeometryToolbar() {
    auto* geomToolbar = new QToolBar("几何操作", this);
    geomToolbar->setMovable(false);
    geomToolbar->setObjectName("geomToolbar");
    geomToolbar->setContentsMargins(4, 0, 4, 0);   // 去掉工具栏自身上下留白, 内容贴近上边界
    this->addToolBar(Qt::TopToolBarArea, geomToolbar);

    // ---- 面板统一样式 ----
    const QString kPanelStyle = R"(
        QFrame#geomPanel { background-color: #FFFFFF; border: 1px solid #DEE2E6;
                           border-radius: 8px; padding: 2px 6px 2px 6px; }
        QFrame#geomPanel QToolButton { background: transparent; border: none;
                                       border-radius: 6px; padding: 0px 5px;
                                       font-size: 13px; color: #212529; }
        QFrame#geomPanel QToolButton:hover { background-color: #E7F1FF; }
        QFrame#geomPanel QToolButton:pressed { background-color: #C5E1FF; }
        QFrame#geomPanel QToolButton::menu-indicator { image: none; }
        QLabel#geomPanelTitle { font-size: 12px; color: #6C757D; font-weight: bold;
                                margin-top: 2px; }
    )";

    // 普通按钮: 图标在左、文字在右 (仅"更多体素"为图标在上文字在下)
    auto makeToolBtn = [](const QString& text, const QString& iconName) {
        auto* btn = new QToolButton();
        btn->setText(" " + text);   // 图标与文字之间留少许间隙
        btn->setIcon(geomIcon(iconName));
        btn->setIconSize(QSize(16, 16));
        btn->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        btn->setAutoRaise(true);
        btn->setMinimumWidth(72);
        btn->setMinimumHeight(34);
        btn->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        return btn;
    };

    // 面板容器: 白底圆角卡片, 底部居中标题
    auto makePanel = [&](const QString& title) {
        auto* panel = new QFrame();
        panel->setObjectName("geomPanel");
        panel->setStyleSheet(kPanelStyle);
        auto* vlay = new QVBoxLayout(panel);
        vlay->setContentsMargins(6, 2, 6, 2);
        vlay->setSpacing(2);
        auto* titleLbl = new QLabel(title);
        titleLbl->setObjectName("geomPanelTitle");
        titleLbl->setAlignment(Qt::AlignCenter);
        vlay->addWidget(titleLbl);
        return panel;
    };

    std::map<QString,int> counters;
    // 新节点: 若无已构建节点则追加到末尾, 否则紧跟最后一个已构建节点
    auto insertAfterLastBuilt = [&](const QString& typeName, const QString& displayName) {
        int lastBuilt = -1;
        for (int i = 0; i < (int)m_geomNodes.size(); i++)
            if (m_geomNodes[i].built) lastBuilt = i;
        int insertIdx = (lastBuilt < 0) ? (int)m_geomNodes.size() : (lastBuilt + 1);
        m_geomNodes.insert(m_geomNodes.begin() + insertIdx,
                           {typeName, nullptr, false});
        // 倒圆角/倒斜角节点: 设置合理的默认参数
        if (typeName == "倒圆角" || typeName == "倒斜角")
            m_geomNodes[insertIdx].w = 5.0;
        auto* item = new QTreeWidgetItem();
        item->setText(0, displayName);
        item->setIcon(0, geomIcon(typeName));
        item->setData(0, Qt::UserRole, insertIdx);
        item->setData(0, Qt::UserRole+1, typeName);
        item->setForeground(0, QBrush(QColor(128,128,128)));
        m_geomNodes[insertIdx].treeItem = item;
        // 同步树节点顺序到 m_geomNodes 顺序
        m_rootGeomNode->insertChild(insertIdx, item);
        for (int i = insertIdx + 1; i < (int)m_geomNodes.size(); i++)
            m_geomNodes[i].treeItem->setData(0, Qt::UserRole, i);
        m_rootGeomNode->setExpanded(true);
        ui->treeGroups->setCurrentItem(item);
        m_projectModified = true;
    };

    QStringList morePrims = {"椭球","棱柱","棱锥","圆台","棱台"};

    // ============================================================
    // 构建面板: 单个大按钮 "全部构建" (图标在上、文字在下, 同"更多体素")
    // 功能 = 构建设置的全部几何 (等价于对最后一个节点执行"构建/更新几何体")
    // ============================================================
    auto* panelBuild = makePanel("构建");
    auto* buildAllBtn = new QToolButton();
    buildAllBtn->setText("全部构建");
    buildAllBtn->setIcon(geomIcon("全部构建"));
    buildAllBtn->setIconSize(QSize(28, 28));
    buildAllBtn->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
    buildAllBtn->setAutoRaise(true);
    buildAllBtn->setMinimumWidth(84);
    buildAllBtn->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    connect(buildAllBtn, &QToolButton::clicked, this, [=]() {
        // 保存最后一个节点的 UI 参数后链式构建全部节点 (空树时自动跳过)
        buildGeometryFromUI((int)m_geomNodes.size() - 1);
    });
    auto* buildLay = new QVBoxLayout();
    buildLay->addWidget(buildAllBtn, 1);   // 拉伸占满两行高度
    qobject_cast<QVBoxLayout*>(panelBuild->layout())->insertLayout(0, buildLay);

    // ============================================================
    // 导入/导出面板: 上下两个横排按钮 (图标在左、文字在右)
    // 功能与 File 菜单的导入/导出几何一致
    // ============================================================
    auto* panelIO = makePanel("导入/导出");
    auto* importBtn = makeToolBtn("导入几何", "导入几何");
    connect(importBtn, &QToolButton::clicked, this, &MainWindow::importGeometry);
    auto* exportBtn = makeToolBtn("导出几何", "导出几何");
    connect(exportBtn, &QToolButton::clicked, this, &MainWindow::exportGeometry);
    auto* ioLay = new QVBoxLayout();
    ioLay->setSpacing(0);
    ioLay->addWidget(importBtn);
    ioLay->addWidget(exportBtn);
    qobject_cast<QVBoxLayout*>(panelIO->layout())->insertLayout(0, ioLay);

    // ============================================================
    // 体素面板: 3 列网格 (左列/中列 6 个 + 右列"更多体素"大按钮)
    // ============================================================
    auto* panelVoxel = makePanel("体素");
    auto* voxelH = new QHBoxLayout();
    voxelH->setSpacing(2);
    qobject_cast<QVBoxLayout*>(panelVoxel->layout())->insertLayout(0, voxelH);

    // 2 行 x 3 列: 上行 立方体/球体/圆柱, 下行 圆锥/圆环/楔形
    QStringList prims = {"立方体","球体","圆柱","圆锥","圆环","楔形"};
    auto* primGrid = new QGridLayout();
    primGrid->setSpacing(0);
    for (int i = 0; i < prims.size(); ++i) {
        auto* btn = makeToolBtn(prims[i], prims[i]);
        connect(btn, &QToolButton::clicked, [=]() mutable {
            counters[prims[i]]++;
            insertAfterLastBuilt(prims[i], QString("%1 %2").arg(prims[i]).arg(counters[prims[i]]));
        });
        primGrid->addWidget(btn, i / 3, i % 3);
    }
    voxelH->addLayout(primGrid);

    // 右列: 放大占满两行高的"更多体素"下拉按钮 (立方体图形 + 下方文字)
    auto* moreBtn = new QToolButton();
    moreBtn->setText("更多体素 ▾");
    moreBtn->setIcon(geomIcon("立方体"));
    moreBtn->setIconSize(QSize(28, 28));
    moreBtn->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
    moreBtn->setAutoRaise(true);
    moreBtn->setMinimumWidth(84);
    moreBtn->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    auto* moreMenu = new QMenu(this);
    for (auto& n : morePrims) {
        auto* act = moreMenu->addAction(geomIcon(n), n, [=]() mutable {
            counters[n]++; insertAfterLastBuilt(n, QString("%1 %2").arg(n).arg(counters[n]));
        });
        act->setIcon(geomIcon(n));
    }
    moreBtn->setMenu(moreMenu);
    moreBtn->setPopupMode(QToolButton::InstantPopup);
    auto* rightLay = new QVBoxLayout();
    rightLay->addWidget(moreBtn, 1);   // 拉伸占满三行高度
    voxelH->addLayout(rightLay, 1);

    // ============================================================
    // 操作面板: 2 列 x 3 行网格
    // ============================================================
    auto* panelOp = makePanel("操作");
    auto* opGrid = new QGridLayout();
    opGrid->setSpacing(0);
    QStringList ops = {"并集","差集","交集","倒圆角","倒斜角","扫掠"};
    for (int i = 0; i < ops.size(); ++i) {
        auto* btn = makeToolBtn(ops[i], ops[i]);
        connect(btn, &QToolButton::clicked, [=]() mutable {
            counters[ops[i]]++;
            insertAfterLastBuilt(ops[i], QString("%1 %2").arg(ops[i]).arg(counters[ops[i]]));
        });
        opGrid->addWidget(btn, i / 3, i % 3);
    }
    qobject_cast<QVBoxLayout*>(panelOp->layout())->insertLayout(0, opGrid);

    // ============================================================
    // 变换面板: 3 列 x 3 行网格 (末行"圆形阵列"居中)
    // ============================================================
    auto* panelXfm = makePanel("变换");
    auto* xfmGrid = new QGridLayout();
    xfmGrid->setSpacing(0);
    QStringList xfms = {"平移","旋转","镜像","缩放","偏移","线性阵列","圆形阵列"};
    for (int i = 0; i < xfms.size(); ++i) {
        auto* btn = makeToolBtn(xfms[i], xfms[i]);
        connect(btn, &QToolButton::clicked, [=]() mutable {
            counters[xfms[i]]++;
            insertAfterLastBuilt(xfms[i], QString("%1 %2").arg(xfms[i]).arg(counters[xfms[i]]));
        });
        int row = i / 4, col = i % 4;
        xfmGrid->addWidget(btn, row, col);
    }
    qobject_cast<QVBoxLayout*>(panelXfm->layout())->insertLayout(0, xfmGrid);

    // ============================================================
    // 容器: 五个面板并排, 窗口变窄时水平滚动而非截断
    // ============================================================
    auto* scroll = new QScrollArea();
    scroll->setWidgetResizable(false);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    scroll->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scroll->setStyleSheet("QScrollArea { background: transparent; border: none; }");
    auto* content = new QWidget();
    auto* contentH = new QHBoxLayout(content);
    contentH->setContentsMargins(0, 0, 0, 0);
    contentH->setSpacing(6);
    contentH->addWidget(panelBuild);
    contentH->addWidget(panelIO);
    contentH->addWidget(panelVoxel);
    contentH->addWidget(panelOp);
    contentH->addWidget(panelXfm);
    scroll->setWidget(content);
    geomToolbar->addWidget(scroll);

    ui->cadViewerWidget->GetRenderer()->GetRenderWindow()->Render();
}

// ====================================================================
// 【核心修改】：重构后处理 UI，支持模式切换与动态显示
// ====================================================================
void MainWindow::setupPostProcessUI() {
    m_pagePostProcess = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_pagePostProcess);
    QScrollArea* scrollArea = new QScrollArea();
    scrollArea->setWidgetResizable(true);
    scrollArea->setFrameShape(QFrame::NoFrame);
    QWidget* scrollContent = new QWidget();
    QVBoxLayout* scrollLayout = new QVBoxLayout(scrollContent);
    scrollLayout->setSpacing(12); 

    // 1. 数据源选择
    QGroupBox* groupData = new QGroupBox("可视化数据源");
    QFormLayout* formMap = new QFormLayout(groupData);
    m_comboPostField = new QComboBox();
    m_comboPostField->addItem("空间电荷密度 rho (C/m^3/ε₀)");
    m_comboPostField->addItem("电场强度幅度 |E| (V/m)");
    m_comboPostField->addItem("电势 φ (V)");
    formMap->addRow("查看物理场:", m_comboPostField);
    scrollLayout->addWidget(groupData);

    // 2. 渲染模式控制
    QGroupBox* groupRender = new QGroupBox("视图渲染模式");
    QVBoxLayout* renderLayout = new QVBoxLayout(groupRender);
    m_radioSurface = new QRadioButton("表面着色云图 (Surface)");
    m_radioVolume = new QRadioButton("三维切面/等值面 (Volume Slice)");
    m_radioSurface->setChecked(true);
    renderLayout->addWidget(m_radioSurface);
    renderLayout->addWidget(m_radioVolume);
    scrollLayout->addWidget(groupRender);

    // 3. 表面模式参数：多选物理组列表
    QGroupBox* groupSelect = new QGroupBox("渲染目标物理组");
    QVBoxLayout* selLayout = new QVBoxLayout(groupSelect);
    m_listPostGroups = new QListWidget();
    m_listPostGroups->setSelectionMode(QAbstractItemView::MultiSelection);
    m_listPostGroups->setStyleSheet("QListWidget { border: 1px solid #CED4DA; padding: 4px; } "
                                    "QListWidget::item { padding: 4px; }");
    selLayout->addWidget(m_listPostGroups);
    scrollLayout->addWidget(groupSelect);

    // 4. 体积模式参数：三维切面位置控制
    QGroupBox* groupVolume = new QGroupBox("三维切面参数");
    QFormLayout* formVol = new QFormLayout(groupVolume);

    m_chkSliceX = new QCheckBox("启用 X 切面");
    m_spinSliceX = new QDoubleSpinBox(); m_spinSliceX->setRange(-1e6, 1e6); m_spinSliceX->setValue(0.0);
    QHBoxLayout* rowX = new QHBoxLayout();
    rowX->addWidget(m_chkSliceX); rowX->addWidget(m_spinSliceX);
    formVol->addRow("X 切面:", rowX);

    m_chkSliceY = new QCheckBox("启用 Y 切面");
    m_spinSliceY = new QDoubleSpinBox(); m_spinSliceY->setRange(-1e6, 1e6); m_spinSliceY->setValue(0.0);
    QHBoxLayout* rowY = new QHBoxLayout();
    rowY->addWidget(m_chkSliceY); rowY->addWidget(m_spinSliceY);
    formVol->addRow("Y 切面:", rowY);

    m_chkSliceZ = new QCheckBox("启用 Z 切面");
    m_spinSliceZ = new QDoubleSpinBox(); m_spinSliceZ->setRange(-1e6, 1e6); m_spinSliceZ->setValue(0.0);
    QHBoxLayout* rowZ = new QHBoxLayout();
    rowZ->addWidget(m_chkSliceZ); rowZ->addWidget(m_spinSliceZ);
    formVol->addRow("Z 切面:", rowZ);

    scrollLayout->addWidget(groupVolume);

    // 初始状态下隐藏体积渲染参数
    groupVolume->setVisible(false);

    // UI 状态机联动：根据单选框动态显示对应面板
    connect(m_radioSurface, &QRadioButton::toggled, [=](bool checked){
        groupSelect->setVisible(checked);
        groupVolume->setVisible(!checked);
    });

    // 5. 颜色范围控制
    QGroupBox* groupColorRange = new QGroupBox("颜色范围");
    QFormLayout* formCR = new QFormLayout(groupColorRange);
    m_chkAutoRange = new QCheckBox("自动范围");
    m_chkAutoRange->setChecked(true);
    m_spinRangeMin = new QDoubleSpinBox(); m_spinRangeMin->setRange(-1e12, 1e12);
    m_spinRangeMin->setDecimals(3); m_spinRangeMin->setValue(0.0); m_spinRangeMin->setEnabled(false);
    m_spinRangeMax = new QDoubleSpinBox(); m_spinRangeMax->setRange(-1e12, 1e12);
    m_spinRangeMax->setDecimals(3); m_spinRangeMax->setValue(1.0); m_spinRangeMax->setEnabled(false);
    formCR->addRow(m_chkAutoRange);
    formCR->addRow("最小值:", m_spinRangeMin);
    formCR->addRow("最大值:", m_spinRangeMax);
    connect(m_chkAutoRange, &QCheckBox::toggled, [=](bool c) {
        m_spinRangeMin->setEnabled(!c); m_spinRangeMax->setEnabled(!c);
    });
    scrollLayout->addWidget(groupColorRange);

    // 6. 操作区
    QGroupBox* groupAction = new QGroupBox("执行操作");
    QVBoxLayout* actionLayout = new QVBoxLayout(groupAction);
    m_btnRenderPost = new QPushButton("加载并执行渲染");
    m_btnRenderPost->setObjectName("btnRenderPost");
    m_btnRenderPost->setMinimumHeight(45);
    actionLayout->addWidget(m_btnRenderPost);
    scrollLayout->addWidget(groupAction);

    scrollLayout->addStretch(); 

    scrollContent->setLayout(scrollLayout);
    scrollArea->setWidget(scrollContent);
    mainLayout->addWidget(scrollArea);

    ui->stackedWidgetSettings->addWidget(m_pagePostProcess);

    // ---------- 收集选中物理组的 Gmsh tag + 获取当前网格文件 ----------
    auto getSelectedTagsAndMesh = [this]() -> std::pair<std::set<int>, QString> {
        std::set<int> tags;
        QString mshPath;
        for (int i = 0; i < m_listPostGroups->count() && i < m_rootGroupNode->childCount(); ++i) {
            if (m_listPostGroups->item(i)->checkState() != Qt::Checked) continue;
            QTreeWidgetItem* treeItem = m_rootGroupNode->child(i);
            if (m_groupTags.count(treeItem))
                tags.insert(m_groupTags[treeItem]);
        }
        for (auto& [node, path] : m_meshFiles)
            if (!path.isEmpty()) { mshPath = path; break; }
        return {tags, mshPath};
    };

    auto executeRender = [this, getSelectedTagsAndMesh]() {
        if (!m_radioSurface->isChecked() || m_lastVtkFilePath.isEmpty() || m_lastFieldName.isEmpty()) return;
        if (!ui->cadViewerWidget->hasResultActor()) return;
        auto [tags, mshPath] = getSelectedTagsAndMesh();
        ui->cadViewerWidget->loadVtkContourSurface(m_lastVtkFilePath, m_lastFieldName, mshPath, tags, 1.0, 0.0);
    };

    // 渲染按钮
    connect(m_btnRenderPost, &QPushButton::clicked, this, [=](){
        logMessage("正在准备映射后处理数据...", "highlight");

        QString workingDir = QApplication::applicationDirPath();
        QString vtkFilePath = workingDir + "/results/Takuma_Results/Cycle000000/data.pvtu";
        if (!QFileInfo(vtkFilePath).exists()) {
            vtkFilePath = workingDir + "/results/Takuma_Results/Takuma_Results_000000/data.pvtu";
            if (!QFileInfo(vtkFilePath).exists()) {
                QMessageBox::warning(this, "警告", "未找到 VTK 结果文件！");
                return;
            }
        }
        QString txt = m_comboPostField->currentText();
        QString fieldName = txt.contains("|E|") ? "E_scalar" : (txt.contains("φ") || txt.contains("电势")) ? "phi" : "rho";

        if (m_radioSurface->isChecked()) {
            auto [selectedTags, mshPath] = getSelectedTagsAndMesh();
            QStringList selectedNames;
            for (int i = 0; i < m_listPostGroups->count() && i < m_rootGroupNode->childCount(); ++i)
                if (m_listPostGroups->item(i)->checkState() == Qt::Checked)
                    selectedNames << m_listPostGroups->item(i)->text();
            if (selectedNames.isEmpty()) {
                QMessageBox::warning(this, "警告", "请至少勾选一个物理组进行表面渲染！");
                return;
            }
            logMessage(QString("已选择的渲染面: %1").arg(selectedNames.join(", ")), "info");
            m_lastVtkFilePath = vtkFilePath;
            m_lastFieldName = fieldName;
            ui->cadViewerWidget->clearMesh();
            // 导入网格: 从 m_groupData 收集 entity tags 直接传给渲染, 避免读取旧 .msh 中的 PG
            std::set<int> directEntityTags;
            bool isImported = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode)
                              && m_isImportedMesh[m_activeMeshNode];
            if (isImported && !m_faceToEntity.empty()) {
                for (int i = 0; i < m_listPostGroups->count() && i < m_rootGroupNode->childCount(); ++i) {
                    if (m_listPostGroups->item(i)->checkState() == Qt::Checked) {
                        for (int fid : m_groupData[m_rootGroupNode->child(i)])
                            if (m_faceToEntity.count(fid))
                                directEntityTags.insert(m_faceToEntity[fid]);
                    }
                }
            }
            ui->cadViewerWidget->loadVtkContourSurface(vtkFilePath, fieldName, mshPath, selectedTags,
                                                        1.0, 0.0, directEntityTags);
            ui->cadViewerWidget->setGeometryHidden(true);
            ui->cadViewerWidget->setResultVisible(true);
            logMessage("已成功加载 VTK 结果并渲染表面云图！", "success");
        } else {
            // 三维切面模式
            bool sx = m_chkSliceX->isChecked(), sy = m_chkSliceY->isChecked(), sz = m_chkSliceZ->isChecked();
            if (!sx && !sy && !sz) {
                QMessageBox::warning(this, "警告", "请至少启用一个方向的切面！");
                return;
            }
            m_lastVtkFilePath = vtkFilePath;
            m_lastFieldName = fieldName;
            ui->cadViewerWidget->clearMesh();  // 清除旧网格显示
            ui->cadViewerWidget->loadVtkSliceView(vtkFilePath, fieldName,
                sx, m_spinSliceX->value(), sy, m_spinSliceY->value(), sz, m_spinSliceZ->value());
            ui->cadViewerWidget->setGeometryHidden(true);
            ui->cadViewerWidget->setResultVisible(true);
            logMessage("已成功渲染三维切面云图！", "success");
        }
    });

    // 物理场切换不自动渲染, 由用户手动点击渲染按钮
}

// 【新增】：刷新后处理界面的物理组列表
void MainWindow::refreshPostProcessGroups() {
    m_listPostGroups->clear();
    for(int i = 0; i < m_rootGroupNode->childCount(); ++i) {
        QTreeWidgetItem* node = m_rootGroupNode->child(i);
        QListWidgetItem* item = new QListWidgetItem(node->text(0), m_listPostGroups);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Checked); // 默认全选
    }
}

void MainWindow::setupBoundaryUI() {
    m_pageBoundary = new QWidget();
    auto* mainLayout = new QVBoxLayout(m_pageBoundary);
    auto* scroll = new QScrollArea(); scroll->setWidgetResizable(true); scroll->setFrameShape(QFrame::NoFrame);
    auto* content = new QWidget(); auto* lay = new QVBoxLayout(content);

    QGroupBox* groupBdr = new QGroupBox("边界条件映射");
    QFormLayout* layBdr = new QFormLayout(groupBdr);
    m_comboBoundaryGroup = new QComboBox();
    layBdr->addRow("目标物理组:", m_comboBoundaryGroup);
    m_chkApplyBoundary = new QCheckBox("施加边界条件");
    layBdr->addRow("", m_chkApplyBoundary);
    // 电压输入: 常数/函数 模式切换
    auto* voltageWrap = new QWidget();
    auto* voltageHLay = new QHBoxLayout(voltageWrap);
    voltageHLay->setContentsMargins(0,0,0,0);
    m_cmbBdrMode = new QComboBox();
    m_cmbBdrMode->addItems({"常数", "函数"});
    voltageHLay->addWidget(m_cmbBdrMode);
    m_stackBdrInput = new QStackedWidget();
    m_spinBdrVoltage = new QDoubleSpinBox();
    m_spinBdrVoltage->setRange(-1e9, 1e9); m_spinBdrVoltage->setDecimals(6); m_spinBdrVoltage->setSuffix(" V");
    m_stackBdrInput->addWidget(m_spinBdrVoltage);
    m_editBdrExpr = new QLineEdit();
    m_editBdrExpr->setPlaceholderText("例: 1000*x + 500*y");
    m_stackBdrInput->addWidget(m_editBdrExpr);
    voltageHLay->addWidget(m_stackBdrInput);
    auto* bindBtn = new QToolButton();
    bindBtn->setText("="); bindBtn->setMaximumSize(22, 22);
    bindBtn->setToolTip("绑定/解绑参数");
    bindBtn->setStyleSheet("QToolButton { border:1px solid #CCC; background:#EEE; font-weight:bold; }");
    if (m_paramBindings.count(m_spinBdrVoltage)) {
        m_spinBdrVoltage->setStyleSheet("QDoubleSpinBox { background-color: #E3F2FD; }");
        m_spinBdrVoltage->setReadOnly(true);
        bindBtn->setStyleSheet("QToolButton { border:1px solid #2196F3; background:#BBDEFB; font-weight:bold; color:#1565C0; }");
    }
    connect(bindBtn, &QToolButton::clicked, [this, bindBtn](){
        showBindDialog(m_spinBdrVoltage, bindBtn);
        // 更新按钮样式
        if (m_paramBindings.count(m_spinBdrVoltage)) {
            bindBtn->setStyleSheet("QToolButton { border:1px solid #2196F3; background:#BBDEFB; font-weight:bold; color:#1565C0; }");
        } else {
            bindBtn->setStyleSheet("QToolButton { border:1px solid #CCC; background:#EEE; font-weight:bold; }");
        }
    });
    voltageHLay->addWidget(bindBtn);
    layBdr->addRow("电压:", voltageWrap);
    auto* hintLabel = new QLabel("变量: x, y, z  支持: + - * / ^ sqrt sin cos abs min max");
    hintLabel->setStyleSheet("color:#888; font-size:13px;");
    hintLabel->setVisible(false);
    layBdr->addRow("", hintLabel);
    connect(m_cmbBdrMode, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int idx){
        m_stackBdrInput->setCurrentIndex(idx);
        hintLabel->setVisible(idx == 1);
        if (!m_isUpdatingUI) onBoundaryParamsChanged();
    });
    m_chkBdrCorona = new QCheckBox("作为电晕极");
    layBdr->addRow("", m_chkBdrCorona);
    lay->addWidget(groupBdr);
    lay->addStretch();
    content->setLayout(lay); scroll->setWidget(content);
    mainLayout->addWidget(scroll);
    ui->stackedWidgetSettings->addWidget(m_pageBoundary);

    connect(m_editBdrExpr, &QLineEdit::textChanged, [this](const QString&){
        if (!m_isUpdatingUI) onBoundaryParamsChanged();
    });
    connect(m_comboBoundaryGroup, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onBoundaryGroupTargetChanged);
    connect(m_chkApplyBoundary, &QCheckBox::toggled, this, &MainWindow::onBoundaryParamsChanged);
    connect(m_spinBdrVoltage, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::onBoundaryParamsChanged);
    connect(m_chkBdrCorona, &QCheckBox::toggled, this, &MainWindow::onBoundaryParamsChanged);
}

void MainWindow::setupSolverUI() {
    m_pageSolver = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_pageSolver);
    QScrollArea* scrollArea = new QScrollArea();
    scrollArea->setWidgetResizable(true);
    scrollArea->setFrameShape(QFrame::NoFrame); 
    QWidget* scrollContent = new QWidget();
    QVBoxLayout* scrollLayout = new QVBoxLayout(scrollContent);
    scrollLayout->setSpacing(12); 

    QGroupBox* groupSolverMesh = new QGroupBox("计算网格绑定");
    QFormLayout* formSolverMesh = new QFormLayout(groupSolverMesh);
    m_comboSolverMesh = new QComboBox();
    formSolverMesh->addRow("参与求解的网格:", m_comboSolverMesh);
    scrollLayout->addWidget(groupSolverMesh);

    QGroupBox* groupPhysics = new QGroupBox("物理场全局设置");
    QFormLayout* formPhysics = new QFormLayout(groupPhysics);
    m_spinE0 = new QDoubleSpinBox(); m_spinE0->setRange(0, 1e9); m_spinE0->setValue(600000.0); m_spinE0->setDecimals(6);
    m_spinRho = new QDoubleSpinBox(); m_spinRho->setRange(0, 1e9); m_spinRho->setValue(10000.0); m_spinRho->setDecimals(6);
    m_spinK = new QDoubleSpinBox(); m_spinK->setRange(0, 1e4); m_spinK->setValue(1.0); m_spinK->setDecimals(6);
    
    QWidget* windWidget = new QWidget();
    QVBoxLayout* windLayout = new QVBoxLayout(windWidget);
    windLayout->setContentsMargins(0, 0, 0, 0);
    m_spinWindX = new QDoubleSpinBox(); m_spinWindX->setRange(-100, 100); m_spinWindX->setDecimals(6); m_spinWindX->setPrefix("X: ");
    m_spinWindY = new QDoubleSpinBox(); m_spinWindY->setRange(-100, 100); m_spinWindY->setDecimals(6); m_spinWindY->setPrefix("Y: ");
    m_spinWindZ = new QDoubleSpinBox(); m_spinWindZ->setRange(-100, 100); m_spinWindZ->setDecimals(6); m_spinWindZ->setPrefix("Z: ");
    windLayout->addWidget(m_spinWindX); windLayout->addWidget(m_spinWindY); windLayout->addWidget(m_spinWindZ);

    formPhysics->addRow("起晕场强 E₀ (V/m):", m_spinE0);
    formPhysics->addRow("表面电荷初值 (C/m^3/ε₀):", m_spinRho);
    formPhysics->addRow("离子迁移率 K:", m_spinK);
    formPhysics->addRow("环境风速 (m/s):", windWidget);
    scrollLayout->addWidget(groupPhysics);


    QGroupBox* groupSolver = new QGroupBox("求解器数值控制");
    QFormLayout* formSolver = new QFormLayout(groupSolver);
    m_spinWRho = new QDoubleSpinBox(); m_spinWRho->setRange(0.01, 1.0); m_spinWRho->setDecimals(6); m_spinWRho->setValue(1.0);
    m_spinGoalConv = new QDoubleSpinBox(); m_spinGoalConv->setRange(0.1, 1.0); m_spinGoalConv->setDecimals(6); m_spinGoalConv->setValue(0.95);
    m_spinTolE = new QDoubleSpinBox(); m_spinTolE->setRange(0, 1); m_spinTolE->setDecimals(6); m_spinTolE->setValue(0.01);
    m_spinTolRho = new QDoubleSpinBox(); m_spinTolRho->setRange(0, 1); m_spinTolRho->setDecimals(6); m_spinTolRho->setValue(0.01);
    QSpinBox* spinIter = new QSpinBox(); spinIter->setRange(1, 10000); spinIter->setValue(100);
    formSolver->addRow("电荷松弛因子 (w_rho):", m_spinWRho);
    formSolver->addRow("目标收敛率 (Goal Conv):", m_spinGoalConv);
    formSolver->addRow("场强收敛容差 (Tol E):", m_spinTolE);
    formSolver->addRow("电荷收敛容差 (Tol Rho):", m_spinTolRho);
    formSolver->addRow("最大更新次数:", spinIter);
    scrollLayout->addWidget(groupSolver);

    QGroupBox* groupHPC = new QGroupBox("并行求解设置");
    QVBoxLayout* layHPC = new QVBoxLayout(groupHPC);
    QFormLayout* formHPC = new QFormLayout();

    m_comboSolverType = new QComboBox();
    m_comboSolverType->addItem("Takuma-电位排序法 CPU", "CPU");
    m_comboSolverType->addItem("Tabata CPU", "TABATA_CPU");
    m_comboSolverType->addItem("Tabata GPU (CUDA)", "GPU");
    m_comboSolverType->addItem("Takuma-双栈法 CPU", "DOUBLESTACK");
    formHPC->addRow("求解器类型:", m_comboSolverType);

    m_spinCores = new QSpinBox();
    int logicalCores = std::thread::hardware_concurrency();
    int physCores = std::max(1, logicalCores / 2);  // 物理核心数 ≈ 逻辑线程数/2
    m_spinCores->setRange(1, logicalCores); m_spinCores->setValue(physCores);
    formHPC->addRow(QString("MPI 并行进程数:"), m_spinCores);
    layHPC->addLayout(formHPC);

    QPushButton* btnCalculate = new QPushButton("生成 JSON 并启动求解器");
    btnCalculate->setObjectName("btnCalculate"); 
    btnCalculate->setMinimumHeight(45);
    layHPC->addWidget(btnCalculate);
    
    QPushButton* btnExportData = new QPushButton("导出计算结果数据");
    btnExportData->setMinimumHeight(35);
    layHPC->addWidget(btnExportData);

    scrollLayout->addWidget(groupHPC);
    scrollLayout->addStretch(); 

    scrollContent->setLayout(scrollLayout);
    scrollArea->setWidget(scrollContent);
    mainLayout->addWidget(scrollArea);

    ui->stackedWidgetSettings->addWidget(m_pageSolver);

    connect(btnExportData, &QPushButton::clicked, this, [=](){
        int index = m_comboSolverMesh->currentIndex();
        if (index < 0 || m_comboSolverMesh->count() == 0 || !m_comboSolverMesh->itemData(0).isValid()) {
            QMessageBox::warning(this, "警告", "没有选择有效的计算网格！"); return;
        }
        QTreeWidgetItem* selectedMeshNode = static_cast<QTreeWidgetItem*>(m_comboSolverMesh->itemData(index).value<void*>());
        if (!m_meshFiles.count(selectedMeshNode)) {
            QMessageBox::warning(this, "警告", "选中的网格文件丢失或无效！"); return;
        }
        QString activeMeshFile = m_meshFiles[selectedMeshNode];

        QString workingDir = QApplication::applicationDirPath();

        QString outputFolder = workingDir + "/results/";
        QFileInfo rhoInfo(outputFolder + "/rho.txt");
        
        if (!rhoInfo.exists()) {
            QMessageBox::warning(this, "提示", "尚未找到计算结果文件，请确认求解是否已成功完成。");
            return;
        }

        QString exportPath = QFileDialog::getExistingDirectory(this, "选择结果导出目录", workingDir);
        if (!exportPath.isEmpty()) {
            QStringList filesToCopy = {"rho.txt", "E.txt", "E0.txt"};
            bool allSuccess = true;
            QString errorDetails;

            for (const QString& fileName : filesToCopy) {
                QString sourceFile = QDir::cleanPath(outputFolder + "/" + fileName);
                QString targetFile = QDir::cleanPath(exportPath + "/" + fileName);
                QFileInfo srcInfo(sourceFile);
                QFileInfo tgtInfo(targetFile);

                if (srcInfo.absoluteFilePath() == tgtInfo.absoluteFilePath()) continue;

                if (!srcInfo.exists()) {
                    allSuccess = false;
                    errorDetails += QString("- 源文件丢失: %1\n").arg(fileName);
                    continue;
                }

                if (tgtInfo.exists()) {
                    if (!QFile::remove(tgtInfo.absoluteFilePath())) {
                        allSuccess = false;
                        errorDetails += QString("- 无法覆盖: %1\n").arg(fileName);
                        continue;
                    }
                }

                if (!QFile::copy(srcInfo.absoluteFilePath(), tgtInfo.absoluteFilePath())) {
                    allSuccess = false;
                    errorDetails += QString("- 无写入权限或磁盘空间: %1\n").arg(fileName);
                }
            }

            if (allSuccess) {
                QMessageBox::information(this, "导出成功", "所有计算结果已成功导出至:\n" + exportPath);
                logMessage("数据已成功导出至: " + exportPath, "success");
            } else {
                QMessageBox::warning(this, "部分失败", "结果文件导出存在异常！\n\n详细原因：\n" + errorDetails);
                logMessage("导出发生异常：\n" + errorDetails, "error");
            }
        }
    });

    connect(btnCalculate, &QPushButton::clicked, this, [=](){
        m_config.physics.E_onset = m_spinE0->value();
        m_config.physics.rho_surface = m_spinRho->value();
        m_config.physics.K_mobility = m_spinK->value();
        m_config.physics.wind_x = m_spinWindX->value();
        m_config.physics.wind_y = m_spinWindY->value();
        m_config.physics.wind_z = m_spinWindZ->value();
        m_config.solver.w_rho = m_spinWRho->value();
        m_config.solver.goal_convergence = m_spinGoalConv->value();
        m_config.solver.tolerance_E = m_spinTolE->value();
        m_config.solver.tolerance_rho = m_spinTolRho->value();
        m_config.solver.max_update_times = spinIter->value();
        m_config.solver.order = 1;
        m_config.solver.num_cores = m_spinCores->value();
        if (m_comboSolverType)
            m_config.solver.solver_type = m_comboSolverType->currentData().toString().toStdString();
        generateJsonAndCalculate();
    });
}

// ============ 参数系统 ============

void MainWindow::setupParamUI() {
    // 初始化默认页
    if (m_config.paramPages.empty()) {
        m_config.paramPages.push_back({"模型参数", {}});
        m_config.paramPages.push_back({"后处理参数", {}});
    }

    m_pageParams = new QWidget();
    auto* mainLayout = new QVBoxLayout(m_pageParams);

    // 工具栏
    auto* toolbar = new QHBoxLayout();
    auto* btnAddPage = new QPushButton("添加页面");
    auto* btnDelPage = new QPushButton("删除页面");
    auto* btnAddParam = new QPushButton("添加参数");
    auto* btnDelParam = new QPushButton("删除参数");
    toolbar->addWidget(btnAddPage);
    toolbar->addWidget(btnDelPage);
    toolbar->addStretch();
    toolbar->addWidget(btnAddParam);
    toolbar->addWidget(btnDelParam);
    mainLayout->addLayout(toolbar);

    // Tab 页
    m_paramTabs = new QTabWidget();
    m_paramTabs->setTabsClosable(false);
    mainLayout->addWidget(m_paramTabs);

    // 状态栏
    m_paramStatus = new QLabel();
    m_paramStatus->setStyleSheet("color:#888; font-size:13px;");
    mainLayout->addWidget(m_paramStatus);

    // 初始化表格
    for (size_t pi = 0; pi < m_config.paramPages.size(); pi++) {
        auto* table = new QTableWidget(0, 5);
        table->setHorizontalHeaderLabels({"名称", "表达式", "值", "单位", "描述"});
        table->horizontalHeader()->setStretchLastSection(true);
        table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
        table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
        table->setSelectionBehavior(QAbstractItemView::SelectRows);
        table->setAlternatingRowColors(true);
        table->verticalHeader()->setVisible(false);

        // 双击修改单元格
        connect(table, &QTableWidget::cellChanged, [this, pi](int row, int col) {
            if (pi >= (int)m_config.paramPages.size()) return;
            auto& entries = m_config.paramPages[pi].entries;
            if (row >= (int)entries.size()) return;
            auto* table = qobject_cast<QTableWidget*>(m_paramTabs->widget(pi));
            if (!table) return;
            m_isUpdatingUI = true;
            if (col == 0) entries[row].name = table->item(row, 0)->text().toStdString();
            else if (col == 1) entries[row].expression = table->item(row, 1)->text().toStdString();
            else if (col == 3) entries[row].unit = table->item(row, 3)->text().toStdString();
            else if (col == 4) entries[row].description = table->item(row, 4)->text().toStdString();
            m_isUpdatingUI = false;
            m_projectModified = true;
            if (col == 1) evalAllParams();  // 表达式变化时重新求值
        });

        m_paramTabs->addTab(table, QString::fromStdString(m_config.paramPages[pi].name));
    }

    // Tab 名称双击修改
    connect(m_paramTabs, &QTabWidget::tabBarDoubleClicked, [this](int idx) {
        bool ok;
        QString newName = QInputDialog::getText(this, "重命名参数页", "页面名称:", QLineEdit::Normal,
                                                m_paramTabs->tabText(idx), &ok);
        if (ok && !newName.isEmpty() && idx < (int)m_config.paramPages.size()) {
            m_config.paramPages[idx].name = newName.toStdString();
            m_paramTabs->setTabText(idx, newName);
        }
    });

    auto rebuildTabs = [this]() {
        m_paramTabs->blockSignals(true);
        while (m_paramTabs->count() > 0) m_paramTabs->removeTab(0);
        for (size_t pi = 0; pi < m_config.paramPages.size(); pi++) {
            auto* table = new QTableWidget(0, 5);
            table->setHorizontalHeaderLabels({"名称", "表达式", "值", "单位", "描述"});
            table->horizontalHeader()->setStretchLastSection(true);
            table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
            table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
            table->setSelectionBehavior(QAbstractItemView::SelectRows);
            table->setAlternatingRowColors(true);
            table->verticalHeader()->setVisible(false);
            connect(table, &QTableWidget::cellChanged, [this, pi](int row, int col) {
                if (pi >= (int)m_config.paramPages.size()) return;
                auto& entries = m_config.paramPages[pi].entries;
                if (row >= (int)entries.size()) return;
                auto* t = qobject_cast<QTableWidget*>(m_paramTabs->widget(pi));
                if (!t) return;
                m_isUpdatingUI = true;
                if (col == 0) entries[row].name = t->item(row, 0)->text().toStdString();
                else if (col == 1) entries[row].expression = t->item(row, 1)->text().toStdString();
                else if (col == 3) entries[row].unit = t->item(row, 3)->text().toStdString();
                else if (col == 4) entries[row].description = t->item(row, 4)->text().toStdString();
                m_isUpdatingUI = false;
                m_projectModified = true;
                if (col == 1) evalAllParams();
            });
            m_paramTabs->addTab(table, QString::fromStdString(m_config.paramPages[pi].name));
        }
        m_paramTabs->blockSignals(false);
        for (size_t pi = 0; pi < m_config.paramPages.size(); pi++)
            refreshParamTable(pi);
        evalAllParams();
    };

    // 添加页面
    connect(btnAddPage, &QPushButton::clicked, [this, rebuildTabs]() {
        bool ok;
        QString name = QInputDialog::getText(this, "新建参数页", "页面名称:", QLineEdit::Normal, "新参数页", &ok);
        if (ok && !name.isEmpty()) {
            m_config.paramPages.push_back({name.toStdString(), {}});
            m_projectModified = true;
            rebuildTabs();
        }
    });

    // 删除页面
    connect(btnDelPage, &QPushButton::clicked, [this, rebuildTabs]() {
        int idx = m_paramTabs->currentIndex();
        if (idx >= 0 && idx < (int)m_config.paramPages.size()) {
            m_config.paramPages.erase(m_config.paramPages.begin() + idx);
            m_projectModified = true;
            rebuildTabs();
        }
    });

    // 添加参数
    connect(btnAddParam, &QPushButton::clicked, [this, rebuildTabs]() {
        int idx = m_paramTabs->currentIndex();
        if (idx < 0 || idx >= (int)m_config.paramPages.size()) return;
        m_config.paramPages[idx].entries.push_back({"param", "0", 0.0, "", ""});
        m_projectModified = true;
        refreshParamTable(idx);
    });

    // 删除参数
    connect(btnDelParam, &QPushButton::clicked, [this, rebuildTabs]() {
        int idx = m_paramTabs->currentIndex();
        if (idx < 0 || idx >= (int)m_config.paramPages.size()) return;
        auto* table = qobject_cast<QTableWidget*>(m_paramTabs->currentWidget());
        if (!table) return;
        int row = table->currentRow();
        if (row >= 0 && row < (int)m_config.paramPages[idx].entries.size()) {
            m_config.paramPages[idx].entries.erase(m_config.paramPages[idx].entries.begin() + row);
            m_projectModified = true;
            refreshParamTable(idx);
        }
    });

    ui->stackedWidgetSettings->addWidget(m_pageParams);

    // 加载默认数据到表格
    for (size_t pi = 0; pi < m_config.paramPages.size(); pi++)
        refreshParamTable(pi);
    evalAllParams();
}

void MainWindow::evalAllParams() {
    // 构建参数名→值映射, 初始用存储的值
    std::map<std::string, double> vars;
    for (auto& pg : m_config.paramPages)
        for (auto& e : pg.entries)
            vars[e.name] = e.value;

    // 拓扑排序: Kahn算法
    std::map<std::string, std::vector<std::string>> deps; // name → 依赖列表
    std::map<std::string, int> inDegree;
    for (auto& pg : m_config.paramPages)
        for (auto& e : pg.entries) {
            auto d = ParamParser::dependencies(e.expression);
            deps[e.name] = std::vector<std::string>(d.begin(), d.end());
            if (inDegree.find(e.name) == inDegree.end()) inDegree[e.name] = 0;
            for (auto& dep : d) inDegree[e.name]++;
            for (auto& dep : d) {
                if (inDegree.find(dep) == inDegree.end()) inDegree[dep] = 0;
            }
        }

    std::vector<std::string> queue;
    for (auto& [name, deg] : inDegree)
        if (deg == 0) queue.push_back(name);

    std::vector<std::string> sorted;
    while (!queue.empty()) {
        std::string name = queue.back(); queue.pop_back();
        sorted.push_back(name);
        for (auto& [n, dlist] : deps) {
            if (std::find(dlist.begin(), dlist.end(), name) != dlist.end()) {
                inDegree[n]--;
                if (inDegree[n] == 0) queue.push_back(n);
            }
        }
    }

    // 按拓扑序求值
    for (auto& name : sorted) {
        for (auto& pg : m_config.paramPages)
            for (auto& e : pg.entries) {
                if (e.name == name) {
                    bool ok;
                    double v = ParamParser::evaluate(e.expression, vars, e.unit, &ok);
                    if (ok) { e.value = v; vars[name] = v; }
                }
            }
    }

    // 检测循环依赖
    bool hasCycle = false;
    std::string cycleNames;
    for (auto& [name, deg] : inDegree)
        if (deg > 0) { hasCycle = true; cycleNames += name + " "; }

    if (hasCycle && m_paramStatus)
        m_paramStatus->setText(QString("循环依赖: %1").arg(QString::fromStdString(cycleNames)));
    else if (m_paramStatus)
        m_paramStatus->setText("参数求值完成");

    // 刷新所有表格的值列
    for (size_t pi = 0; pi < m_config.paramPages.size(); pi++) {
        auto* table = qobject_cast<QTableWidget*>(m_paramTabs->widget(pi));
        if (!table) continue;
        auto& entries = m_config.paramPages[pi].entries;
        for (size_t r = 0; r < entries.size(); r++) {
            auto* item = table->item(r, 2);
            if (item) item->setText(QString::number(entries[r].value, 'g', 6));
        }
    }
    pushParamsToBindings();
}

void MainWindow::refreshParamTable(int pageIdx) {
    if (pageIdx < 0 || pageIdx >= (int)m_config.paramPages.size()) return;
    auto* table = qobject_cast<QTableWidget*>(m_paramTabs->widget(pageIdx));
    if (!table) return;
    auto& entries = m_config.paramPages[pageIdx].entries;

    table->blockSignals(true);
    table->setRowCount(entries.size());
    for (size_t r = 0; r < entries.size(); r++) {
        auto& e = entries[r];
        table->setItem(r, 0, new QTableWidgetItem(QString::fromStdString(e.name)));
        table->setItem(r, 1, new QTableWidgetItem(QString::fromStdString(e.expression)));
        auto* valItem = new QTableWidgetItem(QString::number(e.value, 'g', 6));
        valItem->setFlags(valItem->flags() & ~Qt::ItemIsEditable);
        valItem->setBackground(QColor(240, 240, 240));
        table->setItem(r, 2, valItem);
        auto* unitCombo = new QComboBox();
        unitCombo->addItems({"", "m", "cm", "mm", "km", "deg", "rad", "V", "kV", "MV"});
        unitCombo->setCurrentText(QString::fromStdString(e.unit));
        connect(unitCombo, &QComboBox::currentTextChanged, [this, pageIdx, r](const QString& txt) {
            if (m_isUpdatingUI) return;
            if (pageIdx < (int)m_config.paramPages.size() && r < (int)m_config.paramPages[pageIdx].entries.size())
                m_config.paramPages[pageIdx].entries[r].unit = txt.toStdString();
            m_projectModified = true;
            pushParamsToBindings();
        });
        table->setCellWidget(r, 3, unitCombo);
        table->setItem(r, 4, new QTableWidgetItem(QString::fromStdString(e.description)));
    }
    table->blockSignals(false);
}

void MainWindow::installBindButtons() {
    auto addBtnToSpin = [this](QDoubleSpinBox* spin) {
        if (!spin) return;
        auto* parent = spin->parentWidget();
        if (!parent) return;
        auto makeBtn = [this, spin]() {
            auto* btn = new QToolButton();
            btn->setText("="); btn->setMaximumSize(22, 22);
            btn->setToolTip("绑定/解绑参数");
            btn->setStyleSheet("QToolButton { border:1px solid #CCC; background:#EEE; font-weight:bold; }");
            m_bindBtns[btn] = spin;
            connect(btn, &QToolButton::clicked, [this, spin, btn]() { showBindDialog(spin, btn); });
            if (m_paramBindings.count(spin)) {
                spin->setStyleSheet("QDoubleSpinBox { background-color: #E3F2FD; }");
                spin->setReadOnly(true);
                btn->setStyleSheet("QToolButton { border:1px solid #2196F3; background:#BBDEFB; font-weight:bold; color:#1565C0; }");
            }
            return btn;
        };
        // 递归搜索包含 target 的 QFormLayout (支持嵌套在 QVBoxLayout 等内)
        std::function<QFormLayout*(QLayout*)> findForm = [&](QLayout* lay) -> QFormLayout* {
            if (!lay) return nullptr;
            // 先检查 lay 本身是否就是包含 spin 的 QFormLayout
            if (auto* self = dynamic_cast<QFormLayout*>(lay)) {
                for (int r = 0; r < self->rowCount(); r++)
                    if (self->itemAt(r, QFormLayout::FieldRole) &&
                        self->itemAt(r, QFormLayout::FieldRole)->widget() == spin)
                        return self;
            }
            for (int i = 0; i < lay->count(); i++) {
                auto* item = lay->itemAt(i);
                if (auto* form = dynamic_cast<QFormLayout*>(item->layout())) {
                    for (int r = 0; r < form->rowCount(); r++)
                        if (form->itemAt(r, QFormLayout::FieldRole) &&
                            form->itemAt(r, QFormLayout::FieldRole)->widget() == spin)
                            return form;
                }
                if (item->layout()) {
                    auto* found = findForm(item->layout());
                    if (found) return found;
                }
                if (item->widget() && item->widget()->layout()) {
                    auto* found = findForm(item->widget()->layout());
                    if (found) return found;
                }
            }
            return nullptr;
        };
        // 策略1: QFormLayout → wrapper替换 (递归搜索嵌套layout)
        auto* form = findForm(parent->layout());
        if (form) {
            for (int r = 0; r < form->rowCount(); r++) {
                auto* item = form->itemAt(r, QFormLayout::FieldRole);
                if (item && item->widget() == spin) {
                    auto* wrap = new QWidget();
                    auto* hlay = new QHBoxLayout(wrap);
                    hlay->setContentsMargins(0,0,0,0); hlay->setSpacing(1);
                    spin->setParent(wrap); hlay->addWidget(spin);
                    hlay->addWidget(makeBtn());
                    QLayoutItem* labelItem = form->itemAt(r, QFormLayout::LabelRole);
                    QString labelText;
                    if (labelItem && labelItem->widget())
                        if (auto* lbl = qobject_cast<QLabel*>(labelItem->widget())) labelText = lbl->text();
                    form->removeRow(r);
                    form->insertRow(r, new QLabel(labelText), wrap);
                    return;
                }
            }
        }
        // 策略2: QBoxLayout(QH/QV) → 直接插入按钮
        auto* box = qobject_cast<QBoxLayout*>(parent->layout());
        if (!box) {
            // 如果顶层不是box, 查第一个子layout是不是box
            auto* lay = parent->layout();
            if (lay && lay->count() > 0)
                for (int i = 0; i < lay->count(); i++)
                    if (auto* b = qobject_cast<QBoxLayout*>(lay->itemAt(i)->layout())) { box = b; break; }
        }
        if (box) {
            for (int i = 0; i < box->count(); i++)
                if (box->itemAt(i)->widget() == spin) {
                    box->insertWidget(i + 1, makeBtn());
                    return;
                }
        }
        // 策略3: 其它布局 → 作为sibling放置
        auto* btn = makeBtn();
        btn->setParent(parent);
        QPoint pos = spin->pos();
        btn->move(pos.x() + spin->width() + 2, pos.y());
        btn->show();
    };

    std::vector<QDoubleSpinBox*> spins = {
        m_geomW, m_geomH, m_geomD, m_geomE,
        m_geomPX, m_geomPY, m_geomPZ, m_geomRotAng,
        m_spinFilletR, m_spinChamfD,
        m_spArrSpc, m_spCirCx, m_spCirCy, m_spCirCz, m_spCirAng,
        m_spinE0, m_spinRho, m_spinK,
        m_spinWindX, m_spinWindY, m_spinWindZ,
        m_spinWRho, m_spinGoalConv, m_spinTolE, m_spinTolRho,
        m_spinBdrVoltage
    };
    for (auto* spin : spins) addBtnToSpin(spin);
    for (auto* spin : m_xfmSpins) addBtnToSpin(spin);
    // 注册spinbox稳定名称(持久化绑定用)
    m_spinboxNames.clear();
    m_spinboxNames[m_geomW]="geomW"; m_spinboxNames[m_geomH]="geomH"; m_spinboxNames[m_geomD]="geomD"; m_spinboxNames[m_geomE]="geomE";
    m_spinboxNames[m_geomPX]="geomPX"; m_spinboxNames[m_geomPY]="geomPY"; m_spinboxNames[m_geomPZ]="geomPZ"; m_spinboxNames[m_geomRotAng]="geomRotAng";
    m_spinboxNames[m_spinFilletR]="filletR"; m_spinboxNames[m_spinChamfD]="chamfD";
    m_spinboxNames[m_spArrSpc]="arrSpc"; m_spinboxNames[m_spCirCx]="cirCx"; m_spinboxNames[m_spCirCy]="cirCy"; m_spinboxNames[m_spCirCz]="cirCz"; m_spinboxNames[m_spCirAng]="cirAng";
    m_spinboxNames[m_spinE0]="E0"; m_spinboxNames[m_spinRho]="rho"; m_spinboxNames[m_spinK]="K";
    m_spinboxNames[m_spinWindX]="windX"; m_spinboxNames[m_spinWindY]="windY"; m_spinboxNames[m_spinWindZ]="windZ";
    m_spinboxNames[m_spinWRho]="wRho"; m_spinboxNames[m_spinGoalConv]="goalConv"; m_spinboxNames[m_spinTolE]="tolE"; m_spinboxNames[m_spinTolRho]="tolRho";
    m_spinboxNames[m_spinBdrVoltage]="bdrVoltage";
}

void MainWindow::showBindDialog(QDoubleSpinBox* spin, QToolButton* btn) {
    QDialog dlg(this);
    dlg.setWindowTitle("参数表达式");
    auto* lay = new QVBoxLayout(&dlg);
    // 参数提示
    QStringList paramList;
    for (auto& pg : m_config.paramPages)
        for (auto& e : pg.entries)
            paramList << QString("%1(%2)").arg(QString::fromStdString(e.name)).arg(e.value);
    auto* hint = new QLabel("可用参数: " + paramList.join(", "));
    hint->setStyleSheet("color:#888; font-size:13px;");
    hint->setWordWrap(true);
    lay->addWidget(hint);
    // 表达式输入
    auto* edit = new QLineEdit();
    edit->setPlaceholderText("例: a 或 a + b/2 (清空=解除绑定)");
    if (m_paramBindings.count(spin))
        edit->setText(QString::fromStdString(m_paramBindings[spin]));
    lay->addWidget(edit);
    // 预览
    auto* preview = new QLabel();
    preview->setStyleSheet("color:#1565C0; font-weight:bold;");
    lay->addWidget(preview);
    auto updatePreview = [&]() {
        QString expr = edit->text().trimmed();
        if (expr.isEmpty()) { preview->setText(""); return; }
        std::map<std::string, double> vars;
        for (auto& pg : m_config.paramPages)
            for (auto& e : pg.entries) vars[e.name] = e.value;
        bool ok;
        double v = ParamParser::evaluate(expr.toStdString(), vars, "", &ok);
        preview->setText(ok ? QString("→ 当前值: %1").arg(v, 0, 'g', 6) : "→ 表达式无效");
    };
    connect(edit, &QLineEdit::textChanged, [&](const QString&){ updatePreview(); });
    updatePreview();
    auto* btnBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    lay->addWidget(btnBox);
    connect(btnBox, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(btnBox, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
    if (dlg.exec() == QDialog::Accepted) {
        QString expr = edit->text().trimmed();
        if (expr.isEmpty()) {
            m_paramBindings.erase(spin);
            spin->setStyleSheet("");
            spin->setReadOnly(false);
            btn->setStyleSheet("QToolButton { border:1px solid #CCC; background:#EEE; font-weight:bold; }");
        } else {
            m_paramBindings[spin] = expr.toStdString();
            spin->setStyleSheet("QDoubleSpinBox { background-color: #E3F2FD; }");
            spin->setReadOnly(true);
            btn->setStyleSheet("QToolButton { border:1px solid #2196F3; background:#BBDEFB; font-weight:bold; color:#1565C0; }");
            pushParamsToBindings();
        }
    }
}

void MainWindow::pushParamsToBindings() {
    // 获取 spinbox 的目标基础单位 (几何面板始终用基础单位: 长度=m, 角度=rad, 电压=V)
    auto getTargetUnit = [this](QDoubleSpinBox* spin) -> std::string {
        if (spin == m_spinBdrVoltage) return "V";
        if (spin == m_geomRotAng || spin == m_spCirAng) return "rad";
        if (spin == m_geomW || spin == m_geomH || spin == m_geomD || spin == m_geomE ||
            spin == m_geomPX || spin == m_geomPY || spin == m_geomPZ ||
            spin == m_spinFilletR || spin == m_spinChamfD || spin == m_spArrSpc ||
            spin == m_spCirCx || spin == m_spCirCy || spin == m_spCirCz) return "m";
        return "";
    };
    // 参数名 → unit 映射
    std::map<std::string, std::string> paramUnit;
    for (auto& pg : m_config.paramPages)
        for (auto& e : pg.entries) if (!e.unit.empty()) paramUnit[e.name] = e.unit;

    std::map<std::string, double> vars;
    for (auto& pg : m_config.paramPages)
        for (auto& e : pg.entries) vars[e.name] = e.value;
    for (auto& [spin, expr] : m_paramBindings) {
        if (expr.empty()) continue;
        bool ok;
        double val = ParamParser::evaluate(expr, vars, "", &ok);
        if (!ok) continue;
        // 单位换算: 参数单位 → spinbox 目标单位
        auto itParam = paramUnit.find(expr);
        std::string targetUnit = getTargetUnit(spin);
        if (itParam != paramUnit.end() && !targetUnit.empty())
            val = ParamParser::convertUnit(val, itParam->second, targetUnit);
        // 放宽范围, 允许计算后的小值
        double oldMin = spin->minimum(), oldMax = spin->maximum();
        if (val < oldMin) spin->setMinimum(val);
        if (val > oldMax) spin->setMaximum(val);
        spin->setValue(val);
        spin->setMinimum(oldMin);
        spin->setMaximum(oldMax);
    }
}

void MainWindow::refreshSolverMeshCombo() {
    m_comboSolverMesh->clear();
    for(int i = 0; i < m_rootMeshNode->childCount(); ++i) {
        QTreeWidgetItem* child = m_rootMeshNode->child(i);
        if (m_meshFiles.count(child) && !m_meshFiles[child].isEmpty()) {
            m_comboSolverMesh->addItem(child->text(0), QVariant::fromValue(static_cast<void*>(child)));
        }
    }
    
    if (m_comboSolverMesh->count() == 0) {
        m_comboSolverMesh->addItem("暂无可用网格，请先生成或导入");
        QStandardItemModel* model = qobject_cast<QStandardItemModel*>(m_comboSolverMesh->model());
        if (model && model->item(0)) {
            model->item(0)->setEnabled(false);
        }
    }
}

void MainWindow::refreshBoundaryUI() {
    m_isUpdatingUI = true;
    m_comboBoundaryGroup->clear();

    for (int i = 0; i < m_rootGroupNode->childCount(); ++i) {
        QTreeWidgetItem* groupNode = m_rootGroupNode->child(i);
        QString groupName = groupNode->text(0);
        
        m_comboBoundaryGroup->addItem(groupName, QVariant::fromValue(static_cast<void*>(groupNode)));

        if (m_boundaryConfigs.find(groupNode) == m_boundaryConfigs.end()) {
            BoundaryParams params;
            params.apply = false;
            params.voltage = 0.0;
            params.is_corona = false;
            m_boundaryConfigs[groupNode] = params;
        }
    }
    m_isUpdatingUI = false;
    if (m_comboBoundaryGroup->count() > 0) {
        onBoundaryGroupTargetChanged(0);
    }
}

void MainWindow::onBoundaryGroupTargetChanged(int index) {
    if (m_isUpdatingUI || index < 0) return;
    QTreeWidgetItem* targetGroup = static_cast<QTreeWidgetItem*>(m_comboBoundaryGroup->itemData(index).value<void*>());
    if (!targetGroup || m_boundaryConfigs.find(targetGroup) == m_boundaryConfigs.end()) return;

    m_isUpdatingUI = true;
    const BoundaryParams& params = m_boundaryConfigs[targetGroup];
    m_chkApplyBoundary->setChecked(params.apply);
    m_spinBdrVoltage->setValue(params.voltage);
    m_cmbBdrMode->setCurrentIndex(params.useFunction ? 1 : 0);
    m_editBdrExpr->setText(QString::fromStdString(params.expression));
    m_chkBdrCorona->setChecked(params.is_corona);

    m_spinBdrVoltage->setEnabled(params.apply);
    m_chkBdrCorona->setEnabled(params.apply);
    m_isUpdatingUI = false;
}

void MainWindow::onBoundaryParamsChanged() {
    if (m_isUpdatingUI) return;
    int index = m_comboBoundaryGroup->currentIndex();
    if (index < 0) return;
    QTreeWidgetItem* targetGroup = static_cast<QTreeWidgetItem*>(m_comboBoundaryGroup->itemData(index).value<void*>());
    if (targetGroup) {
        m_projectModified = true;
        BoundaryParams& params = m_boundaryConfigs[targetGroup];
        params.apply = m_chkApplyBoundary->isChecked();
        params.voltage = m_spinBdrVoltage->value();
        params.useFunction = (m_cmbBdrMode->currentIndex() == 1);
        params.expression = m_editBdrExpr->text().toStdString();
        params.is_corona = m_chkBdrCorona->isChecked();
        
        m_spinBdrVoltage->setEnabled(params.apply);
        m_chkBdrCorona->setEnabled(params.apply);
    }
}

// 对导入网格: 基于 m_groupData 同步物理组, 写出统一缓存路径 .msh + .pgroups.json
// 返回同步后的 .msh 路径; 失败返回空串
QString MainWindow::syncImportedMeshToCache() {
    if (!m_activeMeshNode || !m_meshFiles.count(m_activeMeshNode)) return QString();
    QString srcMsh = m_meshFiles[m_activeMeshNode];

    QString cacheDir = QApplication::applicationDirPath() + "/.cache";
    QDir().mkpath(cacheDir);
    QString outMsh = cacheDir + "/imported_synced.msh";
    QString outPg = outMsh + ".pgroups.json";

    try {
        MeshEngine::ensureInit();
        gmsh::clear();
        gmsh::open(srcMsh.toStdString());
        fprintf(stderr, "[DEBUG SYNC] opened %s\n", srcMsh.toStdString().c_str());

        // Step A: 基于 m_groupData 生成 .pgroups.json
        std::ofstream pgf(outPg.toStdString());
        pgf << "{\n  \"mesh_file\": \"" << outMsh.toStdString() << "\",\n  \"groups\": [\n";
        int groupIdx = 0, tagCounter = 1;
        for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
            auto* item = m_rootGroupNode->child(i);
            if (m_groupData[item].empty()) continue;
            int tag = m_groupTags.count(item) ? m_groupTags[item] : tagCounter;
            std::set<int> entityTags;
            for (int fid : m_groupData[item])
                if (m_faceToEntity.count(fid)) entityTags.insert(m_faceToEntity[fid]);
            if (entityTags.empty()) continue;
            std::set<std::size_t> allElems, allNodes;
            for (int et : entityTags) {
                std::vector<int> eTypes;
                std::vector<std::vector<std::size_t>> eTags, nTags;
                gmsh::model::mesh::getElements(eTypes, eTags, nTags, 2, et);
                for (size_t ti = 0; ti < eTags.size(); ti++) {
                    allElems.insert(eTags[ti].begin(), eTags[ti].end());
                    allNodes.insert(nTags[ti].begin(), nTags[ti].end());
                }
            }
            if (allNodes.empty()) continue;
            if (groupIdx > 0) pgf << ",\n";
            pgf << "    {\"name\":\"" << item->text(0).toStdString()
                << "\", \"dim\":2, \"tag\":" << tag;
            pgf << ", \"entities\":[";
            size_t ei = 0; for (int et : entityTags) { if (ei++) pgf << ","; pgf << et; }
            pgf << "], \"elements\":[";
            ei = 0; for (auto e : allElems) { if (ei++) pgf << ","; pgf << e; }
            pgf << "], \"nodes\":[";
            ei = 0; for (auto n : allNodes) { if (ei++) pgf << ","; pgf << n; }
            pgf << "]}";
            m_groupTags[item] = tag;
            fprintf(stderr, "[DEBUG SYNC] group '%s' tag=%d entities=[",
                    item->text(0).toStdString().c_str(), tag);
            for (int et : entityTags) fprintf(stderr, "%d,", et);
            fprintf(stderr, "] nodes=%zu elems=%zu\n", allNodes.size(), allElems.size());
            groupIdx++; tagCounter++;
        }
        pgf << "\n  ]\n}\n";
        pgf.close();
        fprintf(stderr, "[DEBUG SYNC] .pgroups.json written: %d groups\n", groupIdx);

        // Step B: 重建 PG → 写缓存 .msh (修复 element attribute, 主 BC 依赖 .msh 属性)
        std::vector<std::pair<int,int>> oldPGs;
        gmsh::model::getPhysicalGroups(oldPGs);
        for (auto& pg : oldPGs) gmsh::model::removePhysicalGroups({{pg.first, pg.second}});
        for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
            auto* item = m_rootGroupNode->child(i);
            if (!m_groupTags.count(item)) continue;
            int tag = m_groupTags[item];
            std::set<int> ets;
            for (int fid : m_groupData[item])
                if (m_faceToEntity.count(fid)) ets.insert(m_faceToEntity[fid]);
            if (ets.empty()) continue;
            std::vector<int> entVec(ets.begin(), ets.end());
            gmsh::model::addPhysicalGroup(2, entVec, tag);
            gmsh::model::setPhysicalName(2, tag, item->text(0).toStdString());
        }
        // 添加 3D 体积默认组
        std::vector<std::pair<int,int>> volEnts;
        gmsh::model::getEntities(volEnts, 3);
        std::vector<int> volTags;
        for (auto& v : volEnts) volTags.push_back(v.second);
        if (!volTags.empty()) {
            gmsh::model::addPhysicalGroup(3, volTags, 9999);
            gmsh::model::setPhysicalName(3, 9999, "DefaultVolume");
        }
        QFile::remove(outMsh);
        gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
        gmsh::write(outMsh.toStdString());
        gmsh::clear();
        fprintf(stderr, "[DEBUG SYNC] synced mesh written: %s\n", outMsh.toStdString().c_str());

        return QFile::exists(outMsh) ? outMsh : QString();
    } catch (...) {
        try { gmsh::clear(); } catch (...) {}
        return QString();
    }
}

void MainWindow::generateJsonAndCalculate() {
    int index = m_comboSolverMesh->currentIndex();
    if (index < 0 || m_comboSolverMesh->count() == 0 || !m_comboSolverMesh->itemData(0).isValid()) {
        QMessageBox::warning(this, "警告", "没有选择有效的计算网格！"); return;
    }

    QTreeWidgetItem* selectedMeshNode = static_cast<QTreeWidgetItem*>(m_comboSolverMesh->itemData(index).value<void*>());
    if (!m_meshFiles.count(selectedMeshNode)) {
        QMessageBox::warning(this, "警告", "选中的网格文件丢失或无效！"); return;
    }
    QString activeMeshFile = m_meshFiles[selectedMeshNode];

    QString workingDir = QApplication::applicationDirPath();
    QString configPath = workingDir + "/config.json";
    QString outputFolder = workingDir + "/results/";
    QDir().mkpath(outputFolder);

    json rootObj;
    rootObj["mesh_file"] = activeMeshFile.toStdString();

    // Step 1: 对导入网格, 基于 m_groupData 生成 .pgroups.json (放在 boundaries 之前, 确保 tag 一致)
    bool isImported = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode)
                      && m_isImportedMesh[m_activeMeshNode];
    fprintf(stderr, "[DEBUG SYNC] isImported=%d faceToEntity=%zu groups=%d\n",
            (int)isImported, m_faceToEntity.size(), m_rootGroupNode->childCount());
    fflush(stderr);
    if (isImported && !m_faceToEntity.empty()) {
        QString synced = syncImportedMeshToCache();
        if (!synced.isEmpty()) {
            activeMeshFile = synced;
            m_meshFiles[selectedMeshNode] = activeMeshFile;
            rootObj["mesh_file"] = activeMeshFile.toStdString();
            fprintf(stderr, "[DEBUG SYNC] config mesh_file -> %s\n",
                    activeMeshFile.toStdString().c_str());
        }
    } else {
        try {
            MeshEngine::ensureInit();
            gmsh::clear();
            gmsh::open(activeMeshFile.toStdString());
            MeshEngine::savePhysicalGroupsJson(activeMeshFile.toStdString());
            gmsh::clear();
        } catch (...) { try { gmsh::clear(); } catch (...) {} }
    }

    // Step 2: 生成 config.json (此时 m_groupTags 已在 Step 1 中设置)
    m_config.boundaries.clear();
    for (int i = 0; i < m_rootGroupNode->childCount(); ++i) {
        QTreeWidgetItem* groupNode = m_rootGroupNode->child(i);
        if (m_boundaryConfigs.find(groupNode) == m_boundaryConfigs.end()) continue;
        const BoundaryParams& params = m_boundaryConfigs[groupNode];
        if (!params.apply) continue;
        BoundarySetup bdr;
        bdr.name = groupNode->text(0).toStdString();
        bdr.tag = m_groupTags.count(groupNode) ? m_groupTags[groupNode] : (i + 1);
        bdr.voltage = params.voltage;
        bdr.useFunction = params.useFunction;
        bdr.expression = params.expression;
        bdr.is_corona = params.is_corona;
        m_config.boundaries.push_back(bdr);
        fprintf(stderr, "[DEBUG SYNC] boundary '%s' tag=%d\n", bdr.name.c_str(), bdr.tag);
    }
    fprintf(stderr, "[DEBUG SYNC] total boundaries=%zu\n", m_config.boundaries.size());
    fflush(stderr);
    rootObj["output_folder"] = outputFolder.toStdString();
    rootObj["order"] = m_config.solver.order;
    rootObj["export_txt"] = true;
    rootObj["physics"] = json{
        {"E_onset", m_config.physics.E_onset},
        {"rho_surface", m_config.physics.rho_surface},
        {"K_mobility", m_config.physics.K_mobility},
        {"wind_velocity", {m_config.physics.wind_x, m_config.physics.wind_y, m_config.physics.wind_z}}
    };

    rootObj["solver"] = json{
        {"w_rho", m_config.solver.w_rho},
        {"goal_convergence", m_config.solver.goal_convergence},
        {"max_update_times", m_config.solver.max_update_times},
        {"tolerance_E", m_config.solver.tolerance_E},
        {"tolerance_rho", m_config.solver.tolerance_rho}
    };

    json bdrArray = json::array();
    for (const auto& bdr : m_config.boundaries) {
        json bj;
        bj["tag"] = bdr.tag;
        bj["name"] = bdr.name;
        bj["is_corona"] = bdr.is_corona;
        if (bdr.useFunction && !bdr.expression.empty()) {
            bj["use_function"] = true;
            bj["voltage_expression"] = bdr.expression;
            bj["voltage"] = 0.0;
        } else {
            bj["voltage"] = bdr.voltage;
        }
        bdrArray.push_back(bj);
    }
    rootObj["boundaries"] = bdrArray;

    std::ofstream o(configPath.toStdString());
    o << std::setw(4) << rootObj << std::endl;
    o.close();

    if (m_solverProcess && m_solverProcess->state() == QProcess::Running) {
        QMessageBox::warning(this, "提示", "求解器正在后台运行中，请等待其完成。"); return;
    }

    if (!m_solverProcess) {
        m_solverProcess = new QProcess(this);
        m_solverProcess->setProcessChannelMode(QProcess::MergedChannels); 
        connect(m_solverProcess, &QProcess::readyReadStandardOutput, this, [=](){
            QByteArray out = m_solverProcess->readAllStandardOutput();
            m_textConsole->append(QString("<span style='color:#BBBBBB;'>%1</span>")
                                    .arg(QString::fromLocal8Bit(out).trimmed().toHtmlEscaped().replace("\n", "<br>")));
            m_textConsole->verticalScrollBar()->setValue(m_textConsole->verticalScrollBar()->maximum());
        });

        connect(m_solverProcess, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
                this, [=](int exitCode, QProcess::ExitStatus exitStatus) {
            m_btnCancelCompute->setVisible(false);
            if (exitCode == 0 && exitStatus == QProcess::NormalExit) {
                logMessage("MFEM 并行求解顺利收敛并成功落盘！", "success");
                logMessage("求解成功！请点击设置面板下方的【导出计算结果数据】按钮导出结果。", "highlight");
            } else {
                logMessage("MFEM 求解器异常退出，请检查网格质量或参数。", "error");
                QMessageBox::critical(this, "求解失败", "MFEM 求解器发生异常崩溃！\n请仔细检查控制台输出定位错误。");
            }
        });
    }

    logMessage(">>> 正在启动物理计算引擎...", "highlight");
    m_btnCancelCompute->setVisible(true); m_btnCancelCompute->setEnabled(true);
    m_btnCancelCompute->setText("终止计算");
    QString mpiExec = "mpirun";
    QStringList args;
    QString solverExe;
    if (m_config.solver.solver_type == "GPU") solverExe = "Tabata_GPU";
    else if (m_config.solver.solver_type == "TABATA_CPU") solverExe = "Tabata_CPU";
    else if (m_config.solver.solver_type == "DOUBLESTACK") solverExe = "Takuma_DoubleStack_CPU";
    else solverExe = "Takuma_Potential_CPU"; // 默认 "CPU"
    QString solverPath = QApplication::applicationDirPath() + "/" + solverExe;
    args << "-np" << QString::number(m_config.solver.num_cores) << solverPath << "-c" << configPath;
    m_solverProcess->start(mpiExec, args);
}

void MainWindow::updateInteractionModes() {
    ui->cadViewerWidget->setInteractionModes(ui->btnBoxPick->isChecked(), ui->comboSelectMode->currentIndex() == 1, ui->btnHideMode->isChecked());
}

void MainWindow::importGeometry() {
    QString filter = "CAD Files (*.step *.stp *.iges *.igs *.brep *.rle);;All Files (*)";
    QString filePath = QFileDialog::getOpenFileName(this, "导入几何模型", "", filter);
    if (filePath.isEmpty()) return;

    m_currentStepFile = filePath;
    m_projectModified = true;
    MeshEngine::cleanup();

    m_meshFiles.clear();
    m_isImportedMesh.clear();
    m_meshDataMap.clear();
    for(int i = m_rootMeshNode->childCount() - 1; i >= 0; --i) delete m_rootMeshNode->child(i);

    bool ok = false;
    QString suffix = QFileInfo(filePath).suffix().toLower();
    if (suffix == "step" || suffix == "stp")
        ok = ui->cadViewerWidget->loadStepFile(filePath);
    else if (suffix == "iges" || suffix == "igs")
        ok = ui->cadViewerWidget->loadIgesFile(filePath);
    else if (suffix == "brep" || suffix == "rle")
        ok = ui->cadViewerWidget->loadBrepFile(filePath);
    ui->cadViewerWidget->clearMesh();
    ui->cadViewerWidget->setGeometryTransparent(false);

    if (ok)
        logMessage(QString("已成功导入几何模型: %1").arg(QFileInfo(filePath).fileName()), "info");
    else
        logMessage(QString("导入失败: %1").arg(QFileInfo(filePath).fileName()), "error");
}

void MainWindow::exportGeometry() {
    if (m_currentStepFile.isEmpty()) {
        QMessageBox::warning(this, "提示", "当前没有可导出的几何模型。");
        return;
    }

    // 读取几何体
    TopoDS_Shape shape;
    QString suffix = QFileInfo(m_currentStepFile).suffix().toLower();
    if (suffix == "step" || suffix == "stp") {
        STEPControl_Reader reader;
        if (reader.ReadFile(m_currentStepFile.toUtf8().constData()) != IFSelect_RetDone
            || !reader.TransferRoots()) {
            QMessageBox::warning(this, "错误", "无法读取 STEP 几何文件。"); return;
        }
        shape = reader.OneShape();
    } else if (suffix == "iges" || suffix == "igs") {
        IGESControl_Reader reader;
        if (reader.ReadFile(m_currentStepFile.toUtf8().constData()) != IFSelect_RetDone
            || !reader.TransferRoots()) {
            QMessageBox::warning(this, "错误", "无法读取 IGES 几何文件。"); return;
        }
        shape = reader.OneShape();
    } else {
        // BREP / STL cache
        BRep_Builder B;
        if (!BRepTools::Read(shape, m_currentStepFile.toUtf8().constData(), B)) {
            QMessageBox::warning(this, "错误", "无法读取几何文件。"); return;
        }
    }
    if (shape.IsNull()) { QMessageBox::warning(this, "错误", "几何体为空。"); return; }

    // 文件保存对话框
    QString filter = "STEP (*.step);;IGES (*.igs);;BREP (*.brep)";
    QString selectedFilter;
    QString path = QFileDialog::getSaveFileName(this, "导出几何", "", filter, &selectedFilter);
    if (path.isEmpty()) return;

    QByteArray pathBytes = path.toUtf8();
    bool ok = false;
    if (selectedFilter.contains("STEP")) {
        STEPControl_Writer writer;
        if (writer.Transfer(shape, STEPControl_ManifoldSolidBrep) == IFSelect_RetDone) {
            ok = (writer.Write(pathBytes.constData()) == IFSelect_RetDone);
        }
    } else if (selectedFilter.contains("IGES")) {
        IGESControl_Writer writer;
        ok = writer.AddShape(shape) && writer.Write(pathBytes.constData());
    } else if (selectedFilter.contains("BREP")) {
        ok = BRepTools::Write(shape, pathBytes.constData());
    }
    if (ok)
        logMessage(QString("几何已导出: %1").arg(path), "info");
    else
        QMessageBox::warning(this, "错误", "导出失败。");
}

void MainWindow::importExternalMesh() {
    QString filter = "All Mesh Files (*.msh *.unv *.inp *.vtk *.vtu *.mesh *.stl *.bdf *.diff *.dat);;Gmsh Mesh (*.msh);;All Files (*)";
    QString filePath = QFileDialog::getOpenFileName(this, "导入网格文件", "", filter);
    if (filePath.isEmpty()) return;

    // 非 .msh 格式: 用 Gmsh 打开后写为 .msh v2, 供 MFEM 读取和后续渲染
    QString importPath = filePath;
    QString tempMsh;
    if (!filePath.endsWith(".msh", Qt::CaseInsensitive)) {
        try {
            MeshEngine::ensureInit();
            gmsh::clear();
            gmsh::open(filePath.toStdString());
            tempMsh = QDir::tempPath() + "/imported_mesh.msh";
            gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
            gmsh::write(tempMsh.toStdString());
            gmsh::clear();
            importPath = tempMsh;
        } catch (std::exception& e) {
            try { gmsh::clear(); } catch (...) {}
            QMessageBox::warning(this, "导入失败", QString("Gmsh 无法打开此格式: %1").arg(e.what()));
            return;
        }
    }

    std::vector<ImportedGroup> extractedGroups;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    FaceMeshMap meshData = MeshEngine::importMesh(importPath.toStdString(), extractedGroups);
    QApplication::restoreOverrideCursor();

    if (meshData.empty()) {
        QMessageBox::warning(this, "导入失败", "未能从该文件中提取出任何表面网格！");
        logMessage("导入网格失败：无有效表面数据。", "error");
        return;
    }

    QTreeWidgetItem* importedNode = new QTreeWidgetItem(m_rootMeshNode);
    importedNode->setText(0, QFileInfo(filePath).baseName() + " (Imported)");

    m_meshFiles[importedNode] = importPath;
    m_isImportedMesh[importedNode] = true;
    m_meshDataMap[importedNode] = meshData;
    MeshEngine::saveMeshCache((importPath + ".meshcache").toStdString(), meshData);
    // 提取物理组→网格映射, 写入 .pgroups.json
    if (!extractedGroups.empty()) {
        try {
            MeshEngine::ensureInit();
            gmsh::clear();
            gmsh::open(importPath.toStdString());
            MeshEngine::savePhysicalGroupsJson(importPath.toStdString());
            gmsh::clear();
        } catch (...) { try { gmsh::clear(); } catch (...) {} }
    }
    m_projectModified = true;
    m_rootMeshNode->setExpanded(true);

    // 保存当前内部网格的物理组
    if (m_activeMeshNode && m_rootMeshNode->indexOfChild(m_activeMeshNode) >= 0) {
        std::vector<GroupMeta> meta;
        for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
            auto* ci = m_rootGroupNode->child(i);
            GroupMeta gm{ci->text(0), m_groupTags.count(ci) ? m_groupTags[ci] : 0, {}};
            if (m_groupData.count(ci)) gm.faceIds = m_groupData[ci];
            meta.push_back(gm);
        }
        m_importedGroupMeta[m_activeMeshNode] = meta;
        std::map<std::string, BoundaryParams> bdMap;
        for (auto& [node, params] : m_boundaryConfigs)
            bdMap[node->text(0).toStdString()] = params;
        m_importedBoundaryConfigs[m_activeMeshNode] = bdMap;
    }
    for(int i = m_rootGroupNode->childCount() - 1; i >= 0; --i) delete m_rootGroupNode->child(i);
    m_groupData.clear();
    m_boundaryConfigs.clear();
    m_meshConfigs.clear();

    // 构建 face ID ↔ entity tag 映射 (BREP 面序 = Gmsh entity 序)
    m_faceToEntity.clear(); m_entityToFaces.clear(); m_elemToEntity.clear();
    {
        MeshEngine::ensureInit();
        gmsh::clear();
        gmsh::open(importPath.toStdString());
        std::vector<std::pair<int,int>> allEnts;
        gmsh::model::getEntities(allEnts, 2);
        std::sort(allEnts.begin(), allEnts.end(),
                  [](auto& a, auto& b){ return a.second < b.second; });
        int faceId = 1;
        for (auto& ent : allEnts) {
            int entTag = ent.second;
            m_faceToEntity[faceId] = entTag;
            m_entityToFaces[entTag].insert(faceId);
            std::vector<int> eTypes;
            std::vector<std::vector<std::size_t>> eTags, nTags;
            gmsh::model::mesh::getElements(eTypes, eTags, nTags, 2, entTag);
            for (auto& tags : eTags)
                for (auto t : tags) m_elemToEntity[(int)t] = entTag;
            faceId++;
        }
        gmsh::clear();
    }

    // 从 .pgroups.json 读取并映射到 viewer face ID
    std::map<int, std::set<int>> pgFaceMap; // PG tag → viewer face IDs
    {
        std::ifstream pgf((importPath + ".pgroups.json").toStdString());
        if (pgf.is_open()) {
            json pgj = json::parse(pgf);
            for (auto& g : pgj["groups"]) {
                int tag = g["tag"]; int dim = g["dim"];
                if (dim == 2 && g.contains("entities")) {
                    for (auto& e : g["entities"]) {
                        int entTag = e.get<int>();
                        if (m_entityToFaces.count(entTag))
                            pgFaceMap[tag].insert(m_entityToFaces[entTag].begin(),
                                                  m_entityToFaces[entTag].end());
                    }
                }
            }
        }
    }

    std::vector<MainWindow::GroupMeta> groupMeta;
    for (const auto& grp : extractedGroups) {
        if (grp.dimension == 2) {
            QTreeWidgetItem* child = new QTreeWidgetItem(m_rootGroupNode);
            child->setText(0, QString::fromStdString(grp.name));
            child->setData(0, Qt::UserRole, grp.tag);
            child->setFlags(child->flags() | Qt::ItemIsEditable);
            // m_groupData 存 viewer face ID (用于高亮)
            m_groupData[child] = pgFaceMap.count(grp.tag) ? pgFaceMap[grp.tag] : std::set<int>();
            // 存储 entity tag 列表用于面板显示 (纯数字)
            {
                std::set<int> entities;
                for (int fid : m_groupData[child])
                    if (m_faceToEntity.count(fid))
                        entities.insert(m_faceToEntity[fid]);
                QStringList nums;
                for (int e : entities) nums << QString::number(e);
                child->setData(0, Qt::UserRole + 1, nums.join(","));
            }
            m_groupTags[child] = grp.tag;
            groupMeta.push_back({QString::fromStdString(grp.name), grp.tag});
        }
    }
    m_importedGroupMeta[importedNode] = groupMeta;
    m_activeMeshNode = importedNode;
    if (m_btnPartition) m_btnPartition->setEnabled(true);
    m_rootGroupNode->setExpanded(true);
    refreshBoundaryUI();

    ui->cadViewerWidget->setGeometryHidden(true); // 导入网格时隐藏几何体
    ui->treeGroups->setCurrentItem(importedNode);

    QString msg = QString("成功加载外部网格，并解析出 %1 个边界物理组。").arg(m_rootGroupNode->childCount());
    QMessageBox::information(this, "导入成功", msg);
    logMessage(msg, "success");
}

// 【关键修改】：处理后处理节点的显隐切换
void MainWindow::onTreeSelectionChanged() {
    QTreeWidgetItem* current = ui->treeGroups->currentItem();
    if (!current || current == m_rootGroupNode || current == m_rootMeshNode ||
        current == m_rootGeomNode) {
        ui->cadViewerWidget->setResultVisible(false);
        ui->cadViewerWidget->setGeometryHidden(false);
        ui->stackedWidgetSettings->setCurrentIndex(ui->stackedWidgetSettings->count() - 1);
        ui->dockSettings->show();
        return;
    }

    if (current && current->parent() == m_rootGeomNode) {
        // 自动保存上一个节点的参数
        if (m_lastGeomIdx >= 0 && m_lastGeomIdx < (int)m_geomNodes.size())
            saveGeomNodeFromUI(m_lastGeomIdx);
        // 关闭所有打开的添加按钮
        if (m_geomStackPages) {
            for (int p = 0; p < m_geomStackPages->count(); p++) {
                QWidget* pg = m_geomStackPages->widget(p);
                if (!pg) continue;
                QList<QPushButton*> btns = pg->findChildren<QPushButton*>();
                for (auto* b : btns)
                    if (b->isCheckable() && b->isChecked())
                        b->setChecked(false);
            }
        }
        if (m_pageGeomParam && m_geomStackPages) {
            ui->stackedWidgetSettings->setCurrentWidget(m_pageGeomParam);
            int idx = current->data(0, Qt::UserRole).toInt();
            m_lastGeomIdx = idx;
            QString t = current->data(0, Qt::UserRole+1).toString();
            if (t == "并集" || t == "交集") m_geomStackPages->setCurrentIndex(1);
            else if (t == "差集") m_geomStackPages->setCurrentIndex(2);
            else if (t == "平移") m_geomStackPages->setCurrentIndex(3);
            else if (t == "旋转") m_geomStackPages->setCurrentIndex(4);
            else if (t == "镜像") m_geomStackPages->setCurrentIndex(5);
            else if (t == "缩放") m_geomStackPages->setCurrentIndex(6);
            else if (t == "偏移") m_geomStackPages->setCurrentIndex(7);
            else if (t == "倒圆角") m_geomStackPages->setCurrentIndex(8);
            else if (t == "倒斜角") m_geomStackPages->setCurrentIndex(9);
            else if (t == "扫掠") m_geomStackPages->setCurrentIndex(10);
            else if (t == "线性阵列") m_geomStackPages->setCurrentIndex(11);
            else if (t == "圆形阵列") m_geomStackPages->setCurrentIndex(12);
            else m_geomStackPages->setCurrentIndex(0);
            loadGeomNodeToUI(idx);
        } else {
            ui->stackedWidgetSettings->setCurrentIndex(ui->stackedWidgetSettings->count()-1);
        }
        ui->cadViewerWidget->clearMesh();
        ui->cadViewerWidget->setGeometryTransparent(false);
        ui->cadViewerWidget->setGeometryHidden(false);
        ui->dockSettings->show();
        return;
    }

    if (current == m_rootBoundaryNode) {
        refreshBoundaryUI();
        ui->stackedWidgetSettings->setCurrentWidget(m_pageBoundary);
        ui->dockSettings->show();
        return;
    }

    if (current == m_rootSolverNode) {
        ui->stackedWidgetSettings->setCurrentWidget(m_pageSolver);
        refreshSolverMeshCombo();
        // 显示选中的网格 (如果有)
        if (m_activeMeshNode && m_meshDataMap.count(m_activeMeshNode)) {
            ui->cadViewerWidget->loadMesh(m_meshDataMap[m_activeMeshNode]);
            ui->cadViewerWidget->setGeometryTransparent(true);
        }
        ui->cadViewerWidget->setResultVisible(false);
        ui->cadViewerWidget->setGeometryHidden(false);
        ui->dockSettings->show();
        return;
    }

    if (current == m_rootPostNode) {
        ui->stackedWidgetSettings->setCurrentWidget(m_pagePostProcess);
        refreshPostProcessGroups();
        if (ui->cadViewerWidget->hasResultActor()) {
            ui->cadViewerWidget->setGeometryHidden(true);
            ui->cadViewerWidget->setResultVisible(true);
        } else {
            // 显示选中的网格 (如果有)
            if (m_activeMeshNode && m_meshDataMap.count(m_activeMeshNode)) {
                ui->cadViewerWidget->loadMesh(m_meshDataMap[m_activeMeshNode]);
                ui->cadViewerWidget->setGeometryTransparent(true);
            }
        }
        ui->dockSettings->show();
        return;
    }

    if (current == m_rootParamNode) {
        ui->stackedWidgetSettings->setCurrentWidget(m_pageParams);
        ui->dockSettings->show();
        return;
    }

    if (current->parent() == m_rootMeshNode) {
        ui->cadViewerWidget->setResultVisible(false);
        ui->cadViewerWidget->setGeometryHidden(false);
        ui->stackedWidgetSettings->setCurrentIndex(1);
        ui->cadViewerWidget->setSelectedFaces(std::set<int>());

        bool isImported = m_isImportedMesh[current];
        ui->groupMeshTarget->setVisible(!isImported);
        ui->groupBoxMesh->setVisible(!isImported);
        ui->groupGlobalMesh->setVisible(!isImported);
        ui->btnGenerateMesh->setVisible(!isImported);
        ui->btnExportMesh->setText(isImported ? "另存为网格" : "导出当前网格至 .msh 文件");

        // 物理组切换: 内部网格共享组, 导入网格专用组
        if (current != m_activeMeshNode) {
            bool oldImported = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode) && m_isImportedMesh[m_activeMeshNode];
            bool newImported = isImported && m_isImportedMesh[current];

            // 保存当前组到旧网格节点 (离开旧节点前)
            if (m_activeMeshNode) {
                std::vector<GroupMeta> meta;
                for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
                    auto* ci = m_rootGroupNode->child(i);
                    GroupMeta gm{ci->text(0), m_groupTags.count(ci)?m_groupTags[ci]:0, {}};
                    if (m_groupData.count(ci)) gm.faceIds = m_groupData[ci];
                    meta.push_back(gm);
                }
                m_importedGroupMeta[m_activeMeshNode] = meta;
                std::map<std::string, BoundaryParams> bdMap;
                for (auto& [node, params] : m_boundaryConfigs)
                    bdMap[node->text(0).toStdString()] = params;
                m_importedBoundaryConfigs[m_activeMeshNode] = bdMap;
            }

            if (newImported) {
                for (int i = m_rootGroupNode->childCount()-1; i >= 0; i--)
                    delete m_rootGroupNode->child(i);
                m_groupData.clear(); m_boundaryConfigs.clear(); m_groupTags.clear();
                if (m_importedGroupMeta.count(current)) {
                    for (auto& gm : m_importedGroupMeta[current]) {
                        auto* item = new QTreeWidgetItem(m_rootGroupNode);
                        item->setText(0, gm.name);
                        if (gm.tag > 0) { item->setData(0, Qt::UserRole, gm.tag); m_groupTags[item] = gm.tag; }
                        m_groupData[item] = gm.faceIds;
                    }
                }
                if (m_importedBoundaryConfigs.count(current)) {
                    auto& saved = m_importedBoundaryConfigs[current];
                    for (auto& [name, params] : saved)
                        for (int i = 0; i < m_rootGroupNode->childCount(); i++)
                            if (m_rootGroupNode->child(i)->text(0).toStdString() == name)
                                { m_boundaryConfigs[m_rootGroupNode->child(i)] = params; break; }
                }
            } else if (oldImported && !newImported) {
                // 退出导入网格: 删除导入组, 恢复保存的内部物理组
                for (int i = m_rootGroupNode->childCount()-1; i >= 0; i--)
                    delete m_rootGroupNode->child(i);
                m_groupData.clear(); m_boundaryConfigs.clear(); m_groupTags.clear();
                if (m_importedGroupMeta.count(current)) {
                    for (auto& gm : m_importedGroupMeta[current]) {
                        auto* item = new QTreeWidgetItem(m_rootGroupNode);
                        item->setText(0, gm.name);
                        if (gm.tag > 0) { item->setData(0, Qt::UserRole, gm.tag); m_groupTags[item] = gm.tag; }
                        m_groupData[item] = gm.faceIds;
                    }
                }
                if (m_importedBoundaryConfigs.count(current)) {
                    auto& saved = m_importedBoundaryConfigs[current];
                    for (auto& [name, params] : saved)
                        for (int i = 0; i < m_rootGroupNode->childCount(); i++)
                            if (m_rootGroupNode->child(i)->text(0).toStdString() == name)
                                { m_boundaryConfigs[m_rootGroupNode->child(i)] = params; break; }
                }
            }

            m_activeMeshNode = current;
            m_rootGroupNode->setExpanded(true);
            refreshBoundaryUI();
        }

        if (m_meshDataMap.count(current) && !m_meshDataMap[current].empty()) {
            ui->cadViewerWidget->loadMesh(m_meshDataMap[current]);
            ui->cadViewerWidget->setGeometryHidden(true);
            ui->cadViewerWidget->resetCameraView();
        } else {
            ui->cadViewerWidget->clearMesh();
            ui->cadViewerWidget->setGeometryTransparent(false);
        }

        m_isUpdatingUI = true;
        ui->comboMeshGroupTarget->clear();
        for (int i = 0; i < m_rootGroupNode->childCount(); ++i) {
            QTreeWidgetItem* groupNode = m_rootGroupNode->child(i);
            if (m_meshConfigs[current].find(groupNode) == m_meshConfigs[current].end()) {
                m_meshConfigs[current][groupNode] = MeshGroupParams();
            }
            ui->comboMeshGroupTarget->addItem(groupNode->text(0), QVariant::fromValue(static_cast<void*>(groupNode)));
        }
        m_isUpdatingUI = false;

        if (ui->comboMeshGroupTarget->count() > 0) onMeshGroupTargetChanged(0);
        ui->dockSettings->show();

    } else if (current->parent() == m_rootGroupNode) {
        bool lockGroups = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode)
                          && m_isImportedMesh[m_activeMeshNode];
        ui->cadViewerWidget->setResultVisible(false);
        if (!lockGroups) {
            ui->cadViewerWidget->setGeometryHidden(false);
            ui->cadViewerWidget->clearMesh();
            ui->cadViewerWidget->setGeometryTransparent(false);
        }
        if (lockGroups && m_activeMeshNode && m_meshDataMap.count(m_activeMeshNode)) {
            ui->cadViewerWidget->loadMesh(m_meshDataMap[m_activeMeshNode]);
            ui->cadViewerWidget->setGeometryHidden(true);
        }
        ui->stackedWidgetSettings->setCurrentIndex(0);
        ui->cadViewerWidget->setSelectedFaces(m_groupData[current]);
        m_isUpdatingUI = true;
        ui->editGroupName->setEnabled(!lockGroups);
        ui->editGroupName->setText(current->text(0));
        ui->listSelectedFaces->clear();
        if (lockGroups) {
            QString entStr = current->data(0, Qt::UserRole + 1).toString();
            if (!entStr.isEmpty())
                for (auto& s : entStr.split(","))
                    ui->listSelectedFaces->addItem(s.trimmed());
        } else {
            for (int id : m_groupData[current]) ui->listSelectedFaces->addItem(QString::number(id));
        }
        m_isUpdatingUI = false;
        ui->dockSettings->show(); 
    }
}

void MainWindow::showTreeContextMenu(const QPoint& pos) {
    QTreeWidgetItem* item = ui->treeGroups->itemAt(pos);
    if (!item) return;

    QMenu menu(this);
    if (item->parent() == m_rootGeomNode) {
        QAction* actDel = menu.addAction("删除此节点");
        connect(actDel, &QAction::triggered, this, [=](){
            int idx = item->data(0, Qt::UserRole).toInt();
            if (idx >= 0 && idx < (int)m_geomNodes.size()) {
                m_geomNodes.erase(m_geomNodes.begin() + idx);
                delete item;
                for (int i = idx; i < (int)m_geomNodes.size(); i++)
                    m_geomNodes[i].treeItem->setData(0, Qt::UserRole, i);
                m_projectModified = true;
            }
        });
    } else if (item == m_rootGroupNode) {
        QAction* actAdd = menu.addAction("新建物理边界");
        connect(actAdd, &QAction::triggered, this, [=](){
            QTreeWidgetItem* child = new QTreeWidgetItem(m_rootGroupNode);
            child->setText(0, "Boundary_" + QString::number(m_rootGroupNode->childCount()));
            child->setFlags(child->flags() | Qt::ItemIsEditable);
            m_groupData[child] = std::set<int>();
            m_projectModified = true;
            m_rootGroupNode->setExpanded(true);
            ui->treeGroups->setCurrentItem(child);
        });
    } else if (item == m_rootMeshNode) {
        QAction* actAddMesh = menu.addAction("新建网格生成方案");
        connect(actAddMesh, &QAction::triggered, this, [=](){
            QTreeWidgetItem* child = new QTreeWidgetItem(m_rootMeshNode);
            child->setText(0, "Mesh_Scheme_" + QString::number(m_rootMeshNode->childCount()));
            std::map<QTreeWidgetItem*, MeshGroupParams> initialParams;
            for(int i=0; i<m_rootGroupNode->childCount(); ++i) {
                initialParams[m_rootGroupNode->child(i)] = MeshGroupParams();
            }
            m_meshConfigs[child] = initialParams;
            m_isImportedMesh[child] = false;
            m_projectModified = true;
            m_rootMeshNode->setExpanded(true);
            ui->treeGroups->setCurrentItem(child);
        });
    } else if (item->parent() == m_rootGroupNode) {
        QAction* actDelete = menu.addAction("删除此边界");
        connect(actDelete, &QAction::triggered, this, [=](){
            m_projectModified = true;
            m_groupData.erase(item);
            m_boundaryConfigs.erase(item);
            for(auto& pair : m_meshConfigs) pair.second.erase(item);
            if (ui->treeGroups->currentItem() == item) ui->cadViewerWidget->setSelectedFaces(std::set<int>()); 
            delete item;
        });
    } else if (item->parent() == m_rootMeshNode) {
        QAction* actDelete = menu.addAction("删除此网格节点");
        connect(actDelete, &QAction::triggered, this, [=](){
            m_projectModified = true;
            m_meshConfigs.erase(item);
            m_meshDataMap.erase(item);
            m_meshFiles.erase(item);
            m_isImportedMesh.erase(item);
            if (ui->treeGroups->currentItem() == item) {
                ui->cadViewerWidget->clearMesh();
                ui->cadViewerWidget->setGeometryTransparent(false);
            }
            delete item;
        });
    }
    menu.exec(ui->treeGroups->mapToGlobal(pos));
}

void MainWindow::onMeshGroupTargetChanged(int index) {
    if (m_isUpdatingUI || index < 0) return;
    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    if (!currentMesh || currentMesh->parent() != m_rootMeshNode) return;
    QTreeWidgetItem* targetGroup = static_cast<QTreeWidgetItem*>(ui->comboMeshGroupTarget->itemData(index).value<void*>());
    if (!targetGroup) return;

    m_isUpdatingUI = true;
    const MeshGroupParams& params = m_meshConfigs[currentMesh][targetGroup];
    ui->chkEnableLocalField->setChecked(params.enableLocalField);
    ui->groupBoxMesh->setEnabled(params.enableLocalField);
    ui->spinSizeMin->setValue(params.sizeMin);
    ui->spinSizeMax->setValue(params.sizeMax);
    ui->spinDistMin->setValue(params.distMin);
    ui->spinDistMax->setValue(params.distMax);
    m_isUpdatingUI = false;
}

void MainWindow::onEnableLocalFieldToggled(bool checked) {
    ui->groupBoxMesh->setEnabled(checked); 
    if (m_isUpdatingUI) return;
    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    int index = ui->comboMeshGroupTarget->currentIndex();
    if (currentMesh && currentMesh->parent() == m_rootMeshNode && index >= 0) {
        QTreeWidgetItem* targetGroup = static_cast<QTreeWidgetItem*>(ui->comboMeshGroupTarget->itemData(index).value<void*>());
        if (targetGroup) {
            m_meshConfigs[currentMesh][targetGroup].enableLocalField = checked;
            m_meshDataMap.erase(currentMesh);
            ui->cadViewerWidget->clearMesh();
            ui->cadViewerWidget->setGeometryTransparent(false);
        }
    }
}

void MainWindow::onMeshParamsChanged() {
    if (m_isUpdatingUI) return;
    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    int index = ui->comboMeshGroupTarget->currentIndex();
    if (currentMesh && currentMesh->parent() == m_rootMeshNode && index >= 0) {
        QTreeWidgetItem* targetGroup = static_cast<QTreeWidgetItem*>(ui->comboMeshGroupTarget->itemData(index).value<void*>());
        if (targetGroup) {
            MeshGroupParams& params = m_meshConfigs[currentMesh][targetGroup];
            params.sizeMin = ui->spinSizeMin->value();
            params.sizeMax = ui->spinSizeMax->value();
            params.distMin = ui->spinDistMin->value();
            params.distMax = ui->spinDistMax->value();
            m_meshDataMap.erase(currentMesh);
            ui->cadViewerWidget->clearMesh();
            ui->cadViewerWidget->setGeometryTransparent(false);
        }
    }
}

void MainWindow::onFaceSelectionChanged(const std::set<int>& selectedIds) {
    static bool s_updating = false;
    if (s_updating) return;
    QTreeWidgetItem* current = ui->treeGroups->currentItem();
    if (!current || current->parent() != m_rootGroupNode) return;
    bool isImported = m_activeMeshNode && m_isImportedMesh.count(m_activeMeshNode)
                      && m_isImportedMesh[m_activeMeshNode];
    if (isImported && !m_faceToEntity.empty()) {
        // 空集合 → 清空
        if (selectedIds.empty()) {
            m_groupData[current].clear();
            current->setData(0, Qt::UserRole + 1, "");
            ui->listSelectedFaces->clear();
            s_updating = true; ui->cadViewerWidget->setSelectedFaces({}); s_updating = false;
            return;
        }
        // 少量点击 → entity 级 toggle; 批量 → 回显 (size > 3 = click on PG tree node)
        if (selectedIds.size() <= 3) {
            // 维护 entity tag 集合避免每次遍历数千 face
            std::set<int>& faces = m_groupData[current];
            std::set<int> ents; for (int f : faces) if (m_faceToEntity.count(f)) ents.insert(m_faceToEntity[f]);

            // 计算增量: 新增面(在selectedIds但不在faces) vs 移除面(在faces但不在selectedIds)
            std::set<int> addedFaceIds, removedFaceIds;
            for (int fid : selectedIds) if (!faces.count(fid)) addedFaceIds.insert(fid);
            for (int fid : faces)       if (!selectedIds.count(fid)) removedFaceIds.insert(fid);

            // 处理移除: 面被取消选中 → 移除其所属实体全部面
            for (int fid : removedFaceIds) {
                auto it = m_faceToEntity.find(fid);
                if (it == m_faceToEntity.end()) continue;
                int ent = it->second;
                if (ents.count(ent)) {
                    ents.erase(ent);
                    for (int ef : m_entityToFaces[ent]) faces.erase(ef);
                }
            }

            // 处理添加: 新增面 → 添加其所属实体全部面
            for (int fid : addedFaceIds) {
                auto it = m_faceToEntity.find(fid);
                if (it == m_faceToEntity.end()) continue;
                int ent = it->second;
                if (!ents.count(ent)) {
                    ents.insert(ent);
                    faces.insert(m_entityToFaces[ent].begin(), m_entityToFaces[ent].end());
                }
            }

            // 更新显示
            QStringList nums; for (int e : ents) nums << QString::number(e);
            current->setData(0, Qt::UserRole + 1, nums.join(","));
            if(!m_isUpdatingUI) { ui->listSelectedFaces->clear(); for (auto& s : nums) ui->listSelectedFaces->addItem(s); }
            s_updating = true; ui->cadViewerWidget->setSelectedFaces(faces); s_updating = false;
        } else {
            m_groupData[current] = selectedIds;
            std::set<int> ents; for (int f : selectedIds) if (m_faceToEntity.count(f)) ents.insert(m_faceToEntity[f]);
            QStringList nums; for (int e : ents) nums << QString::number(e);
            current->setData(0, Qt::UserRole + 1, nums.join(","));
            if(!m_isUpdatingUI) { ui->listSelectedFaces->clear(); for (auto& s : nums) ui->listSelectedFaces->addItem(s); }
        }
    } else {
        m_groupData[current] = selectedIds;
        if(!m_isUpdatingUI) {
            ui->listSelectedFaces->clear();
            for (int id : selectedIds) ui->listSelectedFaces->addItem(QString::number(id));
        }
    }
}

void MainWindow::onGroupNameEdited(const QString& newName) {
    QTreeWidgetItem* current = ui->treeGroups->currentItem();
    if (current && current->parent() == m_rootGroupNode) {
        current->setText(0, newName); 
        if (ui->stackedWidgetSettings->currentIndex() == 1) {
            int index = ui->comboMeshGroupTarget->currentIndex();
            onTreeSelectionChanged(); 
            ui->comboMeshGroupTarget->setCurrentIndex(index);
        }
    }
}

bool MainWindow::runMeshEngine(const QString& savePath) {
    if (m_isMeshing) {
        QMessageBox::warning(this, "警告", "网格生成正在进行中，请耐心等待！");
        return false;
    }
    if (m_currentStepFile.isEmpty()) return false;

    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    if (!currentMesh || currentMesh->parent() != m_rootMeshNode) {
        QMessageBox::warning(this, "提示", "请在左侧选择一个具体的网格方案 (Mesh_X) 再进行生成。");
        return false;
    }
    
    if (m_isImportedMesh[currentMesh]) {
        QMessageBox::warning(this, "警告", "外部导入的网格无法进行重新生成！");
        return false;
    }

    std::map<std::string, std::pair<std::set<int>, MeshGroupParams>> engineData;
    for (int i=0; i<m_rootGroupNode->childCount(); ++i) {
        QTreeWidgetItem* groupNode = m_rootGroupNode->child(i);
        std::set<int> faceIds = m_groupData[groupNode];
        if (!faceIds.empty()) {
            engineData[groupNode->text(0).toStdString()] = {faceIds, m_meshConfigs[currentMesh][groupNode]};
        }
    }

    int order = 1; 
    bool optimize = ui->checkOptimize->isChecked();
    int numThreads = ui->spinMeshThreads->value(); 

    QString workingDir = QApplication::applicationDirPath();
    QString cacheDir = workingDir + "/.cache";
    QDir().mkpath(cacheDir); 
    QString cacheFile = cacheDir + "/" + currentMesh->text(0) + ".msh";

    m_meshSavePath = savePath;
    m_meshCachePath = cacheFile;
    m_meshingNode = currentMesh;
    m_isMeshing = true;

    m_progressBar->setVisible(true);
    m_simulatedProgress = 0.0;
    m_progressBar->setValue(0);
    m_btnCancelMesh->setText("终止操作");
    m_btnCancelMesh->setVisible(true);
    m_btnCancelMesh->setEnabled(true);
    ui->btnGenerateMesh->setEnabled(false);
    ui->btnExportMesh->setEnabled(false);
    
    m_meshLogTimer->start(100);

    logMessage("正在后台并行生成 3D 网格，请稍候...", "highlight");

    // 计算内部几何首面ID: 取viewer中实际可见的最小面ID
    int firstFaceId = 0;
    for (auto& [fid, oid] : ui->cadViewerWidget->faceObjectMap())
        if (firstFaceId == 0 || fid < firstFaceId) firstFaceId = fid;

    std::string stepPathStr = m_currentStepFile.toStdString();
    std::string cacheFileStr = cacheFile.toStdString();

    m_meshThread = std::thread([this, stepPathStr, engineData, order, optimize, numThreads, cacheFileStr, firstFaceId]() {
        std::string meshLog;
        FaceMeshMap meshData = MeshEngine::generateMesh(
            stepPathStr, engineData, firstFaceId, order, optimize, numThreads, cacheFileStr, meshLog
        );
        this->m_lastMeshData = meshData;
        this->m_lastMeshLog = meshLog;
        QMetaObject::invokeMethod(this, "onMeshFinished", Qt::QueuedConnection);
    });

    return true;
}

void MainWindow::onMeshFinished() {
    fprintf(stderr, "[DEBUG onMeshFinished] called, joinable=%d\n", (int)m_meshThread.joinable());
    if (m_meshThread.joinable()) {
        m_meshThread.join();
    }
    
    m_meshLogTimer->stop();
    m_progressBar->setValue(100);
    m_isMeshing = false;
    m_progressBar->setVisible(false);
    m_btnCancelMesh->setVisible(false);
    ui->btnGenerateMesh->setEnabled(true);
    ui->btnExportMesh->setEnabled(true);
    if (m_btnPartition) m_btnPartition->setEnabled(true);

    if (!m_lastMeshData.empty()) {
        logMessage(QString::fromStdString(m_lastMeshLog), "success");
        m_projectModified = true;
        m_meshFiles[m_meshingNode] = m_meshCachePath;
        m_meshDataMap[m_meshingNode] = m_lastMeshData;
        MeshEngine::saveMeshCache((m_meshCachePath + ".meshcache").toStdString(), m_lastMeshData);

        // 从 .msh 同步 Gmsh 物理组 tag (复用已初始化的 Gmsh)
        try {
            gmsh::open(m_meshCachePath.toStdString());
            std::vector<std::pair<int,int>> dimTags;
            gmsh::model::getPhysicalGroups(dimTags);
            for (auto& dt : dimTags) {
                if (dt.first != 2) continue;
                std::string name;
                gmsh::model::getPhysicalName(dt.first, dt.second, name);
                for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
                    if (m_rootGroupNode->child(i)->text(0).toStdString() == name) {
                        m_groupTags[m_rootGroupNode->child(i)] = dt.second;
                        m_rootGroupNode->child(i)->setData(0, Qt::UserRole, dt.second);
                    }
                }
            }
            gmsh::clear();
        } catch (...) { try { gmsh::clear(); } catch (...) {} }
        
        if (ui->treeGroups->currentItem() == m_meshingNode) {
            ui->cadViewerWidget->loadMesh(m_lastMeshData); 
            ui->cadViewerWidget->setGeometryTransparent(true);
        }
        
        if (!m_meshSavePath.isEmpty()) {
            if (QFile::exists(m_meshSavePath)) QFile::remove(m_meshSavePath);
            if (QFile::copy(m_meshCachePath, m_meshSavePath)) {
                QMessageBox::information(this, "成功", "网格生成并导出成功！\n" + m_meshSavePath); 
            } else {
                QMessageBox::warning(this, "错误", "文件导出失败！");
            }
        }
    } else {
        logMessage("[失败] 网格生成被取消或发生错误。", "error");
        if (!m_meshSavePath.isEmpty()) {
            QMessageBox::warning(this, "失败", "网格未生成且生成操作失败！");
        }
    }
}

void MainWindow::generateMesh() { 
    if (!runMeshEngine("")) QMessageBox::critical(this, "失败", "网格生成失败，请检查模型参数或几何体。"); 
}

void MainWindow::exportMeshToFile(const QString& savePath) {
    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    if (!currentMesh || currentMesh->parent() != m_rootMeshNode) return;
    QString bindFile = m_meshFiles[currentMesh];
    if (bindFile.isEmpty() || !QFile::exists(bindFile)) {
        if (!m_isImportedMesh[currentMesh]) {
            if (!runMeshEngine(savePath))
                QMessageBox::warning(this, "失败", "网格未生成且生成操作失败！");
            else
                QMessageBox::information(this, "成功", "网格生成并导出成功！\n" + savePath);
        } else {
            QMessageBox::warning(this, "错误", "导入的网格文件源已丢失！");
        }
        return;
    }
    QString suffix = QFileInfo(savePath).suffix().toLower();
    if (suffix == "msh") {
        if (QFile::exists(savePath)) QFile::remove(savePath);
        if (QFile::copy(bindFile, savePath))
            QMessageBox::information(this, "成功", "网格文件保存成功！\n" + savePath);
        else
            QMessageBox::warning(this, "错误", "文件写入失败！");
    } else {
        try {
            gmsh::open(bindFile.toStdString());
            gmsh::write(savePath.toStdString());
            gmsh::clear();
            QMessageBox::information(this, "成功", "网格导出成功！\n" + savePath);
        } catch (std::exception& e) {
            try { gmsh::clear(); } catch (...) {}
            QMessageBox::warning(this, "错误", QString("导出失败: %1").arg(e.what()));
        }
    }
}

void MainWindow::exportMesh() {
    QTreeWidgetItem* currentMesh = ui->treeGroups->currentItem();
    if (!currentMesh || currentMesh->parent() != m_rootMeshNode) return;

    QDialog dlg(this);
    dlg.setWindowTitle("导出网格");
    auto* layout = new QVBoxLayout(&dlg);

    layout->addWidget(new QLabel("导出格式:"));
    auto* formatCombo = new QComboBox();
    struct Fmt { QString label; QString ext; bool hasPG; };
    std::vector<Fmt> formats = {
        {"Gmsh v2 ASCII (.msh)",     "msh",  true},
        {"I-DEAS UNV (.unv)",        "unv",  true},
        {"Abaqus INP (.inp)",        "inp",  true},
        {"Tochnog (.dat)",           "dat",  true},
        {"VTK Legacy (.vtk)",        "vtk",  false},
        {"INRIA Medit (.mesh)",      "mesh", false},
        {"STL ASCII (.stl)",         "stl",  false},
        {"Nastran BDF (.bdf)",       "bdf",  false},
        {"Diffpack (.diff)",         "diff", false},
    };
    for (auto& f : formats)
        formatCombo->addItem(f.label);
    layout->addWidget(formatCombo);

    auto* pathLayout = new QHBoxLayout();
    auto* pathEdit = new QLineEdit("IonFlowMesh.msh");
    auto* browseBtn = new QPushButton("浏览...");
    pathLayout->addWidget(pathEdit);
    pathLayout->addWidget(browseBtn);
    layout->addLayout(pathLayout);

    auto* btnLayout = new QHBoxLayout();
    auto* okBtn = new QPushButton("导出");
    auto* cancelBtn = new QPushButton("取消");
    btnLayout->addStretch(); btnLayout->addWidget(okBtn); btnLayout->addWidget(cancelBtn);
    layout->addLayout(btnLayout);

    auto* pgInfo = new QLabel();
    layout->addWidget(pgInfo);

    connect(formatCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int idx) {
        QString base = QFileInfo(pathEdit->text()).completeBaseName();
        pathEdit->setText(base + "." + formats[idx].ext);
        pgInfo->setText(formats[idx].hasPG ? "✓ 此格式支持内嵌物理组" : "✗ 此格式不支持物理组 (仅导出网格几何)");
    });
    // trigger initial update
    pgInfo->setText(formats[0].hasPG ? "✓ 此格式支持内嵌物理组" : "✗ 此格式不支持物理组 (仅导出网格几何)");

    connect(browseBtn, &QPushButton::clicked, [&]() {
        QString filter = "All Mesh Files (*.msh *.vtk *.mesh *.stl *.unv *.inp *.bdf *.diff *.dat);;All Files (*)";
        QString path = QFileDialog::getSaveFileName(&dlg, "选择保存路径", pathEdit->text(), filter);
        if (!path.isEmpty()) pathEdit->setText(path);
    });
    connect(okBtn, &QPushButton::clicked, &dlg, &QDialog::accept);
    connect(cancelBtn, &QPushButton::clicked, &dlg, &QDialog::reject);

    if (dlg.exec() == QDialog::Accepted) {
        QString savePath = pathEdit->text();
        if (!savePath.isEmpty())
            exportMeshToFile(savePath);
    }
}

// ====================================================================
// 从 GeomNode 加载/保存到 UI
// ====================================================================
void MainWindow::loadGeomNodeToUI(int idx) {
    if (idx < 0 || idx >= (int)m_geomNodes.size()) return;
    auto& n = m_geomNodes[idx];
    // 切换节点绑定: 清除全局绑定, 恢复该节点的绑定
    for (auto& [spin, pname] : m_paramBindings) { spin->setStyleSheet(""); spin->setReadOnly(false); }
    for (auto& [btn, spin] : m_bindBtns)
        btn->setStyleSheet("QToolButton { border:1px solid #CCC; background:#EEE; font-weight:bold; }");
    m_paramBindings.clear();
    for (auto& [sname, pname] : n.nodeParamBindings) {
        for (auto& [spin, name] : m_spinboxNames)
            if (name == sname) { m_paramBindings[spin] = pname; spin->setStyleSheet("QDoubleSpinBox { background-color: #E3F2FD; }"); spin->setReadOnly(true); break; }
    }
    m_geomW->setValue(n.w); m_geomH->setValue(n.h); m_geomD->setValue(n.d);
    // 动态更新尺寸标签 + 显隐
    QFormLayout* fl = qobject_cast<QFormLayout*>(m_geomStackPages->widget(0)->layout());
    if (fl) {
        auto setRow = [&](int row, const QString& text, bool vis) {
            auto* labelItem = fl->itemAt(row, QFormLayout::LabelRole);
            auto* fieldItem = fl->itemAt(row, QFormLayout::FieldRole);
            if (labelItem && labelItem->widget())
                if (auto* lb = qobject_cast<QLabel*>(labelItem->widget())) lb->setText(text);
            if (fieldItem && fieldItem->widget()) fieldItem->widget()->setVisible(vis);
        };
        // 第4行默认隐藏
        m_geomE->hide(); setRow(3, "", false);
        auto setDim = [&](int row, const QString& t, bool v){ setRow(row,t,v); };
        auto setDim4 = [&](const QString& t, bool v) {
            setRow(3, t, v); m_geomE->setVisible(v);
        };
        if (n.type=="立方体") { setDim(0,"长 (W):",1); setDim(1,"宽 (H):",1); setDim(2,"高 (D):",1); }
        else if (n.type=="球体") { setDim(0,"半径 (R):",1); setDim(1,"",0); setDim(2,"",0); }
        else if (n.type=="圆柱") { setDim(0,"半径 (R):",1); setDim(1,"高度 (H):",1); setDim(2,"",0); }
        else if (n.type=="圆锥") { setDim(0,"底半径 (R):",1); setDim(1,"",0); setDim(2,"高度 (H):",1); }
        else if (n.type=="圆台") { setDim(0,"底半径 (R1):",1); setDim(1,"顶半径 (R2):",1); setDim(2,"高度 (H):",1); }
        else if (n.type=="圆环") { setDim(0,"主半径 (R1):",1); setDim(1,"副半径 (R2):",1); setDim(2,"",0); }
        else if (n.type=="楔形") { setDim(0,"长 (dx):",1); setDim(1,"宽 (dy):",1); setDim(2,"高 (dz):",1); setDim4("X偏移 (ltx):",1); }
        else if (n.type=="椭球") { setDim(0,"X半径 (Rx):",1); setDim(1,"Y半径 (Ry):",1); setDim(2,"Z半径 (Rz):",1); }
        else if (n.type=="棱柱"||n.type=="棱锥"||n.type=="棱台") {
            if (n.type=="棱台") { setDim(0,"底半径 (R1):",1); setDim(1,"顶半径 (R2):",1); setDim(2,"高度 (H):",1); }
            else if (n.type=="棱柱") { setDim(0,"半径 (R):",1); setDim(1,"高度 (H):",1); setDim(2,"",0); }
            else { setDim(0,"底半径 (R):",1); setDim(1,"高度 (H):",1); setDim(2,"",0); }
            setDim4("棱数:",1);
            m_geomE->setRange(3, 64); m_geomE->setDecimals(0);
        } else {
            setDim(0,"长 (W):",1); setDim(1,"宽 (H):",1); setDim(2,"高 (D):",1);
            m_geomE->setRange(-1e6, 1e6); m_geomE->setDecimals(3);
        }
    }
    m_geomE->setValue(n.e);
    m_geomPX->setValue(n.px); m_geomPY->setValue(n.py); m_geomPZ->setValue(n.pz);
    m_geomRotAng->setValue(n.rotAng); m_geomRotAxis->setCurrentIndex(n.rotAxis);
    if (m_chkKeepInteriorU) m_chkKeepInteriorU->setChecked(n.keepInterior);
    if (m_chkKeepInteriorS) m_chkKeepInteriorS->setChecked(n.keepInterior);
    if (m_chkKeepInputsU) m_chkKeepInputsU->setChecked(n.keepInputs);
    if (m_chkKeepInputsS) m_chkKeepInputsS->setChecked(n.keepInputs);
    if (m_chkKeepToolsS) m_chkKeepToolsS->setChecked(n.keepTools);
    m_listBoolUnion->clear(); m_listBoolSubIn->clear(); m_listBoolSubTool->clear();
    // 变换参数加载
    if (n.type=="平移"||n.type=="旋转"||n.type=="镜像"||n.type=="缩放"||n.type=="偏移") {
        auto* pg = m_geomStackPages->currentWidget();
        if (pg) {
            QList<QDoubleSpinBox*> sps = pg->findChildren<QDoubleSpinBox*>();
            QList<QComboBox*> cmbs = pg->findChildren<QComboBox*>();
            QList<QCheckBox*> chks = pg->findChildren<QCheckBox*>();
            if (!chks.isEmpty()) chks[0]->setChecked(n.keepInputs);
            if (n.type=="平移" && sps.size()>=3) { sps[0]->setValue(n.px); sps[1]->setValue(n.py); sps[2]->setValue(n.pz); }
            else if (n.type=="旋转" && cmbs.size()>=1 && sps.size()>=4) {
                cmbs[0]->setCurrentIndex(n.rotAxis); sps[0]->setValue(n.rotAng);
                sps[1]->setValue(n.w); sps[2]->setValue(n.h); sps[3]->setValue(n.d);
            } else if (n.type=="镜像" && cmbs.size()>=1) { cmbs[0]->setCurrentIndex(n.rotAxis); if (sps.size()>=1) sps[0]->setValue(n.d); }
            else if (n.type=="缩放" && sps.size()>=4) { sps[0]->setValue(std::max(0.01,n.w)); sps[1]->setValue(n.px); sps[2]->setValue(n.py); sps[3]->setValue(n.pz); }
            else if (n.type=="偏移" && sps.size()>=1) sps[0]->setValue(n.w);
            QList<QListWidget*> lws = pg->findChildren<QListWidget*>();
            if (!lws.isEmpty()) { lws[0]->clear();
                for (int ni : n.inputIndices) if (ni>=0&&ni<(int)m_geomNodes.size()&&m_geomNodes[ni].treeItem) {
                    auto* li = new QListWidgetItem(m_geomNodes[ni].treeItem->text(0)); li->setData(Qt::UserRole, ni); lws[0]->addItem(li); }
            }
        }
    }
    auto fill = [&](QListWidget* lw, const std::vector<int>& v) {
        for (int i : v) if (i>=0 && i<(int)m_geomNodes.size()) {
            auto* li = new QListWidgetItem(m_geomNodes[i].treeItem->text(0));
            li->setData(Qt::UserRole, i); lw->addItem(li);
        }
    };
    fill(m_listBoolUnion, n.inputIndices);
    fill(m_listBoolSubIn, n.inputIndices);
    fill(m_listBoolSubTool, n.toolIndices);
    // 加载倒圆角/倒斜角边列表
    m_listFilletEdges->clear(); m_listChamfEdges->clear();
    if (n.type == "倒圆角") {
        m_spinFilletR->setValue(n.w);
        for (auto& [ni, ei] : n.edgeSelections)
            if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].treeItem) {
                auto* li = new QListWidgetItem(QString("%1 - 边%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(ei));
                li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, ei);
                m_listFilletEdges->addItem(li);
            }
    } else if (n.type == "倒斜角") {
        m_spinChamfD->setValue(n.w);
        for (auto& [ni, ei] : n.edgeSelections)
            if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].treeItem) {
                auto* li = new QListWidgetItem(QString("%1 - 边%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(ei));
                li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, ei);
                m_listChamfEdges->addItem(li);
            }
    }
    // 加载扫掠列表
    m_listSwProf->clear(); m_listSwPath->clear();
    if (n.type == "扫掠") {
        // 截面: 存储 (ni, localFaceIdx)
        for (int ni : n.inputIndices)
            if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].treeItem) {
                auto* li = new QListWidgetItem(QString("%1 - 面%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(n.sweepProfileFaceIdx));
                li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, n.sweepProfileFaceIdx);
                m_listSwProf->addItem(li);
            }
        // 路径边
        for (auto& [ni, ei] : n.edgeSelections)
            if (ni >= 0 && ni < (int)m_geomNodes.size() && m_geomNodes[ni].treeItem) {
                auto* li = new QListWidgetItem(QString("%1 - 边%2").arg(m_geomNodes[ni].treeItem->text(0)).arg(ei));
                li->setData(Qt::UserRole, ni); li->setData(Qt::UserRole+1, ei);
                m_listSwPath->addItem(li);
            }
    }
    pushParamsToBindings();
}

void MainWindow::saveGeomNodeFromUI(int idx) {
    if (idx < 0 || idx >= (int)m_geomNodes.size()) return;
    auto& n = m_geomNodes[idx];
    // 倒圆角/倒斜角/扫掠的w等参数不应被通用尺寸微调框覆盖
    if (n.type != "倒圆角" && n.type != "倒斜角" && n.type != "扫掠") {
        n.w=m_geomW->value(); n.h=m_geomH->value(); n.d=m_geomD->value(); n.e=m_geomE->value();
    }
    n.px=m_geomPX->value(); n.py=m_geomPY->value(); n.pz=m_geomPZ->value();
    n.rotAng=m_geomRotAng->value(); n.rotAxis=m_geomRotAxis->currentIndex();
    n.keepInterior = (n.type=="差集") ? m_chkKeepInteriorS->isChecked() : m_chkKeepInteriorU->isChecked();
    n.keepInputs = (n.type=="差集") ? m_chkKeepInputsS->isChecked() : m_chkKeepInputsU->isChecked();
    n.keepTools = m_chkKeepToolsS ? m_chkKeepToolsS->isChecked() : true;

    auto read=[&](QListWidget* lw)->std::vector<int>{
        std::vector<int> v; for(int i=0;i<lw->count();i++) v.push_back(lw->item(i)->data(Qt::UserRole).toInt());
        return v;
    };
    if (n.type == "并集" || n.type == "交集") n.inputIndices = read(m_listBoolUnion);
    else if (n.type == "差集") {
        n.inputIndices = read(m_listBoolSubIn);
        n.toolIndices = read(m_listBoolSubTool);
    } else if (n.type=="平移"||n.type=="旋转"||n.type=="镜像"||n.type=="缩放"||n.type=="偏移"
               ||n.type=="倒圆角"||n.type=="倒斜角"||n.type=="扫掠"
               ||n.type=="线性阵列"||n.type=="圆形阵列") {
        QListWidget* lw = nullptr;
        if (n.type=="平移") lw=m_listXfmT; else if (n.type=="旋转") lw=m_listXfmR;
        else if (n.type=="镜像") lw=m_listXfmM; else if (n.type=="缩放") lw=m_listXfmS;
        else if (n.type=="偏移") lw=m_listXfmO;
        else if (n.type=="倒圆角") {
            lw = m_listFilletEdges;
            n.edgeSelections.clear();
            for (int i=0; i<lw->count(); i++) {
                int ni = lw->item(i)->data(Qt::UserRole).toInt();
                int ei = lw->item(i)->data(Qt::UserRole+1).toInt();
                n.edgeSelections.push_back({ni, ei});
            }
        } else if (n.type=="倒斜角") {
            lw = m_listChamfEdges;
            n.edgeSelections.clear();
            for (int i=0; i<lw->count(); i++) {
                int ni = lw->item(i)->data(Qt::UserRole).toInt();
                int ei = lw->item(i)->data(Qt::UserRole+1).toInt();
                n.edgeSelections.push_back({ni, ei});
            }
        }
        else if (n.type=="扫掠") {
            // 截面: inputIndices存节点, sweepProfileFaceIdx存面索引
            if (m_listSwProf->count() > 0) {
                n.inputIndices = {m_listSwProf->item(0)->data(Qt::UserRole).toInt()};
                n.sweepProfileFaceIdx = m_listSwProf->item(0)->data(Qt::UserRole+1).toInt();
            }
            // 路径边存edgeSelections
            n.edgeSelections.clear();
            for (int i=0; i<m_listSwPath->count(); i++) {
                int ni = m_listSwPath->item(i)->data(Qt::UserRole).toInt();
                int ei = m_listSwPath->item(i)->data(Qt::UserRole+1).toInt();
                n.edgeSelections.push_back({ni, ei});
            }
        } else if (n.type=="线性阵列") lw=m_listArrLin;
        else if (n.type=="圆形阵列") lw=m_listArrCir;
        if (lw) n.inputIndices = read(lw);
        // 保存变换特有参数
        auto* pg = m_geomStackPages->currentWidget();
        if (pg) {
            QList<QDoubleSpinBox*> sps = pg->findChildren<QDoubleSpinBox*>();
            QList<QComboBox*> cmbs = pg->findChildren<QComboBox*>();
            QList<QCheckBox*> chks = pg->findChildren<QCheckBox*>();
            if (!chks.isEmpty()) n.keepInputs = chks[0]->isChecked();
            if (n.type=="平移" && sps.size()>=3) { n.px=sps[0]->value(); n.py=sps[1]->value(); n.pz=sps[2]->value(); }
            else if (n.type=="旋转" && sps.size()>=4) { n.rotAxis=(cmbs.size()>=1?cmbs[0]->currentIndex():0); n.rotAng=sps[0]->value(); n.w=sps[1]->value(); n.h=sps[2]->value(); n.d=sps[3]->value(); }
            else if (n.type=="镜像" && cmbs.size()>=1 && sps.size()>=1) { n.rotAxis=cmbs[0]->currentIndex(); n.d=sps[0]->value(); }
            else if (n.type=="缩放" && sps.size()>=4) { n.w=sps[0]->value(); n.px=sps[1]->value(); n.py=sps[2]->value(); n.pz=sps[3]->value(); }
            else if (n.type=="偏移" && sps.size()>=1) { n.w=sps[0]->value(); }
        }
    }
    // 保存节点级参数绑定
    n.nodeParamBindings.clear();
    for (auto& [spin, pname] : m_paramBindings)
        if (m_spinboxNames.count(spin))
            n.nodeParamBindings[m_spinboxNames[spin]] = pname;
}

// ====================================================================
// 构建几何: 保存指定节点的 UI 参数后链式构建 0..idx
// (btnBuild "构建/更新几何体" 与工具栏"全部构建"共用入口)
// ====================================================================
void MainWindow::buildGeometryFromUI(int idx) {
    if (idx < 0 || idx >= (int)m_geomNodes.size()) return;
    // 变换类型: 先让 saveGeomNodeFromUI 保存通用参数, 再覆盖变换特有参数
    std::vector<int> savedInputIndices;
    auto& node = m_geomNodes[idx];
    bool isXfm = (node.type=="平移"||node.type=="旋转"||node.type=="镜像"||node.type=="缩放"||node.type=="偏移"
                   ||node.type=="倒圆角"||node.type=="倒斜角"||node.type=="扫掠"
                   ||node.type=="线性阵列"||node.type=="圆形阵列");
    if (isXfm) {
        QListWidget* lw = nullptr;
        if (node.type=="平移") lw=m_listXfmT; else if (node.type=="旋转") lw=m_listXfmR;
        else if (node.type=="镜像") lw=m_listXfmM; else if (node.type=="缩放") lw=m_listXfmS;
        else if (node.type=="偏移") lw=m_listXfmO;
        else if (node.type=="倒圆角") {
            lw=m_listFilletEdges; node.w = m_spinFilletR->value();
            node.edgeSelections.clear(); node.inputIndices.clear();
            for (int i=0; i<lw->count(); i++) {
                int ni = lw->item(i)->data(Qt::UserRole).toInt();
                int ei = lw->item(i)->data(Qt::UserRole+1).toInt();
                node.edgeSelections.push_back({ni, ei});
                if (std::find(node.inputIndices.begin(),node.inputIndices.end(),ni)==node.inputIndices.end())
                    node.inputIndices.push_back(ni);
            }
        } else if (node.type=="倒斜角") {
            lw=m_listChamfEdges; node.w = m_spinChamfD->value();
            node.edgeSelections.clear(); node.inputIndices.clear();
            for (int i=0; i<lw->count(); i++) {
                int ni = lw->item(i)->data(Qt::UserRole).toInt();
                int ei = lw->item(i)->data(Qt::UserRole+1).toInt();
                node.edgeSelections.push_back({ni, ei});
                if (std::find(node.inputIndices.begin(),node.inputIndices.end(),ni)==node.inputIndices.end())
                    node.inputIndices.push_back(ni);
            }
        }
        else if (node.type=="扫掠") {
            // 扫掠: 截面存inputIndices+面索引, 路径边存edgeSelections
            node.inputIndices.clear(); node.edgeSelections.clear();
            if (m_listSwProf->count() > 0) {
                node.inputIndices.push_back(m_listSwProf->item(0)->data(Qt::UserRole).toInt());
                node.sweepProfileFaceIdx = m_listSwProf->item(0)->data(Qt::UserRole+1).toInt();
            }
            for (int i=0; i<m_listSwPath->count(); i++) {
                int ni = m_listSwPath->item(i)->data(Qt::UserRole).toInt();
                int ei = m_listSwPath->item(i)->data(Qt::UserRole+1).toInt();
                node.edgeSelections.push_back({ni, ei});
            }
        } else if (node.type=="线性阵列") lw=m_listArrLin;
        else if (node.type=="圆形阵列") lw=m_listArrCir;
        if (lw) for (int i=0; i<lw->count(); i++) savedInputIndices.push_back(lw->item(i)->data(Qt::UserRole).toInt());
    }
    saveGeomNodeFromUI(idx);
    // 变换参数必须在 saveGeomNodeFromUI 之后设置 (避免被覆盖)
    if (isXfm) {
        if (node.type != "扫掠") node.inputIndices = savedInputIndices;
        auto* pg = m_geomStackPages->currentWidget();
        QList<QDoubleSpinBox*> sps = pg ? pg->findChildren<QDoubleSpinBox*>() : QList<QDoubleSpinBox*>();
        QList<QComboBox*> cmbs = pg ? pg->findChildren<QComboBox*>() : QList<QComboBox*>();
        QList<QCheckBox*> chks = pg ? pg->findChildren<QCheckBox*>() : QList<QCheckBox*>();
        if (node.type=="平移" && sps.size()>=3) { node.px=sps[0]->value(); node.py=sps[1]->value(); node.pz=sps[2]->value(); }
        else if (node.type=="旋转" && sps.size()>=4) { node.rotAxis=(cmbs.size()>=1?cmbs[0]->currentIndex():0); node.rotAng=sps[0]->value(); node.w=sps[1]->value(); node.h=sps[2]->value(); node.d=sps[3]->value(); }
        else if (node.type=="镜像" && cmbs.size()>=1 && sps.size()>=1) { node.rotAxis=cmbs[0]->currentIndex(); node.d=sps[0]->value(); }
        else if (node.type=="缩放" && sps.size()>=4) { node.w=sps[0]->value(); node.px=sps[1]->value(); node.py=sps[2]->value(); node.pz=sps[3]->value(); }
        else if (node.type=="偏移" && sps.size()>=1) { node.w=sps[0]->value(); }
        if (!chks.isEmpty()) node.keepInputs = chks[0]->isChecked();
        // 新操作参数
        if (node.type=="倒圆角") node.w = m_spinFilletR->value();
        else if (node.type=="倒斜角") node.w = m_spinChamfD->value();
        else if (node.type=="线性阵列") {
            node.rotAxis = m_cmbArrDir->currentIndex();
            node.px = (node.rotAxis==0) ? m_spArrSpc->value() : 0;
            node.py = (node.rotAxis==1) ? m_spArrSpc->value() : 0;
            node.pz = (node.rotAxis==2) ? m_spArrSpc->value() : 0;
            node.e = m_spArrCnt->value();
        } else if (node.type=="圆形阵列") {
            node.rotAxis = m_cmbCirAx->currentIndex();
            node.w = m_spCirCx->value(); node.h = m_spCirCy->value(); node.d = m_spCirCz->value();
            node.rotAng = m_spCirAng->value(); node.e = m_spCirCnt->value();
        }
    }
    m_btnAddUnion->setChecked(false); m_btnAddSubIn->setChecked(false); m_btnAddSubTool->setChecked(false);
    auto* pg = m_geomStackPages->currentWidget();
    if (pg) { QList<QPushButton*> btns = pg->findChildren<QPushButton*>(); for (auto* b : btns) if (b->isCheckable()) b->setChecked(false); }
    fprintf(stderr, "[DEBUG buildBtn] idx=%d type=%s edgeSel=%zu inputIndices=%zu\n",
            idx, qPrintable(node.type), node.edgeSelections.size(), node.inputIndices.size());
    pushParamsToBindings();
    buildGeometrySequence(idx, idx, &node.inputIndices);
    m_projectModified = true;
}

// ====================================================================
// 几何序列构建 (清空全部, 从节点 0 到 upToIndex 顺序执行)
// ====================================================================
void MainWindow::buildGeometrySequence(int upToIndex, int forceIdx,
                                        const std::vector<int>* forceInputs) {
    ui->cadViewerWidget->clearScene();
    for (auto& nd : m_geomNodes) { nd.built = false; nd.builtFaceId = -1; }

    for (int i = 0; i <= upToIndex && i < (int)m_geomNodes.size(); i++) {
        auto& node = m_geomNodes[i];
        QString t = node.type;
        TopoDS_Shape shape;
        bool ok = true;

        auto getShapes = [&](const std::vector<int>& indices) -> std::vector<TopoDS_Shape> {
            std::vector<TopoDS_Shape> v;
            for (int ni : indices) if (ni>=0&&ni<i&&m_geomNodes[ni].built) v.push_back(m_geomNodes[ni].builtShape);
            return v;
        };

        try {
        if (t == "立方体") shape = BRepPrimAPI_MakeBox(node.w, node.h, node.d).Shape();
        else if (t == "球体") shape = BRepPrimAPI_MakeSphere(node.w).Shape();
        else if (t == "圆柱") shape = BRepPrimAPI_MakeCylinder(node.w, node.h).Shape();
        else if (t == "圆锥") shape = BRepPrimAPI_MakeCone(std::max(0.1,node.w), 0, node.d).Shape();
        else if (t == "圆台") {
            double r1 = std::max(0.1, node.w), r2 = std::max(0.0, node.h);
            if (std::abs(r1 - r2) < 0.01) shape = BRepPrimAPI_MakeCylinder(r1, node.d).Shape();
            else shape = BRepPrimAPI_MakeCone(r1, r2, node.d).Shape();
        }
        else if (t == "圆环") shape = BRepPrimAPI_MakeTorus(node.w, node.h).Shape();
        else if (t == "楔形") shape = BRepPrimAPI_MakeWedge(node.w, node.h, node.d, node.e).Shape();
        else if (t == "椭球") {
            // 用矩阵对角线直接设三轴缩放 (SetAffinity 会互相覆盖)
            shape = BRepPrimAPI_MakeSphere(1.0).Shape();
            if (node.w>0 && node.h>0 && node.d>0) {
                gp_GTrsf gtrsf;
                gtrsf.SetValue(1,1, node.w);
                gtrsf.SetValue(2,2, node.h);
                gtrsf.SetValue(3,3, node.d);
                shape = BRepBuilderAPI_GTransform(shape, gtrsf).Shape();
            }
        }
        else if (t == "棱柱" || t == "棱锥" || t == "棱台") {
            int nSides = std::max(3, std::min(12, node.e > 1 ? (int)node.e : 3));
            double r = std::max(0.1, node.w);
            double r2 = std::max(0.0, t=="棱台" ? node.h : 0.0);
            double ht = std::max(0.1, t=="棱台" ? node.d : node.h);
            BRepBuilderAPI_MakePolygon poly;
            for (int s = 0; s < nSides; s++) {
                double a = 2*M_PI*s/nSides;
                poly.Add(gp_Pnt(r*cos(a), r*sin(a), 0));
            }
            poly.Close();
            BRepBuilderAPI_MakeFace face(poly.Wire());
            if (face.IsDone()) {
                if (t == "棱柱") {
                    BRepPrimAPI_MakePrism prism(face.Shape(), gp_Vec(0,0,ht));
                    prism.Build();
                    if (prism.IsDone()) shape = prism.Shape();
                } else if (t == "棱台") {
                    BRepBuilderAPI_MakePolygon polyTop;
                    for (int s = 0; s < nSides; s++) {
                        double a = 2*M_PI*s/nSides;
                        polyTop.Add(gp_Pnt(r2*cos(a), r2*sin(a), ht));
                    }
                    polyTop.Close();
                    BRepOffsetAPI_ThruSections loft(true, true); // 实体模式 + 盖面
                    loft.AddWire(poly.Wire());
                    loft.AddWire(polyTop.Wire());
                    loft.Build();
                    if (loft.IsDone()) shape = loft.Shape();
                    else shape = BRepPrimAPI_MakeCone(r, r2, ht).Shape();
                } else {
                    BRepBuilderAPI_MakeVertex vtx(gp_Pnt(0,0,ht));
                    BRepOffsetAPI_ThruSections loft(true, true);
                    loft.AddWire(poly.Wire());
                    loft.AddVertex(TopoDS::Vertex(vtx.Vertex()));
                    loft.Build();
                    if (loft.IsDone()) shape = loft.Shape();
                    else shape = BRepPrimAPI_MakeCone(r, 0, ht).Shape();
                }
            }
            if (shape.IsNull()) shape = BRepPrimAPI_MakeCylinder(r, ht).Shape();
        }
        else if (t == "平移" || t == "旋转" || t == "镜像" || t == "缩放" || t == "偏移") {
            const auto& indices = (forceInputs && i == forceIdx) ? *forceInputs : node.inputIndices;
            auto inputs = getShapes(indices);
            if (inputs.empty()) { QMessageBox::warning(this,"构建错误","至少需要1个输入对象"); return; }
            gp_Trsf trsf;
            if (t == "平移") trsf.SetTranslation(gp_Vec(node.px, node.py, node.pz));
            else if (t == "旋转") {
                gp_Ax1 ax(gp_Pnt(node.w,node.h,node.d),
                          node.rotAxis==0?gp_Dir(1,0,0):node.rotAxis==1?gp_Dir(0,1,0):gp_Dir(0,0,1));
                trsf.SetRotation(ax, node.rotAng * M_PI / 180.0);
            } else if (t == "镜像") {
                double pos = node.d;
                gp_Pnt o(0,0,0); gp_Dir n(0,0,1);
                if (node.rotAxis==0) { o.SetZ(pos); n=gp_Dir(0,0,1); }
                else if (node.rotAxis==1) { o.SetX(pos); n=gp_Dir(1,0,0); }
                else { o.SetY(pos); n=gp_Dir(0,1,0); }
                trsf.SetMirror(gp_Ax2(o, n));
            } else if (t == "缩放") trsf.SetScale(gp_Pnt(node.px,node.py,node.pz), std::max(0.01, node.w));
            // 偏移: 使用 BRepOffsetAPI_MakeOffsetShape
            if (t == "偏移") {
                for (auto& in : inputs) {
                    BRepOffsetAPI_MakeOffsetShape off;
                    off.PerformByJoin(in, node.w, 0.001, BRepOffset_Skin, false, false);
                    if (off.IsDone()) { shape = off.Shape(); break; }
                }
            } else {
                // 平移/旋转/镜像/缩放: BRepBuilderAPI_Transform
                for (auto& in : inputs) {
                    BRepBuilderAPI_Transform trans(in, trsf, true); // copy
                    if (trans.IsDone()) { shape = trans.Shape(); break; }
                }
            }
        }
        else if (t == "并集" || t == "交集" || t == "差集") {
            auto inputs = getShapes(node.inputIndices);
            if (t == "并集" || t == "交集") {
                if (inputs.size() < 2) { QMessageBox::warning(this,"构建错误",QString("节点 '%1' 至少需要2个输入对象").arg(node.treeItem->text(0))); return; }
                else {
                    if (node.keepInterior) {
                        // 非破坏性融合: 保留内部交界面
                        BOPAlgo_Builder builder;
                        for (auto& s : inputs) builder.AddArgument(s);
                        builder.SetNonDestructive(true);
                        builder.Perform();
                        if (!builder.HasErrors()) shape = builder.Shape();
                        else { /* fallback */ TopoDS_Compound comp; BRep_Builder B; B.MakeCompound(comp); for (auto& s : inputs) B.Add(comp, s); shape = comp; }
                    } else {
                        shape = inputs[0];
                        for (size_t j=1; j<inputs.size(); j++) {
                            if (t=="并集") { BRepAlgoAPI_Fuse f(shape,inputs[j]); if(f.IsDone()) shape=f.Shape(); }
                            else { BRepAlgoAPI_Common c(shape,inputs[j]); if(c.IsDone()) shape=c.Shape(); }
                        }
                        // 清理融合后的残留边/顶点
                        ShapeUpgrade_UnifySameDomain unify(shape, true, true, true);
                        unify.Build();
                        shape = unify.Shape();
                    }
                }
            } else {
                auto tools = getShapes(node.toolIndices);
                if (inputs.empty()||tools.empty()) { QMessageBox::warning(this,"构建错误",QString("节点 '%1' 需要输入对象和减去对象").arg(node.treeItem->text(0))); return; }
                else {
                    shape = inputs[0];
                    for (size_t j=1; j<inputs.size(); j++) { BRepAlgoAPI_Fuse f(shape,inputs[j]); if(f.IsDone()) shape=f.Shape(); }
                    for (auto& tk : tools) { BRepAlgoAPI_Cut c(shape,tk); if(c.IsDone()) shape=c.Shape(); }
                }
            }
        }
        else if (t == "倒圆角" || t == "倒斜角") {
            fprintf(stderr, "[DEBUG buildGeomSeq] i=%d type=%s inputIndices.size()=%zu edgeSelections.size()=%zu node.w=%f\n",
                    i, qPrintable(t), node.inputIndices.size(), node.edgeSelections.size(), node.w);
            if (node.inputIndices.empty()) { QMessageBox::warning(this,"构建错误","至少需要选择一个边"); return; }
            shape = m_geomNodes[node.inputIndices[0]].builtShape;
            fprintf(stderr, "[DEBUG buildGeomSeq] input[0]=%d, shape.IsNull()=%d\n", node.inputIndices[0], shape.IsNull());
            // 修复基体几何, 确保倒角算法获得干净的边, 避免产生褶皱B样条曲面
            ShapeFix_Shape preFix(shape);
            preFix.SetPrecision(1e-5);
            preFix.Perform();
            shape = preFix.Shape();
            if (t == "倒圆角") {
                BRepFilletAPI_MakeFillet fillet(shape);
                if (node.edgeSelections.empty()) {
                    for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next())
                        fillet.Add(std::max(0.01,node.w), TopoDS::Edge(ex.Current()));
                } else {
                    for (auto& [ni, ei] : node.edgeSelections) {
                        int cnt=0;
                        for (TopExp_Explorer ex(m_geomNodes[ni].builtShape, TopAbs_EDGE); ex.More(); ex.Next(), cnt++)
                            if (cnt==ei) { fillet.Add(std::max(0.01,node.w), TopoDS::Edge(ex.Current())); break; }
                    }
                }
                fillet.Build(); fprintf(stderr, "[DEBUG buildGeomSeq] fillet.IsDone=%d\n", (int)fillet.IsDone()); if (fillet.IsDone()) shape = fillet.Shape();
                else { QMessageBox::warning(this, "构建错误", QString("倒圆角失败: 半径 %1 可能过大, 请减小半径后重试").arg(node.w)); return; }
            } else {
                BRepFilletAPI_MakeChamfer chamfer(shape);
                auto addChamferEdge = [&](const TopoDS_Shape& srcShape, const TopoDS_Edge& E) {
                    // 找到与E相邻的面作为参考面
                    for (TopExp_Explorer fx(srcShape, TopAbs_FACE); fx.More(); fx.Next()) {
                        TopoDS_Face testFace = TopoDS::Face(fx.Current());
                        for (TopExp_Explorer ex2(testFace, TopAbs_EDGE); ex2.More(); ex2.Next())
                            if (ex2.Current().IsSame(E)) {
                                chamfer.Add(std::max(0.01,node.w), std::max(0.01,node.w), E, testFace);
                                return;
                            }
                    }
                };
                if (node.edgeSelections.empty()) {
                    for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next())
                        addChamferEdge(shape, TopoDS::Edge(ex.Current()));
                } else {
                    for (auto& [ni, ei] : node.edgeSelections) {
                        int cnt=0;
                        for (TopExp_Explorer ex(m_geomNodes[ni].builtShape, TopAbs_EDGE); ex.More(); ex.Next(), cnt++)
                            if (cnt==ei) { addChamferEdge(m_geomNodes[ni].builtShape, TopoDS::Edge(ex.Current())); break; }
                    }
                }
                chamfer.Build(); fprintf(stderr, "[DEBUG buildGeomSeq] chamfer.IsDone=%d\n", (int)chamfer.IsDone()); if (chamfer.IsDone()) shape = chamfer.Shape();
                else { QMessageBox::warning(this, "构建错误", QString("倒斜角失败: 距离 %1 可能过大, 请减小距离后重试").arg(node.w)); return; }
            }
        }
        else if (t == "扫掠") {
            if (node.inputIndices.empty() || node.edgeSelections.empty()) { QMessageBox::warning(this,"构建错误","扫掠需要截面和路径边"); return; }
            auto profiles = getShapes(node.inputIndices);
            if (profiles.empty()) { QMessageBox::warning(this,"构建错误","请选择截面对象"); return; }
            // 截面: 取指定的面 (由sweepProfileFaceIdx索引)
            TopoDS_Face profileFace;
            { int cnt=0;
              for (TopExp_Explorer ex(profiles[0], TopAbs_FACE); ex.More(); ex.Next(), cnt++)
                  if (cnt == node.sweepProfileFaceIdx) { profileFace = TopoDS::Face(ex.Current()); break; }
            }
            if (profileFace.IsNull()) { QMessageBox::warning(this,"构建错误","截面对象无有效面"); return; }
            // 路径: 从选中的边构建wire
            BRepBuilderAPI_MakeWire mkWire;
            for (auto& [ni, ei] : node.edgeSelections) {
                if (ni >= 0 && ni < i && m_geomNodes[ni].built) {
                    int cnt=0;
                    for (TopExp_Explorer ex(m_geomNodes[ni].builtShape, TopAbs_EDGE); ex.More(); ex.Next(), cnt++)
                        if (cnt==ei) { mkWire.Add(TopoDS::Edge(ex.Current())); break; }
                }
            }
            TopoDS_Wire pathWire;
            if (mkWire.IsDone()) pathWire = mkWire.Wire();
            if (pathWire.IsNull()) { QMessageBox::warning(this,"构建错误","无法从选中的边构建路径"); return; }
            BRepOffsetAPI_MakePipe pipe(pathWire, profileFace);
            pipe.Build();
            if (pipe.IsDone()) shape = pipe.Shape();
            else { QMessageBox::warning(this,"构建错误","扫掠失败, 请检查:\n1. 截面必须是平面(如圆柱的顶/底面, 不能是侧面)\n2. 路径边需从截面附近出发"); return; }
        }
        else if (t == "线性阵列" || t == "圆形阵列") {
            auto inputs = getShapes(node.inputIndices);
            if (inputs.empty()) { QMessageBox::warning(this,"构建错误","至少需要1个输入对象"); return; }
            int count = std::max(2, std::min(100, node.e > 1 ? (int)node.e : 2));
            TopoDS_Compound comp; BRep_Builder B; B.MakeCompound(comp);
            for (auto& in : inputs) B.Add(comp, in);  // 保留原始
            for (auto& in : inputs) {
                gp_Trsf trsf;
                for (int k = 1; k < count; k++) {
                    if (t == "线性阵列")
                        trsf.SetTranslation(gp_Vec(node.px*k, node.py*k, node.pz*k));
                    else {
                        gp_Ax1 ax(gp_Pnt(node.w,node.h,node.d),
                                  node.rotAxis==0?gp_Dir(1,0,0):node.rotAxis==1?gp_Dir(0,1,0):gp_Dir(0,0,1));
                        trsf.SetRotation(ax, node.rotAng * k * M_PI / 180.0);
                    }
                    BRepBuilderAPI_Transform trans(in, trsf, true);
                    if (trans.IsDone()) B.Add(comp, trans.Shape());
                }
            }
            shape = comp;
        }
        else ok = false;

        if (!ok) continue;
        } catch (const std::exception& ex) {
            qWarning("buildGeometrySequence: node %d (%s) error: %s", i, qPrintable(t), ex.what());
            QMessageBox::warning(this, "构建错误", QString("节点 '%1' 构建失败: %2").arg(node.treeItem->text(0)).arg(ex.what()));
            continue;
        } catch (const Standard_Failure& ex) {
            qWarning("buildGeometrySequence: node %d (%s) OCC error: %s", i, qPrintable(t), ex.GetMessageString());
            QMessageBox::warning(this, "构建错误", QString("节点 '%1' OCC构建失败: %2").arg(node.treeItem->text(0)).arg(ex.GetMessageString()));
            continue;
        }

        // 旋转+平移 (仅对体素, 变换节点已自行处理)
        bool isXfmNode = (t=="平移"||t=="旋转"||t=="镜像"||t=="缩放"||t=="偏移");
        if (!isXfmNode) {
            if (std::abs(node.rotAng) > 0.01) {
                gp_Ax1 axis(gp_Pnt(0,0,0), node.rotAxis==0 ? gp_Dir(1,0,0) :
                                           node.rotAxis==1 ? gp_Dir(0,1,0) : gp_Dir(0,0,1));
                gp_Trsf rot; rot.SetRotation(axis, node.rotAng * M_PI / 180.0);
                shape.Move(TopLoc_Location(rot));
            }
            if (node.px!=0||node.py!=0||node.pz!=0) {
                gp_Trsf trsf; trsf.SetTranslation(gp_Vec(node.px,node.py,node.pz));
                shape.Move(TopLoc_Location(trsf));
            }
        }

        // 表面模式: 从实体 Solid 提取外壳 Shell (布尔操作前再转回 Solid)
        if (!node.asSolid && shape.ShapeType() == TopAbs_SOLID) {
            TopExp_Explorer ex(shape, TopAbs_SHELL);
            if (ex.More()) shape = TopoDS::Shell(ex.Current());
        }
        node.builtShape = shape;
        fprintf(stderr, "[DEBUG buildGeomSeq] node %d type=%s calling addOCCShape, shape.IsNull=%d\n", i, qPrintable(t), shape.IsNull());
        node.builtFaceId = ui->cadViewerWidget->addOCCShape(shape);
        fprintf(stderr, "[DEBUG buildGeomSeq] node %d builtFaceId=%d\n", i, node.builtFaceId);
        node.treeItem->setData(0, Qt::UserRole+2, node.builtFaceId);  // faceId 存 UserRole+2, 不覆写类型
        node.treeItem->setForeground(0, QBrush(QColor(0,0,0)));
        node.built = true;
    }

    // 处理布尔/变换节点的 keepInputs: 不保留时隐藏输入对象
    // 倒圆角/倒斜角/扫掠的结果始终替代原始几何体, 因此总是清除输入
    for (int i = 0; i <= upToIndex && i < (int)m_geomNodes.size(); i++) {
        auto& node = m_geomNodes[i];
        if (!node.built) continue;
        bool isReplaceOp = (node.type == "倒圆角" || node.type == "倒斜角");
        if (!node.keepInputs || isReplaceOp) {
            for (int ni : node.inputIndices)
                if (ni >= 0 && ni < i && m_geomNodes[ni].built && m_geomNodes[ni].builtFaceId > 0)
                    ui->cadViewerWidget->clearFaceGroup(m_geomNodes[ni].builtFaceId);
        }
        if (!node.keepTools || isReplaceOp) {
            for (int ni : node.toolIndices)
                if (ni >= 0 && ni < i && m_geomNodes[ni].built && m_geomNodes[ni].builtFaceId > 0)
                    ui->cadViewerWidget->clearFaceGroup(m_geomNodes[ni].builtFaceId);
        }
    }

    // 标记后续节点为未构建 (灰色)
    for (int i = upToIndex + 1; i < (int)m_geomNodes.size(); i++) {
        m_geomNodes[i].built = false;
        if (m_geomNodes[i].treeItem)
            m_geomNodes[i].treeItem->setForeground(0, QBrush(QColor(128,128,128)));
    }

    // 导出所有已构建几何体为STEP文件, 供网格生成使用
    TopoDS_Compound comp; BRep_Builder B; B.MakeCompound(comp);
    TopoDS_Shape singleShape;
    int exportedCount = 0;
    const auto& fom = ui->cadViewerWidget->faceObjectMap();
    for (int i = 0; i <= upToIndex && i < (int)m_geomNodes.size(); i++) {
        auto& nd = m_geomNodes[i];
        if (!nd.built || nd.builtShape.IsNull()) continue;
        // 只导出viewer中仍可见的节点(未被clearFaceGroup清除的)
        bool visible = false;
        for (auto& [fid, oid] : fom)
            if (oid == nd.builtFaceId) { visible = true; break; }
        if (visible) { B.Add(comp, nd.builtShape); singleShape = nd.builtShape; exportedCount++; }
    }
    if (exportedCount > 0) {
        TopoDS_Shape shapeToExport = (exportedCount == 1) ? singleShape : comp;
        // 修复布尔/倒角产生的B样条曲面容差问题, 避免Gmsh网格生成崩溃
        ShapeFix_Shape fixer(shapeToExport);
        fixer.SetPrecision(1e-5);
        fixer.SetMaxTolerance(1e-3);
        fixer.Perform();
        shapeToExport = fixer.Shape();
        // 修复倒角B样条曲面: 降阶 + 删除退化边 + 归一化参数域
        {
            TopTools_IndexedMapOfShape fmap;
            TopExp::MapShapes(shapeToExport, TopAbs_FACE, fmap);
            BRep_Builder B2;
            int modified = 0;
            for (int fi = 1; fi <= fmap.Extent(); fi++) {
                TopoDS_Face face = TopoDS::Face(fmap(fi));
                Handle(Geom_Surface) surf = BRep_Tool::Surface(face);
                Handle(Geom_BSplineSurface) bs = Handle(Geom_BSplineSurface)::DownCast(surf);
                if (bs.IsNull()) continue;
                double faceTol = BRep_Tool::Tolerance(face);
                double ukr = bs->UKnot(bs->NbUKnots()) - bs->UKnot(1);
                double vkr = bs->VKnot(bs->NbVKnots()) - bs->VKnot(1);
                bool needReduce = (bs->VDegree() > 6);
                bool needNormalize = (ukr > 10.0 || vkr > 10.0);
                if (!needReduce && !needNormalize) continue;
                // 先修复退化边(零长度边/NULL曲线)
                ShapeFix_Face ff(face);
                ff.FixSmallAreaWireMode() = Standard_True;
                ff.Perform();
                face = TopoDS::Face(ff.Result());
                surf = BRep_Tool::Surface(face);
                bs = Handle(Geom_BSplineSurface)::DownCast(surf);
                if (bs.IsNull()) continue;
                // 降阶 + 归一化参数域
                double approxTol = std::max(1e-4, faceTol * 20);
                int maxDeg = needReduce ? 6 : std::max(bs->UDegree(), bs->VDegree());
                GeomConvert_ApproxSurface approx(bs, approxTol, GeomAbs_C1, GeomAbs_C1, maxDeg, maxDeg, 20, 1);
                if (!approx.HasResult()) continue;
                Handle(Geom_BSplineSurface) newSurf = approx.Surface();
                double newTol = std::max(faceTol, approx.MaxError() * 2);
                B2.UpdateFace(face, newSurf, TopLoc_Location(), newTol);
                modified++;
            }
            if (modified > 0) {
                BRepLib::SameParameter(shapeToExport, 1e-4, true);
                ShapeFix_Shape fix2(shapeToExport);
                fix2.SetPrecision(1e-4);
                fix2.SetMaxTolerance(1e-3);
                fix2.Perform();
                shapeToExport = fix2.Shape();
            }
        }
        // 合并C1连续的邻面, 减少冗余B样条曲面
        ShapeUpgrade_UnifySameDomain unifier(shapeToExport, true, true, true);
        unifier.Build();
        shapeToExport = unifier.Shape();
        QString stepPath = QDir::tempPath() + "/IonFlow_geometry_cache.brep";
        BRepTools::Write(shapeToExport, stepPath.toUtf8().constData());
        m_currentStepFile = stepPath;
    }

    ui->cadViewerWidget->GetRenderer()->GetRenderWindow()->Render();
}
void MainWindow::saveProject(const QString& _path) {
    // 提交UI中未保存的编辑
    if (m_paramTabs && m_paramTabs->currentWidget()) {
        auto* table = qobject_cast<QTableWidget*>(m_paramTabs->currentWidget());
        if (table) table->setCurrentCell(-1, -1);
    }
    onBoundaryParamsChanged();
    // 从UI读取求解器/物理场最新值到m_config
    if (m_spinE0) m_config.physics.E_onset = m_spinE0->value();
    if (m_spinRho) m_config.physics.rho_surface = m_spinRho->value();
    if (m_spinK) m_config.physics.K_mobility = m_spinK->value();
    if (m_spinWindX) m_config.physics.wind_x = m_spinWindX->value();
    if (m_spinWindY) m_config.physics.wind_y = m_spinWindY->value();
    if (m_spinWindZ) m_config.physics.wind_z = m_spinWindZ->value();
    if (m_spinWRho) m_config.solver.w_rho = m_spinWRho->value();
    if (m_spinGoalConv) m_config.solver.goal_convergence = m_spinGoalConv->value();
    if (m_spinTolE) m_config.solver.tolerance_E = m_spinTolE->value();
    if (m_spinTolRho) m_config.solver.tolerance_rho = m_spinTolRho->value();
    if (m_spinCores) m_config.solver.num_cores = m_spinCores->value();
    // 强制 .ion 扩展名
    QString path = _path;
    if (!path.endsWith(".ion", Qt::CaseInsensitive)) path += ".ion";

    json j;
    j["version"] = 1;
    j["step_file"] = m_currentStepFile.toStdString();
    j["project_dir"] = QFileInfo(path).absoluteDir().absolutePath().toStdString();

    // 网格方案
    json meshes = json::array();
    for (int i = 0; i < m_rootMeshNode->childCount(); i++) {
        auto* node = m_rootMeshNode->child(i);
        json m;
        m["name"] = node->text(0).toStdString();
        if (m_meshFiles.count(node)) m["file"] = m_meshFiles[node].toStdString();
        if (m_isImportedMesh.count(node)) m["imported"] = m_isImportedMesh[node];
        meshes.push_back(m);
    }
    j["mesh_schemes"] = meshes;

    // 当前活跃网格的物理组(可能还未存入m_importedGroupMeta)
    if (m_activeMeshNode && m_rootMeshNode->indexOfChild(m_activeMeshNode) >= 0) {
        std::vector<GroupMeta> meta;
        for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
            auto* ci = m_rootGroupNode->child(i);
            GroupMeta gm{ci->text(0), m_groupTags.count(ci)?m_groupTags[ci]:0, {}};
            if (m_groupData.count(ci)) gm.faceIds = m_groupData[ci];
            meta.push_back(gm);
        }
        if (!meta.empty()) m_importedGroupMeta[m_activeMeshNode] = meta;
    }

    // 所有网格的物理组元数据
    json impGroups = json::object();
    for (auto& [mnode, meta] : m_importedGroupMeta) {
        if (!mnode) continue;
        bool valid = false;
        for (int ci = 0; ci < m_rootMeshNode->childCount(); ci++)
            if (m_rootMeshNode->child(ci) == mnode) { valid = true; break; }
        if (!valid) continue;
        json mg;
        json grps = json::array();
        for (auto& gm : meta) {
            json g;
            g["name"] = gm.name.toStdString(); g["tag"] = gm.tag;
            json fids = json::array(); for (int fid : gm.faceIds) fids.push_back(fid);
            g["face_ids"] = fids;
            grps.push_back(g);
        }
        mg["groups"] = grps;
        if (m_importedBoundaryConfigs.count(mnode)) {
            json bds = json::object();
            for (auto& [name, bp] : m_importedBoundaryConfigs[mnode]) {
                json b;
                b["apply"] = bp.apply; b["use_function"] = bp.useFunction;
                b["expression"] = bp.expression; b["voltage"] = bp.voltage; b["is_corona"] = bp.is_corona;
                bds[name] = b;
            }
            mg["boundary_configs"] = bds;
        }
        impGroups[mnode->text(0).toStdString()] = mg;
    }
    j["imported_group_meta"] = impGroups;

    // 物理组
    json groups = json::array();
    for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
        auto* g = m_rootGroupNode->child(i);
        json grp;
        grp["name"] = g->text(0).toStdString();
        if (m_groupTags.count(g)) grp["gmsh_tag"] = m_groupTags[g];
        json faceIds = json::array();
        if (m_groupData.count(g)) for (int fid : m_groupData[g]) faceIds.push_back(fid);
        grp["face_ids"] = faceIds;
        groups.push_back(grp);
    }
    j["physical_groups"] = groups;

    // 几何树节点(完整保存)
    json geomNodes = json::array();
    for (auto& nd : m_geomNodes) {
        json gn;
        gn["type"] = nd.type.toStdString();
        gn["name"] = nd.treeItem ? nd.treeItem->text(0).toStdString() : "";
        gn["w"]=nd.w; gn["h"]=nd.h; gn["d"]=nd.d; gn["e"]=nd.e;
        gn["px"]=nd.px; gn["py"]=nd.py; gn["pz"]=nd.pz;
        gn["rotAng"]=nd.rotAng; gn["rotAxis"]=nd.rotAxis;
        gn["asSolid"]=nd.asSolid;
        gn["keepInterior"]=nd.keepInterior; gn["keepInputs"]=nd.keepInputs; gn["keepTools"]=nd.keepTools;
        json inIdx = json::array(); for (int v : nd.inputIndices) inIdx.push_back(v); gn["inputIndices"] = inIdx;
        json tIdx = json::array(); for (int v : nd.toolIndices) tIdx.push_back(v); gn["toolIndices"] = tIdx;
        json es = json::array(); for (auto& [ni,ei] : nd.edgeSelections) es.push_back({ni, ei}); gn["edgeSelections"] = es;
        gn["sweepProfileFaceIdx"] = nd.sweepProfileFaceIdx;
        json pb = json::object();
        for (auto& [sn, pn] : nd.nodeParamBindings) pb[sn] = pn;
        gn["param_bindings"] = pb;
        geomNodes.push_back(gn);
    }
    j["geometry_nodes"] = geomNodes;

    // 求解器/物理场参数
    json solverCfg;
    solverCfg["E_onset"] = m_config.physics.E_onset;
    solverCfg["rho_surface"] = m_config.physics.rho_surface;
    solverCfg["K_mobility"] = m_config.physics.K_mobility;
    solverCfg["wind_x"] = m_config.physics.wind_x;
    solverCfg["wind_y"] = m_config.physics.wind_y;
    solverCfg["wind_z"] = m_config.physics.wind_z;
    solverCfg["w_rho"] = m_config.solver.w_rho;
    solverCfg["goal_convergence"] = m_config.solver.goal_convergence;
    solverCfg["max_update_times"] = m_config.solver.max_update_times;
    solverCfg["tolerance_E"] = m_config.solver.tolerance_E;
    solverCfg["tolerance_rho"] = m_config.solver.tolerance_rho;
    solverCfg["order"] = m_config.solver.order;
    solverCfg["num_cores"] = m_config.solver.num_cores;
    solverCfg["solver_type"] = m_config.solver.solver_type;
    j["solver_config"] = solverCfg;

    // 边界条件(含函数字段)
    json bdrs = json::array();
    for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
        auto* g = m_rootGroupNode->child(i);
        if (m_boundaryConfigs.count(g)) {
            auto& bc = m_boundaryConfigs[g];
            json b;
            b["name"] = g->text(0).toStdString();
            b["apply"] = bc.apply;
            b["use_function"] = bc.useFunction;
            b["expression"] = bc.expression;
            b["voltage"] = bc.voltage;
            b["is_corona"] = bc.is_corona;
            bdrs.push_back(b);
        }
    }
    j["boundary_configs"] = bdrs;

    // 网格参数
    json meshCfgs = json::array();
    for (int mi = 0; mi < m_rootMeshNode->childCount(); mi++) {
        auto* mnode = m_rootMeshNode->child(mi);
        if (m_meshConfigs.count(mnode)) {
            for (auto& [gnode, params] : m_meshConfigs[mnode]) {
                if (!gnode) continue;
                bool gv = false;
                for (int gi = 0; gi < m_rootGroupNode->childCount(); gi++)
                    if (m_rootGroupNode->child(gi) == gnode) { gv = true; break; }
                if (!gv) continue;
                json mc;
                mc["mesh_name"] = mnode->text(0).toStdString();
                mc["group_name"] = gnode->text(0).toStdString();
                mc["enableLocalField"] = params.enableLocalField;
                mc["sizeMin"] = params.sizeMin; mc["sizeMax"] = params.sizeMax;
                mc["distMin"] = params.distMin; mc["distMax"] = params.distMax;
                meshCfgs.push_back(mc);
            }
        }
    }
    j["mesh_params"] = meshCfgs;

    // 参数页
    json paramPages = json::array();
    for (auto& pg : m_config.paramPages) {
        json jp; jp["name"] = pg.name;
        json entries = json::array();
        for (auto& e : pg.entries) {
            json je;
            je["name"] = e.name; je["expression"] = e.expression;
            je["value"] = e.value; je["unit"] = e.unit; je["description"] = e.description;
            entries.push_back(je);
        }
        jp["entries"] = entries; paramPages.push_back(jp);
    }
    j["param_pages"] = paramPages;

    // 参数绑定
    json bindings = json::array();
    for (auto& [spin, pname] : m_paramBindings) {
        if (m_spinboxNames.count(spin)) {
            json b; b["spinbox"] = m_spinboxNames[spin]; b["param"] = pname;
            bindings.push_back(b);
        }
    }
    j["param_bindings"] = bindings;

    // 后处理状态
    j["last_field"] = m_lastFieldName.toStdString();
    j["last_vtk"] = m_lastVtkFilePath.toStdString();
    json post;
    post["surface"] = m_radioSurface ? m_radioSurface->isChecked() : true;
    post["slice_x"] = m_chkSliceX ? m_chkSliceX->isChecked() : false;
    post["slice_x_pos"] = m_spinSliceX ? m_spinSliceX->value() : 0.0;
    post["slice_y"] = m_chkSliceY ? m_chkSliceY->isChecked() : false;
    post["slice_y_pos"] = m_spinSliceY ? m_spinSliceY->value() : 0.0;
    post["slice_z"] = m_chkSliceZ ? m_chkSliceZ->isChecked() : false;
    post["slice_z_pos"] = m_spinSliceZ ? m_spinSliceZ->value() : 0.0;
    post["auto_range"] = m_chkAutoRange ? m_chkAutoRange->isChecked() : true;
    post["range_min"] = m_spinRangeMin ? m_spinRangeMin->value() : 0.0;
    post["range_max"] = m_spinRangeMax ? m_spinRangeMax->value() : 0.0;
    post["selected_groups"] = json::array();
    if (m_listPostGroups) for (int i=0; i<m_listPostGroups->count(); i++)
        if (m_listPostGroups->item(i)->isSelected())
            post["selected_groups"].push_back(m_listPostGroups->item(i)->text().toStdString());
    j["post_process"] = post;

    std::ofstream o(path.toStdString());
    o << std::setw(2) << j << std::endl;
    m_projectModified = false;
    m_projectFilePath = path;
    setWindowTitle(QString("离子流场模拟器 - %1").arg(QFileInfo(path).fileName()));
    logMessage("项目已保存: " + path, "success");
}

void MainWindow::loadProject(const QString& path) {
    std::ifstream f(path.toStdString());
    if (!f.is_open()) { QMessageBox::warning(this, "错误", "无法打开项目文件"); return; }
    json j = json::parse(f, nullptr, false);
    if (j.is_discarded()) { QMessageBox::warning(this, "错误", "项目文件格式无效"); return; }

    m_projectFilePath = path;
    m_projectModified = false;

    // 恢复几何树
    m_geomNodes.clear();
    if (j.contains("geometry_nodes")) {
        for (auto& gn : j["geometry_nodes"]) {
            GeomNode node;
            node.type = QString::fromStdString(gn.value("type", ""));
            node.w=gn.value("w",100.0); node.h=gn.value("h",100.0); node.d=gn.value("d",100.0); node.e=gn.value("e",0.0);
            node.px=gn.value("px",0.0); node.py=gn.value("py",0.0); node.pz=gn.value("pz",0.0);
            node.rotAng=gn.value("rotAng",0.0); node.rotAxis=gn.value("rotAxis",0);
            node.asSolid=gn.value("asSolid",true);
            node.keepInterior=gn.value("keepInterior",false);
            node.keepInputs=gn.value("keepInputs",true); node.keepTools=gn.value("keepTools",true);
            for (auto& v : gn.value("inputIndices", json::array())) node.inputIndices.push_back(v.get<int>());
            for (auto& v : gn.value("toolIndices", json::array())) node.toolIndices.push_back(v.get<int>());
            for (auto& es : gn.value("edgeSelections", json::array()))
                node.edgeSelections.push_back({es[0].get<int>(), es[1].get<int>()});
            node.sweepProfileFaceIdx = gn.value("sweepProfileFaceIdx", 0);
            if (gn.contains("param_bindings"))
                for (auto& [sn, pn] : gn["param_bindings"].items())
                    node.nodeParamBindings[sn] = pn.get<std::string>();
            QString name = QString::fromStdString(gn.value("name", node.type.toStdString()));
            auto* item = new QTreeWidgetItem(m_rootGeomNode);
            item->setText(0, name);
            item->setIcon(0, geomIcon(node.type));
            item->setData(0, Qt::UserRole, (int)m_geomNodes.size());
            item->setData(0, Qt::UserRole+1, node.type);
            item->setForeground(0, QBrush(QColor(128,128,128)));
            node.treeItem = item;
            m_geomNodes.push_back(node);
        }
        m_rootGeomNode->setExpanded(true);
        // 重建几何体
        if (!m_geomNodes.empty())
            buildGeometrySequence((int)m_geomNodes.size() - 1);
    }

    // 恢复STEP文件(外部导入)
    if (j.contains("step_file") && !j["step_file"].get<std::string>().empty() && m_geomNodes.empty()) {
        QString stepPath = QString::fromStdString(j["step_file"].get<std::string>());
        if (QFileInfo::exists(stepPath)) {
            m_currentStepFile = stepPath;
            ui->cadViewerWidget->loadStepFile(stepPath);
        }
    }
    // 恢复网格
    if (j.contains("mesh_schemes")) {
        for (auto& mj : j["mesh_schemes"]) {
            QString name = QString::fromStdString(mj.value("name", ""));
            QString file = QString::fromStdString(mj.value("file", ""));
            bool imported = mj.value("imported", false);
            QTreeWidgetItem* node = new QTreeWidgetItem(m_rootMeshNode);
            node->setText(0, name);
            if (!file.isEmpty() && QFileInfo::exists(file)) {
                m_meshFiles[node] = file;
                m_isImportedMesh[node] = imported;
            }
            m_rootMeshNode->setExpanded(true);
        }
    }
    // 恢复导入网格的物理组元数据
    m_importedGroupMeta.clear(); m_importedBoundaryConfigs.clear();
    if (j.contains("imported_group_meta")) {
        for (auto& [nodeName, mg] : j["imported_group_meta"].items()) {
            for (int i = 0; i < m_rootMeshNode->childCount(); i++) {
                auto* mn = m_rootMeshNode->child(i);
                if (mn->text(0).toStdString() == nodeName) {
                    std::vector<GroupMeta> meta;
                    for (auto& gj : mg["groups"]) {
                        GroupMeta gm;
                        gm.name = QString::fromStdString(gj.value("name", "")); gm.tag = gj.value("tag", 0);
                        for (auto& fid : gj.value("face_ids", json::array())) gm.faceIds.insert(fid.get<int>());
                        meta.push_back(gm);
                    }
                    m_importedGroupMeta[mn] = meta;
                    if (mg.contains("boundary_configs")) {
                        std::map<std::string, BoundaryParams> bdMap;
                        for (auto& [name, bj] : mg["boundary_configs"].items()) {
                            BoundaryParams bp;
                            bp.apply = bj.value("apply", false); bp.useFunction = bj.value("use_function", false);
                            bp.expression = bj.value("expression", ""); bp.voltage = bj.value("voltage", 0.0);
                            bp.is_corona = bj.value("is_corona", false);
                            bdMap[name] = bp;
                        }
                        m_importedBoundaryConfigs[mn] = bdMap;
                    }
                    break;
                }
            }
        }
    }

    // 恢复物理组
    m_groupData.clear(); m_groupTags.clear(); m_boundaryConfigs.clear();
    for (int i = m_rootGroupNode->childCount()-1; i >= 0; i--)
        delete m_rootGroupNode->child(i);
    if (j.contains("physical_groups")) {
        std::vector<QTreeWidgetItem*> groupItems;
        for (auto& gj : j["physical_groups"]) {
            auto* g = new QTreeWidgetItem(m_rootGroupNode);
            g->setText(0, QString::fromStdString(gj.value("name", "")));
            m_groupData[g] = std::set<int>();
            if (gj.contains("gmsh_tag")) m_groupTags[g] = gj["gmsh_tag"].get<int>();
            if (gj.contains("face_ids"))
                for (auto& fid : gj["face_ids"]) m_groupData[g].insert(fid.get<int>());
            groupItems.push_back(g);
        }
        if (!m_rootMeshNode->childCount()) { /* no mesh node yet */ }
        else if (m_importedGroupMeta.count(m_rootMeshNode->child(0)) == 0) {
            for (auto* gi : groupItems)
                m_importedGroupMeta[m_rootMeshNode->child(0)].push_back(
                    {gi->text(0), m_groupTags.count(gi) ? m_groupTags[gi] : 0});
        }
    }
    // 恢复边界配置
    if (j.contains("boundary_configs")) {
        for (auto& bj : j["boundary_configs"]) {
            QString name = QString::fromStdString(bj.value("name", ""));
            for (int i = 0; i < m_rootGroupNode->childCount(); i++) {
                auto* g = m_rootGroupNode->child(i);
                if (g->text(0) == name) {
                    BoundaryParams bp;
                    bp.apply = bj.value("apply", false);
                    bp.useFunction = bj.value("use_function", false);
                    bp.expression = bj.value("expression", "");
                    bp.voltage = bj.value("voltage", 0.0);
                    bp.is_corona = bj.value("is_corona", false);
                    m_boundaryConfigs[g] = bp;
                    break;
                }
            }
        }
    }
    // 恢复参数页
    m_config.paramPages.clear();
    if (j.contains("param_pages")) {
        for (auto& jp : j["param_pages"]) {
            ParameterPage pg;
            pg.name = jp.value("name", "参数页");
            if (jp.contains("entries")) {
                for (auto& je : jp["entries"]) {
                    ParamEntry e;
                    e.name = je.value("name", ""); e.expression = je.value("expression", "0");
                    e.value = je.value("value", 0.0); e.unit = je.value("unit", "");
                    e.description = je.value("description", "");
                    pg.entries.push_back(e);
                }
            }
            m_config.paramPages.push_back(pg);
        }
    }

    // 刷新参数页面UI (setupParamUI创建的初始空tab需重建)
    if (m_paramTabs && !m_config.paramPages.empty()) {
        m_paramTabs->blockSignals(true);
        while (m_paramTabs->count() > 0) m_paramTabs->removeTab(0);
        for (size_t pi = 0; pi < m_config.paramPages.size(); pi++) {
            auto* table = new QTableWidget(0, 5);
            table->setHorizontalHeaderLabels({"名称", "表达式", "值", "单位", "描述"});
            table->horizontalHeader()->setStretchLastSection(true);
            table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
            table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
            table->setSelectionBehavior(QAbstractItemView::SelectRows);
            table->setAlternatingRowColors(true);
            table->verticalHeader()->setVisible(false);
            size_t piCap = pi;
            connect(table, &QTableWidget::cellChanged, [this, piCap](int row, int col) {
                if (piCap >= m_config.paramPages.size()) return;
                auto& entries = m_config.paramPages[piCap].entries;
                if (row >= (int)entries.size()) return;
                auto* t = qobject_cast<QTableWidget*>(m_paramTabs->widget((int)piCap));
                if (!t) return;
                m_isUpdatingUI = true;
                if (col == 0) entries[row].name = t->item(row, 0)->text().toStdString();
                else if (col == 1) entries[row].expression = t->item(row, 1)->text().toStdString();
                else if (col == 3) entries[row].unit = t->item(row, 3)->text().toStdString();
                else if (col == 4) entries[row].description = t->item(row, 4)->text().toStdString();
                m_isUpdatingUI = false;
                m_projectModified = true;
                if (col == 1) evalAllParams();
            });
            m_paramTabs->addTab(table, QString::fromStdString(m_config.paramPages[pi].name));
        }
        m_paramTabs->blockSignals(false);
        for (size_t pi = 0; pi < m_config.paramPages.size(); pi++)
            refreshParamTable(pi);
        evalAllParams();
    }

    // 恢复求解器/物理场参数
    if (j.contains("solver_config")) {
        auto& sc = j["solver_config"];
        m_config.physics.E_onset = sc.value("E_onset", 600000.0);
        m_config.physics.rho_surface = sc.value("rho_surface", 10000.0);
        m_config.physics.K_mobility = sc.value("K_mobility", 1.0);
        m_config.physics.wind_x = sc.value("wind_x", 0.0);
        m_config.physics.wind_y = sc.value("wind_y", 0.0);
        m_config.physics.wind_z = sc.value("wind_z", 0.0);
        m_config.solver.w_rho = sc.value("w_rho", 1.0);
        m_config.solver.goal_convergence = sc.value("goal_convergence", 0.95);
        m_config.solver.max_update_times = sc.value("max_update_times", 100);
        m_config.solver.tolerance_E = sc.value("tolerance_E", 0.01);
        m_config.solver.tolerance_rho = sc.value("tolerance_rho", 0.01);
        m_config.solver.order = sc.value("order", 1);
        m_config.solver.num_cores = sc.value("num_cores", 4);
        m_config.solver.solver_type = sc.value("solver_type", "CPU");
        // 恢复UI控件值
        if (m_spinE0) m_spinE0->setValue(m_config.physics.E_onset);
        if (m_spinRho) m_spinRho->setValue(m_config.physics.rho_surface);
        if (m_spinK) m_spinK->setValue(m_config.physics.K_mobility);
        if (m_spinWindX) m_spinWindX->setValue(m_config.physics.wind_x);
        if (m_spinWindY) m_spinWindY->setValue(m_config.physics.wind_y);
        if (m_spinWindZ) m_spinWindZ->setValue(m_config.physics.wind_z);
        if (m_spinWRho) m_spinWRho->setValue(m_config.solver.w_rho);
        if (m_spinGoalConv) m_spinGoalConv->setValue(m_config.solver.goal_convergence);
        if (m_spinTolE) m_spinTolE->setValue(m_config.solver.tolerance_E);
        if (m_spinTolRho) m_spinTolRho->setValue(m_config.solver.tolerance_rho);
        if (m_spinCores) m_spinCores->setValue(m_config.solver.num_cores);
        if (m_comboSolverType) {
            int idx = m_comboSolverType->findData(QString::fromStdString(m_config.solver.solver_type));
            if (idx >= 0) m_comboSolverType->setCurrentIndex(idx);
        }
    }

    // 恢复网格参数
    if (j.contains("mesh_params")) {
        for (auto& mc : j["mesh_params"]) {
            std::string meshName = mc.value("mesh_name", "");
            std::string groupName = mc.value("group_name", "");
            QTreeWidgetItem* mnode = nullptr, *gnode = nullptr;
            for (int i=0; i<m_rootMeshNode->childCount(); i++)
                if (m_rootMeshNode->child(i)->text(0).toStdString() == meshName)
                    { mnode = m_rootMeshNode->child(i); break; }
            for (int i=0; i<m_rootGroupNode->childCount(); i++)
                if (m_rootGroupNode->child(i)->text(0).toStdString() == groupName)
                    { gnode = m_rootGroupNode->child(i); break; }
            if (mnode && gnode) {
                MeshGroupParams p;
                p.enableLocalField = mc.value("enableLocalField", false);
                p.sizeMin = mc.value("sizeMin", 0.05);
                p.sizeMax = mc.value("sizeMax", 8.0);
                p.distMin = mc.value("distMin", 0.0);
                p.distMax = mc.value("distMax", 30.0);
                m_meshConfigs[mnode][gnode] = p;
            }
        }
    }

    // 恢复后处理状态
    if (j.contains("post_process")) {
        auto& pp = j["post_process"];
        if (m_radioSurface) m_radioSurface->setChecked(pp.value("surface", true));
        if (m_chkSliceX) m_chkSliceX->setChecked(pp.value("slice_x", false));
        if (m_spinSliceX) m_spinSliceX->setValue(pp.value("slice_x_pos", 0.0));
        if (m_chkSliceY) m_chkSliceY->setChecked(pp.value("slice_y", false));
        if (m_spinSliceY) m_spinSliceY->setValue(pp.value("slice_y_pos", 0.0));
        if (m_chkSliceZ) m_chkSliceZ->setChecked(pp.value("slice_z", false));
        if (m_spinSliceZ) m_spinSliceZ->setValue(pp.value("slice_z_pos", 0.0));
        if (m_chkAutoRange) m_chkAutoRange->setChecked(pp.value("auto_range", true));
        if (m_spinRangeMin) m_spinRangeMin->setValue(pp.value("range_min", 0.0));
        if (m_spinRangeMax) m_spinRangeMax->setValue(pp.value("range_max", 0.0));
        if (m_listPostGroups) {
            auto sg = pp.value("selected_groups", json::array());
            for (int i=0; i<m_listPostGroups->count(); i++)
                m_listPostGroups->item(i)->setSelected(false);
            for (auto& s : sg) {
                for (int i=0; i<m_listPostGroups->count(); i++)
                    if (m_listPostGroups->item(i)->text().toStdString() == s.get<std::string>())
                        m_listPostGroups->item(i)->setSelected(true);
            }
        }
    }

    if (j.contains("last_field")) m_lastFieldName = QString::fromStdString(j["last_field"].get<std::string>());
    if (j.contains("last_vtk")) m_lastVtkFilePath = QString::fromStdString(j["last_vtk"].get<std::string>());

    // 重新加载网格数据(从缓存)
    for (int i = 0; i < m_rootMeshNode->childCount(); i++) {
        auto* mn = m_rootMeshNode->child(i);
        if (m_meshFiles.count(mn)) {
            QString cachePath = m_meshFiles[mn] + ".meshcache";
            if (QFileInfo::exists(cachePath))
                m_meshDataMap[mn] = MeshEngine::loadMeshCache(cachePath.toStdString());
        }
    }

    // 恢复参数绑定
    if (j.contains("param_bindings")) {
        std::map<std::string, QDoubleSpinBox*> nameToSpin;
        for (auto& [spin, sname] : m_spinboxNames) nameToSpin[sname] = spin;
        for (auto& bj : j["param_bindings"]) {
            std::string sname = bj.value("spinbox", "");
            std::string pname = bj.value("param", "");
            if (nameToSpin.count(sname)) {
                auto* spin = nameToSpin[sname];
                m_paramBindings[spin] = pname;
                spin->setStyleSheet("QDoubleSpinBox { background-color: #E3F2FD; }");
                spin->setReadOnly(true);
            }
        }
        pushParamsToBindings();
    }

    // 清除导入组残留+恢复内部组(仅当存在网格方案时)
    if (m_rootMeshNode->childCount() > 0) {
    QTreeWidgetItem* firstInternal = nullptr;
    for (int i = 0; i < m_rootMeshNode->childCount(); i++) {
        auto* mn = m_rootMeshNode->child(i);
        if (!m_isImportedMesh.count(mn) || !m_isImportedMesh[mn]) { firstInternal = mn; break; }
    }
    // 暂存边界配置(按组名, 避免清理后指针失效)
    std::map<std::string, BoundaryParams> savedBdrs;
    for (auto& [node, bp] : m_boundaryConfigs)
        savedBdrs[node->text(0).toStdString()] = bp;
    for (int i = m_rootGroupNode->childCount()-1; i >= 0; i--) delete m_rootGroupNode->child(i);
    m_groupData.clear(); m_groupTags.clear(); m_boundaryConfigs.clear();
    if (firstInternal) {
        m_activeMeshNode = firstInternal;
        if (m_importedGroupMeta.count(firstInternal)) {
            for (auto& gm : m_importedGroupMeta[firstInternal]) {
                auto* item = new QTreeWidgetItem(m_rootGroupNode);
                item->setText(0, gm.name);
                if (gm.tag > 0) { item->setData(0, Qt::UserRole, gm.tag); m_groupTags[item] = gm.tag; }
                m_groupData[item] = gm.faceIds;
            }
        }
    }
    // 按名恢复边界配置到新指针
    for (auto& [name, bp] : savedBdrs)
        for (int i = 0; i < m_rootGroupNode->childCount(); i++)
            if (m_rootGroupNode->child(i)->text(0).toStdString() == name)
                { m_boundaryConfigs[m_rootGroupNode->child(i)] = bp; break; }
    m_rootGroupNode->setExpanded(true);
    } // end if mesh node count > 0

    refreshBoundaryUI();
    refreshSolverMeshCombo();
    logMessage("项目已加载: " + path, "success");
}

bool MainWindow::maybeSave() {
    if (!m_projectModified) return true;
    auto ret = QMessageBox::warning(this, "未保存的更改",
        "当前项目有未保存的更改。是否保存？",
        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    if (ret == QMessageBox::Save) {
        if (m_projectFilePath.isEmpty()) {
            QString path = QFileDialog::getSaveFileName(this, "保存项目", "", "IonFlow Project (*.ion)");
            if (path.isEmpty()) return false;
            m_projectFilePath = path;
        }
        saveProject(m_projectFilePath);
        return true;
    }
    return ret == QMessageBox::Discard;
}

void MainWindow::showCaseLibrary() {
    QString casesDir = QCoreApplication::applicationDirPath() + "/../cases/";
    QDir dir(casesDir);
    if (!dir.exists()) dir.mkpath(".");

    // 读取案例列表
    std::vector<CaseInfo> cases;
    QStringList categories = {"全部"};
    QStringList catSet;
    for (QString sub : dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot)) {
        QString infoPath = casesDir + sub + "/info.json";
        QFileInfo fi(infoPath);
        if (!fi.exists()) continue;
        QFile f(infoPath);
        if (f.open(QIODevice::ReadOnly)) {
            auto j = json::parse(f.readAll().toStdString(), nullptr, false);
            f.close();
            if (!j.is_discarded()) {
                CaseInfo ci;
                ci.name = j.value("name", sub.toStdString());
                ci.category = j.value("category", "");
                ci.description = j.value("description", "");
                ci.folderName = sub.toStdString();
                cases.push_back(ci);
                if (!ci.category.empty() && !catSet.contains(QString::fromStdString(ci.category))) {
                    catSet << QString::fromStdString(ci.category);
                    categories << QString::fromStdString(ci.category);
                }
            }
        }
    }

    // 对话框
    QDialog dlg(this);
    dlg.setWindowTitle("案例库");
    dlg.resize(700, 500);
    auto* mainLay = new QHBoxLayout(&dlg);

    // 左侧: 分类筛选
    auto* leftPanel = new QWidget();
    auto* leftLay = new QVBoxLayout(leftPanel);
    leftLay->setContentsMargins(0,0,0,0);
    auto* catLabel = new QLabel("分类筛选");
    catLabel->setStyleSheet("font-weight:bold;");
    leftLay->addWidget(catLabel);
    auto* catList = new QListWidget();
    catList->addItems(categories);
    catList->setCurrentRow(0);
    catList->setMaximumWidth(120);
    leftLay->addWidget(catList);
    mainLay->addWidget(leftPanel);

    // 右侧: 案例列表 + 详情
    auto* rightPanel = new QWidget();
    auto* rightLay = new QVBoxLayout(rightPanel);
    rightLay->setContentsMargins(0,0,0,0);
    auto* searchEdit = new QLineEdit();
    searchEdit->setPlaceholderText("搜索案例...");
    rightLay->addWidget(searchEdit);

    auto* caseList = new QListWidget();
    rightLay->addWidget(caseList);

    auto* detailLabel = new QLabel();
    detailLabel->setWordWrap(true);
    detailLabel->setStyleSheet("color:#555; padding:4px;");
    rightLay->addWidget(detailLabel);
    mainLay->addWidget(rightPanel);

    // 底部按钮
    auto* btnLay = new QHBoxLayout();
    auto* btnImport = new QPushButton("导入案例");
    auto* btnOpen = new QPushButton("打开案例");
    auto* btnDelete = new QPushButton("删除");
    auto* btnClose = new QPushButton("关闭");
    btnLay->addWidget(btnImport);
    btnLay->addWidget(btnOpen);
    btnLay->addWidget(btnDelete);
    btnLay->addStretch();
    btnLay->addWidget(btnClose);

    // 底部按钮放入右侧面板底部
    rightLay->addLayout(btnLay);

    // 刷新列表
    auto refreshList = [&](const QString& filterCat, const QString& search) {
        caseList->clear();
        for (auto& ci : cases) {
            QString cat = QString::fromStdString(ci.category);
            if (filterCat != "全部" && cat != filterCat) continue;
            QString name = QString::fromStdString(ci.name);
            if (!search.isEmpty() && !name.contains(search, Qt::CaseInsensitive)) continue;
            auto* item = new QListWidgetItem(name);
            item->setData(Qt::UserRole, (int)(&ci - &cases[0])); // store index
            caseList->addItem(item);
        }
    };
    refreshList("全部", "");

    connect(catList, &QListWidget::currentTextChanged, [&](const QString& t){ refreshList(t, searchEdit->text()); });
    connect(searchEdit, &QLineEdit::textChanged, [&](const QString& t){ refreshList(catList->currentItem() ? catList->currentItem()->text() : "全部", t); });
    connect(caseList, &QListWidget::currentTextChanged, [&](){
        if (caseList->currentItem()) {
            int idx = caseList->currentItem()->data(Qt::UserRole).toInt();
            detailLabel->setText(QString("名称: %1\n分类: %2\n描述: %3")
                .arg(QString::fromStdString(cases[idx].name))
                .arg(QString::fromStdString(cases[idx].category))
                .arg(QString::fromStdString(cases[idx].description)));
        }
    });

    // 导入案例
    connect(btnImport, &QPushButton::clicked, [&](){
        QString ionFile = QFileDialog::getOpenFileName(&dlg, "导入案例", "", "IonFlow Project (*.ion)");
        if (ionFile.isEmpty()) return;
        QDialog infoDlg(&dlg);
        infoDlg.setWindowTitle("案例信息");
        auto* fl = new QFormLayout(&infoDlg);
        auto* nameEdit = new QLineEdit(QFileInfo(ionFile).baseName());
        auto* catEdit = new QComboBox(); catEdit->setEditable(true);
        for (int i = 1; i < categories.size(); i++) catEdit->addItem(categories[i]);
        catEdit->setCurrentText("基础几何");
        auto* descEdit = new QLineEdit();
        fl->addRow("名称:", nameEdit);
        fl->addRow("分类:", catEdit);
        fl->addRow("描述:", descEdit);
        auto* bb = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        fl->addRow(bb);
        connect(bb, &QDialogButtonBox::accepted, &infoDlg, &QDialog::accept);
        connect(bb, &QDialogButtonBox::rejected, &infoDlg, &QDialog::reject);
        if (infoDlg.exec() == QDialog::Accepted) {
            QString caseDir = casesDir + nameEdit->text();
            QDir().mkpath(caseDir);
            QFile::copy(ionFile, caseDir + "/project.ion");
            json info;
            info["name"] = nameEdit->text().toStdString();
            info["category"] = catEdit->currentText().toStdString();
            info["description"] = descEdit->text().toStdString();
            std::ofstream o((caseDir + "/info.json").toStdString());
            o << std::setw(2) << info << std::endl;
            QMessageBox::information(&dlg, "成功", "案例已导入案例库");
            dlg.accept(); // close and reopen
        }
    });

    // 打开案例
    connect(btnOpen, &QPushButton::clicked, [&](){
        if (!caseList->currentItem()) return;
        if (!maybeSave()) return;
        int idx = caseList->currentItem()->data(Qt::UserRole).toInt();
        QString srcPath = casesDir + QString::fromStdString(cases[idx].folderName) + "/project.ion";
        QString tmpPath = QDir::tempPath() + "/IonFlow_case_open.ion";
        QFile::remove(tmpPath);
        if (!QFile::copy(srcPath, tmpPath)) {
            QMessageBox::warning(&dlg, "错误", "无法打开案例文件"); return;
        }
        // 清理当前项目再加载
        m_geomNodes.clear();
        for (int i = m_rootGeomNode->childCount()-1; i >= 0; i--) delete m_rootGeomNode->child(i);
        for (int i = m_rootMeshNode->childCount()-1; i >= 0; i--) delete m_rootMeshNode->child(i);
        for (int i = m_rootGroupNode->childCount()-1; i >= 0; i--) delete m_rootGroupNode->child(i);
        m_groupData.clear(); m_groupTags.clear(); m_boundaryConfigs.clear();
        m_meshConfigs.clear(); m_meshDataMap.clear(); m_meshFiles.clear();
        m_isImportedMesh.clear(); m_importedGroupMeta.clear(); m_importedBoundaryConfigs.clear();
        m_activeMeshNode = nullptr; m_currentStepFile.clear();
        ui->cadViewerWidget->clearScene();
        loadProject(tmpPath);
        m_projectFilePath.clear();
        m_projectModified = true;
        setWindowTitle(QString("离子流场模拟器 - 案例: %1 (未保存)").arg(QString::fromStdString(cases[idx].name)));
        dlg.accept();
    });

    // 删除案例
    connect(btnDelete, &QPushButton::clicked, [&](){
        if (!caseList->currentItem()) return;
        int idx = caseList->currentItem()->data(Qt::UserRole).toInt();
        auto ret = QMessageBox::question(&dlg, "确认删除", QString("确定要删除案例 '%1' 吗?").arg(QString::fromStdString(cases[idx].name)));
        if (ret == QMessageBox::Yes) {
            QString caseDir = casesDir + QString::fromStdString(cases[idx].folderName);
            QDir(caseDir).removeRecursively();
            QMessageBox::information(&dlg, "已删除", "案例已从案例库中移除");
            dlg.accept();
        }
    });

    connect(btnClose, &QPushButton::clicked, &dlg, &QDialog::reject);
    dlg.exec();
}

void MainWindow::closeEvent(QCloseEvent* event) {
    if (maybeSave()) event->accept();
    else event->ignore();
}