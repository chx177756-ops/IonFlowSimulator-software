#include "mainwindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QDir>
#include <QDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    qputenv("LIBGL_ALWAYS_SOFTWARE", "1");

    QSurfaceFormat format;
    format.setVersion(3, 3);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    QSurfaceFormat::setDefaultFormat(format);

    QApplication a(argc, argv);
    QDir::setCurrent(QCoreApplication::applicationDirPath());

    // ---- 启动前项目选择对话框 (COMSOL 风格) ----
    QDialog dlg;
    dlg.setWindowTitle("离子流场模拟器 - 新建/打开项目");
    dlg.setFixedSize(420, 180);
    auto* layout = new QVBoxLayout(&dlg);
    layout->addWidget(new QLabel("<h3>欢迎使用离子流场模拟器</h3>"));
    layout->addWidget(new QLabel("请新建一个项目或打开已有项目以开始工作："));
    auto* btnLayout = new QHBoxLayout();
    QPushButton* btnNew = new QPushButton("新建项目");
    QPushButton* btnOpen = new QPushButton("打开项目");
    btnNew->setMinimumHeight(36); btnOpen->setMinimumHeight(36);
    btnLayout->addWidget(btnNew); btnLayout->addWidget(btnOpen);
    layout->addLayout(btnLayout);

    QString projectPath;

    QObject::connect(btnNew, &QPushButton::clicked, [&](){
        projectPath = "__untitled__";  // 标记为未保存项目, 直接进入
        dlg.accept();
    });

    QObject::connect(btnOpen, &QPushButton::clicked, [&](){
        QString path = QFileDialog::getOpenFileName(&dlg, "打开项目", QDir::homePath(), "IonFlow Project (*.ion)");
        if (!path.isEmpty()) { projectPath = path; dlg.accept(); }
    });

    dlg.exec();

    if (projectPath.isEmpty()) return 0;  // 用户关闭对话框 → 退出

    MainWindow w;
    w.setProjectPath(projectPath);
    w.show();
    return a.exec();
}