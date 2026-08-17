#include "geom_icons.h"
#include <QPainter>
#include <QPixmap>
#include <QPainterPath>
#include <QMap>
#include <QFileInfo>
#include <QDateTime>
#include <QCoreApplication>
#include <QtMath>

// ====================================================================
// 几何体素/操作线框图标生成器
// 20x20 透明画布, 深灰色线条, CAD 风格
// ====================================================================

static void drawArrowHead(QPainter& p, const QPointF& tip, const QPointF& dir) {
    // 箭头头部: 两条斜线, dir 为箭头指向的单位向量
    double len = 3.2;
    QPointF side(-dir.y(), dir.x());
    p.drawLine(tip, tip - dir * len + side * (len * 0.55));
    p.drawLine(tip, tip - dir * len - side * (len * 0.55));
}

static void drawBoxWire(QPainter& p, QPointF top, QPointF right, QPointF bottom, QPointF left,
                        double depthX, double depthY) {
    // 轴测立方体线框: 顶面 (top-right-bottom-left 菱形) + 拉伸到深度 (depthX,depthY)
    QPointF t2 = top + QPointF(depthX, depthY);
    QPointF r2 = right + QPointF(depthX, depthY);
    QPointF b2 = bottom + QPointF(depthX, depthY);
    QPointF l2 = left + QPointF(depthX, depthY);
    p.drawPolygon(QPolygonF({top, right, bottom, left}));          // 顶面
    p.drawLine(right, r2); p.drawLine(bottom, b2); p.drawLine(left, l2);  // 侧棱
    p.drawPolygon(QPolygonF({right, r2, b2, bottom}));             // 右侧面
    p.drawPolygon(QPolygonF({left, l2, b2, bottom}));              // 左前面
}

static QIcon makeIcon(const QString& type) {
    QPixmap pm(20, 20);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing);
    QPen solid(QColor(70, 70, 70), 1.5);
    p.setPen(solid);

    if (type == "立方体") {
        drawBoxWire(p, {10, 3}, {16.5, 6.5}, {10, 10}, {3.5, 6.5}, 3.2, 6.2);
    } else if (type == "球体") {
        p.drawEllipse(QPointF(10, 10), 6.8, 6.8);
        p.drawEllipse(QPointF(10, 10), 6.8, 2.6);
        p.drawLine(QPointF(10, 3.2), QPointF(10, 16.8));
    } else if (type == "圆柱") {
        p.drawEllipse(QRectF(3.2, 3.0, 13.6, 4.6));   // 上椭圆
        p.drawEllipse(QRectF(3.2, 12.4, 13.6, 4.6));  // 下椭圆
        p.drawLine(QPointF(4.0, 5.3), QPointF(4.0, 14.7));
        p.drawLine(QPointF(16.0, 5.3), QPointF(16.0, 14.7));
    } else if (type == "圆锥") {
        p.drawLine(QPointF(10, 2.6), QPointF(3.2, 14.5));
        p.drawLine(QPointF(10, 2.6), QPointF(16.8, 14.5));
        p.drawEllipse(QRectF(3.2, 12.6, 13.6, 4.4));  // 底椭圆
        p.drawLine(QPointF(10, 12.6), QPointF(10, 17.0));
    } else if (type == "圆环") {
        p.drawEllipse(QPointF(10, 10), 7.2, 7.2);     // 外圆
        p.drawEllipse(QPointF(10, 10), 3.0, 3.0);     // 内圆
        p.drawLine(QPointF(2.8, 10), QPointF(6.0, 10));
        p.drawLine(QPointF(14.0, 10), QPointF(17.2, 10));
    } else if (type == "楔形") {
        // 楔形 = 三角形棱柱
        QPolygonF tri({QPointF(4, 4.5), QPointF(16, 4.5), QPointF(16, 12)});
        p.drawPolygon(tri);
        p.drawLine(QPointF(6, 13.5), QPointF(13, 13.5));
        p.drawLine(QPointF(4, 4.5), QPointF(6, 13.5));
        p.drawLine(QPointF(16, 12), QPointF(13, 13.5));
        p.drawLine(QPointF(16, 4.5), QPointF(13, 13.5));  // 后棱 (可见简化)
    } else if (type == "椭球") {
        p.drawEllipse(QPointF(10, 10), 7.0, 5.0);
        p.drawEllipse(QPointF(10, 10), 7.0, 2.2);
        p.drawLine(QPointF(10, 5), QPointF(10, 15));
    } else if (type == "棱柱") {
        // 三棱柱: 前三角 + 顶/底棱线 + 后三角
        p.drawPolygon(QPolygonF({QPointF(10, 3.2), QPointF(3.5, 11.5), QPointF(16.5, 11.5)}));
        p.drawLine(QPointF(12, 8.5), QPointF(5.5, 16));
        p.drawLine(QPointF(10, 3.2), QPointF(12, 8.5));
        p.drawLine(QPointF(16.5, 11.5), QPointF(12, 16));
        p.drawLine(QPointF(5.5, 16), QPointF(12, 16));
    } else if (type == "棱锥") {
        QPolygonF base({QPointF(4, 13), QPointF(16, 13), QPointF(16, 17), QPointF(4, 17)});
        p.drawPolygon(base);
        p.drawLine(QPointF(10, 3), QPointF(4, 13));
        p.drawLine(QPointF(10, 3), QPointF(16, 13));
        p.drawLine(QPointF(10, 3), QPointF(16, 17));
        p.drawLine(QPointF(10, 3), QPointF(4, 17));
    } else if (type == "圆台") {
        p.drawEllipse(QRectF(4.5, 3.0, 11.0, 4.0));   // 上椭圆
        p.drawEllipse(QRectF(3.0, 13.0, 14.0, 4.0));  // 下椭圆
        p.drawLine(QPointF(5.5, 5.0), QPointF(4.0, 15.0));
        p.drawLine(QPointF(14.5, 5.0), QPointF(16.0, 15.0));
    } else if (type == "棱台") {
        QPolygonF top({QPointF(7, 4.5), QPointF(13.5, 4.5), QPointF(15.5, 8), QPointF(9, 8)});
        p.drawPolygon(top);
        p.drawLine(QPointF(7, 4.5), QPointF(4, 13));
        p.drawLine(QPointF(13.5, 4.5), QPointF(16.5, 13));
        p.drawLine(QPointF(15.5, 8), QPointF(16.5, 13));
        p.drawLine(QPointF(9, 8), QPointF(4, 13));
        p.drawPolygon(QPolygonF({QPointF(4, 13), QPointF(16.5, 13), QPointF(16.5, 16.5), QPointF(4, 16.5)}));
    } else if (type == "并集") {
        p.drawEllipse(QPointF(8.2, 10.5), 5.4, 5.4);
        p.drawEllipse(QPointF(11.8, 10.5), 5.4, 5.4);
    } else if (type == "差集") {
        p.drawEllipse(QPointF(8.2, 10.5), 5.4, 5.4);
        QPen dashed(solid); dashed.setStyle(Qt::DashLine);
        p.setPen(dashed);
        p.drawEllipse(QPointF(11.8, 10.5), 5.4, 5.4);
    } else if (type == "交集") {
        // 透镜形: 两圆公共区域轮廓 (用 QPainterPath 求交太复杂, 手绘透镜)
        QPainterPath lens;
        lens.moveTo(10, 5.2);
        lens.quadTo(12.6, 10.5, 10, 15.8);   // 右弧
        lens.quadTo(7.4, 10.5, 10, 5.2);     // 左弧
        p.drawPath(lens);
        p.drawLine(QPointF(10, 5.2), QPointF(10, 15.8));  // 中心线
    } else if (type == "倒圆角") {
        QPolygonF shape;
        shape << QPointF(4, 13) << QPointF(4, 7) << QPointF(13, 7)
              << QPointF(13, 4) << QPointF(16, 4) << QPointF(16, 16) << QPointF(4, 16);
        p.drawPolyline(shape);
        p.drawArc(QRectF(11.5, 5.5, 4.5, 4.5), 0, 90 * 16);  // 圆角
    } else if (type == "倒斜角") {
        p.drawPolyline(QPolygonF({QPointF(4, 13), QPointF(4, 7), QPointF(11, 7),
                                  QPointF(13.5, 4), QPointF(16, 4), QPointF(16, 16), QPointF(4, 16)}));
    } else if (type == "扫掠") {
        // 曲线路径 + 箭头 + 截面方块
        QPainterPath path;
        path.moveTo(3, 15.5);
        path.cubicTo(7, 12, 13, 12, 17, 5.5);
        p.drawPath(path);
        drawArrowHead(p, {17, 5.5}, {-0.55, -0.84});
        p.drawRect(QRectF(4.5, 12.5, 3.5, 3.5));
    } else if (type == "平移") {
        p.drawLine(QPointF(2.5, 10), QPointF(14, 10));
        drawArrowHead(p, {15.5, 10}, {1, 0});
        p.drawRect(QRectF(4.5, 6, 3, 3));   // 被平移对象
    } else if (type == "旋转") {
        p.drawArc(QRectF(5, 4, 11, 11), 45 * 16, 250 * 16);
        drawArrowHead(p, {10, 15}, {-0.9, 0.45});
        p.drawLine(QPointF(10, 9.5), QPointF(10, 14.5));  // 旋转轴提示
    } else if (type == "镜像") {
        QPen dashed(solid); dashed.setStyle(Qt::DashLine);
        p.setPen(dashed);
        p.drawLine(QPointF(10, 3), QPointF(10, 17));
        p.setPen(solid);
        p.drawPolygon(QPolygonF({QPointF(6.5, 7), QPointF(6.5, 13), QPointF(3, 10)}));
        p.drawPolygon(QPolygonF({QPointF(13.5, 7), QPointF(13.5, 13), QPointF(17, 10)}));
    } else if (type == "缩放") {
        p.drawEllipse(QPointF(8.5, 8.5), 4.8, 4.8);       // 镜片
        p.drawLine(QPointF(12, 12), QPointF(17, 17));     // 柄
        drawArrowHead(p, {17.5, 17.5}, {0.7, 0.7});
    } else if (type == "偏移") {
        p.drawRect(QRectF(4, 4, 12, 12));
        p.drawRect(QRectF(7.5, 7.5, 5, 5));
    } else if (type == "线性阵列") {
        p.drawRect(QRectF(2.5, 7, 4, 6));
        p.drawRect(QRectF(8, 7, 4, 6));
        p.drawRect(QRectF(13.5, 7, 4, 6));
    } else if (type == "圆形阵列") {
        // 圆周上 6 个方块 + 中心点
        p.drawEllipse(QPointF(10, 10), 1.0, 1.0);
        for (int i = 0; i < 6; i++) {
            double ang = M_PI / 3.0 * i;
            QPointF c(10 + 6.2 * qCos(ang), 10 + 6.2 * qSin(ang));
            p.drawRect(QRectF(c.x() - 1.6, c.y() - 1.6, 3.2, 3.2));
        }
    } else {
        // 未知类型: 问号占位
        p.setFont(QFont("Sans", 9, QFont::Bold));
        p.drawText(pm.rect(), Qt::AlignCenter, "?");
    }

    p.end();
    return QIcon(pm);
}

QIcon geomIcon(const QString& type) {
    // 缓存 key 带文件时间戳: 用户替换图标文件后重启即可生效
    struct Entry { QIcon icon; qint64 mtime = 0; };
    static QMap<QString, Entry> cache;

    // 依次尝试用户自定义图标文件: <appDir>/icons/<type>.png / .svg / .jpg
    static const QString iconDir = QCoreApplication::applicationDirPath() + "/icons";
    QString customFile;
    qint64 mtime = 0;
    for (const QString& ext : {QString(".png"), QString(".svg"), QString(".jpg")}) {
        QString cand = iconDir + "/" + type + ext;
        QFileInfo fi(cand);
        if (fi.exists()) { customFile = cand; mtime = fi.lastModified().toMSecsSinceEpoch(); break; }
    }

    auto it = cache.find(type);
    if (it != cache.end() && it->mtime == mtime) return it->icon;  // 缓存命中

    Entry e;
    e.mtime = mtime;
    e.icon = customFile.isEmpty() ? makeIcon(type)            // 无自定义文件 → 自绘回退
                                  : QIcon(customFile);        // 有自定义文件 → 用户图标
    cache[type] = e;
    return e.icon;
}
