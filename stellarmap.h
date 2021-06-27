#ifndef STELLARMAP_H
#define STELLARMAP_H
#include <QWidget>
#include <QObject>
#include <QGridLayout>
#include <QDebug>
#include <QPen>
#include <QStringList>
#include "qcustomplot.h"
#include <synchapi.h>
#include <QDateTime>
#include <QMap>
#include <rinexdata.h>
#include <QThread>
#include <QTimer>
#include <QColor>
#include <QBrush>
#include <QFontDialog>
#include <QFont>
#define PI 3.14159265358979323846  //圆周率

/*******************************
 ***************卫星形状类********
 ******************************/
class SatellitePicture:public QCPItemEllipse
{
     Q_OBJECT
public:
    explicit SatellitePicture(QCustomPlot *parentPlot,const QString &sateName);
    void setPosition(const double &L,const double &B);  //设置卫星在坐标系中的位置，入参纬经纬度(单位为十进制度)
    ~SatellitePicture();
    QString PRN;  //卫星PRN号
    QCPItemText *text;//显示卫星PRN号的文本

};
class StellarMap :public QCustomPlot   //星空图类，用于显示卫星的位置
{
      Q_OBJECT
public:
    explicit StellarMap(QWidget *parent = nullptr);
     ~StellarMap();
    //显示卫星,SatelliteNames是卫星名,L是经度,B是纬度
    void showSatellites();
    double DateToSeconds(const QString &date);//星期几转换成该周的第几秒
    void addSatellite(const GPS &sate); //增加GPS卫星
    void addSatellite(const BDS &sate); //增加BDS卫星
    void addSatellite(const GLONASS &sate); //增加GLONASS卫星
    void addSatellite(const GALILEO &sate); //增加GALILEO卫星
    QMap<QString,QList<double>>  satelliteTable;   //用来存储当前时刻所有的卫星坐标
protected:
    void initStellarMap();
    QCPItemEllipse **latitudeLines;  //纬线
    QCPItemLine **longtitudeLines;//经线
    QThread *calThread;  //计算线程
    QTimer *timeTrigger; //时间触发器;
    //卫星
    QMap<QString,GPS> gps;
    QMap<QString,BDS> bds;
    QMap<QString,SatellitePicture*> SateDraws;  //存储所有卫星的形状
    int interval;  //每次计算卫星位置的时间间隔
signals:
   void tableSignal(QMap<QString,QList<double>>  data);
};
#endif // STELLARMAP_H
