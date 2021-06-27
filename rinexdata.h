#ifndef RINEXDATA_H
#define RINEXDATA_H
#include<QObject>
#include<QStringList>
#include <QVector>
#include <QString>
#include <QFile>
#include <QDebug>
#include <QRegExp>
#include <QList>
#include <qmath.h>
#include <eigen/Eigen/Dense>
using namespace Eigen;
#define PI 3.14159265358979323846  //圆周率

/* class Satellite is a abstract class,it offers the following methods:
 * import ephemeris,calculate the position of satellite and return the coordinate ...
 */
class Satellite
{
public:
    Satellite();
    virtual ~Satellite()=0;
    //return the space rectangular coordinate of a satellite in the Earth-Fixed Coordinate System
    virtual QVector<double> getSatellitePosition(const double &t)=0;
    virtual void printEphemeris()=0;
    virtual QVector<double> getSatelliteLBH(const double &t); //获取卫星的大地坐标(经纬度及高程)
    QVector<double> LBHto3DCoordinate(const QVector<double> &LBH);//大地坐标转空间直角坐标
    struct Ephemeris;
protected:
    virtual void updateEphemeris(const QStringList &ephemeris)=0;

    double GM; //this ia a geocentric gravitational constant,its unit is m^3/s^2
    short int leapSeconds=18;
    double w_e;  //earth's rotation angular velocity, its unit is rad
    double a;//椭球长半轴，单位是m
    double e1s;//第一偏心率的平方
    double e2s;//第二偏心率的平方
};
class GLONASS:public Satellite
{
public:
    GLONASS(const QStringList &ephemerisData);
    virtual ~GLONASS() override;
    virtual QVector<double> getSatellitePosition(const double &t) override;
    virtual void printEphemeris() override;
    struct Ephemeris
    {
        QString PRN; QString epochTime; double r_n;    double tau_n;  double fTime;
                     double x;          double x_v;    double x_a;    double svhlth;
                     double y;          double y_v;    double y_a;    double freNum;
                     double z;          double z_v;    double z_a;    double ageLimit;
    };
     Ephemeris sateEphemeris;
protected:
    virtual void updateEphemeris(const QStringList &List) override; //add or update ephemeris data
};
class GPS:public Satellite
{
public:
    GPS(const QStringList &ephemerisData);//ephemerisData为8行的星历数据
    virtual ~GPS() override;
    virtual QVector<double> getSatellitePosition(const double &t) override;
    virtual void printEphemeris() override;
    struct Ephemeris
    {
        QString PRN; QString epochTime; double a_0;    double a_1;     double a_2;
                     double IODE;       double C_rs;   double delta_n; double M_0;
                     double C_uc;       double e;      double C_us;    double sqrt_a;
                     double t_oe;       double C_ic;   double Omega_0; double C_is;
                     double i_0;        double C_rc;   double omega;   double Omega_p;
                     double i_p;        double cflgl2; double weekno;  double pflgl2;
                     double svacc;      double svhlth; double tgd;     double IODC;
                     double ttm;        double fi;
    };
    Ephemeris sateEphemeris;
protected:
    virtual void updateEphemeris(const QStringList &List) override; //add or update ephemeris data   
};
class BDS:public GPS
{
public:
    BDS();
    BDS(const QStringList &ephemerisData);
    virtual ~BDS() override;
/******************************************
* brief:BDS卫星的MEO和IGSO卫星的坐标计算方式与GPS几乎一样，
*       只有GEO卫星的计算方式稍有不同，所以需定义一个新虚
*       函数来解算GEO卫星的坐标。IGSO与MEO卫星使用继承自
*       GPS卫星计算方式即可。
******************************************/
    virtual QVector<double> getSatellitePosition(const double &t) override;
    bool isValid(); //判断该卫星是否有效
     //其中北斗有一些在轨试验卫星和故障卫星没有加进去
    QList<QString> MEOPRN={"C11","C12","C14","C19","C20","C21","C22","C23","C24",
                           "C25","C26","C27","C28","C29","C30","C32","C33","C34",
                           "C35","C36","C37","C41","C42","C43","C44","C45","C46"};
    QList<QString> IGSOPRN={"C06","C07","C08","C09","C10","C13","C16","C38","C39","C40"};
    QList<QString> GEOPRN={"C01","C02","C03","C04","C05","C59","C60","C61"};
protected:

};
class GALILEO:public GPS
{
public:
    GALILEO(const QStringList &ephemerisData);
    virtual ~GALILEO() override;
};
/*==========================================================
 *======================观测文件数据类========================
 *==========================================================*/
class Observation
{
public:
    Observation(const QStringList &List);
    ~Observation();
    QMap<QString,double> waveFrequency;
private:
    void importFile(const QStringList &List); //输入文件
    struct obsTypes  //该结构体用于存储观测类型
    {
        QString sateName;     //卫星系统名称(一位字母)：G表示GPS，R表示GLONASS，C表示BDS，E表示Galileo
        short typeNum;        //观测类型数量
        QList<QString> types; //观测类型
        QMap<QString,double> phaseShift;//相位改正:key为系统名(1个字母)，值为改正值
    };
    struct epoch   //该结构体用于存储每个历元的数据
    {   /*****************************************************
        * 历元标志：
        *               0表示正常;
        *               1表示在前一历元与当前历元之间发生了电源故障；
        * 大于1为事件标志:2表示天线开始移动;
        *               3表示新设站(动态数据结束)(后面至少需要跟上MARKERNAME记录);
        *               4表示后面紧跟着的是类似于文件头的信息,用于说明观测过程中所发生的一些特殊情况;
        *               5表示外部事件(历元时刻与观测值时标属于相同的时间框架)，
        *               6表示后面为描述所探测出并已被修复周跳的记录(格式与OBSERVATIONS记录相同,不过，
        *                ,用周跳替代了观测值,LLI和信号强度为空格或0)。此项为可选项。
        *******************************************************/
       QString time[6];      //—观测历元时刻:年、月、日、时、分、秒;
       short Marker;         //历元标志
       short sateNum;        //当前历元所观测到的卫星数，被用来说明紧跟在后面的记录数，即后面共有几行用于事件的描述;
       double  clockOffset;   //接收机钟差，单位为s，为可选项
       /*********************************************************
        * 下面的容器用于存储观测值：
        * key(关键字):表示卫星PRN号
        * value(值):为二维数组，每行包括三个值，分别为观测值、LLI、信号强度
        * *******************************************************/
       QMap<QString,QList<QList<double>>> C;  //BDS
       QMap<QString,QList<QList<double>>> G;  //GPS
       QMap<QString,QList<QList<double>>> R;  //GLONASS
       QMap<QString,QList<QList<double>>> E;  //GALILEO
    };
    QString MARKERNAME;              //测站点名
    double APPROXPOSITIONXYZ[3]={0,0,0}; //测站的近似坐标，依次为X,Y,Z,单位为m
    double ANTENNADELTAHEN[3]={0,0,0};   //天线偏差,依次为H、E、N,单位为m
    double INTERVAL;                 //历元间隔，单位为s
    float LEAPSECONDS;               //跳秒数
    bool RCVCLOCKOFFSAPPL;           //是否进行了实时的接收机钟改正：0或空格表示没有；1表示进行了接收机钟改正
    QString TIMEOFFIRSTOBS[7];       //开始观测的时间，依次是：卫星系统、年、月、日、时、分、秒
    QString TIMEOFLASTOBS[7];        //结束观测的时间，依次是：卫星系统、年、月、日、时、分、秒
    QList<epoch> EPOCH;              //存储历元
    QList<obsTypes> OBSTYPES;        //存储观测类型
};
/*==========================================================
 *======================导航电文与观测文件管理==================
 *==========================================================*/
class NavigationObservationData
{
   // constructor is set to private to prevent external instantiation(this class is a singleton pattern)
private:
    NavigationObservationData();
    static NavigationObservationData* m_pInstance;  //class pointer is set to private to prevent to external destruction
public:
   static NavigationObservationData* GetInstance()
    {
        {
            if ( m_pInstance == nullptr )  //creat a object if the class haven't been created,
                m_pInstance = new NavigationObservationData();
            return m_pInstance;         //return the address of this object
        }
    }
    enum satelliteSystem{Gps,Bds, Glonass, Galileo};
    ~NavigationObservationData();
    /*void importFiles(const QStringList &paths)：
     * This function is used to import Rinex data Such as:
     * observation file,
     * navigation message,
     * meteorology file,
     * ...
     */
    void importFiles(const QStringList &paths);
    QVector<double> getSatellitePsition(const QString &satelliteName,const QString &epoch);
    void printEphemeris(satelliteSystem systemName);
    static double radianTodecimaldegree(const double &degree); //弧度转10进制度
    void printObservation();          //打印O文件内容
    //输出卫星
    QList<GPS> getGPSSatellites();
    QList<BDS> getBDSSatellites();
    QList<GALILEO> getGALILEOSatellites();
    QList<GLONASS> getGLONASSSatellites();
protected:
    void importNavigation(const QStringList &paths);   //输入导航电文
    void importObservation(const QStringList &paths); //输入观测文件
    void addSatellite(const QStringList &ephemerisData,const satelliteSystem sysName);
    QList<GPS> GSatellite;            //存储GPS卫星星历
    QList<BDS> CSatellite;            //存储BDS卫星星历
    QList<GALILEO> ESatellite;        //存储Galileo卫星星历
    QList<GLONASS> RSatellite;        //存储Glonass卫星星历
    QList<Observation> OFile;         //O文件
};

//overload operator "<<"  to print the ephemeris of each satellite
//QDebug  operator<<(QDebug out,GLONASS &S)
//{
//    S.printEphemeris();
//    return out;
//}
//QDebug  operator<<(QDebug out, GPS &S)
//{
//    S.printEphemeris();
//    return out;
//}
#endif // RINEXDATA_H
