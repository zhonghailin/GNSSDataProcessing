/******************************************************************************
* @author： Harry Chung
* @date：2021-06-15 08:56:35
* @fileName：rinexdata.cpp
* @brief：
*******************************************************************************/
#include "rinexdata.h"
/*
 * 类内静态成员不能在类内初始化，必须在类外,因为无论你在外部实例化多少个相同的类，
 * 静态成员必须只有一个对象（类内常量成员也一样）。但静态常量成员可以类内初始化
 */
NavigationObservationData *NavigationObservationData::m_pInstance = nullptr;
Satellite::Satellite()
{
}
QVector<double> Satellite::getSatelliteLBH(const double &t)//返回的是弧度
{
    //注意e1s与e2s都已经是平方了
    QVector<double> p=getSatellitePosition(t);
    double X=p[0],Y=p[1],Z=p[2];  //空间坐标
    double B,L,H;                 //纬度、经度、高程
    double N;                     //卯酉圈曲率半径
    //迭代计算纬度
    double B1;
    B=qAtan(Z / qSqrt(X * X + Y * Y)); //设初值
    do{
      B1=B;
      B=qAtan((Z + (a / qSqrt(1 - e1s * qSin(B1) * qSin(B1))) * e1s * qSin(B1)) / qSqrt(X * X + Y * Y));
    }while(qAbs(B-B1)>0.0000000000000001);
    //计算经度
     L=qAtan(Y/X);
     if(X<0&&Y>=0)
    {
      L=L+PI;
    }else if(X<0&&Y<0)
     {
      L=L-PI;
     }
    N=a/qSqrt(1-e1s*qSin(B)*qSin(B)); //计算卯酉圈的曲率半径
    //H=Z/qSin(B)-N*(1-e1s);
    H = qSqrt(X * X + Y * Y) /qCos(B) - N;
    QVector<double> LBH;
    LBH.append(L);
    LBH.append(B);
    LBH.append(H);
    return LBH;
}
QVector<double> Satellite::LBHto3DCoordinate(const QVector<double> &LBH)
{
    double L=LBH[0],B=LBH[1],H=LBH[2];
    double X,Y,Z;
    double N=a/qSqrt(1-e1s*qSin(B)*qSin(B));
    X=(N+H)*qCos(B)*qCos(L);
    Y=(N+H)*qCos(B)*qSin(L);
    Z=(N*(1-e1s)+H)*qSin(B);
    qDebug()<<QString::number(X,'f',4)<<QString::number(Y,'f',4)<<QString::number(Z,'f',4);
    QVector<double> XYZ={X,Y,Z};
    return XYZ;
}
Satellite::~Satellite()
{
    qDebug() << "the Satellite base class was destroyed!";
}
/*==========================================================
 *==========================================================GLONASS
 *==========================================================*/
 GLONASS::GLONASS(const QStringList &ephemerisData):Satellite ()
 {
     //初始化参数列表(PZ90.02椭球参数)
     GM=398600441800000;
     w_e=0.00007292115;
     a=6378136;
     e1s=0.006694366177481926;  //第一偏心率的平方
     e2s=0.006739496742276433;  //第二偏心率的平方
    //导入星历参数
     updateEphemeris(ephemerisData);
 }
void GLONASS::updateEphemeris(const QStringList &List)
{
    if(List.length()!=4) return;
    QString str1=List[0];
    QString str2=List[1];
    QString str3=List[2];
    QString str4=List[3];
    //replace 'D' by 'E'
    str1.replace('D','E');
    str2.replace('D','E');
    str3.replace('D','E');
    str4.replace('D','E');
    sateEphemeris={str1.mid(0,3),str1.mid(4,19),           str1.mid(23,19).toDouble(),str1.mid(42,19).toDouble(),str1.mid(61,19).toDouble(),
                                 str2.mid(4,19).toDouble(),str2.mid(23,19).toDouble(),str2.mid(42,19).toDouble(),str2.mid(61,19).toDouble(),
                                 str3.mid(4,19).toDouble(),str3.mid(23,19).toDouble(),str3.mid(42,19).toDouble(),str3.mid(61,19).toDouble(),
                                 str4.mid(4,19).toDouble(),str4.mid(23,19).toDouble(),str4.mid(42,19).toDouble(),str4.mid(61,19).toDouble()};
}
QVector<double> GLONASS::getSatellitePosition(const double &t)
{
    QVector<double> b;
    b.append(1);
    return b;
}
void GLONASS::printEphemeris()
{
    Ephemeris s=sateEphemeris;
    qDebug()<<s.PRN<<s.epochTime<<s.r_n<<s.tau_n<<s.fTime<<endl
           <<"    "<<s.x<<s.x_v<<s.x_a<<s.svhlth<<endl
           <<"    "<<s.y<<s.y_v<<s.y_a<<s.freNum<<endl
           <<"    "<<s.z<<s.z_v<<s.z_a<<s.ageLimit<<endl;
}
GLONASS::~GLONASS()
{
    qDebug()<<"there is a glonass satellite was destroyed!";
}
/*==========================================================
 *========================GPS================================
 *==========================================================*/
GPS::GPS(const QStringList &ephemerisData):Satellite()
{
    //初始化形参列表(WGS-84椭球参数)
    GM=398600470000000;
    w_e=0.00007292115;
    a=6378137;
    e1s=0.00669437999013; //第一偏心率的平方
    e2s=0.00673949674227; //第二偏心率的平方
    //导入星历
    updateEphemeris(ephemerisData);
}
/******************************************
* brief:
******************************************/
void GPS::updateEphemeris(const QStringList &List)
{
    if(List.length()!=8) return;
    QString str1=List[0];
    QString str2=List[1];
    QString str3=List[2];
    QString str4=List[3];
    QString str5=List[4];
    QString str6=List[5];
    QString str7=List[6];
    QString str8=List[7];
    //replace 'D' by 'E'
    str1.replace('D','E');
    str2.replace('D','E');
    str3.replace('D','E');
    str4.replace('D','E');
    str5.replace('D','E');
    str6.replace('D','E');
    str7.replace('D','E');
    str8.replace('D','E');
    sateEphemeris={str1.mid(0,3),str1.mid(4,19),           str1.mid(23,19).toDouble(),str1.mid(42,19).toDouble(),str1.mid(61,19).toDouble(),
                                 str2.mid(4,19).toDouble(),str2.mid(23,19).toDouble(),str2.mid(42,19).toDouble(),str2.mid(61,19).toDouble(),
                                 str3.mid(4,19).toDouble(),str3.mid(23,19).toDouble(),str3.mid(42,19).toDouble(),str3.mid(61,19).toDouble(),
                                 str4.mid(4,19).toDouble(),str4.mid(23,19).toDouble(),str4.mid(42,19).toDouble(),str4.mid(61,19).toDouble(),
                                 str5.mid(4,19).toDouble(),str5.mid(23,19).toDouble(),str5.mid(42,19).toDouble(),str5.mid(61,19).toDouble(),
                                 str6.mid(4,19).toDouble(),str6.mid(23,19).toDouble(),str6.mid(42,19).toDouble(),str6.mid(61,19).toDouble(),
                                 str7.mid(4,19).toDouble(),str7.mid(23,19).toDouble(),str7.mid(42,19).toDouble(),str7.mid(61,19).toDouble(),
                                 str8.mid(4,19).toDouble(),str8.mid(23,19).toDouble()};
}
/******************************************
* brief:计算并获取卫星的坐标
******************************************/
QVector<double> GPS::getSatellitePosition(const double &t)
{
    Ephemeris &E=sateEphemeris;
    double a_0=E.a_0, a_1=E.a_1, a_2=E.a_2;
    double C_rs=E.C_rs, delta_n=E.delta_n, M_0=E.M_0;
    double C_uc=E.C_uc, e=E.e, C_us=E.C_us,  sqrt_a=E.sqrt_a;
    double t_oe=E.t_oe, C_ic=E.C_ic, Omega_0=E.Omega_0, C_is=E.C_is;
    double i_0=E.i_0, C_rc=E.C_rc, omega=E.omega, Omega_p=E.Omega_p;
    double i_p=E.i_p;
    double a,n_0,n,M_k,E_k, f_k,Fai_k,sigma_u_k,sigma_r_k,sigma_i_k,u_k,r_k,i_k,Omega_k;
    double x_k,y_k;
    double X_k,Y_k,Z_k;
    double t_k;
    // (1)求轨道长半轴
    a=sqrt_a*sqrt_a;
    // (2)calculate the mean angular velocity(计算平均角速度)
    n_0=sqrt(GM/pow(a,3));
    //(3)calculate the time difference between t to t_oe
    t_k=t-t_oe;
    double delta_t=a_0+a_1*t_k+a_2*t_k*t_k;  //计算钟差改正数
    t_k=t_k-delta_t;   //钟差改正
    //(4)correct the mean angular velocity n(平均角速度改正)
    n=n_0+delta_n;
    //(5)calculate the Mean Anomaly(计算平近点角)
    M_k=M_0+n*t_k;
    //(6)calculate the Eccentric Anomaly,its unit is rad(计算偏见近点角，单位为弧度)
    E_k=M_k;
    double Etemp;
    do{
         Etemp=E_k;
         E_k=M_k+e*qSin(Etemp);
    }while(qAbs(E_k-Etemp)>qPow(10,-12)); //iteration
    //(7)calculate the trueanomaly(计算真近点角)
    f_k=acos((cos(E_k)-e)/(1-e*cos(E_k)));
    //(8)calculate the argument of latitude(计算升交角距)
    Fai_k=f_k+omega;
    //(9)Calculate the satellite orbit perturbation correction(计算卫星轨道摄动项改正数)
    sigma_u_k=C_us*sin(2*Fai_k)+C_uc*cos(2*Fai_k);
    sigma_r_k=C_rs*sin(2*Fai_k)+C_rc*cos(2*Fai_k);
    sigma_i_k=C_is*sin(2*Fai_k)+C_ic*cos(2*Fai_k);
    //(10)Calculate the corrected trueanomaly,radius,orbit inclination
    u_k=Fai_k+sigma_u_k;             //改正后的升交角距
    r_k=a*(1-e*cos(E_k))+ sigma_r_k; //改正后的向径
    i_k=i_0+sigma_i_k+i_p*t_k;       //改正后的轨道倾角
    //(11)calculate the instantaneous argument of latitude of satellite(计算观测时刻的升交点精度)
    Omega_k=Omega_0+(Omega_p-w_e)*t_k-w_e*t_oe;
    //(12)calculate the coordinates of a satellite in the orbital plane(计算卫星在轨道平面内的坐标)
    x_k=r_k*cos(u_k);
    y_k=r_k*sin(u_k);
    //(13)calculate the space rectangular coordinate of a satellite in the Earth-Fixed Coordinate System
    X_k=x_k*cos(Omega_k)-y_k*cos(i_k)*sin(Omega_k);
    Y_k=x_k*sin(Omega_k)+y_k*cos(i_k)*cos(Omega_k);
    Z_k=y_k*sin(i_k);
    QVector<double> coordinate;
    coordinate.append(X_k);
    coordinate.append(Y_k);
    coordinate.append(Z_k);
    return coordinate;
}
void GPS::printEphemeris()
{
    Ephemeris s=sateEphemeris;
    qDebug()<<s.PRN<<s.epochTime<<s.a_0<<s.a_1<<s.a_2<<endl
            <<"    "<<s.IODE<<s.C_rs<<s.delta_n<<s.M_0<<endl
            <<"    "<<s.C_uc<<s.e<<s.C_us<<s.sqrt_a<<endl
            <<"    "<<s.t_oe<<s.C_ic<<s.Omega_0<<s.C_is<<endl
            <<"    "<<s.i_0<<s.C_rc<<s.omega<<s.Omega_p<<endl
            <<"    "<<s.i_p<<s.cflgl2<<s.weekno<<s.pflgl2<<endl
            <<"    "<<s.svacc<<s.svhlth<<s.tgd<<s.IODC<<endl
            <<"    "<<s.ttm<<s.fi<<endl;
}
GPS::~GPS()
{
     qDebug()<<"there is a GPS satellite was destroyed!";
}
/*===========================================================
 *==========================BDS==============================
 *==========================================================*/
BDS::BDS(const QStringList &ephemerisData):GPS(ephemerisData)
{
    //初始化形参列表(CGCS2000椭球参数)
    GM=398600441800000;
    w_e=0.00007292115;
    a=6378137;
    e1s=0.00669438002290;
    e2s=0.00673949677548;
}
QVector<double> BDS::getSatellitePosition(const double &t)
{
    Ephemeris &E=sateEphemeris;
    double a_0=E.a_0, a_1=E.a_1, a_2=E.a_2;
    double C_rs=E.C_rs, delta_n=E.delta_n, M_0=E.M_0;
    double C_uc=E.C_uc, e=E.e, C_us=E.C_us,  sqrt_a=E.sqrt_a;
    double t_oe=E.t_oe, C_ic=E.C_ic, Omega_0=E.Omega_0, C_is=E.C_is;
    double i_0=E.i_0, C_rc=E.C_rc, omega=E.omega, Omega_p=E.Omega_p;
    double i_p=E.i_p;
    double a,n_0,n,M_k,E_k, f_k,Fai_k,sigma_u_k,sigma_r_k,sigma_i_k,u_k,r_k,i_k,Omega_k;
    double x_k,y_k,z_k;
    double X_k,Y_k,Z_k;
    double t_k;
    QVector<double> coordinate;  //返回的坐标结果
    // (1)求轨道长半轴
    a=sqrt_a*sqrt_a;
    // (2)calculate the mean angular velocity(计算平均角速度)
    n_0=sqrt(GM/pow(a,3));
    //(3)calculate the time difference between t to t_oe
    t_k=t-t_oe;
    double delta_t=a_0+a_1*t_k+a_2*t_k*t_k;  //计算钟差改正数
    t_k=t_k-delta_t;   //钟差改正
    //(4)correct the mean angular velocity n(平均角速度改正)
    n=n_0+delta_n;
    //(5)calculate the Mean Anomaly(计算平近点角)
    M_k=M_0+n*t_k;
    //(6)calculate the Eccentric Anomaly,its unit is rad(计算偏见近点角，单位为弧度)
    E_k=M_k;
    double Etemp;
    do{
         Etemp=E_k;
         E_k=M_k+e*qSin(Etemp);

    }while(qAbs(E_k-Etemp)>qPow(10,-12)); //iteration
    //(7)calculate the trueanomaly(计算真近点角)
    f_k=acos((cos(E_k)-e)/(1-e*cos(E_k)));
    //(8)calculate the argument of latitude(计算升交角距)
    Fai_k=f_k+omega;
    //(9)Calculate the satellite orbit perturbation correction(计算卫星轨道摄动项改正数)
    sigma_u_k=C_us*sin(2*Fai_k)+C_uc*cos(2*Fai_k);
    sigma_r_k=C_rs*sin(2*Fai_k)+C_rc*cos(2*Fai_k);
    sigma_i_k=C_is*sin(2*Fai_k)+C_ic*cos(2*Fai_k);
    //(10)Calculate the corrected trueanomaly,radius,orbit inclination
    u_k=Fai_k+sigma_u_k;             //改正后的升交角距
    r_k=a*(1-e*cos(E_k))+ sigma_r_k; //改正后的向径
    i_k=i_0+sigma_i_k+i_p*t_k;       //改正后的轨道倾角
    //(11)calculate the coordinates of a satellite in the orbital plane(计算卫星在轨道平面内的坐标)
    x_k=r_k*cos(u_k);
    y_k=r_k*sin(u_k);
    z_k=0;
    if(MEOPRN.contains(sateEphemeris.PRN)||IGSOPRN.contains(sateEphemeris.PRN))  //如果是IGSO或MEO卫星则按如下方法计算
    {
    //(12)calculate the instantaneous argument of latitude of satellite(计算观测时刻的升交点精度)
    Omega_k=Omega_0+(Omega_p-w_e)*t_k-w_e*t_oe;
    //(13)calculate the space rectangular coordinate of a satellite in the Earth-Fixed Coordinate System
    X_k=x_k*cos(Omega_k)-y_k*cos(i_k)*sin(Omega_k);
    Y_k=x_k*sin(Omega_k)+y_k*cos(i_k)*cos(Omega_k);
    Z_k=y_k*sin(i_k);
    coordinate.append(X_k);
    coordinate.append(Y_k);
    coordinate.append(Z_k);
    }else if(GEOPRN.contains(sateEphemeris.PRN))  //如果是GEO卫星则按如下方法计算
    {
        double X_Gk,Y_Gk,Z_Gk;
        Omega_k=Omega_0+Omega_p*t_k-w_e*t_oe;
        //(13)calculate the space rectangular coordinate of a satellite in the Earth-Fixed Coordinate System
        X_Gk=x_k*cos(Omega_k)-y_k*cos(i_k)*sin(Omega_k);
        Y_Gk=x_k*sin(Omega_k)+y_k*cos(i_k)*cos(Omega_k);
        Z_Gk=y_k*sin(i_k);
        Matrix<double,3,1> XYZ_Gk;
        Matrix<double,3,1> XYZ_k;
        Matrix<double,3,3> R_Z;
        Matrix<double,3,3> R_X;
        XYZ_Gk<<X_Gk,
                Y_Gk,
                Z_Gk;
        R_Z<< cos(w_e*t_k) ,  sin(w_e*t_k), 0,
              -sin(w_e*t_k),  cos(w_e*t_k), 0,
              0            ,             0, 1;
        R_X<<1,0            ,0           ,
             0,cos(PI/36) ,sin(PI/36),
             0,-sin(PI/36),cos(PI/36);
        XYZ_k= R_Z* R_X*XYZ_Gk;   //将坐标系旋转回来，使椭圆面与地球赤道面重合
        coordinate.append(XYZ_k[0]);
        coordinate.append(XYZ_k[1]);
        coordinate.append(XYZ_k[2]);
        return coordinate;
    }
    return coordinate;
}
bool BDS::isValid()
{

    if(MEOPRN.contains(sateEphemeris.PRN)
            ||IGSOPRN.contains(sateEphemeris.PRN)
            ||GEOPRN.contains(sateEphemeris.PRN))
    {
        return true;
    }
        return false;
}
BDS::~BDS()
{
    qDebug()<<"there is a BDS satellite was destroyed!";
}
/*==========================================================
 *==========================================================GALILEO
 *==========================================================*/
GALILEO::GALILEO(const QStringList &ephemerisData):GPS(ephemerisData){}
GALILEO::~GALILEO()
{
    qDebug()<<"there is a galileo satellite was destroyed!";
}
/*==========================================================
 *==========================================================NavigationObservationData
 *==========================================================*/
NavigationObservationData::NavigationObservationData()
{

}
void NavigationObservationData::importFiles(const QStringList &paths)
{
     QStringList navigationPaths;  //筛选出的导航电文的路径
     QStringList observationPaths; //筛选出的观测文件的路径
     QRegExp rexObservation(".*\\.[0-9]{2}O");       //匹配观测文件名
     QRegExp rexNavigation(".*\\.[0-9]{2}[G,L,C,N]"); //匹配导航电文
     for(int i=0;i<paths.length();i++)
     {
        rexNavigation.indexIn(paths[i]);
        rexObservation.indexIn(paths[i]);
       if(rexNavigation.matchedLength()!=-1)
        {
         navigationPaths.append(paths[i]);
        }else if(rexObservation.matchedLength()!=-1)
        {
            observationPaths.append(paths[i]);
        }
      }


     importNavigation(navigationPaths);   //输入导航电文
     importObservation(observationPaths); //输入观测文件
}
void NavigationObservationData::addSatellite(const QStringList &ephemerisData, const satelliteSystem sysName)
{
    if(sysName==NavigationObservationData::Gps)
    {
        GPS Sate(ephemerisData);
        GSatellite.append(Sate);
    }else if(sysName==NavigationObservationData::Bds)
    {
        BDS Sate(ephemerisData);
        if(Sate.isValid())  //如果卫星有效，则添加
           CSatellite.append(Sate);
    }else if(sysName==NavigationObservationData::Galileo)
    {
        GALILEO Sate(ephemerisData);
        ESatellite.append(Sate);
    }else if (sysName==NavigationObservationData::Glonass)
    {
        GLONASS Sate(ephemerisData);
        RSatellite.append(Sate);
    }

}
void NavigationObservationData::importNavigation(const QStringList &paths)
{
    for (int i=0;i<paths.length();i++)
    {
        //begin to read text
        QStringList linList;
        QFile file(paths[i]);
        if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            qDebug()<<"Can't open the file!"<<endl;
            return;
        }
        while(!file.atEnd())
        {
            QString line = file.readLine();
            linList.append(line);
        }
        file.close();
        //begin to import ephemeris Data
        if(linList[0].contains("GPS")||linList[0].contains("GALILEO")||linList[0].contains("BDS"))
        {
            QRegExp rx("[C,G,E]{1}[0-9]{2} .+");
            for (int i=0;i<linList.length();i++)
            {
                rx.indexIn(linList[i]);
                if(rx.matchedLength()!=-1)
                {
                    qDebug()<<linList[i].mid(0,3);
                    QStringList temp;
                    temp.append(linList[i]);
                    temp.append(linList[i+1]);
                    temp.append(linList[i+2]);
                    temp.append(linList[i+3]);
                    temp.append(linList[i+4]);
                    temp.append(linList[i+5]);
                    temp.append(linList[i+6]);
                    temp.append(linList[i+7]);
                    if(linList[0].contains("GPS")) addSatellite(temp,NavigationObservationData::Gps);
                    else if(linList[0].contains("GALILEO")) addSatellite(temp,NavigationObservationData::Galileo);
                    else if(linList[0].contains("BDS")) addSatellite(temp,NavigationObservationData::Bds);
                }
            }
        }else if(linList[0].contains("GLONASS"))
        {
            QRegExp rx("R[0-9]{2} .+");
            for (int i=0;i<linList.length();i++)
            {
                rx.indexIn(linList[i]);
                if(rx.matchedLength()!=-1)
                {
                    qDebug()<<linList[i].mid(0,3);
                    QStringList temp;
                    temp.append(linList[i]);
                    temp.append(linList[i+1]);
                    temp.append(linList[i+2]);
                    temp.append(linList[i+3]);
                    addSatellite(temp,NavigationObservationData::Glonass);
                }
            }
        }
    }
}
void NavigationObservationData::importObservation(const QStringList &paths)
{
    for (int i=0;i<paths.length();i++)
    {
        //begin to read text
        QStringList linList;
        QFile file(paths[i]);
        if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            qDebug()<<"Can't open the file!"<<endl;
            return;
        }
        while(!file.atEnd())
        {
            QString line = file.readLine();
            linList.append(line);
        }
        file.close();
        Observation w(linList);
        OFile.append(w);
    }
}
double NavigationObservationData::radianTodecimaldegree(const double &rad)
{
    return 180*rad/PI;
}
QList<GPS> NavigationObservationData::getGPSSatellites()
{
    return GSatellite;
}
QList<BDS> NavigationObservationData::getBDSSatellites()
{
    return CSatellite;
}
QList<GALILEO> NavigationObservationData::getGALILEOSatellites()
{
    return ESatellite;
}
QList<GLONASS> NavigationObservationData::getGLONASSSatellites()
{
    return RSatellite;
}
void NavigationObservationData::printObservation()
{


}
QVector<double> NavigationObservationData::getSatellitePsition(const QString &satelliteName,const QString &epoch)
{
//   if(GSatellite.length()!=0)
//    {
//       GPS sate=GSatellite[0];

//       double t=sate.sateEphemeris.t_oe;
//       for(int i=0;i<100;i++)
//       {
//           QVector<double> p=sate.getSatellitePosition(t+i);
//           QVector<double> LBH=sate.getSatelliteLBH(t+i);
//           QVector<double> XYZ=sate.LBHto3DCoordinate(LBH);
//           qDebug()<<QString::number(p[0],'f',4)
//                   <<QString::number(p[1],'f',4)
//                   <<QString::number(p[2],'f',4)
//                   <<QString::number(radianTodecimaldegree(LBH[0]),'f',4)
//                   <<QString::number(radianTodecimaldegree(LBH[1]),'f',4)
//                   <<QString::number(LBH[2],'f',4)
//                   <<QString::number(XYZ[0],'f',4)
//                   <<QString::number(XYZ[1],'f',4)
//                   <<QString::number(XYZ[2],'f',4);
//       }
//   }
//    QVector<double> p={0,0,0};
//    return p;
}
void NavigationObservationData::printEphemeris(const satelliteSystem systemName)
{
//    if(systemName==NavigationObservationData::Gps)
//    {
//        for (int i=0;i<GSatellite.length();i++)
//        {
//            qDebug()<<GSatellite[i];
//        }
//    }else if(systemName==NavigationObservationData::Bds)
//    {
//        for (int i=0;i<CSatellite.length();i++)
//        {
//            qDebug()<<CSatellite[i];
//        }
//    }else if(systemName==NavigationObservationData::Galileo)
//    {
//        for (int i=0;i<ESatellite.length();i++)
//        {
//          qDebug()<< ESatellite[i];
//        }
//    }else if(systemName==NavigationObservationData::Glonass)
//    {
//        for (int i=0;i<RSatellite.length();i++)
//        {
//           qDebug()<< RSatellite[i];
//        }
//    }
}
NavigationObservationData::~NavigationObservationData()
{
    qDebug()<<"NavigationObservationData class had been deleted!" ;
}

/***************************************************
 ****************************************************
 *******************观测文件类*************************
 ****************************************************/
Observation::Observation(const QStringList &List)
{
    importFile(List);
}
void Observation::importFile(const QStringList &List)
{
    int num=0;  //代表第几行
    int length=List.length();//列表行数
    qDebug()<<"总行数："<<length;
    /***************************************
     * ***************读取表头****************
     * *************************************/
    for(;num<length;num++)
    {
        if(List[num].contains("END OF HEADER")) break;  //到表头时退出循环
        //读取观测类型 (注意：该算法没有考虑到续航情况)
        if(List[num].contains("SYS / # / OBS TYPES"))
        {
            obsTypes type;
            type.sateName=List[num][0];
            type.typeNum=List[num].mid(3,3).toShort();
            QList<QString> list;
            QStringList str=List[num].mid(6,54).split(' '); //空格分割
            for(int i=0;i<str.length();i++)
            {
                if(str[i].length()>=2)
                   list.append(str[i]);
            }
            type.types=list;
            OBSTYPES.append(type);
        }else if(List[num].contains("APPROX POSITION XYZ"))
        {
           //读取接收机近似坐标
          APPROXPOSITIONXYZ[0]=List[num].mid(0,14).toDouble();
          APPROXPOSITIONXYZ[1]=List[num].mid(15,13).toDouble();
          APPROXPOSITIONXYZ[2]=List[num].mid(29,13).toDouble();
        }else if(List[num].contains("MARKER NAME"))
        {
            //读取测站点名
            MARKERNAME=List[num].mid(0,4);
        }else if(List[num].contains("ANTENNA: DELTA H/E/N"))
        {
            //读取接收机天线相位中心偏差
            ANTENNADELTAHEN[0]=List[num].mid(0,14).toDouble();
            ANTENNADELTAHEN[1]=List[num].mid(15,13).toDouble();
            ANTENNADELTAHEN[2]=List[num].mid(29,13).toDouble();
        }else if(List[num].contains("INTERVAL"))
        {
            //读取历元间隔
            INTERVAL= List[num].mid(4,6).toDouble();
        }else if(List[num].contains("INTERVAL"))
        {
            //读取跳秒数
             LEAPSECONDS= List[num].mid(4,2).toFloat();
        }else if(List[num].contains("RCV CLOCK OFFS APPL"))
        {
            //读取跳秒数
             LEAPSECONDS=bool(List[num].mid(5,1).toInt());
        }else if(List[num].contains("TIME OF FIRST OBS"))
        {
             //读取开始观测的时间
             TIMEOFFIRSTOBS[0]=List[num].mid(44,7).replace(" ","");
             TIMEOFFIRSTOBS[1]=List[num].mid(2,4);
             TIMEOFFIRSTOBS[2]=List[num].mid(10,2);
             TIMEOFFIRSTOBS[3]=List[num].mid(16,2);
             TIMEOFFIRSTOBS[4]=List[num].mid(22,2);
             TIMEOFFIRSTOBS[5]=List[num].mid(28,2);
             TIMEOFFIRSTOBS[6]=List[num].mid(33,10);
        }else if(List[num].contains("TIME OF LAST OBS"))
        {
             //读取结束观测的时间
            TIMEOFLASTOBS[0]=List[num].mid(44,7).replace(" ","");
            TIMEOFLASTOBS[1]=List[num].mid(2,4);
            TIMEOFLASTOBS[2]=List[num].mid(10,2);
            TIMEOFLASTOBS[3]=List[num].mid(16,2);
            TIMEOFLASTOBS[4]=List[num].mid(22,2);
            TIMEOFLASTOBS[5]=List[num].mid(28,2);
            TIMEOFLASTOBS[6]=List[num].mid(33,10);
        }else if(List[num].contains("SYS / PHASE SHIFT"))
        {
            //读取相位改正数
            for(int j=0;j<OBSTYPES.length();j++)
            {
                if(OBSTYPES[j].sateName==List[num][0])
                {
                    OBSTYPES[j].phaseShift.insert(List[num].mid(2,3),            //观测类型
                                                  List[num].mid(6,8).toDouble());//相位改正值
                    break;
                }
            }
        }
    }
    qDebug()<<MARKERNAME;
    for (int i=0;i<OBSTYPES.length();i++)
    {
        qDebug()<<OBSTYPES[i].sateName<<OBSTYPES[i].typeNum<<OBSTYPES[i].types;
    }
    /***************************************
     ****************读取历元****************
     ***************************************/
    while(num<length)
    {
        //读取第一行，格式为：A1,1X,I4,4(1X,I2),F11.7,2X,I1,I3,6X,F15.12
        if(List[num][0]==">")
        {
            epoch p;  //创建一个历元
            p.time[0]=List[num].mid(2,4);  //年
            p.time[1]=List[num].mid(7,2);  //月
            p.time[2]=List[num].mid(10,2); //日
            p.time[3]=List[num].mid(13,2); //时
            p.time[4]=List[num].mid(16,2); //分
            p.time[5]=List[num].mid(18,11);//秒
            p.Marker=List[num].mid(31,1).toShort();        //历元标志
            p.sateNum=List[num].mid(32,3).toShort();       //卫星数量
            p.clockOffset=List[num].mid(41,15).toDouble(); //接收机钟差
            ++num;
            while(num<length)
            {
               if(List[num][0]=='>')break;
               QString prn;  //卫星PRN号
               short typenum=0;//观测值数量
               QList<QList<double>> OBS;
               for(int i=0;i<OBSTYPES.length();i++)
               {
                   if(OBSTYPES[i].sateName==List[num][0])
                   {
                       prn=List[num].mid(0,3);
                       typenum=OBSTYPES[i].typeNum;
                       break;
                   }
               }
              int  temp=3;
              for(short i=0;i<typenum;i++)
              {
               QList<double> obs;
               obs.append(List[num].mid(temp,14).toDouble());   //观测值
               obs.append(List[num].mid(temp+14,1).toDouble());//LLI
               obs.append(List[num].mid(temp+15,1).toDouble()); //信号强度
               OBS.append(obs);
               temp+=16;
              }
              if(prn[0]==("G"))
              {
                p.G.insert(prn,OBS);
              }else if(prn[0]==("C"))
              {
                p.C.insert(prn,OBS);
              }else if(prn[0]==("R"))
              {
                 p.R.insert(prn,OBS);
              }else if(prn[0]==("E"))
              {
                 p.E.insert(prn,OBS);
              }
              ++num;
            }
            EPOCH.append(p);//加入该历元
        }else
        {
           ++num;
        }
    }
    qDebug()<<"完成读取O文件";
    qDebug()<<"历元数："<<EPOCH.length();
}
Observation::~Observation()
{
    qDebug()<<"观测值类销毁！";

}
