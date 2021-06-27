#include "stellarmap.h"
StellarMap::StellarMap(QWidget *parent):QCustomPlot(parent)
{
    initStellarMap();
    timeTrigger=new QTimer(); //实例化触发器
    timeTrigger->start(1000); //每秒钟发送一次信号
    qDebug()<<"触发器的线程id："<<timeTrigger->thread()->currentThreadId();
    connect(timeTrigger,&QTimer::timeout,this,&StellarMap::showSatellites);
}
void StellarMap::initStellarMap()
{
    //设置坐标范围、坐标名称
    this->resize(this->width(), this->height());
    this->xAxis->setRange(-100.00, 100.00);
    this->yAxis->setRange(-100.00, 100.00);
    this->xAxis->setVisible(false); // 隐藏坐标轴
    this->yAxis->setVisible(false);
    //设置坐标轴名称
    //this->xAxis->setLabel("latitude");  //纬度
    // this->yAxis->setLabel("longitude"); //经度
    // 增加图层
    //this->addGraph(this->yAxis, this->xAxis);
    //设置坐标轴可见
    //this->xAxis->setVisible(true);
    //this->yAxis->setVisible(true);
    //保存为png图片
    //QString path = QApplication::applicationDirPath() + "\\" + "test123.png";
    //this->savePng(path);
    //实例化纬线
    latitudeLines=new QCPItemEllipse *[3];  //创建三条纬线
    for(int i=0;i<3;i++)
    {
        latitudeLines[i]=new  QCPItemEllipse(this);
        QPen pen(Qt::black,1);
        latitudeLines[i]->setPen(pen);
        latitudeLines[i]->topLeft->setCoords(-30*(i+1), 30*(i+1));      //左上角位置 tL
        latitudeLines[i]->bottomRight->setCoords(30*(i+1), -30*(i+1));  //右下角位置 bR
    }
    //实例化经线
    QPen longtitudePen(Qt::black, 1, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
    longtitudeLines=new QCPItemLine *[4];//创建4条经线
    longtitudeLines[0]=new QCPItemLine(this);
    longtitudeLines[0]->setPen(longtitudePen);
    longtitudeLines[0]->start->setCoords(-90/qSqrt(2),90/qSqrt(2));
    longtitudeLines[0]->end->setCoords(90/qSqrt(2), -90/qSqrt(2));
    longtitudeLines[1]=new QCPItemLine(this);
    longtitudeLines[1]->setPen(longtitudePen);
    longtitudeLines[1]->start->setCoords(-90/qSqrt(2),-90/qSqrt(2));
    longtitudeLines[1]->end->setCoords(90/qSqrt(2), 90/qSqrt(2));
    longtitudeLines[2]=new QCPItemLine(this);
    longtitudeLines[2]->setPen(longtitudePen);
    longtitudeLines[2]->start->setCoords(-90,0);
    longtitudeLines[2]->end->setCoords(90,0);
    longtitudeLines[3]=new QCPItemLine(this);
    longtitudeLines[3]->setPen(longtitudePen);
    longtitudeLines[3]->start->setCoords(0,90);
    longtitudeLines[3]->end->setCoords(0,-90);
    this->setInteractions( QCP::iSelectItems);
    /******************************
     * 若要设置可以拖拽、缩放、多选等，加上如下语句：
     * plot->setInteractions(QCP::iRangeDrag |QCP::iRangeZoom|QCP::iSelectPlottables| QCP::iMultiSelect)
     * **************************/
}
void StellarMap::showSatellites()
{
    satelliteTable.clear();//清空所有保存的坐标
    QDateTime current_date_time =QDateTime::currentDateTime();
    qDebug()<<current_date_time.toString("yyyy.MM.dd hh:mm:ss.zzz ddd");
    double time = DateToSeconds(current_date_time.toString("hh mm ss ddd")); //当前时刻
    /*************************************************************
     *********************开始遍历所有卫星***************************
     *************************************************************/
    //遍历GPS卫星
    if(!gps.empty())
    {
        QMap<QString, GPS>::iterator iter = gps.begin();
        while (iter != gps.end())
        {
            QVector<double> XYZ=iter->getSatellitePosition(time);
            QVector<double> LBH=iter->getSatelliteLBH(time);
           // qDebug()<<iter.key()<<XYZ[0]<<XYZ[1]<<XYZ[2]<<LBH[0]<<LBH[1]<<LBH[2];
            //录入当前时刻的坐标
            QList<double> list={XYZ[0],
                                XYZ[1],
                                XYZ[2],
                                NavigationObservationData::radianTodecimaldegree(LBH[0]),
                                NavigationObservationData::radianTodecimaldegree(LBH[1]),
                                LBH[2]};
            satelliteTable.insert(iter.key(),list);
            //先转换成10进制度再发过去
            SateDraws[iter.key()]->setPosition(NavigationObservationData::radianTodecimaldegree(LBH[0]),
                                               NavigationObservationData::radianTodecimaldegree(LBH[1]));
            iter++;
        }
    }
    if(!bds.empty())
    {
        QMap<QString, BDS>::iterator iter = bds.begin();
        while (iter != bds.end())
        {
            QVector<double> XYZ=iter->getSatellitePosition(time);
            QVector<double> LBH=iter->getSatelliteLBH(time);
           // qDebug()<<iter.key()<<XYZ[0]<<XYZ[1]<<XYZ[2]<<LBH[0]<<LBH[1]<<LBH[2];
            //录入当前时刻的坐标
            QList<double> list={XYZ[0],
                                XYZ[1],
                                XYZ[2],
                                NavigationObservationData::radianTodecimaldegree(LBH[0]),
                                NavigationObservationData::radianTodecimaldegree(LBH[1]),
                                LBH[2]};
            satelliteTable.insert(iter.key(),list);
            //先转换成10进制度再发过去
            SateDraws[iter.key()]->setPosition(NavigationObservationData::radianTodecimaldegree(LBH[0]),
                                               NavigationObservationData::radianTodecimaldegree(LBH[1]));
            iter++;
        }
    }
    replot();
    emit tableSignal(satelliteTable);//发送卫星坐标表
}
double StellarMap::DateToSeconds(const QString &date)//星期几转换成该周的第几秒
{
    QStringList str=date.split(" ");
    double hh,mm,ss,ddd=0.0;   //时，分，秒，星期几
    int multiple=0;
    hh=str[0].toDouble()*3600;
    mm=str[1].toDouble()*60;
    ss=str[2].toDouble();
    if(str[3]=="周一")
    {
        multiple=0;
    }else if(str[3]=="周二")
    {
        multiple=1;
    }
    else if(str[3]=="周三")
    {
        multiple=2;
    }else if(str[3]=="周四")
    {
        multiple=3;
    }else if(str[3]=="周五")
    {
        multiple=4;
    }else if(str[3]=="周六")
    {
        multiple=5;
    }else if(str[3]=="周日")
    {
        multiple=6;
    }
    ddd=multiple*3600*24;
    return hh+mm+ss+ddd;
}
void StellarMap::addSatellite(const GPS &sate) //增加或更新GPS卫星
{
    if(!gps.contains(sate.sateEphemeris.PRN))
    {
        gps.insert(sate.sateEphemeris.PRN,sate);
        SatellitePicture *picture=new SatellitePicture(this,sate.sateEphemeris.PRN);

        SateDraws.insert(sate.sateEphemeris.PRN,picture);  //添加卫星的形状

    }else
    {
        gps.remove(sate.sateEphemeris.PRN);
        gps.insert(sate.sateEphemeris.PRN,sate);
    }
}
void StellarMap::addSatellite(const BDS &sate)//增加BDS卫星
{
        if(!bds.contains(sate.sateEphemeris.PRN))
        {
            bds.insert(sate.sateEphemeris.PRN,sate);
            SatellitePicture *picture=new SatellitePicture(this,sate.sateEphemeris.PRN);
            SateDraws.insert(sate.sateEphemeris.PRN,picture);  //添加卫星的形状
        }else
        {
            bds.remove(sate.sateEphemeris.PRN);
            bds.insert(sate.sateEphemeris.PRN,sate);
        }
}
void StellarMap::addSatellite(const GLONASS &sate) //增加GLONASS卫星
{
}
void StellarMap::addSatellite(const GALILEO &sate) //增加GALILEO卫星
{
}
StellarMap::~StellarMap()
{
    delete latitudeLines;
    delete longtitudeLines;
    timeTrigger->deleteLater();
    qDebug()<<"星空图类被销毁";
}
/*******************************
 ***************卫星形状类********
 ******************************/

SatellitePicture::SatellitePicture(QCustomPlot *parentPlot, const QString &sateName):
    QCPItemEllipse(parentPlot),PRN(sateName)
{
    this->setSelectable(true); //设置该图形可选
    qDebug()<<"添加了"<<PRN;
    text=new QCPItemText(parentPlot);
    //不同卫星设置不同的颜色
    if(PRN.contains("G"))
    {
        this->setBrush(QBrush(QColor(255,100,100)));
    }else if(PRN.contains("C"))
    {
        this->setBrush(QBrush(QColor(100,200,200)));
    }
    text->setText(PRN); //设置卫星名称
    QFont font;//设置字体
    font.setFamily("微软雅黑");//字体
    font.setPixelSize(12);//文字像素大小
    font.setBold(true);//粗体
    // font.setStyle(QFont::StyleOblique);
    // font.setLetterSpacing(QFont::PercentageSpacing,200);//间距
    text->setFont(font);
    text->setVisible(false);
    text->setSelectable(false);//设置文本不可选
    //font.setPointSize(20);//文字大小
    //font.setUnderline(true);//下划线
    //font.setStrikeOut(true);//中划线
    //font.setOverline(true);//上划线
    //font.setItalic(true);//斜体
    // font.setCapitalization(QFont::Capitalize);//首字母大写
}
void SatellitePicture::setPosition(const double &L,const double &B) //十进制度
{
    double x,y;
    x=-(90-qAbs(B))*qCos(PI*L/180);
    y=-(90-qAbs(B))*qSin(PI*L/180);
    //设置卫星在坐标系的位置
    this->topLeft->setCoords(x-3,y+3);
    this->bottomRight->setCoords(x+3,y-3);

    //设置文本在坐标系的位置
    text->setVisible(true);
    text->position->setCoords(x,y);//设置文本位置
    //    textLabel = new QCPItemText(this);//在QCustomplot中新建文字框
    //    textLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);//文字布局：顶、左对齐
    //    textLabel->position->setType(QCPItemPosition::ptAxisRectRatio);//位置类型（当前轴范围的比例为单位/实际坐标为单位）
    //    textLabel->position->setCoords(0.5, 0); //把文字框放在X轴的中间，Y轴的最顶部
    //    textLabel->setText("Text Item Demo");
    //    textLabel->setFont(QFont(font().family(), 16)); //字体大小
    //    textLabel->setPen(QPen(Qt::black)); //字体颜色
    //    textLabel->setPadding(QMargins(2,2,2,2));//文字距离边框几个像素
}
SatellitePicture:: ~SatellitePicture()
{
    qDebug()<<"卫星形状类销毁！";
}

