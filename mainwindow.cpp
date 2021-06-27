#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QTime>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    tableWid=new SatePositionTable();
    ui->SateTableDockWid->setWidget(tableWid);
     connect(ui->stellarMapWid,&StellarMap::tableSignal,tableWid,&SatePositionTable::updateTable);
    qDebug()<<"主窗口的线程id:"<< QThread::currentThreadId();
}
MainWindow::~MainWindow()
{
    delete ui;
}
//打开文件
void MainWindow::on_open_file_action_triggered()
{
    QStringList paths= QFileDialog::getOpenFileNames(this,"open Rinex files","..");
    data=NavigationObservationData::GetInstance();
    data->importFiles(paths);
    data->printEphemeris(NavigationObservationData::Gps);
    QList<GPS> Gsates=data->getGPSSatellites(); //获取GPS卫星
    QList<BDS> Csates=data->getBDSSatellites(); //获取BDS卫星
    for (int i=0;i<Gsates.length();i++)   //将GPS卫星输入星空图
    {
        ui->stellarMapWid->addSatellite(Gsates[i]);
    }
    for (int i=0;i<Csates.length();i++)   //将BDS卫星输入星空图
    {
        ui->stellarMapWid->addSatellite(Csates[i]);
    }
}
