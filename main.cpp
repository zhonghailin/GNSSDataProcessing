#include "mainwindow.h"
#include <QApplication>
#include <QDebug>
#include <iostream>
#include <QList>
#include <QFrame>
#include <QDateTime>
#include <stellarmap.h>
#include <QMap>
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.showMaximized();//最大化显示


    return a.exec();
}
