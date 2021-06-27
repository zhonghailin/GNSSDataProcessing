#ifndef SATEPOSITIONTABLE_H
#define SATEPOSITIONTABLE_H

#include <QWidget>
#include <QTableWidget>
#include <QList>
#include <QMap>
#include <QDebug>
#include <QStandardItemModel>
#include <QListWidget>
#include <QTableWidgetItem>
#include <QHeaderView>
/********************************
 * 该类用于将卫星坐标做成表格并显示
 * *****************************/
class SatePositionTable : public QTableWidget
{
    Q_OBJECT
public:
    explicit SatePositionTable(QTableWidget *parent = nullptr);
    ~SatePositionTable();
    void initialization();//初始化表格
    /******************************************************
     * table的关键字为卫星的PRN号，其值为QList，其中包含的数据依次为：
     * X,Y,Z,L,B,H
     * ****************************************************/
    void clearALLData();   //清除所有行
public slots:
    void updateTable(QMap<QString,QList<double>> table); //更新表格
};

#endif // SATEPOSITIONTABLE_H
