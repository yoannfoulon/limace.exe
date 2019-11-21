#ifndef MAINWI#NDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // les fonctions à compléter
    float faceArea(MyMesh* _mesh, int faceID);
    float baryArea(MyMesh* _mesh, int vertID);
    float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    void H_Curv(MyMesh* _mesh);
    void K_Curv(MyMesh* _mesh);
    OpenMesh::Vec3f LaplaceBeltrami(MyMesh* _mesh, int vertID);
    void LBAll(MyMesh* _mesh, double h, double lambda);
    void matriceLB(MyMesh* _mesh);
    float cot(float angle);

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_LB_clicked();

    void on_h_valueChanged(double arg1);

    void on_lambda_valueChanged(double arg1);

    void on_pushButton_mat_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    double h = 0.25;
    double lambda = 0.50;
    int vertexSelection;
    int edgeSelection;
    int faceSelection;


    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
