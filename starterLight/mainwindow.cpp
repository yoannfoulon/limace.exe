#include "mainwindow.h"
#include "ui_mainwindow.h"

// Calcul de l'aire d'une face du maillage
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );

    std::vector<VertexHandle> vertexes;

    // On récupère les 3 sommets de la face
    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v)
        vertexes.push_back ( *fh_v );

    // On créé deux vecteurs avec le même point de départ
    OpenMesh::Vec3f vectorAB = _mesh->point(vertexes[1]) - _mesh->point(vertexes[0]);
    OpenMesh::Vec3f vectorAC = _mesh->point(vertexes[2]) - _mesh->point(vertexes[0]);

    // On calcule le produit vectoriel de ses deux vecteurs et retourne sa norme divisée par 2
    OpenMesh::Vec3f product = vectorAB % vectorAC;
    float norm = product.norm();

    return norm / 2.0f;
}

// Calcule de l'aire barycentrique sous un sommet
float MainWindow::baryArea(MyMesh* _mesh, int vertID){
    float baryArea = 0;

    VertexHandle vh = _mesh->vertex_handle ( vertID );

    // On somme toutes les aires des faces voisines au sommet
    MyMesh::VertexFaceIter vf = _mesh->vf_iter ( vh );
    for ( ; vf.is_valid ( ) ; ++vf ) {
        FaceHandle current = *vf;
        baryArea += faceArea ( _mesh , current.idx( ) );
    }

    // On retourne cette somme divisée par 3
    return baryArea / 3.0f;
}

// Calcul de l'angle entre deux faces
float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    int sign = 0;

    VertexHandle vh0 = _mesh->vertex_handle ( vertID0 );
    FaceHandle fh0 = _mesh->face_handle ( faceID0 );


    // On créé un itérateur clockwise sur les sommets de la face d'handle faceID0
    // et on tourne jusqu'à arriver sur le sommet d'handle vertID0
    MyMesh::FaceVertexCWIter fh_cwv = _mesh->fv_cwiter ( fh0 );
    while ( fh_cwv.is_valid ( ) && *fh_cwv != vh0 ) ++fh_cwv;

    VertexHandle next = *++fh_cwv;

    // Si le suivant du sommet d'handle vertID0 sur la face d'id faceID0 est le sommet
    // d'handle vertID1, on est dans le sens négatif des faces, le signe de l'angle
    // est négatif. Sinon, il est positif.
    if ( next.idx ( ) == vertID1 ) sign = -1;
    else sign = 1;

    // On récupère les normales, calcule leur produit scalaire et on renvoie l'angle
    // signé avec le produit scalaire et le signe trouvé plus tot
    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( _mesh->face_handle ( faceID1 ) ) );
    float scalar = normal0 | normal1;

    return sign * acos ( scalar );
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );
    VertexHandle vh = _mesh->vertex_handle ( vertexID );
    std::vector<VertexHandle> vertexes;

    // On récupère les sommets de la faces qui ne sont pas le sommet d'ID vertexID
    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v) {
        VertexHandle current = *fh_v;
        if( current.idx() != vertexID )
            vertexes.push_back ( current );
    }

    // On retourne l'angle calculé depuis le produit scalaire
    OpenMesh::Vec3f vectorAB = _mesh->point(vertexes[0]) - _mesh->point(vh);
    OpenMesh::Vec3f vectorAC = _mesh->point(vertexes[1]) - _mesh->point(vh);
    return acos ( vectorAB.normalize() | vectorAC.normalize() );
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        MyMesh::VertexHandle current = *curVert;
        float val = 0.0f;

        for (MyMesh::VertexEdgeIter currentEdge = _mesh->ve_iter ( current ); currentEdge.is_valid(); currentEdge++)
        {
            MyMesh::EdgeHandle eh = *currentEdge;
            MyMesh::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);
            MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1);

            FaceHandle fh0 = _mesh->face_handle(heh0);
            FaceHandle fh1 = _mesh->face_handle(heh1);

            // Si l'arête est en bordure, on ne traite qu'une face
            if ( fh1.idx ( ) > _mesh->n_faces ( ) )
                fh1 = fh0;

            // On détermine le sommet opposé au sommet courant sur l'arête courante
            int vertex2ID = _mesh->to_vertex_handle(heh1).idx();
            if (vertex2ID == current.idx ( ) )
                vertex2ID = _mesh->to_vertex_handle(heh0).idx();

            // Vecteur entre le sommet courant et son opposé sur l'arête courante
            OpenMesh::Vec3f currentOppVector = _mesh->point ( _mesh->vertex_handle ( vertex2ID ) ) - _mesh->point ( current );

            OpenMesh::Vec3f normal0 ( _mesh->normal ( fh0 ) );
            OpenMesh::Vec3f normal1 ( _mesh->normal ( fh1 ) );

            // On ordonne les faces
            if ( ( ( normal0 % normal1 ) | currentOppVector ) < 0 )
            {
                FaceHandle tempF = fh0;
                fh0 = fh1;
                fh1 = tempF;
            }

            // On somme les valeurs des angles entre les faces * le vecteur entre le sommet et son opposé
            val += currentOppVector.norm ( ) * angleFF ( _mesh , fh0.idx ( ) , fh1.idx ( ) , current.idx() , vertex2ID );
        }

        // On calcule le résultat sur le sommet et on le lui attribue dans les propriétés du maillage
        val /= ( 4 * baryArea ( _mesh , current.idx ( ) ) );
        _mesh->data ( current ).value = val;
    }
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;

        // On calcule l'angle barycentrique sous le sommet courant
        float area = baryArea ( _mesh , current.idx() );
        float angleEESum = 0;

        //On itère les faces autour du sommet courant et somme les angles entre leurs arêtes
        MyMesh::VertexFaceIter vf = _mesh->vf_iter ( current );
        for ( ; vf.is_valid ( ) ; ++vf ) {
            FaceHandle currentFace = *vf;
            angleEESum += angleEE ( _mesh , current.idx ( ) , currentFace.idx ( ) );
        }

        // On calcule le résultat sur le sommet et on le lui attribue dans les propriétés du maillage
        _mesh->data ( current ).value = ( 1 / area ) * ( 2 * M_PI - angleEESum );
    }
}

float MainWindow::cot(float angle){
    return cos(angle) / sin(angle);
}

void MainWindow::matriceLB(MyMesh *_mesh){
    //std::map<std::pair<int, int>, float> matriceLaplaceBeltrami;
    //size_t nb_vertices = _mesh->n_vertices;
    std::vector<std::vector<float>> matriceLaplaceBeltrami(_mesh->n_vertices(), std::vector<float>(_mesh->n_vertices(), 0));

    for ( MyMesh::EdgeIter curEdge = _mesh->edges_begin() ; curEdge!=_mesh->edges_end() ; ++curEdge ) {
        HalfedgeHandle heh0 = _mesh->halfedge_handle(curEdge, 0);
        HalfedgeHandle heh1 = _mesh->halfedge_handle(curEdge, 1);

        FaceHandle fh0 = _mesh->face_handle(heh0);
        FaceHandle fh1 = _mesh->face_handle(heh1);

        if(!_mesh->is_valid_handle(fh0) || !_mesh->is_valid_handle(fh1)) continue;

        VertexHandle vh0 = _mesh->to_vertex_handle(heh0);
        VertexHandle vh1 = _mesh->to_vertex_handle(heh1);

        OpenMesh::Vec3f coordsVh0 = _mesh->point(vh0);
        OpenMesh::Vec3f coordsVh1 = _mesh->point(vh1);

        float angle0;
        float angle1;
        VertexHandle thirdPoint0;
        OpenMesh::Vec3f vec0_0;
        OpenMesh::Vec3f vec0_1;
        VertexHandle thirdPoint1;
        OpenMesh::Vec3f vec1_0;
        OpenMesh::Vec3f vec1_1;

        MyMesh::FaceVertexIter fv = _mesh->fv_iter(fh0);
        for ( ; fv.is_valid ( ) ; ++fv ) {
            VertexHandle currentVert = *fv;
            if(currentVert != vh0  && currentVert != vh1){
                thirdPoint0 = currentVert;
            }
        }

        vec0_0 = coordsVh0 - _mesh->point ( thirdPoint0 );
        vec0_1 = coordsVh1 - _mesh->point ( thirdPoint0 );

        angle0 = acos(vec0_0 | vec0_1);

        fv = _mesh->fv_iter(fh1);
        for ( ; fv.is_valid ( ) ; ++fv ) {
            VertexHandle currentVert = *fv;
            if(currentVert != vh0 && currentVert != vh1){
                thirdPoint1 = currentVert;
            }
        }

        vec1_0 = coordsVh0 - _mesh->point ( thirdPoint1 );
        vec1_1 = coordsVh1 - _mesh->point ( thirdPoint1 );

        angle1 = acos(vec1_0 | vec1_1);

        float coef = cot(angle0) + cot(angle1);
        matriceLaplaceBeltrami[vh0.idx()][vh1.idx()] = coef;
        matriceLaplaceBeltrami[vh0.idx()][vh0.idx()] += coef;
        matriceLaplaceBeltrami[vh1.idx()][vh1.idx()] += coef;
    }

    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;
        matriceLaplaceBeltrami[current.idx()][current.idx()] = -matriceLaplaceBeltrami[current.idx()][current.idx()];
    }

    for ( int i = 0 ; i < matriceLaplaceBeltrami.size() ; ++i ){
        for ( int j = 0 ; j < matriceLaplaceBeltrami[i].size() ; ++j ){
            std::cout << std::setprecision(5) << "Matrice[" << i << "][" << j << "]=" << matriceLaplaceBeltrami[i][j] << std::endl;
        }
    }
}

OpenMesh::Vec3f MainWindow::LaplaceBeltrami(MyMesh* _mesh, int vertID){
    float area = baryArea(_mesh, vertID);
    VertexHandle current = _mesh->vertex_handle(vertID);
    MyMesh::VertexEdgeIter ve = _mesh->ve_iter ( current );
    OpenMesh::Vec3f coordsV = _mesh->point(current);
    OpenMesh::Vec3f sum;
    for ( ; ve.is_valid ( ) ; ++ve ) {
        EdgeHandle currentEdge = *ve;
        HalfedgeHandle heh0 = _mesh->halfedge_handle(currentEdge, 0);
        HalfedgeHandle heh1 = _mesh->halfedge_handle(currentEdge, 1);

        FaceHandle fh0 = _mesh->face_handle(heh0);
        FaceHandle fh1 = _mesh->face_handle(heh1);

        if(!_mesh->is_valid_handle(fh0) || !_mesh->is_valid_handle(fh1)) continue;

        VertexHandle opp = _mesh->to_vertex_handle(heh0);
        if(opp == current) opp = _mesh->from_vertex_handle(heh0);

        OpenMesh::Vec3f coordsVOpp = _mesh->point(opp);

        float angle0;
        float angle1;
        VertexHandle thirdPoint0;
        OpenMesh::Vec3f vec0_0;
        OpenMesh::Vec3f vec0_1;
        VertexHandle thirdPoint1;
        OpenMesh::Vec3f vec1_0;
        OpenMesh::Vec3f vec1_1;

        MyMesh::FaceVertexIter fv = _mesh->fv_iter(fh0);
        for ( ; fv.is_valid ( ) ; ++fv ) {
            VertexHandle currentVert = *fv;
            if(currentVert != current  && currentVert != opp){
                thirdPoint0 = currentVert;
            }
        }

        vec0_0 = coordsV - _mesh->point ( thirdPoint0 );
        vec0_1 = coordsVOpp - _mesh->point ( thirdPoint0 );

        angle0 = acos(vec0_0 | vec0_1);

        fv = _mesh->fv_iter(fh1);
        for ( ; fv.is_valid ( ) ; ++fv ) {
            VertexHandle currentVert = *fv;
            if(currentVert != current && currentVert != opp){
                thirdPoint1 = currentVert;
            }
        }

        vec1_0 = coordsV - _mesh->point ( thirdPoint1 );
        vec1_1 = coordsVOpp - _mesh->point ( thirdPoint1 );

        angle1 = acos(vec1_0 | vec1_1);

        OpenMesh::Vec3f vToOpp = coordsVOpp - coordsV;
        float coef = cot(angle0) + cot(angle1);
        sum += coef * vToOpp;
    }

    return sum / (2.0 * area);
}

void MainWindow::LBAll(MyMesh* _mesh, double h, double lambda){
    resetAllColorsAndThickness(_mesh);
    std::map<int, OpenMesh::Vec3f> newCoords;
    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;
        OpenMesh::Vec3f newPosition = LaplaceBeltrami(_mesh, current.idx());
        newCoords[current.idx()] = _mesh->point(current) + h * lambda * newPosition;
    }

    for (std::map<int, OpenMesh::Vec3f>::iterator it= newCoords.begin(); it!=newCoords.end(); ++it){
        MyMesh::Point newPos = it->second;
        _mesh->set_point ( _mesh->vertex_handle(it->first) , newPos );
        qDebug() << it->second[0] << it->second[1] << it->second[2];
    }


    displayMesh(_mesh);
}


void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_LB_clicked()
{
    LBAll(&mesh, h, lambda);
}
void MainWindow::on_h_valueChanged(double arg1)
{
    h = arg1;
}

void MainWindow::on_lambda_valueChanged(double arg1)
{
    lambda = arg1;
}

void MainWindow::on_pushButton_mat_clicked()
{
    matriceLB(&mesh);
}
