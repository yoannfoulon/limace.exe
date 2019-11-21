short sign = 0;

for (MyMesh::VertexFaceCWIter curFace = _mesh->vf_cwiter(_mesh->vertex_handle(vertID0)) ; curFace.is_valid() && sign == 0; curFace ++)
{
    MyMesh::FaceHandle fh = *curFace;
    if (fh.idx() == faceID0)
    {
        curFace++;
        fh = *curFace;
        if (fh.idx() == faceID1)
            sign = 1;
        else
        {
            curFace--;
            fh = *curFace;
        }
    }
    if (fh.idx() == faceID1)
    {
        curFace++;
        fh = *curFace;
        if (fh.idx() == faceID0)
            sign = -1;
    }
}

VectorT<float,3> n1(_mesh->normal(_mesh->face_handle(faceID0)));
VectorT<float,3> n2(_mesh->normal(_mesh->face_handle(faceID1)));

float dp = n1|n2;

float normal = acos(dp);

return normal*sign;
