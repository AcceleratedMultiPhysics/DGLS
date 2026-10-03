#!/usr/bin/env python3
"""Structured quadrilateral mesh of [x0,x1]x[y0,y1] in Gmsh 2.2 ASCII format, with the
boundary condition id of each side as the physical tag of its boundary lines
(ins/mps convention: 1 wall, 2 inflow, 3 outflow, 4 x-slip (u=0), 5 y-slip (v=0)).

usage: genRectQuad.py out.msh x0 x1 y0 y1 NX NY bcBottom bcRight bcTop bcLeft
"""
import sys
out = sys.argv[1]
x0, x1, y0, y1 = map(float, sys.argv[2:6])
NX, NY = map(int, sys.argv[6:8])
bcB, bcR, bcT, bcL = map(int, sys.argv[8:12])

def nid(i, j): return j*(NX+1) + i + 1
lines = []
for i in range(NX):   lines.append((bcB, nid(i,0),   nid(i+1,0)))
for j in range(NY):   lines.append((bcR, nid(NX,j),  nid(NX,j+1)))
for i in range(NX):   lines.append((bcT, nid(i+1,NY),nid(i,NY)))
for j in range(NY):   lines.append((bcL, nid(0,j+1), nid(0,j)))

with open(out, 'w') as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n" % ((NX+1)*(NY+1)))
    for j in range(NY+1):
        for i in range(NX+1):
            f.write("%d %.16g %.16g 0\n" % (nid(i,j), x0+(x1-x0)*i/NX, y0+(y1-y0)*j/NY))
    f.write("$EndNodes\n$Elements\n%d\n" % (len(lines) + NX*NY))
    k = 1
    for (bc, a, b) in lines:
        f.write("%d 1 2 %d 1 %d %d\n" % (k, bc, a, b)); k += 1
    for j in range(NY):
        for i in range(NX):
            f.write("%d 3 2 9 1 %d %d %d %d\n" % (k, nid(i,j), nid(i+1,j), nid(i+1,j+1), nid(i,j+1))); k += 1
    f.write("$EndElements\n")
print("wrote", out, NX*NY, "quads")
