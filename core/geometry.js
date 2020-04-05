"use strict";

class Geometry {
    /**
     * This class represents the geometry of a {@link module:Core.Mesh Mesh}. This includes information such
     * as the position of vertices as well as methods to compute edge lengths, corner
     * angles, face area, normals, discrete curvatures etc.
     * @constructor module:Core.Geometry
     * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
     * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
     * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is true.
     * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
     * @property {Object} positions A dictionary mapping each vertex to a normalized position.
     */
    constructor(mesh, positions, normalizePositions = true) {
        this.mesh = mesh;
        this.positions = {};
        for (let i = 0; i < positions.length; i++) {
            let v = this.mesh.vertices[i];
            let p = positions[i];

            this.positions[v] = p;
        }

        if (normalizePositions) {
            normalize(this.positions, mesh.vertices);
        }
    }

    /**
     * Computes the vector along a halfedge.
     * @method module:Core.Geometry#vector
     * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vector(h) {
        let a = this.positions[h.vertex];
        let b = this.positions[h.next.vertex];

        return b.minus(a);
    }

    /**
     * Computes the length of an edge.
     * @method module:Core.Geometry#length
     * @param {module:Core.Edge} e The edge whose length needs to be computed.
     * @returns {number}
     */
    length(e) {
        return this.vector(e.halfedge).norm();
    }

    /**
     * Computes the midpoint of an edge.
     * @method module:Core.Geometry#midpoint
     * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
     * @returns {number}
     */
    midpoint(e) {
        let h = e.halfedge;
        let a = this.positions[h.vertex];
        let b = this.positions[h.twin.vertex];

        return (a.plus(b)).over(2);
    }

    /**
     * Computes the mean edge length of all the edges in a mesh.
     * @method module:Core.Geometry#meanEdgeLength
     * @returns {number}
     */
    meanEdgeLength() {
        let sum = 0;
        let edges = this.mesh.edges;
        for (let e of edges) {
            sum += this.length(e);
        }

        return sum / edges.length;
    }

    /**
     * Computes the area of a face.
     * @method module:Core.Geometry#area
     * @param {module:Core.Face} f The face whose area needs to be computed.
     * @returns {number}
     */
    area(f) {
        if (f.isBoundaryLoop()) return 0.0;

        let u = this.vector(f.halfedge);
        let v = this.vector(f.halfedge.prev).negated();

        return 0.5 * u.cross(v).norm();
    }

    /**
     * Computes the total surface area of a mesh.
     * @method module:Core.Geometry#totalArea
     * @returns {number}
     */
    totalArea() {
        let sum = 0.0;
        for (let f of this.mesh.faces) {
            sum += this.area(f);
        }

        return sum;
    }

    /**
     * Computes the normal of a face.
     * @method module:Core.Geometry#faceNormal
     * @param {module:Core.Face} f The face whose normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    faceNormal(f) {
        if (f.isBoundaryLoop()) return undefined;

        let u = this.vector(f.halfedge);
        let v = this.vector(f.halfedge.prev).negated();

        return u.cross(v).unit();
    }

    /**
     * Computes the centroid of a face.
     * @method module:Core.Geometry#centroid
     * @param {module:Core.Face} f The face whose centroid needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    centroid(f) {
        let h = f.halfedge;
        let a = this.positions[h.vertex];
        let b = this.positions[h.next.vertex];
        let c = this.positions[h.prev.vertex];

        if (f.isBoundaryLoop()) return a.plus(b).over(2);

        return a.plus(b).plus(c).over(3);
    }

    /**
     * Computes the circumcenter of a face.
     * @method module:Core.Geometry#circumcenter
     * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    circumcenter(f) {
        let h = f.halfedge;
        let a = this.positions[h.vertex];
        let b = this.positions[h.next.vertex];
        let c = this.positions[h.prev.vertex];

        if (f.isBoundaryLoop()) return a.plus(b).over(2);

        let ac = c.minus(a);
        let ab = b.minus(a);
        let w = ab.cross(ac);

        let u = (w.cross(ab)).times(ac.norm2());
        let v = (ac.cross(w)).times(ab.norm2());
        let x = (u.plus(v)).over(2 * w.norm2());

        return x.plus(a);
    }

    /**
     * Computes an orthonormal bases for a face.
     * @method module:Core.Geometry#orthonormalBases
     * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
     * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
     */
    orthonormalBases(f) {
        let e1 = this.vector(f.halfedge).unit();

        let normal = this.faceNormal(f);
        let e2 = normal.cross(e1);

        return [e1, e2];
    }

    /**
     * Computes the angle (in radians) at a corner.
     * @method module:Core.Geometry#angle
     * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
     * @returns {number} The angle clamped between 0 and π.
     */
    angle(c) {

        let edge1 = this.vector(c.halfedge.next);
        let edge2 = this.vector(c.halfedge.prev);

        let dot = edge1.dot(edge2.negated());
        return Math.acos(dot / (edge1.norm() * edge2.norm()));
    }

    /**
     * Computes the cotangent of the angle opposite to a halfedge.
     * @method module:Core.Geometry#cotan
     * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
     * @returns {number}
     */
    cotan(h) {
        if (h.onBoundary) return 0; // if halfedge on boundary face, angle doesn't make sense, return 0 for Laplace

        let corner = h.corner;

        let edge1 = this.vector(corner.halfedge.next);
        let edge2 = this.vector(corner.halfedge.prev);

        let dot = edge1.dot(edge2.negated());
        let cross = edge1.cross(edge2).norm();

        return dot / cross;
    }

    /**
     * Computes the signed angle (in radians) between two adjacent faces.
     * @method module:Core.Geometry#dihedralAngle
     * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
     * the dihedral angle is computed.
     * @returns {number} The dihedral angle.
     */
    dihedralAngle(h) {
        let e = this.vector(h);  // e_ij

        let f1 = h.face;
        let f2 = h.twin.face;

        let N1 = this.faceNormal(f1);
        let N2 = this.faceNormal(f2);

        return Math.atan2(e.unit().dot(N1.cross(N2)), N1.dot(N2));
    }

    /**
     * Computes the barycentric dual area of a vertex.
     * @method module:Core.Geometry#barycentricDualArea
     * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
     * @returns {number}
     */
    barycentricDualArea(v) {

        let area = 0;
        for (let face of v.adjacentFaces()) {
            area += this.area(face) / 3;
        }

        return area;
    }

    /**
     * Computes the circumcentric dual area of a vertex.
     * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
     * @method module:Core.Geometry#circumcentricDualArea
     * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
     * @returns {number}
     */
    circumcentricDualArea(v) {

        let A = 0;
        for (let h of v.adjacentHalfedges()) {

            A += this.vector(h).norm2() * this.cotan(h);
            A += this.vector(h.prev).norm2() * this.cotan(h.prev);
        }

        return A / 8;
    }

    /**
     * Computes the normal at a vertex using the "equally weighted" method.
     * @method module:Core.Geometry#vertexNormalEquallyWeighted
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalEquallyWeighted(v) {
        let n = new Vector();
        for (let f of v.adjacentFaces()) {
            let normal = this.faceNormal(f);

            n.incrementBy(normal);
        }

        n.normalize();

        return n;
    }

    /**
     * Computes the normal at a vertex using the "face area weights" method.
     * @method module:Core.Geometry#vertexNormalAreaWeighted
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalAreaWeighted(v) {
        let N = new Vector(0, 0, 0);

        for (let f of v.adjacentFaces()) {
            let N_ijk = this.faceNormal(f);
            let A_ijk = this.area(f);
            N.incrementBy(N_ijk.times(A_ijk));
        }

        return N.unit();
    }

    /**
     * Computes the normal at a vertex using the "tip angle weights" method.
     * @method module:Core.Geometry#vertexNormalAngleWeighted
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalAngleWeighted(v) {
        let N = new Vector();

        for (let f of v.adjacentFaces()) {
            for (let c of f.adjacentCorners()) {
                if (c.vertex === v) {
                    N.incrementBy(this.faceNormal(f).times(this.angle(c)));
                    break;
                }
            }
        }

        return N.unit();
    }

    /**
     * Computes the normal at a vertex using the "gauss curvature" method.
     * @method module:Core.Geometry#vertexNormalGaussCurvature
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalGaussCurvature(v) {
        let N = new Vector();

        for (let h of v.adjacentHalfedges()) {
            let theta = this.dihedralAngle(h);
            let e = this.vector(h).unit();
            N.incrementBy(e.times(theta));
        }

        return N.over(2).unit();
    }

    /**
     * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
     * @method module:Core.Geometry#vertexNormalMeanCurvature
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalMeanCurvature(v) {
        let N = new Vector();

        for (let h of v.adjacentHalfedges()) {

            let cotA = this.cotan(h);
            let cotB = this.cotan(h.twin);

            N.incrementBy(this.vector(h).times(cotA + cotB));

        }
        return N.unit();
    }

    /**
     * Computes the normal at a vertex using the "inscribed sphere" method.
     * @method module:Core.Geometry#vertexNormalSphereInscribed
     * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
     * @returns {module:LinearAlgebra.Vector}
     */
    vertexNormalSphereInscribed(v) {
        let N = new Vector();

        for (let h of v.adjacentHalfedges()) {

            let e_ij = this.vector(h);
            let e_ik = this.vector(h.prev.twin);

            N.incrementBy(e_ij.cross(e_ik).over(e_ij.norm2() * e_ik.norm2()))
        }

        return N.unit();
    }

    /**
     * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
     * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
     * @method module:Core.Geometry#angleDefect
     * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
     * @returns {number}
     */
    angleDefect(v) {

        let incidentAngles = 0;
        let isBoundary = v.onBoundary();

        for (let h of v.adjacentHalfedges()) {
            let c = h.next.corner;
            incidentAngles += this.angle(c);
        }

        if (isBoundary) return Math.PI - incidentAngles;
        else return 2 * Math.PI - incidentAngles;
    }

    /**
     * Computes the (integrated) scalar gauss curvature at a vertex.
     * @method module:Core.Geometry#scalarGaussCurvature
     * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
     * @returns {number}
     */
    scalarGaussCurvature(v) {
        return this.angleDefect(v);
    }

    /**
     * Computes the (integrated) scalar mean curvature at a vertex.
     * @method module:Core.Geometry#scalarMeanCurvature
     * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
     * @returns {number}
     */
    scalarMeanCurvature(v) {
        let sum = 0;

        for (let h of v.adjacentHalfedges()) {
            sum += this.length(h.edge) * this.dihedralAngle(h);
        }

        return sum / 2;
    }

    /**
     * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
     * @method module:Core.Geometry#totalAngleDefect
     * @returns {number}
     */
    totalAngleDefect() {
        let V = this.mesh.vertices.length;
        let E = this.mesh.edges.length;
        let F = this.mesh.faces.length;
        return 2 * Math.PI * (V - E + F)
    }

    /**
     * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
     * @method module:Core.Geometry#principalCurvatures
     * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
     * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
     */
    principalCurvatures(v) {

        let H = this.scalarMeanCurvature(v) / this.circumcentricDualArea(v);
        let K = this.scalarGaussCurvature(v) / this.circumcentricDualArea(v);

        return [H - Math.sqrt(H * H - K), H + Math.sqrt(H * H - K)];
    }

    /**
     * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
     * instead we build a positive definite matrix by multiplying the entries of the
     * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
     * @method module:Core.Geometry#laplaceMatrix
     * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    laplaceMatrix(vertexIndex) {
        // TODO

        return SparseMatrix.identity(1, 1); // placeholder
    }

    /**
     * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
     * of a mesh.
     * @method module:Core.Geometry#massMatrix
     * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    massMatrix(vertexIndex) {
        // TODO

        return SparseMatrix.identity(1, 1); // placeholder
    }

    /**
     * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
     * instead we build a positive definite matrix by multiplying the entries of the
     * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
     * @method module:Core.Geometry#complexLaplaceMatrix
     * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
     * @returns {module:LinearAlgebra.ComplexSparseMatrix}
     */
    complexLaplaceMatrix(vertexIndex) {
        // TODO

        return ComplexSparseMatrix.identity(1, 1); // placeholder
    }
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
    // compute center of mass
    let N = vertices.length;
    let cm = new Vector();
    for (let v of vertices) {
        let p = positions[v];

        cm.incrementBy(p);
    }
    cm.divideBy(N);

    // translate to origin and determine radius
    let radius = -1;
    for (let v of vertices) {
        let p = positions[v];

        p.decrementBy(cm);
        radius = Math.max(radius, p.norm());
    }

    // rescale to unit radius
    if (rescale) {
        for (let v of vertices) {
            let p = positions[v];

            p.divideBy(radius);
        }
    }
}