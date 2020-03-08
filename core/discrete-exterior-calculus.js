"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
    /**
     * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
     * By convention, the area of a vertex is 1.
     * @static
     * @param {module:Core.Geometry} geometry The geometry of a mesh.
     * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    static buildHodgeStar0Form(geometry, vertexIndex) {
        let T = new Triplet(geometry.mesh.vertices.length, geometry.mesh.vertices.length);

        for (let vertex of geometry.mesh.vertices) {
            let v = vertexIndex[vertex];

            let dual = geometry.barycentricDualArea(vertex);
            T.addEntry(dual, v, v);
        }

        return SparseMatrix.fromTriplet(T);
    }

    /**
     * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
     * @static
     * @param {module:Core.Geometry} geometry The geometry of a mesh.
     * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    static buildHodgeStar1Form(geometry, edgeIndex) {
        let T = new Triplet(geometry.mesh.edges.length, geometry.mesh.edges.length);

        for (let edge of geometry.mesh.edges) {
            let e = edgeIndex[edge];

            let dual = (geometry.cotan(edge.halfedge) + geometry.cotan(edge.halfedge.twin)) / 2;
            T.addEntry(dual, e, e);
        }

        return SparseMatrix.fromTriplet(T);
    }

    /**
     * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
     * By convention, the area of a vertex is 1.
     * @static
     * @param {module:Core.Geometry} geometry The geometry of a mesh.
     * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    static buildHodgeStar2Form(geometry, faceIndex) {
        let T = new Triplet(geometry.mesh.faces.length, geometry.mesh.faces.length);

        for (let face of geometry.mesh.faces) {
            let f = faceIndex[face];

            let edges = [...face.adjacentEdges()];


            let s = (geometry.length(edges[0]) + geometry.length(edges[1]) + geometry.length(edges[2])) / 2;

            let dual = 1 / Math.sqrt(s * (s - geometry.length(edges[0])) * (s - geometry.length(edges[1])) * (s - geometry.length(edges[2])));
            T.addEntry(dual, f, f);
        }

        return SparseMatrix.fromTriplet(T);
    }

    /**
     * Builds a sparse matrix encoding the exterior derivative on 0-forms.
     * @static
     * @param {module:Core.Geometry} geometry The geometry of a mesh.
     * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
     * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
        let mesh = geometry.mesh;
        let T = new Triplet(mesh.edges.length, mesh.vertices.length);

        for (let i = 0; i < mesh.edges.length; i++) {
            let e = edgeIndex[i];

            let v1 = vertexIndex[mesh.edges[e].halfedge.vertex.index];
            let v2 = vertexIndex[mesh.edges[e].halfedge.twin.vertex.index];
            //vert -> edge (oriented)
            T.addEntry(-1, e, v1); // orientation The vertex at the base of this halfedge.
            T.addEntry(1, e, v2);
        }
        return SparseMatrix.fromTriplet(T);
    }


    /**
     * Builds a sparse matrix encoding the exterior derivative on 1-forms.
     * @static
     * @param {module:Core.Geometry} geometry The geometry of a mesh.
     * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
     * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
     * @returns {module:LinearAlgebra.SparseMatrix}
     */
    static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
        //I copied this from a fork, I don't understand this either

        const edges = geometry.mesh.edges;
        const faces = geometry.mesh.faces;
        let T = new Triplet(faces.length, edges.length);
        for (let f of faces) {
            const row = faceIndex[f];
            for (let e of f.adjacentEdges()) {
                // get edge orientation here. Not sure if I understand this fully.
                // see http://brickisland.net/DDGSpring2019/2019/02/13/1669/#comment-107
                const entry = e.halfedge.face === f ? 1 : -1;
                T.addEntry(entry, row, edgeIndex[e]);
            }
        }
        return SparseMatrix.fromTriplet(T);
    }
}
