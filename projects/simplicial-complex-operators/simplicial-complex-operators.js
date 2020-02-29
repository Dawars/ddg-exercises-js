"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

    /** This class implements various operators (e.g. boundary, star, link) on a mesh.
     * @constructor module:Projects.SimplicialComplexOperators
     * @param {module:Core.Mesh} mesh The input mesh this class acts on.
     * @property {module:Core.Mesh} mesh The input mesh this class acts on.
     * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
     * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
     */
    constructor(mesh) {
        this.mesh = mesh;
        this.assignElementIndices(this.mesh);

        this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
        this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
    }

    /** Assigns indices to the input mesh's vertices, edges, and faces
     * @method module:Projects.SimplicialComplexOperators#assignElementIndices
     * @param {module:Core.Mesh} mesh The input mesh which we index.
     */
    assignElementIndices(mesh) {

        let set_index = function (item, index, arr) {
            item.index = index;
        };

        mesh.vertices.forEach(set_index);
        mesh.edges.forEach(set_index);
        mesh.faces.forEach(set_index);

    }

    /** Returns the vertex-edge adjacency matrix of the given mesh.
     * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
     * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
     * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
     */
    buildVertexEdgeAdjacencyMatrix(mesh) {
        let T = new Triplet(mesh.edges.length, mesh.vertices.length);

        for (let i = 0; i < mesh.edges.length; i++) {
            T.addEntry(1, mesh.edges[i].index, mesh.edges[i].halfedge.vertex.index);
            T.addEntry(1, mesh.edges[i].index, mesh.edges[i].halfedge.twin.vertex.index);
        }
        return SparseMatrix.fromTriplet(T);
    }

    /** Returns the edge-face adjacency matrix.
     * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
     * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
     * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
     */
    buildEdgeFaceAdjacencyMatrix(mesh) {
        let T = new Triplet(mesh.faces.length, mesh.edges.length);
        for (let i = 0; i < mesh.faces.length; i++) {

            let it = mesh.faces[i].adjacentEdges();
            let edges = [...it];

            for (let j = 0; j < edges.length; j++) {
                T.addEntry(1, mesh.faces[i].index, edges[j].index);
            }
        }
        return SparseMatrix.fromTriplet(T);
    }

    /** Returns a column vector representing the vertices of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildVertexVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
     *  vertex i is in the given subset and 0 otherwise
     */
    buildVertexVector(subset) {
        let vector = DenseMatrix.zeros(this.mesh.vertices.length);
        for (let vertexIndex of subset.vertices) {
            vector.set(1, vertexIndex);
        }

        return vector;
    }

    /** Returns a column vector representing the edges of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
     *  edge i is in the given subset and 0 otherwise
     */
    buildEdgeVector(subset) {
        let vector = DenseMatrix.zeros(this.mesh.edges.length);
        for (let edgeIndex of subset.edges) {
            vector.set(1, edgeIndex);
        }

        return vector;
    }

    /** Returns a column vector representing the faces of the
     * given subset.
     * @method module:Projects.SimplicialComplexOperators#buildFaceVector
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
     *  face i is in the given subset and 0 otherwise
     */
    buildFaceVector(subset) {
        let vector = DenseMatrix.zeros(this.mesh.faces.length);
        for (let faceIndex of subset.faces) {
            vector.set(1, faceIndex);
        }

        return vector;
    }

    /** Returns the star of a subset.
     * @method module:Projects.SimplicialComplexOperators#star
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The star of the given subset.
     */
    star(subset) {
        let star = MeshSubset.deepCopy(subset);

        // vector containing original vertices of subset
        let vertexVector = this.buildVertexVector(subset);
        let edgeVector = this.buildEdgeVector(subset);

        // vector containing edges containing vertexVector
        let newEdges = this.A0.timesDense(vertexVector);

        // faces containing original edges
        let faces1 = this.A1.timesDense(edgeVector);
        // faces containing new edges
        let faces2 = this.A1.timesDense(newEdges);


        for (let i = 0; i < this.A0.nRows(); i++) {
            if (newEdges.get(i, 0) > 0) {
                star.addEdge(i);
            }
        }
        for (let i = 0; i < this.A1.nRows(); i++) {
            if (faces1.get(i, 0) > 0 || faces2.get(i, 0) > 0) {
                star.addFace(i);
            }
        }

        return star; // placeholder
    }

    /**
     * Print matrix elements helper function
     * @method module:Projects.SimplicialComplexOperators#printMatrix
     * @param {module:LinearAlgebra.DenseMatrix} mat A subset of our mesh.
     */
    printMatrix(mat) {

        let text = "";

        text += "[" + mat.nRows() + ", " + mat.nCols() + "]\n";
        for (let i = 0; i < mat.nRows(); i++) {
            let row = "";
            for (let j = 0; j < mat.nCols(); j++) {
                row += mat.get(i, j) + " ";
            }
            text += row + "\n";
        }
        console.log(text);
    }

    /** Returns the closure of a subset.
     * @method module:Projects.SimplicialComplexOperators#closure
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The closure of the given subset.
     */
    closure(subset) {
        let closure = MeshSubset.deepCopy(subset);

        // faces -> edges
        let faces = this.buildFaceVector(closure);
        let connectedEdges = this.A1.transpose().timesDense(faces);
        //add edges
        for (let i = 0; i < connectedEdges.nRows(); i++) {
            if (connectedEdges.get(i, 0) !== 0) {
                closure.addEdge(i);
            }
        }

        //edges -> verts
        let edges = this.buildEdgeVector(closure);

        let connectedVerts = this.A0.transpose().timesDense(edges);
        for (let i = 0; i < connectedVerts.nRows(); i++) {
            if (connectedVerts.get(i, 0) !== 0) {
                closure.addVertex(i);
            }
        }

        return closure; // placeholder
    }

    /** Returns the link of a subset.
     * @method module:Projects.SimplicialComplexOperators#link
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {module:Core.MeshSubset} The link of the given subset.
     */
    link(subset) {
        let closure_star = this.closure(this.star(subset));
        let star_closure = this.star(this.closure(subset));
        closure_star.deleteSubset(star_closure);

        return closure_star;
    }

    /** Returns true if the given subset is a subcomplex and false otherwise.
     * @method module:Projects.SimplicialComplexOperators#isComplex
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
     */
    isComplex(subset) {
        let closure = this.closure(subset);

        return closure.equals(subset);
    }

    /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
     * @method module:Projects.SimplicialComplexOperators#isPureComplex
     * @param {module:Core.MeshSubset} subset A subset of our mesh.
     * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
     */
    isPureComplex(subset) {
        if (!this.isComplex(subset))
            return -1;  // if not simplical complex, cannot be pure

        //empty subset
        if (subset.equals(new MeshSubset()))
            return -1;

        // if no edges and faces
        if (subset.edges.size === 0 && subset.faces.size === 0)
            return 0;

        // check number of edges containing the vertices in this set
        let verts = this.buildVertexVector(subset);

        for (let i = 0; i < verts.nRows(); i++) {
            // if vertex is in subset
            if (verts.get(i) !== 0) {
                // get column of A0
                let col = this.A0.subMatrix(0, this.A0.nRows(), i, i + 1).toDense();  // todo maybe not a good idea

                let connectedEdges = 0;  // vertex has to be connected to an edge in this subset, not the mesh

                for (let edgeIndex of subset.edges) {
                    connectedEdges += col.get(edgeIndex);  // I assume elements are 1 or 0
                }

                //    this must be > 0 for every vertex *in the subset*
                if (connectedEdges === 0)
                    return -1;
            }
        }

        // if no faces
        if (subset.faces.size === 0)
            return 1;

        // check number of faces containing the edges in this set
        let edges = this.buildEdgeVector(subset);

        for (let i = 0; i < edges.nRows(); i++) {
            // if vertex is in subset
            if (edges.get(i) !== 0) {
                // get column of A1
                let col = this.A1.subMatrix(0, this.A1.nRows(), i, i + 1).toDense();  // todo reconsider

                let connectedFaces = 0;  // edge has to be connected to a face in this subset, not the mesh

                for (let faceIndex of subset.faces) {
                    connectedFaces += col.get(faceIndex);  // I assume elements are 1 or 0
                }

                //  this must be > 0 for every edge *in the subset*
                if (connectedFaces === 0)
                    return -1;
            }
        }
        return 2;
    }

    /** Returns the boundary of a subset.
     * @method module:Projects.SimplicialComplexOperators#boundary
     * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
     * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
     */
    boundary(subset) {
        let bounds = new MeshSubset();

        let k_complex = this.isPureComplex(subset);
        if (k_complex === -1)
            return bounds; // empty

        if (k_complex === 2) {
            // faces, edges, verts

            //vector to count occurances of edges (verts)
            let numEdges = DenseMatrix.zeros(this.mesh.edges.length);
            // loop faces
            for (let face of subset.faces) {
                //get edges connected to face

                let set = new MeshSubset([], [], [face]);
                let subsets = this.closure(set);  // creates *faces* of set (incl itself)

                // loop edges
                for (let e of subsets.edges) {
                    numEdges.set(1 + numEdges.get(e), e);
                }
            }

            for (let i = 0; i < this.mesh.edges.length; i++) {
                if (numEdges.get(i) === 1) {
                    bounds.addEdge(i);
                }
            }
        }

        if (k_complex === 1) {
            // only lines and verts

            //vector to count occurances of edges (verts)
            let numVertices = DenseMatrix.zeros(this.mesh.vertices.length);
            // loop faces
            for (let edge of subset.edges) {
                //get edges connected to face

                let set = new MeshSubset([], [edge], []);
                let subsets = this.closure(set);  // creates *faces* of set (incl itself)

                // loop verts
                for (let v of subsets.vertices) {
                    numVertices.set(1 + numVertices.get(v), v);
                }
            }

            for (let i = 0; i < this.mesh.vertices.length; i++) {
                if (numVertices.get(i) === 1) {
                    bounds.addVertex(i);
                }
            }
        }

        return this.closure(bounds); // placeholder
    }
}
