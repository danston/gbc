// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is the halfedge data structure.

    This class does not have any external dependencies.

*/

#ifndef GBC_HALFEDGE_HPP
#define GBC_HALFEDGE_HPP

namespace gbc {

    // Halfedge data structure.
    class Halfedge {

    public:
        // Constructor.
        Halfedge() : prev(-1), next(-1), neigh(-1), dest(-1) { }

        // Flags.
        int prev;  // index of the previous halfedge
        int next;  // index of the next halfedge
        int neigh; // neighboring halfedge
        int dest;  // index of the vertex at the end of the halfedge
    };

} // namespace gbc

#endif // GBC_HALFEDGE_HPP
