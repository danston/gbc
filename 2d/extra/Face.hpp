// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is the face data structure.

    This class does not have any external dependencies.

*/

#ifndef GBC_FACE_HPP
#define GBC_FACE_HPP

namespace gbc {

    // Face data structure.
    class Face {

    public:
        // Constructors.
        Face() { v[0] = v[1] = v[2] = v[3] = f[0] = f[1] = f[2] = f[3] = -1; }

        Face(const int v0, const int v1, const int v2) {
            v[0] = v0;
            v[1] = v1;
            v[2] = v2;
            v[3] = f[0] = f[1] = f[2] = f[3] = -1;
        }

        Face(const int v0, const int v1, const int v2, const int v3) {
            v[0] = v0;
            v[1] = v1;
            v[2] = v2;
            v[3] = v3;
            f[0] = f[1] = f[2] = f[3] = -1;
        }

        // Elements.
        int v[4]; // indices of the vertices
        int f[4]; // indices of the neighboring faces
    };

} // namespace gbc

#endif // GBC_FACE_HPP
