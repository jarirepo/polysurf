# polysurf
Bilinear interpolation of four connected polylines

Performs a bilinear interpolation between four connected 3-D parametric functions in the XYZ-space. The parametric functions are first discretized into polylines (piecewise linear functions) and then interpolated using the Coon's bilinear interpolation method which is implemented for polyline objects. This implementation allows each boundary function to have an arbitrary number of sampling points. The bilinear interpolation method itself allows advanced surfaces to be generated from arbitrary analytical functions and point sequences.
