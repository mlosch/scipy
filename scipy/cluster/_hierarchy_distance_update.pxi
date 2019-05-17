"""
A `linkage_distance_update` function calculates the distance from cluster i
to the new cluster xy after merging cluster x and cluster y

Parameters
----------
d_xi : float
    Distance from cluster x to cluster i
d_yi : float
    Distance from cluster y to cluster i
d_xy : float
    Distance from cluster x to cluster y
size_x : int
    Size of cluster x
size_y : int
    Size of cluster y
size_i : int
    Size of cluster i

Returns
-------
d_xyi : float
    Distance from the new cluster xy to cluster i
"""
ctypedef float (*linkage_distance_update)(float d_xi, float d_yi,
                                           float d_xy, int size_x,
                                           int size_y, int size_i)


cdef float _single(float d_xi, float d_yi, float d_xy,
                    int size_x, int size_y, int size_i):
    return min(d_xi, d_yi)


cdef float _complete(float d_xi, float d_yi, float d_xy,
                      int size_x, int size_y, int size_i):
    return max(d_xi, d_yi)


cdef float _average(float d_xi, float d_yi, float d_xy,
                     int size_x, int size_y, int size_i):
    return (size_x * d_xi + size_y * d_yi) / (size_x + size_y)


cdef float _centroid(float d_xi, float d_yi, float d_xy,
                      int size_x, int size_y, int size_i):
    return sqrt((((size_x * d_xi * d_xi) + (size_y * d_yi * d_yi)) -
                 (size_x * size_y * d_xy * d_xy) / (size_x + size_y)) /
                (size_x + size_y))


cdef float _median(float d_xi, float d_yi, float d_xy,
                    int size_x, int size_y, int size_i):
    return sqrt(0.5 * (d_xi * d_xi + d_yi * d_yi) - 0.25 * d_xy * d_xy)


cdef float _ward(float d_xi, float d_yi, float d_xy,
                  int size_x, int size_y, int size_i):
    cdef float t = 1.0 / (size_x + size_y + size_i)
    return sqrt((size_i + size_x) * t * d_xi * d_xi +
                (size_i + size_y) * t * d_yi * d_yi -
                size_i * t * d_xy * d_xy)


cdef float _weighted(float d_xi, float d_yi, float d_xy,
                      int size_x, int size_y, int size_i):
    return 0.5 * (d_xi + d_yi)
