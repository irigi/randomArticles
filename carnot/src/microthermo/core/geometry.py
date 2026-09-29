from __future__ import annotations

import math
import numpy as np


def cross2(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def rotate(points: np.ndarray, angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return points @ np.array([[c, s], [-s, c]])


def world_polygon(local: np.ndarray, position: np.ndarray, angle: float) -> np.ndarray:
    return rotate(local, angle) + position


def polygon_signed_area(vertices: np.ndarray) -> float:
    return 0.5 * float(np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1) -
                              vertices[:, 1] * np.roll(vertices[:, 0], -1)))


def ensure_ccw(vertices: np.ndarray) -> np.ndarray:
    return vertices if polygon_signed_area(vertices) > 0 else vertices[::-1].copy()


def closest_point_segment(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, float]:
    d = b - a
    dd = float(d @ d)
    if dd == 0.0:
        return a.copy(), 0.0
    t = max(0.0, min(1.0, float((p - a) @ d) / dd))
    return a + t * d, t


def point_in_polygon(point: np.ndarray, polygon: np.ndarray) -> bool:
    inside = False
    j = len(polygon) - 1
    for i in range(len(polygon)):
        pi, pj = polygon[i], polygon[j]
        if ((pi[1] > point[1]) != (pj[1] > point[1])):
            x = (pj[0] - pi[0]) * (point[1] - pi[1]) / (pj[1] - pi[1]) + pi[0]
            if point[0] < x:
                inside = not inside
        j = i
    return inside


def polygon_axes(poly: np.ndarray) -> list[np.ndarray]:
    edges = np.roll(poly, -1, axis=0) - poly
    axes = np.column_stack((-edges[:, 1], edges[:, 0]))
    norms = np.linalg.norm(axes, axis=1)
    return [axes[i] / norms[i] for i in range(len(poly)) if norms[i] > 0]


def convex_separation(a: np.ndarray, b: np.ndarray) -> tuple[float, np.ndarray]:
    """SAT separation; positive means separated, nonpositive means overlap."""
    best_sep = -math.inf
    best_axis = np.array([1.0, 0.0])
    ca, cb = np.mean(a, axis=0), np.mean(b, axis=0)
    for axis in polygon_axes(a) + polygon_axes(b):
        if (cb - ca) @ axis < 0:
            axis = -axis
        pa, pb = a @ axis, b @ axis
        sep = float(np.min(pb) - np.max(pa))
        if sep > best_sep:
            best_sep, best_axis = sep, axis
    return best_sep, best_axis


def convex_witnesses(a: np.ndarray, b: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Separation lower bound and closest physical features of two convex polygons.

    The SAT value is a conservative lower bound while separated. Witnesses are
    obtained from vertex-to-edge queries, including endpoint/vertex cases.
    Feature IDs encode polygon vertex or edge by ``2*index`` or ``2*index+1``.
    """
    separation, axis = convex_separation(a, b)
    best = math.inf
    wa = wb = None
    fa = fb = -1
    for i, vertex in enumerate(a):
        for j, start in enumerate(b):
            candidate, u = closest_point_segment(vertex, start, b[(j+1) % len(b)])
            d2 = float(np.dot(candidate-vertex, candidate-vertex))
            if d2 < best:
                best, wa, wb, fa = d2, vertex.copy(), candidate, 2*i
                fb = 2*j if u == 0.0 else (2*((j+1) % len(b)) if u == 1.0 else 2*j+1)
    for j, vertex in enumerate(b):
        for i, start in enumerate(a):
            candidate, u = closest_point_segment(vertex, start, a[(i+1) % len(a)])
            d2 = float(np.dot(candidate-vertex, candidate-vertex))
            if d2 < best:
                best, wa, wb, fb = d2, candidate, vertex.copy(), 2*j
                fa = 2*i if u == 0.0 else (2*((i+1) % len(a)) if u == 1.0 else 2*i+1)
    delta = wb-wa
    length = math.sqrt(best)
    if separation > 0.0 and fb % 2 == 1:
        edge=b[(fb//2+1)%len(b)]-b[fb//2]
        orientation=1.0 if polygon_signed_area(b)>0 else -1.0
        normal=orientation*np.array([-edge[1],edge[0]])/np.linalg.norm(edge)
    elif separation > 0.0 and fa % 2 == 1:
        edge=a[(fa//2+1)%len(a)]-a[fa//2]
        orientation=1.0 if polygon_signed_area(a)>0 else -1.0
        normal=orientation*np.array([edge[1],-edge[0]])/np.linalg.norm(edge)
    else:
        normal = delta/length if separation > 0.0 and length > 0.0 else axis
    return separation, normal, wa, wb, fa, fb


def polygon_segment_witnesses(polygon: np.ndarray, start: np.ndarray,
                              end: np.ndarray) -> tuple[float, np.ndarray, np.ndarray,
                                                        np.ndarray, int, int]:
    """Closest points on a convex polygon and finite segment.

    The distance is zero when the segment intersects or lies inside the
    polygon. Feature IDs use even numbers for vertices/endpoints and odd
    numbers for edge interiors.
    """
    start, end = np.asarray(start, float), np.asarray(end, float)
    segment = end-start
    best = math.inf
    wp = ws = None
    fp = fs = -1

    def offer(poly_point, segment_point, poly_feature, segment_feature):
        nonlocal best, wp, ws, fp, fs
        d2 = float(np.dot(segment_point-poly_point,segment_point-poly_point))
        if d2 < best:
            best, wp, ws = d2, poly_point.copy(), segment_point.copy()
            fp, fs = poly_feature, segment_feature

    for i, vertex in enumerate(polygon):
        point, u = closest_point_segment(vertex,start,end)
        offer(vertex,point,2*i,0 if u == 0.0 else 2 if u == 1.0 else 1)
    for j, endpoint in enumerate((start,end)):
        for i, vertex in enumerate(polygon):
            point, u = closest_point_segment(endpoint,vertex,polygon[(i+1)%len(polygon)])
            feature = 2*i if u == 0.0 else 2*((i+1)%len(polygon)) if u == 1.0 else 2*i+1
            offer(point,endpoint,feature,2*j)

    # Endpoint-edge queries alone miss a segment crossing the polygon while
    # both endpoints lie outside it. Detect that case before using the gap as
    # a conservative CCD distance.
    if point_in_polygon(start,polygon):
        offer(start,start,-1,0)
    if point_in_polygon(end,polygon):
        offer(end,end,-1,2)
    for i, vertex in enumerate(polygon):
        edge=polygon[(i+1)%len(polygon)]-vertex
        denominator=cross2(segment,edge)
        if abs(denominator) <= 1e-15*np.linalg.norm(segment)*np.linalg.norm(edge):
            continue
        delta=vertex-start
        u=cross2(delta,edge)/denominator
        v=cross2(delta,segment)/denominator
        if 0.0 <= u <= 1.0 and 0.0 <= v <= 1.0:
            point=start+u*segment
            offer(point,point,2*i+1,1)

    distance=math.sqrt(best)
    if distance > 0.0:
        normal=(ws-wp)/distance
    else:
        tangent=segment/np.linalg.norm(segment)
        normal=np.array([-tangent[1],tangent[0]])
        if (np.mean(polygon,axis=0)-.5*(start+end))@normal > 0:
            normal=-normal
    return distance,normal,wp,ws,fp,fs


def disc_polygon_separation(center: np.ndarray, radius: float, polygon: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, int]:
    """Signed separation from a disc to a convex polygon with outward normal."""
    polygon = ensure_ccw(polygon)
    best_d2 = math.inf
    best_point = polygon[0]
    best_edge = 0
    for i, a in enumerate(polygon):
        q, _ = closest_point_segment(center, a, polygon[(i + 1) % len(polygon)])
        d2 = float((center - q) @ (center - q))
        if d2 < best_d2:
            best_d2, best_point, best_edge = d2, q, i
    delta = best_point - center
    dist = math.sqrt(best_d2)
    if dist > 0:
        normal = delta / dist
    else:
        edge = polygon[(best_edge + 1) % len(polygon)] - polygon[best_edge]
        normal = np.array([edge[1], -edge[0]]) / np.linalg.norm(edge)
    inside = point_in_polygon(center, polygon)
    signed = -dist - radius if inside else dist - radius
    return signed, normal, best_point, best_edge


def polygon_inertia(vertices: np.ndarray, mass: float) -> tuple[np.ndarray, float]:
    """Return centroid and central inertia of a uniform simple polygon."""
    v = ensure_ccw(np.asarray(vertices, float))
    w = np.roll(v, -1, axis=0)
    cr = v[:, 0] * w[:, 1] - w[:, 0] * v[:, 1]
    area2 = float(np.sum(cr))
    centroid = np.sum((v + w) * cr[:, None], axis=0) / (3.0 * area2)
    second = np.sum(cr * (np.sum(v*v, axis=1) + np.sum(v*w, axis=1) + np.sum(w*w, axis=1))) / 12.0
    inertia_origin = mass * second / (area2 / 2.0)
    return centroid, float(inertia_origin - mass * (centroid @ centroid))


def aperture_clearance(width: float, disc_radius: float) -> float:
    return width - 2.0 * disc_radius
