#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math


def minimum_distance(left_atoms, right_atoms):
    output = atom_to_atom_distance(left_atoms[0], right_atoms[0])
    for left in left_atoms:
        for right in right_atoms:
            output = min(output, atom_to_atom_distance(left, right))
    return output


def atom_to_atom_distance(left_atom, right_atom):
    distances = [left_atom.x - right_atom.x,
                 left_atom.y - right_atom.y,
                 left_atom.z - right_atom.z]
    return math.sqrt(sum([x * x for x in distances]))


def point_to_line_distance(source, target, point):
    source_target = vector_diff(source, target)
    source_point = vector_diff(source, point)
    product = vector_product(source_point, source_target)
    return vector_size(product) / vector_size(source_target)


def atom_as_vector(atom):
    return [atom.x, atom.y, atom.z]


def vector_dot_product(left, right):
    return (left[0] * right[0]) + (left[1] * right[1]) + (left[2] * right[2])


def vector_product(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
    ]


def vector_diff(left, right):
    return [l - r for l, r in zip(left, right)]


def vector_normalized(vector):
    size = vector_size(vector)
    return [value / size for value in vector]


def vector_size(vector):
    return math.sqrt((vector[0] * vector[0]) +
                     (vector[1] * vector[1]) +
                     (vector[2] * vector[2]))


def do_sphere_obstruct_connection(
        start, end, sphere_center, sphere_radius, ignore_in_radius=None):
    """
    :return: (is there collision, distance to collision).
    """
    sphere_center = vector_diff(sphere_center, start)
    end = vector_diff(end, start)
    start = vector_diff(start, start)

    direction = vector_normalized(vector_diff(end, start))
    projection = multiply_by_scalar(
        direction, vector_dot_product(direction, sphere_center))

    # If the center of sphere is of a greater distance then radius
    # then there is no obstruction.
    distance_to_line = vector_size(vector_diff(projection, sphere_center))
    if distance_to_line > sphere_radius:
        return False, None

    # We could use vector size as the start is shifted to the origin.
    projection_distance = vector_size(projection)

    # We need to calculate where does the ray enter the sphere,
    # call this point X.
    projection_to_x = math.sqrt((sphere_radius * sphere_radius) -
                                (distance_to_line * distance_to_line))
    distance_to_x = projection_distance - projection_to_x

    distance_to_end = vector_size(end)
    if distance_to_end < distance_to_x:
        # The collision is after the end of line.
        return False, distance_to_x

    if ignore_in_radius is not None:
        if distance_to_end > ignore_in_radius:
            # Outside ignore radius, this is OK.
            return True, distance_to_x
        # Inside ignore radius, try the other collision point.
        distance_to_x = projection_distance + projection_to_x
        if distance_to_end > ignore_in_radius:
            # Nope still in radius.
            return True, distance_to_x
        if distance_to_end < distance_to_x:
            # The collision is after the end of line.
            return False, distance_to_x

    return True, distance_to_x


def multiply_by_scalar(vector, scalar):
    return [value * scalar for value in vector]


def vector_add(left, right):
    return [l + r for l, r in zip(left, right)]


def rotation_matrix_from_vector_and_angle(v, angle_in_radian):
    cos = math.cos(angle_in_radian)
    sin = math.sin(angle_in_radian)
    x = v[0]
    y = v[1]
    z = v[2]
    return [
        [cos + (math.pow(x, 2) * (1 - cos)),
         (y * x * (1 - cos)) + z * sin,
         (z * x * (1 - cos)) - y * sin
         ], [
            (x * y * (1 - cos)) - z * sin,
            cos + (math.pow(y, 2) * (1 - cos)),
            (z * y * (1 - cos)) + x * sin
        ], [
            (x * z * (1 - cos)) + y * sin,
            (y * z * (1 - cos)) - x * sin,
            cos + (math.pow(z, 2) * (1 - cos))
        ]
    ]


def multiply_by_matrix(vector, matrix):
    return [
        sum([v * c for v, c in zip(vector, column)])
        for column in matrix
    ]
