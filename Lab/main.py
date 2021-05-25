# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

import math


class Point:
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return "(" + str(self.x) + " ; " + str(self.y) + ")"

class HermiteSpline:
    class KeyPoint:
        Input = 0.0
        Output = 0.0
        M = 0.0

        def __str__(self):
            return '[' + str(self.Input) + ', ' + str(self.Output) + ', ' + str(self.M) + ']'

    class Param:
        pass

    def __init__(self):
        self.index_previous = 0
        self.Param = self.Param()

    def FindIndex(self, t, index_previous=0):
        index = index_previous
        if index >= len(self.KeyPoints): index = len(self.KeyPoints) - 1
        while index + 1 < len(self.KeyPoints) and t > self.KeyPoints[index + 1].Input:  index += 1
        while index >= 0 and t < self.KeyPoints[index].Input:  index -= 1
        return index

    def Evaluate(self, t):
        index = self.FindIndex(t, self.index_previous)
        if abs(t - self.KeyPoints[-1].Input) < 1.0e-6:  index = len(self.KeyPoints) - 2
        if index < 0 or index >= len(self.KeyPoints) - 1:
            if index < 0:
                index = 0
                t = self.KeyPoints[0].Input
            else:
                index = len(self.KeyPoints) - 2
                t = self.KeyPoints[-1].Input

        h00 = lambda t: t * t * (2.0 * t - 3.0) + 1.0
        h10 = lambda t: t * (t * (t - 2.0) + 1.0)
        h01 = lambda t: t * t * (-2.0 * t + 3.0)
        h11 = lambda t: t * t * (t - 1.0)

        self.index_previous = index
        p0 = self.KeyPoints[index]
        p1 = self.KeyPoints[index + 1]
        tr = (t - p0.Input) / (p1.Input - p0.Input)
        return h00(tr) * p0.Output + h10(tr) * (p1.Input - p0.Input) * p0.M + h01(tr) * p1.Output + h11(tr) * (p1.Input - p0.Input) * p1.M

    CARDINAL = 1
    GRAD = 1

    def StartCardinal(self, data, tan_method=CARDINAL, end_tan=GRAD, c=0.0, m=1.0):
        if data != None:
            self.KeyPoints = [self.KeyPoint() for i in range(len(data))]
            for idx in range(len(data)):
                self.KeyPoints[idx].Input = data[idx][0]
                self.KeyPoints[idx].Output = data[idx][1]

        if tan_method == None: tan_method = self.Param.TanMethod
        else: self.Param.TanMethod = tan_method
        if end_tan == None: end_tan = self.Param.EndTan
        else: self.Param.EndTan = end_tan
        if c == None: c = self.Param.C
        else: self.Param.C = c
        if m == None: c = self.Param.M
        else: self.Param.M = m

        grad = lambda idx1, idx2: (self.KeyPoints[idx2].Output - self.KeyPoints[idx1].Output) / (
                self.KeyPoints[idx2].Input - self.KeyPoints[idx1].Input)

        if tan_method == self.CARDINAL:
            for idx in range(1, len(self.KeyPoints) - 1):
                self.KeyPoints[idx].M = (1.0 - c) * grad(idx - 1, idx + 1)

        if end_tan == self.GRAD:
            self.KeyPoints[0].M = m * grad(0, 1)
            self.KeyPoints[-1].M = m * grad(-2, -1)
        # if end_tan == self.ZERO:
        #     self.KeyPoints[0].M = 0.0
        #     self.KeyPoints[-1].M = 0.0
        # elif end_tan == self.CYCLIC:
        #     if tan_method == self.CARDINAL:
        #         T = self.KeyPoints[-1].Input - self.KeyPoints[0].Input
        #         X = self.KeyPoints[-1].Output - self.KeyPoints[0].Output
        #         grad_2 = (X + self.KeyPoints[1].Output - self.KeyPoints[-2].Output) / (T + self.KeyPoints[1].Input - self.KeyPoints[-2].Input)
        #         M = (1.0 - c) * grad_2
        #         self.KeyPoints[0].M = M
        #         self.KeyPoints[-1].M = M
        # ZERO = 0
        # CYCLIC = 2

def GetPoints(file_name):
    points_list = []
    input_array = open(file_name).read().split()

    i = 0
    while i < len(input_array):
        new_point = Point(float(input_array[i]), float(input_array[i + 1]))
        points_list.append(new_point)
        i += 2
    return points_list


def FindOrientation(a, b, c):
    value = (b.y - a.y) * (c.x - b.x) - (c.y - b.y) * (b.x - a.x)
    if value == 0:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1


def JarvisHull(points, n):
    if n <= 2:
        return

    result = []

    leftmost = 0
    p = leftmost
    q = -1
    while q != leftmost:
        result.append(points[p])
        q = GetNextIndex(points, p)
        for i in range(n):
            if FindOrientation(points[p], points[i], points[q]) == -1:
                q = i
        p = q

    return result


def GetNextIndex(array, current_index):
    n = len(array)
    if current_index == n - 1:
        return 0
    return current_index + 1


def FindHullApprox(points):
    n = len(points)

    left = points[0]
    right = points[n - 1]
    distance = right.x - left.x

    k = int(distance / 2)

    width = distance / k
    stripes = []
    plt.plot([left.x, left.x], [-10, 10])

    for i in range(k):
        stripe = []
        for point in points:
            if left.x + (i + 1) * width > point.x >= left.x + i * width:
                stripe.append(point)
        stripes.append(stripe)
        if i == k - 1:
            stripes[i].append(right)
        plt.plot([left.x + (i + 1) * width, left.x + (i + 1) * width], [-10, 10], '-b')

    plt.draw()
    potential_hull = []
    for stripe in stripes:
        if len(stripe) == 0:
            continue
        low = min(stripe, key=lambda point: point.y)
        high = max(stripe, key=lambda point: point.y)
        if potential_hull.count(low) == 0:
            potential_hull.append(low)
        if potential_hull.count(high) == 0:
            potential_hull.append(high)
    plt.draw()

    if potential_hull.count(left) == 0:
        potential_hull.append(left)
    if potential_hull.count(right) == 0:
        potential_hull.append(right)

    return JarvisHull(potential_hull, len(potential_hull))


def main():
    points = GetPoints("points.txt")
    points = sorted(points, key=lambda point: point.x)
    print(points)

    xs = []
    ys = []
    for point in points:
        xs.append(point.x)
        ys.append(point.y)

    plt.plot(xs, ys, 'ko')
    plt.draw()

    hull = FindHullApprox(points)
    xs = []
    ys = []
    key_points = []
    for point in hull:
        key_points.append([point.x, point.y])
        xs.append(point.x)
        ys.append(point.y)

    xs.append(xs[0])
    ys.append(ys[0])

    flag = True
    points1 = []
    points2 = []
    for i in range(len(xs) - 1):
        if (xs[i] < xs[i + 1] and flag):
            points1.append([xs[i], ys[i]])
        elif (xs[i] < xs[i + 1] and not flag):
            if (xs[i - 1] > xs[i]):
                points2.append([xs[i], ys[i]])
            points1.insert(0, [xs[i], ys[i]])
        else:
            if (flag):
                points1.append([xs[i], ys[i]])
            points2.append([xs[i], ys[i]])
            flag = False

    spline1 = HermiteSpline()
    spline1.StartCardinal(points1, tan_method=spline1.CARDINAL, c=0.0)

    X = []
    Y = []
    import numpy as np
    for i in range(len(points1) - 1):
        for t in np.arange(points1[i][0], points1[i + 1][0], 0.001):
            x = spline1.Evaluate(t)
            X.append(t)
            Y.append(x)

    plt.plot(X, Y, '-m')

    X = []
    Y = []

    points2.sort()
    spline2 = HermiteSpline()
    spline2.StartCardinal(points2, tan_method=spline2.CARDINAL, c=0.0)

    for i in range(len(points2) - 1):
        for t in np.arange(points2[i][0], points2[i + 1][0], 0.001):
            x = spline2.Evaluate(t)
            X.append(t)
            Y.append(x)

    plt.plot(X, Y, '-m')

    plt.axis([-10, 10, -10, 10])
    plt.show()


if __name__ == '__main__':
    main()