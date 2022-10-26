import numpy as np
import scipy.interpolate


class QuadrilateralTranform:
    """Finds a bilinear transformation from an
    arbitrary quadrilateral to a unit square"""

    """wish list: test and fix quadrilaterals
    presented in the wrong order."""
    """wish list: test to make sure points inside
    a quadrilateral never give an error
                  test non-convex quadrilaterals"""

    def __init__(
        self, x_vertices, y_vertices  # array/list of 4 x locations
    ):  # array/list of 4 y locations

        # A = np.zeros((4,4),dtype = np.float64)
        # A[0,:] = [1.0,0.0,0.0,0.0]
        # A[1,:] = [1.0,1.0,0.0,0.0]
        # A[2,:] = [1.0,1.0,1.0,1.0]
        # A[3,:] = [1.0,0.0,1.0,0.0]

        # we seek a bilinear tranformation from x,y to l,m,
        # where l ranges from 0 to 1
        # and m ranges from 0 to 1

        # The transform from l,m to x,y is given by

        # x = a0 + a1*l + a2*m + a4*l*m
        # y = b0 + b1*l + b2*m + b4*l*m

        # X = A*a and Y = A*b
        # A_inv*X = a and A_inv*Y = b

        # A_inv tranforms X and Y to alpha and beta

        # We don't need to recalculate A_inv because is is constant

        A_inv = np.zeros((4, 4), dtype=np.float64)
        A_inv[0, :] = [1.0, 0.0, 0.0, 0.0]
        A_inv[1, :] = [-1.0, 1.0, 0.0, 0.0]
        A_inv[2, :] = [-1.0, 0.0, 0.0, 1.0]
        A_inv[3, :] = [1.0, -1.0, 1.0, -1.0]

        self.a = np.matmul(A_inv, x_vertices)
        self.b = np.matmul(A_inv, y_vertices)

        # Now the we know the a's and b's we need to go back to solve
        #
        # x = a0 + a1*l + a2*m + a4*l*m
        # y = b0 + b1*l + b2*m + b4*l*m
        #
        # for l and m given x and y.
        #
        # This can be rearranged into a quadratic for m
        #
        # aa = a3*b2 - a2*b3 (independent of x, y)
        # bb = a3*b0 - a0*b3  + a1*b2 - a2*b1 + b3*x - a3*y
        # cc = a1*b0 - a0*b1 + b1*x -a1*y
        #
        # here we compute aa and the x adn y independent parts of bb and cc

        self.aa = self.a[3] * self.b[2] - self.a[2] * self.b[3]
        # bb = a(4)*b(1) -a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + x*b(4) - y*a(4)
        self.bb_part = (
            self.a[3] * self.b[0]
            - self.a[0] * self.b[3]
            + self.a[1] * self.b[2]
            - self.a[2] * self.b[1]
        )
        self.cc_part = self.a[1] * self.b[0] - self.a[0] * self.b[1]

    def __call__(self, x, y):

        bb = self.bb_part + x * self.b[3] - y * self.a[3]
        cc = self.cc_part + x * self.b[1] - y * self.a[1]

        if self.aa == 0.0:
            try:
                m = -cc / bb
            except ZeroDivisionError:
                raise ValueError("x,y not inside quadrilateral")
        else:
            z = bb * bb - 4 * self.aa * cc
            if z >= 0.0:
                det = np.sqrt(z)
                m = (-bb + det) / (2 * self.aa)
            else:
                raise ValueError("x,y not inside quadrilateral")

        ll = (x - self.a[0] - self.a[2] * m) / (self.a[1] + self.a[3] * m)

        if (ll < 0.0) or (ll > 1.0) or (m < 0.0) or (m > 1.0):
            raise ValueError("x,y not inside quadrilateral")

        return ll, m


class InterpBilinearQuadrilateral:

    """defines and performs bilinear interpolation on an arbitrary quadrilateral"""

    def __init__(
        self,
        x_vertices,  # array/list of 4 x locations
        y_vertices,  # array/list of 4 y locations
        z,
    ):  # array/list of 4 z locations)

        self.tranform = QuadrilateralTranform(x_vertices, y_vertices)

        self.a00 = z[0]
        self.a10 = z[1] - z[0]
        self.a01 = z[3] - z[0]
        self.a11 = z[2] - z[1] - z[3] + z[0]

    def __call__(self, x, y):

        # transform to l,m
        l, m = self.transform(x, y)
        return self.a00 + self.a10 * l + self.a01 * m + self.a11 * m * l


class InterpParallelogram:
    """Interpolation inside an arbitratry parallelogram
    using a remapped bilinear interpolator"""

    def __init__(
        self, x, y, z  # array/list of 4 x locations  # array/list of 4 y locations
    ):  # array/list of 4 z locations)

        try:
            assert len(x) == 4
            assert len(y) == 4
            assert len(z) == 4
        except AssertionError:
            raise ValueError("x, y, and z must each have 4 elements")

        self.x = np.array(x, dtype=np.float64)
        self.y = np.array(y, dtype=np.float64)
        self.z = np.array(z, dtype=np.float64)

        self.vec01 = np.array([self.x[1] - self.x[0], self.y[1] - self.y[0]])
        self.vec32 = np.array([self.x[2] - self.x[3], self.y[2] - self.y[3]])
        self.vec03 = np.array([self.x[3] - self.x[0], self.y[3] - self.y[0]])
        self.vec12 = np.array([self.x[2] - self.x[1], self.y[2] - self.y[1]])

        self.len_01_sqr = self.vec01[0] ** 2 + self.vec01[1] ** 2
        self.len_03_sqr = self.vec03[0] ** 2 + self.vec03[1] ** 2

        self.a = np.array([z[0], z[1] - z[0], z[3] - z[0],
                           z[2] - z[1] - z[3] + z[0]])
        # the interpolation is done on the unit square, where the
        # interpolated value is given by
        # a00 + a10*x + a01*y + a11*x*y
        # a00 = z00 = z[0]
        # a10 = z10 - z00 = z[1]-z[0]
        # a01 = z01 - z00 = z[3]-z[0]
        # a11 = z11 - z10 - z01 +z00 = z[2] - z[1] - z[3] + z[0]
        # https://en.wikipedia.org/wiki/Bilinear_interpolation

        # need to find the transform that transforms the
        # parallelogram onto the unit square
        #
        #  | a b | |vec01x|    | 1 |
        #  |     |x|      | =  |   |
        #  | c d | |vec01y|    | 0 |
        #
        # and
        #
        #  | a b | |vec03x|    | 0 |
        #  |     |x|      | =  |   |
        #  | c d | |vec03y|    | 1 |

        A = np.zeros((4, 4), dtype=np.float64)
        A[0:2, 0] = self.vec01
        A[0:2, 1] = self.vec03
        A[2:4, 2] = self.vec01
        A[2:4, 3] = self.vec03

        B = np.array([1.0, 0.0, 0.0, 1.0], dtype=np.float64)

        solution = np.linalg.solve(A, B)

        self.transform_matrix = np.transpose(np.reshape(solution, (2, 2)))

        if not self.check_if_parallelogram():
            raise ValueError("points (x,y) are not a parallelogram")

    def check_if_parallelogram(self, rel_threshold=0.005):
        """Checks to see if the points in self represent a parallelogram"""
        """Points need to be in either clockwise or counterclockwise order"""

        threshold = (
            rel_threshold * 0.5 * (np.sqrt(self.len_01_sqr) + np.sqrt(self.len_03_sqr))
        )

        is_parallelogram = True

        # 01 and 32 should be the same
        max_abs_diff = np.max(np.abs(self.vec01 - self.vec32))
        if max_abs_diff > threshold:
            is_parallelogram = False

        # 03 and 12 should be the same
        max_abs_diff = np.max(np.abs(self.vec03 - self.vec12))
        if max_abs_diff > threshold:
            is_parallelogram = False

        return is_parallelogram

    def __call__(self, x0, y0, check_for_inside=True):

        dx_y = np.array([x0 - self.x[0], y0 - self.y[0]])

        # This transforms the locations dx,dy (relative to x[0],y[0])
        # onto the unit square
        uv = np.dot(self.transform_matrix, dx_y)

        if check_for_inside:
            if (np.max(uv) > 1.000000001) or (np.min(uv) < -0.000000001):
                raise ValueError("point is not inside the parallelogram")

        # this does the bilinear interpolation
        uv0 = np.array([1, uv[0], uv[1], uv[0] * uv[1]])
        return np.dot(self.a, uv0)


class InterpQuadrilateralFit:
    """Interpolation in an arbitratry quadrilateral using a minimum curvature bi-quadratic
    interpolator"""

    """The function is in the form ax^2 + bxy + cy^2 + dx + ey + 1.
       The function is equal to the values z and the four corners, and the
       quadratic terms are minimized using a lagrange multiplier method"""

    """Based on the discussion in
       https://math.stackexchange.com/questions/828392/
       spatial-interpolation-for-irregular-grid"""

    def __init__(
        self, x, y, z  # array/list of 4 x locations  # array/list of 4 y locations
    ):  # array/list of 4 z locations
        """defines the interpolation function"""

        try:
            assert len(x) == 4
            assert len(y) == 4
            assert len(z) == 4
        except AssertionError:
            raise ValueError("x, y, and z must each have 4 elements")

        self.x = np.array(x, dtype=np.float64)
        self.y = np.array(y, dtype=np.float64)
        self.z = np.array(z, dtype=np.float64)

        X = np.zeros((6, 4), dtype=np.float64)

        X[0, :] = np.square(self.x)
        X[1, :] = self.x * self.y
        X[2, :] = np.square(self.y)
        X[3, :] = self.x
        X[4, :] = self.y
        X[5, :] = 1.0

        E = np.diag([1, 1, 1, 0, 0, 0]).astype(np.float64)

        A = np.zeros((10, 10), dtype=np.float64)

        A[0:6, 0:6] = E
        A[6:10, 0:6] = np.transpose(X)
        A[0:6, 6:10] = X

        B = np.zeros((10), dtype=np.float64)
        B[6:10] = z

        solution = np.linalg.solve(A, B)

        self.parameters = solution[0:6]
        # the other 4 values are the lagrange multipliers
        self.lagrange = solution[6:10]

    def __call__(self, x, y):
        """performs the interpolation"""
        return (
            self.parameters[0] * np.square(x)
            + self.parameters[1] * x * y
            + self.parameters[2] * np.square(y)
            + self.parameters[3] * x
            + self.parameters[4] * y
            + self.parameters[5]
        )

    def check_inside(self, x, y):
        """checks to see if x,y is inside the quadrilateral
        defined for the interpolator"""

        # Not implemented yet
        # I know how to do it for a convex quadrilateral,
        # but not for one that is not convex
        # probably it is not very good to do this for a
        # non convex quadrilateral anyway....
        pass


if __name__ == "__main__":

    a_small_number = 0.00000001
    x = [1.0, 3.0, 3.5, 1.5]
    y = [1.0, 2.0, 3.0, 3.0]
    z = [2.0, 3.0, 6.0, 3.0]

    qt = QuadrilateralTranform(x, y)

    x0 = 1.0
    y0 = 1.0
    print(qt(x0, y0))

    x0 = 3.0
    y0 = 2.0
    print(qt(x0, y0))

    x0 = 1.5
    y0 = 3.0
    print(qt(x0, y0))

    x0 = 3.5
    y0 = 3.0
    print(qt(x0, y0))

    x0 = 3.25
    y0 = 2.5
    print(qt(x0, y0))

    qi = InterpQuadrilateralFit(x, y, z)
    pi = InterpParallelogram(x, y, z)

    sci = scipy.interpolate.interp2d(x, y, z)

    print("General Quadrilateral Method on a parallelogram")

    for i in range(4):
        print(f"Interpolated Value ax {x[i]},{y[i]} = {qi(x[i],y[i])} should be {z[i]}")
        if np.abs(qi(x[i], y[i]) - z[i]) > a_small_number:
            raise Exception(
                "Fail - error in interploated value at a corner of the quadrilateral"
            )
    inside_test_point = [2.0, 1.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{qi(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [2.25, 2.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{qi(inside_test_point[0],inside_test_point[1])}"
    )

    print()
    print("Parallelogram Method on a parallelogram")

    for i in range(4):
        print(f"Interpolated Value ax {x[i]},{y[i]} = {pi(x[i],y[i])} should be {z[i]}")
        if np.abs(pi(x[i], y[i]) - z[i]) > a_small_number:
            raise Exception(
                "Fail - error in interploated value at a corner of the quadrilateral"
            )
    inside_test_point = [2.0, 1.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{pi(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [2.5, 3.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{pi(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [2.5, 1.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{pi(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [3.25, 2.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{pi(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [1.25, 2.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{pi(inside_test_point[0],inside_test_point[1])}"
    )

    print()
    print("SciPy Method on a parallelogram")
    for i in range(4):
        print(
            f"Interpolated Value ax {x[i]},{y[i]} = " +
            f"{sci(x[i],y[i])} should be {z[i]}"
        )
        if np.abs(sci(x[i], y[i]) - z[i]) > a_small_number:
            raise Exception(
                "Fail - error in interploated value at a corner " +
                "of the quadrilateral"
            )

    inside_test_point = [2.0, 1.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{sci(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [2.5, 3.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{sci(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [2.5, 1.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{sci(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [3.25, 2.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} = " +
        f"{sci(inside_test_point[0],inside_test_point[1])}"
    )
    inside_test_point = [1.25, 2.0]
    print(
        f"Interpolated Value {inside_test_point[0]:.2f}," +
        f"{inside_test_point[1]:.2f} =" +
        f" {sci(inside_test_point[0],inside_test_point[1])}"
    )

    print()

    outside_test_point = [3.5, 2.7]

    try:
        print(
            f"Interpolated Value {outside_test_point[0]:.2f}," +
            f"{outside_test_point[1]:.2f} = " +
            f"{pi(outside_test_point[0],inside_test_point[1])}"
        )
        raise Exception("Fail - should have caused an error")
    except ValueError as e:
        print(e)
        print("Pass - this should be an error")

    # NOT a parallelogram
    print()
    print("General Quadrilateral Method on a NOT a parallelogram")
    x = [1.0, 3.0, 3.5, 1.5]
    y = [1.0, 1.0, 3.0, 4.0]
    z = [2.0, 3.0, 4.0, 3.0]

    qi = InterpQuadrilateralFit(x, y, z)
    print("General Quaudrilateral")

    for i in range(4):
        print(f"Interpolated Value ax {x[i]},{y[i]} = " +
              f"{qi(x[i],y[i])} should be {z[i]}")
    print(f"Interpolated Value = {qi(2.0,2.0)}")

    # this should cause an error
    print()
    print("Parallelogram Method on a NOT a parallelogram")
    try:
        pi = InterpParallelogram(x, y, z)
        raise Exception("Fail - should have caused an error")
    except ValueError as e:
        print(e)
        print("Pass - Error should have happened")

    print("All tests passed")
