from cryptography318.linear_algebra import augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals
from cryptography318.matrix import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
import textwrap


def testLinearSolve(it=500, count=False):
    acceptable = 0
    for _ in range(it):
        m1 = RandomMatrix()
        m2 = RandomMatrix(rows=len(m1[0]), cols=1)
        m3 = MultiplyMatrix(m1, m2)
        matrix = augRREF(augmentMatrix(m1, m3))
        if IsConsistent(matrix) and IsSolvable(matrix, aug=True):
            acceptable += 1
            a1, a2 = separateMatrix(matrix)
            res = Solve(a1, sol=a2)
            assert MatrixEquals(res, m2)

    if count:
        print(f"\nLinearSolve test tried to solve {it} random linear systems.")
        print(textwrap.indent(f"Consistent and solvable systems: {acceptable}\n"
              f"Solved: {acceptable}\n", ' ' * 5))


def testAllTests():
    testLinearSolve(count=True)


if __name__ == '__main__':
    testAllTests()
