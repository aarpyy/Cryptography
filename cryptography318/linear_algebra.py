import numpy
from matrix import RREF, IsRREF, augmentMatrix, augRREF, augIsRREF


def IsConsistent(matrix):
    """Function that takes as input an augmented coefficient matrix and returns
    True if the matrix is consistent, False if otherwise"""

    c = len(matrix[0]) - 1
    for i in range(len(matrix)):
        non_zero = 0
        for j in range(c):
            if matrix[i][j] != 0:
                non_zero += 1
        if non_zero == 0 and matrix[i][c] != 0:
            print(non_zero, i)
            return False
    return True


def _solveSystem(matrix, sol=None, aug=False):
    if sol is None and not aug:
        sol = numpy.array([[0]] * len(matrix[0]))

    if not aug:
        matrix = augmentMatrix(matrix, sol)

    if not augIsRREF(matrix):
        matrix = augRREF(matrix)

    r = len(matrix[0]) - 1

    _vars = {}
    if len(matrix) == 1:
        for i in range(r):
            if matrix[0][i] != 0:
                _vars[i] = None
        if len(_vars) > 1:
            return False

        for x in _vars:
            _vars[x] = matrix[0][r]
        return _vars

    result = _solveSystem(matrix[1:], sol=None, aug=True)
    if result is False:
        return False

    unknowns = []
    for i in range(r):
        if matrix[0][i] != 0:
            if i not in result:
                unknowns.append(i)

    if len(unknowns) > 1:
        return False
    if len(unknowns) == 0:
        return result

    solution = 0
    index = None
    for i in unknowns:
        if i in result:
            solution += result[i] * matrix[0][1]
        else:
            index = i

    result[index] = matrix[0][r] - solution
    return result


def Solve(matrix, sol=None, aug=False):
    if not IsSolveable(matrix, aug):
        raise ArithmeticError("System could not be solved")

    result = _solveSystem(matrix, sol, aug)

    solution = []
    for var in sorted(result):
        solution.append([result[var]])

    return numpy.array(solution)


def IsSolveable(matrix, aug=False):
    if aug:
        if not augIsRREF(matrix):
            matrix = augRREF(matrix)
        if not IsConsistent(matrix):
            raise ArithmeticError("System is inconsistent")
    else:
        if not IsRREF(matrix):
            matrix = RREF(matrix)

    r = len(matrix[0])
    if aug:
        r -= 1

    for j in range(r):
        for i in range(len(matrix)):
            if matrix[i][j] not in [0, 1]:
                raise ArithmeticError("System has free variables")
    return True
