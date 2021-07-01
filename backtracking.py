

def PlaceQueens(size, col, sofar):
    if len(sofar) == size:
        return True, sofar

    # iterate through all possible values in current column
    for row in range(size):
        curr = (row, col)
        invalid = False
        for pair in sofar:

            # if current row value of current column conflict with a previous placement, it is not a valid spot
            if pair[0] == curr[0] or pair[1] == curr[1] or (abs(pair[0] - curr[0]) == abs(pair[1] - curr[1])):
                invalid = True

        # if no current conflict exist, try placing a queen there
        if not invalid:
            sofar_copy = sofar[:]
            sofar_copy.append(curr)
            check = PlaceQueens(size, col+1, sofar_copy)
            # if all iterations return True, return output
            if check[0]:
                return check

    # if never found possible setup, returns false
    return False, []


def Draw(size, lst):
    for i in range(size):
        for j in range(size):
            if (i, j) in lst:
                print("Q", end=" ")
            else:
                print("_", end=" ")
        print()


n = 5

possible, lst = PlaceQueens(n, 0, [])
Draw(n, lst)
