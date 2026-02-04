'''
[[  x,  x,  x,  x, 15,  x,  x,  x,  x,  x,  x,  x,  x],
 [  x,  x,  x,  x,  x,  x,  x, 11,  x,  x,  x,  x,  x],
 [  x, 15,  x, 05,  x,  x, 15,  x, 11,  x, 11,  x,  x],
 [  x,  x,  x,  x, 15,  x,  x, 08,  x, 12,  x, 12,  x],
 [  x, 16,  x,  x,  x, 08,  x,  x,  x,  x, 06,  x,  x],
 [  x,  x,  x, 16,  x,  x,  x,  x,  x,  x,  x,  x, 06],
 [  x,  x, 16,  x, 03,  x, 16,  x, 01,  x, 12,  x,  x],
 [ 13,  x,  x,  x,  x,  x,  x,  x,  x, 04,  x,  x,  x],
 [  x,  x, 07,  x,  x,  x,  x, 12,  x,  x,  x, 10,  x],
 [  x, 02,  x, 13,  x, 16,  x,  x, 14,  x,  x,  x,  x],
 [  x,  x, 13,  x, 14,  x, 14,  x,  x, 14,  x, 10,  x],
 [  x,  x,  x,  x,  x, 09,  x,  x,  x,  x,  x,  x,  x],
 [  x,  x,  x,  x,  x,  x,  x,  x, 09,  x,  x,  x,  x]]
Strategy:
xSolve for a,b,c tiles
-12 and 13 seem somewhat restricted... how many possible versions are there to have for 13 to contain 12? 12 is seemingly more restricted

after:
[[  x, 05, 05, 05, 15, 15,  x, 11,  x,  x,  x,  x,  x],
 [  x,  x,  x, 05,  x, 15,  x, 11,  x,  x, 11, 11, 11],
 [ 15, 15, 15, 05,  x, 15, 15, 11, 11, 11, 11,  x, 11],
 [ 15, 16, 15, 15, 15, 15, 08, 08, 08, 12, 12, 12, 11],
 [ 15, 16,  x,  x, 08, 08, 08,  x, 08, 12, 06, 06, 06],
 [ 15, 16,  x, 16, 16, 16, 16,  x, 08, 12, 12,  x, 06],
 [  x, 16, 16, 16, 03, 03, 16, 16, 01, 04, 12, 06, 06],
 [ 13, 13, 13, 13, 14, 03, 16, 04, 04, 04, 12, 10, 10],
 [ 07, 07, 07, 13, 14, 16, 16, 12, 12, 12, 12, 10,  x],
 [ 07, 02, 02, 13, 14, 16, 14, 14, 14, 14,  x, 10,  x],
 [ 07, 07, 13, 13, 14, 14, 14,  x, 09, 14, 14, 10, 10],
 [  x, 07, 13, 09, 09, 09, 09,  x, 09, 14,  x,  x, 10],
 [  x,  x, 13, 13, 13, 13, 09, 09, 09, 14, 10, 10, 10]]
 
row sums: [56, 64, 135, 162, 93, 133, 115, 129, 138, 120, 139, 89, 123]
min: 56, max: 162

ANSWER: 9072
'''

from numpy import sqrt, log, arange
from copy import deepcopy as dc

class Grid:
    def __init__(self, r=13, c=13):
        self.r = r
        self.c = c
        
        self.grid = [[0]*c for _ in range(r)]
    
    def fillGrid(self, r, c, n):
        if n < 1 or n > 17:
            return -1
        self.grid[r][c] = n
    
    def getGridSpace(self, r, c):
        if r < 0 or r >= self.r or c < 0 or c >= self.c:
            #print("Warning: Grabbing invalid grid space")
            return -2
        n = self.grid[r][c]
        return n
    
    def printGrid(self):
        s = "\n[["
        for r in range(self.r):
            if r > 0:
                s += ",\n ["
            for c in range(self.c):
                n = self.getGridSpace(r, c)
                if n == 0:
                    s += "  x"
                elif 1 <= n <= 9:
                    s += " 0" + str(n)
                else:
                    s += " " + str(n)
                if c+1 < self.c:
                    s += ","
            s += "]"
        s += "]"
        
        print(s)
    
    def getFinalAns(self):
        def getRowSums():
            rs = [0]*self.r
            for r in range(self.r):
                for c in range(self.c):
                    n = self.getGridSpace(r, c)
                    if n > 0:
                        rs[r] += n
            
            return rs
        
        rowSums = getRowSums()
        mn, mx = min(rowSums), max(rowSums)
        ans = mn*mx
        print(f"\nrow sums: {rowSums}\nmin: {mn}, max: {mx}, ANSWER: {ans}\n")
        
        return ans
        
    def getNominoOfN(self, n, extras=True):
        spaces = []
        for r in range(self.r):
            for c in range(self.c):
                space = self.getGridSpace(r, c)
                if space == n:
                    spaces.append((r,c))
        
        nomino = Nomino(spaces, offset=extras, doVariations=extras)
        return nomino
    
    def isSameGrid(self, otherGrid):
        if otherGrid.r != self.r or otherGrid.c != otherGrid.c:
            return False
        
        for r in range(self.r):
            for c in range(self.c):
                n1 = self.getGridSpace(r, c)
                n2 = otherGrid.getGridSpace(r, c)
                if n1 != n2:
                    return False
        
        return True
        
        
        
class Nomino:
    def __init__(self, spaces, offset=True, doVariations=True):
        self.spaces = dc(spaces)
        self.n = len(spaces)

        self.variations = []

        if offset:
            self.offset()
        
        self.isConnected = self.connectedCheck()
        if not self.isConnected:
            return
        
        if doVariations:
            self.rotFlip()
    
    def offset(self):
        space0 = self.spaces[0]
        r0, c0 = space0
        if r0 == 0 and c0 == 0:
            return
        newSpaces = []
        for space in self.spaces:
            r, c = space
            newSpace = (r-r0, c-c0)
            newSpaces.append(newSpace)
        
        self.spaces = newSpaces
    
    def connectedCheck(self, spaces=None):
        if spaces is None:
            spaces = self.spaces
        clearedList = [0]*self.n
        clearedList[0] = 1
        
        prevChecked = [spaces[0]]
        oDirs = [(-1, 0), (0, 1), (1, 0), (0, -1)]  # nesw
        while len(prevChecked) > 0:
            nextChecked = []
            for prevSpace in prevChecked:
                r1, c1 = prevSpace
                for oDir in oDirs:
                    ro, co = oDir
                    nextSpace = (r1+ro, c1+co)
                    index = -1
                    for i, space in enumerate(spaces):
                        if space == nextSpace:
                            index = i
                            break
                    if index == -1:
                        continue
                    if clearedList[index] == 1:
                        continue
                    
                    clearedList[index] = 1
                    
                    nextChecked.append(nextSpace)
            
            prevChecked = nextChecked
        
        return 0 not in clearedList
    
    def rotFlip(self):
        self.variations = [[],[],[],[],[],[],[],[]]
        for space in self.spaces:
            r, c = space
            self.variations[0].append(space)
            self.variations[1].append((r, -c))
            self.variations[2].append((-r, c))
            self.variations[3].append((-r, -c))
            self.variations[4].append((c, r))
            self.variations[5].append((c, -r))
            self.variations[6].append((-c, r))
            self.variations[7].append((-c, -r))
            
        

class Solver:
    def __init__(self):
        self.eqsList = [
            lambda a, b, c: 6 * c - 4 * b, # 0
            lambda a, b, c: 8 - b,
            lambda a, b, c: (a ** b - 4) / (6 * c + 1),
            lambda a, b, c: (b + c) / (c - 1),
            lambda a, b, c: b * b - b / c,
            lambda a, b, c: sqrt(30 + a) / c, # 5
            lambda a, b, c: (a + b) / (c - 3 * a),
            lambda a, b, c: (b - 3 * a) / (a - c),
            lambda a, b, c: 8 * a - 2 * b,
            lambda a, b, c: b / (a-c),
            lambda a, b, c: (b+9)/sqrt(c-a), # 10
            lambda a, b, c: 18 / (a*c+1),
            lambda a, b, c: c**b,
            lambda a, b, c: (3+b*b)/sqrt(3+2*c),
            lambda a, b, c: b / (a*a-c*c),
            lambda a, b, c: sqrt(a+2) / a, # 15
            lambda a, b, c: a**b - 12/a,
            lambda a, b, c: 2*c + c/a,
            lambda a, b, c: 4*a - 5*b,
            lambda a, b, c: c+2*a,
            lambda a, b, c: b/(9*a-5*c), # 20
            lambda a, b, c: (b**3+2*c)/(b+2*c),
            lambda a, b, c: b/(a-1),
            lambda a, b, c: (c-b)/(2*a),
            lambda a, b, c: b/(a-c),
            lambda a, b, c: (b+c)/(a-c), # 25
            lambda a, b, c: log(a)/log(c),
            lambda a, b, c: (c*c-b)/a,
            lambda a, b, c: (b-1)**2,
            lambda a, b, c: ((43-a*c)**(1/3.))/a,
            lambda a, b, c: (b-a)/(a-c), # 30
            lambda a, b, c: 11-b,
            lambda a, b, c: (b-2*a)/(a-c),
            lambda a, b, c: (c+3)/a,
            lambda a, b, c: 8*c - b/c,
            lambda a, b, c: b**2, # 35
            lambda a, b, c: (2**b + 1) / (a*c) # 36
        ]
        self.eqsSols = [15, 11, 15, 5, 15, 11, 11, 15, 8, 12, 12, 16, 8, 6, 16, 6, 16, 3, 16, 1, 12, 13, 4, 7, 12, 10, 2, 13, 16, 14, 13, 14, 14, 14, 10, 9, 9]
        # Sorted: [1, 2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16, 16]
        self.eqPos = [
            (0,4), # 0
            (1,7),
            (2,1),
            (2,3),
            (2,6),
            (2,8), # 5
            (2,10),
            (3,4),
            (3,7),
            (3,9),
            (3,11), # 10
            (4,1),
            (4,5),
            (4,10),
            (5,3),
            (5,12), # 15
            (6,2),
            (6,4),
            (6,6),
            (6,8),
            (6,10), # 20
            (7,0),
            (7,9),
            (8,2),
            (8,7),
            (8,11), # 25
            (9,1),
            (9,3),
            (9,5),
            (9,8),
            (10,2), # 30
            (10,4),
            (10,6),
            (10,9),
            (10,11),
            (11,5), # 35
            (12,8) # 36
        ]
        
        self.grid = Grid(13, 13)
    
    def solve_equations(self):
        # Narrow down a
        def test1():
            eq15 = self.eqsList[15]
            for a in arange(-2, 3, 0.25):
                if a == 0:
                    continue
                x = eq15(a,0,0)
                print(a, x)
        
        # Narrow down c
        def test2():
            eq26 = self.eqsList[26]
            for a in [0.25, 2]:
                for x in range(1,18):
                    c = x-2*a
                    y = eq26(a, 0, c)
                    s_extra = "" if a != c else "... invalid because a cannot = c"
                    print(f"a: {a}, c: {c}, log_c(a): {y}{s_extra}")
        
        # Narrow down b
        def test3():
            a = 0.25
            c = 0.5
            for b in range(-4, 5):
                if b in [-1]:
                    continue
                eqs_sols = []
                for eqi, eq in enumerate(self.eqsList):
                    eq_sol = eq(a, b, c)
                    eqs_sols.append(eq_sol)
                print(f"a: {a}, b: {b}, c: {c}, eq solutions: {eqs_sols}")
        
        # a, b, c = 0.25, -3, 0.5
        def compute():
            a, b, c = 0.25, -3, 0.5
            err = 1e-7 # For float errors to be ironed out when converted to ints
            eqs_sols = []
            self.eqsSols = []
            for eqi, eq in enumerate(self.eqsList):
                eq_sol = eq(a, b, c)
                eqs_sols.append(eq_sol)
                self.eqsSols.append(int(eq_sol + err))
            print(f"a: {a}, b: {b}, c: {c}, eq solutions: {eqs_sols}")
            print(f"True eq sols: {self.eqsSols}")
            print(f"and sorted: {sorted(self.eqsSols)}")
        
        #compute()
        # Pre-computed...
        
        for i in range(len(self.eqsSols)):
            sol = self.eqsSols[i]
            pos = self.eqPos[i]
            r, c = pos
            self.grid.fillGrid(r, c, sol)
        
        
    def main(self):
        self.solve_equations()
        self.grid.printGrid()
        
        realGrid = self.narrowStart()
        realGrid.getFinalAns()
    
    def compareAllGrids(self, grids):
        uniqueGrids = []
        for grid in grids:
            same = False
            for grid0 in uniqueGrids:
                same = grid.isSameGrid(grid0)
                if same:
                    break
            if same:
                continue
            uniqueGrids.append(grid)
        
        return uniqueGrids
    
    def narrowStart(self):
        nominos12to13 = [
            Nomino([
                (0,0),
                (0,1),
                (0,2),
                (1,0),
                (2,0),
                (2,1), #
                (3,1),
                (2,-1),
                (2,-2),
                (3,-2),
                (4,-2),
                (5,-2)
            ]),
            Nomino([
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (2, 0),
                (3,0),  #
                (3, 1),
                (2, -1),
                (2, -2),
                (3, -2),
                (4, -2),
                (5, -2)
            ]),
            Nomino([
                (0, 0),
                (0, 1),
                (0, 2),
                (1,0), #
                (2,0), #
                (3,0), #
                (3,1),
                (4,1),
                (5,1),
                (5,0),
                (5,-1),
                (5,-2)
            ]),
            Nomino([
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),  #
                (2, 0),  #
                (2,1),  #
                (3, 1),
                (4, 1),
                (5, 1),
                (5, 0),
                (5, -1),
                (5, -2)
            ]),
        ]
        validGrids = []
        print("narrow start...")
        
        for nomino in nominos12to13:
            validGrids1 = self.placeNominoOnGrid(self.grid, nomino, 12)
            validGrids = validGrids + validGrids1
        validGrids = self.narrowMtoNUp(12, 13, validGrids)
        
        print(f"num valid grids now: {len(validGrids)} with 12 & 13")
        
        validGrids = self.narrowMtoNUp(13, 14, validGrids)
        print(f"num valid grids now: {len(validGrids)} with 12-14")
        
        validGrids = self.narrowMtoNUp(14, 15, validGrids)
        print(f"num valid grids now: {len(validGrids)} with 12-15")
        validGrids = self.narrowMtoNUp(15, 16, validGrids)
        print(f"num valid grids now: {len(validGrids)} with 12-16")
        
        validGrids = validGrids[1:]
        
        validGrids = self.narrowMtoNDown(12, 11, validGrids)
        print(f"num valid grids now: {len(validGrids)} with 11-16")
        validGrids = self.narrowMtoNDown(11,10, validGrids)
        print(f"num valid grids now: {len(validGrids)} with 10-16")
        
        for i in range(10, 2, -1):
            validGrids = self.narrowMtoNDown(i, i-1, validGrids)
            print(f"num valid grids now: {len(validGrids)} with {i-1}-16")
        
        validGrids = self.compareAllGrids(validGrids)
        print(f"After eliminating same grids, final grid count: {len(validGrids)}")
        
        for i, grid in enumerate(validGrids):
            grid.printGrid()
        
        print("ending the narrowing...")
        
        return validGrids[0]
    
    def narrowMtoNUp(self, n1, n2, grids):
        validGrids = []
        for grid in grids:
            nomino = grid.getNominoOfN(n1)
            
            ''' Add 1 to nomino in any way '''
            validSpaces = []
            newominos = []
            for alreadySpace in nomino.spaces:
                r, c = alreadySpace
                oDirs = [(-1, 0), (0, 1), (1, 0), (0, -1)]  # nesw
                for oDir in oDirs:
                    ro, co = oDir
                    potSpace = (r+ro, c+co)
                    if potSpace in nomino.spaces:
                        continue
                    if potSpace in validSpaces:
                        continue
                    validSpaces.append(potSpace)
                    
                    newSpaces = nomino.spaces + [potSpace]
                    newomino = Nomino(newSpaces)
                    newominos.append(newomino)
            ''''''
            validGrids2 = self.narrowMtoNUp1Grid(newominos, n1, n2, grid)
            validGrids = validGrids + validGrids2
        
        return validGrids
    
    def narrowMtoNDown(self, n1, n2, grids):
        validGrids = []
        for grid in grids:
            nomino = grid.getNominoOfN(n1, extras=False)
            # Can delete one if it doesn't break the connect and is not already a given space (original abc formulas)
            nominoSpaces = nomino.spaces
            for i, nominoSpace in enumerate(nominoSpaces):
                newSpaces = nominoSpaces[:i] + nominoSpaces[i + 1:]
                newomino = Nomino(newSpaces)
                isConnected = newomino.isConnected
                if not isConnected:
                    continue
                
                validGrids2 = self.placeNominoOnGrid(grid, newomino, n2)
                validGrids = validGrids + validGrids2
        
        return validGrids
            
    
    def narrowMtoNUp1Grid(self, nominos, n1, n2, grid):
        validGrids = []
        
        for nomino in nominos:
            someValidGrids = self.placeNominoOnGrid(grid, nomino, n2)
            validGrids = validGrids + someValidGrids
            
        return validGrids
    
    
    def placeNominoOnGrid(self, grid, nomino, placeN):
        validGrids = []
        
        if abs(placeN - nomino.n) > 1:
            return -1
        
        R, C = -1, -1
        for r in range(grid.r):
            for c in range(grid.c):
                n = grid.getGridSpace(r, c)
                if n == placeN:
                    R, C = r, c
                    break
            if R > -1:
                break
        
        for vari in nomino.variations:
            for originSpace in vari:
                newGrid = dc(grid)
                r0, c0 = originSpace
                fine = True
                for space in vari:
                    r1, c1 = space
                    r1, c1 = r1 + R - r0, c1 + C - c0
                    gridSpace = newGrid.getGridSpace(r1, c1)
                    
                    if gridSpace == 0 or gridSpace == placeN:
                        newGrid.fillGrid(r1, c1, placeN)
                        continue
                    # origin space is flawed, move to next
                    fine = False
                    break
                
                if not fine:
                    continue
                
                newomino = newGrid.getNominoOfN(placeN)
                if newomino.n > placeN:
                    continue
                if not newomino.isConnected:
                    continue
                
                validGrids.append(newGrid)
        
        return validGrids


s = Solver()
s.main()

''' Limit is 17 '''
def checkLimit():
    r, c = 13, 13
    nTiles = r*c
    print(f"nTiles: {nTiles}")
    
    sum = 0
    i = 1
    while True:
        if sum+i > nTiles:
            i -= 1
            break
        sum += i
        i += 1
        print(sum)
    
    print(f"limit: {i}; sum: {sum}")

#checkLimit()

testL = [(0,0), (0,1), (0,2), (1,1), (2,2)]
nom = Nomino(testL)
cc = nom.connectedCheck()
print(f"test: {cc}")
