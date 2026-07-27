'''
Author: Logan Richan
Completed Date: 2026-07-26
Look out for my internship application soon!

Final Output (before Solution class):
Grid 0:
Towers: 0,1,2,3,4,5,6,7,8,9,10,11,12,
[(0, 3), (2, 7), (1, 1), (2, 4), (1, 5), (4, 2), (3, 0), (3, 4), (5, 4), (5, 7), (7, 0), (6, 2), (7, 4)]
Matching Stamps: [0, 1, 16, 23, 528, 37, 88, 138, 272, 449, 750, 1100]
current position: (2, 7)
   x,   x,  33, 555,  14,  37,   x,1100,
  44, 541, 530,  53, 372,2712, 299, 113,
 489,   x, 410,  23, 112, 138,1047,59400,
 528, 572,  70, 164,2689, 335,   x, 264,
 127, 449,2667, 797,  16, 894,   x, 995,
 615, 750, 191,  88,   3, 272,   1,8976,
 219, 107,   1, 704, 845,  10, 944,   x,
   0, 659, 248,   x,7440,   x, 240,   5,
1 towerLists = [[(0, 3), (2, 7), (1, 1), (2, 4), (1, 5), (4, 2), (3, 0), (3, 4), (5, 4), (5, 7), (7, 0), (6, 2), (7, 4)]]
1 stampLists = [[0, 1, 16, 23, 528, 37, 88, 138, 272, 449, 750, 1100]]

...
Final Answer: 33609
'''
from copy import deepcopy as dcp

regions = {
    0: {(0,0), (0,1), (0,2), (0,3), (0,4)},
    1: {(0,5), (0,6), (0,7), (1,7), (2,7)},
    2: {(1,0), (1,1), (1,2), (2,0), (2,2)},
    3: {(1,3), (1,4), (2,3), (2,4), (2,5)},
    4: {(1,5), (1,6), (2,6), (3,6), (3,7)},
    5: {(2,1), (3,1), (3,2), (4,2), (5,2)},
    6: {(3,0), (4,0), (4,1), (5,0), (6,0)},
    7: {(3,3), (3,4), (4,3), (4,4)},
    8: {(3,5), (4,5), (4,6), (5,4), (5,5)},
    9: {(4,7), (5,7), (6,7), (7,7), (7,6)},
    10: {(5,1), (6,1), (7,1), (7,0), (7,2)},
    11: {(5,3), (6,3), (6,2), (6,4), (7,3)},
    12: {(5,6), (6,6), (6,5), (7,5), (7,4)}
}

coord_to_val = {
    (0,5): 37,
    (0,7): 1100,
    (2,3): 23,
    (2,5): 138,
    (3,0): 528,
    (4,1): 449,
    (4,4): 16,
    (5,1): 750,
    (5,3): 88,
    (5,5): 272,
    (5,6): 1,
    (7,0): 0
}

coord_to_region = {}
for region in regions:
    for coord in regions[region]:
        coord_to_region[coord] = region

def debug_check():
    print(coord_to_region)
    for i in range(8):
        for j in range(8):
            coords = (i, j)
            print(f"{coords} is in region {coord_to_region[coords]}" if coords in coord_to_region else f"{coords} not in region")
        print("")

basegrid = [[-1 for x in range(8)] for y in range(8)]
basegrid[7][0] = 0

class Grid:
    def __init__(self, startTop=False):
        self.grid = dcp(basegrid)
        self.top = startTop
        self.r = 7
        self.c = 0
        self.score = 0

        self.towers = [None for i in range(13)]
        self.nTowers = 0
        if self.top:
            self.towers[10] = (7,0)
            self.nTowers += 1

        self.scoreStamps = [0]
        self.stop = False

    def move(self, r, c, score, stamp=False):
        self.r, self.c, self.score = r, c, score
        self.grid[r][c] = score
        if stamp:
            self.scoreStamps.append(score)

    def printGrid(self, index=-1, printStamp=False):
        print(f"Grid {index}:")
        print(f"Towers: {''.join([f'{i},' if self.towers[i] else '' for i in range(len(self.towers))])}")
        if self.nTowers == 13:
            print(self.towers)
        if printStamp:
            print(f"Matching Stamps: {self.scoreStamps}")
        print(f"current position: {self.r, self.c}")
        for row in self.grid:
            print(''.join([f"{x:>4}," if x > -1 else "   x," for x in row]))


class Solver:
    def __init__(self):
        self.firstGrid = Grid()
        self.validGrids = [Grid(startTop=True), Grid(startTop=False)]
        self.newGrids = []

        self.movesBase = ((1, 2), (2, 1), (2, -1), (1, -2), (-1, -2), (-2, -1), (-2, 1), (-1, 2))
        self.movesSpec = ((0, 2), (2, 0), (0, -2), (-2, 0))

        self.N = 1
        self.K = 3 # Initial jump value

        self.k_ok = False
        self.KPhase = 0

        self.stoppedGrids = []

    def printValidGrids(self, printStamp=False, stopped=False):
        towerLists, stampLists = [], []
        grids = self.validGrids if not stopped else self.stoppedGrids
        for i, grid in enumerate(grids):
            grid.printGrid(index=i, printStamp=printStamp)
            if grid.towers not in towerLists:
                towerLists.append(grid.towers)
            if grid.scoreStamps not in stampLists:
                stampLists.append(grid.scoreStamps)
        print(f"{len(towerLists)} {towerLists = }")
        print(f"{len(stampLists)} {stampLists = }")

    def solve(self):
        for i in range(6):
            self.jumpK()
        print(f"After 18 moves, num valid grids: {len(self.validGrids)}; finding true K...\n")
        self.printValidGrids(printStamp=True)

        print("Now jumping by K=7...")
        self.K = 7
        for i in range(6):
            self.KPhase += 1
            self.jumpK()
            if self.KPhase >= 5:
                self.printValidGrids(printStamp=True)
            print(f"Jump {i+1}/6: {len(self.validGrids)} valid grids now")

        self.printValidGrids(printStamp=True, stopped=True)

    def jumpK(self):
        for i in range(self.K):
            self.k_ok = i + 1 == self.K
            if self.KPhase > 0:
                print(f"{self.KPhase = }; Jump #{i+1}/{self.K}")
            self.newGrids = []
            for grid in self.validGrids:
                if self.KPhase == 6 and grid.stop:
                    self.newGrids.append(grid)
                    continue
                self.tryAllMoves(grid)

            self.validGrids = self.newGrids
            self.N += 1

        self.newGrids = []
        for grid in self.validGrids:
            coord = (grid.r, grid.c)
            if coord not in coord_to_val:
                continue
            if coord_to_val[coord] != grid.score:
                continue
            self.newGrids.append(grid)
        self.validGrids = self.newGrids


    def tryAllMoves(self, grid):
        for move in self.movesBase:
            self.tryMove(grid, move, False)
        for move in self.movesSpec:
            self.tryMove(grid, move, True)

    def tryMove(self, grid, move, flip):
        dr, dc = move
        r2, c2 = grid.r + dr, grid.c + dc
        coord2 = (r2, c2)
        if r2 < 0 or r2 >= 8 or c2 < 0 or c2 >= 8:
            return False
        if grid.grid[r2][c2] > -1:
            return False
        if flip^grid.top:
            region = coord_to_region[coord2]
            if grid.towers[region]:
                return False

        if not self.k_ok and coord2 in coord_to_val:
            return False

        grid = dcp(grid)
        if not flip:
            grid.move(r2, c2, grid.score + self.N, stamp=self.k_ok)
            if grid.top:
                self.addTower(grid, region, coord2)
            self.newGrids.append(grid)
            return True
        if grid.top:
            if grid.score % self.N != 0:
                return False
            grid.move(r2, c2, grid.score // self.N, stamp=self.k_ok)
        else:
            grid.move(r2, c2, grid.score * self.N, stamp=self.k_ok)
            self.addTower(grid, region, coord2)

        grid.top = not grid.top
        self.newGrids.append(grid)
        return True

    def addTower(self, grid, region, coord):
        grid.towers[region] = coord
        grid.nTowers += 1
        if self.KPhase == 6 and grid.nTowers == 13 and not grid.stop:
            grid.stop = True
            self.stoppedGrids.append(grid)


solver = Solver(); solver.solve()

'''
Use this for analysis to help find K...
K found! 
K=7
'''
class Counting:
    def __init__(self, val=88, nStart=19, top=False, maxK=11):
        self.val = val
        self.ns = nStart
        self.top = top

        self.maxK = maxK


        self.cur = [(val, top, [val])]
        self.lookfor = {138: [], 272: [], 449: [], 750: [], 1100: []}

    def solve(self):
        k = 1
        for n in range(self.ns, self.ns+self.maxK):
            newcur = []
            for val, top, history in self.cur:
                toadd1a = val+n
                history.append(toadd1a)
                toadd1 = (toadd1a, top, history.copy())
                newcur.append(toadd1)
                if toadd1a in self.lookfor:
                    self.lookfor[toadd1a].append((k, top, history.copy()))

                if top and val % n != 0:
                    continue
                toadd2a = val * n if not top else val / n
                history[-1] = toadd2a
                toadd2 = (toadd2a, not top, history.copy())
                newcur.append(toadd2)
                if toadd2a in self.lookfor:
                    self.lookfor[toadd2a].append((k, not top, history.copy()))

            self.cur = newcur
            k += 1

        print(f"{self.val}:", self.lookfor)
        '''
        88: {138: [(7, False)], 272: [], 449: [], 750: [], 1100: []}
        [88, 107, 127, 2667, 2689, 2712, 113, 138] 3T
        
        138: {138: [], 272: [(7, False)], 449: [], 750: [], 1100: []}
        [138, 164, 191, 219, 248, 7440, 240.0, 272.0] 1T
        [272, 8976, 264.0, 299.0, 335.0, 372.0, 410.0, 449.0] 1T
        [449, 489, 530, 572, 615, 659, 704, 750]
        [750, 797, 845, 894, 944, 995, 1047, 1100]
        
        1 tower left after all numbers filled...
        '''
"""
counting = Counting(88, 19, False); counting.solve()
counting = Counting(138, 26, False); counting.solve()
counting = Counting(272, 26+7, False); counting.solve()
counting = Counting(449, 26+7*2); counting.solve()
counting = Counting(750, 26+7*3); counting.solve()
"""

class Solution:
    def __init__(self):
        # Taken from Solver output:
        self.grid = [
            [-1, -1, 33, 555, 14, 37, -1, 1100],
            [44, 541, 530, 53, 372, 2712, 299, 113],
            [489, -1, 410, 23, 112, 138, 1047, 59400],
            [528, 572, 70, 164, 2689, 335, -1, 264],
            [127, 449, 2667, 797, 16, 894, -1, 995],
            [615, 750, 191, 88, 3, 272, 1, 8976],
            [219, 107, 1, 704, 845, 10, 944, -1],
            [0, 659, 248, -1, 7440, -1, 240, 5]
        ]
        self.towers = [(0, 3), (2, 7), (1, 1), (2, 4), (1, 5), (4, 2), (3, 0), (3, 4), (5, 4), (5, 7), (7, 0), (6, 2), (7, 4)]

        # Checking stuff:
        self.moves = [(7, 0)]
        self.moveValues = [0]
        '''
        Final:
        self.moves = [(7, 0), (6, 2), (5, 4), (5, 6), (7, 7), (6, 5), (4, 4), (2, 4), (0, 4), (2, 3), (0, 2), (1, 0), (3, 0), (1, 1), (0, 3), (0, 5), (1, 3), (3, 2), (5, 3), (6, 1), (4, 0), (4, 2), (3, 4), (1, 5), (1, 7), (2, 5), (3, 3), (5, 2), (6, 0), (7, 2), (7, 4), (7, 6), (5, 5), (5, 7), (3, 7), (1, 6), (3, 5), (1, 4), (2, 2), (4, 1), (2, 0), (1, 2), (3, 1), (5, 0), (7, 1), (6, 3), (5, 1), (4, 3), (6, 4), (4, 5), (6, 6), (4, 7), (2, 6), (0, 7), (2, 7)]
        self.moveValues = [0, 1, 3, 1, 5, 10, 16, 112, 14, 23, 33, 44, 528, 541, 555, 37, 53, 70, 88, 107, 127, 2667, 2689, 2712, 113, 138, 164, 191, 219, 248, 7440, 240, 272, 8976, 264, 299, 335, 372, 410, 449, 489, 530, 572, 615, 659, 704, 750, 797, 845, 894, 944, 995, 1047, 1100, 59400]
        '''

        self.k1 = 3
        self.k2 = 7

        self.sum = 0

    def check(self):
        checkGrid = [[-1]*8 for _ in range(8)]
        checkGrid[7][0] = 0
        pos = (7, 0)
        top = True
        checkTowers = []

        countE, countF = 0, 0
        for r in range(8):
            for c in range(8):
                if self.grid[r][c] >= 0:
                    countF += 1
                    continue
                countE += 1
        print(f"Num empty spaces: {countE}; num filled spaces: {countF}")

        for n in range(1, countF):
            R, C = pos
            prev = self.moveValues[-1]
            valReg = prev + n
            valSpec = prev * n if not top else (prev / n if prev % n == 0 else None)
            for move in ((1, 2), (2, 1), (2, -1), (1, -2), (-1, -2), (-2, -1), (-2, 1), (-1, 2), (0, 2), (2, 0), (0, -2), (-2, 0)):
                dr, dc = move
                r2, c2 = R+dr, C+dc
                spec = dr == 0 or dc == 0
                val = valSpec if spec else valReg
                #print(f"trying move {n} at ({r2}, {c2}) with value {val}; {top = }; {spec = }; {prev = }")

                if r2 < 0 or r2 >= 8 or c2 < 0 or c2 >= 8 or checkGrid[r2][c2] > -1:
                    continue
                if top and spec and prev % n != 0:
                    continue

                if val == self.grid[r2][c2]:
                    if spec^top and (r2,c2) not in self.towers:
                        continue
                    if not spec^top and (r2,c2) in self.towers:
                        continue
                    if (r2, c2) in coord_to_val and ((n <= 18 and n % 3 != 0) or (n > 18 and (n-18)%7 != 0)):
                        continue
                    pos = (r2, c2)
                    if spec:
                        top = not top
                    if top:
                        checkTowers.append(pos)
                    self.moves.append(pos)
                    self.moveValues.append(int(val))
                    checkGrid[r2][c2] = val
                    break
            else:
                print(f"NOT GOOD? {self.moveValues = }, {n = }")
                return

        print(f"GOOD BOARD!")
        print(f"{self.moves = }\n{self.moveValues = }")

    def solve(self):
        print("\nFinal count...")
        for r, row in enumerate(self.grid):
            for c, val in enumerate(row):
                if val >= 0:
                    continue
                for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    r2, c2 = r+dr, c+dc
                    if r2 < 0 or r2 >= len(self.grid) or c2 < 0 or c2 >= len(self.grid[0]):
                        continue
                    toadd = self.grid[r2][c2]
                    if toadd < 0:
                        continue
                    print(f"from ({r}, {c}) adding {toadd}")
                    self.sum += toadd

        print(f"Final Answer: {self.sum}")

solution = Solution()
solution.check()
solution.solve()
