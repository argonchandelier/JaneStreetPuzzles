"""
Author:    Logan Richan
Completed: 2025-11-10
Made for solving the Jane Street Puzzle, Nov 2025
"""

'''
Things to consider:
-Arrow tiles are never box tiles
-Number tiles are always box tiles
-Tiles directly orthogonal to arrows that the arrow does not face are never box tiles
-Arrows can point to a tile past another arrow
-Number tiles within a row or column of an arrow give the upper bound on how far that arrow can point
-"Squeeze Theorem": arrows that point toward each other always have their closest box tile in between them,
 which also creates an upper bound
-Arrows pointing from different axes do not necessarily point to the same tile nor are the tiles necessarily closer than the
 intersecting tile
-Invalid tiles (never-box-tiles) may create a lower bound if orthogonal an arrow pointed at it, and can even create a bigger
 lower bound by chaining more invalid tiles in that direction
-Boundary is also the upper bound
-If number=#tiles valid, then all valid tiles are box tiles
-Numbers include their own cell + kings move
-Invalid can be determined by arrows and their lower bound
-Make arrows in binary: 0000=urbl (up, right, down, left), so 0 is no arrow and 1111=16 is arrow in every direction
-Base off of tile in the corner of the box (say, lower-bottom-left), and have the dimensions of the box (in the right
 order), then it folds only one possible way
-Number of box tiles has to factor into 2 and after that is a sum of 3 rectangular numbers (which still could be anything since
 1 is included)

Elimination done, now folding consideration:
-False: It is possible that all the uncertain tiles can be ignored for a check of using circles and squares to pick a corner & dimensions,
 then either they line up or they don't, so we can eliminate many corner/dimensions in that way
-Unknown: dimensions, corners, 0s on main grid
-Consider the minimum viable path; going thru minimum number of 0s to use as connecting paths, then maybe brute forcing to check
 for any overlaps, unpluggable holes, misaligned squares or circles that are inevitable on any folding
-...or at least which 0s can totally be avoided (before a full check, not that they are necessarily invalid) like spurs
-Find any paths that connect circles and/or squares first
-Remember that even 1 flip to certain can change a lot, as it anchors arrows which may even start a chain reaction with
 anchoring other arrows
-Go thru all possible paths before folding consideration? Might be a bit strenuous, even with chain reactions, although
 maybe there is another way to eliminate, like num tiles being specific "magic" numbers
 -Limit for a dimension is... 18

Method:
-Main grid should have a baseline for invalid/valid/certain tiles for being box tiles
-Other "arrow" grid keeping track of upper and lower bounds (or just store them as lists of x,y,lower,upper);
 store pointers
-Grid or list or some other way to track numbers with their constrictions as well

-For each arrow, find upper/lower bounds based on closest number in row/column, then squeeze/main grid boundary
 (remember to offset bounds based on if number is in the pointed direction or not)
-Mark all < lower bounds as invalid (including self), then also = lower bounds as invalid iff not in pointed direction
-Check each number to see if # valid tiles = number, and turn into certain tiles if so (if < number, throw an error)
-Repeat process
'''
from copy import deepcopy as dc

# From validCounts()
dims = {134: [[1, 3, 16]], 138: [[1, 4, 13], [1, 6, 9], [3, 3, 10]], 136: [[2, 2, 16], [2, 4, 10], [2, 6, 7], [3, 4, 8]], 132: [[2, 3, 12], [2, 5, 8]]}


class Arrow:
    def __init__(self, row, col, direction, lower=1, upper=19):
        self.r = row
        self.c = col
        self.dr = direction
        self.lower = lower
        self.upper = upper
        
        self.left = self.dr & 1 == 1
        self.down = (self.dr >> 1) & 1 == 1
        self.right = (self.dr >> 2) & 1 == 1
        self.up = (self.dr >> 3) & 1 == 1
        
        self.completed = False
    
    def SetAsComplete(self):
        self.completed = True
    
    def GetRowCol(self):
        return self.r, self.c
    
    def RestrictRangeFromNum(self, numInGrid):
        r2, c2, n = numInGrid
        if r2 != self.r or c2 != self.c:
            return
        dv = [r2-self.r, c2-self.c]
        d = dv[0]+dv[1]
        isPointingHere = self.IsPointingHere(dv)
        
        newMax = abs(d)
        if not isPointingHere:
            newMax -= 1
        
        if newMax < self.upper:
            self.upper = newMax
    
    def RestrictRangeFromBoundary(self):
        m1, m2, m3, m4 = [self.upper]*4
        if self.up:
            m1 = self.r
        if self.right:
            m2 = 19 - self.c
        if self.down:
            m3 = 19 - self.r
        if self.left:
            m4 = self.c
        
        self.upper = min(m1, m2, m3, m4)
        
    def IsPointingHere(self, displacement):
        dr, dc = displacement
        if self.up and dr < 0 and dc == 0:
            return True
        if self.right and dr == 0 and dc > 0:
            return True
        if self.down and dr > 0 and dc == 0:
            return True
        if self.left and dr == 0 and dc < 0:
            return True
        
        return False
        

class Grid:
    def __init__(self, setup=True):
        self.mainGrid = [[0] * 20 for _ in range(20)] # -1 invalid, 1 certain, 0 undetermined
        self.arrowPointerGrid = [[-1] * 20 for _ in range(20)]
        self.numsInGrid = [[-1] * 20 for _ in range(20)]
        
        self.countMin = 0
        self.countMax = 400
        
        self.invalidGrid = False
        
        # Each arrow: [r, c, dir (binary: urdl), lower, upper]
        self.arrows = [
            Arrow( 0, 1,0b0100,1,19),
            Arrow( 0, 8,0b0011,1,19),
            Arrow( 0,12,0b0010,1,19),
            Arrow( 0,15,0b0010,1,19),
            Arrow( 0,19,0b0011,1,19),
            Arrow( 1, 4,0b0100,1,19),
            Arrow( 1,13,0b0111,1,19),
            Arrow( 2, 6,0b0100,1,19),
            Arrow( 2, 9,0b0111,1,19),
            Arrow( 2,18,0b1011,1,19),
            Arrow( 3, 1,0b0010,1,19),
            Arrow( 3, 3,0b0110,1,19),
            Arrow( 4, 0,0b0100,1,19),
            Arrow( 4, 5,0b0110,1,19),
            Arrow( 4,15,0b1101,1,19),
            Arrow( 5,11,0b0110,1,19),
            Arrow( 6, 3,0b0010,1,19),
            Arrow( 6,10,0b0101,1,19),
            Arrow( 6,13,0b1011,1,19),
            Arrow( 6,16,0b0100,1,19),
            Arrow( 6,19,0b0011,1,19),
            Arrow( 7, 1,0b0110,1,19),
            Arrow( 7,14,0b0011,1,19),
            Arrow( 8, 5,0b1110,1,19),
            Arrow( 8, 9,0b1111,1,19),
            Arrow( 9, 1,0b1100,1,19),
            Arrow( 9,12,0b1101,1,19),
            Arrow( 9,16,0b0110,1,19),
            Arrow( 9,18,0b1011,1,19),
            Arrow(10, 1,0b1000,1,19),
            Arrow(10, 3,0b1010,1,19),
            Arrow(11, 8,0b1101,1,19),
            Arrow(11,10,0b0101,1,19),
            Arrow(11,14,0b1101,1,19),
            Arrow(11,17,0b1001,1,19),
            Arrow(12, 2,0b1100,1,19),
            Arrow(12, 5,0b1111,1,19),
            Arrow(12,18,0b1000,1,19),
            Arrow(13, 0,0b0100,1,19),
            Arrow(13,12,0b1001,1,19),
            Arrow(13,16,0b1000,1,19),
            Arrow(14, 8,0b1011,1,19),
            Arrow(14,11,0b1010,1,19),
            Arrow(14,13,0b0010,1,19),
            Arrow(15, 2,0b1100,1,19),
            Arrow(15, 4,0b1000,1,19),
            Arrow(15, 9,0b0011,1,19),
            Arrow(15,14,0b0011,1,19),
            Arrow(15,19,0b0001,1,19),
            Arrow(16, 5,0b1000,1,19),
            Arrow(16, 7,0b1100,1,19),
            Arrow(16,16,0b0001,1,19),
            Arrow(16,18,0b1001,1,19),
            Arrow(17, 1,0b0100,1,19),
            Arrow(17,10,0b1011,1,19),
            Arrow(18, 4,0b1100,1,19),
            Arrow(18, 6,0b0100,1,19),
            Arrow(18,14,0b1001,1,19),
            Arrow(19, 0,0b0100,1,19),
            Arrow(19, 4,0b1100,1,19),
            Arrow(19, 7,0b0100,1,19),
            Arrow(19,11,0b0001,1,19),
            Arrow(19,18,0b1001,1,19)
        ]
        
        # [row, col, num]
        self.numsInGridList = [
            [ 1,11, 4],
            [ 1,15, 4],
            [ 2, 7, 5],
            [ 3,12, 7],
            [ 3,14, 5],
            [ 4,10, 4],
            [ 4,13, 7],
            [ 4,17, 4],
            [ 5, 6, 4],
            [ 5, 8, 7],
            [ 6, 7, 9],
            [ 7,17, 6],
            [ 8, 2, 7],
            [ 8,11, 5],
            [ 9,14, 5],
            [10, 5, 4],
            [10, 7, 7],
            [10,18, 3],
            [13, 3, 5],
            [13, 6, 6],
            [13, 9, 2],
            [15, 6, 5],
            [17,12, 5],
            [17,13, 5],
            [18, 8, 4]
        ]
        
        self.circleTiles = [
            [1, 15],
            [3, 14],
            [4, 10],
            [4, 17],
            [8, 2],
            [17, 12]
        ]
        self.squareTiles = [
            [2, 7],
            [5, 8],
            [7, 17],
            [8, 11],
            [10, 18],
            [13, 3],
            [18, 8]
        ]
        
        self.showErrors = False
        
        if setup:
            for i, arrow in enumerate(self.arrows):
                r, c = arrow.GetRowCol()
                self.arrowPointerGrid[r][c] = i
            
            for numInCell in self.numsInGridList:
                r, c, n = numInCell[0], numInCell[1], numInCell[2]
                self.numsInGrid[r][c] = n
        
            print(f"Num arrows: {len(self.arrows)}")
            print(f"Num nums: {len(self.numsInGridList)}")
            
            self.changeInLoop = True
            
            # These were found manually
            self.setMainGrid(8, 6, 1)
            self.setMainGrid(1, 7, 1)
            self.setMainGrid(10, 13, -1)
            self.setMainGrid(8, 14, 1)
            self.setMainGrid(11, 4, -1)
            self.setMainGrid(9, 10, -1)
            
            self.FindInitBaseGrid()
            
            self.mainGridCounts()
            
    
    def GetNum(self, r, c):
        for numTile in self.numsInGridList:
            r2, c2, n = numTile
            if r != r2:
                continue
            if c != c2:
                continue
            
            return n
        return 0
    
    def compareGrid(self, otherGrid):
        for r in range(20):
            for c in range(20):
                vThis = self.getMainGridTile(r, c)
                vThat = otherGrid.getMainGridTile(r, c)
                if vThis != vThat:
                    return False
        return True
    
    def setGridAsInvalid(self):
        self.invalidGrid = True
    
    def mainGridCounts(self, prnt = True):
        counts = [0, 0, 0]
        for r in range(20):
            for c in range(20):
                val = self.getMainGridTile(r, c)
                counts[val+1] += 1
        
        cN1, c0, c1 = counts
        self.countMin = c1
        self.countMax = c0 + c1
        if prnt:
            print(f"\nCounts on box: invalid: {cN1}, valid range: {self.countMin}-{self.countMax}, unknown: {c0}")
        
        return {-1: cN1, 0: c0, 1: c1}
    
    def printGrid(self, grid=None):
        if grid is None:
            grid = self.mainGrid
        s = "\n"
        for r in range(len(grid)):
            s += "["
            for c in range(len(grid[0])):
                if c > 0:
                    s += ","
                val = grid[r][c]
                if val == 1:
                    s += "  1"
                elif val == 0:
                    s += "  ."
                elif val == -1:
                    s += "XXX"
                else:
                    s += "  " + str(val)
            s += "]\n"
        
        print(s)
        
    def getArrow(self, row, col):
        arrowIndex = self.arrowPointerGrid[row][col]
        if arrowIndex == -1:
            return None
        
        arrow = self.arrows[arrowIndex]
        return arrow
    
    def setMainGrid(self, row, col, val):
        if row < 0 or row > 19 or col < 0 or col > 19:
            return
        if val*-1 == self.getMainGridTile(row, col) and abs(val) > 0:
            if self.showErrors:
                print("ERROR, OVERWRITING VALIDITY THAT WAS ALREADY DETERMINED!")
            self.setGridAsInvalid()
            
        if self.getMainGridTile(row, col) != val:
            self.changeInLoop = True
        
        self.mainGrid[row][col] = val
    
    def getMainGridTile(self, row, col):
        if row < 0 or row > 19 or col < 0 or col > 19:
            if self.showErrors:
                print(f"WARNING: getting main grid tile from invalid location [{row},{col}]")
            return -2
        return self.mainGrid[row][col]
        
    def FindInitBaseGrid(self):
        ''' Num tiles immediately in box'''
        for numInCell in self.numsInGridList:
            r, c, n = numInCell
            self.setMainGrid(r, c, 1)
        ''' Arrow tiles (and orthogonal not-pointed-to tiles) immediately out of box'''
        for arrow in self.arrows:
            r, c = arrow.GetRowCol()
            self.setMainGrid(r, c, -1)
            if not arrow.up:
                self.setMainGrid(r-1, c, -1)
            if not arrow.right:
                self.setMainGrid(r, c+1, -1)
            if not arrow.down:
                self.setMainGrid(r+1, c, -1)
            if not arrow.left:
                self.setMainGrid(r, c-1, -1)
        
        # Restrict arrow range from numbers
        for arrow in self.arrows:
            for numInCell in self.numsInGridList:
                arrow.RestrictRangeFromNum(numInCell)
        
        # Restrict arrow range based on boundary
        for i, arrow in enumerate(self.arrows):
            arrow.RestrictRangeFromBoundary()
        
        # Restrict arrow range based on squeeze
        for i, arrow, in enumerate(self.arrows):
            #print(f"\nArrow {i} ([{arrow.r}, {arrow.c}]) bounds: upper: {arrow.upper}, lower: {arrow.lower}")
            for j, arrow2 in enumerate(self.arrows):
                if arrow is arrow2:
                    continue
                dv1 = [arrow2.r-arrow.r, arrow2.c-arrow.c]
                dv2 = [arrow.r-arrow2.r, arrow.c-arrow2.c]
                p1 = arrow.IsPointingHere(dv1)
                p2 = arrow2.IsPointingHere(dv2)
                if p1 and p2:
                    #print(f"SQUEEZE: Arrow2 {j} ([{arrow2.r}, {arrow2.c}])")
                    #print(f"dv1: {dv1}, dv2: {dv2}")
                    dist = abs(dv1[0] + dv1[1])
                    newMax = dist - 1
                    arrow.upper = min(arrow.upper, newMax)
                    arrow2.upper = min(arrow2.upper, newMax)
        
        self.RepeatGridRestriction()
    
    def RepeatGridRestriction(self):
        self.changeInLoop = True
        while self.changeInLoop:
            self.changeInLoop = False
            
            # Restrict range based on invalid tiles
            self.RestrictArrowRangeFromInvalidTiles()
            
            # Guaranteed tiles based on arrows
            self.ArrowGuaranteedTiles()
            if self.invalidGrid:
                return False
            
            # Guaranteed AND invalidates tiles based on nums
            self.LockNumCertainTiles()
            if self.invalidGrid:
                return False
            
            # Makes tiles invalid
            self.MazeCheck()
            if self.invalidGrid:
                return False
            
            # Confirms tiles
            self.MazeBlockingCheck()
            
            # Restrict arrow range
            self.RestrictArrowRangeFromConfirmedTiles()
        return True
    
    
    
    def RestrictArrowRangeFromInvalidTiles(self):
        iteration = 0
        while True:
            change = False
            for i, arrow in enumerate(self.arrows):
                '''
                Check the max squares and restrict max until all are valid
                Repeat for the mins and mark invalid in all dirs if moving min up
                '''
                r, c = arrow.GetRowCol()
                minTiles = []
                maxTiles = []
                minim = arrow.lower
                maxim = arrow.upper
                #print(f"Arrow {i} ([{arrow.r}, {arrow.c}]) bounds: upper: {arrow.upper}, lower: {arrow.lower}")
                arrowDirs = []
                if arrow.up:
                    minTiles.append([r - minim, c])
                    maxTiles.append([r - maxim, c])
                    arrowDirs.append([-1, 0])
                if arrow.right:
                    minTiles.append([r, c + minim])
                    maxTiles.append([r, c + maxim])
                    arrowDirs.append([0, 1])
                if arrow.down:
                    minTiles.append([r + minim, c])
                    maxTiles.append([r + maxim, c])
                    arrowDirs.append([1, 0])
                if arrow.left:
                    minTiles.append([r, c - minim])
                    maxTiles.append([r, c - maxim])
                    arrowDirs.append([0, -1])
                
                # Restrict upper
                while True:
                    restrict = False
                    for j, tile in enumerate(maxTiles):
                        r2, c2 = tile
                        if self.getMainGridTile(r2, c2) == -1:
                            restrict = True
                            change = True
                            self.changeInLoop = True
                            break
                    if not restrict:
                        break
                    
                    for j, tile in enumerate(maxTiles):
                        maxTiles[j] = [tile[0] - arrowDirs[j][0], tile[1] - arrowDirs[j][1]]
                    
                    arrow.upper -= 1
                
                # Restrict lower
                while True:
                    restrict = False
                    for j, tile in enumerate(minTiles):
                        r2, c2 = tile
                        if self.getMainGridTile(r2, c2) == -1:
                            restrict = True
                            change = True
                            self.changeInLoop = True
                            break
                    if not restrict:
                        break
                    
                    for j, tile in enumerate(minTiles):
                        r2, c2 = tile
                        self.setMainGrid(r2, c2, -1)  # Set invalid tile here (in arrow direction)
                        minTiles[j] = [tile[0] + arrowDirs[j][0], tile[1] + arrowDirs[j][1]]
                    
                    arrow.lower += 1
                    minim += 1
                    
                    # Also set invalid tile here (in not arrow direction)
                    if not arrow.up:
                        self.setMainGrid(r - minim, c, -1)
                    if not arrow.right:
                        self.setMainGrid(r, c + minim, -1)
                    if not arrow.down:
                        self.setMainGrid(r + minim, c, -1)
                    if not arrow.left:
                        self.setMainGrid(r, c - minim, -1)
            
            iteration += 1
            if not change:
                break
    
    def ArrowGuaranteedTiles(self):
        for i, arrow in enumerate(self.arrows):
            if arrow.upper < arrow.lower:
                if self.showErrors:
                    print(f"ERROR: LOWER AND UPPER BOUNDS YIELD NO SOLUTION... Arrow {i}: [{arrow.r}, {arrow.c}], "
                          f"lower: {arrow.lower}, upper: {arrow.upper}")
                self.setGridAsInvalid()
                return
            if arrow.upper == arrow.lower and not arrow.completed:
                arrow.SetAsComplete()
                d = arrow.upper
                r, c = arrow.GetRowCol()
                if arrow.up:
                    self.setMainGrid(r - d, c, 1)
                if arrow.right:
                    self.setMainGrid(r, c + d, 1)
                if arrow.down:
                    self.setMainGrid(r + d, c, 1)
                if arrow.left:
                    self.setMainGrid(r, c - d, 1)
                    
    def LockNumCertainTiles(self):
        disps = []
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                disps.append([dr, dc])
        for numInCell in self.numsInGridList:
            r, c, n = numInCell
            count = 0
            count2 = 0
            validCheckTiles = []
            uncCheckTiles = []
            for disp in disps:
                checkTile = [r+disp[0], c+disp[1]]
                r2, c2 = checkTile
                val = self.getMainGridTile(r2, c2)
                if val >= 0:
                    validCheckTiles.append(checkTile)
                    count += 1
                    if val == 0:
                        uncCheckTiles.append(checkTile)
                    if val == 1:
                        count2 += 1
            
            if count < n:
                if self.showErrors:
                    print(f"ERROR: NOT ENOUGH TILES IN NUMBER VICINITY (n={n}, num tiles={count}): [{r}, {c}]")
                self.setGridAsInvalid()
                return
            elif count == n:
                #print(f"Number tile is at VALID capacity, setting surrounding tiles to certain... (n={n}): [{r}, {c}]")
                for validCheckTile in validCheckTiles:
                    r2, c2 = validCheckTile
                    self.setMainGrid(r2, c2, 1)
            elif count2 == n:
                #print(f"Number tile is at ABSOLUTE capacity, setting surrounding uncertain to invalid (n={n}): [{r}, {c}]")
                for uncCheckTile in uncCheckTiles:
                    r2, c2 = uncCheckTile
                    self.setMainGrid(r2, c2, -1)
                
                    
    def MazeCheck(self, blockTile = (0, 0)): # [0,0] already is invalid, so default doesn't affect anything
        mazeGrid = [[0] * 20 for _ in range(20)]  # -1 wall (invalid), 0 untouched, 1 reached
        for r in range(20):
            for c in range(20):
                mainGridVal = self.getMainGridTile(r, c)
                if mainGridVal == -1:
                    mazeGrid[r][c] = -1
        
        sr, sc = 2, 7 # easy number tile to start off with
        mazeGrid[sr][sc] = 1
        blockR, blockC = blockTile
        mazeGrid[blockR][blockC] = -1
        
        oldTiles = [[sr, sc]]
        orthoDirs = [[-1, 0], [0, 1], [1, 0], [0, -1]]
        while True:
            newTiles = []
            for oldTile in oldTiles:
                for odir in orthoDirs:
                    potTile = [oldTile[0]+odir[0], oldTile[1]+odir[1]]
                    ptr, ptc = potTile
                    if ptr < 0 or ptr > 19 or ptc < 0 or ptc > 19:
                        continue
                    if mazeGrid[ptr][ptc] == 0:
                        mazeGrid[ptr][ptc] = 1
                        newTiles.append(potTile)
            
            if len(newTiles) == 0:
                break
            
            oldTiles = newTiles
        
        # All tiles connected to mainland are accounted for now...
        
        for r in range(20):
            for c in range(20):
                mainGridVal = self.getMainGridTile(r, c)
                mazeGridVal = mazeGrid[r][c]
                
                if mazeGridVal == 0: # couldn't be reached
                    if blockTile == (0, 0): # first maze
                        if mainGridVal == 1:
                            if self.showErrors:
                                print(f"ERROR: COULDN'T CONNECT ALL CONFIRMED TILES [{r}, {c}]")
                            self.setGridAsInvalid()
                            return
                            #continue
                        self.setMainGrid(r, c, -1) # Set to invalid tile
                    else: # block check maze
                        if mainGridVal == 1: # block causes disconnect in confirmed tiles, must be confirmed also
                            return False
        
        return True
    
    def MazeBlockingCheck(self):
        for r in range(20):
            for c in range(20):
                mainGridVal = self.getMainGridTile(r, c)
                if mainGridVal == 0:
                    blockTile = (r, c)
                    mazeBlockCheck = self.MazeCheck(blockTile)
                    if not mazeBlockCheck:
                        self.setMainGrid(r, c, 1)
                        
    def RestrictArrowRangeFromConfirmedTiles(self):
        # Restricting only maximum, so no making invalid tiles
        # outward
        orthoDirs = [[-1, 0], [0, 1], [1, 0], [0, -1]]
        for i, arrow in enumerate(self.arrows):
            checkTiles = []
            moveDirs = []
            if arrow.up:
                checkTiles.append([arrow.r-arrow.lower, arrow.c])
                moveDirs.append(orthoDirs[0])
            if arrow.right:
                checkTiles.append([arrow.r, arrow.c+arrow.lower])
                moveDirs.append(orthoDirs[1])
            if arrow.down:
                checkTiles.append([arrow.r+arrow.lower, arrow.c])
                moveDirs.append(orthoDirs[2])
            if arrow.left:
                checkTiles.append([arrow.r, arrow.c-arrow.lower])
                moveDirs.append(orthoDirs[3])
            
            def discoveryLoop(arrow):
                for dist in range(arrow.lower, arrow.upper):
                    for j, checkTile in enumerate(checkTiles):
                        r, c = checkTile
                        if self.getMainGridTile(r, c) == 1:
                            arrow.upper = dist
                            self.changeInLoop = True
                            return
                        checkTiles[j] = [checkTile[0] + moveDirs[j][0], checkTile[1] + moveDirs[j][1]]
            
            discoveryLoop(arrow)
            
    
    def DeadMazeCheck(self):
        mazeGrid = [[0] * 22 for _ in range(22)]  # -1 wall (invalid), 0 untouched, 1 reached
        for r in range(20):
            for c in range(20):
                mainGridVal = self.getMainGridTile(r, c)
                if mainGridVal != -1:
                    mazeGrid[r+1][c+1] = -1
        
        sr, sc = 0, 0  # easy number tile to start off with
        mazeGrid[sr][sc] = 1
        
        oldTiles = [[sr, sc]]
        orthoDirs = [[-1, 0], [0, 1], [1, 0], [0, -1]]
        while True:
            newTiles = []
            for oldTile in oldTiles:
                for odir in orthoDirs:
                    potTile = [oldTile[0] + odir[0], oldTile[1] + odir[1]]
                    ptr, ptc = potTile
                    if ptr < 0 or ptr > 21 or ptc < 0 or ptc > 21:
                        continue
                    if mazeGrid[ptr][ptc] == 0:
                        mazeGrid[ptr][ptc] = 1
                        newTiles.append(potTile)
            
            if len(newTiles) == 0:
                break
            
            oldTiles = newTiles
        
        for r in range(1, 21):
            for c in range(1, 21):
                if mazeGrid[r][c] == 0:
                    return False
        
        return True
    

class MeshGrid:
    def __init__(self, nRows, nCols, nLayers, sr, sc, gridI):
        self.nr = nRows
        self.nc = nCols
        self.nl = nLayers
        
        # Data storage purposes
        self.sr = sr
        self.sc = sc
        self.gridI = gridI
        
        self.frontGrid = [[0] * self.nc for _ in range(self.nr)]
        self.backGrid = [[0] * self.nc for _ in range(self.nr)]
        self.topGrid = [[0] * self.nc for _ in range(self.nl)]
        self.botGrid = [[0] * self.nc for _ in range(self.nl)]
        self.leftGrid = [[0] * self.nl for _ in range(self.nr)]
        self.rightGrid = [[0] * self.nl for _ in range(self.nr)]
        
        self.circles = []
        self.squares = []
        
        self.sideToGrid = {"front": self.frontGrid, "back": self.backGrid, "top": self.topGrid,
                           "bot": self.botGrid, "left": self.leftGrid, "right": self.rightGrid}
        
        self.sideToPlaceFunc = {"front": self.PlaceOnFrontGrid, "back": self.PlaceOnBackGrid, "top": self.PlaceOnTopGrid,
                                "bot": self.PlaceOnBotGrid, "left": self.PlaceOnLeftGrid, "right": self.PlaceOnRightGrid}
        
        self.nums = [0, 0, 0, 0, 0, 0]
        self.currentN = 0
    
    def PrintMeshGrid(self):
        print(f"\nGrid index: {self.gridI}\n"
              f"-------------------\n"
              f"Dimensions: [{self.nr}, {self.nc}, {self.nl}]\n"
              f"(19, 9) start at front pos: [{self.sr}, {self.sc}]")
        self.PrintFace("front")
        self.PrintFace("top")
        self.PrintFace("back")
        self.PrintFace("bot")
        self.PrintFace("left")
        self.PrintFace("right")
    
    def PrintFace(self, face):
        faceGrid = self.sideToGrid[face]
        s = f"{face}:\n"
        for r in range(len(faceGrid)):
            s += "["
            for c in range(len(faceGrid[0])):
                if c > 0:
                    s += ","
                val = faceGrid[r][c]
                if val == 1:
                    s += "  1"
                elif val == 2:
                    s += "op2"
                elif val == 3:
                    s += "nx3"
                else:
                    s += "  ."
            s += "]\n"
        print(s)
    
    def Place(self, side, r, c, t, o, n):
        self.currentN = n
        return self.sideToPlaceFunc[side](r, c, t, o)
    
    def PlaceType(self, side, r, c, t):
        if t == 2:
            self.circles.append([side, r, c])
        if t == 3:
            self.squares.append([side, r, c])
    
    def PlaceOnFrontGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnTopGrid(self.nl - 1, c, t, o)
        if r == self.nr:
            return self.PlaceOnBotGrid(0, c, t, o)
        if c == -1:
            return self.PlaceOnLeftGrid(r, self.nl - 1, t, o)
        if c == self.nc:
            return self.PlaceOnRightGrid(r, 0, t, o)
        
        if self.frontGrid[r][c] > 0:
            return None
        
        self.nums[0] += self.currentN
        self.frontGrid[r][c] = t
        self.PlaceType("front", r, c, t)
        return ["front", r, c, o%4]
    
    def PlaceOnBackGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnBotGrid(self.nl - 1, c, t, o)
        if r == self.nr:
            return self.PlaceOnTopGrid(0, c, t, o)
        if c == -1:
            return self.PlaceOnLeftGrid(self.nr - r - 1, 0, t, o+2)
        if c == self.nc:
            return self.PlaceOnRightGrid(self.nr - r - 1, self.nl - 1, t, o+2)
        
        if self.backGrid[r][c] > 0:
            return None
        
        self.nums[1] += self.currentN
        self.backGrid[r][c] = t
        self.PlaceType("back", r, c, t)
        return ["back", r, c, o%4]
    
    def PlaceOnTopGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnBackGrid(self.nr - 1, c, t, o)
        if r == self.nl:
            return self.PlaceOnFrontGrid(0, c, t, o)
        if c == -1:
            return self.PlaceOnLeftGrid(0, r, t, o-1)
        if c == self.nc:
            return self.PlaceOnRightGrid(0, self.nl - r - 1, t, o+1)
        
        if self.topGrid[r][c] > 0:
            return None
        
        self.nums[2] += self.currentN
        self.topGrid[r][c] = t
        self.PlaceType("top", r, c, t)
        return ["top", r, c, o%4]
    
    def PlaceOnBotGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnFrontGrid(self.nr - 1, c, t, o)
        if r == self.nl:
            return self.PlaceOnBackGrid(0, c, t, o)
        if c == -1:
            return self.PlaceOnLeftGrid(self.nr - 1, self.nl - r - 1, t, o+1)
        if c == self.nc:
            return self.PlaceOnRightGrid(self.nr - 1, r, t, o-1)
        
        if self.botGrid[r][c] > 0:
            return None
        
        self.nums[3] += self.currentN
        self.botGrid[r][c] = t
        self.PlaceType("bot", r, c, t)
        return ["bot", r, c, o%4]
    
    def PlaceOnLeftGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnTopGrid(c, 0, t, o+1)
        if r == self.nr:
            return self.PlaceOnBotGrid(self.nl - c - 1, 0, t, o-1)
        if c == -1:
            return self.PlaceOnBackGrid(self.nr - r - 1, 0, t, o+2)
        if c == self.nl:
            return self.PlaceOnFrontGrid(r, 0, t, o)
        
        if self.leftGrid[r][c] > 0:
            return None
        
        self.nums[4] += self.currentN
        self.leftGrid[r][c] = t
        self.PlaceType("left", r, c, t)
        return ["left", r, c, o%4]
    
    def PlaceOnRightGrid(self, r, c, t, o):
        if r == -1:
            return self.PlaceOnTopGrid(self.nl - c - 1, self.nc - 1, t, o-1)
        if r == self.nr:
            return self.PlaceOnBotGrid(c, self.nc - 1, t, o+1)
        if c == -1:
            return self.PlaceOnFrontGrid(r, self.nc - 1, t, o)
        if c == self.nl:
            return self.PlaceOnBackGrid(self.nr - r - 1, self.nc - 1, t, o+2)
        
        if self.rightGrid[r][c] > 0:
            return None
        
        self.nums[5] += self.currentN
        self.rightGrid[r][c] = t
        self.PlaceType("right", r, c, t)
        return ["right", r, c, o%4]
    
    
    def CheckCircles(self):
        for circle in self.circles:
            side, r, c = circle
            
            if side == "front":
                opposite = ["back", self.nr - r - 1, c]
            elif side == "top":
                opposite = ["bot", self.nr - r - 1, c]
            elif side == "left":
                opposite = ["right", r, self.nl - c - 1]
            else:
                continue
            # We don't need to check the other sides, because they already pair with these ones
            if opposite not in self.circles:
                return False
        
        return True
    
    def CheckSquares(self):
        for square in self.squares:
            side, r, c, = square
            grid = self.sideToGrid[side]
            nr, nc = len(grid), len(grid[0])
            checkTiles = [[r-1, c], [r+1, c], [r, c-1], [r, c+1]]
            goodTile = False
            for checkTile in checkTiles:
                r2, c2 = checkTile
                if r2 < 0 or c2 < 0 or r2 == nr or c2 == nc:
                    continue
                if [side, r2, c2] in self.squares:
                    goodTile = True
                    break
            if not goodTile:
                return False
        
        return True
    
    def CalcAnswer(self):
        ans = 1
        for numSum in self.nums:
            ans *= numSum
        
        print(f"Face sums: {self.nums}")
        print(f"Final answer: {ans}")
        return ans
        
            
class ShutTheBox:
    def __init__(self):
        self.grid = Grid()
        self.grid.printGrid()
        self.grids = [self.grid]
        
        self.fullDims = {}
        
        self.goodGridsToCheck = []
        self.goodMeshes = []
        
    def MainProcess1(self):
        oldCounts = self.grid.mainGridCounts()
        while True:
            newCounts = self.GridSplitBinaryTree()
            if oldCounts[0] == newCounts[0]:
                print("Counts didn't change, finally BREAKING THE LOOP")
                break
            oldCounts = newCounts
        
        self.GridSplitBinaryTree(split=True)
        
        self.CalcFullDims()
        
        self.SaveGoodGridsToFile()
    
    def CalcFullDims(self):
        self.fullDims = {}
        for nTiles in dims:
            newListForN = []
            dsListForN = dims[nTiles]
            for ds in dsListForN:
                d1, d2, d3 = ds
                newListForN += [[d1, d2, d3]]
                if d1 == d2 and d2 == d3:
                    continue
                newListForN += [[d2, d3, d1]]
                newListForN += [[d3, d1, d2]]
                if d2 != d3 and d1 != d2 and d1 != d3:
                    newListForN += [[d1, d3, d2]]
                    newListForN += [[d2, d1, d3]]
                    newListForN += [[d3, d2, d1]]
            
            self.fullDims[nTiles] = newListForN
        
        # For tests:
        self.fullDims[54] = [[3,3,3]]
        
        print(f"Full dims: {self.fullDims}")
        
    
    def GridSplitBinaryTree(self, split=False):
        oldGrids = self.grids
        ogGrid = self.grid
        
        for r in range(20):
            for c in range(20):
                if ogGrid.getMainGridTile(r, c) != 0:
                    continue
                newGrids = []
                invalidReached = False
                for grid in oldGrids:
                    if grid.getMainGridTile(r, c) != 0:
                        newGrids.append(grid)
                        continue
                    for val in [-1, 1]:
                        newGrid = dc(grid)
                        newGrid.setMainGrid(r, c, val)
                        
                        # Determine validity HERE; most intensive part probably
                        newGrid.RepeatGridRestriction()
                        
                        if not newGrid.invalidGrid:
                            newGrids.append(newGrid)
                        else:
                            invalidReached = True
                            
                if not invalidReached and not split:
                    pass
                    #print(f"Grid space [{r}, {c}] is still unknown for now...")
                else:
                    oldGrids = newGrids
                    print(f"Filled unknown grid space [{r}, {c}]; number of grids: {len(oldGrids)}")
                
        print(f"\nDONE with split; number of grids FINAL1: {len(oldGrids)}\n")
        if split:
            removeGrids = []
            for i, g in enumerate(oldGrids):
                count = g.mainGridCounts(prnt=False)
                if count[1] not in dims: # Note: I changed dims later because of the already-ran code after this
                    removeGrids.append(g)
            for g in removeGrids:
                oldGrids.remove(g)
            
            print(f"num of grids final2: {len(oldGrids)}\n")
            
            removeGrids = []
            for i, g in enumerate(oldGrids):
                valid = g.DeadMazeCheck()
                if not valid:
                    removeGrids.append(g)
            for g in removeGrids:
                oldGrids.remove(g)
            
            print(f"num of grids final3: {len(oldGrids)}\n")
            
            counts = []
            for i, g in enumerate(oldGrids):
                c = g.mainGridCounts(prnt=False)
                if c[1] not in counts:
                    counts.append(c[1])
                print(f"i: {i}")
                g.printGrid()
            counts.sort()
            print(f"FINAL: all possible num tiles: {counts}")
            
        self.grids = oldGrids
        counts = self.grids[0].mainGridCounts()
        return counts
    
    
    def SaveGoodGridsToFile(self):
        grids = self.goodGridsToCheck = self.grids
        with open("JSP-2025-11-saved-grids-LR.txt", "w") as file:
            for i, grid in enumerate(grids):
                file.write(f"Grid #{i}:\n")
                gridMeat = grid.mainGrid
                for r in range(len(gridMeat)):
                    for c in range(len(gridMeat[0])):
                        val = grid.getMainGridTile(r, c)
                        file.write(str(val+1))
                    file.write("\n")
            file.write(f"End; num grids: {len(grids)}")
    
    def LoadGoodGridsFromFile(self, prnt=False):
        self.goodGridsToCheck = []
        with open("JSP-2025-11-saved-grids-LR.txt", "r") as file:
            while True:
                line = file.readline()
                if line[0] == "G":
                    newGrid = Grid(setup=False)
                    for r in range(20):
                        line = file.readline()
                        for c in range(20):
                            char = line[c]
                            val = int(char)-1
                            newGrid.setMainGrid(r, c, val)
                    self.goodGridsToCheck.append(newGrid)
                    
                if line[0] == "E":
                    break
        
        if prnt:
            for i, grid in enumerate(self.goodGridsToCheck):
                print(f"Grid #{i}:")
                grid.printGrid()
                
    
    
    def BoxAGrid(self, grid, gridI, startTile=(19, 9)):
        nGoodMeshes = 0
        nTiles = grid.mainGridCounts(prnt=False)[1]
        possibleDims = self.fullDims[nTiles]
        for ds in possibleDims:
            # Number of rows, columns, and layers
            nr, nc, nl = ds
            for sr in range(nr):
                for sc in range(nc):
                    meshGrid = MeshGrid(nr, nc, nl, sr, sc, gridI)
                    meshGrid.PlaceOnFrontGrid(sr, sc, 1, 0)
                    mazeGrid = dc(grid.mainGrid)
                    
                    oldTilesMesh = [["front", sr, sc, 0]]
                    oldTilesFlat = [[startTile[0], startTile[1]]]
                    mazeGrid[startTile[0]][startTile[1]] = 2
                    badMesh = False
                    while True:
                        newTilesMesh = []
                        newTilesFlat = []
                        for i, oldTileMesh in enumerate(oldTilesMesh):
                            oldTileFlat = oldTilesFlat[i]
                            orthoDirs = [[-1, 0], [0, 1], [1, 0], [0, -1]]
                            for j, odir in enumerate(orthoDirs):
                                potTileFlat = [oldTileFlat[0] + odir[0], oldTileFlat[1] + odir[1]]
                                ptrF, ptcF = potTileFlat
                                
                                if ptrF < 0 or ptrF > 19 or ptcF < 0 or ptcF > 19:
                                    continue
                                if mazeGrid[ptrF][ptcF] != 1:
                                    continue
                                
                                n = grid.GetNum(ptrF, ptcF)

                                o = oldTileMesh[3]
                                odir2 = orthoDirs[(j+o)%4]
                                potTileMesh = [oldTileMesh[0], oldTileMesh[1] + odir2[0], oldTileMesh[2] + odir2[1], oldTileMesh[3]]
                                side, ptrM, ptcM, o = potTileMesh
                                
                                mazeGrid[ptrF][ptcF] = 2
                                
                                t = 1
                                baseTile = [ptrF, ptcF]
                                if baseTile in grid.circleTiles:
                                    t = 2
                                if baseTile in grid.squareTiles:
                                    t = 3
                                newTileMesh = meshGrid.Place(side, ptrM, ptcM, t, o, n)
                                if newTileMesh is None:
                                    badMesh = True
                                    break
                                newTilesMesh.append(newTileMesh)
                                newTilesFlat.append(potTileFlat)
                            if badMesh:
                                break

                        if len(newTilesMesh) == 0 or badMesh:
                            break
                        
                        oldTilesMesh = newTilesMesh
                        oldTilesFlat = newTilesFlat
                    
                    if badMesh:
                        continue
                        
                    if not meshGrid.CheckCircles():
                        continue
                    if not meshGrid.CheckSquares():
                        continue
                    
                    print("good mesh")
                    nGoodMeshes += 1
                    goodMesh = meshGrid
                    self.goodMeshes.append(goodMesh)
        
        print(f"nGoodMeshes in this grid net: {nGoodMeshes}")
    
    def BoxGrids(self):
        for i, grid in enumerate(self.goodGridsToCheck):
            print(f"Checking Grid #{i}")
            self.BoxAGrid(grid, i)
        
        print(f"\n\n\n\n!!!!!!!!!\nNum good meshes: {len(self.goodMeshes)}\n!!!!!!!!!")
        for goodMesh in self.goodMeshes:
            gridI = goodMesh.gridI
            print(f"Mesh good from grid #{gridI}")
            goodMesh.PrintMeshGrid()
            self.goodGridsToCheck[gridI].printGrid()
            goodMesh.CalcAnswer()
    
    def CompareAllGoodGrids(self):
        repeatGridIs = []
        for i, grid1 in enumerate(self.goodGridsToCheck):
            if i in repeatGridIs:
                continue
            for j, grid2 in enumerate(self.goodGridsToCheck[i+1:]):
                j += i+1
                if j in repeatGridIs:
                    continue
                same = grid1.compareGrid(grid2)
                if same:
                    repeatGridIs.append(j)
        if len(repeatGridIs) > 0:
            repeatGridIs.sort()
            print(f"Repeated grids is: {repeatGridIs}")
        else:
            print("No repeat grids")
            
    
    def MainProcess2(self):
        self.CalcFullDims()
        self.LoadGoodGridsFromFile()
        self.BoxGrids()
            
            
                    
    
        
        

stb = ShutTheBox()
# Only needs to be ran once:
stb.MainProcess1()
# Checks:
#stb.LoadGoodGridsFromFile(prnt=True)
#stb.CompareAllGoodGrids()
stb.MainProcess2()


def validCounts(mn, mx=None):
    if mx is None:
        mx = mn
    validDimensions = {}
    validNumbersOfTiles = []
    
    def nTilesBasedOnDim(d1, d2, d3):
        # exs: 3,4,5: 94; 4,4,4: 96
        nTiles = 2 * (d1*d2 + d2*d3 + d1*d3)
        return nTiles
    
    for d1 in range(1, 100):
        c2 = 0
        for d2 in range(d1, 100):
            c3 = 0
            #for d3 in range(d2, mx//4):
            for d3 in range(d2, 19):
                n = nTilesBasedOnDim(d1, d2, d3)
                
                if n > mx:
                    break
                c2 += 1
                c3 += 1
                if n < mn:
                    continue

                dimensions = [d1, d2, d3]
                if n not in validNumbersOfTiles:
                    validNumbersOfTiles.append(n)
                    validDimensions[n] = [dimensions]
                else:
                    validDimensions[n] = validDimensions[n] + [dimensions]
            if c3 == 0:
                break
        if c2 == 0:
            break
    
    validNumbersOfTiles.sort()
    print(validNumbersOfTiles)
    print(len(validDimensions), validDimensions)

#validCounts(132, 138)

'''
Face sums: [57, 28, 5, 11, 11, 17]
Final answer: 16414860
'''

